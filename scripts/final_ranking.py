#!/usr/bin/env python3
"""
final_ranking.py - 候选分子终选排序脚本
综合所有维度，输出最终候选分子优先级榜单
"""

import os
import sys
import json
import random
import numpy as np
import pandas as pd
import yaml
from pathlib import Path
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent.parent))

import warnings
warnings.filterwarnings('ignore')


class FinalRanking:
    """终选排序与报告生成"""

    def __init__(self, config_path=None):
        if config_path is None:
            config_path = Path(__file__).parent.parent / "config" / "config.yaml"
        
        with open(config_path, 'r', encoding='utf-8') as f:
            self.config = yaml.safe_load(f)
        
        self.targets = self.config['targets']
        self.weights = self.config['screening']['final_ranking']['weights']
        self.base_dir = Path(__file__).parent.parent
        
        self.target_colors = self.config['visualization']['target_colors']
    
    def load_all_data(self):
        """加载所有分析结果"""
        all_data = {}
        
        for target_name in self.targets:
            # 加载相互作用分析数据
            ia_path = self.base_dir / "results" / "alphafold3" / target_name / "interaction_analysis.csv"
            # 加载复合物预测数据
            af3_path = self.base_dir / "results" / "alphafold3" / target_name / "complex_predictions.csv"
            # 加载初始对接数据
            dd_path = self.base_dir / "results" / "diffdock" / target_name / "docking_results.csv"
            # 加载化合物性质
            cmpd_path = self.base_dir / "data" / "library" / "compounds_standardized.csv"
            
            data_sources = {}
            if ia_path.exists():
                data_sources['interaction'] = pd.read_csv(ia_path)
            if af3_path.exists():
                data_sources['af3'] = pd.read_csv(af3_path)
            if dd_path.exists():
                data_sources['docking'] = pd.read_csv(dd_path)
            if cmpd_path.exists():
                data_sources['compounds'] = pd.read_csv(cmpd_path)
            
            if data_sources:
                all_data[target_name] = data_sources
        
        return all_data
    
    def calculate_multi_dimension_score(self, row):
        """计算七维综合得分（含功能预测）"""
        np.random.seed(hash(row.get('mol_id', 'unknown')) % 2**32)
        
        # 1. 对接分数 (0-1)
        docking_score = row.get('best_confidence', row.get('composite_score', 0.5))
        docking_score_normalized = docking_score if isinstance(docking_score, (int, float)) and docking_score > 0 else 0.5
        
        # 2. 结构置信度 (基于pLDDT和ipTM)
        plddt = row.get('plddt', 70)
        iptm = row.get('iptm', 0.5)
        structure_confidence = (plddt / 100 * 0.5 + iptm * 0.5)
        
        # 3. 相互作用强度
        estimated_dG = row.get('estimated_dG', -5)
        dG_normalized = np.clip(abs(estimated_dG) / 15, 0, 1)
        key_residue_contacts = row.get('key_residue_contacts', 5)
        contact_normalized = min(key_residue_contacts / 15, 1.0)
        interaction_strength = dG_normalized * 0.6 + contact_normalized * 0.4
        
        # 4. 成药性
        qed = row.get('QED', row.get('drug_likeness_score', 0.5))
        drug_likeness = qed if isinstance(qed, (int, float)) else 0.5
        
        # 5. 结合位点匹配度
        binding_site_match = row.get('binding_efficiency_index', 1.0)
        site_score = min(binding_site_match / 3.0, 1.0) if isinstance(binding_site_match, (int, float)) else 0.5
        
        # 6. 🆕 竞争性抑制概率
        competitive_prob = row.get('competitive_probability', 0.5)
        competitive_score = competitive_prob if isinstance(competitive_prob, (int, float)) else 0.5
        
        # 7. 🆕 选择性指数（归一化）
        selectivity_index = row.get('selectivity_index', 1.5)
        selectivity_score = min(selectivity_index / 4.0, 1.0) if isinstance(selectivity_index, (int, float)) else 0.4
        
        # 综合七维加权得分
        final_score = (
            self.weights.get('docking_score', 0.20) * docking_score_normalized +
            self.weights.get('structure_confidence', 0.15) * structure_confidence +
            self.weights.get('interaction_strength', 0.20) * interaction_strength +
            self.weights.get('drug_likeness', 0.10) * drug_likeness +
            self.weights.get('binding_site_match', 0.10) * site_score +
            self.weights.get('competitive_inhibition', 0.15) * competitive_score +
            self.weights.get('selectivity_index', 0.10) * selectivity_score
        )
        
        return round(final_score, 4), {
            'docking_score': round(docking_score_normalized, 4),
            'structure_confidence': round(structure_confidence, 4),
            'interaction_strength': round(interaction_strength, 4),
            'drug_likeness': round(drug_likeness, 4),
            'binding_site_match': round(site_score, 4),
            'competitive_inhibition': round(competitive_score, 4),
            'selectivity_index': round(selectivity_score, 4),
        }
    
    def rank_target(self, target_name, data_sources):
        """对单个靶点进行终选排序"""
        print(f"\n{'=' * 60}")
        print(f"终选排序：{target_name}")
        print(f"{'=' * 60}")
        
        # 合并所有数据源
        df_main = None
        
        if 'interaction' in data_sources:
            df_main = data_sources['interaction'].copy()
        elif 'af3' in data_sources:
            df_main = data_sources['af3'].copy()
        
        if df_main is None:
            print(f"  警告：{target_name} 没有可用的分析数据")
            return None, None
        
        # 计算综合得分
        final_scores = []
        score_breakdowns = []
        
        for _, row in df_main.iterrows():
            final_score, breakdown = self.calculate_multi_dimension_score(row)
            final_scores.append(final_score)
            score_breakdowns.append(breakdown)
        
        df_main['final_score'] = final_scores
        df_main['score_breakdown'] = score_breakdowns
        
        # 排序
        df_ranked = df_main.sort_values('final_score', ascending=False)
        df_ranked['rank'] = range(1, len(df_ranked) + 1)
        df_ranked['target'] = target_name
        
        # 分级
        def classify(score):
            if score >= 0.75:
                return 'A (高优先级)'
            elif score >= 0.6:
                return 'B (中优先级)'
            elif score >= 0.45:
                return 'C (低优先级)'
            else:
                return 'D (不推荐)'
        
        df_ranked['priority_class'] = df_ranked['final_score'].apply(classify)
        
        # 统计
        print(f"排名分子数：{len(df_ranked)}")
        print(f"最高得分：{df_ranked['final_score'].max():.4f}")
        print(f"平均得分：{df_ranked['final_score'].mean():.4f}")
        print(f"优先级分布：")
        for cls in ['A (高优先级)', 'B (中优先级)', 'C (低优先级)', 'D (不推荐)']:
            cnt = (df_ranked['priority_class'] == cls).sum()
            print(f"  {cls}: {cnt} 个")
        
        # 保存
        output_dir = self.base_dir / "results" / "final_report" / target_name
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 保存不含breakdown的CSV
        df_save = df_ranked.drop(columns=['score_breakdown'], errors='ignore')
        df_save.to_csv(output_dir / "final_ranking.csv", index=False)
        
        # 保存完整JSON
        ranked_data = df_ranked.to_dict('records')
        with open(output_dir / "final_ranking.json", 'w', encoding='utf-8') as f:
            json.dump(ranked_data, f, indent=2, ensure_ascii=False, default=str)
        
        return df_ranked, ranked_data
    
    def generate_final_report(self, all_rankings):
        """生成最终汇总报告"""
        print(f"\n{'=' * 60}")
        print("生成终选汇总报告")
        print(f"{'=' * 60}")
        
        all_top = []
        report_sections = []
        
        for target_name, (df_ranked, _) in all_rankings.items():
            if df_ranked is None:
                continue
            
            # 取每个靶点的Top10
            df_top10 = df_ranked.head(10)
            all_top.append(df_top10)
            
            # 生成靶点摘要
            a_count = (df_ranked['priority_class'] == 'A (高优先级)').sum()
            b_count = (df_ranked['priority_class'] == 'B (中优先级)').sum()
            
            section = {
                'target': target_name,
                'description': self.targets[target_name]['description'],
                'total_ranked': len(df_ranked),
                'class_a': int(a_count),
                'class_b': int(b_count),
                'top_score': float(df_ranked['final_score'].max()),
                'top_candidates': df_top10[['rank', 'mol_id', 'final_score', 'priority_class']].to_dict('records')
            }
            report_sections.append(section)
            
            print(f"  {target_name}: 共{len(df_ranked)}个, "
                  f"A类{a_count}个, B类{b_count}个, "
                  f"最高得分={df_ranked['final_score'].max():.4f}")
        
        # 合并所有Top候选
        if all_top:
            df_all_top = pd.concat(all_top, ignore_index=True)
            output_dir = self.base_dir / "results" / "final_report"
            df_all_top.to_csv(output_dir / "all_top_candidates.csv", index=False)
        
        # 保存报告JSON
        final_report = {
            'generated_at': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'total_candidates_across_targets': sum(s['total_ranked'] for s in report_sections),
            'target_summaries': report_sections,
            'screening_weights': self.weights,
            'recommendations': self.generate_recommendations(report_sections)
        }
        
        with open(output_dir / "final_report.json", 'w', encoding='utf-8') as f:
            json.dump(final_report, f, indent=2, ensure_ascii=False)
        
        # 可视化数据
        viz_data = {
            'targets': list(all_rankings.keys()),
            'final_scores': {},
            'priority_distributions': {},
            'score_breakdowns': {},
            'top_candidates': {}
        }
        
        for target_name, (df_ranked, _) in all_rankings.items():
            if df_ranked is None:
                continue
            viz_data['final_scores'][target_name] = df_ranked['final_score'].tolist()
            viz_data['priority_distributions'][target_name] = df_ranked['priority_class'].value_counts().to_dict()
            
            # 各维度平均得分
            breakdowns = df_ranked['score_breakdown'].tolist()
            if breakdowns:
                avg_breakdown = {}
                for key in breakdowns[0].keys():
                    avg_breakdown[key] = round(np.mean([b[key] for b in breakdowns]), 4)
                viz_data['score_breakdowns'][target_name] = avg_breakdown
            
            # Top5 候选
            viz_data['top_candidates'][target_name] = df_ranked.head(5)[['rank', 'mol_id', 'final_score', 'priority_class']].to_dict('records')
        
        viz_path = self.base_dir / "results" / "final_report" / "visualization_data.json"
        with open(viz_path, 'w', encoding='utf-8') as f:
            json.dump(viz_data, f, indent=2, ensure_ascii=False)
        
        print(f"\n最终报告已保存：{self.base_dir / 'results' / 'final_report' / 'final_report.json'}")
        
        return final_report
    
    def generate_recommendations(self, report_sections):
        """生成后续实验建议"""
        recommendations = []
        
        for section in report_sections:
            target = section['target']
            if section['class_a'] > 0:
                recommendations.append({
                    'target': target,
                    'priority': '高',
                    'action': f'推荐对{target}的{section["class_a"]}个A类候选分子进行体外结合实验（SPR/BLI）',
                    'assay': '表面等离子体共振 (SPR) 或 生物层干涉 (BLI)',
                    'expected_Kd_range': '< 1 μM'
                })
            
            if section['class_b'] > 0:
                recommendations.append({
                    'target': target,
                    'priority': '中',
                    'action': f'对{target}的{section["class_b"]}个B类候选进行细胞水平功能验证',
                    'assay': 'NFAT-luciferase 报告基因实验 或 混合淋巴细胞反应 (MLR)',
                    'expected_IC50_range': '1-10 μM'
                })
        
        if not recommendations:
            recommendations.append({
                'target': '所有靶点',
                'priority': '一般',
                'action': '当前候选分子评分偏低，建议扩大化合物库或优化筛选条件后重新筛选',
                'note': '可尝试路线二进行全新分子生成'
            })
        
        return recommendations
    
    def run(self):
        """运行完整终选排序流程"""
        print("\n" + "█" * 60)
        print("█" + "候选分子终选排序与报告生成模块".center(58) + "█")
        print("█" + f"启动时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}".center(58) + "█")
        print("█" * 60)
        
        print(f"\n排序权重配置：")
        for key, val in self.weights.items():
            print(f"  {key}: {val}")
        
        # 加载所有数据（如果有的话），否则生成模拟数据
        all_data = self.load_all_data()
        
        if not all_data:
            print("\n未找到已有数据，运行完整模拟管线...")
            from scripts.interaction_analysis import InteractionAnalyzer
            analyzer = InteractionAnalyzer()
            ia_results = analyzer.run()
            all_data = self.load_all_data()
        
        # 对每个靶点排序
        all_rankings = {}
        for target_name in self.targets:
            if target_name in all_data:
                df_ranked, ranked_data = self.rank_target(target_name, all_data[target_name])
                all_rankings[target_name] = (df_ranked, ranked_data)
            else:
                # 模拟排序（如果数据不足）
                print(f"\n{target_name}: 使用模拟数据排序...")
                df_ranked, ranked_data = self.rank_target(target_name, {})
                all_rankings[target_name] = (df_ranked, ranked_data)
        
        # 生成最终报告
        final_report = self.generate_final_report(all_rankings)
        
        print("\n" + "█" * 60)
        print("█" + "终选排序与报告生成全部完成！".center(58) + "█")
        print("█" * 60)
        
        return {
            'rankings': all_rankings,
            'final_report': final_report
        }


if __name__ == "__main__":
    ranker = FinalRanking()
    results = ranker.run()