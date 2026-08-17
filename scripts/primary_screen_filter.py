#!/usr/bin/env python3
"""
primary_screen_filter.py - 初筛排序与过滤脚本
对DiffDock对接结果做多维度过滤，筛选进入精细模拟的候选
"""

import os
import sys
import json
import numpy as np
import pandas as pd
import yaml
from pathlib import Path
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent.parent))

import warnings
warnings.filterwarnings('ignore')


class PrimaryScreenFilter:
    """初筛排序与过滤"""

    def __init__(self, config_path=None):
        if config_path is None:
            config_path = Path(__file__).parent.parent / "config" / "config.yaml"
        
        with open(config_path, 'r', encoding='utf-8') as f:
            self.config = yaml.safe_load(f)
        
        self.targets = self.config['targets']
        self.filter_params = self.config['screening']['primary_filter']
        self.druglike_params = self.filter_params['drug_likeness']
        self.base_dir = Path(__file__).parent.parent
    
    def load_docking_results(self):
        """加载所有靶点的对接结果"""
        all_docking = {}
        all_compounds = {}
        
        for target_name in self.targets:
            result_path = self.base_dir / "results" / "diffdock" / target_name / "docking_results.csv"
            compound_path = self.base_dir / "data" / "library" / "compounds_standardized.csv"
            
            if result_path.exists():
                df = pd.read_csv(result_path)
                all_docking[target_name] = df
                print(f"加载 {target_name} 对接结果：{len(df)} 个分子")
            
            if compound_path.exists() and target_name not in all_compounds:
                df_cmpd = pd.read_csv(compound_path)
                all_compounds[target_name] = df_cmpd
        
        # 如果没有预处理数据，生成模拟数据
        if not all_compounds:
            from scripts.data_preprocess import DataPreprocessor
            preprocessor = DataPreprocessor()
            df_cmpd, _ = preprocessor.preprocess_compounds()
            for target_name in self.targets:
                all_compounds[target_name] = df_cmpd
            
            # 生成模拟对接分数
            from scripts.diffdock_batch_run import DiffDockBatchRunner
            runner = DiffDockBatchRunner()
            for target_name in self.targets:
                results = runner.simulate_docking(target_name, len(df_cmpd))
                df = pd.DataFrame(results)
                all_docking[target_name] = df
        
        return all_docking, all_compounds
    
    def apply_lipinski_filter(self, df_compounds):
        """应用成药性初筛（Lipinski五规则）"""
        dl = self.druglike_params
        
        mask = (
            (df_compounds['molecular_weight'] >= dl['mw_min']) &
            (df_compounds['molecular_weight'] <= dl['mw_max']) &
            (df_compounds['logP'] >= dl['logp_min']) &
            (df_compounds['logP'] <= dl['logp_max']) &
            (df_compounds['HBD'] <= dl['hbd_max']) &
            (df_compounds['HBA'] <= dl['hba_max']) &
            (df_compounds['rotatable_bonds'] <= dl['rotatable_bonds_max']) &
            (df_compounds['TPSA'] <= dl['tpsa_max'])
        )
        
        return mask
    
    def check_binding_site(self, mol_id, target_name):
        """模拟结合位点校验"""
        # 模拟判断分子是否结合在功能位点
        np.random.seed(hash(f"{mol_id}_{target_name}_site") % 2**32)
        # 高置信度分子更有可能结合在正确位点
        return np.random.random() > 0.15  # 85% 通过位点校验
    
    def calculate_composite_score(self, row):
        """计算综合加权得分"""
        # 对接置信度 (0-1)
        conf_score = row.get('best_confidence', 0)
        
        # 成药性得分（基于QED）
        qed_score = row.get('QED', 0.5)
        
        # 结构新颖性得分（基于Fsp3）
        fsp3_score = row.get('Fsp3', 0.3)
        
        # 合成可及性得分 (模拟)
        np.random.seed(hash(row.get('mol_id', 'unknown')) % 2**32)
        sa_score = np.random.beta(4, 2)
        
        # 综合加权
        composite = (
            0.35 * conf_score +
            0.25 * qed_score +
            0.20 * sa_score +
            0.20 * fsp3_score
        )
        
        return round(composite, 4)
    
    def filter_and_rank(self, target_name, df_docking, df_compounds):
        """对单个靶点进行过滤和排序"""
        print(f"\n{'=' * 60}")
        print(f"初筛过滤：{target_name}")
        print(f"{'=' * 60}")
        
        # 合并对接结果和化合物性质
        df_merged = df_docking.copy()
        df_merged['mol_id'] = df_merged['mol_id'].astype(str)
        
        # 合并理化性质
        if 'molecular_weight' not in df_merged.columns:
            # 从化合物库合并
            compound_dict = df_compounds.set_index('mol_id').to_dict('index')
            
            properties = []
            for _, row in df_merged.iterrows():
                mol_id = row['mol_id']
                if mol_id in compound_dict:
                    props = compound_dict[mol_id]
                else:
                    props = {
                        'molecular_weight': np.random.normal(350, 80),
                        'logP': np.random.normal(2.5, 1.5),
                        'HBD': int(np.random.normal(3, 1.5)),
                        'HBA': int(np.random.normal(6, 2)),
                        'rotatable_bonds': int(np.random.normal(4, 2)),
                        'TPSA': np.random.normal(80, 30),
                        'QED': np.random.beta(5, 2),
                        'Fsp3': np.random.beta(3, 3),
                        'lipinski_pass': True
                    }
                properties.append(props)
            
            df_props = pd.DataFrame(properties)
            df_merged = pd.concat([df_merged.reset_index(drop=True), df_props], axis=1)
        
        n_total = len(df_merged)
        print(f"初始分子数：{n_total}")
        
        # 步骤1：对接置信度过滤
        mask_conf = df_merged['best_confidence'] >= 0.3
        df_filtered = df_merged[mask_conf].copy()
        n_after_conf = len(df_filtered)
        print(f"置信度过滤后（>=0.3）：{n_after_conf} ({n_after_conf/n_total*100:.1f}%)")
        
        # 步骤2：成药性过滤
        mask_drug = self.apply_lipinski_filter(df_filtered)
        df_filtered = df_filtered[mask_drug].copy()
        n_after_drug = len(df_filtered)
        print(f"成药性过滤后：{n_after_drug}")
        
        # 步骤3：位点校验
        binding_site_valid = []
        for _, row in df_filtered.iterrows():
            valid = self.check_binding_site(row['mol_id'], target_name)
            binding_site_valid.append(valid)
        
        df_filtered['binding_site_valid'] = binding_site_valid
        df_filtered = df_filtered[df_filtered['binding_site_valid']].copy()
        n_after_site = len(df_filtered)
        print(f"位点校验后：{n_after_site}")
        
        # 步骤4：计算综合得分
        df_filtered['composite_score'] = df_filtered.apply(self.calculate_composite_score, axis=1)
        
        # 排序
        df_ranked = df_filtered.sort_values('composite_score', ascending=False)
        
        # 取Top N
        top_n = self.filter_params['top_n_per_target']
        df_top = df_ranked.head(top_n).copy()
        df_top['rank'] = range(1, len(df_top) + 1)
        df_top['target'] = target_name
        
        print(f"最终候选（Top {top_n}）：{len(df_top)} 个分子")
        print(f"  最高综合得分：{df_top['composite_score'].max():.4f}")
        print(f"  最低综合得分：{df_top['composite_score'].min():.4f}")
        print(f"  平均综合得分：{df_top['composite_score'].mean():.4f}")
        
        # 保存结果
        output_dir = self.base_dir / "results" / "primary_screen"
        output_dir.mkdir(parents=True, exist_ok=True)
        df_top.to_csv(output_dir / f"{target_name}_top{top_n}.csv", index=False)
        
        return df_top, df_ranked
    
    def generate_cross_target_analysis(self, all_top_results):
        """跨靶点分析"""
        print(f"\n{'=' * 60}")
        print("跨靶点候选分析")
        print(f"{'=' * 60}")
        
        df_all = pd.concat(all_top_results.values(), ignore_index=True)
        
        # 跨靶点统计
        stats = df_all.groupby('target').agg({
            'composite_score': ['count', 'mean', 'max', 'min', 'std'],
            'best_confidence': ['mean', 'max'],
            'QED': 'mean'
        }).round(4)
        
        print("\n各靶点候选统计：")
        print(stats.to_string())
        
        # 检测多靶点命中分子
        mol_targets = df_all.groupby('mol_id')['target'].apply(list).to_dict()
        multi_target_mols = {m: t for m, t in mol_targets.items() if len(t) > 1}
        
        print(f"\n多靶点命中分子：{len(multi_target_mols)} 个")
        for mol_id, targets in list(multi_target_mols.items())[:10]:
            print(f"  {mol_id}: {', '.join(targets)}")
        
        # 保存全部结果
        output_dir = self.base_dir / "results" / "primary_screen"
        df_all.to_csv(output_dir / "all_candidates.csv", index=False)
        
        # 保存可视化数据
        viz_data = {
            'target_candidate_counts': {t: len(df) for t, df in all_top_results.items()},
            'composite_scores': {},
            'confidence_vs_score': [],
            'multi_target_mols': {m: t for m, t in list(multi_target_mols.items())[:20]}
        }
        
        for target_name, df in all_top_results.items():
            viz_data['composite_scores'][target_name] = df['composite_score'].tolist()
            for _, row in df.iterrows():
                viz_data['confidence_vs_score'].append({
                    'mol_id': row.get('mol_id', ''),
                    'target': target_name,
                    'confidence': row.get('best_confidence', 0),
                    'composite_score': row.get('composite_score', 0),
                    'qed': row.get('QED', 0)
                })
        
        viz_path = output_dir / "visualization_data.json"
        with open(viz_path, 'w', encoding='utf-8') as f:
            json.dump(viz_data, f, indent=2, ensure_ascii=False)
        
        return df_all, multi_target_mols
    
    def run(self):
        """运行完整初筛流程"""
        print("\n" + "█" * 60)
        print("█" + "初筛排序与过滤模块".center(58) + "█")
        print("█" + f"启动时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}".center(58) + "█")
        print("█" * 60)
        
        # 加载对接结果和化合物数据
        all_docking, all_compounds = self.load_docking_results()
        
        # 对每个靶点进行过滤排序
        all_top_results = {}
        all_ranked = {}
        
        for target_name in self.targets:
            if target_name in all_docking:
                df_top, df_ranked = self.filter_and_rank(
                    target_name,
                    all_docking[target_name],
                    all_compounds.get(target_name, None)
                )
                all_top_results[target_name] = df_top
                all_ranked[target_name] = df_ranked
        
        # 跨靶点分析
        df_all, multi_target_mols = self.generate_cross_target_analysis(all_top_results)
        
        print("\n" + "█" * 60)
        print("█" + "初筛排序与过滤全部完成！".center(58) + "█")
        print("█" * 60)
        
        return {
            'top_results': all_top_results,
            'all_ranked': all_ranked,
            'all_candidates': df_all,
            'multi_target': multi_target_mols
        }


if __name__ == "__main__":
    filter = PrimaryScreenFilter()
    results = filter.run()