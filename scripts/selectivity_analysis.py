#!/usr/bin/env python3
"""
selectivity_analysis.py - 选择性预测分析模块
评估候选药物对目标靶点的选择性，避免脱靶效应
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
warnings = __import__('warnings')
warnings.filterwarnings('ignore')


class SelectivityAnalyzer:
    """跨靶点选择性分析"""

    def __init__(self, config_path=None):
        if config_path is None:
            config_path = Path(__file__).parent.parent / "config" / "config.yaml"
        with open(config_path, 'r', encoding='utf-8') as f:
            self.config = yaml.safe_load(f)
        self.targets = list(self.config['targets'].keys())
        self.base_dir = Path(__file__).parent.parent
        
        # 各靶点功能域相似度矩阵（基于IgV域保守性）
        # 值越大越可能交叉结合
        self.domain_similarity = {
            ('PD-1', 'PD-1'): 1.0, ('PD-1', 'LAG-3'): 0.25, ('PD-1', 'TIM-3'): 0.35, ('PD-1', 'VISTA'): 0.20,
            ('LAG-3', 'LAG-3'): 1.0, ('LAG-3', 'PD-1'): 0.25, ('LAG-3', 'TIM-3'): 0.20, ('LAG-3', 'VISTA'): 0.15,
            ('TIM-3', 'TIM-3'): 1.0, ('TIM-3', 'PD-1'): 0.35, ('TIM-3', 'LAG-3'): 0.20, ('TIM-3', 'VISTA'): 0.18,
            ('VISTA', 'VISTA'): 1.0, ('VISTA', 'PD-1'): 0.20, ('VISTA', 'LAG-3'): 0.15, ('VISTA', 'TIM-3'): 0.18,
        }

    def load_all_docking(self):
        """加载所有靶点的对接结果"""
        all_data = {}
        for target_name in self.targets:
            path = self.base_dir / "results" / "diffdock" / target_name / "docking_results.csv"
            if path.exists():
                all_data[target_name] = pd.read_csv(path)
        return all_data

    def predict_off_target_binding(self, mol_id, target_name, target_affinity, all_docking):
        """预测分子对其他靶点的结合亲和力"""
        off_target_affinities = {}
        
        for off_target in self.targets:
            if off_target == target_name:
                continue
            
            sim = self.domain_similarity.get((target_name, off_target), 0.2)
            
            # 基于结构相似度 + 对接结果推断脱靶亲和力
            if off_target in all_docking:
                df = all_docking[off_target]
                matched = df[df['mol_id'] == mol_id]
                if len(matched) > 0:
                    off_affinity = float(matched['estimated_binding_energy'].iloc[0])
                else:
                    # 分子未对该靶点对接，根据相似度估算
                    off_affinity = target_affinity + 2.0 + np.random.uniform(-1, 1) * (1 - sim)
            else:
                off_affinity = target_affinity + 2.5 + np.random.uniform(-1, 1) * (1 - sim)
            
            off_target_affinities[off_target] = round(off_affinity, 2)
        
        return off_target_affinities

    def calculate_selectivity_index(self, target_name, mol_id, target_affinity, off_affinities):
        """计算选择性指数 SI = ΔG_target / min(|ΔG_off-target|)"""
        target_dG = abs(target_affinity)
        off_values = [abs(v) for v in off_affinities.values()]
        
        if len(off_values) == 0:
            return 2.0, 'Selective'
        
        min_off_dG = min(off_values)
        
        if min_off_dG < 0.5:
            return 5.0, 'Highly Selective'
        
        si = target_dG / min_off_dG
        
        # 选择性分类
        if si >= 3.0:
            selectivity_class = 'Highly Selective'
        elif si >= 1.5:
            selectivity_class = 'Selective'
        elif si >= 1.0:
            selectivity_class = 'Moderately Selective'
        else:
            selectivity_class = 'Non-selective'
        
        return round(si, 2), selectivity_class

    def predict_pan_target_potential(self, target_name, mol_id, off_affinities):
        """评估多靶点潜力（有的免疫疗法需要多靶点抑制）"""
        strong_binders = sum(1 for v in off_affinities.values() if v <= -6.5)
        moderate_binders = sum(1 for v in off_affinities.values() if v <= -5.0)
        
        if strong_binders >= 2:
            return 1.0, 'Strong Pan-inhibitor'
        elif moderate_binders >= 2:
            return 0.6, 'Moderate Pan-inhibitor'
        elif strong_binders >= 1:
            return 0.3, 'Low Pan-inhibitor'
        else:
            return 0.0, 'Target-specific'

    def run_for_target(self, target_name, df_competitive, all_docking):
        """对单个靶点运行选择性分析"""
        if df_competitive is None or len(df_competitive) == 0:
            print(f"  {target_name}: 无竞争性结合数据，跳过")
            return []
        
        print(f"\n选择性分析 — {target_name}")
        
        results = []
        for _, row in df_competitive.iterrows():
            mol_id = row.get('mol_id', 'unknown')
            drug_dG = row.get('drug_dG', -7.0)
            
            off_affinities = self.predict_off_target_binding(
                mol_id, target_name, drug_dG, all_docking
            )
            
            si, si_class = self.calculate_selectivity_index(
                target_name, mol_id, drug_dG, off_affinities
            )
            
            pan_score, pan_class = self.predict_pan_target_potential(
                target_name, mol_id, off_affinities
            )
            
            min_off_dG = min(off_affinities.values()) if off_affinities else 0
            
            results.append({
                'mol_id': mol_id,
                'target': target_name,
                'selectivity_index': si,
                'selectivity_class': si_class,
                'pan_target_score': pan_score,
                'pan_target_class': pan_class,
                'min_off_target_dG': min_off_dG,
                'off_target_affinities': json.dumps(off_affinities),
            })
        
        df_results = pd.DataFrame(results)
        
        highly_sel = (df_results['selectivity_class'] == 'Highly Selective').sum()
        print(f"  高选择性分子: {highly_sel}/{len(df_results)}")
        print(f"  平均选择性指数: {df_results['selectivity_index'].mean():.2f}")
        
        output_dir = self.base_dir / "results" / "alphafold3" / target_name
        df_results.to_csv(output_dir / "selectivity_analysis.csv", index=False)
        
        return results

    def run(self):
        print("\n" + "█" * 60)
        print("█" + "选择性预测分析模块".center(58) + "█")
        print("█" * 60)
        
        all_docking = self.load_all_docking()
        
        all_results = {}
        for target_name in self.targets:
            # 加载竞争性结合结果
            cb_path = self.base_dir / "results" / "alphafold3" / target_name / "competitive_binding.csv"
            df_cb = pd.read_csv(cb_path) if cb_path.exists() else None
            
            results = self.run_for_target(target_name, df_cb, all_docking)
            if results:
                all_results[target_name] = results
        
        # 可视化数据
        viz_data = {}
        for target_name, results in all_results.items():
            df = pd.DataFrame(results)
            viz_data[target_name] = {
                'selectivity_indices': df['selectivity_index'].tolist(),
                'selectivity_distribution': df['selectivity_class'].value_counts().to_dict(),
                'mean_si': round(df['selectivity_index'].mean(), 2),
                'highly_selective': int((df['selectivity_class'] == 'Highly Selective').sum()),
            }
        
        viz_path = self.base_dir / "results" / "alphafold3" / "selectivity_viz.json"
        with open(viz_path, 'w', encoding='utf-8') as f:
            json.dump(viz_data, f, indent=2, ensure_ascii=False)
        
        print(f"\n可视化数据已保存: {viz_path}")
        return all_results


if __name__ == "__main__":
    analyzer = SelectivityAnalyzer()
    analyzer.run()