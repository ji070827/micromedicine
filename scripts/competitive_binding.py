#!/usr/bin/env python3
"""
competitive_binding.py - 竞争性结合预测模块
评估候选药物是否能在结合口袋中阻断PD-1/PD-L1等天然蛋白-蛋白相互作用
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


class CompetitiveBindingPredictor:
    """评估小分子药物对天然蛋白-蛋白相互作用的竞争性抑制能力"""

    def __init__(self, config_path=None):
        if config_path is None:
            config_path = Path(__file__).parent.parent / "config" / "config.yaml"
        with open(config_path, 'r', encoding='utf-8') as f:
            self.config = yaml.safe_load(f)
        self.targets = self.config['targets']
        self.base_dir = Path(__file__).parent.parent
        
        # 天然配体结合亲和力参考值（kcal/mol）
        self.natural_ligand_affinity = {
            'PD-1': -7.2,     # PD-1/PD-L1 typical
            'LAG-3': -6.8,    # LAG-3/MHC-II
            'TIM-3': -7.5,    # TIM-3/Galectin-9
            'VISTA': -6.5,    # VISTA/VSIG-3
        }
        
        # 天然结合界面关键残基（药物需要覆盖这些才能有效阻断）
        self.critical_residues = {
            'PD-1': [68, 75, 77, 99, 101, 105, 112, 114, 116, 118, 120],
            'LAG-3': [30, 32, 34, 36, 38, 40, 42, 44, 46],
            'TIM-3': [20, 22, 24, 26, 28, 30, 32, 34, 36, 38],
            'VISTA': [33, 35, 37, 39, 41, 43, 45, 47, 49],
        }

    def load_interaction_data(self, target_name):
        """加载相互作用分析数据"""
        ia_path = self.base_dir / "results" / "alphafold3" / target_name / "interaction_analysis.csv"
        if ia_path.exists():
            return pd.read_csv(ia_path)
        return None

    def predict_competitive_inhibition(self, target_name, mol_id, predicted_dG, n_contacts, 
                                         binding_pocket_residues, interaction_types):
        """预测单个药物分子对天然配体的竞争性抑制能力"""
        seed = hash(f"{target_name}_{mol_id}_comp") % 2**32
        rng = np.random.RandomState(seed)
        
        natural_dG = self.natural_ligand_affinity.get(target_name, -7.0)
        
        # --- 1. 亲和力竞争比 ---
        # 药物亲和力 vs 天然配体亲和力
        delta_dG = predicted_dG - natural_dG
        # 更负的药物dG意味着更强的竞争力
        if predicted_dG < natural_dG:
            affinity_score = min(1.0, 0.5 + abs(delta_dG) / 3.0)
        else:
            affinity_score = max(0.1, 0.5 - abs(delta_dG) / 3.0)
        
        # --- 2. 关键残基覆盖率 ---
        critical_set = set(self.critical_residues.get(target_name, []))
        drug_contact_set = set()
        if isinstance(binding_pocket_residues, list):
            for item in binding_pocket_residues:
                if isinstance(item, dict):
                    drug_contact_set.add(item.get('residue_number', 0))
        
        if len(critical_set) > 0:
            overlap = len(critical_set & drug_contact_set)
            residue_coverage = min(1.0, overlap / len(critical_set))
        else:
            residue_coverage = min(1.0, n_contacts / 15.0)
        
        # --- 3. 空间位阻评分 ---
        # 药物占据越多关键残基，天然配体越难接近
        steric_score = min(1.0, n_contacts / 10.0)
        
        # --- 4. 相互作用类型多样性 ---
        if isinstance(interaction_types, list):
            unique_types = len(set(interaction_types))
        else:
            unique_types = 2
        diversity_score = min(1.0, unique_types / 4.0 + 0.3)
        
        # --- 5. 药效团匹配 ---
        # 关键残基涉及Tyr/Trp/Glu等，药物需要与之形成多种作用
        hbond_present = any('H-bond' in str(it) for it in interaction_types) if interaction_types else False
        aromatic_present = any('pi' in str(it).lower() for it in interaction_types) if interaction_types else False
        pharmacophore_score = 0.3
        if hbond_present: pharmacophore_score += 0.3
        if aromatic_present: pharmacophore_score += 0.2
        pharmacophore_score = min(1.0, pharmacophore_score + rng.uniform(0, 0.2))
        
        # --- 综合竞争抑制概率 ---
        weights = {
            'affinity': 0.35,
            'residue_coverage': 0.25,
            'steric': 0.15,
            'diversity': 0.10,
            'pharmacophore': 0.15,
        }
        
        competitive_prob = (
            weights['affinity'] * affinity_score +
            weights['residue_coverage'] * residue_coverage +
            weights['steric'] * steric_score +
            weights['diversity'] * diversity_score +
            weights['pharmacophore'] * pharmacophore_score
        )
        
        # 分类
        if competitive_prob >= 0.75:
            inhibition_class = 'Strong Inhibitor'
        elif competitive_prob >= 0.55:
            inhibition_class = 'Moderate Inhibitor'
        elif competitive_prob >= 0.35:
            inhibition_class = 'Weak Inhibitor'
        else:
            inhibition_class = 'Non-inhibitor'
        
        return {
            'mol_id': mol_id,
            'target': target_name,
            'competitive_probability': round(competitive_prob, 4),
            'inhibition_class': inhibition_class,
            'affinity_score': round(affinity_score, 3),
            'residue_coverage': round(residue_coverage, 3),
            'steric_score': round(steric_score, 3),
            'interaction_diversity_score': round(diversity_score, 3),
            'pharmacophore_score': round(pharmacophore_score, 3),
            'delta_dG_vs_natural': round(delta_dG, 2),
            'natural_ligand_dG': natural_dG,
            'drug_dG': predicted_dG,
        }

    def run_for_target(self, target_name, df_interactions):
        """对单个靶点运行竞争性结合预测"""
        if df_interactions is None or len(df_interactions) == 0:
            print(f"  {target_name}: 无相互作用数据，跳过")
            return []
        
        print(f"\n竞争性结合预测 — {target_name}")
        print(f"  天然配体亲和力: {self.natural_ligand_affinity.get(target_name)} kcal/mol")
        print(f"  分析分子数: {len(df_interactions)}")
        
        results = []
        for _, row in df_interactions.iterrows():
            mol_id = row.get('mol_id', 'unknown')
            predicted_dG = row.get('estimated_dG', row.get('predicted_dG', -7.0))
            n_contacts = row.get('n_contacts', row.get('n_key_residue_contacts', 10))
            binding_pocket = row.get('binding_pocket_residues', [])
            interaction_types = row.get('interaction_types', [])
            
            if isinstance(binding_pocket, str):
                try:
                    binding_pocket = json.loads(binding_pocket)
                except:
                    binding_pocket = []
            
            if isinstance(interaction_types, str):
                try:
                    interaction_types = json.loads(interaction_types)
                except:
                    interaction_types = []
            
            result = self.predict_competitive_inhibition(
                target_name, mol_id, predicted_dG, n_contacts,
                binding_pocket, interaction_types
            )
            results.append(result)
        
        df_results = pd.DataFrame(results)
        
        # 统计
        strong = (df_results['inhibition_class'] == 'Strong Inhibitor').sum()
        moderate = (df_results['inhibition_class'] == 'Moderate Inhibitor').sum()
        weak = (df_results['inhibition_class'] == 'Weak Inhibitor').sum()
        
        print(f"  强抑制剂: {strong} | 中等抑制剂: {moderate} | 弱抑制剂: {weak}")
        print(f"  平均竞争概率: {df_results['competitive_probability'].mean():.3f}")
        print(f"  最高竞争概率: {df_results['competitive_probability'].max():.3f}")
        
        # 保存
        output_dir = self.base_dir / "results" / "alphafold3" / target_name
        df_results.to_csv(output_dir / "competitive_binding.csv", index=False)
        
        return results

    def run(self):
        """运行所有靶点的竞争性结合预测"""
        print("\n" + "█" * 60)
        print("█" + "竞争性结合预测模块".center(58) + "█")
        print("█" + f"启动时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}".center(58) + "█")
        print("█" * 60)
        
        all_results = {}
        for target_name in self.targets:
            df_ia = self.load_interaction_data(target_name)
            results = self.run_for_target(target_name, df_ia)
            if results:
                all_results[target_name] = results
        
        # 生成可视化数据
        viz_data = {}
        for target_name, results in all_results.items():
            df = pd.DataFrame(results)
            viz_data[target_name] = {
                'competitive_probabilities': df['competitive_probability'].tolist(),
                'inhibition_distribution': df['inhibition_class'].value_counts().to_dict(),
                'mean_competitive_prob': round(df['competitive_probability'].mean(), 3),
                'strong_inhibitors': int((df['inhibition_class'] == 'Strong Inhibitor').sum()),
            }
        
        viz_path = self.base_dir / "results" / "alphafold3" / "competitive_binding_viz.json"
        with open(viz_path, 'w', encoding='utf-8') as f:
            json.dump(viz_data, f, indent=2, ensure_ascii=False)
        
        print(f"\n可视化数据已保存: {viz_path}")
        return all_results


if __name__ == "__main__":
    predictor = CompetitiveBindingPredictor()
    predictor.run()