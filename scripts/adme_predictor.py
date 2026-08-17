#!/usr/bin/env python3
"""
adme_predictor.py - ADME/Tox 类药性预测（真实规则版）
使用真实的药物化学经验规则（Lipinski/Veber/Egan/Muegge/GSK）+ RDKit 真实性质，
对候选分子的类药性做真实评估。不再使用随机分布。
"""

import os
import sys
import json
import numpy as np
import pandas as pd
import yaml
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.real_chemistry import parse_molecule, compute_properties, check_pains_brenk

import warnings
warnings.filterwarnings('ignore')


class ADMEPredictor:
    """真实类药性预测器（药物化学经验规则）"""

    def __init__(self, config_path=None):
        if config_path is None:
            config_path = Path(__file__).parent.parent / "config" / "config.yaml"
        with open(config_path, 'r', encoding='utf-8') as f:
            self.config = yaml.safe_load(f)
        self.targets = list(self.config['targets'].keys())
        self.base_dir = Path(__file__).parent.parent

    def _evaluate_rules(self, props):
        """用真实药物化学规则评估类药性"""
        mw = props['molecular_weight']
        logp = props['logP']
        tpsa = props['TPSA']
        hbd = props['HBD']
        hba = props['HBA']
        rot = props['rotatable_bonds']
        rings = props['num_rings']

        # 1. Lipinski 五规则 (Ro5)
        lipinski = {
            'mw': mw <= 500,
            'logp': logp <= 5,
            'hbd': hbd <= 5,
            'hba': hba <= 10,
        }
        n_lipinski_pass = sum(lipinski.values())
        lipinski_violations = 4 - n_lipinski_pass

        # 2. Veber 规则
        veber_pass = (rot <= 10) and (tpsa <= 140)

        # 3. Egan 规则
        egan_pass = (logp <= 5.88) and (tpsa <= 131.6)

        # 4. Muegge 规则（简化）
        muegge_conditions = [
            mw >= 200 and mw <= 600,
            logp >= -2 and logp <= 5,
            tpsa <= 150,
            rings <= 7,
        ]
        muegge_pass = all(muegge_conditions)

        # 5. GSK 4/400 规则
        gsk_pass = (mw <= 400) and (logp <= 4)

        return {
            'lipinski_violations': lipinski_violations,
            'lipinski_pass': n_lipinski_pass >= 4,
            'veber_pass': veber_pass,
            'egan_pass': egan_pass,
            'muegge_pass': muegge_pass,
            'gsk_4_400_pass': gsk_pass,
        }

    def _estimate_admet_scores(self, props):
        """基于真实性质的 ADMET 简化评分（启发式但输入真实）"""
        mw = props['molecular_weight']
        logp = props['logP']
        tpsa = props['TPSA']
        rot = props['rotatable_bonds']

        # 吸收：TPSA < 140 且 rot < 10 有利于口服吸收
        absorption = 1.0
        if tpsa > 140:
            absorption -= (tpsa - 140) / 100
        if rot > 10:
            absorption -= (rot - 10) * 0.05
        absorption = float(np.clip(absorption, 0, 1))

        # 分布：MW 适中、logP 适中有利于分布
        distribution = 1.0 - abs(mw - 350) / 400 - abs(logp - 3) / 8
        distribution = float(np.clip(distribution, 0, 1))

        # 代谢：避免过高的 logP（预示快速代谢）
        metabolism = 1.0 - max(0, logp - 4) / 6
        metabolism = float(np.clip(metabolism, 0, 1))

        # 排泄：MW < 500 利于肾排泄
        excretion = float(np.clip(1.0 - mw / 800, 0.1, 1.0))

        # 毒性：无 PAINS（外部传入）
        return {
            'absorption_score': round(absorption, 3),
            'distribution_score': round(distribution, 3),
            'metabolism_score': round(metabolism, 3),
            'excretion_score': round(excretion, 3),
        }

    def predict_one(self, target_name, row):
        """预测单个分子的类药性"""
        smiles = row.get('smiles', '')
        mol = parse_molecule(smiles)
        if mol is None:
            return None

        props = compute_properties(smiles)
        rules = self._evaluate_rules(props)
        scores = self._estimate_admet_scores(props)

        # PAINS/Brenk（真实子结构过滤）
        is_pains, _ = check_pains_brenk(mol)
        toxicity_score = 0.1 if is_pains else 0.9  # PAINS命中 → 高毒性风险

        # 综合类药性评分（真实规则加权）
        admet_composite = (
            0.25 * scores['absorption_score'] +
            0.15 * scores['distribution_score'] +
            0.20 * scores['metabolism_score'] +
            0.10 * scores['excretion_score'] +
            0.30 * toxicity_score
        )
        admet_composite = float(np.clip(admet_composite, 0, 1))

        # 分级
        if admet_composite >= 0.75:
            admet_class = 'Excellent'
        elif admet_composite >= 0.55:
            admet_class = 'Good'
        elif admet_composite >= 0.35:
            admet_class = 'Moderate'
        else:
            admet_class = 'Poor'

        return {
            'mol_id': row.get('mol_id', 'unknown'),
            'target': target_name,
            'smiles': smiles,
            'admet_score': round(admet_composite, 4),
            'admet_class': admet_class,
            'absorption_score': scores['absorption_score'],
            'distribution_score': scores['distribution_score'],
            'metabolism_score': scores['metabolism_score'],
            'excretion_score': scores['excretion_score'],
            'toxicity_score': round(toxicity_score, 3),
            'lipinski_violations': rules['lipinski_violations'],
            'veber_pass': rules['veber_pass'],
            'egan_pass': rules['egan_pass'],
            'muegge_pass': rules['muegge_pass'],
            'gsk_4_400_pass': rules['gsk_4_400_pass'],
            'pains_flag': is_pains,
            'molecular_weight': props['molecular_weight'],
            'logP': props['logP'],
            'TPSA': props['TPSA'],
        }

    def run_for_target(self, target_name, df_molecules):
        if df_molecules is None or len(df_molecules) == 0:
            print(f"  {target_name}: 无分子数据，跳过")
            return []

        print(f"\n真实类药性评估 — {target_name} ({len(df_molecules)}个分子)")
        results = []
        for _, row in df_molecules.iterrows():
            r = self.predict_one(target_name, row)
            if r:
                results.append(r)

        df_results = pd.DataFrame(results)
        if len(df_results) == 0:
            return []

        print(f"  Excellent: {(df_results['admet_class']=='Excellent').sum()} | "
              f"Good: {(df_results['admet_class']=='Good').sum()} | "
              f"PAINS: {df_results['pains_flag'].sum()}")
        print(f"  平均ADMET评分: {df_results['admet_score'].mean():.3f}")

        out_dir = self.base_dir / "results" / "alphafold3" / target_name
        out_dir.mkdir(parents=True, exist_ok=True)
        df_results.to_csv(out_dir / "admet_prediction.csv", index=False)

        return results

    def run(self):
        print("\n" + "█" * 60)
        print("█" + "ADME/Tox 类药性预测 (真实规则)".center(58) + "█")
        print("█" * 60)

        all_results = {}
        # 从标准化化合物库读取
        lib_path = self.base_dir / "data" / "library" / "compounds_standardized.csv"
        df_lib = pd.read_csv(lib_path) if lib_path.exists() else None

        if df_lib is None:
            print("未找到标准化化合物库，跳过")
            return all_results

        for target_name in self.targets:
            df_t = df_lib[df_lib['target'] == target_name] if 'target' in df_lib.columns else df_lib
            results = self.run_for_target(target_name, df_t)
            if results:
                all_results[target_name] = results

        # 可视化
        viz_data = {}
        for t, results in all_results.items():
            df = pd.DataFrame(results)
            viz_data[t] = {
                'admet_scores': df['admet_score'].tolist(),
                'admet_distribution': df['admet_class'].value_counts().to_dict(),
                'mean_admet': round(df['admet_score'].mean(), 3),
                'pains_flags': int(df['pains_flag'].sum()),
            }
        with open(self.base_dir / "results" / "alphafold3" / "admet_viz.json", 'w') as f:
            json.dump(viz_data, f, indent=2, ensure_ascii=False)

        print("\n真实类药性预测完成！")
        return all_results


if __name__ == "__main__":
    predictor = ADMEPredictor()
    predictor.run()