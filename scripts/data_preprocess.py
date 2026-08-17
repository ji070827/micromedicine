#!/usr/bin/env python3
"""
data_preprocess.py - 数据预处理脚本（真实化学计算版）
统一完成靶点信息整理与小分子库标准化
小分子理化性质使用 RDKit 真实计算（real_chemistry 模块）。
"""

import os
import sys
import json
import hashlib
import numpy as np
import pandas as pd
import yaml
from pathlib import Path
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.real_chemistry import compute_properties, parse_molecule

import warnings
warnings.filterwarnings('ignore')


class DataPreprocessor:
    """数据预处理器：靶点信息 + 小分子库真实理化性质计算"""

    def __init__(self, config_path=None):
        if config_path is None:
            config_path = Path(__file__).parent.parent / "config" / "config.yaml"

        with open(config_path, 'r', encoding='utf-8') as f:
            self.config = yaml.safe_load(f)

        self.targets = self.config['targets']
        self.screening_params = self.config['screening']
        self.base_dir = Path(__file__).parent.parent

    def preprocess_proteins(self):
        """
        整理靶点蛋白信息（从 config 读取真实的 PDB ID、结合位点残基、功能域）。
        不做随机坐标模拟——真实结构由 generate_3d_complex.py 用拓扑模板生成，
        这里只负责汇总生物学注释信息。
        """
        print("=" * 60)
        print("靶点蛋白信息整理...")
        print("=" * 60)

        protein_results = {}
        for target_name in self.targets:
            t = self.targets[target_name]
            info = {
                'target': target_name,
                'pdb_id': t['pdb_id'],
                'chain': t['chain'],
                'functional_domain': t['functional_domain'],
                'description': t['description'],
                'binding_site_residues': t['binding_site_residues'],
                'n_binding_residues': len(t['binding_site_residues']),
            }
            protein_results[target_name] = info

            # 保存
            out = self.base_dir / "data" / "targets" / f"{target_name}_preprocessed.json"
            with open(out, 'w', encoding='utf-8') as f:
                json.dump(info, f, indent=2, ensure_ascii=False)
            print(f"  {target_name}: PDB={t['pdb_id']}, "
                  f"{info['n_binding_residues']}个结合位点残基")

        print("\n靶点信息整理完成！")
        return protein_results

    def load_library(self):
        """
        加载小分子化合物库（通用：自动识别 SMILES 列 + 多格式支持）。
        优先读取 config 里配置的 active_library，失配时回退到默认库。
        """
        from scripts.library_loader import load_library_file

        # 1. 优先使用 config 里配置的激活库
        active_lib = self.config.get('data', {}).get('active_library', 'pubchem_all_targets.csv')
        lib_path = self.base_dir / "data" / "library" / active_lib
        if lib_path.exists():
            df = load_library_file(lib_path)
            print(f"  加载化合物库: {lib_path.name}（{len(df)} 个分子，SMILES列自动识别）")
            return df

        # 2. 回退到默认库
        default_path = self.base_dir / "data" / "library" / "pubchem_all_targets.csv"
        if default_path.exists():
            df = load_library_file(default_path)
            print(f"  加载默认库: {default_path.name}（{len(df)} 个分子）")
            return df

        # 3. 都没有则尝试生成
        print("  警告: 未找到化合物库文件")
        print("  请先运行: python scripts/generate_real_library.py")
        return None

    def preprocess_compounds(self):
        """
        用 RDKit 真实计算每个分子的理化性质、QED、SA score、PAINS 过滤。
        替代原先的 np.random 随机模拟。
        """
        print("\n" + "=" * 60)
        print("小分子库真实理化性质计算 (RDKit)...")
        print("=" * 60)

        df = self.load_library()
        if df is None:
            return None, None

        records = []
        n_parse_fail = 0
        n_pains = 0

        for _, row in df.iterrows():
            smiles = row.get('smiles', '')
            mol = parse_molecule(smiles)

            if mol is None:
                n_parse_fail += 1
                continue

            # 真实计算理化性质
            props = compute_properties(smiles)

            # 真实 PAINS/Brenk 过滤
            from scripts.real_chemistry import check_pains_brenk
            is_pains, pains_desc = check_pains_brenk(mol)
            if is_pains:
                n_pains += 1

            rec = {
                'mol_id': row.get('mol_id', f"CMPD_{len(records)+1:06d}"),
                'source': row.get('source', 'library'),
                'target': row.get('target', 'all'),
                'smiles': smiles,
                'iupac_name': row.get('iupac_name', ''),
                'hash_id': hashlib.md5(smiles.encode()).hexdigest()[:8],
                **props,
                'pains_brenk_flag': is_pains,
                'pains_description': '; '.join(pains_desc) if pains_desc else '',
            }
            records.append(rec)

        df_clean = pd.DataFrame(records)
        print(f"  解析成功: {len(df_clean)} 个")
        print(f"  解析失败(SMILES非法): {n_parse_fail} 个")
        print(f"  PAINS/Brenk 警示: {n_pains} 个")

        if len(df_clean) == 0:
            print("  错误: 没有有效分子")
            return None, None

        # 按 SMILES 去重
        n_before = len(df_clean)
        df_clean = df_clean.drop_duplicates(subset=['hash_id'])
        n_dup = n_before - len(df_clean)
        if n_dup > 0:
            print(f"  去除重复: {n_dup} 个")

        # 统计
        physchem_stats = {
            'data_source': 'RDKit真实计算',
            'total_compounds': len(df_clean),
            'parse_failed': n_parse_fail,
            'duplicates_removed': n_dup,
            'pains_flags': n_pains,
            'lipinski_pass': int(df_clean['lipinski_pass'].sum()),
            'lipinski_fail': int((~df_clean['lipinski_pass']).sum()),
            'pass_rate': round(df_clean['lipinski_pass'].mean() * 100, 1),
            'mean_mw': round(df_clean['molecular_weight'].mean(), 2),
            'mean_logp': round(df_clean['logP'].mean(), 2),
            'mean_tpsa': round(df_clean['TPSA'].mean(), 2),
            'mean_qed': round(df_clean['QED'].mean(), 3),
            'mean_sa': round(df_clean['sa_score'].mean(), 3),
        }

        # 保存标准化库
        output_path = self.base_dir / "data" / "library" / "compounds_standardized.csv"
        df_clean.to_csv(output_path, index=False)

        # 保存统计
        stats_path = self.base_dir / "data" / "library" / "physchem_statistics.json"
        with open(stats_path, 'w', encoding='utf-8') as f:
            json.dump(physchem_stats, f, indent=2, ensure_ascii=False)

        print(f"\n真实理化性质统计:")
        print(f"  总化合物: {physchem_stats['total_compounds']}")
        print(f"  Lipinski通过率: {physchem_stats['pass_rate']}%")
        print(f"  平均MW: {physchem_stats['mean_mw']}")
        print(f"  平均LogP: {physchem_stats['mean_logp']}")
        print(f"  平均TPSA: {physchem_stats['mean_tpsa']}")
        print(f"  平均QED: {physchem_stats['mean_qed']}")
        print(f"  平均SA score: {physchem_stats['mean_sa']}")

        return df_clean, physchem_stats

    def generate_visualization_data(self, df, physchem_stats):
        """生成可视化数据（真实分布）"""
        if df is None:
            return {}

        viz_data = {
            'physchem_stats': physchem_stats,
            'mw_distribution': df['molecular_weight'].tolist(),
            'logp_distribution': df['logP'].tolist(),
            'tpsa_distribution': df['TPSA'].tolist(),
            'qed_distribution': df['QED'].tolist(),
            'sa_distribution': df['sa_score'].tolist(),
            'lipinski_counts': {
                'pass': int(df['lipinski_pass'].sum()),
                'fail': int((~df['lipinski_pass']).sum()),
            },
            'property_ranges': {
                'mw': {'min': float(df['molecular_weight'].min()), 'max': float(df['molecular_weight'].max())},
                'logp': {'min': float(df['logP'].min()), 'max': float(df['logP'].max())},
                'tpsa': {'min': float(df['TPSA'].min()), 'max': float(df['TPSA'].max())},
                'hbd': {'min': int(df['HBD'].min()), 'max': int(df['HBD'].max())},
                'hba': {'min': int(df['HBA'].min()), 'max': int(df['HBA'].max())},
                'qed': {'min': float(df['QED'].min()), 'max': float(df['QED'].max())},
            }
        }

        viz_path = self.base_dir / "data" / "library" / "visualization_data.json"
        with open(viz_path, 'w', encoding='utf-8') as f:
            json.dump(viz_data, f, indent=2, ensure_ascii=False)

        return viz_data

    def run(self):
        print("\n" + "█" * 60)
        print("█" + "数据预处理模块 (RDKit真实计算)".center(58) + "█")
        print("█" * 60)

        protein_results = self.preprocess_proteins()
        df_compounds, physchem_stats = self.preprocess_compounds()
        viz_data = self.generate_visualization_data(df_compounds, physchem_stats)

        print("\n" + "█" * 60)
        print("█" + "数据预处理全部完成！".center(58) + "█")
        print("█" * 60)

        return {
            'protein_results': protein_results,
            'compounds_df': df_compounds,
            'physchem_stats': physchem_stats,
            'viz_data': viz_data
        }


if __name__ == "__main__":
    preprocessor = DataPreprocessor()
    preprocessor.run()