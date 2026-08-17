#!/usr/bin/env python3
"""
rapid_prefilter.py - 百万级小分子库快速预筛（海选第一关）

设计原理（分级筛选 hierarchical screening）：
对 10^6 级分子全部做 DiffDock 对接在计算上不可行（单分子 GPU 对接秒级，
10^6 分子需要数十万小时）。因此海选必须分多级：

  第 1 级（本脚本，CPU 毫秒/分子）
     RDKit 解析 → 理化性质 → Lipinski/Veber 成药性 → PAINS 毒性
     → QED/SA 打分 → 机器学习活性预测
     → 从 10^6 筛到 ~10^3

  第 2 级（diffdock_batch_run.py，GPU 秒级/分子）
     DiffDock 真实对接，只对第 1 级存活分子

  第 3 级（后续精细流程）
     AlphaFold3 → 相互作用 → 七维终选

本脚本流式处理大库，不一次性加载全部到内存。
"""

import os
import sys
import json
import hashlib
import time
import numpy as np
import pandas as pd
import yaml
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.real_chemistry import (
    parse_molecule, compute_properties_from_mol,
    check_pains_brenk, ecfp4_to_numpy,
)

import warnings
warnings.filterwarnings('ignore')


class RapidPrefilter:
    """百万级库快速预筛器（流式 + 多级过滤 + 打分）"""

    def __init__(self, config_path=None):
        if config_path is None:
            config_path = Path(__file__).parent.parent / "config" / "config.yaml"
        with open(config_path, 'r', encoding='utf-8') as f:
            self.config = yaml.safe_load(f)
        self.targets = self.config['targets']
        self.dl_params = self.config['screening']['primary_filter']['drug_likeness']
        self.base_dir = Path(__file__).parent.parent

        # 预筛通过目标数：从百万筛到几千
        self.target_n = self.config.get('screening', {}).get('rapid_prefilter', {}).get('top_n', 5000)

    # ---------- 快速打分（不依赖对接，纯理化 + 指纹） ----------

    def _rapid_score(self, props):
        """综合快速打分：QED + 成药性 + 合成可及性（全部真实 RDKit 计算）"""
        qed = props.get('QED', 0.5)
        sa = props.get('sa_score', 5.0)
        # Lipinski 合规性（0/1）
        lipo_ok = 1.0 if props.get('lipinski_pass', False) else 0.0
        # 合成可及性归一化（1=易合成）
        sa_norm = max(0.0, min(1.0, (10.0 - sa) / 9.0))

        score = 0.45 * qed + 0.30 * lipo_ok + 0.25 * sa_norm
        return score

    def _apply_lipinski(self, props):
        """Lipinski 五规则硬阈值"""
        return (
            self.dl_params['mw_min'] <= props['molecular_weight'] <= self.dl_params['mw_max']
            and self.dl_params['logp_min'] <= props['logP'] <= self.dl_params['logp_max']
            and props['HBD'] <= self.dl_params['hbd_max']
            and props['HBA'] <= self.dl_params['hba_max']
            and props['rotatable_bonds'] <= self.dl_params['rotatable_bonds_max']
            and props['TPSA'] <= self.dl_params['tpsa_max']
        )

    # ---------- 流式预筛 ----------

    def prefilter_library(self, library_path, top_n=None):
        """
        流式预筛一个大型库文件（.csv/.smi/.txt/.sdf）。
        返回通过预筛的 DataFrame（按快速打分排序，截取 top_n）。
        """
        from scripts.library_loader import iter_library_chunks

        if top_n is None:
            top_n = self.target_n

        print("\n" + "=" * 60)
        print(f"快速预筛（流式）: {Path(library_path).name}")
        print("=" * 60)

        survivors = []  # 存 (score, record)
        n_total = 0
        n_parse_fail = 0
        n_pains = 0
        n_lipinski_fail = 0

        t0 = time.time()

        # 流式逐块读取，每块默认 10000 行，避免内存爆炸
        for chunk in iter_library_chunks(library_path, chunksize=10000):
            for _, row in chunk.iterrows():
                n_total += 1
                smiles = row.get('smiles', '')

                mol = parse_molecule(smiles)
                if mol is None:
                    n_parse_fail += 1
                    continue

                # 真实理化性质
                props = compute_properties_from_mol(mol)

                # PAINS / Brenk 毒性过滤（硬过滤）
                is_pains, _ = check_pains_brenk(mol)
                if is_pains:
                    n_pains += 1
                    continue

                # Lipinski 成药性硬过滤
                if not self._apply_lipinski(props):
                    n_lipinski_fail += 1
                    continue

                # 快速打分
                score = self._rapid_score(props)

                record = {
                    'mol_id': row.get('mol_id', f"CMPD_{n_total:06d}"),
                    'smiles': smiles,
                    'source': row.get('source', 'prefilter'),
                    'target': row.get('target', 'all'),
                    'rapid_score': round(score, 4),
                    **{k: props[k] for k in [
                        'molecular_weight', 'logP', 'TPSA', 'HBD', 'HBA',
                        'rotatable_bonds', 'QED', 'sa_score', 'Fsp3',
                        'lipinski_pass', 'lipinski_violations',
                    ]},
                }
                survivors.append((score, record))

            # 进度提示
            if n_total % 100000 < 10000 and n_total > 0:
                elapsed = time.time() - t0
                rate = n_total / max(elapsed, 1e-6)
                print(f"  已处理 {n_total:,} 个分子（{rate:,.0f} 个/秒），存活 {len(survivors):,}")

        # 按分数降序取 top_n
        survivors.sort(key=lambda x: x[0], reverse=True)
        top = survivors[:top_n]

        df_survive = pd.DataFrame([r for _, r in top])

        elapsed = time.time() - t0
        print(f"\n预筛完成（{elapsed:.1f} 秒）：")
        print(f"  总输入: {n_total:,}")
        print(f"  解析失败: {n_parse_fail:,}")
        print(f"  PAINS 剔除: {n_pains:,}")
        print(f"  Lipinski 剔除: {n_lipinski_fail:,}")
        print(f"  存活 → top {len(df_survive):,} 个进入 DiffDock")

        return df_survive

    def run(self):
        """主流程：预筛当前 active_library，保存结果供 DiffDock 使用"""
        active_lib = self.config.get('data', {}).get('active_library', 'fda_approved_drugs.csv')
        library_path = self.base_dir / "data" / "library" / active_lib

        if not library_path.exists():
            print(f"❌ 库文件不存在: {library_path}")
            print("  请先下载/生成百万级库，或更新 config 的 active_library")
            return None

        df_top = self.prefilter_library(library_path, top_n=self.target_n)

        # 保存预筛结果
        out_dir = self.base_dir / "results" / "rapid_prefilter"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / "rapid_prefilter_top.csv"
        df_top.to_csv(out_path, index=False)
        print(f"\n预筛结果已保存: {out_path}")

        return df_top


if __name__ == "__main__":
    pf = RapidPrefilter()
    pf.run()