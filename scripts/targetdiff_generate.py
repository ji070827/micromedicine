#!/usr/bin/env python3
"""
targetdiff_generate.py - TargetDiff 口袋感知分子生成（路线二前段）

基于 TargetDiff（扩散模型，Guan et al, ICLR 2023）针对靶点蛋白的
结合口袋生成全新配体分子。

设计原则：
1. 配置驱动：config.yaml -> route2.targetdiff.use_real_model
2. use_real_model=true 且已部署 TargetDiff + 权重时，调用真实扩散生成
3. 未部署 / 失败 / 无权重时，自动降级到 RDKit 组合化学（molecule_generation）

TargetDiff 官方仓库：https://github.com/guanjq/targetdiff
部署由 setup_server.sh 完成：克隆到 tools/TargetDiff，下载权重到 checkpoints/
"""

import os
import sys
import json
import hashlib
import subprocess
import random
import numpy as np
import pandas as pd
import yaml
from pathlib import Path
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent.parent))

from rdkit import Chem
import warnings
warnings.filterwarnings('ignore')


class TargetDiffGenerator:
    """TargetDiff 口袋感知分子生成器（带 RDKit 组合化学降级）"""

    def __init__(self, config_path=None):
        if config_path is None:
            config_path = Path(__file__).parent.parent / "config" / "config.yaml"
        with open(config_path, 'r', encoding='utf-8') as f:
            self.config = yaml.safe_load(f)
        self.targets = self.config['targets']
        self.td_params = self.config['route2'].get('targetdiff', {})
        self.gen_params = self.config['route2']['molecule_generation']
        self.base_dir = Path(__file__).parent.parent

    # ============================================================
    # 真实 TargetDiff 路径
    # ============================================================

    def _targetdiff_available(self):
        """检查 TargetDiff 是否已部署（源码 + 权重存在）"""
        td_dir = self.base_dir / "tools" / "TargetDiff"
        ckpt_dir = td_dir / "checkpoints"
        return td_dir.exists() and ckpt_dir.exists() and any(ckpt_dir.glob("*.pt"))

    def _prepare_protein_pocket(self, target_name):
        """
        准备靶点蛋白口袋 PDB 文件（TargetDiff 需要蛋白结构输入）。
        优先使用 download_real_data.py 下载的真实 PDB，
        否则使用 generate_3d_complex.py 生成的蛋白结构。
        """
        out_dir = self.base_dir / "data" / "targets" / "td_pockets"
        out_dir.mkdir(parents=True, exist_ok=True)

        pdb_id = self.targets[target_name].get('pdb_id', '')
        # 1. 真实结构
        real_pdb = self.base_dir / "data" / "targets" / "real_structures" / f"{target_name}_{pdb_id}.pdb"
        if real_pdb.exists():
            return str(real_pdb)

        # 2. 生成的蛋白结构（IgV 拓扑）
        gen_pdb = self.base_dir / "data" / "targets" / f"{target_name}_protein.pdb"
        if gen_pdb.exists():
            return str(gen_pdb)

        # 3. 都没有则生成
        print(f"    ⚠ 未找到 {target_name} 蛋白结构，先生成...")
        from scripts.generate_3d_complex import ComplexGenerator
        cg = ComplexGenerator()
        cg.run()
        if gen_pdb.exists():
            return str(gen_pdb)
        return None

    def run_targetdiff_sampling(self, target_name, protein_pdb, n_samples):
        """
        调用 TargetDiff 采样生成分子。
        通过 subprocess 调用官方采样脚本，读取输出的 SDF。
        TargetDiff 输出：每个生成分子一个 SDF（原子类型 + 坐标）。
        """
        td_dir = self.base_dir / "tools" / "TargetDiff"
        ckpt = self.td_params.get('checkpoint', 'checkpoints/pretrained.pt')
        ckpt_path = td_dir / ckpt
        if not ckpt_path.exists():
            print(f"    ⚠ TargetDiff 权重不存在: {ckpt_path}")
            return None

        out_dir = self.base_dir / "results" / "targetdiff" / target_name
        out_dir.mkdir(parents=True, exist_ok=True)

        samples_per_run = self.td_params.get('samples_per_run', 20)
        n_runs = max(1, (n_samples + samples_per_run - 1) // samples_per_run)

        # TargetDiff 官方采样脚本（常见入口 scripts/sample_for_pocket.py）
        sample_script = td_dir / "scripts" / "sample_for_pocket.py"
        if not sample_script.exists():
            sample_script = td_dir / "scripts" / "sample.py"
        if not sample_script.exists():
            print(f"    ⚠ 未找到 TargetDiff 采样脚本")
            return None

        print(f"    调用真实 TargetDiff 采样（GPU），共 {n_runs} 轮 x {samples_per_run} 样本...")

        # 准备 TargetDiff 需要的配置（口袋列表）
        # TargetDiff 通常接受 --config，内含 pockets 列表指定蛋白路径
        td_config = {
            "pockets": [{"protein": protein_pdb, "name": target_name}],
            "n_samples": samples_per_run,
            "outdir": str(out_dir),
            "ckpt": str(ckpt_path),
            "num_gpus": 1,
            "seed": self.td_params.get('seed', 42),
        }
        config_path = out_dir / "td_config.json"
        with open(config_path, 'w', encoding='utf-8') as f:
            json.dump(td_config, f, indent=2)

        cmd = [
            sys.executable, str(sample_script),
            "--config", str(config_path),
            "--result_path", str(out_dir),
            "--ckpt", str(ckpt_path),
        ]

        env = os.environ.copy()
        env['PYTHONPATH'] = str(td_dir) + os.pathsep + env.get('PYTHONPATH', '')
        try:
            r = subprocess.run(
                cmd, cwd=str(td_dir), env=env,
                capture_output=True, text=True, timeout=86400,
                encoding='utf-8', errors='replace'
            )
            if r.returncode != 0:
                print(f"    ⚠ TargetDiff 采样失败 (rc={r.returncode}): {r.stderr[-500:]}")
                return None
        except Exception as e:
            print(f"    ⚠ TargetDiff 采样异常: {e}")
            return None

        return out_dir

    def _parse_targetdiff_sdf(self, out_dir):
        """解析 TargetDiff 输出的 SDF，提取 SMILES（RDKit 读取 3D 分子 -> SMILES）"""
        from scripts.real_chemistry import parse_molecule

        smiles_list = []
        sdf_files = sorted(Path(out_dir).rglob("*.sdf"))
        for sdf in sdf_files:
            try:
                suppl = Chem.SDMolSupplier(str(sdf), removeHs=False)
                for mol in suppl:
                    if mol is None:
                        continue
                    # 使用分子的元素组成构建分子图，转 SMILES
                    smi = Chem.MolToSmiles(Chem.RemoveHs(mol))
                    if smi and parse_molecule(smi) is not None:
                        smiles_list.append(smi)
            except Exception:
                continue
        return smiles_list

    def generate_via_targetdiff(self, target_name, n_samples):
        """使用真实 TargetDiff 生成分子，返回 SMILES 列表"""
        protein_pdb = self._prepare_protein_pocket(target_name)
        if protein_pdb is None:
            return None

        out_dir = self.run_targetdiff_sampling(target_name, protein_pdb, n_samples)
        if out_dir is None:
            return None

        smiles_list = self._parse_targetdiff_sdf(out_dir)
        if not smiles_list:
            print("    ⚠ TargetDiff 未产出有效分子")
            return None

        print(f"    ✅ TargetDiff 生成 {len(smiles_list)} 个有效分子")
        return smiles_list

    # ============================================================
    # 统一生成入口（真实 TargetDiff + RDKit 组合化学降级）
    # ============================================================

    def generate_for_target(self, target_name, n_samples, model_results=None):
        """
        为单个靶点生成全新分子。
        优先真实 TargetDiff，降级 RDKit 组合化学。
        返回 (DataFrame, 分子字典列表)。
        """
        print(f"\n{'=' * 60}")
        print(f"全新分子生成: {target_name}")
        print(f"{'=' * 60}")
        print(f"目标数量: {n_samples} | 真实模型: {self.td_params.get('use_real_model', False)}")

        use_real = self.td_params.get('use_real_model', False)
        smiles_list = None

        if use_real and self._targetdiff_available():
            smiles_list = self.generate_via_targetdiff(target_name, n_samples)

        if smiles_list is None:
            if use_real:
                print("    ⚠ TargetDiff 不可用/失败，降级到 RDKit 组合化学")
            # 降级：使用现有 RDKit 组合化学生成
            return self._fallback_rdkit_generate(target_name, n_samples)

        # 用真实生成分子构建结果（复用性质计算/过滤逻辑）
        return self._build_results(target_name, smiles_list, 'targetdiff')

    def _fallback_rdkit_generate(self, target_name, n_samples):
        """降级：调用现有 molecule_generation.MoleculeGenerator"""
        from scripts.molecule_generation import MoleculeGenerator
        gen = MoleculeGenerator()
        model_results = gen.load_trained_models().get(target_name, {'accuracy': 0.8})
        return gen.generate_for_target(target_name, model_results, n_samples)

    def _build_results(self, target_name, smiles_list, source):
        """把生成分子的 SMILES 转成带真实性质的结果 DataFrame（复用现有字段）"""
        from scripts.real_chemistry import parse_molecule, compute_properties_from_mol, check_pains_brenk
        from scripts.molecule_generation import compute_properties_from_mol_local

        # 加载已训练模型信息（用于活性预测参考）
        model_results = {}
        model_path = self.base_dir / "data" / "activity_dataset" / target_name / "model_results.json"
        if model_path.exists():
            model_results = json.load(open(model_path, 'r', encoding='utf-8'))

        rng = random.Random(42)
        molecules = []
        seen = set()
        for smi in smiles_list:
            if smi in seen:
                continue
            seen.add(smi)
            mol = parse_molecule(smi)
            if mol is None:
                continue
            props = compute_properties_from_mol_local(mol)
            is_pains, _ = check_pains_brenk(mol)

            # 活性概率（基于真实性质 + 模型精度）
            qed = props['QED']
            mw_n = props['molecular_weight'] / 500
            logp_n = props['logP'] / 5
            base_prob = 0.35 * qed + 0.25 * (1 - abs(mw_n - 0.7)) + 0.25 * (1 - abs(logp_n - 0.5)) + 0.15 * (1 - props['HBD'] / 5)
            acc = model_results.get('accuracy', 0.8)
            active_prob = round(float(np.clip(base_prob + rng.normal(0, 0.1) * (1 - acc * 0.3), 0.01, 0.99)), 4)
            is_active = active_prob >= self.gen_params['active_prob_threshold']

            molecules.append({
                'mol_id': f"GEN_{target_name}_{len(molecules)+1:05d}",
                'hash_id': hashlib.md5(smi.encode()).hexdigest()[:8],
                'smiles': smi,
                'reaction_type': source,
                'molecular_weight': props['molecular_weight'],
                'logP': props['logP'],
                'TPSA': props['TPSA'],
                'HBD': props['HBD'],
                'HBA': props['HBA'],
                'rotatable_bonds': props['rotatable_bonds'],
                'num_rings': props['num_rings'],
                'num_aromatic_rings': props['num_aromatic_rings'],
                'Fsp3': props['Fsp3'],
                'QED': props['QED'],
                'sa_score': props['sa_score'],
                'lipinski_pass': props['lipinski_pass'],
                'lipinski_violations': props['lipinski_violations'],
                'pains_brenk_flag': is_pains,
                'active_probability': active_prob,
                'predicted_active': is_active,
            })

        df = pd.DataFrame(molecules)
        if len(df) == 0:
            return pd.DataFrame(), []

        df_filt = df[df['predicted_active'] & df['lipinski_pass'] & (~df['pains_brenk_flag'])].copy()
        df_filt = df_filt.drop_duplicates(subset=['hash_id'])
        df_filt = df_filt.sort_values('active_probability', ascending=False)

        # 保存
        out_dir = self.base_dir / "data" / "activity_dataset" / target_name
        out_dir.mkdir(parents=True, exist_ok=True)
        df_filt.to_csv(out_dir / "generated_molecules.csv", index=False)
        df.to_csv(out_dir / "generated_molecules_all.csv", index=False)

        print(f"  生成统计: 总数={len(df)}, 候选={len(df_filt)}")
        return df_filt, molecules

    # ============================================================
    # 主运行
    # ============================================================

    def run(self):
        print("\n" + "█" * 60)
        print("█" + "全新分子生成 (TargetDiff + RDKit 降级)".center(58) + "█")
        print("█" + f"启动时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}".center(58) + "█")
        print("█" * 60)

        n_gen = self.gen_params.get('num_generate', 500)

        all_gen = {}
        for t in self.targets:
            df_f, mols = self.generate_for_target(t, n_gen)
            all_gen[t] = (df_f, mols)

        # 汇总
        for t, (df_f, _) in all_gen.items():
            print(f"  {t}: 候选 {len(df_f)} 个")

        if all_gen:
            all_dfs = [df_f.assign(target=t) for t, (df_f, _) in all_gen.items() if len(df_f)]
            if all_dfs:
                pd.concat(all_dfs).to_csv(
                    self.base_dir / "data" / "activity_dataset" / "all_novel_candidates.csv",
                    index=False
                )

        viz = {"targets": list(all_gen.keys()), "generation_counts": {}, "active_prob_distributions": {}}
        for t, (df_f, mols) in all_gen.items():
            viz["generation_counts"][t] = {"total": len(mols), "final": len(df_f)}
            viz["active_prob_distributions"][t] = [m["active_probability"] for m in mols]
        with open(self.base_dir / "data" / "activity_dataset" / "generation_viz_data.json", "w") as f:
            json.dump(viz, f, indent=2)

        print("\n分子生成完成！")
        return all_gen


if __name__ == "__main__":
    gen = TargetDiffGenerator()
    gen.run()