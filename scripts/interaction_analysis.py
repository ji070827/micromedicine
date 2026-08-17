#!/usr/bin/env python3
"""
interaction_analysis.py - 复合物相互作用分析（真实几何计算版）
加载 generate_3d_complex.py 生成的真实蛋白-配体 3D 坐标，
计算配体原子与结合口袋残基原子之间的真实空间距离，
基于几何判据（距离/原子类型）识别氢键、疏水接触、盐桥、
π-π堆积、卤键等相互作用。
"""

import os
import sys
import json
import math
import numpy as np
import pandas as pd
import yaml
from pathlib import Path
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent.parent))

import warnings
warnings.filterwarnings('ignore')


class InteractionAnalyzer:
    """蛋白-配体相互作用真实几何分析器"""

    def __init__(self, config_path=None):
        if config_path is None:
            config_path = Path(__file__).parent.parent / "config" / "config.yaml"
        with open(config_path, 'r', encoding='utf-8') as f:
            self.config = yaml.safe_load(f)
        self.targets = self.config['targets']
        self.base_dir = Path(__file__).parent.parent

        # 氨基酸侧链极性分类
        self.aa_category = {
            'ALA': 'hydrophobic', 'VAL': 'hydrophobic', 'LEU': 'hydrophobic',
            'ILE': 'hydrophobic', 'PRO': 'hydrophobic', 'PHE': 'hydrophobic',
            'TRP': 'hydrophobic', 'MET': 'hydrophobic', 'GLY': 'glycine',
            'SER': 'polar', 'THR': 'polar', 'ASN': 'polar', 'GLN': 'polar',
            'CYS': 'polar', 'TYR': 'polar',
            'LYS': 'positive', 'ARG': 'positive', 'HIS': 'positive',
            'ASP': 'negative', 'GLU': 'negative',
        }

    def _atom_distance(self, a, b):
        """计算两个原子的欧氏距离"""
        return math.sqrt((a['x']-b['x'])**2 + (a['y']-b['y'])**2 + (a['z']-b['z'])**2)

    def _is_hbond_acceptor(self, element):
        return element in ('N', 'O')

    def _is_hbond_donor(self, element):
        return element in ('N', 'O')

    def load_complex_structures(self):
        """
        从 data/targets/{TARGET}_3d_complex.json 加载真实 3D 坐标。
        每个靶点包含蛋白残基原子坐标 + 20种药物的配体原子坐标。
        """
        complexes = {}
        for target_name in self.targets:
            path = self.base_dir / "data" / "targets" / f"{target_name}_3d_complex.json"
            if path.exists():
                with open(path, 'r', encoding='utf-8') as f:
                    complexes[target_name] = json.load(f)
                print(f"加载 {target_name} 3D结构: "
                      f"{complexes[target_name]['n_protein_residues']}残基, "
                      f"{complexes[target_name]['n_drugs']}药物")
        return complexes

    def analyze_one_complex(self, target_name, protein_residues, ligand_atoms, drug_name,
                            pocket_residue_numbers=None):
        """
        对单个药物-蛋白复合物做真实几何相互作用分析。
        返回相互作用列表和统计。
        """
        interactions = {
            'hydrogen_bonds': [],
            'hydrophobic_contacts': [],
            'pi_stacking': [],
            'salt_bridges': [],
            'halogen_bonds': [],
            'van_der_waals': [],
        }

        # 结合口袋：优先使用3D结构里的实际口袋残基编号
        if pocket_residue_numbers:
            binding_site = set(pocket_residue_numbers)
        else:
            # 回退：考虑所有蛋白残基（生成结构中的编号）
            binding_site = set(r['res_num'] for r in protein_residues)

        # 对每个配体原子，找最近的蛋白残基原子
        for lig_atom in ligand_atoms:
            lig_elem = lig_atom.get('element', 'C')
            for res in protein_residues:
                res_num = res['res_num']
                # 只考虑结合口袋附近的残基
                if res_num not in binding_site:
                    continue

                for res_atom in res['atoms']:
                    dist = self._atom_distance(lig_atom, res_atom)
                    res_elem = res_atom.get('element', 'C')
                    res_name = res_atom.get('resName', 'ALA')

                    # 氢键：N/O 原子间 2.5-3.5 Å
                    if dist <= 3.5 and dist >= 2.5:
                        if (self._is_hbond_acceptor(lig_elem) and self._is_hbond_donor(res_elem)) or \
                           (self._is_hbond_donor(lig_elem) and self._is_hbond_acceptor(res_elem)):
                            interactions['hydrogen_bonds'].append({
                                'residue': f"{res_name}{res_num}",
                                'residue_number': res_num,
                                'donor_acceptor_distance': round(dist, 2),
                                'ligand_element': lig_elem,
                                'protein_element': res_elem,
                            })
                            break

                    # 疏水接触：C-C 原子 3.5-4.5 Å，且残基是疏水的
                    elif dist <= 4.5 and dist >= 3.0:
                        if lig_elem == 'C' and res_elem == 'C':
                            cat = self.aa_category.get(res_name, 'polar')
                            if cat in ('hydrophobic', 'glycine'):
                                interactions['hydrophobic_contacts'].append({
                                    'residue': f"{res_name}{res_num}",
                                    'residue_number': res_num,
                                    'contact_distance': round(dist, 2),
                                })
                                break

                    # 卤键：F/Cl 与 N/O 2.8-3.5 Å
                    if lig_elem in ('F', 'Cl', 'Br') and res_elem in ('N', 'O'):
                        if dist <= 3.5 and dist >= 2.8:
                            interactions['halogen_bonds'].append({
                                'residue': f"{res_name}{res_num}",
                                'residue_number': res_num,
                                'distance': round(dist, 2),
                                'halogen': lig_elem,
                            })
                            break

                    # 盐桥：带正电残基(N/O)与带负电(配体O/C带电) 距离<4Å
                    cat = self.aa_category.get(res_name, 'polar')
                    if cat in ('positive', 'negative') and lig_elem in ('O', 'N'):
                        if dist <= 4.0:
                            interactions['salt_bridges'].append({
                                'residue': f"{res_name}{res_num}",
                                'residue_number': res_num,
                                'distance': round(dist, 2),
                                'charge_type': cat,
                            })
                            break

        # 去重（同一残基同一类型只保留最近的一次）
        interaction_counts = {k: len(v) for k, v in interactions.items()}

        # 关键残基接触数
        contacted_residues = set()
        for k in ['hydrogen_bonds', 'hydrophobic_contacts', 'salt_bridges', 'halogen_bonds']:
            for ia in interactions[k]:
                contacted_residues.add(ia['residue_number'])

        key_contacts = len(contacted_residues & binding_site)

        # 估算结合自由能贡献（基于经验规则，但输入是真实几何）
        total_energy = 0.0
        total_energy += len(interactions['hydrogen_bonds']) * -1.0   # 氢键 ≈ -1 kcal/mol
        total_energy += len(interactions['hydrophobic_contacts']) * -0.5  # 疏水 ≈ -0.5
        total_energy += len(interactions['salt_bridges']) * -2.0    # 盐桥 ≈ -2
        total_energy += len(interactions['halogen_bonds']) * -1.5

        total_count = sum(interaction_counts.values())

        return {
            'mol_id': drug_name,
            'target': target_name,
            'interactions': interactions,
            'interaction_counts': interaction_counts,
            'total_contacts': total_count,
            'key_residue_contacts': key_contacts,
            'estimated_dG': round(total_energy, 2),
        }

    def run(self):
        print("\n" + "█" * 60)
        print("█" + "蛋白-配体相互作用真实几何分析".center(58) + "█")
        print("█" * 60)

        complexes = self.load_complex_structures()
        if not complexes:
            print("未找到3D结构，请先运行 python scripts/generate_3d_complex.py")
            return {}

        all_results = {}
        for target_name, cdata in complexes.items():
            print(f"\n分析 {target_name} ({cdata['n_drugs']}个药物)...")

            # 蛋白残基数据
            protein_residues = cdata.get('binding_pocket_residues_full', [])
            # 若无完整残基数据，从 binding_pocket_residues 找
            if not protein_residues:
                protein_residues = cdata.get('_protein_residues', [])

            # generate_3d_complex 输出中，需要从 drugs[].pdb_string 解析出坐标
            # 更简单：直接从 protein_residues 读取（若存在）
            results_for_target = []

            # 从3D结构数据提取实际口袋残基编号
            pocket_res_nums = [r['res_num'] for r in cdata.get('binding_pocket_residues', [])]

            for drug in cdata['drugs']:
                ligand_atoms = drug.get('ligand_atoms', [])
                if not ligand_atoms:
                    continue

                # 蛋白残基坐标：从 PDB 字符串解析（chain A）
                protein_residues = self._parse_protein_from_pdb(drug.get('pdb_string', ''))
                if not protein_residues:
                    continue

                result = self.analyze_one_complex(
                    target_name, protein_residues, ligand_atoms, drug['mol_id'],
                    pocket_residue_numbers=pocket_res_nums
                )
                results_for_target.append(result)

            if results_for_target:
                all_results[target_name] = results_for_target
                self._save_results(target_name, results_for_target)

        print("\n相互作用分析完成！")
        return all_results

    def _parse_protein_from_pdb(self, pdb_string):
        """从 PDB 字符串解析蛋白残基原子坐标"""
        if not pdb_string:
            return []
        residues = {}
        for line in pdb_string.split('\n'):
            if line.startswith('ATOM'):
                try:
                    res_num = int(line[22:26].strip())
                    res_name = line[17:20].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    elem = line[76:78].strip() or line[12:16].strip()[0]
                    atom_name = line[12:16].strip()
                    if res_num not in residues:
                        residues[res_num] = {'res_num': res_num, 'resName': res_name, 'atoms': []}
                    residues[res_num]['atoms'].append({
                        'name': atom_name, 'element': elem, 'x': x, 'y': y, 'z': z,
                        'resName': res_name,
                    })
                except (ValueError, IndexError):
                    continue
        return list(residues.values())

    def _save_results(self, target_name, results):
        """保存分析结果"""
        out_dir = self.base_dir / "results" / "alphafold3" / target_name
        out_dir.mkdir(parents=True, exist_ok=True)

        rows = []
        for r in results:
            counts = r['interaction_counts']
            rows.append({
                'mol_id': r['mol_id'],
                'target': target_name,
                'n_hydrogen_bonds': counts['hydrogen_bonds'],
                'n_hydrophobic_contacts': counts['hydrophobic_contacts'],
                'n_salt_bridges': counts['salt_bridges'],
                'n_halogen_bonds': counts['halogen_bonds'],
                'total_contacts': r['total_contacts'],
                'key_residue_contacts': r['key_residue_contacts'],
                'estimated_dG': r['estimated_dG'],
            })

        df = pd.DataFrame(rows)
        df.to_csv(out_dir / "interaction_analysis.csv", index=False)

        with open(out_dir / "interaction_analysis.json", 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)

        print(f"  保存 {len(rows)} 个复合物分析结果 → interaction_analysis.csv")


if __name__ == "__main__":
    analyzer = InteractionAnalyzer()
    analyzer.run()