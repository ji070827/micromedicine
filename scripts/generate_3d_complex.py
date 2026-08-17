#!/usr/bin/env python3
"""
generate_3d_complex.py - 生成蛋白-配体复合物的3D结构数据
用于3Dmol.js前端可视化，输出PDB格式坐标
使用 IgV β-三明治拓扑模板生成紧凑可靠的蛋白结构
支持为每个靶点生成多个候选小分子药物的模拟结合构象
"""

import os
import sys
import json
import math
import random
import hashlib
import numpy as np
from pathlib import Path
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent.parent))

import warnings
warnings.filterwarnings('ignore')

# ========== 氨基酸映射 ==========
AA_MAP = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
}

AA_COLORS = {
    'hydrophobic': '#FFD700',
    'polar': '#00FF88',
    'positive': '#4488FF',
    'negative': '#FF4444',
    'glycine': '#AAAAAA',
}

def get_residue_category(aa_code):
    if aa_code in 'ALVIMFWPG': return 'hydrophobic'
    elif aa_code in 'STNQCY': return 'polar'
    elif aa_code in 'KRH': return 'positive'
    elif aa_code in 'DE': return 'negative'
    return 'polar'

# ========== 各靶点序列（IgV域 β-三明治） ==========
TARGET_POCKET_SEQUENCES = {
    "PD-1": "LNWYRMSPSNQTDKLAAFPEDRSQPGQDCRFRVTQLPNGRDFHMSVVRARRNDSGTYLCGAISLAPKAQIKESLRAELRVTERRAE",
    "LAG-3": "LSLRRAGVTWQHQPDSGPPAAAPGHPLAPGPHPAAPSSWGPRPRRYTVLSVGPGGLRSGRLPLQPRVQLDERGRQRGDF",
    "TIM-3": "KGACPVFECGNVVLRTDERDVNYWTSRYWLNGDFRKGDVSLTIENVTLADSGIYCCRIQIPGIMNDEKFNLK",
    "VISTA": "WYRSSRGEVQTCSERRPIRNLTFQDLHLHHGGHQAANTSHDLAQRHGLESASDHHGNFSITMRNLTLLDSGLY",
}

# 结合口袋残基索引（偏移量，从START_RES开始算）
TARGET_BINDING_SITE_OFFSETS = {
    "PD-1": [4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21,22],
    "LAG-3": [4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],
    "TIM-3": [4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],
    "VISTA": [4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22],
}

START_RES = 1  # 蛋白残基起始编号

# ========== IgV β-三明治 拓扑坐标模板 ==========
# 基于典型 IgV 域的原子级坐标模板（归一化到≈25x35x40Å范围）
# 9条β链: A', G, F, C, C', C" (前片) | B, E, D (后片)
# 残基对标: 从N端到C端约80残基

# 每条β链的 CA 原子坐标（Å）— 9条链，各含3~6个CA位置
# 前片 (front sheet) A'-G-F-C-C'
IGV_FRONT_STRANDS = {
    "A":  [(0.5, -15.0, 4.0), (0.8, -15.0, 7.5), (0.3, -15.0, 11.0), (0.0, -15.0, 14.5)],  # A' strand
    "G":  [(9.5, -15.0, 2.0), (9.2, -15.0, 5.5), (9.5, -15.0, 9.0), (10.0, -15.0, 12.5), (10.3, -15.0, 16.0)],  # G strand
    "F":  [(15.5, -14.0, 2.0), (15.2, -14.0, 5.5), (15.0, -14.0, 9.0), (15.5, -14.0, 12.5)],  # F strand
    "CC": [(21.0, -13.0, 3.0), (20.5, -13.0, 6.5), (20.2, -13.0, 10.0), (20.8, -13.0, 13.5), (21.5, -13.0, 17.0)],  # C' strand
    "C2": [(25.0, -11.0, 5.0), (24.5, -11.0, 8.5), (24.8, -11.0, 12.0)],  # C strand
}

# 后片 (back sheet) B-E-D
IGV_BACK_STRANDS = {
    "B":  [(2.0, 2.0, 2.0), (2.5, 2.0, 5.5), (2.8, 2.0, 9.0), (2.3, 2.0, 12.5)],  # B strand
    "E":  [(13.0, 3.0, 1.0), (12.5, 3.0, 4.5), (12.8, 3.0, 8.0), (13.3, 3.0, 11.5), (13.0, 3.0, 15.0)],  # E strand
    "D":  [(18.5, 4.0, 2.5), (18.0, 4.0, 6.0), (18.3, 4.0, 9.5), (18.8, 4.0, 13.0)],  # D strand
}

# 连接环 (loops) — 桥接各条β链
IGV_LOOPS = [
    [(1.0, -6.0, 16.0), (2.0, -5.0, 17.5), (3.5, -8.0, 18.0)],  # A'-B loop
    [(4.0, -8.0, 14.0), (5.5, -9.0, 15.0), (6.0, -10.0, 13.0)],  # B-E loop
    [(9.0, -10.0, 17.0), (8.5, -12.0, 18.0), (7.0, -13.0, 17.5)],  # E-D loop
    [(15.0, -10.0, 14.0), (16.5, -11.0, 15.5), (17.0, -13.0, 16.0)],  # D-F loop
    [(18.5, -12.0, 18.0), (19.0, -15.0, 18.5), (20.5, -14.0, 20.0)],  # F-G loop
    [(23.0, -11.0, 15.0), (24.0, -13.0, 15.5), (25.5, -14.0, 14.0)],  # G-C' loop
    [(26.5, -9.0, 18.0), (27.0, -10.0, 19.5), (27.5, -12.0, 20.0)],  # C'-C loop  ** ← 此处形成 CDR-like 结合口袋 **
]

# 将拓扑序列映射为连续的CA坐标路径
# 残基顺序: A'(4) → loop(3) → B(4) → loop(3) → E(5) → loop(3) → D(4) → loop(3) → F(4) → loop(3) → G(5) → loop(3) → C'(5) → loop(3) → C(3)
IGV_TOPOLOGY_ORDER = [
    ("A", 4), ("loop", 0),  # A' strand + AB loop
    ("B", 4), ("loop", 1),  # B strand + BE loop
    ("E", 5), ("loop", 2),  # E strand + ED loop  
    ("D", 4), ("loop", 3),  # D strand + DF loop
    ("F", 4), ("loop", 4),  # F strand + FG loop
    ("G", 5), ("loop", 5),  # G strand + GC' loop
    ("CC", 5), ("loop", 6),  # C' strand + C'C loop (口袋区域)
    ("C2", 3),                # C strand
]

def build_igv_ca_coords():
    """构建 IgV 域的完整 CA 坐标列表"""
    coords = []
    strand_idx = {"A": 0, "B": 0, "E": 0, "D": 0, "F": 0, "G": 0, "CC": 0, "C2": 0}
    loop_idx = 0
    
    for item_type, item_id in IGV_TOPOLOGY_ORDER:
        if item_type == "loop":
            loops = IGV_LOOPS
            idx = item_id if isinstance(item_id, int) else loop_idx
            idx = min(idx, len(loops) - 1)
            for pt in loops[idx]:
                coords.append(pt)
            loop_idx += 1
        else:
            strands = {**IGV_FRONT_STRANDS, **IGV_BACK_STRANDS}
            pts = strands[item_type]
            n = item_id
            if isinstance(n, int):
                pts = pts[:n]
            for pt in pts:
                coords.append(pt)
    
    return coords


def generate_protein_coords_from_topology(sequence, target_name, start_res=1):
    """
    从 IgV 拓扑模板分配残基到 CA 坐标位置
    每个 CA 位置分配序列中的一个残基
    """
    ca_coords = build_igv_ca_coords()
    residues = []
    
    # 对靶点做微小扰动以区分
    seed_val = hash(target_name + "_igv") % 2**32
    rng = random.Random(seed_val)
    
    # 处理序列长度与CA数量的差异
    n_ca = len(ca_coords)
    n_seq = len(sequence)
    
    if n_seq < n_ca:
        # 序列较短：截断CA坐标
        ca_coords = ca_coords[:n_seq]
    elif n_seq > n_ca:
        # 序列较长：在最后几个loop处插入额外CA（C端延伸）
        extra = n_seq - n_ca
        last = ca_coords[-1] if ca_coords else (30.0, -12.0, 18.0)
        for k in range(extra):
            ca_coords.append((last[0] + 1.5, last[1] + k * 1.2, last[2] + rng.uniform(-1, 1)))
    
    n_atoms = len(ca_coords)
    
    # 为每个CA坐标位置生成主链+CB原子
    n_offset_v = np.array([-0.52, 0.30, 0.88])
    c_offset_v = np.array([0.52, -0.30, -0.44])
    o_offset_v = np.array([1.24, -0.40, -0.66])
    cb_offset_v = np.array([-0.54, -0.69, -0.50])  # 用于非Gly
    
    for i in range(n_atoms):
        ca = np.array(ca_coords[i]) + np.array([rng.uniform(-0.2, 0.2) for _ in range(3)])
        
        res_num = start_res + i
        aa_one = sequence[i] if i < n_seq else 'A'
        aa_three = AA_MAP.get(aa_one, 'ALA')
        category = get_residue_category(aa_one)
        
        n_atom = ca + n_offset_v
        c_atom = ca + c_offset_v
        o_atom = ca + o_offset_v
        cb_atom = ca + cb_offset_v if aa_one != 'G' else ca + np.array([0.0, 0.0, 0.0])
        
        atoms = [
            {"serial": i * 5 + 1, "name": "N", "alt": " ", "resName": aa_three, "chain": "A",
             "resSeq": res_num, "x": n_atom[0], "y": n_atom[1], "z": n_atom[2], "element": "N"},
            {"serial": i * 5 + 2, "name": "CA", "alt": " ", "resName": aa_three, "chain": "A",
             "resSeq": res_num, "x": ca[0], "y": ca[1], "z": ca[2], "element": "C"},
            {"serial": i * 5 + 3, "name": "C", "alt": " ", "resName": aa_three, "chain": "A",
             "resSeq": res_num, "x": c_atom[0], "y": c_atom[1], "z": c_atom[2], "element": "C"},
            {"serial": i * 5 + 4, "name": "O", "alt": " ", "resName": aa_three, "chain": "A",
             "resSeq": res_num, "x": o_atom[0], "y": o_atom[1], "z": o_atom[2], "element": "O"},
            {"serial": i * 5 + 5, "name": "CB", "alt": " ", "resName": aa_three, "chain": "A",
             "resSeq": res_num, "x": cb_atom[0], "y": cb_atom[1], "z": cb_atom[2], "element": "C"},
        ]
        
        residue = {
            "res_num": res_num,
            "aa_one": aa_one,
            "aa_three": aa_three,
            "category": category,
            "color": AA_COLORS.get(category, '#CCCCCC'),
            "atoms": atoms,
        }
        residues.append(residue)
    
    return residues


# ========== 3D 药物模板 (带真实 z 坐标) ==========
DRUG_MOLECULE_TEMPLATES = {
    "mol_001": {
        "name": "BMS-936558 Analog (Aryl Ether)",
        "mol_type": "Aryl Ether",
        "smiles_like": "CCOC1=CC=C2NC3=CC=CC=C3C2=C1",
        "atoms": [
            ("C1", "C", -1.5, 0.0, 0.2), ("C2", "C", -0.6, 1.2, 0.5),
            ("C3", "C", 0.8, 1.2, 0.3), ("C4", "C", 1.6, 0.0, -0.2),
            ("C5", "C", 0.8, -1.2, -0.5), ("C6", "C", -0.6, -1.2, -0.3),
            ("C7", "C", -2.9, 0.0, 0.4), ("O1", "O", -3.5, 1.2, 0.7),
            ("C8", "C", -4.9, 1.2, 0.9), ("C9", "C", -5.5, 0.0, 0.8),
            ("N1", "N", -1.3, 2.5, 0.9), ("C10", "C", 2.9, 2.3, 0.5),
            ("C11", "C", 3.6, 3.5, 0.7), ("C12", "C", 4.9, 3.5, 0.5),
            ("O2", "O", 5.6, 2.3, 0.1), ("C13", "C", 3.0, -2.3, -0.9),
            ("F1", "F", 3.5, -3.3, -1.3),
        ],
        "bonds": [(0,1),(1,2),(2,3),(3,4),(4,5),(5,0),(0,6),(6,7),(7,8),(8,9),(1,10),(2,11),(11,12),(12,13),(4,14),(14,15)],
        "binding_affinity": -9.2, "mw": 380.4, "drug_likeness": 0.82,
        "predicted_interactions": {"hydrogen_bonds": 3, "pi_stacking": 2, "hydrophobic": 5}
    },
    "mol_002": {
        "name": "LAG-525 Analog (Benzamide)",
        "mol_type": "Benzamide",
        "smiles_like": "O=C(NCC1=CC=CC=C1)C2=CC=C(OCCN3CCCC3)C=C2",
        "atoms": [
            ("C1", "C", -2.0, 0.0, 0.15), ("C2", "C", -0.5, 0.0, -0.15),
            ("C3", "C", 1.0, 0.0, 0.1), ("C4", "C", 2.5, 0.0, -0.1),
            ("C5", "C", 3.2, 1.2, 0.2), ("C6", "C", 4.6, 1.2, -0.15),
            ("C7", "C", 5.3, 0.0, 0.3), ("C8", "C", 4.6, -1.2, -0.3),
            ("C9", "C", 3.2, -1.2, 0.2), ("O1", "O", -3.5, 0.0, -0.5),
            ("N1", "N", 1.7, -1.3, -0.5), ("C10", "C", 3.1, -1.8, 0.5),
            ("C11", "C", 3.8, -3.1, 0.3), ("N2", "N", 5.2, -3.1, -0.4),
            ("C12", "C", 5.9, -1.8, -0.8), ("C13", "C", 5.9, -4.4, 0.6),
            ("C14", "C", 7.3, -4.4, -0.3),
        ],
        "bonds": [(0,1),(0,9),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,3),(2,10),(10,11),(11,12),(12,13),(12,14)],
        "binding_affinity": -8.7, "mw": 345.4, "drug_likeness": 0.78,
        "predicted_interactions": {"hydrogen_bonds": 2, "pi_stacking": 1, "hydrophobic": 4}
    },
    "mol_003": {
        "name": "TSR-042 Analog (Benzimidazole)",
        "mol_type": "Benzimidazole",
        "smiles_like": "CC1=NC2=CC=CC=C2N1CC3=CC=C(C=C3)C(=O)NO",
        "atoms": [
            ("C1", "C", -2.5, 0.0, 0.2), ("C2", "C", -1.0, 0.0, -0.1),
            ("C3", "C", 0.5, 0.0, 0.3), ("C4", "C", 2.0, 0.0, -0.2),
            ("C5", "C", 2.7, 1.2, 0.4), ("C6", "C", 4.1, 1.2, -0.2),
            ("C7", "C", 2.7, -1.2, -0.5), ("C8", "C", 4.1, -1.2, 0.3),
            ("N1", "N", 4.8, 0.0, -0.5), ("N2", "N", 0.5, 1.3, -0.4),
            ("C9", "C", -0.2, 2.5, 0.6), ("C10", "C", -1.6, 2.5, -0.3),
            ("C11", "C", -2.3, 1.3, 0.4), ("C12", "C", -1.6, 3.7, -0.7),
            ("C13", "C", -3.0, 3.7, 0.5), ("C14", "C", -3.7, 2.5, -0.4),
            ("O1", "O", -3.9, 0.0, 0.7), ("F1", "F", -4.2, 1.3, -0.5),
            ("C15", "C", -5.0, 0.0, -0.3),
        ],
        "bonds": [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(3,7),(7,8),(8,4),(2,9),(9,10),(10,11),(11,12),(12,13),(13,14),(0,15),(15,17)],
        "binding_affinity": -9.5, "mw": 310.3, "drug_likeness": 0.85,
        "predicted_interactions": {"hydrogen_bonds": 3, "pi_stacking": 2, "hydrophobic": 3}
    },
    "mol_004": {
        "name": "CA-170 Analog (Biphenyl Urea)",
        "mol_type": "Biphenyl Urea",
        "smiles_like": "CN(C)CCOC1=CC=C(C=C1)C(=O)NC2=CC=CC=C2",
        "atoms": [
            ("C1", "C", -3.0, 0.0, 0.3), ("C2", "C", -1.5, 0.0, -0.2),
            ("C3", "C", 0.0, 0.0, 0.1), ("C4", "C", 1.5, 0.0, -0.3),
            ("C5", "C", 3.0, 0.0, 0.2), ("C6", "C", 3.7, 1.2, -0.3),
            ("C7", "C", 5.1, 1.2, 0.2), ("C8", "C", 5.8, 0.0, -0.4),
            ("C9", "C", 5.1, -1.2, 0.3), ("C10","C", 3.7, -1.2, -0.2),
            ("N1", "N", 1.5, -1.3, 0.5), ("O1", "O", 0.0, 1.3, -0.5),
            ("C11","C", 0.0, 2.7, 0.4), ("C12","C", 1.0, 3.5, -0.3),
            ("C13","C", 1.0, 4.9, 0.5), ("C14","C", 0.0, 5.5, -0.4),
            ("N2", "N", 2.0, 5.7, 0.6), ("C15","C", 3.3, 5.1, -0.3),
            ("C16","C", 4.2, 5.9, 0.5),
        ],
        "bonds": [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,9),(9,4),(3,10),(2,11),(11,12),(12,13),(12,14),(14,15),(15,16)],
        "binding_affinity": -8.4, "mw": 355.4, "drug_likeness": 0.76,
        "predicted_interactions": {"hydrogen_bonds": 2, "pi_stacking": 2, "hydrophobic": 4}
    },
    "mol_005": {
        "name": "SHR-1210 Analog (Quinoline)",
        "mol_type": "Quinoline",
        "smiles_like": "CN1CCC(CC1)OC2=CC=C3C(=C2)C=CC=N3",
        "atoms": [
            ("C1", "C", -3.5, 0.0, 0.1), ("C2", "C", -2.0, 0.0, -0.2),
            ("C3", "C", -0.5, 0.0, 0.3), ("C4", "C", 1.0, 0.0, -0.1),
            ("C5", "C", 2.5, 0.0, 0.2), ("C6", "C", 3.2, 1.2, -0.3),
            ("C7", "C", 4.6, 1.2, 0.1), ("C8", "C", 5.3, 0.0, -0.2),
            ("C9", "C", 4.6, -1.2, 0.3), ("C10","C", 3.2, -1.2, -0.1),
            ("N1", "N", 6.7, 0.0, 0.5), ("C11","C", 7.4, 1.2, -0.3),
            ("C12","C", 8.8, 1.2, 0.2), ("C13","C", 9.5, 0.0, -0.4),
            ("C14","C", 8.8, -1.2, 0.3), ("C15","C", 7.4, -1.2, -0.2),
            ("O1", "O", -4.2, 1.2, -0.5), ("C16","C", -5.6, 1.2, 0.4),
            ("N2", "N", -6.3, 0.0, -0.4), ("C17","C", -7.7, 0.0, 0.3),
        ],
        "bonds": [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,9),(9,4),(7,10),(10,11),(11,12),(12,13),(13,14),(14,15),(0,16),(16,17),(17,18)],
        "binding_affinity": -8.9, "mw": 345.4, "drug_likeness": 0.73,
        "predicted_interactions": {"hydrogen_bonds": 2, "pi_stacking": 3, "hydrophobic": 5}
    },
    "mol_006": {
        "name": "Pembrolizumab Small Molecule Mimetic",
        "mol_type": "Indole-piperidine",
        "smiles_like": "COC1=CC2=C(C=C1)N(C=C2)CCN3CCC(CC3)O",
        "atoms": [
            ("C1", "C", -2.5, 0.0, 0.3), ("C2", "C", -1.0, 0.0, -0.2),
            ("C3", "C", 0.5, 0.0, 0.1), ("C4", "C", 2.0, 0.0, -0.4),
            ("C5", "C", 2.7, 1.2, 0.3), ("C6", "C", 4.1, 1.2, -0.2),
            ("C7", "C", 4.8, 0.0, 0.5), ("C8", "C", 4.1, -1.2, -0.3),
            ("C9", "C", 2.7, -1.2, 0.2), ("N1", "N", 0.5, 1.3, -0.5),
            ("C10","C", -0.2, 2.5, 0.4), ("C11","C", -1.6, 2.5, -0.3),
            ("C12","C", -2.3, 1.3, 0.5), ("C13","C", -2.3, 3.7, -0.6),
            ("C14","C", -3.7, 3.7, 0.3), ("C15","C", -4.4, 2.5, -0.4),
            ("O1", "O", -5.8, 2.5, 0.5), ("C16","C", -6.5, 1.3, -0.6),
            ("O2", "O", 6.2, 0.0, -0.5),
        ],
        "bonds": [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,3),(2,9),(9,10),(10,11),(11,12),(11,13),(13,14),(14,15),(15,16),(16,17),(6,18)],
        "binding_affinity": -10.1, "mw": 326.4, "drug_likeness": 0.81,
        "predicted_interactions": {"hydrogen_bonds": 4, "pi_stacking": 2, "hydrophobic": 3}
    },
    "mol_007": {
        "name": "MBG453 Analog (Thiadiazole)",
        "mol_type": "Thiadiazole",
        "smiles_like": "CC1=NN=C(S1)NCC2=CC=C(C=C2)OC(F)(F)F",
        "atoms": [
            ("C1", "C", -3.0, 0.5, 0.2), ("C2", "C", -1.5, 0.0, -0.3),
            ("C3", "C", 0.0, 0.0, 0.1), ("C4", "C", 1.5, 0.0, -0.4),
            ("C5", "C", 3.0, 0.0, 0.3), ("C6", "C", 3.7, 1.2, -0.2),
            ("C7", "C", 5.1, 1.2, 0.5), ("C8", "C", 5.8, 0.0, -0.3),
            ("C9", "C", 5.1, -1.2, 0.2), ("C10","C", 3.7, -1.2, -0.4),
            ("S1", "S", 0.0, 1.5, -0.5), ("N1", "N", 1.2, 1.8, 0.3),
            ("N2", "N", 2.0, 0.8, -0.4), ("O1", "O", 7.2, 0.0, 0.5),
            ("F1", "F", 7.8, 0.8, -0.6), ("F2", "F", 7.8, -0.8, 0.6),
            ("F3", "F", 8.0, 0.0, -0.3), ("C11","C", -4.5, 0.0, 0.4),
            ("N3", "N", -5.2, 1.2, -0.5),
        ],
        "bonds": [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,9),(9,4),(2,10),(10,11),(11,12),(7,13),(13,14),(13,15),(0,16),(16,17)],
        "binding_affinity": -9.3, "mw": 358.3, "drug_likeness": 0.79,
        "predicted_interactions": {"hydrogen_bonds": 3, "pi_stacking": 1, "hydrophobic": 5, "halogen": 3}
    },
    "mol_008": {
        "name": "BMS-986207 Analog (Triazolopyrimidine)",
        "mol_type": "Triazolopyrimidine",
        "smiles_like": "CN1N=NC2=C1N=C(N=C2)NC3=CC=C(C=C3)C#N",
        "atoms": [
            ("C1", "C", -2.5, 0.0, 0.2), ("C2", "C", -1.0, 0.0, -0.3),
            ("C3", "C", 0.5, 0.0, 0.4), ("C4", "C", 2.0, 0.0, -0.2),
            ("C5", "C", 3.5, 0.0, 0.3), ("C6", "C", 4.2, 1.2, -0.4),
            ("C7", "C", 5.6, 1.2, 0.2), ("C8", "C", 6.3, 0.0, -0.5),
            ("C9", "C", 5.6, -1.2, 0.3), ("C10","C", 4.2, -1.2, -0.2),
            ("N1", "N", 2.0, 1.3, -0.5), ("N2", "N", 3.3, 1.3, 0.4),
            ("N3", "N", 4.0, 0.0, -0.4), ("N4", "N", 6.8, 0.0, 0.6),
            ("C11","C", 7.5, -1.2, -0.8), ("C12","C", 8.0, 0.0, 0.5),
            ("C13","C", 7.5, 1.2, -0.6), ("C14","C", 6.8, 0.0, -1.2),
        ],
        "bonds": [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,9),(9,4),(3,10),(10,11),(11,12),(12,13),(13,14),(6,15),(15,16)],
        "binding_affinity": -9.8, "mw": 370.4, "drug_likeness": 0.83,
        "predicted_interactions": {"hydrogen_bonds": 4, "pi_stacking": 3, "hydrophobic": 4}
    },
    "mol_009": {
        "name": "TSR-033 Analog (Pyrazolopyrimidine)",
        "mol_type": "Pyrazolopyrimidine",
        "smiles_like": "CN1C=CC2=C1N=CN=C2NCC3=CC=C(C=C3)F",
        "atoms": [
            ("C1", "C", -2.0, 0.0, 0.2), ("C2", "C", -0.5, 0.0, -0.3),
            ("C3", "C", 1.0, 0.0, 0.4), ("C4", "C", 2.5, 0.0, -0.2),
            ("C5", "C", 4.0, 0.0, 0.3), ("C6", "C", 4.7, 1.2, -0.4),
            ("C7", "C", 6.1, 1.2, 0.2), ("C8", "C", 6.8, 0.0, -0.5),
            ("C9", "C", 6.1, -1.2, 0.3), ("C10","C", 4.7, -1.2, -0.2),
            ("N1", "N", 2.5, 1.3, -0.6), ("N2", "N", 3.8, 1.3, 0.4),
            ("N3", "N", 4.5, 0.0, -0.5), ("C11","C", 3.8, -1.3, 0.4),
            ("C12","C", 2.5, -1.3, -0.3), ("F1", "F", 8.2, 0.0, 0.6),
            ("N4", "N", 0.0, 1.3, -0.5), ("O1", "O", -2.7, 1.2, 0.4),
            ("C13","C", -4.1, 1.2, -0.3),
        ],
        "bonds": [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,9),(9,4),(3,10),(10,11),(11,12),(12,13),(4,13),(0,16),(16,17)],
        "binding_affinity": -8.6, "mw": 340.3, "drug_likeness": 0.75,
        "predicted_interactions": {"hydrogen_bonds": 2, "pi_stacking": 2, "hydrophobic": 3, "halogen": 1}
    },
    "mol_010": {
        "name": "AB122 Analog (Pyridine-Benzamide)",
        "mol_type": "Pyridine-Benzamide",
        "smiles_like": "CN1CCN(CC1)C2=CC=C(C=C2)C(=O)NC3=CC=NC=C3",
        "atoms": [
            ("C1", "C", -2.5, 0.0, 0.3), ("C2", "C", -1.0, 0.0, -0.2),
            ("C3", "C", 0.5, 0.0, 0.1), ("C4", "C", 2.0, 0.0, -0.4),
            ("C5", "C", 3.5, 0.0, 0.2), ("C6", "C", 4.2, 1.2, -0.3),
            ("C7", "C", 5.6, 1.2, 0.5), ("C8", "C", 6.3, 0.0, -0.4),
            ("C9", "C", 5.6, -1.2, 0.3), ("C10","C", 4.2, -1.2, -0.2),
            ("N1", "N", 7.7, 0.0, -0.5), ("N2", "N", 2.0, 1.3, 0.5),
            ("O1", "O", 2.0, -1.3, -0.5), ("C11","C", 3.3, 1.8, 0.4),
            ("C12","C", 3.3, 3.2, -0.3), ("N3", "N", 4.6, 3.7, 0.5),
            ("C13","C", 5.8, 2.7, -0.6), ("C14","C", 5.8, 1.3, 0.4),
            ("O2", "O", -3.2, 1.2, -0.5), ("C15","C", -4.6, 1.2, 0.4),
        ],
        "bonds": [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,9),(9,4),(7,10),(3,11),(3,12),(11,13),(13,14),(14,15),(0,16),(16,17)],
        "binding_affinity": -9.0, "mw": 325.4, "drug_likeness": 0.88,
        "predicted_interactions": {"hydrogen_bonds": 3, "pi_stacking": 2, "hydrophobic": 3}
    },
}

MORE_SIMPLE_LIGANDS = [
    {"name": "4-(4-Fluorobenzyl)piperidine", "smiles": "Fc1ccc(CC2CCNCC2)cc1"},
    {"name": "Dopamine Analog", "smiles": "CC(N)Cc1ccc(O)c(O)c1"},
    {"name": "Sulfonamide-thiazole", "smiles": "Nc1ccc(S(=O)(=O)Nc2nccs2)cc1"},
    {"name": "Cyclopropyl-carboxamide", "smiles": "O=C(NC1CC1)c2cccc3c2CCN(C)C3"},
    {"name": "t-Butyl benzamide", "smiles": "CC(C)(C)c1ccc(C(=O)NCc2ccccc2)cc1"},
    {"name": "N-methyl piperidine amide", "smiles": "CNC(=O)C1CCN(CC1)c2ccccc2"},
    {"name": "Diethoxy benzamide", "smiles": "CCOc1cc(OCC)cc(C(=O)Nc2ncccn2)c1"},
    {"name": "Thiadiazole thioether", "smiles": "Cc1nnc(SCC(=O)Nc2cccc(C)c2)s1"},
    {"name": "Tetrahydroisoquinoline furan", "smiles": "COc1ccc2c(c1)CCN(C2)C(=O)c3ccco3"},
    {"name": "Morpholine benzamide", "smiles": "NC(=O)c1ccc(CN2CCOCC2)cc1"},
]


class Complex3DGenerator:
    def __init__(self):
        self.base_dir = Path(__file__).parent.parent
        self.output_dir = self.base_dir / "data" / "targets"

    def generate_protein_coords(self, sequence, target_name, start_res=1):
        """使用 IgV β-三明治拓扑模板生成蛋白主链坐标"""
        return generate_protein_coords_from_topology(sequence, target_name, start_res)

    def generate_ligand_from_template(self, template_data, center, scale=1.0, mol_id="mol_001"):
        """从模板生成配体3D坐标，放置在口袋中心"""
        atoms_list = template_data.get("atoms", [])
        bonds_list = template_data.get("bonds", [])

        cx, cy, cz = center
        atoms = []

        seed_val = hash(mol_id + str(center)) % 2**32
        rng = random.Random(seed_val)
        rot_angle = rng.uniform(0, 2 * math.pi)
        
        # 小幅偏移使不同药物在口袋中有不同位置
        offset_x = rng.uniform(-0.5, 0.5)
        offset_y = rng.uniform(-0.5, 0.5)
        offset_z = rng.uniform(-0.5, 0.5)

        for i, atom_info in enumerate(atoms_list):
            name, elem = atom_info[0], atom_info[1]
            tx, ty = atom_info[2], atom_info[3]
            tz = atom_info[4] if len(atom_info) >= 5 else 0.0

            rx = tx * math.cos(rot_angle) - ty * math.sin(rot_angle)
            ry = tx * math.sin(rot_angle) + ty * math.cos(rot_angle)

            atom = {
                "serial": i + 1,
                "name": name,
                "alt": " ",
                "resName": "LIG",
                "chain": "B",
                "resSeq": 1,
                "x": round(cx + rx * scale + offset_x, 3),
                "y": round(cy + ry * scale + offset_y, 3),
                "z": round(cz + tz * scale + offset_z, 3),
                "element": elem,
            }
            atoms.append(atom)

        return {"atoms": atoms, "bonds": bonds_list}

    def generate_simple_ligand(self, center, mol_id, seed_offset=0):
        """为简单分子生成3D坐标"""
        cx, cy, cz = center
        seed = hash(mol_id + str(seed_offset)) % 2**32
        rng = random.Random(seed)
        natoms = rng.randint(8, 14)

        atoms = []
        for i in range(natoms):
            angle = (i / natoms) * 2 * math.pi + rng.uniform(-0.3, 0.3)
            radius = rng.uniform(1.5, 3.5)
            elem = rng.choice(['C', 'C', 'C', 'N', 'O', 'C', 'C', 'S', 'F'])
            atoms.append({
                "serial": i + 1,
                "name": f"{elem}{i+1}",
                "alt": " ",
                "resName": "LIG",
                "chain": "B",
                "resSeq": 1,
                "x": round(cx + radius * math.cos(angle), 3),
                "y": round(cy + rng.uniform(-2.5, 2.5), 3),
                "z": round(cz + radius * math.sin(angle), 3),
                "element": elem,
            })

        bonds = []
        for j in range(natoms - 1):
            bonds.append((j, j + 1))
        if natoms > 4:
            bonds.append((0, natoms - 1))

        return {"atoms": atoms, "bonds": bonds}

    def to_pdb_single_ligand(self, protein_residues, ligand_data):
        """转换为含单个配体的PDB格式文本"""
        lines = []
        lines.append("HEADER    PROTEIN-LIGAND COMPLEX")
        lines.append("TITLE     IgV Domain Drug-Binding Complex")
        lines.append("REMARK    Simulated structure based on IgV beta-sandwich topology")
        lines.append("REMARK    Protein chain A, Ligand chain B")

        # 蛋白原子
        for res in protein_residues:
            for atom in res["atoms"]:
                line = self._format_atom_line(atom, is_hetatm=False)
                lines.append(line)

        # 配体原子 (HETATM, chain B)
        for atom in ligand_data["atoms"]:
            a = atom.copy()
            a["chain"] = "B"
            line = self._format_atom_line(a, is_hetatm=True)
            lines.append(line)

        # 连接记录
        for bond in ligand_data["bonds"]:
            a1 = ligand_data["atoms"][bond[0]]["serial"]
            a2 = ligand_data["atoms"][bond[1]]["serial"]
            lines.append(f"CONECT{a1:5d}{a2:5d}")

        lines.append("END")
        return "\n".join(lines)

    def _format_atom_line(self, atom, is_hetatm=False):
        """Format a single ATOM/HETATM line complying with PDB v3 format"""
        record = "HETATM" if is_hetatm else "ATOM  "
        serial = atom.get("serial", 1)
        name = atom.get("name", "CA")
        alt = atom.get("alt", " ")
        res = atom.get("resName", "ALA")
        chain = atom.get("chain", "A")
        res_seq = atom.get("resSeq", 1)
        x = atom.get("x", 0.0)
        y = atom.get("y", 0.0)
        z = atom.get("z", 0.0)
        elem = atom.get("element", "C")

        # PDB格式: ATOM/HETATM(1-6) serial(7-11) 空格(12) name(13-16) altLoc(17) resName(18-20) 空格(21) chainID(22) resSeq(23-26) iCode(27) 空格(28-30) x(31-38) y(39-46) z(47-54) occupancy(55-60) tempFactor(61-66) 空格(67-76) element(77-78) charge(79-80)
        # 简化版
        line = (
            f"{record:<6s}"
            f"{serial:5d} "
            f"{name:<4s}"
            f"{alt}"
            f"{res:<3s} "
            f"{chain}"
            f"{res_seq:4d}    "
            f"{x:8.3f}"
            f"{y:8.3f}"
            f"{z:8.3f}"
            f"  1.00  0.00          "
            f"{elem:<2s}"
        )
        return line

    def generate_complex_with_drugs(self, target_name):
        """为单个靶点生成蛋白+所有候选药物的复合物数据"""
        sequence = TARGET_POCKET_SEQUENCES.get(target_name, "")
        if not sequence:
            return None

        # 生成蛋白结构
        residues = self.generate_protein_coords(sequence, target_name, start_res=START_RES)
        
        # 结合口袋中心：使用 CC'-C loop 区域（CDR-like）
        # 在拓扑中 C'→C loop 对应残基范围 ≈ 60-72
        n_res = len(residues)
        pocket_start = max(0, n_res - 25)
        pocket_end = min(n_res - 1, n_res - 10)
        
        pocket_x = 0
        pocket_y = 0
        pocket_z = 0
        count = 0
        for r in residues[pocket_start:pocket_end]:
            for atom in r["atoms"]:
                if atom["name"] == "CA":
                    pocket_x += atom["x"]
                    pocket_y += atom["y"]
                    pocket_z += atom["z"]
                    count += 1
        
        if count == 0:
            pocket_center = (25.0, -11.0, 18.0)
        else:
            pocket_center = (pocket_x/count, pocket_y/count, pocket_z/count)

        # 为每个药物模板生成配体
        drug_complexes = []
        for mol_id, mol_template in DRUG_MOLECULE_TEMPLATES.items():
            ligand_data = self.generate_ligand_from_template(
                mol_template, pocket_center, scale=1.0, mol_id=mol_id
            )
            pdb_single = self.to_pdb_single_ligand(residues, ligand_data)

            drug_complexes.append({
                "mol_id": mol_id,
                "name": mol_template.get("name", mol_id),
                "mol_type": mol_template.get("mol_type", "Unknown"),
                "smiles_like": mol_template.get("smiles_like", ""),
                "binding_affinity": mol_template.get("binding_affinity", -7.0),
                "mw": mol_template.get("mw", 350),
                "drug_likeness": mol_template.get("drug_likeness", 0.7),
                "predicted_interactions": mol_template.get("predicted_interactions", {}),
                "pdb_string": pdb_single,
                "ligand_atoms": [
                    {"name": a["name"], "element": a["element"],
                     "x": a["x"], "y": a["y"], "z": a["z"]}
                    for a in ligand_data["atoms"]
                ],
                "ligand_bonds": ligand_data["bonds"],
                "n_ligand_atoms": len(ligand_data["atoms"]),
            })

        # 简单分子
        for idx, simple_mol in enumerate(MORE_SIMPLE_LIGANDS):
            mol_id = f"simple_{idx+1:03d}"
            ligand_data = self.generate_simple_ligand(pocket_center, mol_id, idx)
            pdb_single = self.to_pdb_single_ligand(residues, ligand_data)

            drug_complexes.append({
                "mol_id": mol_id,
                "name": simple_mol.get("name", mol_id),
                "mol_type": "Fragment",
                "smiles_like": simple_mol.get("smiles", ""),
                "binding_affinity": round(-6.0 - random.uniform(0, 2.5), 1),
                "mw": round(220 + random.uniform(50, 180), 1),
                "drug_likeness": round(0.5 + random.uniform(0, 0.4), 2),
                "predicted_interactions": {
                    "hydrogen_bonds": random.randint(1, 3),
                    "pi_stacking": random.randint(0, 2),
                    "hydrophobic": random.randint(2, 5)
                },
                "pdb_string": pdb_single,
                "ligand_atoms": [
                    {"name": a["name"], "element": a["element"],
                     "x": a["x"], "y": a["y"], "z": a["z"]}
                    for a in ligand_data["atoms"]
                ],
                "ligand_bonds": ligand_data["bonds"],
                "n_ligand_atoms": len(ligand_data["atoms"]),
            })

        drug_complexes.sort(key=lambda x: x["binding_affinity"])

        # 结合口袋残基
        binding_pocket_residues = []
        offsets = TARGET_BINDING_SITE_OFFSETS.get(target_name, [])
        for r in residues:
            rn = r["res_num"] - START_RES
            if rn in offsets or (pocket_start <= rn <= pocket_end):
                binding_pocket_residues.append({
                    "res_num": r["res_num"], "aa": r["aa_one"], "aa3": r["aa_three"]
                })

        result = {
            "target": target_name,
            "n_protein_residues": len(residues),
            "n_protein_atoms": len(residues) * 5,
            "binding_pocket_center": list(pocket_center),
            "binding_pocket_residues": binding_pocket_residues[:25],
            "n_drugs": len(drug_complexes),
            "drugs": drug_complexes,
        }
        return result

    def run(self):
        targets = ["PD-1", "LAG-3", "TIM-3", "VISTA"]
        all_data = {}

        for target_name in targets:
            print(f"\n生成 {target_name} 复合物结构（IgV β-三明治 + 候选药物）...")
            result = self.generate_complex_with_drugs(target_name)
            if result:
                all_data[target_name] = result

                out_path = self.output_dir / f"{target_name}_3d_complex.json"
                with open(out_path, 'w', encoding='utf-8') as f:
                    json.dump(result, f, indent=2, ensure_ascii=False)
                print(f"  {result['n_protein_residues']} 残基, {result['n_drugs']} 候选药物 → {out_path}")

                for drug in result["drugs"]:
                    safe_name = drug["mol_id"]
                    drug_pdb_path = self.output_dir / f"{target_name}_{safe_name}.pdb"
                    with open(drug_pdb_path, 'w', encoding='utf-8') as f:
                        f.write(drug["pdb_string"])

        # 保存轻量药物列表
        drug_list_data = {}
        for target_name in targets:
            if target_name in all_data:
                drug_list_data[target_name] = [
                    {
                        'mol_id': d['mol_id'],
                        'name': d['name'],
                        'mol_type': d.get('mol_type', ''),
                        'binding_affinity': d.get('binding_affinity', 0),
                        'mw': d.get('mw', 0),
                        'drug_likeness': d.get('drug_likeness', 0),
                        'predicted_interactions': d.get('predicted_interactions', {}),
                        'n_ligand_atoms': d.get('n_ligand_atoms', 0),
                    }
                    for d in all_data[target_name]["drugs"]
                ]

        drug_list_path = self.output_dir / "all_drug_candidates.json"
        with open(drug_list_path, 'w', encoding='utf-8') as f:
            json.dump(drug_list_data, f, indent=2, ensure_ascii=False)

        print(f"\n全部3D结构已生成 → {drug_list_path}")
        return all_data


if __name__ == "__main__":
    generator = Complex3DGenerator()
    generator.run()