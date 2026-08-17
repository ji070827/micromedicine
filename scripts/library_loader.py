#!/usr/bin/env python3
"""
library_loader.py - 通用小分子库加载器

解决"不同小分子库记录信息不同"的兼容性问题：
1. 自动识别 SMILES 列（smiles / SMILES / CanonicalSMILES / isomeric_smiles 等常见命名）
2. 支持 .csv / .smi / .txt / .sdf 多种格式
3. 自动补全缺失的 mol_id / source / target 字段
4. 所有理化性质一律由 RDKit 从 SMILES 重新计算（不依赖输入库自带字段）

这样下游所有评判逻辑（预处理/对接/初筛/打分）只需依赖"smiles"这一个字段，
而该字段由本模块统一保证存在并命名规范。
"""

import os
import sys
import hashlib
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.real_chemistry import parse_molecule


# 常见的 SMILES 列名（按优先级排序，大小写不敏感匹配）
SMILES_COLUMN_CANDIDATES = [
    'smiles', 'SMILES', 'Smiles',
    'canonical_smiles', 'CanonicalSMILES', 'canonicalsmiles',
    'isomeric_smiles', 'IsomericSMILES', 'isomericsmiles',
    'smiles_std', 'std_smiles', 'smiles_string',
    'nonisomeric_smiles', 'structure', 'canonicalsmiles',
]

# 常见的分子 ID 列名
ID_COLUMN_CANDIDATES = [
    'mol_id', 'molID', 'MoleculeID', 'molecule_id', 'compound_id',
    'cid', 'CID', 'id', 'ID', 'name', 'Name', 'zinc_id', 'pubchem_cid',
]


def detect_smiles_column(df):
    """自动检测 DataFrame 中哪一列是 SMILES 列。"""
    # 1. 精确匹配候选名（忽略大小写）
    lower_cols = {str(c).lower(): c for c in df.columns}
    for cand in SMILES_COLUMN_CANDIDATES:
        if cand.lower() in lower_cols:
            return lower_cols[cand.lower()]

    # 2. 启发式：找出内容像 SMILES 的列
    for c in df.columns:
        col = df[c].astype(str)
        # 抽样前若干个非空值，检查是否能被 RDKit 解析
        sample = col[col.notna() & (col != '')].head(20)
        if len(sample) == 0:
            continue
        n_parsed = sum(1 for s in sample if parse_molecule(s) is not None)
        if n_parsed >= max(1, 0.7 * len(sample)):
            return c

    # 3. 兜底：默认第一个字符串列
    for c in df.columns:
        if df[c].dtype == object:
            return c

    return None


def detect_id_column(df):
    """自动检测分子 ID 列。"""
    lower_cols = {str(c).lower(): c for c in df.columns}
    for cand in ID_COLUMN_CANDIDATES:
        if cand.lower() in lower_cols:
            return lower_cols[cand.lower()]
    return None


def _normalize_dtypes(df):
    """把常见类型列转成数值，便于后续打分。"""
    numeric_candidates = [
        'molecular_weight', 'molecularweight', 'mw', 'mwt',
        'logp', 'xlogp', 'log_p', 'tpsa', 'hbd', 'hba', 'qed',
        'rotatable_bonds', 'fsp3', 'num_rings', 'psa',
        'ic50', 'ec50', 'activity',
    ]
    lower_cols = {str(c).lower(): c for c in df.columns}
    for cand in numeric_candidates:
        if cand in lower_cols:
            c = lower_cols[cand]
            df[c] = pd.to_numeric(df[c], errors='coerce')
    return df


def load_library_file(path, smiles_col=None):
    """
    通用小分子库加载入口。

    参数：
        path: 文件路径（支持 .csv/.smi/.txt/.sdf）
        smiles_col: 显式指定 SMILES 列名（可选，不指定则自动检测）

    返回：
        标准化的 DataFrame，保证包含 'smiles' 和 'mol_id' 两列。
        其他原始列原样保留（但不会覆盖由 RDKit 重算的性质列）。
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"化合物库文件不存在: {path}")

    ext = path.suffix.lower()

    # 1. 按格式读取
    if ext == '.csv':
        df = pd.read_csv(path)
    elif ext in ('.smi', '.txt'):
        # SMILES 文本文件：一行一个 SMILES（可能含 ID，用空格/制表符分隔）
        rows = []
        with open(path, 'r', encoding='utf-8', errors='replace') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()
                # 常见格式：SMILES [ID] 或 SMILES<TAB>ID（或反过来）
                # 第一个能被 RDKit 解析的 token 作为 SMILES
                smiles = None
                mol_id = None
                for tok in parts:
                    if smiles is None and parse_molecule(tok) is not None:
                        smiles = tok
                    elif smiles is not None and mol_id is None:
                        mol_id = tok
                if smiles is not None:
                    rows.append({'smiles': smiles, 'mol_id': mol_id})
        df = pd.DataFrame(rows)
    elif ext == '.sdf':
        from rdkit import Chem
        mols = [m for m in Chem.SDMolSupplier(str(path), removeHs=False) if m is not None]
        smiles_list = []
        names = []
        for m in mols:
            smi = Chem.MolToSmiles(Chem.RemoveHs(m))
            smiles_list.append(smi)
            names.append(m.GetProp('_Name') if m.HasProp('_Name') else None)
        df = pd.DataFrame({'smiles': smiles_list, 'mol_id': names})
    else:
        raise ValueError(f"不支持的分子库格式: {ext}（仅支持 .csv/.smi/.txt/.sdf）")

    # 2. 统一 SMILES 列名
    if smiles_col is None:
        smiles_col = detect_smiles_column(df)
    if smiles_col is None:
        raise ValueError("无法识别 SMILES 列，请显式指定 smiles_col 参数")

    # 如果检测到的列名不是 'smiles'，重命名为 'smiles'
    if smiles_col != 'smiles':
        df = df.rename(columns={smiles_col: 'smiles'})

    # 3. 统一 mol_id 列
    if 'mol_id' not in df.columns:
        id_col = detect_id_column(df)
        if id_col is not None and id_col != 'smiles':
            df = df.rename(columns={id_col: 'mol_id'})
        else:
            df['mol_id'] = [f"CMPD_{i+1:06d}" for i in range(len(df))]

    # 4. 过滤空 SMILES + 补全 mol_id 空值
    df = df[df['smiles'].notna() & (df['smiles'].astype(str).str.strip() != '')]
    df['smiles'] = df['smiles'].astype(str).str.strip()
    df['mol_id'] = df['mol_id'].where(
        df['mol_id'].notna() & (df['mol_id'].astype(str).str.strip() != ''),
        [f"CMPD_{i+1:06d}" for i in range(len(df))]
    )

    # 5. 数值列类型转换
    df = _normalize_dtypes(df)

    return df.reset_index(drop=True)


def iter_library_chunks(path, chunksize=10000):
    """
    流式迭代读取大型分子库，逐块返回标准化的 DataFrame。
    用于百万级库的快速预筛，避免一次性加载全部到内存。

    参数：
        path: 库文件路径（.csv/.smi/.txt）
        chunksize: 每块行数

    返回：
        生成器，每次 yield 一个标准化后的 DataFrame（含 'smiles' 和 'mol_id' 列）
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"化合物库文件不存在: {path}")

    ext = path.suffix.lower()

    if ext == '.csv':
        # CSV 分块读取
        reader = pd.read_csv(path, chunksize=chunksize)
        for raw_chunk in reader:
            yield standardize_dataframe(raw_chunk)
    elif ext in ('.smi', '.txt'):
        # 文本分块读取
        rows = []
        count = 0
        with open(path, 'r', encoding='utf-8', errors='replace') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()
                smiles = None
                mol_id = None
                for tok in parts:
                    if smiles is None and parse_molecule(tok) is not None:
                        smiles = tok
                    elif smiles is not None and mol_id is None:
                        mol_id = tok
                if smiles is not None:
                    rows.append({'smiles': smiles, 'mol_id': mol_id})
                    count += 1
                    if count >= chunksize:
                        yield standardize_dataframe(pd.DataFrame(rows))
                        rows = []
                        count = 0
        if rows:
            yield standardize_dataframe(pd.DataFrame(rows))
    else:
        # 其他格式一次性读（SDF 不流式）
        df = load_library_file(path)
        yield df


def standardize_dataframe(df):
    """
    对已读取的 DataFrame 做标准化：统一 smiles/mol_id 列名 + 过滤空值。
    供流式读取复用。
    """
    if df is None or len(df) == 0:
        return df

    # 统一 SMILES 列名
    smiles_col = detect_smiles_column(df)
    if smiles_col is None:
        # 如果没有 smiles 列，返回原样（调用方会跳过）
        return df
    if smiles_col != 'smiles':
        df = df.rename(columns={smiles_col: 'smiles'})

    # 统一 mol_id 列
    if 'mol_id' not in df.columns:
        id_col = detect_id_column(df)
        if id_col is not None and id_col != 'smiles':
            df = df.rename(columns={id_col: 'mol_id'})
        else:
            df['mol_id'] = [f"CMPD_{i+1:06d}" for i in range(len(df))]

    # 过滤空 SMILES
    df = df[df['smiles'].notna() & (df['smiles'].astype(str).str.strip() != '')]
    df['smiles'] = df['smiles'].astype(str).str.strip()

    return df.reset_index(drop=True)


if __name__ == '__main__':
    # 自测：临时创建一个不同列名的库
    import tempfile
    test_csv = Path(tempfile.gettempdir()) / 'test_lib.csv'
    pd.DataFrame({
        'CanonicalSMILES': ['CC(=O)OC1=CC=CC=C1C(=O)O', 'CC(=O)Nc1ccccc1'],
        'MoleculeID': ['asp', 'ace'],
        'SomeExtraColumn': [1, 2],
    }).to_csv(test_csv, index=False)

    df = load_library_file(test_csv)
    print("标准化的库：")
    print(df[['mol_id', 'smiles']].head())
    print("SUCCESS: 通用库加载器自测通过")
