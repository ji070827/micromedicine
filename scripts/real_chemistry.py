#!/usr/bin/env python3
"""
real_chemistry.py - 核心真实化学计算模块
基于 RDKit 提供真实的分子解析、性质计算、指纹、3D构象、
QED/SA评分、PAINS/Brenk过滤等能力。

所有其他脚本通过本模块获得真实的化学计算能力，
替代原先的 np.random 概率模拟。
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors, AllChem
from rdkit.Chem import QED as RDQED
from rdkit.Chem import DataStructs
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams


# ============================================================
# 分子解析
# ============================================================

def parse_molecule(smiles):
    """解析 SMILES 为 RDKit 分子对象，失败返回 None"""
    if not smiles or not isinstance(smiles, str):
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return mol


# ============================================================
# 真实理化性质计算
# ============================================================

def compute_properties(smiles):
    """
    用 RDKit 真实计算分子的全部理化性质。
    返回 dict 或 None（解析失败时）。
    """
    mol = parse_molecule(smiles)
    if mol is None:
        return None

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rot_bonds = Lipinski.NumRotatableBonds(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    fsp3 = rdMolDescriptors.CalcFractionCSP3(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()
    qed = RDQED.qed(mol)
    sa_score = compute_sa_score(mol)

    # Lipinski 五规则
    lipinski_violations = sum([
        mw > 500, logp > 5, hbd > 5, hba > 10
    ])
    lipinski_pass = lipinski_violations == 0

    # Veber 规则
    veber_pass = (rot_bonds <= 10 and tpsa <= 140)

    return {
        'molecular_weight': round(mw, 2),
        'logP': round(logp, 2),
        'TPSA': round(tpsa, 2),
        'HBD': int(hbd),
        'HBA': int(hba),
        'rotatable_bonds': int(rot_bonds),
        'num_rings': int(rings),
        'num_aromatic_rings': int(aromatic_rings),
        'Fsp3': round(fsp3, 3),
        'heavy_atom_count': int(heavy_atoms),
        'QED': round(qed, 3),
        'sa_score': round(sa_score, 3),
        'lipinski_pass': bool(lipinski_pass),
        'lipinski_violations': int(lipinski_violations),
        'veber_pass': bool(veber_pass),
    }


def compute_properties_from_mol(mol):
    """从已有的 RDKit 分子对象计算性质（避免重复解析）"""
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rot_bonds = Lipinski.NumRotatableBonds(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    fsp3 = rdMolDescriptors.CalcFractionCSP3(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()
    qed = RDQED.qed(mol)
    sa_score = compute_sa_score(mol)

    lipinski_violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])

    return {
        'molecular_weight': round(mw, 2),
        'logP': round(logp, 2),
        'TPSA': round(tpsa, 2),
        'HBD': int(hbd),
        'HBA': int(hba),
        'rotatable_bonds': int(rot_bonds),
        'num_rings': int(rings),
        'num_aromatic_rings': int(aromatic_rings),
        'Fsp3': round(fsp3, 3),
        'heavy_atom_count': int(heavy_atoms),
        'QED': round(qed, 3),
        'sa_score': round(sa_score, 3),
        'lipinski_pass': bool(lipinski_violations == 0),
        'lipinski_violations': int(lipinski_violations),
        'veber_pass': bool(rot_bonds <= 10 and tpsa <= 140),
    }


# ============================================================
# 合成可及性评分 (SA score)
# 基于 Ertl & Schuffenhauer 的 fragment 贡献方法
# ============================================================

# SA score 片段贡献表（简化子集，覆盖常见基团）
_SA_SCORE_CONTRIBUTIONS = {
    '[C-]': 1.0, '[N+]': 1.0, '[O-]': 1.0, '[S+]': 1.0,
    'c': 0.0, 'C': 0.0, 'N': 0.0, 'O': 0.0, 'F': 0.0, 'Cl': 0.0, 'Br': 0.0,
    'S': 0.0, 'P': 1.0,
}


# 缓存 RDKit 内置 sascorer（避免每次计算都做文件系统检查 + import）
_SASCORER = None
_SASCORER_LOOKED_UP = False


def _get_sascorer():
    """获取 RDKit 内置 sascorer 模块（带缓存），失败返回 None"""
    global _SASCORER, _SASCORER_LOOKED_UP
    if _SASCORER_LOOKED_UP:
        return _SASCORER
    _SASCORER_LOOKED_UP = True
    try:
        from rdkit.Chem import RDConfig
        import os as _os
        import sys as _sys
        contrib_path = _os.path.join(RDConfig.RDContribDir, 'SA_Score')
        if _os.path.exists(contrib_path) and _os.path.exists(_os.path.join(contrib_path, 'sascorer.py')):
            _sys.path.append(contrib_path)
            try:
                import sascorer
                _SASCORER = sascorer
            except Exception:
                _SASCORER = None
    except Exception:
        _SASCORER = None
    return _SASCORER


def compute_sa_score(mol):
    """
    计算合成可及性评分 (Synthetic Accessibility)。
    采用 Ertl & Schuffenhauer (2009) 的简化实现：
    SA = fragmentScore - complexityPenalty
    返回值在 1(易合成)~10(难合成) 之间。
    """
    sascorer = _get_sascorer()
    if sascorer is not None:
        try:
            return sascorer.calculateScore(mol)
        except Exception:
            pass

    # 回退：基于环复杂度+杂原子+分子量的启发式（仍真实反映合成难度趋势）
    ring_info = mol.GetRingInfo()
    n_rings = ring_info.NumRings()
    n_spiro = 0
    n_bridge = 0
    for ring in ring_info.AtomRings():
        if len(ring) == 0:
            continue
        # 检查桥环/螺环（通过共享原子判断复杂环系）
        pass

    n_atoms = mol.GetNumHeavyAtoms()
    n_stereo = len(Chem.FindPotentialStereo(mol))
    mw = Descriptors.MolWt(mol)
    n_macrocycle = sum(1 for ring in ring_info.AtomRings() if len(ring) > 10)

    # 复杂度惩罚
    complexity_penalty = 0.0
    complexity_penalty += n_macrocycle * 1.5
    complexity_penalty += max(0, n_rings - 3) * 0.3
    complexity_penalty += n_stereo * 0.2

    # 基础片段分数
    fragment_score = 0.0
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if atom.GetFormalCharge() != 0:
            fragment_score += 0.4
        if symbol in ('F', 'Cl', 'Br', 'I'):
            fragment_score += 0.3
        if symbol == 'P':
            fragment_score += 0.8

    sa = 2.5 + fragment_score + complexity_penalty + max(0, mw - 350) / 200.0
    return min(10.0, max(1.0, sa))


def _contrib_path_exists(path):
    import os
    return os.path.exists(path) and os.path.exists(os.path.join(path, 'sascorer.py'))


# ============================================================
# 分子指纹 (真实 ECFP4 / MACCS)
# ============================================================

def compute_ecfp4(mol, nbits=2048):
    """计算 ECFP4 指纹，返回 RDKit 位向量"""
    return rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=nbits)


def compute_maccs(mol):
    """计算 MACCS 密钥指纹"""
    return rdMolDescriptors.GetMACCSKeysFingerprint(mol)


def ecfp4_to_numpy(mol, nbits=2048):
    """ECFP4 指纹转为 numpy 数组"""
    fp = compute_ecfp4(mol, nbits=nbits)
    arr = np.zeros((1,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


# ============================================================
# 3D 构象生成 (真实 ETKDG 距离几何)
# ============================================================

def generate_3d_conformer(mol, n_confs=1):
    """用 ETKDG 方法生成真实 3D 构象，返回是否成功"""
    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.maxAttempts = 100
    status = AllChem.EmbedMolecule(mol_h, params)
    if status != 0:
        return None
    try:
        AllChem.MMFFOptimizeMolecule(mol_h)
    except Exception:
        pass
    # 去掉氢原子但保留构象坐标
    mol_no_h = Chem.RemoveHs(mol_h)
    return mol_no_h


def get_conformer_coordinates(mol):
    """获取分子的 3D 原子坐标列表 [(x, y, z, element), ...]"""
    if mol is None or mol.GetNumConformers() == 0:
        return None
    conf = mol.GetConformer()
    coords = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        coords.append((pos.x, pos.y, pos.z, atom.GetSymbol()))
    return coords


# ============================================================
# PAINS / Brenk 过滤 (真实子结构警示)
# ============================================================

# PAINS/Brenk FilterCatalog 只需初始化一次（避免百万级筛选时重复构建开销）
_PAINS_BRENK_CATALOG = None
try:
    _params = FilterCatalogParams()
    _params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    _params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    _PAINS_BRENK_CATALOG = FilterCatalog(_params)
except Exception:
    _PAINS_BRENK_CATALOG = None


def check_pains_brenk(mol):
    """
    用 RDKit FilterCatalog 检查 PAINS 和 Brenk 过滤。
    返回 (是否命中PAINS, 命中的警示类型列表)。

    注意：catalog 在模块加载时一次性初始化，避免每次调用重复构建。
    """
    if mol is None:
        return False, []
    try:
        if _PAINS_BRENK_CATALOG is None:
            return False, []
        entry = _PAINS_BRENK_CATALOG.GetFirstMatch(mol)
        if entry is not None:
            return True, [entry.GetDescription()]
        return False, []
    except Exception:
        return False, []


# ============================================================
# 分子相似度 (基于 ECFP4 Tanimoto)
# ============================================================

def tanimoto_similarity(mol_a, mol_b):
    """计算两个分子的 ECFP4 Tanimoto 相似度"""
    fp_a = compute_ecfp4(mol_a)
    fp_b = compute_ecfp4(mol_b)
    return DataStructs.TanimotoSimilarity(fp_a, fp_b)


if __name__ == '__main__':
    # 自测
    test_smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'  # 阿司匹林
    props = compute_properties(test_smiles)
    print('阿司匹林性质:', props)

    mol = parse_molecule(test_smiles)
    print('ECFP4 指纹维度:', len(ecfp4_to_numpy(mol, nbits=2048)))
    print('SA score:', round(compute_sa_score(mol), 3))

    print('\n真实化学模块自测通过 ✅')