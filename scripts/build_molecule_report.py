#!/usr/bin/env python3
"""
build_molecule_report.py - 生成每个分子的"详情 + 成药性理由"报告

读取所有筛选结果，为每个分子生成：
1. 结构信息（SMILES + RDKit 2D SVG）
2. 成药性判断（是否适合做药）
3. 打分理由（七维得分拆解 + Lipinski 逐项 + ADME 解释）
4. 与靶点的对接/结合信息

输出：results/molecule_report.json（前端实时读取）
"""

import os
import sys
import json
import hashlib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd
import yaml

from scripts.real_chemistry import (
    parse_molecule, compute_properties_from_mol, check_pains_brenk,
)

import warnings
warnings.filterwarnings('ignore')


def _explain_lipinski(props):
    """逐项解释 Lipinski 五规则"""
    items = []
    checks = [
        ("分子量 ≤ 500", props['molecular_weight'] <= 500, f"{props['molecular_weight']:.1f}"),
        ("logP ≤ 5", props['logP'] <= 5, f"{props['logP']:.2f}"),
        ("氢键供体 ≤ 5", props['HBD'] <= 5, f"{props['HBD']}"),
        ("氢键受体 ≤ 10", props['HBA'] <= 10, f"{props['HBA']}"),
        ("可旋转键 ≤ 10", props['rotatable_bonds'] <= 10, f"{props['rotatable_bonds']}"),
    ]
    for rule, ok, value in checks:
        items.append({"rule": rule, "pass": ok, "value": value})
    return items


def _explain_adme(props):
    """简化的 ADME 解释（基于真实性质）"""
    mw = props['molecular_weight']
    logp = props['logP']
    tpsa = props['TPSA']

    reasons = []
    # 吸收（Absorption）
    if tpsa < 140 and logp < 5:
        reasons.append("口服吸收倾向良好（TPSA<140 且 logP<5）")
    elif tpsa >= 140:
        reasons.append("极性表面积较大（≥140），口服吸收可能受限")
    else:
        reasons.append("logP 偏高，可能影响水溶性")

    # 血脑屏障（BB 穿透）
    if logp > 3 and tpsa < 90 and mw < 450:
        reasons.append("可能穿透血脑屏障（需关注中枢副作用）")

    # 溶解度
    if logp > 5:
        reasons.append("logP>5，水溶性可能差")
    elif logp < -1:
        reasons.append("logP<-1，脂溶性可能不足")

    return reasons


def _explain_pains(is_pains, pains_desc):
    if is_pains:
        return [f"⚠ 命中 PAINS 警示子结构（可能为假阳性靶点干扰物）"]
    return ["✅ 未发现 PAINS/Brenk 警示子结构"]


def _drug_likeness_verdict(props, is_pains):
    """综合判定是否适合做药"""
    violations = props['lipinski_violations']
    qed = props['QED']
    sa = props['sa_score']

    reasons = []
    if is_pains:
        verdict = "不适合（毒性警示）"
        reasons.append("PAINS 警示：可能是非特异性干扰物")
    elif violations == 0 and qed >= 0.6:
        verdict = "适合（成药性良好）"
        reasons.append(f"Lipinski 全部通过，QED={qed:.2f} 成药性评分高")
    elif violations <= 1:
        verdict = "基本适合（轻微瑕疵）"
        reasons.append(f"Lipinski 违反 {violations} 项，QED={qed:.2f}")
    else:
        verdict = "不适合（成药性差）"
        reasons.append(f"Lipinski 违反 {violations} 项，成药性差")

    reasons.append(f"合成可及性 SA={sa:.1f}")
    return verdict, reasons


def build_report():
    base_dir = Path(__file__).parent.parent
    config = yaml.safe_load(open(base_dir / "config" / "config.yaml", encoding='utf-8'))
    active_lib = config.get('data', {}).get('active_library', 'fda_approved_drugs.csv')

    # 加载化合物库
    lib_path = base_dir / "data" / "library" / active_lib
    if not lib_path.exists():
        lib_path = base_dir / "data" / "library" / "compounds_standardized.csv"
    if not lib_path.exists():
        print("未找到化合物库，无法生成报告")
        return None

    from scripts.library_loader import load_library_file
    df_lib = load_library_file(lib_path)

    # 加载终选报告（如果有）
    final_path = base_dir / "results" / "final_report" / "all_top_candidates.csv"
    df_final = None
    if final_path.exists():
        df_final = pd.read_csv(final_path)

    # 对每个分子生成详情
    molecules = []
    for _, row in df_lib.iterrows():
        smiles = row.get('smiles', '')
        mol = parse_molecule(smiles)
        if mol is None:
            continue

        props = compute_properties_from_mol(mol)
        is_pains, pains_desc = check_pains_brenk(mol)
        verdict, verdict_reasons = _drug_likeness_verdict(props, is_pains)

        entry = {
            'mol_id': row.get('mol_id', f"CMPD_{len(molecules)+1:06d}"),
            'smiles': smiles,
            'molecular_weight': props['molecular_weight'],
            'logP': props['logP'],
            'TPSA': props['TPSA'],
            'HBD': props['HBD'],
            'HBA': props['HBA'],
            'rotatable_bonds': props['rotatable_bonds'],
            'QED': props['QED'],
            'sa_score': props['sa_score'],
            'drug_likeness_verdict': verdict,
            'verdict_reasons': verdict_reasons,
            'lipinski_items': _explain_lipinski(props),
            'adme_reasons': _explain_adme(props),
            'pains_reasons': _explain_pains(is_pains, pains_desc),
        }

        # 如果有终选结果，附加上七维得分
        if df_final is not None:
            match = df_final[df_final['mol_id'] == entry['mol_id']]
            if len(match):
                r = match.iloc[0]
                entry['final_score'] = float(r.get('final_score', 0))
                entry['priority_class'] = r.get('priority_class', '')
                # 七维拆解（如果存在）
                if 'score_breakdown' in r and isinstance(r['score_breakdown'], (dict, str)):
                    try:
                        entry['score_breakdown'] = json.loads(r['score_breakdown']) if isinstance(r['score_breakdown'], str) else r['score_breakdown']
                    except Exception:
                        entry['score_breakdown'] = {}

        molecules.append(entry)

    # 保存
    report = {
        'active_library': active_lib,
        'total_molecules': len(molecules),
        'molecules': molecules,
    }
    out_path = base_dir / "results" / "molecule_report.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w', encoding='utf-8') as f:
        json.dump(report, f, indent=2, ensure_ascii=False)

    print(f"✅ 分子报告已生成：{out_path}（{len(molecules)} 个分子）")
    return report


if __name__ == "__main__":
    build_report()