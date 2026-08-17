#!/usr/bin/env python3
"""
download_real_data.py - 下载真实蛋白结构 + 真实小分子库 (服务器运行时)

在具有正常国际网络访问的 Linux 服务器上运行，获取：
1. 四个免疫检查点靶点的真实 PDB 晶体结构 (RCSB PDB)
2. 已知的免疫检查点小分子抑制剂 (PubChem)

注意：本机 Windows 环境访问 RCSB/PubChem 受限，
此脚本设计为在服务器上执行。
"""

import os
import sys
import json
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

BASE_DIR = Path(__file__).parent.parent

# 四个靶点的真实 PDB ID
TARGET_PDB_IDS = {
    "PD-1": "4ZQK",    # PD-1/PD-L1 复合物晶体结构
    "LAG-3": "7TZH",   # LAG-3
    "TIM-3": "5F71",   # TIM-3
    "VISTA": "6OIL",   # VISTA
}

# 已知免疫检查点小分子抑制剂（PubChem CID，用于构建真实阳性分子库）
KNOWN_INHIBITOR_CIDS = {
    "PD-1": [23629198, 25023587, 447290, 49867904],       # BMS 系列等
    "LAG-3": [137347988, 72721919],
    "TIM-3": [444899, 16074],
    "VISTA": [54694254],
}


def download_with_retry(url, output_path, retries=3, timeout=60):
    """带重试的下载"""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                content = resp.read()
            with open(output_path, "wb") as f:
                f.write(content)
            return True
        except Exception as e:
            print(f"    尝试 {attempt+1}/{retries} 失败: {e}")
            if attempt < retries - 1:
                time.sleep(5)
    return False


def download_real_pdb_structures():
    """下载真实 PDB 结构"""
    print("\n=== 下载真实蛋白质结构 (RCSB PDB) ===")
    out_dir = BASE_DIR / "data" / "targets" / "real_structures"
    out_dir.mkdir(parents=True, exist_ok=True)

    results = {}
    for target, pdb_id in TARGET_PDB_IDS.items():
        print(f"\n{target}: PDB {pdb_id}")
        pdb_path = out_dir / f"{target}_{pdb_id}.pdb"
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

        if download_with_retry(url, pdb_path):
            size_kb = pdb_path.stat().st_size / 1024
            print(f"  ✅ {pdb_id}.pdb ({size_kb:.0f} KB)")
            results[target] = {"pdb_id": pdb_id, "file": str(pdb_path)}
        else:
            print(f"  ❌ {pdb_id} 下载失败")
            results[target] = {"pdb_id": pdb_id, "file": None}

    # 保存摘要
    with open(out_dir / "download_summary.json", "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    return results


def fetch_pubchem_smiles(cid):
    """从 PubChem 获取化合物的 SMILES"""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES,MolecularWeight,XLogP/JSON"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read())
        props = data["PropertyTable"]["Properties"][0]
        return {
            "cid": cid,
            "smiles": props.get("CanonicalSMILES", ""),
            "mw": props.get("MolecularWeight"),
            "logp": props.get("XLogP"),
        }
    except Exception as e:
        print(f"    CID {cid} 获取失败: {e}")
        return None


def build_real_compound_library():
    """构建已知抑制剂的小分子库"""
    print("\n=== 构建真实小分子抑制剂库 (PubChem) ===")
    all_compounds = []

    for target, cids in KNOWN_INHIBITOR_CIDS.items():
        print(f"\n{target}: {len(cids)} 个已知抑制剂")
        for cid in cids:
            info = fetch_pubchem_smiles(cid)
            if info and info["smiles"]:
                info["target"] = target
                all_compounds.append(info)
                print(f"  ✅ CID {cid}: {info['smiles'][:50]}...")
            time.sleep(1)  # 避免请求过快

    # 保存
    if all_compounds:
        import pandas as pd
        df = pd.DataFrame(all_compounds)
        out_path = BASE_DIR / "data" / "library" / "real_inhibitors.csv"
        df.to_csv(out_path, index=False)
        print(f"\n✅ 真实抑制剂库已保存: {out_path} ({len(df)} 个化合物)")

    return all_compounds


if __name__ == "__main__":
    print("=" * 60)
    print("真实数据下载脚本 (需服务器环境 + 正常网络)")
    print("=" * 60)
    download_real_pdb_structures()
    build_real_compound_library()
    print("\n✅ 全部真实数据下载完成")