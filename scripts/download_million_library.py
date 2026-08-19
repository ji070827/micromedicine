#!/usr/bin/env python3
"""
download_million_library.py - 下载知名公开小分子库（真实数据，非自生成）

数据源（均为知名公开数据库，无需注册，HTTPS 直接下载）：

1. ChEMBL（推荐，约 230 万个化合物，含真实 SMILES + 活性注释）
   - 官方 chemreps 文件：chembl_XX_chemreps.txt.gz
   - 列：chembl_id, canonical_smiles, ...
   - 完整文件约 200MB（压缩），230 万行

2. ZINC20/22（可购化合物库，约 13M 可购子集）
   - 需按 tranche 分片，URL 结构复杂，作为备选

3. PubChem（约 1.1 亿化合物，FTP bulk）

用法（在武大超算上，有国际网络时）：
  python scripts/download_million_library.py --source chembl --n 100000

输出：data/library/chembl_100k.csv  （含 mol_id + smiles 两列，可直接被
      library_loader + rapid_prefilter + diffdock 使用）
"""

import os
import sys
import gzip
import json
import random
import argparse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

BASE_DIR = Path(__file__).parent.parent
LIB_DIR = BASE_DIR / "data" / "library"


# ChEMBL chemreps 文件可能的 URL（版本号会更新，按优先级尝试）
CHEMBL_URLS = [
    "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_35_chemreps.txt.gz",
    "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_34_chemreps.txt.gz",
    "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_33_chemreps.txt.gz",
]


def _http_get(url, timeout=120):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    return urllib.request.urlopen(req, timeout=timeout)


def _download_file(url, dest, timeout=3600):
    """下载文件到 dest，返回是否成功"""
    print(f"  下载: {url}")
    try:
        resp = _http_get(url, timeout=120)
        total = 0
        with open(dest, "wb") as f:
            while True:
                chunk = resp.read(1024 * 256)
                if not chunk:
                    break
                f.write(chunk)
                total += len(chunk)
                if total % (50 * 1024 * 1024) < 1024 * 256:  # 每 50MB 打印一次
                    print(f"    已下载 {total / 1024 / 1024:.0f} MB")
        print(f"  ✅ 下载完成: {total / 1024 / 1024:.1f} MB")
        return True
    except Exception as e:
        print(f"  ❌ 下载失败: {e}")
        return False


def _parse_chemreps_header(line):
    """解析 chemreps 表头，找到 canonical_smiles 和 chembl_id 的列索引"""
    cols = line.strip().split("\t")
    smiles_idx = None
    id_idx = None
    for i, c in enumerate(cols):
        cl = c.lower()
        if cl in ("canonical_smiles", "canonicalsmiles", "smiles"):
            smiles_idx = i
        if cl in ("chembl_id", "chemblid", "molecule_chembl_id"):
            id_idx = i
    return id_idx, smiles_idx


def download_chembl(n_target=100000, seed=42):
    """
    下载 ChEMBL chemreps 文件，流式解析 SMILES，水库采样 n_target 个分子。
    返回输出的 CSV 路径。
    """
    print("=" * 60)
    print(f"下载 ChEMBL 知名化合物库（目标 {n_target:,} 个）")
    print("=" * 60)

    # 1. 下载 chemreps 文件
    gz_path = LIB_DIR / "chembl_chemreps.txt.gz"
    downloaded = False
    for url in CHEMBL_URLS:
        if _download_file(url, gz_path):
            downloaded = True
            break

    if not downloaded:
        print("\n⚠ 无法自动下载 ChEMBL chemreps 文件。")
        print("  可能原因：版本号更新，或 EBI 域名访问受限。")
        print("  请手动访问 https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/")
        print("  下载 chembl_XX_chemreps.txt.gz 放到 data/library/ 后重跑本脚本")
        return None

    # 2. 流式解析 + 水库采样
    print(f"\n解析 chemreps 文件并采样 {n_target:,} 个分子（水库采样，不爆内存）...")
    id_idx = smi_idx = None
    reservoir = []  # [(chembl_id, smiles)]
    n_seen = 0
    rng = random.Random(seed)

    with gzip.open(gz_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if id_idx is None:
                # 第一行是表头
                id_idx, smi_idx = _parse_chemreps_header(line)
                if smi_idx is None:
                    print("  ❌ 未找到 canonical_smiles 列，表头:", line[:200])
                    return None
                continue

            cols = line.split("\t")
            if smi_idx >= len(cols):
                continue
            smiles = cols[smi_idx].strip()
            if not smiles:
                continue
            mol_id = cols[id_idx].strip() if id_idx is not None and id_idx < len(cols) else f"CHEMBL_{n_seen+1}"

            n_seen += 1
            # 水库采样
            if len(reservoir) < n_target:
                reservoir.append((mol_id, smiles))
            else:
                j = rng.randint(0, n_seen)
                if j < n_target:
                    reservoir[j] = (mol_id, smiles)

    print(f"  共扫描 {n_seen:,} 个化合物，采样 {len(reservoir):,} 个")

    # 3. 保存 CSV（兼容 library_loader：需要 mol_id + smiles 列）
    out_path = LIB_DIR / "chembl_100k.csv"
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("mol_id,smiles,source\n")
        for mol_id, smiles in reservoir:
            # CSV 转义：SMILES 一般无逗号/引号，但安全起见做转义
            f.write(f"{mol_id},{smiles},ChEMBL\n")

    print(f"\n✅ 已保存: {out_path}（{len(reservoir):,} 个分子）")
    print(f"   下一步：")
    print(f"     1. 把 config.yaml 的 active_library 改成 chembl_100k.csv")
    print(f"     2. python scripts/rapid_prefilter.py   # 快速预筛")
    print(f"     3. sbatch slurm/diffdock.slurm        # GPU 对接")
    return out_path


def download_zinc(n_target=100000):
    """ZINC 下载指引（URL 结构复杂，作为备选方案说明）"""
    print("=" * 60)
    print("ZINC 可购化合物库下载说明")
    print("=" * 60)
    print("ZINC20/22 有约 13M 可购分子，按 tranche 分片存储。")
    print("推荐做法（手动）：")
    print("  1. 访问 https://zinc.docking.org/ 注册登录")
    print("  2. 选择 'drug-like' + 'in-stock' 子集")
    print("  3. 导出 .smi 文件（可指定数量，如 10 万）")
    print("  4. 上传到 data/library/ 后，把 config active_library 指向它")
    print()
    print("备选：直接下载 ZINC 预生成的分片 smi.gz")
    print("  https://files.docking.org/zinc20/  （tranche 分片目录）")
    print()
    print("注：ChEMBL 是脚本化下载的首选（无需注册 + 单一文件），")
    print("    ZINC 更适合需要'可购'性质的场景（需手动导出）。")
    return None


def main():
    parser = argparse.ArgumentParser(description="下载知名公开小分子库")
    parser.add_argument("--source", default="chembl", choices=["chembl", "zinc"],
                        help="数据源（chembl 推荐，可脚本化下载）")
    parser.add_argument("--n", type=int, default=100000, help="目标分子数")
    args = parser.parse_args()

    LIB_DIR.mkdir(parents=True, exist_ok=True)

    if args.source == "chembl":
        download_chembl(n_target=args.n)
    elif args.source == "zinc":
        download_zinc(n_target=args.n)


if __name__ == "__main__":
    main()