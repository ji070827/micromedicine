#!/usr/bin/env python3
"""
download_million_library.py - 下载百万级小分子库（用于海选）

真实数据源：
1. ZINC20/22（可购化合物库，~2B 分子，支持按子集下载）
2. PubChem（~110M 化合物，可通过 FTP bulk 下载）

用法（在服务器上、有正常网络时运行）：
  python scripts/download_million_library.py --source zinc --size 1000000
"""

import os
import sys
import argparse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

BASE_DIR = Path(__file__).parent.parent


def download_zinc_subset(output_path, n_target=1000000):
    """下载 ZINC 库子集（说明 + 指引）"""
    print("=" * 60)
    print("下载 ZINC 百万级小分子库")
    print("=" * 60)
    print(f"目标: 约 {n_target} 个小分子")
    print("数据源: ZINC (可购化合物库)")
    print()
    print("ZINC 库非常大（~2B），需按 tranche 分片下载：")
    print("  1. 访问 https://zinc.docking.org/ 选择可购性子集")
    print("  2. 下载 .smi.gz 分片文件到 data/library/")
    print("  3. 解压后运行: python scripts/rapid_prefilter.py")
    print()
    print("常用公开百万级数据集（可直接下载 SMILES）：")
    print("  - ZINC Gold 子集（数百万分子）")
    print("  - Enamine / ChemBridge 可购库")
    print("  - PubChem bulk SDF 分片（每个 ~2GB，10万+ 分子）")
    return None


def main():
    parser = argparse.ArgumentParser(description="下载百万级小分子库")
    parser.add_argument("--source", default="zinc", choices=["zinc", "pubchem"], help="数据源")
    parser.add_argument("--size", type=int, default=1000000, help="目标分子数")
    parser.add_argument("--output", default=None, help="输出文件路径")
    args = parser.parse_args()

    if args.output is None:
        args.output = BASE_DIR / "data" / "library" / f"{args.source}_million.smi"

    print("=" * 60)
    print("百万级小分子库下载")
    print("=" * 60)
    print(f"数据源: {args.source}")
    print(f"目标规模: {args.size:,} 分子")
    print(f"输出: {args.output}")
    print()

    if args.source == "zinc":
        download_zinc_subset(args.output, args.size)
    elif args.source == "pubchem":
        print("PubChem bulk 下载方案：")
        print("  1. 访问 https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/")
        print("  2. 下载 SDF 分片文件（每个 ~2GB，含 10 万+ 分子）")
        print("  3. 用 RDKit 转成 SMILES 后放入 data/library/")
        print("  PubChem 完整库 ~110M 分子，按需下载分片")

    print()
    print("下载完成后，海选流程：")
    print("  1. python scripts/rapid_prefilter.py     # 快速预筛（百万→几千）")
    print("  2. python scripts/diffdock_batch_run.py  # DiffDock 对接（几千→几百）")
    print("  3. AlphaFold3 → 相互作用 → 终选")


if __name__ == "__main__":
    main()