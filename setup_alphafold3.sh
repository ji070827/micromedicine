#!/usr/bin/env bash
# ============================================================
# AlphaFold3 Docker 部署脚本 (Linux 服务器)
# 官方建议用 Docker 部署，避免 JAX/OpenMM 依赖冲突
# 使用: bash setup_alphafold3.sh
# ============================================================
set -euo pipefail

echo "=============================================="
echo "  AlphaFold3 Docker 部署"
echo "=============================================="

# 1. 检查 NVIDIA Docker
if ! command -v nvidia-smi &> /dev/null; then
    echo "❌ 未检测到 NVIDIA GPU，请先安装驱动"
    exit 1
fi

if ! docker info 2>/dev/null | grep -q "Runtimes.*nvidia"; then
    echo "⚠️  未检测到 nvidia-container-runtime"
    echo "    请先安装: https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/"
    echo "    安装后重启 Docker"
    exit 1
fi

# 2. 下载官方数据（注意：完整数据库 ~2TB，小批量测试可精简）
echo "[1/4] 下载 AlphaFold3 数据库（可配置精简）..."
mkdir -p data/alphafold3_databases
cd data/alphafold3_databases

# 精简测试集（约 100GB）
# 完整版见: https://github.com/google-deepmind/alphafold3#databases
if [ ! -f "uniref90.fasta" ]; then
    echo "    下载 uniref90（精简版）..."
    wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
    gunzip uniref90.fasta.gz
fi
cd ../..

# 3. 构建或拉取 Docker 镜像
echo "[2/4] 拉取 AlphaFold3 镜像..."
docker pull ghcr.io/google-deepmind/alphafold3:latest

# 4. 运行测试
echo "[3/4] 运行 AlphaFold3 测试..."
docker run --gpus all \
    -v $(pwd)/data:/data \
    -v $(pwd)/results/alphafold3:/output \
    ghcr.io/google-deepmind/alphafold3:latest \
    python /app/alphafold/run_alphafold.py \
    --input_dir /data/af3_input \
    --output_dir /output \
    --model_dir /data/alphafold3_databases/models

echo "[4/4] ✅ AlphaFold3 部署完成！"
echo "  输入文件放在 data/af3_input/"
echo "  输出结果在 results/alphafold3/"