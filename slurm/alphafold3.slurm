#!/bin/bash
# ============================================================
# AlphaFold3 复合物预测作业脚本（武大超算 Slurm + Singularity）
# 提交方式：sbatch slurm/alphafold3.slurm
#
# 武大超算不支持 Docker，用 Singularity 替代。
# 注意：AF3 数据库巨大（数百GB），需要先确认存储空间，
#       且 Singularity 镜像需提前拉取。
# ============================================================
#SBATCH --job-name=alphafold3
#SBATCH --output=logs/af3_%j.log
#SBATCH --error=logs/af3_%j.err
#SBATCH --partition=gpu          # ← 改成实际 GPU 分区名
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00

CONDA_MODULE="anaconda"
CONDA_ENV="immuno_drug_screen"
AF3_SIF="$HOME/alphafold3.sif"   # Singularity 镜像路径（需提前拉取）
AF3_DB_DIR="$HOME/af3_databases" # AF3 数据库目录

echo "=== AlphaFold3 作业开始 ==="
date
hostname

module load ${CONDA_MODULE}
source activate ${CONDA_ENV}

# 进入项目目录
cd ~/project/micromedicine

# 检查 Singularity 镜像
if [ ! -f "${AF3_SIF}" ]; then
    echo "⚠ 未找到 AlphaFold3 Singularity 镜像: ${AF3_SIF}"
    echo "  请先拉取镜像："
    echo "  singularity pull alphafold3.sif docker://ghcr.io/google-deepmind/alphafold3:latest"
    exit 1
fi

# 运行 AlphaFold3（通过 Singularity，替代 docker run）
python scripts/af3_complex_prediction.py \
  --use_singularity \
  --sif_path "${AF3_SIF}" \
  --db_dir "${AF3_DB_DIR}"

echo "=== AlphaFold3 作业完成 ==="
date