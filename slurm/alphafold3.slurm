#!/bin/bash
# ============================================================
# AlphaFold3 复合物预测作业脚本（武大超算 Slurm + Singularity）
# 提交方式：sbatch slurm/alphafold3.slurm
#
# 武大超算不支持 Docker，用 Singularity 替代。
# 注意：AF3 数据库巨大，先确认 /project 存储配额（200GB）。
#       镜像需先在登录节点拉取（singularity pull）。
# ============================================================
#SBATCH --job-name=alphafold3
#SBATCH --output=logs/af3_%j.log
#SBATCH --error=logs/af3_%j.err
#SBATCH --partition=gpu          # 武大 GPU 分区（建议 A100 节点）
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00

ENV_PATH="${HOME}/project/conda_envs/immuno_drug_screen"
AF3_SIF="${HOME}/project/alphafold3.sif"      # Singularity 镜像放 /project
AF3_DB_DIR="${HOME}/project/af3_databases"   # AF3 数据库放 /project

# 临时目录指到 /project（避免 /home 爆盘）
export TMPDIR="${HOME}/project/tmp"
export TMP="${HOME}/project/tmp"
export TEMP="${HOME}/project/tmp"
mkdir -p "${TMPDIR}"

echo "=== AlphaFold3 作业开始 ==="
date
hostname

# 激活 Python 环境（conda 或 venv 都能激活）
if [ -f "${ENV_PATH}/bin/activate" ]; then
    source "${ENV_PATH}/bin/activate"      # venv
elif [ -d "${ENV_PATH}" ]; then
    source activate "${ENV_PATH}"           # conda
else
    echo "❌ 未找到环境: ${ENV_PATH}，请先运行 bash setup_whu.sh"
    exit 1
fi

cd ~/project/micromedicine
mkdir -p logs

# 检查 Singularity 镜像
if [ ! -f "${AF3_SIF}" ]; then
    echo "⚠ 未找到 Singularity 镜像: ${AF3_SIF}"
    echo "  请先在登录节点拉取："
    echo "  singularity pull ${AF3_SIF} docker://ghcr.io/google-deepmind/alphafold3:latest"
    exit 1
fi

# 运行 AlphaFold3（Singularity 后端，命令行参数覆盖 config）
python scripts/af3_complex_prediction.py \
  --use_singularity \
  --sif_path "${AF3_SIF}" \
  --db_dir "${AF3_DB_DIR}"

echo "=== AlphaFold3 作业完成 ==="
date