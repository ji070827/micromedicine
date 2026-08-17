#!/bin/bash
# ============================================================
# DiffDock 真实对接作业脚本（武大超算 Slurm）
# 提交方式：sbatch slurm/diffdock.slurm
#
# 武大 GPU 分区：gpu（V100: g0001-g0135 / A100: g0136-g0179）
# 注意：先运行 sinfo 确认分区名和可用 GPU 类型
# ============================================================
#SBATCH --job-name=diffdock
#SBATCH --output=logs/diffdock_%j.log
#SBATCH --error=logs/diffdock_%j.err
#SBATCH --partition=gpu          # 武大 GPU 分区名
#SBATCH --gres=gpu:1             # 请求 1 块 GPU
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00

# 环境路径（与 setup_whu.sh 一致）
ENV_PATH="${HOME}/project/conda_envs/immuno_drug_screen"

# 临时目录指到 /project（避免运行时 /home 爆盘）
export TMPDIR="${HOME}/project/tmp"
export TMP="${HOME}/project/tmp"
export TEMP="${HOME}/project/tmp"
mkdir -p "${TMPDIR}"

echo "=== DiffDock 作业开始 ==="
date
hostname

# 激活 Python 环境（conda 或 venv 都能激活）
if [ -f "${ENV_PATH}/bin/activate" ]; then
    # venv 环境
    source "${ENV_PATH}/bin/activate"
elif [ -d "${ENV_PATH}" ]; then
    # conda 环境
    source activate "${ENV_PATH}"
else
    echo "❌ 未找到环境: ${ENV_PATH}"
    echo "  请先运行: bash setup_whu.sh"
    exit 1
fi

# 确认 GPU 可见
nvidia-smi
python -c "import torch; print('CUDA available:', torch.cuda.is_available()); print('Device:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'N/A')"

# 进入项目目录
cd ~/project/micromedicine
mkdir -p logs

# 真实 DiffDock 对接（GPU 重计算，必须作业提交）
python scripts/diffdock_batch_run.py

# 后续精细流程
python scripts/primary_screen_filter.py
python scripts/interaction_analysis.py
python scripts/final_ranking.py

echo "=== DiffDock 作业完成 ==="
date