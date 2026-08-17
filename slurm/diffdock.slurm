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

# 可配置项（按需改）
CONDA_MODULE="anaconda"
ENV_PATH="${HOME}/project/conda_envs/immuno_drug_screen"   # conda 环境装在 /project

echo "=== DiffDock 作业开始 ==="
date
hostname

# 加载 module + 激活 conda 环境（--prefix 路径）
module load ${CONDA_MODULE} 2>/dev/null || module load anaconda/3.7
source activate "${ENV_PATH}"

# 确认 GPU 可见
nvidia-smi
python -c "import torch; print('CUDA available:', torch.cuda.is_available()); print('Device:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'N/A')"

# 进入项目目录（在 /project 下）
cd ~/project/micromedicine

# 建日志目录
mkdir -p logs

# 1. 快速预筛（CPU 即可，登录节点可先做，也可放这里）
# python scripts/rapid_prefilter.py

# 2. 真实 DiffDock 对接（GPU 重计算，必须作业提交）
python scripts/diffdock_batch_run.py

# 3. 后续精细流程（轻量，可放同一作业）
python scripts/primary_screen_filter.py
python scripts/interaction_analysis.py
python scripts/final_ranking.py

echo "=== DiffDock 作业完成 ==="
date