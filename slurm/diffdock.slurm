#!/bin/bash
# ============================================================
# DiffDock 真实对接作业脚本（武大超算 Slurm）
# 提交方式：sbatch slurm/diffdock.slurm
#
# 注意：需要先确认你的 GPU 分区名（用 sinfo 查看），
#       把 --partition= 和 --gres= 改成实际值。
# ============================================================
#SBATCH --job-name=diffdock
#SBATCH --output=logs/diffdock_%j.log
#SBATCH --error=logs/diffdock_%j.err
#SBATCH --partition=gpu          # ← 改成实际 GPU 分区名
#SBATCH --gres=gpu:1             # 请求 1 块 GPU
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00

# 可配置项（按需改）
CONDA_MODULE="anaconda"
CONDA_ENV="immuno_drug_screen"

echo "=== DiffDock 作业开始 ==="
date
hostname

# 加载 module + 激活 conda 环境
module load ${CONDA_MODULE}
source activate ${CONDA_ENV}

# 确认 GPU 可见
nvidia-smi
python -c "import torch; print('CUDA available:', torch.cuda.is_available()); print('Device:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'N/A')"

# 进入项目目录（按实际路径改）
cd ~/project/micromedicine

# 1. 快速预筛（CPU 即可，但如果已在登录节点做过可跳过）
# python scripts/rapid_prefilter.py

# 2. 真实 DiffDock 对接（GPU）
python scripts/diffdock_batch_run.py

# 3. 后续精细流程（CPU/轻量）
python scripts/primary_screen_filter.py
python scripts/interaction_analysis.py
python scripts/final_ranking.py

echo "=== DiffDock 作业完成 ==="
date