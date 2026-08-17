#!/usr/bin/env bash
# ============================================================
# 武大超算 Swarm 集群部署脚本（Slurm + module + conda，无 sudo）
#
# 适用：武汉大学超算 swarm01/02 登录节点
# 特点：
#   1. 无 root 权限，不用 sudo，依赖装用户级 conda 环境
#   2. 用 module 加载系统预装的 Python/Anaconda/CUDA
#   3. 重计算（DiffDock/AF3）用 sbatch 提交，不在登录节点直接跑
#
# 使用：bash setup_whu.sh
# ============================================================
set -euo pipefail

# 可配置项
CONDA_ENV_NAME="immuno_drug_screen"
CONDA_MODULE="anaconda"        # 武大预装 conda 的 module 名（可能为 anaconda/miniconda3）
PYTHON_VERSION="3.10"          # 超算建议用 3.10（兼容 DiffDock 生态）

echo "=============================================="
echo "  武大超算 Swarm 集群部署（用户级，无 sudo）"
echo "=============================================="

# ---------- 1. 加载 module ----------
echo "[1/5] 加载系统预装模块..."
module load ${CONDA_MODULE} 2>/dev/null || {
    echo "  ⚠ module load ${CONDA_MODULE} 失败"
    echo "    请先运行 'module avail' 查看可用的 conda/anaconda 模块名"
    echo "    然后修改本脚本开头的 CONDA_MODULE 变量"
    exit 1
}
echo "  ✅ 已加载 ${CONDA_MODULE}"

# ---------- 2. 创建 conda 环境 ----------
echo "[2/5] 创建 conda 环境 ${CONDA_ENV_NAME} (Python ${PYTHON_VERSION})..."
if conda env list | grep -q "^${CONDA_ENV_NAME} "; then
    echo "  ℹ 环境已存在，跳过创建"
else
    # 接受 conda 服务条款（若需要）
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main 2>/dev/null || true
    conda create -n ${CONDA_ENV_NAME} python=${PYTHON_VERSION} pip -y
fi

# 激活环境
source activate ${CONDA_ENV_NAME}
echo "  ✅ conda 环境就绪"

# ---------- 3. 安装依赖（pip） ----------
echo "[3/5] 安装 Python 依赖..."
pip install --upgrade pip

# 基础依赖（纯 Python，无 C++ 编译）
pip install numpy pandas scipy pyyaml scikit-learn rdkit flask tqdm requests

# 深度学习（CUDA 版 torch，超算用模块加载的 CUDA，需匹配版本）
# 注：武大 GPU 分区的 CUDA 版本需先确认（module avail cuda），一般 cu118/cu121
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu118

# DiffDock 依赖（PyG 生态，需从专门源装 C++ 扩展）
pip install e3nn fair-esm ema-pytorch torchmetrics prody
pip install torch-geometric
pip install torch-scatter torch-sparse torch-cluster torch-spline-conv \
    -f https://pytorch-geometric.com/whl/torch-2.2.0+cu118.html

echo "  ✅ 依赖安装完成"

# ---------- 4. 开启真实模型开关 ----------
echo "[4/5] 开启真实模型开关..."
python - <<'PYEOF'
import yaml
from pathlib import Path
cfg_path = Path("config/config.yaml")
with open(cfg_path, 'r', encoding='utf-8') as f:
    cfg = yaml.safe_load(f)
cfg['screening']['diffdock']['use_real_model'] = True
with open(cfg_path, 'w', encoding='utf-8') as f:
    yaml.safe_dump(cfg, f, allow_unicode=True, sort_keys=False)
print("  ✅ use_real_model (DiffDock) 已设为 true")
PYEOF

# ---------- 5. 下载 DiffDock 权重 + 真实数据 ----------
echo "[5/5] 下载 DiffDock 权重 + 真实数据..."
mkdir -p tools/DiffDock/models
cd tools/DiffDock
python -c "
from utils.download import download_and_extract
download_and_extract('https://github.com/gcorso/DiffDock/releases/latest/download/diffdock_models.zip', 'models/')
" || echo "  ⚠ DiffDock 权重下载失败（可能网络受限），可稍后手动下载"
cd ../..

python scripts/download_real_data.py || echo "  ⚠ 真实数据下载失败（可能网络受限）"

echo "=============================================="
echo "  ✅ 武大超算部署完成！"
echo ""
echo "  重要：重计算必须用 sbatch 提交，不要在登录节点直接跑！"
echo ""
echo "  运行 DiffDock 对接（提交作业）："
echo "    sbatch slurm/diffdock.slurm"
echo ""
echo "  运行 AlphaFold3（Singularity）："
echo "    sbatch slurm/alphafold3.slurm"
echo ""
echo "  激活环境（交互式调试用）："
echo "    module load ${CONDA_MODULE}"
echo "    source activate ${CONDA_ENV_NAME}"
echo "=============================================="