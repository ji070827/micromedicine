#!/usr/bin/env bash
# ============================================================
# 武大超算 Swarm 集群部署脚本（Slurm + module + conda，无 sudo）
#
# 适用：武汉大学超算 登录节点 swarm01/02/03
# 关键约束（来自武大官方说明）：
#   1. 无 root 权限，不用 sudo
#   2. /home 配额仅 20GB（试用），/project 200GB
#      → conda 环境必须装到 /project，否则爆 /home！
#   3. 重计算必须 sbatch 提交，禁止登录节点直接跑
#   4. 不支持 Docker，用 Singularity
#
# 使用：bash setup_whu.sh
# ============================================================
set -euo pipefail

# 可配置项
ENV_PREFIX="${HOME}/project/conda_envs"   # conda 环境装在 /project（关键！避免 /home 爆盘）
ENV_NAME="immuno_drug_screen"
CONDA_MODULE="anaconda"                   # 官方示例：module load anaconda/3.7
PYTHON_VERSION="3.10"

echo "=============================================="
echo "  武大超算 Swarm 集群部署（conda→/project）"
echo "=============================================="

# ---------- 1. 加载 module ----------
echo "[1/5] 加载 anaconda 模块..."
module load ${CONDA_MODULE} 2>/dev/null || module load anaconda/3.7 2>/dev/null || {
    echo "  ❌ module load anaconda 失败"
    echo "     请先运行 'module avail' 确认可用的 anaconda 模块名"
    echo "     然后修改本脚本 CONDA_MODULE 变量"
    exit 1
}
echo "  ✅ anaconda 已加载"

# ---------- 2. 创建 conda 环境（在 /project 下，关键！） ----------
echo "[2/5] 创建 conda 环境到 /project（避免 /home 20GB 爆盘）..."
mkdir -p "${ENV_PREFIX}"

ENV_PATH="${ENV_PREFIX}/${ENV_NAME}"
if [ -d "${ENV_PATH}" ]; then
    echo "  ℹ 环境已存在: ${ENV_PATH}"
else
    # 接受 conda ToS（若需要）
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main 2>/dev/null || true
    echo "  （首次创建环境需几分钟，请耐心等待）"
    conda create --prefix "${ENV_PATH}" python=${PYTHON_VERSION} pip -y
fi

# 激活环境（用 --prefix 路径激活）
source activate "${ENV_PATH}"
echo "  ✅ conda 环境就绪: ${ENV_PATH}"

# ---------- 3. 安装依赖（pip） ----------
echo "[3/5] 安装 Python 依赖..."
pip install --upgrade pip

# 基础依赖（纯 Python）
pip install numpy pandas scipy pyyaml scikit-learn rdkit flask tqdm requests

# 深度学习（CUDA 版 torch）
# ⚠ 武大 GPU 分区有 V100/A100，CUDA 版本需先确认（module avail cuda 或 nvidia-smi）
#    V100 兼容 cu118；A100 可用 cu121
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu118

# DiffDock 依赖（PyG C++ 扩展）
pip install e3nn fair-esm ema-pytorch torchmetrics prody
pip install torch-geometric
pip install torch-scatter torch-sparse torch-cluster torch-spline-conv \
    -f https://pytorch-geometric.com/whl/torch-2.2.0+cu118.html

echo "  ✅ 依赖安装完成"

# ---------- 4. 开启真实模型开关 ----------
echo "[4/5] 开启真实 DiffDock 开关..."
python - <<'PYEOF'
import yaml
from pathlib import Path
cfg_path = Path("config/config.yaml")
with open(cfg_path, 'r', encoding='utf-8') as f:
    cfg = yaml.safe_load(f)
cfg['screening']['diffdock']['use_real_model'] = True
cfg['screening']['af3']['backend'] = 'singularity'   # 武大无 Docker，用 Singularity
with open(cfg_path, 'w', encoding='utf-8') as f:
    yaml.safe_dump(cfg, f, allow_unicode=True, sort_keys=False)
print("  ✅ use_real_model(DiffDock)=true, af3.backend=singularity")
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
echo "  ⚠ 重要：重计算必须 sbatch 提交，禁止登录节点直接跑！"
echo ""
echo "  运行 DiffDock 对接："
echo "    sbatch slurm/diffdock.slurm"
echo ""
echo "  运行 AlphaFold3（Singularity，需先拉镜像）："
echo "    singularity pull alphafold3.sif docker://ghcr.io/google-deepmind/alphafold3:latest"
echo "    sbatch slurm/alphafold3.slurm"
echo ""
echo "  交互式调试（srun 小任务）："
echo "    srun -p gpu -n 1 --gres=gpu:1 --pty bash"
echo "    source activate ${ENV_PATH}"
echo "=============================================="