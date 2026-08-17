#!/usr/bin/env bash
# ============================================================
# 武大超算 Swarm 集群部署脚本（Slurm + module + conda，无 sudo）
#
# 关键约束（武大官方说明）：
#   1. /home 配额仅 20GB，/project 200GB → 缓存/临时/环境全部放 /project
#   2. 重计算必须 sbatch 提交，禁止登录节点直接跑
#   3. 无 Docker，用 Singularity
#
# 使用：bash setup_whu.sh
# ============================================================
set -euo pipefail

ENV_PREFIX="${HOME}/project/conda_envs"
ENV_NAME="immuno_drug_screen"
PYTHON_VERSION="3.10"

# ===== 关键修复：把 pip/TMP 缓存全部指到 /project（避免 /home 爆盘） =====
export TMPDIR="${HOME}/project/tmp"
export TMP="${HOME}/project/tmp"
export TEMP="${HOME}/project/tmp"
export PIP_CACHE_DIR="${HOME}/project/pip_cache"
export CONDA_PKGS_DIRS="${HOME}/project/conda_pkgs"
mkdir -p "${TMPDIR}" "${PIP_CACHE_DIR}" "${CONDA_PKGS_DIRS}"

echo "=============================================="
echo "  武大超算部署（缓存/临时全部→/project）"
echo "=============================================="

# ---------- 1. 自动探测 conda / 用 venv 兜底 ----------
echo "[1/5] 探测 Python 环境方案..."
CONDA_MODULE=""
for m in miniconda3 anaconda3 anaconda python/anaconda miniconda conda; do
    if module load "$m" 2>/dev/null; then
        CONDA_MODULE="$m"
        break
    fi
done

if [ -z "${CONDA_MODULE}" ]; then
    CONDA_MODULE=$(module avail 2>&1 | grep -iE 'conda|anaconda|miniconda' | head -1 | awk '{print $1}' || true)
    [ -n "${CONDA_MODULE}" ] && module load "${CONDA_MODULE}" 2>/dev/null || CONDA_MODULE=""
fi

if [ -z "${CONDA_MODULE}" ] && command -v conda >/dev/null 2>&1; then
    CONDA_MODE="conda"
elif [ -n "${CONDA_MODULE}" ]; then
    CONDA_MODE="conda"
else
    CONDA_MODE="venv"
fi
echo "  ✅ 方案: ${CONDA_MODE}"

# ---------- 2. 创建 Python 环境 ----------
echo "[2/5] 创建 Python 环境到 /project..."
mkdir -p "${ENV_PREFIX}"
ENV_PATH="${ENV_PREFIX}/${ENV_NAME}"

if [ "${CONDA_MODE}" = "conda" ]; then
    if [ ! -d "${ENV_PATH}" ]; then
        conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main 2>/dev/null || true
        conda create --prefix "${ENV_PATH}" python=${PYTHON_VERSION} pip -y
    fi
    source activate "${ENV_PATH}"
else
    if [ ! -d "${ENV_PATH}" ]; then
        python3 -m venv "${ENV_PATH}"
    fi
    source "${ENV_PATH}/bin/activate"
fi
echo "  ✅ 环境就绪: ${ENV_PATH}"

# ---------- 3. 安装依赖（pip，缓存已在 /project） ----------
echo "[3/5] 安装依赖（ProDy 可选，失败不中断）..."
pip install --upgrade pip

# 基础依赖
pip install numpy pandas scipy pyyaml scikit-learn rdkit flask tqdm requests

# 深度学习（CUDA 版 torch）
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu118

# DiffDock 依赖（PyG C++ 扩展）
pip install e3nn fair-esm ema-pytorch torchmetrics
pip install torch-geometric
pip install torch-scatter torch-sparse torch-cluster torch-spline-conv \
    -f https://pytorch-geometric.com/whl/torch-2.2.0+cu118.html

# ProDy 是可选依赖（全原子模型才用，粗粒度模型不需要），编译失败不中断
pip install prody || echo "  ⚠ ProDy 编译失败（非核心依赖，跳过，不影响粗粒度 DiffDock）"

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
cfg['screening']['af3']['backend'] = 'singularity'
with open(cfg_path, 'w', encoding='utf-8') as f:
    yaml.safe_dump(cfg, f, allow_unicode=True, sort_keys=False)
print("  ✅ use_real_model(DiffDock)=true")
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
echo "  ✅ 部署完成！"
echo ""
echo "  运行 DiffDock 对接："
echo "    sbatch slurm/diffdock.slurm"
echo ""
echo "  AlphaFold3（Singularity，需先拉镜像）："
echo "    singularity pull ~/project/alphafold3.sif docker://ghcr.io/google-deepmind/alphafold3:latest"
echo "    sbatch slurm/alphafold3.slurm"
echo "=============================================="