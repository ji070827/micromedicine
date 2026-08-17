#!/usr/bin/env bash
# ============================================================
# DiffDock + AlphaFold3 + TargetDiff 服务器部署脚本 (Linux)
# 在配置好 GPU + 正常网络的 Linux 服务器上运行
# 使用: bash setup_server.sh
#
# 注意：本仓库已包含 tools/DiffDock 源码（含 torch_scatter/torch_cluster shim），
#       无需重复 clone，只需装依赖 + 下载权重 + 开启真实模型开关。
#       TargetDiff 需要单独 clone + 下载权重（权重下载失败时路线二会自动降级到 RDKit 组合化学）。
# ============================================================
set -euo pipefail

echo "=============================================="
echo "  多免疫检查点小分子AI筛选 - 服务器部署"
echo "=============================================="

# ---------- 1. 系统依赖 ----------
echo "[1/8] 安装系统依赖..."
sudo apt-get update -y
sudo apt-get install -y \
    git wget curl \
    openbabel \
    python3 python3-pip python3-venv

# ---------- 2. Python 环境 ----------
echo "[2/8] 创建 Python 环境..."
python3 -m venv venv_diffdock
source venv_diffdock/bin/activate

pip install --upgrade pip

# PyTorch (CUDA 版，按服务器 GPU 的 CUDA 版本调整，默认 cu118)
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu118

# DiffDock 依赖（普通纯 Python 包）
pip install e3nn==0.5.0
pip install torch-geometric==2.2.0
pip install torch-scatter torch-sparse torch-cluster \
    -f https://pytorch-geometric.com/whl/torch-2.2.0+cu118.html
pip install fair-esm[esmfold]==2.0.0
pip install rdkit pandas scipy scikit-learn prody torchmetrics openbabel pyyaml flask

# ---------- 3. 下载 DiffDock 预训练权重 ----------
echo "[3/8] 下载 DiffDock 预训练权重 (~4.5GB)..."
mkdir -p tools/DiffDock/models
cd tools/DiffDock
python -c "
from utils.download import download_and_extract
download_and_extract('https://github.com/gcorso/DiffDock/releases/latest/download/diffdock_models.zip', 'models/')
"
cd ../..

# ---------- 4. 部署 TargetDiff (路线二分子生成) ----------
echo "[4/8] 部署 TargetDiff（口袋感知分子生成）..."
if [ ! -d "tools/TargetDiff" ]; then
    git clone https://github.com/guanjq/targetdiff.git tools/TargetDiff || \
        echo "  ⚠ TargetDiff 克隆失败，路线二将使用 RDKit 组合化学降级"
fi

if [ -d "tools/TargetDiff" ]; then
    # 安装 TargetDiff 依赖（torchdrug 较重，失败不中断）
    (pip install torchdrug hydra-core omegaconf || echo "  ⚠ TargetDiff 依赖安装失败")

    # 下载预训练权重（TargetDiff 权重在 Google Drive / HuggingFace，网络受限可能失败）
    mkdir -p tools/TargetDiff/checkpoints
    echo "  ℹ TargetDiff 预训练权重需手动下载（Google Drive / HuggingFace），"
    echo "    下载后放到 tools/TargetDiff/checkpoints/ 下。"
    echo "    若未下载权重，路线二会自动降级到 RDKit 组合化学。"
fi

# ---------- 5. 开启真实模型开关 ----------
echo "[5/8] 开启真实模型开关..."
python - <<'PYEOF'
import yaml
from pathlib import Path

cfg_path = Path("config/config.yaml")
with open(cfg_path, 'r', encoding='utf-8') as f:
    cfg = yaml.safe_load(f)

# 开启真实 DiffDock
cfg['screening']['diffdock']['use_real_model'] = True
# 若 TargetDiff 目录存在，则开启真实 TargetDiff
if Path("tools/TargetDiff").exists():
    cfg['route2']['targetdiff']['use_real_model'] = True
# 默认保持 AF3 为模拟（需先完成 setup_alphafold3.sh + 数据库）

with open(cfg_path, 'w', encoding='utf-8') as f:
    yaml.safe_dump(cfg, f, allow_unicode=True, sort_keys=False)

print("  ✅ use_real_model（DiffDock）已设为 true")
print("  ✅ use_real_model（TargetDiff）已设为 true（若已部署）")
PYEOF

# ---------- 6. 下载真实数据 ----------
echo "[6/8] 下载真实蛋白质结构 + 真实小分子抑制剂库..."
python scripts/download_real_data.py

# ---------- 7. 完成 ----------
echo "[7/8] 完成"
echo "=============================================="
echo "  ✅ 服务器部署完成！"
echo ""
echo "  当前状态："
echo "    - DiffDock 真实模型：已启用 (use_real_model=true)"
echo "    - TargetDiff 真实模型：已启用（若已部署 + 权重下载成功）"
echo "    - AlphaFold3：仍为模拟（要启用需先跑 setup_alphafold3.sh）"
echo ""
echo "  运行完整路线一："
echo "    source venv_diffdock/bin/activate"
echo "    python scripts/diffdock_batch_run.py"
echo ""
echo "  运行完整路线二（TargetDiff 生成 + 共享筛选）："
echo "    source venv_diffdock/bin/activate"
echo "    python scripts/targetdiff_generate.py"
echo "    python scripts/diffdock_batch_run.py"
echo ""
echo "  启用真实 AlphaFold3（需 Docker + 数据库）："
echo "    bash setup_alphafold3.sh"
echo "    然后手动把 config/config.yaml 里 screening.af3.use_real_model 改为 true"
echo ""
echo "  ⚠ TargetDiff 预训练权重说明："
echo "     Official: https://github.com/guanjq/targetdiff"
echo "     权重下载后放 tools/TargetDiff/checkpoints/，否则路线二自动降级 RDKit 组合化学"
echo "=============================================="