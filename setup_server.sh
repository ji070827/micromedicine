#!/usr/bin/env bash
# ============================================================
# DiffDock + AlphaFold3 服务器部署脚本 (Linux)
# 在配置好 GPU + 正常网络的 Linux 服务器上运行
# 使用: bash setup_server.sh
#
# 注意：本仓库已包含 tools/DiffDock 源码（含 torch_scatter/torch_cluster shim），
#       无需重复 clone，只需装依赖 + 下载权重 + 开启真实模型开关。
# ============================================================
set -euo pipefail

echo "=============================================="
echo "  多免疫检查点小分子AI筛选 - 服务器部署"
echo "=============================================="

# ---------- 1. 系统依赖 ----------
echo "[1/6] 安装系统依赖..."
sudo apt-get update -y
sudo apt-get install -y \
    git wget curl \
    openbabel \
    python3 python3-pip python3-venv

# ---------- 2. Python 环境 ----------
echo "[2/6] 创建 Python 环境..."
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
echo "[3/6] 下载 DiffDock 预训练权重 (~4.5GB)..."
mkdir -p tools/DiffDock/models
cd tools/DiffDock
python -c "
from utils.download import download_and_extract
download_and_extract('https://github.com/gcorso/DiffDock/releases/latest/download/diffdock_models.zip', 'models/')
"
cd ../..

# ---------- 4. 开启真实模型开关 ----------
echo "[4/6] 开启真实 DiffDock 模型开关..."
python - <<'PYEOF'
import yaml
from pathlib import Path

cfg_path = Path("config/config.yaml")
with open(cfg_path, 'r', encoding='utf-8') as f:
    cfg = yaml.safe_load(f)

# 开启真实 DiffDock
cfg['screening']['diffdock']['use_real_model'] = True
# 默认保持 AF3 为模拟（需先完成 setup_alphafold3.sh + 数据库）
# 用户如需开启真实 AF3，手动改为 true 即可

with open(cfg_path, 'w', encoding='utf-8') as f:
    yaml.safe_dump(cfg, f, allow_unicode=True, sort_keys=False)

print("  ✅ use_real_model（DiffDock）已设为 true")
PYEOF

# ---------- 5. 下载真实数据 ----------
echo "[5/6] 下载真实蛋白质结构 + 真实小分子抑制剂库..."
python scripts/download_real_data.py

# ---------- 6. 完成 ----------
echo "[6/6] 完成"
echo "=============================================="
echo "  ✅ 服务器部署完成！"
echo ""
echo "  当前状态："
echo "    - DiffDock 真实模型：已启用 (use_real_model=true)"
echo "    - AlphaFold3：仍为模拟（要启用需先跑 setup_alphafold3.sh）"
echo ""
echo "  运行真实对接："
echo "    source venv_diffdock/bin/activate"
echo "    python scripts/diffdock_batch_run.py"
echo ""
echo "  启用真实 AlphaFold3（需 Docker + 数据库）："
echo "    bash setup_alphafold3.sh"
echo "    然后手动把 config/config.yaml 里 screening.af3.use_real_model 改为 true"
echo "=============================================="