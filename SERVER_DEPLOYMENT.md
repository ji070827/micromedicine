# 服务器部署指南

> 本文档说明如何在 **Linux 服务器**上部署真实的 DiffDock 和 AlphaFold3 模型，以及下载真实数据。

---

## 一、部署条件

| 需求 | DiffDock | AlphaFold3 |
|------|----------|-----------|
| GPU | NVIDIA GPU ≥8GB | NVIDIA GPU ≥16GB（推荐 A100/4090） |
| 系统 | Linux | Linux |
| 网络 | 可访问 GitHub / RCSB / PubChem | 可访问 Google / UniProt |
| 磁盘 | ~10GB（权重+模型） | ~2TB（完整数据库）或 ~100GB（精简） |

---

## 二、部署步骤

### 步骤 1：DiffDock 部署（约 30 分钟）

```bash
# 在服务器项目根目录执行
bash setup_server.sh
```

这个脚本会自动：
1. 安装系统依赖（git、openbabel、cuda-toolkit）
2. 创建 Python 虚拟环境
3. 安装 PyTorch (CUDA) + PyG + e3nn + esm + 所有依赖
4. 克隆 DiffDock 仓库
5. 下载预训练权重
6. 下载真实蛋白结构 + 真实小分子库

### 步骤 2：AlphaFold3 部署（约 2-4 小时）

```bash
# 完整数据库很大，先确保磁盘空间足够
bash setup_alphafold3.sh
```

这个脚本会：
1. 检查 NVIDIA Docker 环境
2. 下载 AlphaFold3 数据库（可配置精简 vs 完整）
3. 拉取官方 Docker 镜像
4. 运行测试

### 步骤 3：下载真实数据

```bash
# 下载真实 PDB 结构 + 已知抑制剂
python scripts/download_real_data.py
```

---

## 三、真实数据来源

### 3.1 真实蛋白质结构（RCSB PDB）

| 靶点 | PDB ID | 说明 |
|------|--------|------|
| PD-1 | 4ZQK | PD-1/PD-L1 复合物晶体结构 |
| LAG-3 | 7TZH | LAG-3 |
| TIM-3 | 5F71 | TIM-3 |
| VISTA | 6OIL | VISTA |

### 3.2 真实小分子库

**已知免疫检查点抑制剂（PubChem）：**

| 靶点 | 代表性化合物 | PubChem CID |
|------|-------------|-------------|
| PD-1 | BMS-936558 相关 | 23629198, 25023587 |
| PD-1 | 其他 | 447290, 49867904 |
| LAG-3 | LAG525 相关 | 137347988, 72721919 |
| TIM-3 | 相关抑制剂 | 444899, 16074 |
| VISTA | 相关抑制剂 | 54694254 |

**可扩展的数据源：**
- PubChem（~110M 化合物，关键词搜索）
- ZINC20/22（~2B 可购化合物）
- ChEMBL（~2.4M 带活性数据）
- BindingDB（~2.4M 结合数据）

---

## 四、部署后的运行流程

```bash
# 1. 激活环境
source venv_diffdock/bin/activate

# 2. 运行完整路线一（真实 DiffDock + AF3）
python scripts/data_preprocess.py      # 预处理（RDKit 真实计算）
python scripts/diffdock_batch_run.py   # 真实 DiffDock 对接
python scripts/primary_screen_filter.py # 初筛过滤
python scripts/af3_complex_prediction.py # 真实 AlphaFold3
python scripts/interaction_analysis.py  # 相互作用分析
python scripts/final_ranking.py        # 终选排序

# 3. 运行路线二
python scripts/activity_model_train.py  # 活性模型训练（真实 ML）
python scripts/molecule_generation.py   # 分子生成（RDKit 组合化学）

# 4. 启动 Web
python app.py
```

---

## 五、关键注意事项

1. **CUDA 版本匹配**：服务器 GPU 的 CUDA 版本需与 PyTorch 安装命令匹配（setup_server.sh 中默认 cu118，需按需调整）
2. **网络**：所有真实数据下载需在国际网络正常的环境执行
3. **磁盘空间**：AlphaFold3 数据库巨大，先检查磁盘
4. **GPU 内存**：DiffDock 单个复合物 ~2-4GB，AF3 需要更大
5. **Windows 限制**：本项目的真实 DiffDock 在 Windows 上无法完整运行（C++ 扩展 + 网络限制），必须 Linux

---

## 六、当前本地环境状态（供参考）

| 组件 | 状态 |
|------|------|
| 真实化学计算（RDKit） | ✅ 已实现并验证 |
| 真实机器学习（scikit-learn） | ✅ 已实现并验证 |
| 真实相互作用分析（几何） | ✅ 已实现 |
| 真实 DiffDock GPU | ⚠️ 环境已就绪，权重下载受阻于网络 |
| 真实 AlphaFold3 | ❌ 需服务器 Docker |
| 真实蛋白结构 | ❌ 需服务器网络 |
| 真实小分子库 | ⚠️ 已有 200 个真实类药分子 + 需补充已知抑制剂 |