# 服务器迁移清单

> 本文件说明从本地 Windows 电脑迁移到 Linux 服务器时，**哪些文件需要拷贝**、哪些不需要、以及拷贝后如何运行。

---

## 一、需要拷贝的文件（核心源码）

### 必拷目录

```
项目根目录/
├── app.py                          # Flask 主应用
├── config/
│   └── config.yaml                 # 配置文件
├── scripts/                        # 全部 Python 脚本（17个）
├── templates/                      # Web 前端页面（5个）
├── static/                         # 静态资源（css/js）
├── tools/
│   └── DiffDock/                   # DiffDock 源码（含兼容 shim）
├── setup_server.sh                 # 服务器部署脚本
├── setup_alphafold3.sh             # AlphaFold3 部署脚本
├── requirements_server.txt         # 服务器依赖清单
├── SERVER_DEPLOYMENT.md            # 部署指南
├── IMPLEMENTATION_STATUS.md        # 实现状态说明
├── README.md                       # 项目说明
├── Process.md                      # 设计方案
└── MIGRATION.md                    # 本文件
```

### 需要拷贝的数据文件

```
data/
├── targets/
│   ├── {TARGET}_preprocessed.json  # 蛋白预处理信息
│   ├── {TARGET}_3d_complex.json    # 3D 复合物数据
│   └── all_drug_candidates.json    # 候选药物列表
├── library/
│   ├── pubchem_all_targets.csv     # 真实化合物库（200个）
│   └── compounds_standardized.csv  # 标准化库
└── activity_dataset/
    └── （路线二模型数据，可选）
```

> **注意**：`data/targets/{TARGET}_*.pdb` 和 `data/library/pubchem_*.csv` 是本地生成的可再生文件，可不拷（服务器上可重新生成）。但 `pubchem_all_targets.csv` 建议拷，避免重复生成。

---

## 二、不需要拷贝的文件

| 目录/文件 | 原因 |
|-----------|------|
| `venv/` | 本地虚拟环境，服务器需重新创建（Linux 与 Windows 不兼容） |
| `results/` | 本地运行结果，服务器会重新生成 |
| `logs/` | 本地日志 |
| `tools/DiffDock/__pycache__/` | Python 缓存 |
| `*.pyc` | 编译缓存 |
| `data/library/user_uploads/` | 用户上传的数据集（按需拷贝） |

---

## 三、拷贝方式

### 方式 1：Git（推荐，如果本地已有 git 仓库）

```bash
# 在服务器上
git clone <你的仓库地址>
```

### 方式 2：压缩打包拷贝

在本地 Windows 执行：
```bash
# 打包核心文件（排除 venv、results、logs 等）
tar -czf project_source.tar.gz \
    app.py config scripts templates static tools \
    data/targets/*.json data/library/pubchem_all_targets.csv \
    data/library/compounds_standardized.csv \
    setup_server.sh setup_alphafold3.sh \
    requirements_server.txt *.md
```

然后 scp 到服务器：
```bash
scp project_source.tar.gz user@server:/opt/project/
```

### 方式 3：直接拷贝整个目录（简单但体积大）

```
整个项目目录拷到服务器 → 删除 venv → 重新创建环境
```

---

## 四、服务器上的运行步骤（拷贝完成后）

```bash
# 1. 解压（如果用了打包）
cd /opt/project
tar -xzf project_source.tar.gz

# 2. 部署 DiffDock（自动安装依赖+克隆+下载权重）
bash setup_server.sh

# 3. 部署 AlphaFold3（可选，需要 GPU≥16GB + 空间）
bash setup_alphafold3.sh

# 4. 下载真实数据
python scripts/download_real_data.py

# 5. 运行完整管线
python scripts/data_preprocess.py
python scripts/diffdock_batch_run.py   # 真实 GPU 对接
python scripts/primary_screen_filter.py
python scripts/af3_complex_prediction.py
python scripts/interaction_analysis.py
python scripts/final_ranking.py

# 6. 启动 Web
python app.py
```

---

## 五、重要提醒

1. **CUDA 版本**：`setup_server.sh` 默认 cu118，需根据服务器 GPU 调整 `requirements_server.txt` 中的 torch 版本
2. **torch-geometric 扩展**：需要从专门源安装（`pytorch-geometric.com`），版本需与 torch 严格匹配
3. **DiffDock 权重**：`setup_server.sh` 会自动下载，但在服务器上需能访问 GitHub releases
4. **权限**：服务器上可能需要 sudo 权限安装系统依赖
5. **磁盘**：AlphaFold3 完整数据库 ~2TB，DiffDock 权重 ~4.5GB

---

## 六、快速验证清单

服务器上部署完成后，运行以下命令验证：

```bash
# 验证 GPU
python -c "import torch; print(torch.cuda.is_available())"  # True

# 验证 DiffDock 依赖
python -c "import torch_geometric, e3nn, esm; print('OK')"

# 验证 RDKit
python -c "from rdkit import Chem; print('OK')"

# 验证 OpenBabel
python -c "import openbabel; print('OK')"

# 验证 Web
python app.py  # 浏览器访问 http://<服务器IP>:5050
```

全部通过即部署成功。