# 项目部署到服务器完整流程（Git 方案）

> 本文档说明如何把项目推送到你自己的 Git 仓库，再从服务器克隆并在服务器上运行真实模型。

---

## 一、整体流程概览

```
本地电脑（开发）          Git 仓库          Linux 服务器（SSH）
     │                      │                    │
     │ ① git push           │                    │
     ├─────────────────────→│                    │
     │                      │ ② git clone        │
     │                      │←───────────────────┤
     │                      │                    │ ③ 部署 + 运行
```

---

## 二、阶段 1：把项目推送到 Git 仓库（本地电脑）

### 前提
在 GitHub（或 Gitee）网页上创建一个**空仓库**，不要勾选"用 README 初始化"。

### 操作步骤
```bash
# ① .gitignore 已创建，排除了 venv/、results/、*.pdb、缓存等

# ② 配置提交身份（改成你自己的）
git config user.name "你的名字"
git config user.email "你的邮箱"

# ③ 添加远程仓库（改成你的仓库地址）
git remote add origin https://github.com/用户名/仓库名.git

# ④ 提交所有代码
git add .
git commit -m "初始提交：多免疫检查点小分子AI筛选系统"

# ⑤ 推送到远程
git push -u origin main
```

---

## 三、会提交 / 不会提交的文件

### ✅ 会提交到 Git 仓库

| 内容 | 说明 |
|------|------|
| `app.py` | Flask 主应用 |
| `config/config.yaml` | 配置文件 |
| `scripts/*.py` | 18 个核心脚本 |
| `templates/*.html` | 5 个前端页面 |
| `static/` | 静态资源 |
| `tools/DiffDock/` | DiffDock 源码 + 兼容 shim |
| `setup_server.sh`、`setup_alphafold3.sh` | 部署脚本 |
| `requirements_server.txt` | 服务器依赖清单 |
| `*.md` 全部文档 | 说明文档 |
| `start_server.bat` | 本地一键启动 |
| `data/library/*.csv` | 化合物库（200个真实分子） |
| `data/library/*.json` | 理化性质统计 |
| `data/targets/*.json` | 靶点配置 + 3D 复合物数据 |
| `data/activity_dataset/` | 路线二数据 |

### ❌ 不提交（被 .gitignore 排除）

| 内容 | 原因 |
|------|------|
| `venv/` | 本地虚拟环境（服务器重新建） |
| `results/` | 运行结果（可重新生成） |
| `data/targets/*.pdb` | 3D 结构 PDB（可重新生成） |
| `data/targets/real_structures/` | 下载的真实结构（体积大） |
| `logs/`、`__pycache__/`、`*.log` | 缓存和日志 |

---

## 四、阶段 2：服务器克隆项目（SSH）

```bash
# 登录服务器
ssh 用户名@服务器IP

# 克隆代码
git clone https://github.com/用户名/仓库名.git project
cd project
```

---

## 五、阶段 3：服务器部署 + 运行

```bash
# 1. 部署 DiffDock（装依赖 + 下载权重 ~4.5GB）
bash setup_server.sh

# 2. 部署 AlphaFold3（可选，需 GPU≥16GB + 数百GB磁盘）
bash setup_alphafold3.sh

# 3. 下载真实蛋白结构 + 真实抑制剂
python scripts/download_real_data.py

# 4. 用 tmux 跑管线（防 SSH 断线）
tmux new -s run
python scripts/data_preprocess.py
python scripts/diffdock_batch_run.py
python scripts/primary_screen_filter.py
python scripts/af3_complex_prediction.py
python app.py
```

> **tmux 用法**：`Ctrl+B` 然后 `D` 脱离（任务继续跑），`tmux attach -s run` 重新挂回。

---

## 六、重要提醒

1. **网络前提**：服务器必须能访问 GitHub（克隆代码）+ RCSB/PubChem（下载数据）。如果服务器也在受限网络，需先配代理或走镜像。
2. **仓库设为私有**：如果项目含未公开的数据，建议在 GitHub 上把仓库设为 Private。
3. **大文件**：DiffDock 权重（4.5GB）和 AF3 数据库（数百GB）**不要**提交到 git，它们由 `setup_server.sh` 和 `setup_alphafold3.sh` 在服务器上直接下载。