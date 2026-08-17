# 武汉大学超算（Swarm 集群）部署指南

> 本文档说明如何在武汉大学超算集群上部署和运行本项目。
> 武大超算使用 **Slurm 作业调度 + module 环境 + Singularity**，与普通独立服务器有本质区别。

---

## 一、武大超算关键信息

从官方文档（https://docs.hpc.whu.edu.cn/）确认：

| 项目 | 内容 |
|------|------|
| **登录节点** | `swarm01.whu.edu.cn`（202.114.96.179）、`swarm02.whu.edu.cn`（202.114.96.180） |
| **专用传输节点** | `swarm-xfe.whu.edu.cn`（202.114.96.177） |
| **登录方式** | 仅支持 SSH（Windows 自带 OpenSSH / Xshell / MobaXterm） |
| **校外访问** | 必须先连武大 SSLVPN（vpn.whu.edu.cn） |
| **权限** | 无 root，不能用 sudo |
| **容器** | 不支持 Docker，用 **Singularity** |
| **IDE** | 不支持 VSCode 服务端，仅命令行 |

---

## 二、关键差异：为什么不能用普通 setup_server.sh

| 普通服务器 | 武大超算 |
|-----------|---------|
| `sudo apt-get install` | ❌ 无 root，改用 `module load` |
| Docker | ❌ 改用 Singularity |
| 直接在终端跑重计算 | ❌ 必须 `sbatch` 提交作业 |
| pip 装系统级 | ✅ 装用户级 conda 环境 |

---

## 三、部署步骤

### 步骤 1：登录

```bash
# 本机先连 VPN（校外），然后 SSH
ssh 你的用户名@swarm01.whu.edu.cn
```

### 步骤 2：传代码到 /project

**方式 A：git clone（若登录节点能访问 GitHub）**
```bash
cd ~/project
git clone https://github.com/ji070827/micromedicine.git
cd micromedicine
```

**方式 B：本地打包 + scp（若不能访问 GitHub）**
```bash
# 本机打包（排除 venv/results/.git）
tar -czf micromedicine.tar.gz app.py config scripts templates static tools data/library data/activity_dataset setup_whu.sh slurm requirements_server.txt *.md

# 传到专用传输节点
scp micromedicine.tar.gz 用户名@swarm-xfe.whu.edu.cn:~/project/

# 登录后解压
ssh 用户名@swarm01.whu.edu.cn
cd ~/project && tar -xzf micromedicine.tar.gz
```

> XFTP 8：主机填 `swarm-xfe.whu.edu.cn`，协议 SFTP，端口 22。

### 步骤 3：运行部署脚本

```bash
cd ~/project/micromedicine
bash setup_whu.sh
```

**这个脚本会：**
1. `module load anaconda` 加载系统 conda
2. 创建用户级 conda 环境 `immuno_drug_screen`
3. 装所有 Python 依赖（torch GPU 版 + PyG 生态 + rdkit）
4. 开启真实 DiffDock 模型开关
5. 下载 DiffDock 权重 + 真实数据

### 步骤 4：先确认环境

```bash
# 查看 GPU 分区名和可用状态
sinfo

# 查看你的账号可用的分区
sacctmgr show assoc format=Account,User,Partition,QOS
```

**把 `slurm/diffdock.slurm` 和 `slurm/alphafold3.slurm` 里的 `--partition=` 和 `--gres=` 改成你实际可用的分区名/GPU 类型**（如 `--partition=GPU --gres=gpu:1`）。

---

## 四、提交作业运行

### 运行 DiffDock 对接（GPU）

```bash
sbatch slurm/diffdock.slurm
```

### 运行 AlphaFold3（需先拉 Singularity 镜像）

```bash
# 先拉镜像（在登录节点执行，只拉一次）
singularity pull alphafold3.sif docker://ghcr.io/google-deepmind/alphafold3:latest

# 提交作业
sbatch slurm/alphafold3.slurm
```

### 查看作业

```bash
squeue -u 你的用户名   # 查看排队/运行状态
sacct                  # 查看历史作业
tail -f logs/*.log     # 查看输出
scancel 作业ID         # 取消作业
```

---

## 五、重要提醒

1. **登录节点禁重计算**：DiffDock/AF3 的 GPU 推理必须 `sbatch` 提交，不要直接在 swarm01/02 上跑，否则会被杀或封号。
2. **先确认 GPU 分区**：不同账号分区不同，用 `sinfo` 和 `sacctmgr show assoc` 查你的可用分区，再改 slurm 脚本。
3. **数据及时下载**：超算存储不宜长期存放结果，算完尽快下载到本地。
4. **AlphaFold3 数据库巨大**：完整 AF3 数据库数百 GB，需确认存储配额，或只用精简版。
5. **CUDA 版本匹配**：`setup_whu.sh` 默认 cu118，需用 `module avail cuda` 查集群 CUDA 版本后调整 torch 安装命令。

---

## 六、完整命令速查

```bash
# 登录
ssh 用户名@swarm01.whu.edu.cn

# 部署
cd ~/project/micromedicine && bash setup_whu.sh

# 查分区
sinfo

# 提交 DiffDock 作业
sbatch slurm/diffdock.slurm

# 提交 AF3 作业（需先 singularity pull）
sbatch slurm/alphafold3.slurm

# 查看作业
squeue -u 用户名