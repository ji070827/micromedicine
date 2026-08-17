# 武汉大学超算（Swarm 集群）部署指南

> 本文档依据武大超算中心官方说明（http://hpc.whu.edu.cn）整理，说明如何在武大 Swarm 集群上部署和运行本项目。
> 武大超算使用 **SLURM 作业调度 + module 环境 + Singularity**，与普通独立服务器有本质区别。

---

## 一、武大超算关键信息（官方确认）

| 项目 | 内容 |
|------|------|
| **登录地址** | `swarm.whu.edu.cn` 或 `202.114.96.180`，端口 22 |
| **登录节点** | swarm01、swarm02、swarm03（仅用于编译/提交/文件管理，**禁止运行计算**） |
| **文件传输** | SFTP，校内地址 `202.114.96.177`（校外需先连 VPN） |
| **校外访问** | 必须先连武大 VPN（https://vpn.whu.edu.cn） |
| **操作系统** | Rocky 9.4 Linux |
| **调度系统** | SLURM |
| **文件系统** | Lustre 并行文件系统 |
| **权限** | 无 root，不能用 sudo |
| **容器** | Singularity（**不支持 Docker**） |
| **软件管理** | module（`module load anaconda/3.7` 等） |

---

## 二、存储配额（关键！）

武大超算的分区配额（**决定了 conda 环境必须装哪里**）：

| 分区 | 用途 | 试用用户配额 | 说明 |
|------|------|-------------|------|
| `/home` | 家目录（环境变量/代码/编译产物） | 20GB | **容量小，不要装 conda 环境！** |
| `/project` | 项目数据区（计算数据/作业） | 200GB | **~ `~/project`，conda 环境和权重都放这里** |
| `/scratch` | 临时数据（3个月清理） | 无 | 仅临时文件 |

> ⚠️ **重要**：PyTorch GPU 版 + PyG 生态 + esm 全套依赖约 10-20GB，conda 环境如果装在默认的 `/home` 会立即爆盘。**必须把 conda 环境装到 `/project`**，本项目 `setup_whu.sh` 已自动处理。

---

## 三、GPU 分区（官方确认）

| GPU 类型 | 节点范围 | 数量 | 说明 |
|---------|---------|------|------|
| Tesla V100 | g0001-g0135 | 135 台（540 卡） | 第一代 GPU |
| Nvidia A100 | g0136-g0179 | 44 台（216 卡） | 第二代 GPU（**DiffDock 推荐用这个**） |
| Hopper 高端 GPU | 13 台 | 120 卡 | 第三代 |

GPU 分区名统一为 **`gpu`**，通过 `--gres=gpu:1` 请求。

---

## 四、部署步骤

### 步骤 1：登录（先连 VPN）

```bash
# 校外先连 VPN（vpn.whu.edu.cn），然后 SSH
ssh 你的用户名@swarm.whu.edu.cn
```

### 步骤 2：传代码到 /project

**方式 A：git clone（若登录节点能访问 GitHub）**
```bash
cd ~/project
git clone https://github.com/ji070827/micromedicine.git
cd micromedicine
```

**方式 B：本地打包 + SFTP（若不能访问 GitHub）**
```bash
# 本机打包（排除 venv/results/.git）
tar -czf micromedicine.tar.gz app.py config scripts templates static tools data/library data/activity_dataset setup_whu.sh slurm requirements_server.txt *.md

# XFTP/FileZilla 传到 202.114.96.177 的 ~/project/ 下
# 登录后解压
cd ~/project && tar -xzf micromedicine.tar.gz
```

> XFTP 8 / FileZilla：主机填 `202.114.96.177`（或 swarm.whu.edu.cn），协议 SFTP，端口 22。

### 步骤 3：运行部署脚本

```bash
cd ~/project/micromedicine
bash setup_whu.sh
```

**这个脚本自动完成：**
1. `module load anaconda` 加载系统 conda
2. 在 **`~/project/conda_envs/`** 下创建 conda 环境（避免 /home 爆盘）
3. 装所有 Python 依赖（torch GPU 版 + PyG + rdkit + esm）
4. 开启 DiffDock 真实模型 + af3 切换 Singularity 后端
5. 下载 DiffDock 权重 + 真实数据

### 步骤 4：确认环境

```bash
# 查看 GPU 分区状态
sinfo

# 查看可用模块
module avail | grep -i anaconda
module avail | grep -i cuda

# 查看磁盘配额
/bin/myDiskQuota
```

---

## 五、提交作业运行

### 运行 DiffDock 对接（GPU）

```bash
cd ~/project/micromedicine
sbatch slurm/diffdock.slurm
```

### 运行 AlphaFold3（Singularity）

```bash
# 先在登录节点拉 Singularity 镜像（只拉一次，放 /project）
singularity pull ~/project/alphafold3.sif docker://ghcr.io/google-deepmind/alphafold3:latest

# 提交作业
sbatch slurm/alphafold3.slurm
```

### 作业管理命令

```bash
squeue -u 你的用户名              # 查看排队/运行状态
sacct                             # 查看历史作业
tail -f logs/*.log                # 查看输出
scancel 作业ID                    # 取消作业
scontrol show job 作业ID          # 查看作业详情
```

---

## 六、重要提醒（官方约束）

1. **登录节点禁重计算**：DiffDock/AF3 的 GPU 推理必须 `sbatch` 提交，禁止在 swarm01/02/03 直接跑，否则账号可能被限制。
2. **conda 环境装 /project**：这是最容易踩的坑，/home 只有 20GB，装 GPU 版 torch + PyG 必爆盘。
3. **CUDA 版本**：武大有 V100 和 A100。V100 支持 cu11.x，A100 支持 cu11.x/cu12.x。`setup_whu.sh` 默认 cu118（两者都兼容），如需调整用 `module avail cuda` 查看。
4. **数据及时下载**：超算不提供备份，算完尽快下载结果到本地。
5. **成果致谢**（官方要求）：
   - 中文：本论文的数值计算得到了武汉大学超级计算中心的计算支持和帮助。
   - 英文：The numerical calculations in this paper have been done on the supercomputing system in the Supercomputing Center of Wuhan University.

---

## 七、完整命令速查

```bash
# 登录（先连 VPN）
ssh 用户名@swarm.whu.edu.cn

# 部署
cd ~/project/micromedicine && bash setup_whu.sh

# 查分区/配额
sinfo
/bin/myDiskQuota

# 提交 DiffDock 作业
sbatch slurm/diffdock.slurm

# AF3（先拉镜像）
singularity pull ~/project/alphafold3.sif docker://ghcr.io/google-deepmind/alphafold3:latest
sbatch slurm/alphafold3.slurm

# 查看作业
squeue -u 用户名

# 交互式调试（小任务，srun）
srun -p gpu -n 1 --gres=gpu:1 --pty bash
source activate ~/project/conda_envs/immuno_drug_screen