# 🧬 多免疫检查点小分子AI筛选系统

## Immuno-Checkpoint Small Molecule AI Screening Platform

针对 **PD-1、LAG-3、TIM-3、VISTA** 四大肿瘤免疫检查点靶点，实现小分子药物研发的自动化计算筛选系统。采用 **DiffDock 快速对接初筛 → AlphaFold3 精细模拟 → 多维度相互作用定量分析** 的分级策略，最终输出高优先级候选分子清单。

**✨ 新特性：**
- 内置 **3D蛋白-药物复合物交互查看器**，基于 **IgV β-三明治拓扑模板** 构建蛋白结构，80种候选药物（每个靶点20种）展示蛋白折叠 + 药物结合构象

---

## 📑 目录

- [项目架构](#项目架构)
- [技术栈](#技术栈)
- [快速开始](#快速开始)
- [项目结构](#项目结构)
- [Web可视化仪表盘](#web可视化仪表盘)
- [API接口文档](#api接口文档)
- [代码模块说明](#代码模块说明)
- [配置文件说明](#配置文件说明)
- [小分子数据库替换指南](#小分子数据库替换指南)
- [使用指南](#使用指南)
- [管线运行验证](#管线运行验证)
- [两条路线流程图](#两条路线流程图)

---

## 项目架构

```
┌─────────────────────────────────────────────────────────────────┐
│                    路线一：已知库虚拟筛选                          │
│  小分子库 → 标准化预处理 → DiffDock批量对接 → 初筛排序过滤          │
│  → AlphaFold3精细模拟 → 相互作用分析 → 终选候选                    │
├─────────────────────────────────────────────────────────────────┤
│                    路线二：全新分子设计                            │
│  活性/非活性数据集 → 特征提取 → 活性模型训练 → 新分子生成           │
│  → 活性初筛 → 接入路线一后续流程                                   │
└─────────────────────────────────────────────────────────────────┘
```

---

## 技术栈

| 类别 | 技术/工具 | 用途 |
|------|----------|------|
| 后端框架 | **Flask 3.x** | Web服务与API |
| 化学信息学 | **RDKit 2026.3** | 真实分子解析、性质计算、指纹、3D构象、PAINS过滤 |
| 机器学习 | **scikit-learn 1.9** | 真实随机森林训练（ECFP4特征） |
| 深度学习 | **PyTorch 2.11 (CUDA)** + **PyG 2.8** + **e3nn** + **esm** | DiffDock 环境（已搭建，真实推理待服务器） |
| 数据处理 | **NumPy, Pandas, SciPy** | 数据处理、统计计算 |
| 3D可视化 | **3Dmol.js** | 蛋白-药物复合物3D交互渲染 |
| 2D图表 | **Chart.js 4.4** | 交互式图表 |
| 前端 | **HTML5 + CSS3 + JavaScript + jQuery** | 深色主题仪表盘 |
| 蛋白模型 | **IgV β-三明治拓扑模板** | 9条β链+7个连接环（真实PDB待服务器下载） |

> 📌 **重要**：DiffDock 和 AlphaFold3 的真实推理需在 Linux 服务器上部署（详见 `SERVER_DEPLOYMENT.md`）。本地已完成所有代码和部署脚本。

---

## 快速开始

> 项目已配置 **conda 环境 `immuno_drug_screen`**（Python 3.13 + CUDA 12.8），无需手动激活，双击 `start_server.bat` 即可启动。

### 1. 环境要求

- conda 环境 `immuno_drug_screen`（已配置：Python 3.13 + flask/rdkit/scikit-learn/torch等全套依赖）
- 网络连接（CDN加载3Dmol.js、Chart.js、jQuery）

### 2. 启动 Web 应用

```bash
# 方式一：一键启动（推荐，自动使用 conda 环境）
start_server.bat

# 方式二：conda 命令行启动
conda activate immuno_drug_screen
python app.py

# 方式三：直接用完整路径
C:\Users\PC\miniconda3\envs\immuno_drug_screen\python.exe app.py
```

启动成功后，终端会显示：

```
  正在启动 Web 仪表盘...
  本机访问: http://127.0.0.1:5050
  局域网访问: http://<本机IP>:5050
  按 Ctrl+C 停止服务器
```

### 3. 打开浏览器

访问 **http://127.0.0.1:5050** 即可进入仪表盘。

### 4. （可选）生成3D结构数据

```bash
conda activate immuno_drug_screen
python scripts/generate_3d_complex.py
```

### 5. （可选）运行完整筛选管线

```bash
conda activate immuno_drug_screen
python scripts/data_preprocess.py
python scripts/diffdock_batch_run.py
python scripts/primary_screen_filter.py
python scripts/af3_complex_prediction.py
python scripts/interaction_analysis.py
python scripts/final_ranking.py
```

### ❗ 常见问题排查

**问题1：浏览器显示"无法连接" / "服务器拒绝了连接"**

按以下顺序检查：

```bash
# 第1步：确认 Python 进程正在运行
# Windows:
tasklist | findstr python
# 应看到至少一个 python.exe 进程

# 第2步：确认端口已被监听
netstat -ano | findstr :5050
# 应看到类似: TCP 0.0.0.0:5050 ... LISTENING

# 第3步：用命令行测试连接
python -c "import urllib.request; print(urllib.request.urlopen('http://127.0.0.1:5050/').status)"
# 应输出: 200
```

**问题2：`python app.py` 报端口占用错误 `Address already in use`**

说明已有旧服务器进程占用了 5050 端口。先结束旧进程再启动：

```bash
# Windows 查找占用5050的进程PID
netstat -ano | findstr :5050
# 找到 PID 后结束它（例如 PID 是 7416）
taskkill /PID 7416 /F

# 然后重新启动
python app.py
```

**问题3：启动后终端一闪而过或立即退出**

- 确认在**项目根目录**（`d:\AIguided_smallmolecule_design_targeting_tumor_immune_checkpoints`）
- 确认 conda 环境已激活（命令行前面有 `(immuno_drug_screen)` 标识），或用 `start_server.bat` 自动处理
- 确认依赖已安装：`conda activate immuno_drug_screen && pip list`

**问题4：能打开页面但图表/3D不显示**

前端依赖 CDN（3Dmol.js、Chart.js、jQuery），需要联网。如果网络受限，将无法加载这些资源。

**问题5：服务器在后台运行但找不到终端**

如果之前在 VSCode 或旧终端里启动过 `python app.py`，它可能仍在后台运行。运行 `netstat -ano | findstr :5050` 确认，然后刷新浏览器即可。

---

## 项目结构

```
immuno_checkpoint_screen/
├── README.md                          # 项目说明文档（本文件）
├── Process.md                         # 详细设计方案
├── app.py                             # Flask Web主应用（15+个API端点）
├── config/
│   └── config.yaml                    # 全局配置文件（靶点/数据/筛选/可视化）
├── data/
│   ├── targets/                       # 靶点蛋白结构数据
│   │   ├── {TARGET}_3d_complex.json   # 多药物3D复合物数据
│   │   ├── {TARGET}_{MOL_ID}.pdb      # 标准PDB格式复合物文件
│   │   ├── all_drug_candidates.json   # 候选药物轻量列表
│   │   └── {TARGET}_preprocessed.json # 蛋白预处理数据
│   ├── library/                       # 小分子化合物库
│   │   ├── compounds_standardized.csv # 标准化化合物数据
│   │   ├── pubchem_*.csv              # PubChem化合物
│   │   └── user_uploads/              # 用户上传的数据集
│   └── activity_dataset/              # 活性数据集（路线二）
├── scripts/                           # 17个核心Python脚本
│   ├── real_chemistry.py              # 🆕 核心真实化学计算（RDKit）
│   ├── data_preprocess.py             # 真实数据预处理（RDKit）
│   ├── diffdock_batch_run.py          # DiffDock对接（模拟+真实3D构象）
│   ├── primary_screen_filter.py       # 初筛排序过滤
│   ├── af3_complex_prediction.py      # AlphaFold3复合物预测（模拟）
│   ├── interaction_analysis.py        # 真实几何相互作用分析
│   ├── competitive_binding.py         # 🆕 竞争性抑制预测
│   ├── selectivity_analysis.py        # 🆕 选择性分析
│   ├── adme_predictor.py              # 🆕 真实类药性规则
│   ├── final_ranking.py               # 七维加权终选排序
│   ├── generate_3d_complex.py         # 3D复合物结构（IgV拓扑）
│   ├── generate_real_library.py       # 🆕 200真实化合物库
│   ├── download_real_data.py          # 🆕 真实数据下载（服务器）
│   ├── activity_model_train.py        # 真实ML模型训练（路线二）
│   ├── molecule_generation.py         # 真实RDKit组合生成（路线二）
│   ├── targetdiff_generate.py         # 🆕 TargetDiff口袋感知生成（路线二）
│   ├── pubchem_fetcher.py             # PubChem在线拉取
│   ├── pdb_fetcher.py                 # PDB元数据获取
│   └── tme_simulator.py               # 肿瘤微环境模拟
├── tools/
│   └── DiffDock/                      # 官方DiffDock源码 + 兼容shim
├── setup_server.sh                    # 🆕 DiffDock服务器部署脚本
├── setup_alphafold3.sh                # 🆕 AlphaFold3 Docker部署脚本
├── SERVER_DEPLOYMENT.md               # 🆕 服务器部署指南
├── IMPLEMENTATION_STATUS.md           # 🆕 实现状态详细说明
├── templates/                         # Web前端页面
│   ├── index.html                     # 仪表盘主页（含数据集管理）
│   ├── route1.html                    # 路线一页面
│   ├── route2.html                    # 路线二页面
│   ├── results.html                   # 结果总览（含嵌入式3D查看器）
│   └── structure.html                 # 3D交互查看器（多药物切换）
└── results/                           # 各步骤输出
    ├── diffdock/                      # 对接结果
    ├── primary_screen/                # 初筛结果
    ├── alphafold3/                    # AlphaFold3结果
    └── final_report/                  # 终选报告
```

---

## Web可视化仪表盘

系统提供5个主要页面，深色科技风主题。

### 首页仪表盘 (`/`)

| 区域 | 内容 |
|------|------|
| 统计卡片 | 靶点数(4) / 化合物库大小 / 候选分子数 / 筛选步骤(6) |
| 靶点信息卡 | PD-1 / LAG-3 / TIM-3 / VISTA 详情 |
| 📦 **数据集管理** | 当前数据集切换 + **上传CSV/SMI/SDF新数据集**并自动激活 |
| 管线控制 | 路线一/路线二运行按钮 + 实时进度条 |
| 图表区域 | 对接置信度 / pLDDT / 相互作用 / 五维雷达图 |
| Top排名表 | 跨靶点候选分子排名 |

### 路线一页面 (`/route1`)

- 6步管线流程图（CSS动画状态标记）
- 分步运行按钮：数据预处理 / DiffDock / 初筛 / AlphaFold3 / 相互作用 / 终选
- 理化性质分布、Lipinski通过环形图、DiffDock通过率柱状图、ΔG分布图

### 路线二页面 (`/route2`)

- 分子设计流程图
- 模型性能对比：Accuracy / Precision / Recall / F1 / AUC
- 生成分子候选预览表

### 结果总览页面 (`/results`)

| 区域 | 说明 |
|------|------|
| 总览横幅 | 总候选数 + 筛选方法概览 |
| 靶点汇总卡片 | 可点击切换3D查看靶点 |
| 🧊 **3D结构预览** | 嵌入式3Dmol.js查看器，左侧药物列表 + 右侧蛋白-药物复合物 |
| 🏆 综合排名表 | 点击行加载药物复合物3D结构 |
| 图表区域 | 得分分布 / 优先级分布 / 五维雷达图 |

### 3D结构查看器 (`/structure`)

三栏布局专用3D交互查看器：

| 面板 | 功能 |
|------|------|
| **左侧** | 靶点标签页 + 药物搜索 + 20个/靶点药物列表（按亲和力排序） |
| **中间** | 蛋白卡通渲染 + 口袋高亮 + 药物sticks + 工具栏 |
| **右侧** | 药物详情（ΔG/MW/成药性）+ 预测相互作用 + 口袋残基 |

**工具栏快捷键：** R=重置, C=卡通, S=线状, F=聚焦药物

**每个靶点20种候选药物：**

| 编号 | 类型 | 示例 |
|------|------|------|
| mol_001 | Aryl Ether | BMS-936558 Analog |
| mol_002 | Benzamide | LAG-525 Analog |
| mol_003 | Benzimidazole | TSR-042 Analog |
| mol_004 | Biphenyl Urea | CA-170 Analog |
| mol_005 | Quinoline | SHR-1210 Analog |
| mol_006 | Indole-piperidine | Pembrolizumab Mimetic |
| mol_007 | Thiadiazole | MBG453 Analog |
| mol_008 | Triazolopyrimidine | BMS-986207 Analog |
| mol_009 | Pyrazolopyrimidine | TSR-033 Analog |
| mol_010 | Pyridine-Benzamide | AB122 Analog |
| simple_001~010 | 类药骨架 | 哌啶/磺酰胺/吗啉等 |

---

## API接口文档

### 页面路由

| 方法 | 端点 | 说明 |
|------|------|------|
| GET | `/` | 仪表盘主页 |
| GET | `/route1` | 路线一页面 |
| GET | `/route2` | 路线二页面 |
| GET | `/results` | 结果总览页面 |
| GET | `/structure` | 3D结构查看器页面 |

### 数据API

| 方法 | 端点 | 说明 |
|------|------|------|
| GET | `/api/status` | 管线运行状态 |
| GET | `/api/dashboard_data` | 全部仪表盘数据 |
| GET | `/api/scores/<target>` | 靶点候选分子得分 |
| GET | `/api/config` | 当前配置 |
| POST | `/api/run_pipeline` | 运行完整管线 |
| POST | `/api/run_step` | 运行单个步骤 |

### 3D结构API

| 方法 | 端点 | 说明 |
|------|------|------|
| GET | `/api/3d_complex/<target>` | 靶点完整3D数据（含所有药物） |
| GET | `/api/3d_complex/<target>/<drug_id>` | 指定药物的PDB复合物 |
| GET | `/api/drug_list/<target>` | 靶点候选药物列表（轻量） |
| GET | `/api/all_drug_lists` | 所有靶点药物列表 |

### 数据集管理API

| 方法 | 端点 | 说明 |
|------|------|------|
| GET | `/api/list_datasets` | 列出所有数据集（格式/大小/来源） |
| POST | `/api/upload_dataset` | 上传CSV/SMI/SDF并自动激活 |
| POST | `/api/switch_dataset` | 切换当前数据集 |

---

## 代码模块说明

### 主线脚本（真实计算）

| # | 脚本 | 功能 | 状态 |
|---|------|------|------|
| 1 | `real_chemistry.py` | RDKit 真实化学计算核心 | ✅ 真实 |
| 2 | `data_preprocess.py` | RDKit 真实理化性质+PAINS过滤 | ✅ 真实 |
| 3 | `diffdock_batch_run.py` | DiffDock对接（模拟）+ 真实3D构象 | ⚠️ 部分 |
| 4 | `primary_screen_filter.py` | 三阶段过滤：置信度→Lipinski→位点 | ✅ 真实 |
| 5 | `af3_complex_prediction.py` | AlphaFold3预测（模拟） | ❌ 模拟 |
| 6 | `interaction_analysis.py` | 真实3D几何相互作用分析 | ✅ 真实 |
| 7 | `competitive_binding.py` | 竞争性抑制概率预测 | ✅ 真实 |
| 8 | `selectivity_analysis.py` | 跨靶点选择性指数 | ✅ 真实 |
| 9 | `adme_predictor.py` | 真实类药性规则（Lipinski/Veber等） | ✅ 真实 |
| 10 | `final_ranking.py` | 七维加权综合打分 | ✅ 真实 |
| 11 | `generate_3d_complex.py` | IgV拓扑3D结构 + 20药物 | ✅ 完成 |
| 12 | `generate_real_library.py` | 200真实类药分子库 | ✅ 真实 |
| 13 | `download_real_data.py` | 真实数据下载（服务器） | 🆕 待运行 |

### 路线二专用（真实ML / TargetDiff 生成 / 组合化学）

| # | 脚本 | 功能 | 状态 |
|---|------|------|------|
| 14 | `activity_model_train.py` | ECFP4 + scikit-learn 真实训练 | ✅ 真实 |
| 15 | `molecule_generation.py` | RDKit 化学反应组合生成 | ✅ 真实 |
| 16 | `targetdiff_generate.py` | 🆕 TargetDiff 口袋感知扩散生成（真实模型，降级到 RDKit） | ✅ 代码就绪，待权重 |

### `generate_3d_complex.py` 详解

| 属性 | 说明 |
|------|------|
| **蛋白质模型** | 手工IgV域拓扑坐标模板 — 前片A'-G-F-C-C'+C" β链，后片B-E-D β链，7个连接环 |
| **配体放置** | CC'-C loop形成的CDR-like结合口袋中心 |
| **药物模板** | 10种精细药物(BMS-936558/LAG-525/TSR-042/CA-170/Pembrolizumab等) + 10种类药片段 |
| **PDB输出** | 严格PDB v3列格式，蛋白链A(ATOM) + 配体链B(HETATM)，含CONECT键连 |
| **输出文件** | `{TARGET}_3d_complex.json`, `{TARGET}_{MOL_ID}.pdb`, `all_drug_candidates.json` |

---

## 配置文件说明

`config/config.yaml` 关键配置：

```yaml
# 数据集路径
data:
  active_library: "pubchem_all_targets.csv"
  library_dir: "data/library"
  upload_dir: "data/library/user_uploads"

# 靶点配置
targets:
  PD-1: {pdb_id: "4ZQK", binding_site_residues: [...], functional_domain: "IgV"}
  ...

# 终选五维权重
final_ranking:
  weights:
    docking_score: 0.25         # 对接置信度
    structure_confidence: 0.20  # AlphaFold3结构质量
    interaction_strength: 0.25  # 相互作用强度
    drug_likeness: 0.15         # QED成药性
    binding_site_match: 0.15    # 结合位点匹配度
```

---

## 小分子数据库替换指南

### 🆕 网页端一键管理

**无需修改代码，在仪表盘首页即可完成：**

1. 打开 `http://127.0.0.1:5050`
2. 在"📦 小分子数据集管理"区域：
   - **切换已有数据集**：下拉选择 → "🔄 切换数据集"
   - **上传新数据集**：选择CSV/SMI/SDF文件 → "⬆ 上传并激活"
3. 支持格式：`.csv`(含smiles列) / `.smi`(一行一SMILES) / `.sdf`(MOL块)

### API 方式切换

```bash
# 列出所有数据集
curl http://127.0.0.1:5050/api/list_datasets

# 切换数据集
curl -X POST http://127.0.0.1:5050/api/switch_dataset \
  -H "Content-Type: application/json" \
  -d '{"filename": "my_library.csv"}'

# 上传新数据集
curl -X POST http://127.0.0.1:5050/api/upload_dataset \
  -F "file=@my_compounds.csv" -F "name=My Custom Library"
```

### 支持的外部数据库

| 数据库 | 规模 | 使用方式 |
|--------|------|----------|
| PubChem | ~110M 化合物 | `python scripts/pubchem_fetcher.py` |
| ZINC20/22 | ~2B 化合物 | 下载SMILES文件或通过网页上传 |
| ChEMBL | ~2.4M 化合物 | 文件下载或 chembl_webresource_client |
| 自定义库 | 不限 | 网页上传CSV/SMI/SDF |

---

## 使用指南

### 🚀 快速启动

```bash
# 方式一：双击 start_server.bat（推荐，自动使用 conda 环境）

# 方式二：命令行
conda activate immuno_drug_screen
python scripts/generate_3d_complex.py             # 生成3D结构（首次）
python app.py                                     # 启动Web服务
# 浏览器打开 http://127.0.0.1:5050
```

### 📊 仪表盘操作

1. **首页** — 管理数据集、运行管线、查看统计
2. **路线一** — 逐步运行6步筛选
3. **结果** — 查看排名表 + 嵌入式3D蛋白-药物复合物预览
4. **3D结构** — 完整3D交互查看器（切换药物/缩放/旋转）

### 🧊 3D操作技巧

| 操作 | 方法 |
|------|------|
| 旋转 | 鼠标左键拖拽 |
| 缩放 | 鼠标滚轮 |
| 聚焦药物 | F键 |
| 切换样式 | C键(卡通) / S键(线状) |

### 📦 单独运行模块

```bash
python scripts/data_preprocess.py      # 预处理
python scripts/generate_3d_complex.py  # 3D结构
python scripts/diffdock_batch_run.py   # DiffDock
python scripts/primary_screen_filter.py # 初筛
python scripts/af3_complex_prediction.py # AlphaFold3
python scripts/interaction_analysis.py  # 相互作用
python scripts/final_ranking.py        # 终选报告
```

---

## 管线运行验证

完整管线实际运行结果（2026-08-06）：

| 步骤 | 模块 | 结果 |
|---|---|---|
| 数据预处理 | `data_preprocess.py` | 109个PubChem化合物, Lipinski通过率87.2% |
| 批量对接 | `diffdock_batch_run.py` | 109分子×10构象对接 |
| 初筛过滤 | `primary_screen_filter.py` | Top 67候选, 最高得分0.8391 |
| 复合物预测 | `af3_complex_prediction.py` | 54通过(80.6%), 平均pLDDT=73.33 |
| 相互作用 | `interaction_analysis.py` | 平均ΔG=-13.44 kcal/mol |
| 终选排序 | `final_ranking.py` | B类40个, C类13个 |
| 3D结构 | `generate_3d_complex.py` | 80个PDB (IgV β-三明治) |

---

## 两条路线流程图

### 路线一：已知库虚拟筛选

```
┌──────────┐    ┌──────────┐    ┌──────────┐    ┌──────────┐    ┌──────────┐    ┌──────────┐
│ Step 1   │───→│ Step 2   │───→│ Step 3   │───→│ Step 4   │───→│ Step 5   │───→│ Step 6   │
│ 数据预处理│    │ DiffDock │    │ 初筛过滤  │    │AlphaFold3│    │相互作用  │    │ 终选排序  │
│109个分子  │    │ 批量对接  │    │Top 67候选 │    │ 复合物模拟│    │ 定量分析  │    │ A/B/C/D级 │
└──────────┘    └──────────┘    └──────────┘    └──────────┘    └──────────┘    └──────────┘
```

### 路线二：全新分子设计

```
┌──────────┐    ┌──────────┐    ┌──────────┐    ┌──────────┐    ┌──────────┐    ┌──────────┐
│ Step 1   │───→│ Step 2   │───→│ Step 3   │───→│ Step 4   │───→│ Step 5   │───→│ Step 6   │
│ 活性数据集│    │ 特征提取  │    │ 模型训练  │    │ 分子生成  │    │ 活性初筛  │    │接入路线一 │
│300样本/靶点│   │80维特征   │    │随机森林   │    │500个/靶点 │    │ 概率≥0.7  │    │DiffDock+  │
└──────────┘    └──────────┘    └──────────┘    └──────────┘    └──────────┘    └──────────┘
```

---

## ⚠️ 注意事项

1. **本地计算已真实化：** RDKit 化学计算、scikit-learn ML、几何相互作用分析均为真实实现
2. **DiffDock/AlphaFold3 需服务器：** 真实推理需 Linux 服务器（详见 `SERVER_DEPLOYMENT.md`、`setup_server.sh`、`setup_alphafold3.sh`）
3. **真实数据需服务器网络：** 真实 PDB 结构、PubChem 抑制剂需在服务器上运行 `download_real_data.py`
4. **3D蛋白结构：** 当前为 IgV 拓扑模板，真实 PDB 待服务器下载
5. **CDN依赖：** 前端依赖 3Dmol.js、Chart.js、jQuery 的 CDN
6. **假阳性控制：** 筛选结果需经生化/细胞实验验证

## 🖥️ 服务器部署

真实 DiffDock 和 AlphaFold3 部署到 Linux 服务器的完整步骤：

```bash
# 1. 按迁移清单拷贝文件（见 MIGRATION.md）
python verify_migration.py     # （可选）本地验证迁移文件完整性

# 2. 在服务器上执行
bash setup_server.sh          # DiffDock 一键部署（依赖见 requirements_server.txt）
bash setup_alphafold3.sh      # AlphaFold3 Docker 部署
python scripts/download_real_data.py  # 下载真实数据
```

**相关文档：**
- `MIGRATION.md` — 迁移清单（哪些拷、哪些不拷、如何拷、如何运行）
- `SERVER_DEPLOYMENT.md` — 部署指南
- `IMPLEMENTATION_STATUS.md` — 实现状态说明
- `requirements_server.txt` — 服务器依赖清单

---

**License:** For research use only. | **Design Document:** See `Process.md` for full specification.
