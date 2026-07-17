# 单细胞 RNA-seq 分析项目模板

这是从实际项目中整理出的 Seurat 单细胞分析代码模板。仓库只保存可复用代码和项目配置框架，不保存原始数据、运行结果或某个项目的临时分析。

核心约定：

- `src/`：通用函数和可复用分析方法。
- `config/config.R`：当前项目需要调整的参数。
- `workflow/`：项目开始时为空；每实际运行一个分析点，新建一个 `stepXX.名称.R`。
- 已经运行过的 workflow 不回头改；参数变化或补充分析使用新的步骤文件记录。
- workflow 只负责加载配置、加载方法和调用函数，保持简短。

## 开始一个新项目

推荐在 GitHub 中把本仓库设为 **Template repository**，之后通过 **Use this template** 创建项目仓库。

克隆新项目后：

1. 把 Cell Ranger 输出放到 `data/<sample>/outs/`。
2. 修改 `config/config.R` 中的分组、QC、归一化、聚类和子聚类参数。
3. 在项目根目录启动 R。
4. 按分析进度在 `workflow/` 新建步骤文件。

第一个分析脚本可以写成：

```r
rm(list = ls())
setwd("你的项目根目录")
source("config/config.R")
source("src/load_all.R")

step_qc_before_filter()
```

保存为 `workflow/step01.QC_before_filter.R`。运行后检查并修改自动生成的 `config/qc_thresholds.csv`，再新建：

```r
rm(list = ls())
setwd("你的项目根目录")
source("config/config.R")
source("src/load_all.R")

step_qc_filter()
```

后续可调用的标准分析点见 [新项目指南](docs/new_project.md)。

## 目录

```text
.
├── config/          # 项目参数与人工表格
├── data/            # 原始矩阵和中间对象（不提交）
├── docs/            # 使用说明和项目报告
├── logs/            # 运行日志（不提交）
├── results/         # 图表和结果表（不提交）
├── src/             # 通用函数与可复用分析方法
└── workflow/        # 当前项目实际运行的 stepXX 脚本
```

## 环境

Conda 环境：

```bash
conda env create -f environment.yml
conda activate scrna-workflow
```

代码以 Seurat 5 为主要运行环境；Harmony 和 ggalluvial 仅在对应方法启用时需要。

## 发布前

本仓库没有替你选择开源许可证。公开发布前，请根据代码复用范围补充合适的 `LICENSE`。
