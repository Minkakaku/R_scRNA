# R-scrna-workflow-template

这是可移植的 Seurat 单细胞主流程模板。复用函数、公共配置和具体项目分开保存；历史参考代码仍保留在 `source_catalog/`，不参与正式运行。

```text
.
├── src/R/                         # 通用函数和复用分析步骤
├── config/                        # 公共默认参数和 marker 参考
├── projects/project1/
│   ├── config/                    # 项目分组、QC 阈值和人工注释
│   ├── workflow/                  # 按实际运行顺序编号的短入口
│   ├── data/{raw,processed,reference}/
│   ├── results/                   # 按分析点编号的图表
│   ├── logs/
│   └── docs/reports/
├── envs/
├── tests/
└── source_catalog/                # 只读历史实现目录
```

## 开始项目

1. 复制 `projects/project1/` 并改成项目名。
2. 把 Cell Ranger 输出放到 `data/raw/<sample>/outs/`。
3. 修改项目 `config/config.yaml`；公共默认值在根目录 `config/default.yaml`。
4. 从仓库根目录按编号运行项目 `workflow/`。
5. 在 QC 阈值、主注释和子注释处先人工确认对应 TSV，再继续下一步。

单个入口保持简短：

```r
source("src/R/project_config.R")
load_project("projects/project1")
step_global_clustering()
```

已经执行过的 workflow 不覆盖；参数变化或补充分析新建下一个编号脚本。模板不提供跨越人工检查点的 `run_all`。

## 环境与测试

```bash
conda env create -f envs/environment.yaml
conda activate scrna-workflow
Rscript tests/test_project_config.R
```

代码面向 Seurat 5。PCA 完成后会同时调用 `PCAtools::findElbowPoint()` 和 `nFactors::nScree()` 输出自动拐点参考，但不会覆盖 YAML 中人工确认的聚类维度。

## 输出约定

- 中间对象保存到项目 `data/processed/`，按步骤编号。
- 最终图表保存到项目 `results/01_qc` 至 `06_subcluster`。
- 图默认同名输出 PNG/PDF，重要结果同时保留支撑表。
- 原始数据、中间对象、results 和 logs 不提交到模板仓库。

## 历史代码

`source_catalog/` 收录 609 份 R/R Markdown/Shell 语义副本；`SOURCE_CATALOG.tsv` 记录场景、角色、owner 和 SHA-256。正式分析优先使用 `src/R/`，不要把历史脚本中的内部路径和固定生物学判断复制到复用层。
