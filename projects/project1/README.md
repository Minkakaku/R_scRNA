# project1

这是模板中的示例项目层。复制为新项目名后，只修改本目录中的配置和 workflow。

- `config/config.yaml`：项目分组、比较和子聚类目标；覆盖根目录公共参数。
- `config/samples.tsv`：QC 步骤生成的样本和分组表，人工确认后继续分析。
- `config/qc_thresholds.tsv`：逐样本 QC 阈值，运行 02 前人工检查。
- `config/cluster_annotation.tsv`：主聚类人工注释。
- `config/subcluster_annotation.tsv`：子聚类人工注释。
- `data/raw/`：Cell Ranger 输出，推荐布局为 `raw/<sample>/outs/`。
- `data/processed/`：按步骤编号的 QS2 中间对象。
- `results/`：按分析点编号的最终图和表。
- `logs/`：运行日志。

从仓库根目录按编号运行 workflow；QC 阈值和人工注释确认后再进入下一步。不要用一键脚本跨过这些人工检查点。
