library(tidyverse)
library(Seurat)

library(hdWGCNA)


# https://mp.weixin.qq.com/s/qxgFsP2WF7s-cUKHaDx9gw

set.seed(12345)

enableWGCNAThreads(nThreads = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb

outdir <- "wgcna_output"
dir.create(outdir)

## read data  -------
# all_data <- readRDS("all_data_subcluster.rds")
# annotation_cols <- "annotation3"
# 
# cds <- readRDS("monocle_project_batch_3/cds.rds")
# meta_data_cds <- Biobase::pData(cds)[c("annotation_full2", "State")]

# epi
sp_epi <- readRDS("sp_subcluster_LUSC/sp_epi.rds")

# wgcna 
all_data_wgcna <- sp_epi[, sp_epi$annotation_epi == 'ALG11 Epithelial']
# all_data_wgcna <- all_data_wgcna %>% AddMetaData(meta_data_cds)
all_data_wgcna <- all_data_wgcna %>% JoinLayers()

# add useful columns
all_data_wgcna$cell_type <- all_data_wgcna$annotation_epi
all_data_wgcna$Sample <- str_c(all_data_wgcna$project, all_data_wgcna$group, sep = "_")
  
  
# data load  -------
DimPlot(sp_epi, group.by = "annotation_epi", label = T) # + umap_theme() + NoLegend()
DimPlot(all_data_wgcna, group.by = "cell_type", label = T) # + umap_theme() + NoLegend()


# seurat为V5 ---------
seurat_obj_version <- SeuratObject::Version(all_data_wgcna) %>% package_version() %>% .$major
seurat_version <- packageDescription("Seurat")$Version %>% package_version() %>% .$major

if(seurat_obj_version < seurat_version){
  all_data_wgcna <- UpdateSeuratObject(all_data_wgcna)
}

# construct wgcna --------
# gene_select: fraction variable all custom
# fraction需传入fraction 
# custom需传入gene_list 
all_data_wgcna <- SetupForWGCNA(
  all_data_wgcna,
  gene_select = "variable", # "fraction",
  fraction    = 0.05,
  wgcna_name  = "ALG11 Epithelial" # hdWGCNA实验名称
)

all_data_wgcna %>% GetWGCNAGenes() # 12639


# construct meta cells ---------

# 构造meta cells，将相似的细胞表达量进行合并，减少过于稀疏的单细胞数据带来的问题
# 一般以分组的方式进行
undebug(MetacellsByGroups)
all_data_wgcna <- MetacellsByGroups(
  seurat_obj = all_data_wgcna,
  group.by = c("cell_type", "Sample"),  # 将来自同一Sample的同种cell_type构造成metacell
  reduction = "harmony", # 用于执行KNN算法，此处harmony为通过harmony校正过的PCA
  k = 25, # KNN参数
  min_cells = 100,
  max_shared = 10, # 两个metacell可共享的最大细胞数
  ident.group = "cell_type"  # Idents名
)


# 标准化metacell表达矩阵
all_data_wgcna <- NormalizeMetacells(all_data_wgcna)

# 36671 cells -> 4000 meta cells
all_data_wgcna@misc$`ALG11 Epithelial`$wgcna_metacell_obj


# construct WGCNA --------------

# set expression
all_data_wgcna <- SetDatExpr(
  all_data_wgcna,
  group_name = "ALG11 Epithelial", # 挑选感兴趣的细胞类型 
  group.by = "cell_type", # TFF3 Goblet所在的metadata列名
  assay = "RNA",
  slot = "data"
)

# select soft power
all_data_wgcna <- TestSoftPowers(
  all_data_wgcna,
  networkType = "signed" # 此外，还可选择“unsigned”或“signed hybrid”参数
)

plot_list <- PlotSoftPowers(all_data_wgcna)

patchwork::wrap_plots(plot_list, ncol=2)

# export soft power
pdf(
  file.path(outdir, "1 soft power selection.pdf"),
  width  = 15,
  height = 12
)
patchwork::wrap_plots(plot_list, ncol=2)
dev.off()

# construct co-expression network
all_data_wgcna <- ConstructNetwork(
  all_data_wgcna, 
  soft_power = 7, # 上图中黑色圆圈中的最佳软阈值
  setDatExpr = FALSE,
  deepSplit = 3,
  detectCutHeight = 0.995,
  minModuleSize = 30,
  mergeCutHeight = 0.2,
  tom_name = 'ALG11 Epithelial', # TOM矩阵名称，运行后会在当前目录生成一个tom_outdir指定的文件夹保存TOM矩阵
  overwrite_tom = TRUE,
  tom_outdir = outdir
)

PlotDendrogram(all_data_wgcna, main='ALG11 Epithelial hdWGCNA Dendrogram')

# export nentwork
pdf(
  file.path(outdir, "2 wgcna network.pdf"),
  width  = 10,
  height = 6
)
PlotDendrogram(all_data_wgcna, main='hdWGCNA Dendrogram')
dev.off()


# 读入TOM矩阵
TOM <- GetTOM(all_data_wgcna)

# 查看TOM矩阵元素
TOM[1:4, 1:4]


# compute ME and kME ----
# 若根据sample去批次，则需要先进行scale，harmony前需要进行正态化（Scale）
all_data_wgcna <- ScaleData(all_data_wgcna, features = VariableFeatures(all_data_wgcna))

# 计算所有的MEs，比较耗时
all_data_wgcna <- ModuleEigengenes(
  all_data_wgcna,
  group.by.vars = "Sample" # 根据Sample去批次
)


# compute eigengene-based connectivity (kME):
all_data_wgcna <- ModuleConnectivity(
  all_data_wgcna,
  group.by = 'cell_type', 
  group_name = 'ALG11 Epithelial'  # 这里继续计算（感兴趣的细胞类型）的kME
)


# 重命名module，以表明这些module由" "计算所得
all_data_wgcna <- ResetModuleNames(
  all_data_wgcna,
  new_name = "ALG11 Epithelial"
)

# 可视化每个模块中的基因的kME
PlotKMEs(all_data_wgcna, ncol = 5, n_hubs = 10)


# export kME
pdf(
  file.path(outdir, "3 wgcna kME plot.pdf"),
  width  = 15,
  height = 12
)
PlotKMEs(all_data_wgcna, ncol = 5, n_hubs = 10)
dev.off()


# get the module assignment table and add module score ----------------
modules <- GetModules(all_data_wgcna)

# 通过GetHubGenes函数抽取每个模块的topN hub节点，即每个模块中kME高的topN个基因
hub_df <- GetHubGenes(seurat_obj = all_data_wgcna, n_hubs = 20)
saveRDS(hub_df, file.path(outdir, "hub_df.rds"))


# 根据每个模块的topN个基因对细胞进行打分
all_data_wgcna <- ModuleExprScore(
  all_data_wgcna,
  n_genes = 25, # topN hub genes
  method = "Seurat" # 打分方法，分为Seurat、UCell
)

# 至此，基本完成了hdWGCNA分析
saveRDS(all_data_wgcna, file = file.path(outdir, 'hdWGCNA_object.rds'))

# visualization -------------

## FeaturePlot - hMEs -------------
all_data_wgcna[['umap2']] <- 
  CreateDimReducObject(
    embeddings = sp_epi[['harmonyumap']] %>% Embeddings() %>% .[colnames(all_data_wgcna), ],
    assay = "RNA"
  )

plot_list <- ModuleFeaturePlot(
  all_data_wgcna,
  reduction = "umap2",
  features = 'hMEs', # 可选择MEs、hMEs、scores、average
  order = TRUE # order so the points with highest hMEs are on top
)

#  stitch together with patchwork
patchwork::wrap_plots(plot_list, ncol = 6)


# export feature plot
pdf(
  file.path(outdir, "4 wgcna module feature plot.pdf"),
  width  = 18,
  height = 12
)
patchwork::wrap_plots(plot_list, ncol = 6)
dev.off()


## module correlation -------------
ModuleCorrelogram(all_data_wgcna)

pdf(
  file.path(outdir, "4 module correlation.pdf"),
  width  = 12,
  height = 12
)
ModuleCorrelogram(all_data_wgcna)
dev.off()



## module feature plot ------
undebug(ModuleFeaturePlot)
plot_list <- ModuleFeaturePlot(
  all_data_wgcna,
  features = 'hMEs', # 可选择MEs、hMEs、scores、average
  reduction = "full.harmonyumap",
  order = TRUE # order so the points with highest hMEs are on top
)
#  stitch together with patchwork
patchwork::wrap_plots(plot_list, ncol = 6)

## dotplot -------------

# get hMEs from seurat object
MEs <- GetMEs(all_data_wgcna, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
all_data_wgcna@meta.data <- cbind(all_data_wgcna@meta.data, MEs)

p_dot <-
  DotPlot(
    all_data_wgcna[, all_data_wgcna$cell_type == "ALG11 Epithelial"], 
    features = mods, 
    group.by = "group"
  ) + 
  RotatedAxis() + # 旋转坐标轴标签
  scale_color_gradient2(high = "red", mid = "grey95", low = "blue")  # 改颜色
p_dot2 <-
  DotPlot(
    all_data_wgcna[, all_data_wgcna$cell_type == "ALG11 Epithelial"], 
    features = mods, 
    group.by = "cell_type"
  ) + 
  RotatedAxis() + # 旋转坐标轴标签
  scale_color_gradient2(high = "red", mid = "grey95", low = "blue")  # 改颜色

# export
pdf(
  file.path(outdir, "4 dotplot.pdf"),
  width  = 12,
  height = 5
)
p_dot
p_dot2
dev.off()

## vlnplot -------------
p_vln <- VlnPlot(
  all_data_wgcna[, all_data_wgcna$cell_type == "ALG11 Epithelial"],
  features = c("brown", "blue", "green"),
  group.by = "group",
  pt.size = 0
)

p_vln <-
  p_vln &
  geom_boxplot(width = .25, fill = "white", outliers = FALSE) &
  labs(x = "", y = "hME") &
  NoLegend()

# export
pdf(
  file.path(outdir, "4 vlnplot.pdf"),
  width  = 12,
  height = 12
)
p_vln
dev.off()



## network  --------------

mod_name <- hdWGCNA::GetModules(all_data_wgcna)$module %>% unique() %>% as.character()
ModuleNetworkPlot(seurat_obj = all_data_wgcna, mods = mod_name)

HubGeneNetworkPlot(
  all_data_wgcna[, all_data_wgcna$cell_type == "ALG11 Epithelial"],
  n_hubs = 3, # 用于可视化的hub gene
  n_other = 5, # 随机选取的gene
  edge_prop = 0.75, # 采样的边数
  mods = "all"
)

# export
pdf(
  file.path(outdir, "4 hub gene network.pdf"),
  width  = 12,
  height = 12
)
HubGeneNetworkPlot(
  all_data_wgcna[, all_data_wgcna$cell_type == "ALG11 Epithelial"],
  n_hubs = 5, # 用于可视化的hub gene
  n_other = 10, # 随机选取的gene
  edge_prop = 0.75, # 采样的边数
  mods = "all"
)
dev.off()
##  umap -----------

# 将基因降维成UMAP图，该函数会在Seurat_obj@misc$wgcna_name中生成一个module_umap矩阵
all_data_wgcna <- RunModuleUMAP(
  all_data_wgcna, 
  n_hubs = 10, # 用于UMAP嵌入的hub genes
  n_neighbors = 15, # UMAP参数
  min_dist = 0.1 # 两个点（基因）在UMAP空间中的最短距离
)
head(all_data_wgcna@misc$`ALG11 Epithelial`$module_umap, 2)


# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(all_data_wgcna)

p_uamp <-
  ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) + # 设置坐标轴
  geom_point(
    color = umap_df$color, # 按module对点（gene）着色
    size = umap_df$kME * 2 # 点的大小与其kME值相关
  ) + 
  umap_theme()

ModuleUMAPPlot(
  all_data_wgcna,
  edge.alpha = 0.25,
  sample_edges = T,
  edge_prop = 0.1, # 采样10%的边用于可视化
  label_hubs = 2, # 每个module用于展示的hub genes
  keep_grey_edges = F
)


## module trait correlation ---------------
# convert sex to factor
all_data_wgcna$group2 <- as.factor(all_data_wgcna$group)
all_data_wgcna$project <- as.factor(all_data_wgcna$project)




# 需要与模块进行关联的特征，该特征从Seurat_obj@meta.data抽取
cur_traits <- c('group2', 'project', 'nCount_RNA', 'nFeature_RNA')

# 将基因模块与特征进行关联
all_data_wgcna <- ModuleTraitCorrelation(
  all_data_wgcna,
  traits = cur_traits,  # 需要与模块进行关联的特征
  group.by = 'cell_type'  # 该参数表示以cell_type进行分组，以cell_type为单位计算每个cell_type的模块-特征相关性矩阵
)

# 查看相关性矩阵
mt_cor <- GetModuleTraitCorrelation(all_data_wgcna)

names(mt_cor)
# 'cor''pval''fdr'

names(mt_cor$cor)
# 'all_cells''INH''EX''OPC''ODC''ASC''MG'

head(mt_cor$cor$all_cells[, 1:5])

PlotModuleTraitCorrelation(
  all_data_wgcna,
  label = 'fdr',  # 可选pval、fdr作为显著性分数
  label_symbol = 'stars',  # 以*号作为显著性标记，numeric则显示数值
  text_size = 2,
  text_digits = 2,
  text_color = 'white',
  high_color = 'yellow',
  mid_color = 'black',
  low_color = 'purple',
  plot_max = 0.2,
  combine=TRUE
)





