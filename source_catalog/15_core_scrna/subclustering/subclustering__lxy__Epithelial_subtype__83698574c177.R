library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Mm.eg.db)


all_data <- readRDS('/mnt/NexData_file/home/lxy/case/XE2603162231-02/20260616_XE2603162231-02/R_OUTPUT/all_data.rds')

#all_data <- readRDS(file.path("R_OUTPUT", "all_data.rds"))

# markers: https://cloud.tencent.com/developer/article/2115987
major_marker_ls <- c(
  "B cell" = "Cd19",
  "B cell" = "Cd79a",
  "B cell" = "Ms4a1",
  
  "Plasma" = "Mzb1",
  "Plasma" = "Jchain",
  
  "T cell" = "Cd3d",
  "T cell" = "Cd3e",
  "T cell" = "Cd4",
  "T cell" = "Foxp3",
  "T cell" = "Cd8a",
  "T cell" = "Gzmb",
  
  "NK" = "Ncam1",
  
  "Proliferating" = "Mki67",
  "Proliferating" = "Top2a",
  
  "Myeloid" = "Itgam",
  
  "Macrophage" = "Cd68",
  "Macrophage" = "Cd163",
  
  "Mast cell" = "Ms4a2",
  "Mast cell" = "Kit",
  "Mast cell" = "Hdc",
  "Mast cell" = "Tpsab1",
  
  "Monocytes" = "Cd14",
  
  "Dendritic" = "Cd83",
  "Dendritic" = "Clec9a",
  "Dendritic" = "Ido1",
  
  "Neutrophi" = "Csf3r",
  "Neutrophi" = "Cxcr2",
  
  "Endothelial" = "Pecam1",
  "Endothelial" = "Vwf",
  
  "Epithelial" = "Epcam",
  "Epithelial" = "Krt20",
  
  "Enterocyte" = "Slc26a3",
  "Enterocyte" = "Ca1",
  "Enterocyte" = "Guca2a",
  "Enterocyte" = "Guca2b",
  "Enterocyte" = "Fabp2",
  "Enterocyte" = "Rnase1",
  "Enterocyte" = "Spink1",
  
  "Goblet" = "Muc2",
  "Goblet" = "Tff3",
  
  "Enteroendocrine" = "Chga",
  "Enteroendocrine" = "Chgb",
  
  "Neuron" = "Snap25",
  "Neuron" = "Tubb3",
  "Neuron" = "Map2",
  
  "Fibroblasts" = "Pdgfra",
  "Fibroblasts" = "Col1a1",
  
  "Smooth muscle cells" = "Pdgfrb",
  "Smooth muscle cells" = "Pln",
  "Smooth muscle cells" = "Myh11",
  "Smooth muscle cells" = "Tagln",
  
  "Mesothelial" = "Msln",
  
  "Paneth-like cell" = "Defa5",
  "Paneth-like cell" = "Defa6",
  "Paneth-like cell" = "Reg3a",
  
  "Tuft cell" = "Trpm5",
  "Tuft cell" = "Sh2d6",
  "Tuft cell" = "Pou2f3",
  
  "Colonocyte" = "Ceacam1"
  #"Colonocyte" = "Guca2a"
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )



# Epithelial -------------------------------------------------------------------
all_data$annotation %>% table()
all_data$annotation <- str_replace(all_data$annotation, 'Enterocyte', 'Epithelial')
all_data$annotation <- str_replace(all_data$annotation, 'Goblet', 'Epithelial')
all_data$annotation %>% table()


Epithelial_out_dir <- "Epithelial_subtype" 
dir.create(Epithelial_out_dir)
Epithelial <- subset(all_data, annotation == 'Epithelial')

all_data$annotation %>% table()
all_data$celltype %>% table()
## PCA --------
Epithelial <- FindVariableFeatures(Epithelial)

Epithelial <- ScaleData(Epithelial, vars.to.regress = 'nCount_RNA')
Epithelial <- RunPCA(Epithelial)

### pca inspector ------
pdf(file.path(Epithelial_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Epithelial, ndims = 50)
dev.off()


## dimension reduction --------
dims <- 17 # 18

#library(future)
#options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
#plan(multisession, workers = 10)

Epithelial <- RunUMAP(Epithelial, dims = 1:dims,return.model = TRUE)

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

Epithelial <- Epithelial %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Epithelial <- RunUMAP(
  Epithelial, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Epithelial  <-
  FindNeighbors(
    Epithelial ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Epithelial <-
  FindClusters(
    Epithelial,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Epithelial, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Epithelial, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Epithelial, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

Idents(Epithelial) <- "RNA_snn_res.0.8"

### dotplot -------
Epithelial_major_marker_ls <- c(
  "Epithelial" = "Epcam",
  "Epithelial" = "Krt20",
  
  'Enterocyte' = 'Slc26a3', 
  'Enterocyte' = 'Ca1', 
  'Enterocyte' = 'Guca2a', 
  'Enterocyte' = 'Guca2b', 
  'Enterocyte' = 'Fabp2', 
  'Enterocyte' = 'Rnase1', 
  'Enterocyte' = 'Spink1', 
  
  'Goblet' = 'Muc2', 
  'Goblet' = 'Tff3', 
  
  'Enteroendocrine' = 'Chga', 
  "Enteroendocrine" = "Chgb",
  
  'Paneth cell' = 'Defa5', 
  'Paneth cell' = 'Defa6', 
  'Paneth cell' = 'Reg3a',
  
  'Mesothelial' = 'Msln', 
  
  'Neuron' = 'Tubb3', 
  'Neuron' = 'Snap25',
  'Neuron' =  'Map2',
  
  'Tuft cell' = 'Trpm5',
  'Tuft cell' = 'Sh2d6',
  'Tuft cell' = 'Pou2f3',
  
  #'Colonocyte' = 'Ceacam1',
  #'Colonocyte' = 'Guca2a'
  
  'Proliferating cell' = 'Mki67',
  'Proliferating cell' = 'Top2a',
  
  #'Intestinal stem cell' = 'Lgr5', # 也是癌干细胞,Wnt/β-catenin通路基因,ref：https://pmc.ncbi.nlm.nih.gov/articles/PMC8941949/和https://pmc.ncbi.nlm.nih.gov/articles/PMC7616986
  #'Intestinal stem cell' = 'Olfm4', #也是癌干细胞
  'Intestinal stem cell' = 'Ascl2',
  'Intestinal stem cell' = 'Axin2',
  'Intestinal stem cell' = 'Rnf43', #待定
  'Intestinal stem cell' = 'Ephb2',
  #'Intestinal stem cell' = 'Sox9',
  
  'Cancer stem cell' = 'Lgr5',
  'Cancer stem cell' = 'Olfm4',
  'Cancer stem cell' = 'Sox9',
  
  'Cancer cell' = 'Muc1', #ERO1A-positive tumor epithelial cells in colorectal cancer progression: a multi-omics perspective
  'Cancer cell' = 'Mmp7',
  'Cancer cell' = 'Ero1a'
  
)


Idents(Epithelial)<-'RNA_snn_res.0.8'
Epithelial_markers <-
  FindAllMarkers(
    JoinLayers(Epithelial),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Epithelial_markers, file = "Epithelial_subtype/Epithelial_all_markers.csv", row.names = TRUE)

Epithelial_top5_markers <- 
  Epithelial_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 5) %>%
  nest() %>%
  deframe() %>% 
  map(pull, 'gene') %>%
  imap(function(genes, name){
    re <- genes %>%
      map(function(gene){
        # browser()
        names(gene) = str_glue("cluster_{name}")
        gene
      })
    re
  }) %>%
  unlist() %>% 
  .[!duplicated(.)]

names(Epithelial_top5_markers) <- sub('.*\\.','',names(Epithelial_top5_markers))


# dotplot
# 如果已做其他细胞类型的亚群分析，这里就可以顺便画一下其他类型的细胞亚群基因，比如Fibro_major_marker_ls2，否则只用画T细胞
pdf(file.path(Epithelial_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = Epithelial_major_marker_ls)+RotatedAxis()+theme_90angle
#DotPlot(Epithelial, group.by = "RNA_snn_res.0.4", features = PT_major_marker_ls_1)+RotatedAxis()+theme_90angle
DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = Epithelial_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()

pdf(file.path(Epithelial_out_dir, str_glue("dotplot dims{dims} major_marker_0.8.pdf")), width = 25, height = 6)
DotPlot(Epithelial, group.by = "RNA_snn_res.0.8", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Epithelial, group.by = "RNA_snn_res.0.8", features = Epithelial_major_marker_ls)+RotatedAxis()+theme_90angle
#DotPlot(Epithelial, group.by = "RNA_snn_res.0.4", features = PT_major_marker_ls_1)+RotatedAxis()+theme_90angle
DotPlot(Epithelial, group.by = "RNA_snn_res.0.8", features = Epithelial_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()


saveRDS(Epithelial ,file = "Epithelial_subtype/Epithelial_cluster.rds")
Epithelial <- readRDS("Epithelial_subtype/Epithelial_cluster.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "Cancer stem cell", # 
  '1'  = "Cancer stem cell", #
  '2'  = "Goblet", # 
  '3'  = "Cancer stem cell", # 
  '4'  = "Cancer stem cell", # 
  '5'  = "Enterocyte",#  
  '6'  = "Goblet", #  
  '7'  = "Cancer stem cell", #  
  '8'  = "Cancer stem cell", #  
  '9'  = "Enterocyte",  
  '10'  = "Enterocyte", 
  '11'  = "Cancer stem cell", #
  '12'  = "Paneth cell", # 
  '13'  = "Enterocyte", # 
  '14'  = "Cancer stem cell", # 
  '15'  = "Fibroblast-like_Epithelial",#  
  '16'  = "Goblet", 
  '17'  = "Cancer stem cell",  
  '18'  = "Fibroblast-like_Epithelial",  
  '19'  = "Enterocyte",  
  '20'  = "Monocyte-like_Epithelial"
)

Idents(Epithelial) <- "RNA_snn_res.0.8"




Epithelial <- RenameIdents(Epithelial, cell_cluster)
Epithelial$subtype <- Idents(Epithelial)


pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 15)
DimPlot(Epithelial, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


keep <- c(
  "Enterocyte",
  "Goblet",
  "Paneth cell",
  "Cancer stem cell"
)

Epithelial_clean <- subset(Epithelial, subtype %in% keep)

Epithelial_clean$subtype <- droplevels(Epithelial_clean$subtype)
Idents(Epithelial_clean) <- Epithelial_clean$subtype

table(Epithelial_clean$subtype)
saveRDS(Epithelial_clean, "Epithelial_subtype/Epithelial_clean.rds")

Epithelial <- Epithelial_clean

#Epithelial <- RunUMAP(Epithelial, dims = 1:{dims},n.neighbors = 50,return.model = TRUE)

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} Epithelial annotation.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} Epithelial annotation split.pdf")), height = 6, width = 15)
DimPlot(Epithelial, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

Epithelial <- RunUMAP(Epithelial, dims = 1:{dims},n.neighbors = 50,return.model = TRUE)

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} Epithelial annotation_1.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} Epithelial annotation split_1.pdf")), height = 6, width = 15)
DimPlot(Epithelial, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()



Epithelial_1 <- readRDS("Epithelial_subtype/Epithelial_cluster.rds")
cell_cluster <- c(
  '0'  = "Cancer stem cell", # 
  '1'  = "Cancer stem cell", #
  '2'  = "Goblet", # 
  '3'  = "Cancer stem cell", # 
  '4'  = "Cancer stem cell", # 
  '5'  = "Enterocyte",#  
  '6'  = "Goblet", #  
  '7'  = "Cancer stem cell", #  
  '8'  = "Cancer stem cell", #  
  '9'  = "Enterocyte",  
  '10'  = "Enterocyte", 
  '11'  = "Cancer stem cell", #
  '12'  = "Paneth cell", # 
  '13'  = "Enterocyte", # 
  '14'  = "Cancer stem cell", # 
  '15'  = "Fibroblast-like_Epithelial",#  
  '16'  = "Goblet", 
  '17'  = "Cancer stem cell",  
  '18'  = "Fibroblast-like_Epithelial",  
  '19'  = "Enterocyte",  
  '20'  = "Monocyte-like_Epithelial"
)

Idents(Epithelial_1) <- "RNA_snn_res.0.8"

Epithelial_1 <- RenameIdents(Epithelial_1, cell_cluster)
Epithelial_1$subtype <- Idents(Epithelial_1)

saveRDS(Epithelial_1, "Epithelial_subtype/Epithelial_allannotation.rds")

######下面合并到alldata的注释
all_data$annotation %>% table()
#all_data <- readRDS("output/data/all_data.rds")
# 1. 先把 all_data$annotation 复制到新列 T_annotation
all_data$annotation1 <- as.character(all_data$annotation)

# 2. 创建映射：细胞名 → 精细T亚型
cell_subtype <- Epithelial_1$subtype
names(cell_subtype) <- colnames(Epithelial_1)

# 3. 只把 "T cell" 替换成精细亚型，其他细胞保持原样
all_data$annotation1[colnames(Epithelial_1)] <- as.character(cell_subtype)

# 查看结果（正确）
table(all_data$annotation1)
table(all_data$annotation)
table(Epithelial_1$subtype)
saveRDS(all_data, "Epithelial_subtype/all_data_Epithelial_annotation.rds")




Epithelial <- readRDS("Epithelial_subtype/Epithelial_clean.rds")
#画一个注释完的点图
Epithelial_major_marker_ls_2 <- c(
  "Epithelial" = "Epcam",
  "Epithelial" = "Krt20",
  
  'Enterocyte' = 'Slc26a3', 
  'Enterocyte' = 'Spink1', 
  
  'Goblet' = 'Muc2', 
  
  'Paneth cell' = 'Defa5', 

  
  'Cancer stem cell' = 'Mki67',
  'Cancer stem cell' = 'Top2a',
  'Cancer stem cell' = 'Lgr5',
  'Cancer stem cell' = 'Sox9',
  'Cancer stem cell' = 'Mmp7',
  'Cancer stem cell' = 'Ascl2',
  'Cancer stem cell' = 'Axin2',
  'Cancer stem cell' = 'Ephb2'
)


pdf(file.path(Epithelial_out_dir, str_glue("dotplot Epithelial subtype major_marker.pdf")), width = 12, height = 6)
DotPlot(Epithelial, group.by = "subtype", features = Epithelial_major_marker_ls_2)+RotatedAxis()+theme_90angle
dev.off()


Epithelial$subtype %>% table()
Epithelial$sample %>% table()
Epithelial$group %>% table()
## cell props barplot  ----------------------------------------------------------
Epithelial_meta <- Epithelial[[]]

Epithelial_meta$sample <- factor(Epithelial_meta$sample,levels=c(Epithelial_meta$sample %>% unique()))

p_col_proj <-
  Epithelial_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  scale_fill_brewer(palette = "Set2", name = "Epithelial Subtype") +
  #ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(Epithelial_out_dir, "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  Epithelial_meta %>%
  group_by(group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path(Epithelial_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)

#改图
# 提取meta
Epithelial_meta <- Epithelial[[]]

# 新分组显示顺序（左→右：Normal, Mild, Severe）
group_levels <- c("CON", "HDFMT", "CRCFMT")
# 对应配色（保持你指定色值不变）
group_col <- c(
  "CON"  = "#2E86AB",
  "HDFMT"    = "#A23B72",
  "CRCFMT"  = "#F18F01"
)


# 确保meta里的group是因子，顺序和颜色对应
#Epithelial_meta$Group <- recode(Epithelial_meta$group, !!!group_map)
Epithelial_meta$group <- factor(Epithelial_meta$group, levels = group_levels)


# ========== 2. 计算样本水平百分比（用于箱线图 & 均值堆叠图） ==========
# 每个样本内各细胞亚群百分比
data_long <- Epithelial_meta %>%
  group_by(sample, group, subtype) %>%
  summarise(Count = n(), .groups = "drop_last") %>%
  group_by(sample) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

# 保存原始百分比数据
write.csv(data_long, 
          file.path("Epithelial_subtype", "Epithelial_sample_percentage.csv"), 
          row.names = FALSE)

# 按Group + CellType 计算均值（用于堆叠柱状图）
category_avg <- data_long %>%
  group_by(group, subtype) %>%
  summarise(MeanPercentage = mean(Percentage), .groups = "drop")

# 保存均值统计数据
write.csv(category_avg, 
          file.path("Epithelial_subtype", "Epithelial_group_mean_percentage.csv"), 
          row.names = FALSE)

# ========== 堆叠柱状图 ==========
p_category_stack <- ggplot(category_avg, 
                           aes(x = group, y = MeanPercentage, fill = subtype)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_brewer(palette = "Set2", name = "Epithelial Subtype") +
  labs(title = "Epithelial Subtype Composition by Group",
       x = "Group", y = "Percentage (%)") +
  coord_cartesian(ylim = c(0, 100)) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

ggsave(file.path("Epithelial_subtype", "Epithelial_StackedBar.pdf"), 
       p_category_stack, width = 8, height = 6, dpi = 300)
print(p_category_stack)


# ==================== 按细胞类型分组的箱线图 ====================
comp_list <- list(c("CON","HDFMT"), c("CON","CRCFMT"), c("HDFMT","CRCFMT"))
p <- ggplot(data_long, aes(x = group, y = Percentage, fill = group)) +
  geom_boxplot(width = 0.7, alpha = 0.8, outlier.shape = 19, outlier.size = 1.5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7, aes(color = group)) +
  facet_wrap(~subtype, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = group_col) +
  scale_color_manual(values = group_col) +
  labs(title = "Epithelial - Cell Type Percentages",
       x = "Group", y = "Percentage (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))+
  stat_compare_means(
    comparisons = comp_list,
    method = "wilcox.test",
    label = "p.format",
    size = 3
  )

# 保存+展示
ggsave(file.path("Epithelial_subtype","Epithelial_Boxplot.pdf"), p, width = 10, height = 8, dpi = 300)
print(p)


# 映射到空间切片
Epithelial_spatial_out_dir <- "/mnt/NexData_file/home/lxy/case/XE2603162231-02/20260616_XE2603162231-02/Epithelial_subtype/spatial" 
dir.create(Epithelial_spatial_out_dir)

cells_boundary_ls <- readRDS("R_OUTPUT/data_for_analysis/cells_boundary_ls.rds")

head(cells_boundary_ls)
tail(cells_boundary_ls)


all_data_meta <- 
  all_data[[]] %>% 
  rownames_to_column('cell_id') 
head(all_data_meta)
tail(all_data_meta)


cn_sp_p_ls <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){ # browser()
    anno_data   <- all_data_meta[c("sample", "group")] %>% distinct()
    group_info  <- anno_data$Group[anno_data$sample == sample]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    dat <- dat %>%
      mutate(plot_group = ifelse(annotation1 %in% keep, annotation1, "Others"))
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation1, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() + 
      # 自定义颜色：为keep列表中的细胞分配颜色，"Others"设为灰色
      scale_fill_manual(
        values = c(
          "Enterocyte" = "#1f77b4",    # 蓝色
          "Goblet" = "#ff7f0e",         # 橙色
          "Paneth cell" = "#2ca02c",    # 绿色
          "Cancer stem cell" = "#d62728",# 红色
          "Others" = "#dcdcdc"          # 浅灰色
        ),
        limits = c(keep, "Others")) +
      theme_void() +
      labs(fill = NULL, title = paste(sample, group_info, sep = "_") )
    p
  })

cn_sp_p_ls %>% 
  iwalk(function(p, name){
    if(all(is.na(p$data$vertex_x))) return(NULL)
    if(any(is.na(p$data$vertex_x))) warning("vertex_x has NAs ...", .immediate = T)
    pdf(
      file.path(Epithelial_spatial_out_dir, str_glue("Subtype spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


