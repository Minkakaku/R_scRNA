library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Hs.eg.db)

#all_data <- all_data2
all_data <- readRDS('R_OUTPUT/all_data.rds')

all_data$annotation %>% table()

major_marker_ls <- c(
  # 少突胶质谱系
  "OPC" = "Pdgfra",
  "OPC" = "Sox10",
  "OPC" = "Cspg4",
  "OPC" = "Enpp6",
  
  "Oligodendrocyte" = "Mbp",
  "Oligodendrocyte" = "Mog",
  "Oligodendrocyte" = "Mag",
  #"Oligodendrocyte" = "Sox10",
  
  # 血管壁细胞
  #"VSMC" = "Pdgfrb",
  "VSMC" = "Myh11",
  "VSMC" = "Mylk",
  "VSMC" = "Acta2",
  
  "Pericyte" = "Pdgfrb",
  "Pericyte" = "Rgs5",
  #"Pericyte" = "Cspg4",
  "Pericyte" = "Kcnj8",
  
  # 免疫细胞
  "Immune cell" = "Ptprc",
  "Immune cell" = "Cd68",
  "Immune cell" = "Cd44",
  
  "Border macrophage" = "Lyve1",
  "Border macrophage" = "Mrc1",
  "Border macrophage" = "Dab2",
  "Border macrophage" = "Clec10a",
  
  "Dural macrophage" = "Clec4n",
  #"Dural macrophage" = "Folr1",
  "Dural macrophage" = "Ccr2",
  "Dural macrophage" = "Lgals3",
  
  "Neutrophil" = "Cxcr2",
  
  "B cell" = "Cd19",
  
  "T cell" = "Nkg7",
  "T cell" = "Cd3e",
  "T cell" = "Cd4",
  "T cell" = "Cd8a",
  
  "Microglia" = "Aif1",
  "Microglia" = "Sall1",
  "Microglia" = "Tmem119",
  "Microglia" = "Cx3cr1",
  "Microglia" = "Iba1",
  
  # 星形胶质
  "Astrocyte" = "Aqp4",
  "Astrocyte" = "Slc1a3",
  "Astrocyte" = "Arxes2",
  "Astrocyte" = "Cxcl14",
  
  "Border astrocyte" = "Myoc",
  "Border astrocyte" = "Gfap",
  
  # 施万细胞
  "Schwann cell" = "Plp1",
  "Schwann cell" = "Sox2",
  "Schwann cell" = "Ngfr",
  #"Schwann cell" = "Sox10",
  "Schwann cell" = "Sema3d",
  
  # 血管内皮
  "Endothelial cell" = "Pecam1",
  "Endothelial cell" = "Plvap",
  "Endothelial cell" = "Ly6c1",
  "Endothelial cell" = "Cdh5",
  
  # 脑膜/硬膜/脉络丛
  "Ependymal cell" = "Cd24a",
  "Ependymal cell" = "Foxj1",
  "Ependymal cell" = "Calml4",
  "Ependymal cell" = "Folr1",
  "Ependymal cell" = "Ccdc153",
  "Ependymal cell" = "Tmem212",
  "Ependymal cell" = "Dynlrb2",
  "Ependymal cell" = "Rsph4a",
  
  "Choroid plexus" = "Kcnj13",
  "Choroid plexus" = "Flor1",
  "Choroid plexus" = "Clic6",
  "Choroid plexus" = "Kl",
  
  #"Leptomeningeal cell" = "Pdgfra",
  #"Leptomeningeal cell" = "Col1a1",
  "Leptomeningeal cell" = "Gjb6",
  "Leptomeningeal cell" = "Gjb2",
  "Leptomeningeal cell" = "Bmp7",
  "Leptomeningeal cell" = "Bmp6",
  "Leptomeningeal cell" = "Wnt5a",
  "Leptomeningeal cell" = "Lama1",
  "Leptomeningeal cell" = "Ptgds",
  
  #"Dural cell" = "Pdgfra",
  "Dural cell" = "Col1a1",
  "Dural cell" = "H19",
  "Dural cell" = "Matn4",
  "Dural cell" = "Clec11a",
  "Dural cell" = "Shisa3",
  
  # 神经元
  "Interneuron" = "Gad1",
  "Interneuron" = "Gad2",
  
  "Excitatory neuron" = "Slc17a7",
  "Excitatory neuron" = "Slc17a6"
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

Ependymal_cell_out_dir <- "/mnt/NexData_file/home/lxy/case/XE2601300820-01/stRNA_2/analysis/20260623_analysis/Ependymal_cell_subtype" 
dir.create(Ependymal_cell_out_dir)
Ependymal_cell <- subset(all_data, annotation == 'Ependymal cell')

all_data$annotation %>% table()
all_data$celltype %>% table()
## PCA --------
Ependymal_cell <- FindVariableFeatures(Ependymal_cell)

Ependymal_cell <- ScaleData(Ependymal_cell, vars.to.regress = 'nCount_RNA')
Ependymal_cell <- RunPCA(Ependymal_cell)

### pca inspector ------
pdf(file.path(Ependymal_cell_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Ependymal_cell, ndims = 50)
dev.off()


## dimension reduction --------
dims <- 18 # 18

#library(future)
#options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
#plan(multisession, workers = 10)

Ependymal_cell <- RunUMAP(Ependymal_cell, dims = 1:dims,return.model = TRUE)

pdf(file.path(Ependymal_cell_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Ependymal_cell, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

Ependymal_cell <- Ependymal_cell %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Ependymal_cell <- RunUMAP(
  Ependymal_cell, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(Ependymal_cell_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Ependymal_cell, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Ependymal_cell  <-
  FindNeighbors(
    Ependymal_cell ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Ependymal_cell <-
  FindClusters(
    Ependymal_cell,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(Ependymal_cell_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Ependymal_cell, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Ependymal_cell, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Ependymal_cell, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Ependymal_cell, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Ependymal_cell, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Ependymal_cell, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Ependymal_cell, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


### dotplot -------
Idents(Ependymal_cell)<-'RNA_snn_res.0.6'
Ependymal_cell_markers <-
  FindAllMarkers(
    JoinLayers(Ependymal_cell),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Ependymal_cell_markers, file = "Ependymal_cell_subtype/Ependymal_cell_all_markers.csv", row.names = TRUE)

Ependymal_cell_top5_markers <- 
  Ependymal_cell_markers %>% 
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

names(Ependymal_cell_top5_markers) <- sub('.*\\.','',names(Ependymal_cell_top5_markers))


# dotplot
# 如果已做其他细胞类型的亚群分析，这里就可以顺便画一下其他类型的细胞亚群基因，比如Fibro_major_marker_ls2，否则只用画T细胞
pdf(file.path(Ependymal_cell_out_dir, str_glue("dotplot dims{dims} major_marker_0.6.pdf")), width = 25, height = 6)
DotPlot(Ependymal_cell, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()+theme_90angle
#DotPlot(Ependymal_cell, group.by = "RNA_snn_res.0.4", features = PT_major_marker_ls)+RotatedAxis()+theme_90angle
#DotPlot(Ependymal_cell, group.by = "RNA_snn_res.0.4", features = PT_major_marker_ls_1)+RotatedAxis()+theme_90angle
DotPlot(Ependymal_cell, group.by = "RNA_snn_res.0.6", features = Ependymal_cell_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(Ependymal_cell, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()


saveRDS(Ependymal_cell ,file = "Ependymal_cell_subtype/Ependymal_cell_cluster.rds")
Ependymal_cell <- readRDS("Ependymal_cell_subtype/Ependymal_cell_cluster.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "Ependymal_cell_1", # 
  '1'  = "Ependymal_cell_2", #
  '2'  = "Ependymal_cell_3", # 
  '3'  = "Ependymal_cell_4", # 
  '4'  = "Ependymal_cell_5", # 
  '5'  = "Ependymal_cell_6",#  
  '6'  = "Astrocyte-like_Ependymal_cell", #  
  '7'  = "Astrocyte-like_Ependymal_cell", #  
  '8'  = "mix", # 
  '9'  = "Choroid_plexus-like_Ependymal_cell", # 
  '10'  = "Oligodendrocyte-like_Ependymal_cell", # 
  '11'  = "Neuron-like_Ependymal_cell", # 
  '12'  = "Leptomeningeal_cell-like_Ependymal_cell", # 
  '13'  = "Neuron-like_Ependymal_cell", # 
  '14'  = "mix"
)

Idents(Ependymal_cell) <- "RNA_snn_res.0.6"




Ependymal_cell <- RenameIdents(Ependymal_cell, cell_cluster)
Ependymal_cell$subtype <- Idents(Ependymal_cell)


pdf(file.path(Ependymal_cell_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Ependymal_cell, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Ependymal_cell_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 30)
DimPlot(Ependymal_cell, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Ependymal_cell_out_dir, str_glue("Dimplot dim{dims} annotation split by group.pdf")), height = 6, width = 15)
DimPlot(Ependymal_cell, group.by = "subtype",split.by = 'Group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

keep <- c(
  "Ependymal_cell_1", # 
  "Ependymal_cell_2", #
  "Ependymal_cell_3", # 
  "Ependymal_cell_4", # 
  "Ependymal_cell_5", # 
  "Ependymal_cell_6"
)

Ependymal_cell_clean <- subset(Ependymal_cell, subtype %in% keep)

Ependymal_cell_clean$subtype <- droplevels(Ependymal_cell_clean$subtype)
Idents(Ependymal_cell_clean) <- Ependymal_cell_clean$subtype

table(Ependymal_cell_clean$subtype)
saveRDS(Ependymal_cell_clean, "Ependymal_cell_subtype/Ependymal_cell_clean.rds")

Ependymal_cell <- Ependymal_cell_clean

Ependymal_cell <- RunUMAP(Ependymal_cell, dims = 1:{dims},n.neighbors = 50,return.model = TRUE)

pdf(file.path(Ependymal_cell_out_dir, str_glue("Dimplot dim{dims} Ependymal_cell annotation.pdf")), height = 6, width = 8)
DimPlot(Ependymal_cell, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Ependymal_cell_out_dir, str_glue("Dimplot dim{dims} Ependymal_cell annotation split.pdf")), height = 6, width = 30)
DimPlot(Ependymal_cell, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Ependymal_cell_out_dir, str_glue("Dimplot dim{dims} Ependymal_cell annotation split by group.pdf")), height = 6, width = 15)
DimPlot(Ependymal_cell, group.by = "subtype",split.by = 'Group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()


pdf(file.path(Ependymal_cell_out_dir, str_glue("Dimplot dim{dims} Ependymal_cell annotation harmony.pdf")), height = 6, width = 8)
DimPlot(Ependymal_cell, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()



# 改组名
group_map <- c(
  "CON1"  = "CON",
  "CON2"  = "CON",
  "CON3"  = "CON",
  "CRS4"  = "CRS",
  "CRS5"  = "CRS",
  "CRS6"  = "CRS"
)


Ependymal_cell@meta.data <- Ependymal_cell@meta.data %>%
  mutate(
    Group = recode(group, !!!group_map)
  )

table(Ependymal_cell$Group)
table(Ependymal_cell$group)
table(Ependymal_cell$subtype)
saveRDS(Ependymal_cell, "Ependymal_cell_subtype/Ependymal_cell_annotation.rds")
#nohup Rscript /mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/deg_gokegg/deg_enrich.R -i /mnt/NexData_file/home/lxy/case/XE2601300820-01/stRNA_2/analysis/20260623_analysis/Ependymal_cell_subtype/Ependymal_cell_annotation.rds --species mouse -g Group -c subtype --ident1 CRS --ident2 CON -o CRS_vs_CON_DEG &

Ependymal_cell_1 <- readRDS("Ependymal_cell_subtype/Ependymal_cell_cluster.rds")
cell_cluster <- c(
  '0'  = "Ependymal_cell_1", # 
  '1'  = "Ependymal_cell_2", #
  '2'  = "Ependymal_cell_3", # 
  '3'  = "Ependymal_cell_4", # 
  '4'  = "Ependymal_cell_5", # 
  '5'  = "Ependymal_cell_6",#  
  '6'  = "Astrocyte-like_Ependymal_cell", #  
  '7'  = "Astrocyte-like_Ependymal_cell", #  
  '8'  = "mix", # 
  '9'  = "Choroid_plexus-like_Ependymal_cell", # 
  '10'  = "Oligodendrocyte-like_Ependymal_cell", # 
  '11'  = "Neuron-like_Ependymal_cell", # 
  '12'  = "Leptomeningeal_cell-like_Ependymal_cell", # 
  '13'  = "Neuron-like_Ependymal_cell", # 
  '14'  = "mix"
)

Idents(Ependymal_cell_1) <- "RNA_snn_res.0.6"
Ependymal_cell_1 <- RenameIdents(Ependymal_cell_1, cell_cluster)
Ependymal_cell_1$subtype <- Idents(Ependymal_cell_1)

# 改组名
group_map <- c(
  "CON1"  = "CON",
  "CON2"  = "CON",
  "CON3"  = "CON",
  "CRS4"  = "CRS",
  "CRS5"  = "CRS",
  "CRS6"  = "CRS"
)


Ependymal_cell_1@meta.data <- Ependymal_cell_1@meta.data %>%
  mutate(
    Group = recode(group, !!!group_map)
  )

table(Ependymal_cell_1$Group)
table(Ependymal_cell_1$group)
table(Ependymal_cell_1$subtype)
saveRDS(Ependymal_cell_1, "Ependymal_cell_subtype/Ependymal_cell_allannotation.rds")


######下面合并到alldata的注释先不运行
all_data$annotation %>% table()
#all_data <- readRDS("output/data/all_data.rds")
# 1. 先把 all_data$annotation 复制到新列 T_annotation
all_data$annotation1 <- as.character(all_data$annotation)

# 2. 创建映射：细胞名 → 精细T亚型
cell_subtype <- Ependymal_cell_1$subtype
names(cell_subtype) <- colnames(Ependymal_cell_1)

# 3. 只把 "T cell" 替换成精细亚型，其他细胞保持原样
all_data$annotation1[colnames(Ependymal_cell_1)] <- as.character(cell_subtype)

# 查看结果（正确）
table(all_data$annotation1)
table(all_data$annotation)
table(Ependymal_cell_1$subtype)


saveRDS(all_data, "Ependymal_cell_subtype/all_data_T_annotation.rds")
saveRDS(Ependymal_cell_1 ,file = "Ependymal_cell_subtype/Ependymal_cell_annotation.rds")
##########不运行


Ependymal_cell <- readRDS("Ependymal_cell_subtype/Ependymal_cell_clean.rds")
#画一个注释完的点图，没有marker基因，先不运行
PT_major_marker_ls_2 <- c(
  "PT" = "LRP2",
  "PT" = "CUBN",
  "PTS1" = "SLC5A2",
  "PTS1" = "SLC16A1",
  "PTS2" = "SLC22A12",
  "PTS2" = "SLC22A8",
  "PTS3" = "ABCC3",
  "PTS3" = "SATB2",
  "PTS3" = "GPT"
)

pdf(file.path(Ependymal_cell_out_dir, str_glue("dotplot Ependymal_cell subtype major_marker.pdf")), width = 12, height = 6)
DotPlot(Ependymal_cell, group.by = "subtype", features = PT_major_marker_ls_2)+RotatedAxis()+theme_90angle
dev.off()


Ependymal_cell$subtype %>% table()
Ependymal_cell$sample %>% table()
Ependymal_cell$group %>% table()
## cell props barplot  ----------------------------------------------------------
Ependymal_cell_meta <- Ependymal_cell[[]]

Ependymal_cell_meta$sample <- factor(Ependymal_cell_meta$sample,levels=c(Ependymal_cell_meta$sample %>% unique()))

p_col_proj <-
  Ependymal_cell_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  scale_fill_brewer(palette = "Set2", name = "Ependymal_cell Subtype") +
  #ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(Ependymal_cell_out_dir, "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  Ependymal_cell_meta %>%
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
  filename = file.path(Ependymal_cell_out_dir, "cell_freq_group_1.pdf"),
  height = 6,
  width  = 4
)

p_col_proj <-
  Ependymal_cell_meta %>%
  group_by(Group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = Group, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path(Ependymal_cell_out_dir, "cell_freq_group_2.pdf"),
  height = 6,
  width  = 4
)


#改图
# 提取meta
Ependymal_cell_meta <- Ependymal_cell[[]]

# 对应配色（保持你指定色值不变）
group_col <- c(
  "CON"  = "#2E86AB",
  "CRS"    = "#A23B72"
)

# ========== 2. 计算样本水平百分比（用于箱线图 & 均值堆叠图） ==========
# 每个样本内各细胞亚群百分比
data_long <- Ependymal_cell_meta %>%
  group_by(sample, Group, subtype) %>%
  summarise(Count = n(), .groups = "drop_last") %>%
  group_by(sample) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

# 保存原始百分比数据
write.csv(data_long, 
         file.path("Ependymal_cell_subtype", "Ependymal_cell_sample_percentage.csv"), 
         row.names = FALSE)


# ==================== 按细胞类型分组的箱线图 ====================
comp_list <- list(c("CON", "CRS"))
p <- ggplot(data_long, aes(x = Group, y = Percentage, fill = Group)) +
  geom_boxplot(width = 0.7, alpha = 0.8, outlier.shape = 19, outlier.size = 1.5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7, aes(color = Group)) +
  facet_wrap(~subtype, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_col) +
  scale_color_manual(values = group_col) +
  labs(title = "Ependymal_cell - Cell Type Percentages",
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
ggsave(file.path("Ependymal_cell_subtype","Ependymal_cell_Boxplot.pdf"), p, width = 11, height = 8, dpi = 300)
print(p)

names(Ependymal_cell)

# 映射到空间切片
Ependymal_cell_spatial_out_dir <- "/mnt/NexData_file/home/lxy/case/XE2601300820-01/stRNA_2/analysis/20260623_analysis/Ependymal_cell_subtype/spatial" 
dir.create(Ependymal_cell_spatial_out_dir)

cells_boundary_ls <- readRDS("R_OUTPUT/data_for_analysis/cells_boundary_ls.rds")

head(cells_boundary_ls)
tail(cells_boundary_ls)


Ependymal_cell_meta <- 
  Ependymal_cell[[]] %>% 
  rownames_to_column('cell_id') 
head(Ependymal_cell_meta)
tail(Ependymal_cell_meta)


cn_sp_p_ls <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){ # browser()
    anno_data   <- Ependymal_cell_meta[c("sample", "Group")] %>% distinct()
    group_info  <- anno_data$Group[anno_data$sample == sample]
    
    dat <-
      cells_boundary %>%
      right_join(Ependymal_cell_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() + 
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL, title = paste(sample, group_info, sep = "_") )
    p
  })

cn_sp_p_ls %>% 
  iwalk(function(p, name){
    if(all(is.na(p$data$vertex_x))) return(NULL)
    if(any(is.na(p$data$vertex_x))) warning("vertex_x has NAs ...", .immediate = T)
    pdf(
      file.path(Ependymal_cell_spatial_out_dir, str_glue("Subtype spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })



