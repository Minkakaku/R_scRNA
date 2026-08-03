library(Seurat)
library(tidyverse)
all_data <- readRDS('/mnt/NexData_file/home/lzl/case/XE2510316140-03/output/data/adata_filtered_add_celltype.rds')
unique(all_data$celltype)
all_data[['RNA']] <- as(all_data[['RNA']], "Assay5")
all_data@assays$RNA$data <- NULL
all_data[['RNA']] <- split(all_data[['RNA']], all_data$project, layers = "counts")
all_data$annotation <- as.character(all_data$annotation)
all_data$annotation[all_data$annotation == "CD4⁺ T cell"] <- "CD4+ T cell"
all_data$annotation <- as.factor(all_data$annotation)
# normalize
future::plan("sequential")
all_data <- all_data %>% NormalizeData()
Idents(all_data) <- "annotation"
saveRDS(all_data ,file = "/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/data/all_data.rds")

all_data <- readRDS("/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/data/all_data.rds")
all_data@meta.data <- all_data@meta.data  %>%
  mutate(
    group_2 = case_when(
      grepl("S", sample) ~ "S",
      grepl("P", sample) ~ "P",
      grepl("R", sample) ~ "R",
      TRUE ~ NA_character_  # 都不匹配时设为 NA
    )
  )

all_data@meta.data <- all_data@meta.data  %>%
  mutate(
    group_3 = case_when(
      grepl("A", sample)&grepl("S", sample) ~ "A_S",
      grepl("A", sample)&grepl("P", sample) ~ "A_P",
      grepl("A", sample)&grepl("R", sample) ~ "A_R",
      grepl("C", sample)&grepl("S", sample) ~ "C_S",
      grepl("C", sample)&grepl("P", sample) ~ "C_P",
      TRUE ~ NA_character_  # 都不匹配时设为 NA
    )
  )
table(all_data$sample)
saveRDS(all_data ,file = "/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/data/all_data.rds")
all_data <- subset(all_data,!sample%in%c())
####
# Epi -------------------------------------------------------------------

Epi_out_dir <- "Epi_subtype" 
dir.create(Epi_out_dir)

Epi <- subset(all_data, annotation%in% c('Epithelial',"Columnar cell","Basal cell"))


## PCA --------
Epi <- FindVariableFeatures(Epi)

options(future.globals.maxSize = 50 * 1024^3) # 把并行全局变量上限改成 50GB（足够用）
Epi <- ScaleData(Epi, vars.to.regress = 'nCount_RNA')
Epi <- RunPCA(Epi)

### pca inspector ------
pdf(file.path(Epi_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Epi, ndims = 50)
dev.off()


Epi$annotation %>% table()
#Epi$celltype %>% table()
## dimension reduction --------
dims <- 20 # 18

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 20 * 1024 ^ 3) # 20Gb

Epi <- RunUMAP(Epi, dims = 1:dims,return.model = TRUE)

pdf(file.path(Epi_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Epi, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

Epi <- Epi %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Epi <- RunUMAP(
  Epi, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(Epi_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Epi, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Epi  <-
  FindNeighbors(
    Epi ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Epi <-
  FindClusters(
    Epi,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(Epi_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Epi, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epi, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epi, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epi, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Epi, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Epi, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Epi, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()
###########
Idents(Epi) <- "RNA_snn_res.0.4"
###marker_ls
Epi_major_marker_ls <- c(
  "AT1 cell" = "AGER",
  "AT1 cell" = "SEMA3B",
  "AT1 cell" = "PDPN",
  "AT2 cell" = "ABCA3",
  "AT2 cell" = "SFTPA1",
  "AT2 cell" = "SFTPB",
  "AT2 cell" = "SFTPC",
  "AT2 cell" = "HHIP",
  "Basal cell" = "KRT5",
  "Basal cell" = "NGFR",
  "Basal cell" = "TP63",
  "Basal cell" = "KRT15",
  "Deuterosomal cell" = "CDC20B",
  "Goblet cell" = "MUC5AC",
  "Goblet cell" = "CEACAM5",
  "Pulmonary neuroendocrine cell" = "GRP",
  "Pulmonary neuroendocrine cell" = "CALCA",
  "Tuft cell" = "POU2F3",
  "Tuft cell" = "ASCL2",
  "Lonocyte" = "ASCL3",
  "Lonocyte" = "BSND",
  "Lonocyte" = "FOXI1",
  "Hillock cell" = "KRT14",
  "Hillock cell" = "KRT6A",
  "Columnar cell" = "FOXJ1",
  "Columnar cell" = "RSPH1",
  "Columnar cell" = "CDHR3",
  "Columnar cell" = "CDHR4",
  "Mucous cell" = "MUC5B",
  "Mucous cell" = "SPDEF",
  "Mucous cell" = "BPIFB2",
  "Serous cell" = "LYZ",
  "Serous cell" = "LTF",
  "Secretory cell" = "SCGB1A1",
  "Secretory cell" = "SCGB3A2",
  "Ciliated cell" = "AKAP14",
  "Ciliated cell" = "ALDH3B1",
  "Ciliated cell" = "ANKRD66",
  "Ciliated cell" = "C11ORF97",
  "Ciliated cell" = "DNAI1",
  "Ciliated cell" = "DLEC1",
  "Club cell" = "SCGB3A1",
  "Club cell" = "PIGR",
  "Lymphatic endothelial cell" = "LYVE1",
  "Lymphatic endothelial cell" = "CCL21A",
  "Lymphatic endothelial cell" = "CCL21B"
)
existing_genes <- rownames(all_data)

# 2. 过滤 marker 列表，只保留存在的基因
Epi_major_marker_ls <- Epi_major_marker_ls[
  Epi_major_marker_ls %in% existing_genes
]


Epi_markers <-
  FindAllMarkers(
    JoinLayers(Epi),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Epi_markers, file = "Epi_subtype/Epi_all_markers.csv", row.names = TRUE)


Epi_top5_markers <- 
  Epi_markers %>% 
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

names(Epi_top5_markers) <- sub('.*\\.','',names(Epi_top5_markers))


# dotplot
# 如果已做其他细胞类型的亚群分析，这里就可以顺便画一下其他类型的细胞亚群基因，比如Fibro_major_marker_ls2，否则只用画T细胞
pdf(file.path(Epi_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(Epi, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Epi, group.by = "RNA_snn_res.0.4", features = Epi_major_marker_ls)+RotatedAxis()+theme_90angle
#DotPlot(Macrophage, group.by = "RNA_snn_res.0.4", features = Macrophage_major_marker_ls_1)+RotatedAxis()+theme_90angle
#DotPlot(Macrophage, group.by = "RNA_snn_res.0.6", features = Macrophage_markers)+RotatedAxis()
DotPlot(Epi, group.by = "RNA_snn_res.0.4", features = Epi_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(Macrophage, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()

## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "AT2 cell", # 
  '1'  = "Basal cell", #
  '2'  = "Fibroblast-like", #
  '3'  = "Unknown", # 
  '4'  = "Mucous cell", # 
  '5'  = "Proliferating basal cell", # 
  '6'  = "AT1 cell", # 
  '7'  = "Columnar cell", # 
  '8'  = "Mac-like",# 
  '9'  = "Basal cell", # 
  '10' = "Unknown", # 
  '11' = "Interferon-stimulated epithelial"
)

Idents(Epi) <- "RNA_snn_res.0.4"



Epi <- RenameIdents(Epi, cell_cluster)
Epi$subtype <- Idents(Epi)


pdf(file.path(Epi_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Epi, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Epi_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 15)
DimPlot(Epi, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

saveRDS(Epi, file.path(Epi_out_dir, 'Epi_subtype.rds'))
###
Epi_clean <- Epi[, !Epi$subtype %in% c('Fibroblast-like',"Unknown","Mac-like")]
Epi_clean <- RunUMAP(Epi_clean, dims = 1:dims, return.model = TRUE)


pdf(file.path(Epi_out_dir,str_glue("Dimplot dim{dims} annotation - pure.pdf")), height = 6, width = 8)
DimPlot(Epi_clean, group.by = "subtype", label = TRUE, repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Epi_out_dir, str_glue("Dimplot dim{dims} annotation split by group - pure.pdf")), height = 7, width = 18)
DimPlot(Epi_clean, group.by = "subtype",split.by = 'group', label = TRUE, repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Epi_out_dir, str_glue("Dimplot dim{dims} annotation split by group_2 - pure.pdf")), height = 7, width = 24)
DimPlot(Epi_clean, group.by = "subtype",split.by = 'group_2', label = TRUE, repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Epi_out_dir, str_glue("Dimplot dim{dims} annotation split by group_3 - pure.pdf")), height = 7, width = 30)
DimPlot(Epi_clean, group.by = "subtype",split.by = 'group_3', label = TRUE, repel = TRUE) + NoLegend()
dev.off()


saveRDS(Epi_clean, file.path(Epi_out_dir, 'Epi_subtype_clean.rds'))
#######
library(msigdbr)
# 获取人类 Hallmark 基因集
hallmark <- msigdbr(species = "Homo sapiens", category = "H")

# 提取目标通路
target_pathways <- c(
  "HALLMARK_HYPOXIA",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_GLYCOLYSIS"
)

# 筛选并整理
pathway_genes <- hallmark %>%
  filter(gs_name %in% target_pathways) %>%
  select(pathway = gs_name, gene = gene_symbol) %>%
  arrange(pathway, gene)
pathway_genes %>% count(pathway)

pathway_list <- pathway_genes %>%
  group_by(pathway) %>%
  summarise(genes = list(gene)) %>%
  deframe() 
genes_in_expr <- rownames(Epi_clean)
pathway_list <- lapply(pathway_list, function(genes) {
  intersect(genes, genes_in_expr)
})
###
Epi_clean <- JoinLayers(Epi_clean)
Epi_clean <- AddModuleScore(
  object = Epi_clean,
  features = list(
    Glycolysis = pathway_list[["HALLMARK_GLYCOLYSIS"]],
    Hypoxia = pathway_list[["HALLMARK_HYPOXIA"]],
    OXPHOS = pathway_list[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]]
  ),
  name = c("Glycolysis_", "Hypoxia_", "OXPHOS_")
)

pdf(file.path(Epi_out_dir, "score_vlnplot.pdf"), width = 8, height = 4)
VlnPlot(
  object = Epi_clean,
  features = c("Glycolysis_1", "Hypoxia_2", "OXPHOS_3"),
  group.by = "subtype",
  pt.size = 0,           # 不显示散点
  ncol = 3
) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
dev.off()

pdf(file.path(Epi_out_dir,str_glue("Feature tumor marker.pdf")), height = 6, width = 30)
FeaturePlot(Epi_clean, features = "MUC1",split.by = "group_3", label = FALSE, repel = TRUE)
FeaturePlot(Epi_clean, features = "TP63",split.by = "group_3", label = FALSE, repel = TRUE)
FeaturePlot(Epi_clean, features = "SOX2",split.by = "group_3", label = FALSE, repel = TRUE)
dev.off()


pdf(file.path(Epi_out_dir,str_glue("Feature score.pdf")), height = 6, width = 30)
FeaturePlot(Epi_clean, features = "Glycolysis_1",split.by = "group_3", label = FALSE, repel = TRUE)
FeaturePlot(Epi_clean, features = "Hypoxia_2",split.by = "group_3", label = FALSE, repel = TRUE)
FeaturePlot(Epi_clean, features = "OXPHOS_3",split.by = "group_3", label = FALSE, repel = TRUE)
dev.off()


#####
Epi_meta <-Epi_clean[[]]
Epi_meta$sample <- factor(Epi_meta$sample,levels=c(Epi_meta$sample %>% unique()))
Epi_meta$subtype <- factor(Epi_meta$subtype,levels=c(Epi_meta$subtype %>% unique()))
Epi_meta <- Epi_meta %>%
  mutate(
    group_2 = case_when(
      grepl("S", sample) ~ "S",
      grepl("P", sample) ~ "P",
      grepl("R", sample) ~ "R",
      TRUE ~ NA_character_  # 都不匹配时设为 NA
    )
  )



df = Epi_meta
df_p = df[grepl("P",df$project),]
summary <- table(df_p[,c('subtype','group')])
roe <- as.data.frame(ROIE(summary))
long_df <- tidyr::gather(roe, key = "group", value = "Ro/e", colnames(roe)[1:ncol(roe)-1])
write.csv(roe, file = "/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Epi_subtype/Roe_P.csv",quote=F)


p = pheatmap::pheatmap(roe, 
                       color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(10),
                       cluster_rows = F, 
                       cluster_cols = F, 
                       display_numbers = TRUE,
                       number_color = "black",
                       cellwidth = 30, cellheight = 16, fontsize = 10,                    # fontsize字号大小
                       border_color = '#ffffff',
                       angle_col='45')
dev.off()
p$gtable$heights[2] = unit(15,"bigpts")

p$gtable[[2]][p$gtable[[2]]$name =="legend",1:4] = c(4,5,4,6)
p$gtable <- gtable_add_grob(p$gtable, textGrob("Ro/e",x = unit(0.32,"npc"),y = unit(0.7,"npc")), t = 2, l = 5, r= 5, b=2)
#pdf("heatmap.pdf",width = ncol(roe)*1.5,height = nrow(roe)*0.3)
pdf("/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Epi_subtype/heatmap_P.pdf",width = ncol(roe)*2,height = nrow(roe)*0.3)
p
dev.off()

###S,P,R
df = Epi_meta
summary <- table(df[,c('subtype','group_2')])
roe <- as.data.frame(ROIE(summary))
long_df <- tidyr::gather(roe, key = "group_2", value = "Ro/e", colnames(roe)[1:ncol(roe)-1])
write.csv(roe, file = "/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Epi_subtype/Roe_S_P_R.csv",quote=F)


p = pheatmap::pheatmap(roe, 
                       color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(10),
                       cluster_rows = F, 
                       cluster_cols = F, 
                       display_numbers = TRUE,
                       number_color = "black",
                       cellwidth = 30, cellheight = 16, fontsize = 10,                    # fontsize字号大小
                       border_color = '#ffffff',
                       angle_col='45')
p$gtable$heights[2] = unit(15,"bigpts")

p$gtable[[2]][p$gtable[[2]]$name =="legend",1:4] = c(4,5,4,6)
p$gtable <- gtable_add_grob(p$gtable, textGrob("Ro/e",x = unit(0.32,"npc"),y = unit(0.7,"npc")), t = 2, l = 5, r= 5, b=2)
#pdf("heatmap.pdf",width = ncol(roe)*1.5,height = nrow(roe)*0.3)
pdf("/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Epi_subtype/heatmap_S_P_R.pdf",width = ncol(roe)*2,height = nrow(roe)*0.3)
p
dev.off()


###group_3
df = Epi_meta
df <- df[df$group_3%in%c("A_P","A_R","A_S"),]
summary <- table(df[,c('subtype','group_3')])
roe <- as.data.frame(ROIE(summary))
long_df <- tidyr::gather(roe, key = "group_3", value = "Ro/e", colnames(roe)[1:ncol(roe)-1])
write.csv(roe, file = "/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Epi_subtype/Roe_group_3.csv",quote=F)


p = pheatmap::pheatmap(roe, 
                       color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(10),
                       cluster_rows = F, 
                       cluster_cols = F, 
                       display_numbers = TRUE,
                       number_color = "black",
                       cellwidth = 30, cellheight = 16, fontsize = 10,                    # fontsize字号大小
                       border_color = '#ffffff',
                       angle_col='45')
p$gtable$heights[2] = unit(15,"bigpts")

p$gtable[[2]][p$gtable[[2]]$name =="legend",1:4] = c(4,5,4,6)
p$gtable <- gtable_add_grob(p$gtable, textGrob("Ro/e",x = unit(0.32,"npc"),y = unit(0.7,"npc")), t = 2, l = 5, r= 5, b=2)
#pdf("heatmap.pdf",width = ncol(roe)*1.5,height = nrow(roe)*0.3)
pdf("/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Epi_subtype/heatmap_group_3.pdf",width = ncol(roe)*2,height = nrow(roe)*0.3)
p
dev.off()


#########Basal sub
Basal_data <- subset(Epi_clean,subtype%in%c("Basal cell","Proliferating basal cell"))
Basal_out_dir <- "Basal_subtype" 
dir.create(Basal_out_dir)

## PCA --------
Basal_data <- FindVariableFeatures(Basal_data)

options(future.globals.maxSize = 20 * 1024 ^ 3) # 20Gb
Basal_data <- ScaleData(Basal_data, vars.to.regress = 'nCount_RNA')
Basal_data <- RunPCA(Basal_data)

### pca inspector ------
pdf(file.path(Basal_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Basal_data, ndims = 50)
dev.off()


Basal_data$annotation %>% table()
#Basal_data$celltype %>% table()
## dimension reduction --------
dims <- 15 # 18

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 20 * 1024 ^ 3) # 20Gb

Basal_data <- RunUMAP(Basal_data, dims = 1:dims,return.model = TRUE)

pdf(file.path(Basal_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Basal_data, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

Basal_data <- Basal_data %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Basal_data <- RunUMAP(
  Basal_data, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(Basal_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Basal_data, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(Basal_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Basal_data, group.by = "sample", label = TRUE, reduction = "umap") + NoLegend()
dev.off()

reduction.use <- 'pca' # 'harmony' 'pca'
### cluster -------剔除部分复发细胞
Basal_data  <-
  FindNeighbors(
    Basal_data ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Basal_data <-
  FindClusters(
    Basal_data,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(Basal_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Basal_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(Basal_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(Basal_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(Basal_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'umap') + NoLegend()
# DimPlot(Basal_data, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Basal_data, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Basal_data, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

saveRDS(Basal_data, file.path(Basal_out_dir, 'Basal_data_all.rds'))
Basal_data <- readRDS(file.path(Basal_out_dir, 'Basal_data_all.rds'))
#####剔除两个复发的大部分细胞

Basal_data_use <- subset(Basal_data,`RNA_snn_res.0.6`%in%c("12","11","8"))
table(Basal_data_use$sample)

## PCA --------
Basal_data_use <- FindVariableFeatures(Basal_data_use)

options(future.globals.maxSize = 20 * 1024 ^ 3) # 20Gb
Basal_data_use <- ScaleData(Basal_data_use, vars.to.regress = 'nCount_RNA')
Basal_data_use <- RunPCA(Basal_data_use)

### pca inspector ------
pdf(file.path(Basal_out_dir, 'pca plot use.pdf'), height = 6, width = 8)
ElbowPlot(Basal_data_use, ndims = 50)
dev.off()


Basal_data_use$subtype %>% table()
#Basal_data$celltype %>% table()
## dimension reduction --------
dims <- 10 # 18

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 20 * 1024 ^ 3) # 20Gb

Basal_data_use <- RunUMAP(Basal_data_use, dims = 1:dims,return.model = TRUE)

Basal_data_use <- Basal_data_use %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Basal_data_use <- RunUMAP(
  Basal_data_use, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(Basal_out_dir, str_glue("Dimplot dim{dims} sample  - harmony use.pdf")), height = 6, width = 8)
DimPlot(Basal_data_use, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

### cluster -------
Basal_data_use  <-
  FindNeighbors(
    Basal_data_use ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Basal_data_use <-
  FindClusters(
    Basal_data_use,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(Basal_out_dir, str_glue("Dimplot dim{dims} cluster use.pdf")), height = 6, width = 8)
DimPlot(Basal_data_use, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Basal_data_use, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Basal_data_use, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Basal_data_use, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Basal_data, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Basal_data, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Basal_data, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()
###########
Idents(Basal_data_use) <- "RNA_snn_res.0.4"

###marker_ls
# 2. 过滤 marker 列表，只保留存在的基因
Basal_markers <-
  FindAllMarkers(
    JoinLayers(Basal_data_use),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Basal_markers, file = "Basal_subtype/Basal_all_markers.csv", row.names = TRUE)
Basal_markers <- read.csv("Basal_subtype/Basal_all_markers.csv")

Basal_top5_markers <- 
  Basal_markers %>% 
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

names(Basal_top5_markers) <- sub('.*\\.','',names(Basal_top5_markers))


# dotplot
# 如果已做其他细胞类型的亚群分析，这里就可以顺便画一下其他类型的细胞亚群基因，比如Fibro_major_marker_ls2，否则只用画T细胞
pdf(file.path(Basal_out_dir, str_glue("dotplot dims{dims} major_marker use.pdf")), width = 25, height = 6)
DotPlot(Basal_data_use, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Basal_data_use, group.by = "RNA_snn_res.0.4", features = Epi_major_marker_ls)+RotatedAxis()+theme_90angle
#DotPlot(Macrophage, group.by = "RNA_snn_res.0.4", features = Macrophage_major_marker_ls_1)+RotatedAxis()+theme_90angle
#DotPlot(Macrophage, group.by = "RNA_snn_res.0.6", features = Macrophage_markers)+RotatedAxis()
DotPlot(Basal_data_use, group.by = "RNA_snn_res.0.4", features = Basal_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(Macrophage, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()


pdf(file.path(Basal_out_dir, "ncount_vlnplot.pdf"), width = 8, height = 4)
VlnPlot(Basal_data_use, features = c("nCount_RNA", "nFeature_RNA"), group.by = "RNA_snn_res.0.2", pt.size = 0)

VlnPlot(Basal_data_use, features = c("TP63"), group.by = "subtype", pt.size = 0)
dev.off()
#####
cell_cluster <- c(
  '0'  = "Basal", # 
  '1'  = "Unknown", #
  '2'  = "AKT3+ basal", #
  '3'  = "AT1/AT2-like", # 
  '4'  = "DSG3+ basal", # 
  '5'  = "Pulmonary neuroendocrine cell", # 
  '6'  = "Unknown"
)

Idents(Basal_data_use) <- "RNA_snn_res.0.4"



Basal_data_use <- RenameIdents(Basal_data_use, cell_cluster)
Basal_data_use$subtype <- Idents(Basal_data_use)


pdf(file.path(Basal_out_dir, str_glue("Dimplot dim{dims} annotation use 2.pdf")), height = 6, width = 8)
DimPlot(Basal_data_use, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Basal_out_dir, str_glue("Dimplot dim{dims} annotation split use 2.pdf")), height = 6, width = 15)
DimPlot(Basal_data_use, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

saveRDS(Basal_data_use, file.path(Basal_out_dir, 'Basal_data_use.rds'))


Basal_data_use <- readRDS(file.path(Basal_out_dir, 'Basal_data_use.rds'))

####
Basal_clean <- Basal_data_use[, !Basal_data_use$subtype %in% c('Unknown',"AT1/AT2-like","Pulmonary neuroendocrine cell")]
Basal_clean <- RunUMAP(Basal_clean, dims = 1:dims, return.model = TRUE)


pdf(file.path(Basal_out_dir,str_glue("Dimplot dim{dims} annotation - pure use 2.pdf")), height = 6, width = 8)
DimPlot(Basal_clean, group.by = "subtype", label = TRUE, repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Basal_out_dir, str_glue("Dimplot dim{dims} annotation split by group - pure use 2.pdf")), height = 7, width = 18)
DimPlot(Basal_clean, group.by = "subtype",split.by = 'group', label = TRUE, repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Basal_out_dir, str_glue("Dimplot dim{dims} annotation split by group_2 - pure use 2.pdf")), height = 7, width = 24)
DimPlot(Basal_clean, group.by = "subtype",split.by = 'group_2', label = TRUE, repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Basal_out_dir, str_glue("Dimplot dim{dims} annotation split by group_3 - pure use 2.pdf")), height = 7, width = 40)
DimPlot(Basal_clean, group.by = "subtype",split.by = 'group_3', label = TRUE, repel = TRUE) + NoLegend()
dev.off()

saveRDS(Basal_clean, file.path(Basal_out_dir, 'Basal_data_clean.rds'))

pdf(file.path(Basal_out_dir, str_glue("Dimplot dim{dims} annotation split by sample - pure 2.pdf")), 
    height = 14, width = 30)

# 2. 绘制基础图形（注意：这里不要加 NoLegend()，因为分面后图例可能丢失）
p <- DimPlot(Basal_clean, group.by = "subtype", split.by = 'sample', 
             label = TRUE, repel = TRUE)

# 3. 使用 facet_wrap 将图形强制排列为两行
p + facet_wrap(~sample, nrow = 2) + NoLegend()

# 4. 关闭PDF设备
dev.off()

####compare
df = Basal_clean[[]]
df$subtype <- as.character(df$subtype)
df$subtype <- as.factor(df$subtype)

#df <- df[df$group_3%in%c("A_P","A_R","A_S"),]
summary <- table(df[,c('subtype','group_3')])
roe <- as.data.frame(ROIE(summary))
long_df <- tidyr::gather(roe, key = "group_3", value = "Ro/e", colnames(roe)[1:ncol(roe)-1])
write.csv(roe, file = "/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Basal_subtype/Roe_group_3.csv",quote=F)


p = pheatmap::pheatmap(roe, 
                       color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(10),
                       cluster_rows = F, 
                       cluster_cols = F, 
                       display_numbers = TRUE,
                       number_color = "black",
                       cellwidth = 30, cellheight = 16, fontsize = 10,                    # fontsize字号大小
                       border_color = '#ffffff',
                       angle_col='45')
p$gtable$heights[2] = unit(15,"bigpts")

p$gtable[[2]][p$gtable[[2]]$name =="legend",1:4] = c(4,5,4,6)
p$gtable <- gtable_add_grob(p$gtable, textGrob("Ro/e",x = unit(0.32,"npc"),y = unit(0.7,"npc")), t = 2, l = 5, r= 5, b=2)
#pdf("heatmap.pdf",width = ncol(roe)*1.5,height = nrow(roe)*0.3)
pdf("/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Basal_subtype/heatmap_group_3.pdf",width = ncol(roe)*1.5,height = nrow(roe)*0.6)
p
dev.off()


summary_t <- as.data.frame(summary)
summary_prop <- summary_t %>%
  group_by(sample) %>%             # 按分组计算
  mutate(total = sum(Freq), prop = Freq/total) %>%  # 计算组内总数、占比
  ungroup()


basal_colors <- c(
  "Basal"         = "#4A7C59",   # 深绿 - 基础/正常
  "AKT3+ basal"       = "#F77F00",   # 橙红 - 侵袭/危险
  "DSG3+ basal"      = "#9B5DE5"   # 紫色 - 淋巴/特殊微环境
)
p2 <- ggplot(summary_prop, aes(x = group_3, y = prop, fill = subtype)) +
  geom_col(position = "stack", width = 0.7) +  # 堆叠柱状图
  scale_y_continuous(labels = scales::percent_format(), expand = c(0,0)) + # y轴转为百分比
  # 配色
  scale_fill_manual(values = basal_colors) +
  # 坐标轴、标签主题设置
  labs(
    x = "Group",
    y = "Cell proportion",
    fill = "Cell type"
  ) +
  # 论文级主题美化
  theme_classic() +
  theme(
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 13, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    legend.position = "right",
    plot.margin = margin(10,20,10,10)
  )

pdf("/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/Basal_subtype/subtypetype_prop.pdf",width = 8,height = 8)
p2
dev.off()

##
p_col_proj <-
  df %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  scale_fill_manual(values = basal_colors) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path("Basal_subtype/cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  df %>%
  group_by(group_3, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group_3, y = freq, fill = subtype)) + 
  geom_col() +
  scale_fill_manual(values = basal_colors) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


ggsave(
  p_col_proj,
  filename = file.path("Basal_subtype/cell_freq_group3.pdf"),
  height = 6,
  width  = 4
)

####
Basal_clean <- AddModuleScore(
  object = Basal_clean,
  features = list(
    Glycolysis = pathway_list[["HALLMARK_GLYCOLYSIS"]],
    Hypoxia = pathway_list[["HALLMARK_HYPOXIA"]],
    OXPHOS = pathway_list[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]]
  ),
  name = c("Glycolysis_", "Hypoxia_", "OXPHOS_")
)
pdf(file.path(Basal_out_dir,str_glue("Feature score Basal.pdf")), height = 6, width = 30)
FeaturePlot(Basal_clean, features = "Glycolysis_1",split.by = "group_3", label = FALSE, repel = TRUE)
FeaturePlot(Basal_clean, features = "Hypoxia_2",split.by = "group_3", label = FALSE, repel = TRUE)
FeaturePlot(Basal_clean, features = "OXPHOS_3",split.by = "group_3", label = FALSE, repel = TRUE)
dev.off()


pdf(file.path(Basal_out_dir,str_glue("Feature tumor marker Basal.pdf")), height = 6, width = 30)
FeaturePlot(Basal_clean, features = "MUC1",split.by = "group_3", label = FALSE, repel = TRUE)
FeaturePlot(Basal_clean, features = "TP63",split.by = "group_3", label = FALSE, repel = TRUE)
FeaturePlot(Basal_clean, features = "SOX2",split.by = "group_3", label = FALSE, repel = TRUE)
FeaturePlot(Basal_clean, features = "DSG3",split.by = "group_3", label = FALSE, repel = TRUE)
FeaturePlot(Basal_clean, features = "AKT3",split.by = "group_3", label = FALSE, repel = TRUE)
dev.off()

pdf(file.path(Basal_out_dir,str_glue("Feature TP63 Basal.pdf")), height = 6, width = 6)
FeaturePlot(Basal_clean, features = "TP63", label = FALSE, repel = TRUE)
dev.off()
#####
Basal_clean <- readRDS("Basal_subtype/Basal_data_clean.rds")
head(Idents(Basal_clean))
Basal_markers <-
  FindAllMarkers(
    JoinLayers(Basal_clean),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Basal_markers, file = "Basal_subtype/Basal_all_markers_celltype.csv", row.names = TRUE)

Basal_top5_markers <- 
  Basal_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 30) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%

  nest() %>%
  deframe() %>% 
  map(pull, 'gene') %>%
  imap(function(genes, name){
    re <- genes %>%
      map(function(gene){
        # browser()
        names(gene) = str_glue("{name}")
        gene
      })
    re
  }) %>%
  unlist() %>% 
  .[!duplicated(.)]

names(Basal_top5_markers) <- sub('.*\\.','',names(Basal_top5_markers))

pdf(file.path(Basal_out_dir, str_glue("dotplottop5_marker subtype.pdf")), width = 25, height = 6)
DotPlot(Basal_clean, group.by = "subtype", features = major_marker_ls)+RotatedAxis()+scale_size(range = c(3, 6))+theme_90angle
DotPlot(Basal_clean, group.by = "subtype", features = Epi_major_marker_ls)+RotatedAxis()+scale_size(range = c(3, 6))+theme_90angle
DotPlot(Basal_clean, group.by = "subtype", features = Basal_top5_markers)+RotatedAxis()+scale_size(range = c(3, 6))+theme_90angle
#DotPlot(Macrophage, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()



####
out_dir <- "spatial plot"
dir.create(out_dir)

all_data <- readRDS("/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/data/all_data.rds")
all_data$subtype <- as.character(all_data$annotation)  # 或 NA
# 提取 KC_data 的细胞名和对应注释
Basal_cells <- colnames(Basal_clean)
Basal_subtype <- Basal_clean$subtype  # 如 KC_data$cell_type 或 KC_data$RNA_snn_res.0.4

unique(Basal_clean$subtype)
# 创建命名向量
subtype_map <- setNames(as.character(Basal_subtype), Basal_cells)

# 赋值到 all_data
all_data$subtype[match(Basal_cells, colnames(all_data))] <- subtype_map[Basal_cells]

# 非 KC 细胞保留原注释
non_kc <- setdiff(colnames(all_data), Basal_cells)
all_data$subtype[non_kc] <- as.character(all_data$annotation[non_kc])
unique(all_data$subtype)



annotation_col <- 'subtype'

cells_position_ls <-
  cbind(
    all_data[['sample']],
    Embeddings(all_data[['spatial']])
  ) %>%
  rownames_to_column('cell_id') %>%
  group_by(sample) %>%
  nest() %>%
  deframe()


plot_dat <-
  FetchData(
    all_data, 
    vars = c(annotation_col)
  ) %>%
  rownames_to_column('cell_id')


basal_colors <- c(
  "Basal"         = "#4A7C59",   # 深绿 - 基础/正常
  "AKT3+ basal"       = "#F77F00",   # 橙红 - 侵袭/危险
  "DSG3+ basal"      = "#9B5DE5"   # 紫色 - 淋巴/特殊微环境
)
all_cell_types <- names(basal_colors)

# 遍历每种细胞类型
for(target_cell_type in all_cell_types){
  
  message("正在绘制：", target_cell_type)
  
  # 强制锁定当前细胞类型和颜色
  current_cell <- target_cell_type
  current_color <- basal_colors[[current_cell]]
  
  p_sp_ct_ls <- cells_position_ls %>%
    imap(function(coords, sam){
      
      ct <- current_cell
      
      dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
      
      # ========== 关键修复1：检查该样本是否包含目标细胞类型 ==========
      cell_types_vec <- dat2 %>% pull(!!sym(annotation_col))
      
      if (!ct %in% cell_types_vec) {
        message("  跳过样本 ", sam, "：无 ", ct, " 细胞")
        return(NULL)
      }
      
      # ========== 关键修复2：统一字符型比较，避免 factor 问题 ==========
      cell_types_char <- as.character(cell_types_vec)
      ct_char <- as.character(ct)
      
      highlight_vec <- ifelse(cell_types_char == ct_char, ct_char, "Other")
      highlight_vec[is.na(highlight_vec)] <- "Other"
      
      dat2$highlight <- highlight_vec
      
      dat2$highlight <- factor(
        dat2$highlight,
        levels = c(ct_char, "Other")
      )
      
      # ========== 绘图 ==========
      fill_values <- setNames(c(current_color, "gray60"), c(ct_char, "Other"))
      alpha_values <- setNames(c(1, 0.3), c(ct_char, "Other"))
      
      p <- ggplot(dat2, aes(x = SPATIAL_1, y = SPATIAL_2)) +
        geom_tile(aes(fill = highlight, alpha = highlight), height = 10, width = 10) +
        theme_void() +
        coord_equal() +
        scale_y_reverse() +
        scale_alpha_manual(values = alpha_values) +
        scale_fill_manual(
          values = fill_values,
          guide = guide_legend(override.aes = list(alpha = c(1, 1)))
        ) +
        theme(
          panel.background = element_rect(fill = NA, color = NA),
          plot.title = element_text(hjust = 0.5)
        ) +
        ggtitle(paste0(sam, " (", ct, ")"))
      
      return(p)
    }) %>%
    # ========== 关键修复3：过滤掉 NULL（跳过的样本）==========
  compact()
  
  # 如果所有样本都被跳过，给出提示并继续下一个细胞类型
  if (length(p_sp_ct_ls) == 0) {
    message("  警告：", target_cell_type, " 在所有样本中均未找到，跳过输出")
    next
  }
  
  # 输出PDF（根据样本数调整宽度）
  pdf(
    file.path(out_dir, paste0("spatial plot of ", target_cell_type, ".pdf")),
    height = 8, width = 8 
  )
  walk(p_sp_ct_ls, print)
  dev.off()
  
  message("  完成：", target_cell_type, "，共 ", length(p_sp_ct_ls), " 个样本")
}
#####功能

library(clusterProfiler)
library(org.Hs.eg.db)

# 按 cluster 分组做富集
marker_list <- Basal_markers %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  split(.$cluster)

enrichGO_results <- imap(marker_list, function(df, cluster_name) {
  message("正在富集：", cluster_name, " (", nrow(df), " genes)")
  
  entrez <- bitr(df$gene, "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID
  
  if (length(entrez) < 10) return(NULL)
  
  enrichGO(
    gene = entrez,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    readable = TRUE
  )
}) %>% compact()

names(enrichGO_results)
enrichGO_results
# 可视化
iwalk(enrichGO_results, function(ego, cluster_name) {
  
  # 检查是否有显著结果
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    message(cluster_name, " 无显著富集结果")
    return(NULL)
  }
  
  # 点图
  p <- dotplot(ego, showCategory = 15) +
    ggtitle(paste0(cluster_name, " - GO BP Enrichment")) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 10)
    )
  
  # 保存
  ggsave(
    filename = file.path(out_dir, paste0("GO_dotplot_", cluster_name, ".pdf")),
    plot = p,
    width = 10,
    height = 8
  )
  
  message("已保存: ", cluster_name)
})
