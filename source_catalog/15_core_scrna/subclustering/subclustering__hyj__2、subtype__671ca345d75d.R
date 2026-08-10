result_dir <- 'result' 
dir.create(result_dir)

#read RDS----
library(Seurat)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Hs.eg.db)

all_data <- readRDS(file.path("R_OUTPUT", "all_data.rds"))
all_data2 <- all_data
#xenium <- readRDS("xenium_FOV.rds")
#xenium <- NormalizeData(xenium)

xenium$annotation %>% table()

result_dir <- "result"
################Add_umap#####################
# new_umap <- umap[row_cell, ] %>% as.matrix()
# colnames(new_umap) <- c("UMAP_1", "UMAP_2")
harmonyumap <- all_data@reductions$harmonyumap@cell.embeddings
rownames(harmonyumap) <- rownames(xenium@meta.data)

harmonyumap_reduc <- CreateDimReducObject(
  embeddings = harmonyumap, 
  key = "UMAP_",
  assay = "Xenium"
)
xenium@reductions$harmonyumap <- harmonyumap_reduc

pdf(file.path(result_dir, str_glue("Dimplot sample - harmony.pdf")), height = 6, width = 8)
DimPlot(xenium, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(result_dir, str_glue("Dimplot group - harmony.pdf")), height = 6, width = 8)
DimPlot(xenium, group.by = "group", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(result_dir, str_glue("Dimplot annotation - harmony.pdf")), height = 6, width = 8)
DimPlot(xenium, group.by = "annotation", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


major_marker <- c(
  # 1. B细胞
  "B cell" = "CD19",
  "B cell" = "CD79A", 
  "B cell" = "CD79B",
  "B cell" = "MS4A1",
  
  # 2. 浆细胞
  "Plasma" = "TNFRSF17",
  "Plasma" = "MZB1",
  
  # 3. T细胞
  "T cell" = "CD3E",
  "T cell" = "CD3G",
  "T cell" = "CD4",
  "T cell" = "FOXP3",
  "T cell" = "CD8A", 
  "T cell" = "GZMB",
  
  # 4. NK细胞
  "NK" = "NCAM1",
  
  # 5. 增殖细胞
  "Proliferating" = "MKI67",
  "Proliferating" = "TOP2A",
  
  # 6. 髓系细胞
  "Myeloid" = "ITGAM",
  
  # 7. 巨噬细胞
  "Macrophage" = "CD68",
  "Macrophage" = "CD163",
  "Macrophage" = "CSF1R",
  
  # 8. 肥大细胞
  "Mast cell" = "CTSG",
  "Mast cell" = "MS4A2",
  "Mast cell" = "KIT",
  "Mast cell" = "HDC",
  
  # 9. 单核细胞
  "Monocytes" = "CD14",
  "Monocytes" = "FCGR3A",
  
  # 10. 树突状细胞
  "Dendritic" = "FCER1A",
  "Dendritic" = "CD1C",
  "Dendritic" = "CD83",
  "Dendritic" = "LILRA4",
  "Dendritic" = "CLEC9A",
  "Dendritic" = "IDO1",
  
  # 11. 中性粒细胞
  "Neutrophil" = "CSF3R",  # 注意: 修正了拼写 Neutrophi -> Neutrophil
  "Neutrophil" = "FCGR3B",
  
  # 12. 内皮细胞
  "Endothelial" = "PECAM1",
  "Endothelial" = "CLDN5",
  
  # 13. 上皮细胞
  "Epithelial" = "EPCAM",
  
  # 14. 成纤维细胞
  "Fibroblasts" = "PDGFRA",
  "Fibroblasts" = "CCDC80",
  
  # 15. 平滑肌细胞
  "Smooth muscle cells" = "PDGFRB",
  "Smooth muscle cells" = "PLN",
  "Smooth muscle cells" = "MYLK"
)


valid_genes <- major_marker[major_marker %in% rownames(all_data2)]
cat("有效基因:", length(valid_genes), "/", length(major_marker))

missing_gene <- setdiff(major_marker, valid_genes)
missing_gene

pdf(file.path(result_dir, str_glue("dotplot_main_marker.pdf")), width = 25, height = 10)
DotPlot(xenium, group.by = "res_0.8", features = valid_genes, dot.scale = 10) + theme_90angle
DotPlot(xenium, group.by = "annotation", features = valid_genes, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "sketch_snn_res.0.8", features = major_marker_ls1, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "sketch_snn_res.0.8", features = major_marker_ls2, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "sketch_snn_res.0.8", features = major_marker_ls3, dot.scale = 10) + theme_90angle
dev.off()


#映射add亚类注释到all_data----
library(tibble)
library(data.table)
subtype <- rbind(ENDO[['subtype']] %>% rownames_to_column('cells'), 
                 EPI[['subtype']] %>% rownames_to_column('cells'),
                 LYMPHOID[['subtype']]%>% rownames_to_column('cells'),
                 MYELOID[['subtype']]%>% rownames_to_column('cells'),
                 MESENCHYMAL[['subtype']]%>% rownames_to_column('cells')
)

meta.data <-
  subtype%>%
  right_join(
    xenium[['annotation']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    subtype = subtype
    # subtype = ifelse((as.character(annotation) %in% c('T cells','Fibroblasts',"Myeloids")),
    #                      as.character(subtype), 
    #                      as.character(annotation))
  )%>%
  dplyr::select(cells, subtype) %>%
  column_to_rownames('cells')


xenium<-AddMetaData(xenium,meta.data['subtype'])
xenium$annotation %>% table()
xenium$subtype %>% table()

# xenium$annotation2 <- NULL

# 确认是否加入
table(xenium$subtype, xenium$annotation, useNA = "ifany")

#1、亚类注释------
# xenium@meta.data$subtype %>% table()
# 
# group_mapping <- c(
#   "Arteriole" = "Endothelial",
#   "Monocyte" = "Immune",
#   "Venous" = "Endothelial",
#   "AT2" = "Epithelial",
#   "AT1" = "Epithelial",
#   "Activated fibroblasts FBs" = "Mesenchymal",
#   "Capillary" = "Endothelial",
#   "CD4 T" = "Immune",
#   "NK" = "Immune",
#   "Neutrophils" = "Immune",
#   "Spp1 macrophage" = "Immune",
#   "Alveolar FBs" = "Mesenchymal",
#   "SMCs or fibrocytes FBs" = "Mesenchymal",
#   "Subpleural mesothelial" = "Mesenchymal",
#   "cDCs" = "Immune",
#   "B cell" = "Immune",
#   "Alveolar macrophage" = "Immune",
#   "Proliferating myeloid" = "Immune",
#   "Myofibroblasts FBs" = "Mesenchymal",
#   "Lymphatic" = "Endothelial",
#   "Proliferating T cell" = "Immune",
#   "Tregs" = "Immune",
#   "Interstitial macrophage" = "Immune",
#   "Ciliated" = "Epithelial",
#   "Plasma" = "Immune",
#   "Migratory DCs" = "Immune",
#   "Inf_Epi" = "Epithelial",
#   "delete" = "delete"
# )
# 
# xenium@meta.data$cell_group <- group_mapping[match(
#   xenium@meta.data$subtype,  # 待匹配的subtype列
#   names(group_mapping)       # 映射关系中的亚型名称
# )]
# 
# table(xenium@meta.data$subtype, xenium@meta.data$cell_group)
# 
# saveRDS(xenium, 'result/xenium_subtype_cellgroup.rds')
# 
# 
table(all_data2@meta.data$celltype, all_data2@meta.data$sample)


#ENDO--------
ENDO <- subset(all_data2, annotation == 'Endothelial')
ENDO@images <- list()

## PCA --------
DefaultAssay(ENDO) <- "RNA"
ENDO <- FindVariableFeatures(ENDO)
ENDO <- ScaleData(ENDO, vars.to.regress = 'nCount_RNA')
ENDO <- RunPCA(ENDO)

ENDO_out_dir<-'result/ENDO_subtype/' 
dir.create(ENDO_out_dir)


### pca inspector ------
pdf(file.path(ENDO_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(ENDO, ndims = 50)
dev.off()

# pdf(file.path(ENDO_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
# DimHeatmap(ENDO, dims = 1:5)
# dev.off()
# pdf(file.path(ENDO_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
# DimHeatmap(ENDO, dims = 6:10)
# dev.off()
# pdf(file.path(ENDO_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
# DimHeatmap(ENDO, dims = 11:15)
# dev.off()
# pdf(file.path(ENDO_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
# DimHeatmap(ENDO, dims = 16:20)
# dev.off()


## dimension reduction --------
dims <- 15

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


ENDO <- RunUMAP(ENDO, dims = 1:dims,return.model = TRUE)

pdf(file.path(ENDO_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(ENDO, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

ENDO <- ENDO %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

ENDO <- RunUMAP(
  ENDO, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(ENDO_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(ENDO, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'

### cluster -------
ENDO  <-
  FindNeighbors(
    ENDO ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


ENDO <-
  FindClusters(
    ENDO,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(ENDO_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(ENDO, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(ENDO, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(ENDO, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(ENDO, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


### dotplot -------
endo_major_marker_ls <- c(
  # 补体产生细胞 (巨噬细胞/肝细胞)
  "Complement_producing" = "C1QA",
  "Complement_producing" = "C1QB", 
  "Complement_producing" = "C1QC",
  "Complement_producing" = "C2",
  "Complement_producing" = "C3",
  "Complement_producing" = "C3AR1",
  
  # DAMP信号细胞 (应激/损伤细胞)
  "DAMP_signaling" = "HMGB1",
  "DAMP_signaling" = "HSPA1A",
  "DAMP_signaling" = "S100A8",
  "DAMP_signaling" = "S100A9", 
  "DAMP_signaling" = "S100A12",
  
  # 模式识别受体细胞 (先天免疫细胞)
  "PRR_expressing" = "TLR1",
  "PRR_expressing" = "TLR5",
  "PRR_expressing" = "TLR7",
  "PRR_expressing" = "IFIH1",
  "PRR_expressing" = "DHX58",
  

  "IGHG1",
  "IGHG2",
  "IGHG3", 
  "IGHG4",
  "IGHM",
  "IGKC",
  "IGLC2",
  "IGHA1",
  "IGHA2",
  
  "FCGR1A",
  "FCGR2A", 
  "FCGR3A"
)

valid_endogenes <- endo_major_marker_ls[endo_major_marker_ls %in% rownames(ENDO)]
cat("有效基因:", length(valid_endogenes), "/", length(endo_major_marker_ls))

# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_gene <- setdiff(endo_major_marker_ls, valid_endogenes)
missing_gene


# find markers
Idents(ENDO) <- "RNA_snn_res.0.4"
endo_markers <-
  FindAllMarkers(
    JoinLayers(ENDO),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

endo_top10_markers <- 
  endo_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10) %>%
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

names(endo_top10_markers) <- sub('.*\\.','',names(endo_top10_markers))



#提取deepseek
endo_top10_markers1 <-
  endo_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10)



# dotplot
pdf(file.path(ENDO_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 35, height = 6)
DotPlot(ENDO, group.by = "RNA_snn_res.0.4", features = major_marker)+RotatedAxis()
DotPlot(ENDO, group.by = "RNA_snn_res.0.4", features = endo_major_marker_ls)+RotatedAxis()
DotPlot(ENDO, group.by = "RNA_snn_res.0.4", features = endo_top10_markers)+RotatedAxis()
dev.off()


for (gene in valid_endogenes) {
  p <- FeaturePlot(ENDO, features = gene, reduction = 'harmonyumap')+RotatedAxis()
  ggsave(
    filename = file.path(ENDO_out_dir, paste0("featureplot/",gene, "_feature.png")),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}


##空间imagefeatureplot-------
SegImagefeature_dir <- 'result/ENDO_subtype/SegImagefeature/'
dir.create(SegImagefeature_dir)

source("ImageFeature-plot.R") 

all_colors <- c(
  "#A6CEE3", "#5B9EC9", "#2D82AF", "#7EBA98", "#98D277",
  "#52AF43", "#6F9E4C", "#DD9A88", "#FDBB69", "#FE982C",
  "#F78620", "#D9A295", "#feee7d", "#02c9c9", "#600473",
  "#B294C7", "#7D54A5", "#9E8099", "#F0EB99", "#DBB466",
  "#B15928", # 原始21色
  "#FF6B6B", "#4ECDC4", "#FF9F1C", "#7209B7", "#4361EE", 
  "#F15BB5", "#00BBF9", "#00F5D4", "#9B5DE5", "#00F5A7",
  "#FEE440", "#9C89B8", "#06D6A0"  # 新增13色
)

all_colors10 <- c(
  "#FF6B6B","#F15BB5","#A6CEE3","#00F5A7", "#FDBB69","#feee7d", "#9E8099",
  "#00BBF9","#7209B7","#9B5DE5"
)

# all_colors10 <- c(
#   "gray90","#F15BB5","gray90","gray90","gray90","gray90","gray90","gray90","gray90","gray90"
# )

Idents(xenium) <- "annotation"
names(all_colors10) <- levels(xenium)


###所有细胞类型------
for (sample_id in unique(xenium@meta.data$sample)){
  for (gene in valid_endogenes) {
  
    P2 <- ImageFeaturePlot_wyt(xenium, 
                               features = gene, 
                               size = 0.75, 
                               cols = c('#003366', '#FFD700'),
                               fov =sample_id,
                               nmols = 20000,
                               border.color = all_colors10, 
                               border.size = 0.3,
                               border.by = "celltype",
                               dark.background = T
                               )
    
    ggsave(P2,filename = paste0(SegImagefeature_dir,sample_id,"_",gene,"_Imagefeature.pdf"),width = 25,height = 25,dpi = 300)
    ggsave(P2,filename = paste0(SegImagefeature_dir,sample_id,"_",gene,"_Imagefeature.png"),width = 25,height = 25,dpi = 300)
    
    message("成功保存: ", sample_id, "_", gene)
    }
}

###只画内皮的--------
for (sample_id in unique(xenium@meta.data$sample)){
  for (gene in valid_endogenes) {
    
    P4 <- ImageFeaturePlot_wyt(xenium, 
                               features = gene, 
                               size = 0.5, 
                               cols = c('#003366', '#FFD700'),
                               fov =sample_id,
                               nmols = 20000,
                               border.color = all_colors10, 
                               border.size = 0.2,
                               border.by = "celltype",
                               cells.use = c("Endothelial")
                               )
    ggsave(P4,filename = paste0(SegImagefeature_dir,sample_id,"_",gene,"_Imagefeature_ENDO.pdf"),width = 25,height = 25,dpi = 300)
    ggsave(P4,filename = paste0(SegImagefeature_dir,sample_id,"_",gene,"_Imagefeature_ENDO.png"),width = 25,height = 25,dpi = 300)
    
  }
}






# annotation -----------
cell_cluster <- c(
  '0' = "0", # HEY1 VEGFA HERPUD1
  '1' = "2", # IL7R FCN3 CA4 CLDN5 APLN APLNR
  '2' = "2", # PLVAP ACKR1 POSTN COL15A1
  '3' = "3", # PLVAP ACKR1 POSTN COL15A1
  '4' = "4", # CCL21 FABP4 MRC1 GNG11
  '5' = "5", # PLVAP ACKR1 POSTN COL15A1
  '6' = "6", # HEY1 VEGFA HERPUD1
  '7' = "7", # no marker
  '8' = "8", # no marker
  '9' = '9'
)

Idents(ENDO) <- "RNA_snn_res.0.4"

ENDO <- RenameIdents(ENDO, cell_cluster)
ENDO$subtype <- Idents(ENDO)

pdf(file.path(ENDO_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(ENDO, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(ENDO_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 7, width = 18)
DimPlot(ENDO, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


## cell props barplot  --------
ENDO_meta <- ENDO[[]]

ENDO_meta$sample <- factor(ENDO_meta$sample, 
                           levels=c(ENDO_meta$sample %>% unique())
)
p_col_proj <-
  ENDO_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path("result/ENDO_subtype/cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  ENDO_meta %>%
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
  filename = file.path("result/ENDO_subtype/cell_freq_group.pdf"),
  height = 6,
  width  = 4
)

#ENDO <- readRDS("result/ENDO_subtype/Endothelial.rds")
## volcano-------
endo_markers <-
  FindAllMarkers(
    JoinLayers(ENDO),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(endo_markers,'result/ENDO_subtype/endo_diff_gene.csv')

library(EnhancedVolcano)
type<-unique(endo_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- endo_markers[endo_markers$cluster == type, ]
    p<-EnhancedVolcano(markers,
                       lab = rownames(markers),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = str_glue('{name} volcano plot'),
                       pCutoff = 10e-2,
                       FCcutoff = 0,
                       pointSize = 3.0,
                       labSize = 6.0)
    p
  })
volcano_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path('result/ENDO_subtype', str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })

endo_markers_ls <- list(endo_markers)
names(endo_markers_ls) <- 'Endothelial'

# scatter
p_volcanoplot_major_ls <- 
  endo_markers_ls %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 0.5
    dodge_width <- 0.1
    
    df_scatter <-
      dat %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      slice_sample(prop = 1) %>%
      group_by(cluster) %>%
      
      mutate(
        x_delta = 1:n(),
        x_delta = scales::rescale(x_delta, to = c(dodge_width, 1 - dodge_width)),
        xx = x_base - 1 + x_delta,
        y = ifelse(avg_log2FC > 0, avg_log2FC + rect_height, avg_log2FC - rect_height)
      ) %>%
      mutate(
        group = ifelse(p_val_adj > 0.01, "p > 0.01", "p < 0.01"),
        group = factor(group, levels = c("p > 0.01", "p < 0.01"))
      ) %>%
      dplyr::arrange(group)
    
    # df_label 
    df_label_up <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_max(avg_log2FC, n = 10)
    
    df_label_down <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_min(avg_log2FC, n = 10)
    
    
    
    # grp rect  
    df_grp_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise() %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1, xmax = x_base) %>%
      mutate(ymin = -rect_height, ymax = rect_height)
    
    
    # background rect
    df_bg_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise(ymax = max(avg_log2FC) + rect_height, ymin = min(avg_log2FC) - rect_height) %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1 + dodge_width, xmax = x_base - dodge_width)
    
    
    
    df_scatter %>%
      ggplot(aes(x = xx, y = y)) +
      geom_rect(
        data = df_bg_rect, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        inherit.aes = FALSE
      ) +
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_point(aes(color = group), size = 1) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text.x = element_blank()
      ) +
      labs(
        x = NULL,
        y = "Log2FC",
        title = name
      )
    
  })

p_volcanoplot_major_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path('result/ENDO_subtype', str_glue("{name} diffgene dot.pdf")),
      height = 10,
      width = 20
    )
  })


## enrichment ----------------------------
endo_de_genes_for_enrichmenet <- 
  endo_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.15) %>%
  dplyr::filter(abs(avg_log2FC) > 0.5) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")


enrich_ls <- 
  endo_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichGO(
      x,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
      ont='ALL',
      pAdjustMethod = 'BH',
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.05,
      keyType = 'SYMBOL'
    )
  })


er_output<-'result/ENDO_subtype/enrichment' 
dir.create(er_output)

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(er_output, name))
    
    write.csv(
      x@result,
      file = file.path(er_output, name, "Go enrichment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
# export dotplot
enrich_dot_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "go dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })

endo_de_genes_id_for_enrichmenet <-
  endo_de_genes_for_enrichmenet %>%
  map(function(x){ # browser()
    clusterProfiler::bitr(
      x,
      OrgDb = "org.Hs.eg.db",
      fromType = "SYMBOL",
      toType = "ENTREZID"
    ) %>%
      .$ENTREZID
  }) %>%
  purrr::discard(function(x) length(x) < 1)


enrich_ls <-
  endo_de_genes_id_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichKEGG(
      x,
      organism = "hsa", 
      keyType = "ncbi-geneid", 
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
    )
  })

# saveRDS(enrich_ls, file.path(er_output, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(er_output, name))
    
    write.csv(
      x@result,
      file = file.path(er_output, name, "KEGG enriment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
# export dotplot
enrich_dot_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "kegg dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })

saveRDS(ENDO, 'result/ENDO_subtype/Endothelial.rds')


#EPI--------
EPI <- subset(xenium, annotation %in% c('AT1','AT2','Ciliated'))
EPI@images <- list()

## PCA --------
DefaultAssay(EPI) <- "Xenium"
EPI <- FindVariableFeatures(EPI)
EPI <- ScaleData(EPI, vars.to.regress = 'nCount_RNA')
EPI <- RunPCA(EPI)

EPI_out_dir<-'result/EPI_subtype/' 
dir.create(EPI_out_dir)


### pca inspector ------
pdf(file.path(EPI_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(EPI, ndims = 50)
dev.off()

pdf(file.path(EPI_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(EPI, dims = 1:5)
dev.off()
pdf(file.path(EPI_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(EPI, dims = 6:10)
dev.off()
pdf(file.path(EPI_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(EPI, dims = 11:15)
dev.off()
pdf(file.path(EPI_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(EPI, dims = 16:20)
dev.off()


## dimension reduction --------
dims <- 10

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


EPI <- RunUMAP(EPI, dims = 1:dims,return.model = TRUE)

pdf(file.path(EPI_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(EPI, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

EPI <- EPI %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

EPI <- RunUMAP(
  EPI, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(EPI_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(EPI, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'

### cluster -------
EPI  <-
  FindNeighbors(
    EPI ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


EPI <-
  FindClusters(
    EPI,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(EPI_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(EPI, group.by = "Xenium_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(EPI, group.by = "Xenium_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(EPI, group.by = "Xenium_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(EPI, group.by = "Xenium_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


### dotplot -------
EPI_major_marker_ls <- c(
  # AT2：肺泡 2 型细胞
  "AT2" = "Sftpc",
  "AT2" = "Sftpd",
  "AT2" = "Pgc",
  "AT2" = "Napsa",
  "AT2" ="Wif1", 
  "AT2" ="Etv5", 
  "AT2" ="Abca3",
  
  # Basal：基底细胞
  "Basal" = "Krt5",
  "Basal" = "Krt15",
  "Basal" = "Krt17",
  "Basal" = "Trp63",
  
  # AT1：肺泡 1 型细胞
  "AT1" = "Rtkn2",
  "AT1" = "Ager",
  "AT1" = "Col4a3",
  
  # Transitional AT2：过渡态肺泡 2 型细胞
  "Transitional AT2" = "Lamp3",
  
  # Ciliated：纤毛细胞
  "Ciliated" = "C20orf85",  # 保留原名（无标准小鼠同源基因）
  "Ciliated" = "Foxj1",
  "Ciliated" = "Scgb1a1",
  "Ciliated" = "Trp73",
  
  # Secretory：分泌细胞
  "Secretory" = "Muc5b",
  "Secretory" = "Wfdc2",
  "Secretory" = "Mmp7",
  "Secretory" = "Scgb3a2",
  
  # PNEC：肺神经内分泌细胞
  "PNEC" = "Scg2",
  "PNEC" = "Chgb",
  "PNEC" = "Calca",
  "PNEC" = "Ascl1",
  
  # Mast：肥大细胞
  "Mast" = "Cpa3",
  "Mast" = "Tpsab1",
  "Mast" = "Kit",
  
  # Proliferating：增殖细胞
  "Proliferating" = "Top2a",
  "Proliferating" = "Mki67",
  "Proliferating" = "Cdk1",
  
  # 未明确分类细胞
  "Unclassified" = "Nkx2-1",
  "Unclassified" = "Bmp4",
  "Unclassified" = "Spink1"
)
valid_EPIgenes <- EPI_major_marker_ls[EPI_major_marker_ls %in% rownames(EPI)]
cat("有效基因:", length(valid_EPIgenes), "/", length(EPI_major_marker_ls))

# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_gene <- setdiff(EPI_major_marker_ls, valid_EPIgenes)
missing_gene


# find markers
Idents(EPI) <- "Xenium_snn_res.0.4"
EPI_markers <-
  FindAllMarkers(
    JoinLayers(EPI),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

EPI_top10_markers <- 
  EPI_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10) %>%
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

names(EPI_top10_markers) <- sub('.*\\.','',names(EPI_top10_markers))


#提取deepseek
EPI_top10_markers1 <-
  EPI_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10)


# dotplot
pdf(file.path(EPI_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 35, height = 6)
DotPlot(EPI, group.by = "Xenium_snn_res.0.4", features = major_marker)+RotatedAxis()
DotPlot(EPI, group.by = "Xenium_snn_res.0.4", features = EPI_major_marker_ls)+RotatedAxis()
DotPlot(EPI, group.by = "Xenium_snn_res.0.4", features = EPI_top10_markers)+RotatedAxis()
dev.off()


for (gene in valid_EPIgenes) {
  p <- FeaturePlot(EPI, features = gene, reduction = 'harmonyumap')+RotatedAxis()
  ggsave(
    filename = file.path(EPI_out_dir, paste0("featureplot/",gene, "_feature.png")),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}


## annotation -----------
cell_cluster <- c(
  '0' = "AT2", # SFTPC SFTPB SFTPD
  '1' = "Ciliated", # FOXJ1 TPPP3 DNAH5
  '2' = "AT1", # AGER PDPN CLDN4
  '3' = "AT2", # SFTPC SFTPB SFTPD
  '4' = "Ciliated", # FOXJ1 TPPP3 DNAH5
  '5' = "AT1", # AGER PDPN CLDN4
  '6' = "AT2", # SFTPC SFTPB SFTPD
  '7' = "AT1", # AGER PDPN CLDN4
  '8' = "AT1", # AGER PDPN CLDN4
  '9' = "AT2", # SFTPC SFTPB SFTPD
  '10' = "Inf_Epi", # IFITM3 ISG15 IFI44L
  '11' = "delete" # no marker
)

Idents(EPI) <- "Xenium_snn_res.0.4"

EPI <- RenameIdents(EPI, cell_cluster)
EPI$subtype <- Idents(EPI)

pdf(file.path(EPI_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(EPI, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(EPI_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 7, width = 18)
DimPlot(EPI, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


## cell props barplot  --------
EPI_meta <- EPI[[]]

EPI_meta$sample <- factor(EPI_meta$sample, 
                          levels=c(EPI_meta$sample %>% unique())
)
p_col_proj <-
  EPI_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path("result/EPI_subtype/cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  EPI_meta %>%
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
  filename = file.path("result/EPI_subtype/cell_freq_group.pdf"),
  height = 6,
  width  = 4
)


## volcano-------
EPI_markers <-
  FindAllMarkers(
    JoinLayers(EPI),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(EPI_markers,'result/EPI_subtype/EPI_diff_gene.csv')

library(EnhancedVolcano)
type<-unique(EPI_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- EPI_markers[EPI_markers$cluster == type, ]
    p<-EnhancedVolcano(markers,
                       lab = rownames(markers),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = str_glue('{name} volcano plot'),
                       pCutoff = 10e-2,
                       FCcutoff = 0,
                       pointSize = 3.0,
                       labSize = 6.0)
    p
  })
volcano_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path('result/EPI_subtype', str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })

EPI_markers_ls <- list(EPI_markers)
names(EPI_markers_ls) <- 'Epithelial'

# scatter
p_volcanoplot_major_ls <- 
  EPI_markers_ls %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 0.5
    dodge_width <- 0.1
    
    df_scatter <-
      dat %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      slice_sample(prop = 1) %>%
      group_by(cluster) %>%
      
      mutate(
        x_delta = 1:n(),
        x_delta = scales::rescale(x_delta, to = c(dodge_width, 1 - dodge_width)),
        xx = x_base - 1 + x_delta,
        y = ifelse(avg_log2FC > 0, avg_log2FC + rect_height, avg_log2FC - rect_height)
      ) %>%
      mutate(
        group = ifelse(p_val_adj > 0.01, "p > 0.01", "p < 0.01"),
        group = factor(group, levels = c("p > 0.01", "p < 0.01"))
      ) %>%
      dplyr::arrange(group)
    
    # df_label 
    df_label_up <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_max(avg_log2FC, n = 10)
    
    df_label_down <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_min(avg_log2FC, n = 10)
    
    
    
    # grp rect  
    df_grp_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise() %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1, xmax = x_base) %>%
      mutate(ymin = -rect_height, ymax = rect_height)
    
    
    # background rect
    df_bg_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise(ymax = max(avg_log2FC) + rect_height, ymin = min(avg_log2FC) - rect_height) %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1 + dodge_width, xmax = x_base - dodge_width)
    
    
    
    df_scatter %>%
      ggplot(aes(x = xx, y = y)) +
      geom_rect(
        data = df_bg_rect, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        inherit.aes = FALSE
      ) +
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_point(aes(color = group), size = 1) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text.x = element_blank()
      ) +
      labs(
        x = NULL,
        y = "Log2FC",
        title = name
      )
    
  })

p_volcanoplot_major_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path('result/EPI_subtype', str_glue("{name} diffgene dot.pdf")),
      height = 10,
      width = 20
    )
  })


## enrichment ----------------------------
EPI_de_genes_for_enrichmenet <- 
  EPI_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.15) %>%
  dplyr::filter(abs(avg_log2FC) > 0.5) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")


enrich_ls <- 
  EPI_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichGO(
      x,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
      ont='ALL',
      pAdjustMethod = 'BH',
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.05,
      keyType = 'SYMBOL'
    )
  })


er_output<-'result/EPI_subtype/enrichment' 
dir.create(er_output)

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(er_output, name))
    
    write.csv(
      x@result,
      file = file.path(er_output, name, "Go enrichment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
# export dotplot
enrich_dot_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "go dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })

EPI_de_genes_id_for_enrichmenet <-
  EPI_de_genes_for_enrichmenet %>%
  map(function(x){ # browser()
    clusterProfiler::bitr(
      x,
      OrgDb = "org.Hs.eg.db",
      fromType = "SYMBOL",
      toType = "ENTREZID"
    ) %>%
      .$ENTREZID
  }) %>%
  purrr::discard(function(x) length(x) < 1)


enrich_ls <-
  EPI_de_genes_id_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichKEGG(
      x,
      organism = "hsa", 
      keyType = "ncbi-geneid", 
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
    )
  })

# saveRDS(enrich_ls, file.path(er_output, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(er_output, name))
    
    write.csv(
      x@result,
      file = file.path(er_output, name, "KEGG enriment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
# export dotplot
enrich_dot_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "kegg dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })

saveRDS(EPI, 'result/EPI_subtype/Epithelial.rds')



#T_cell--------
all_data2$annotation %>% table()
T_cell <- subset(all_data2, annotation %in% c('T_cell'))
T_cell@images <- list()

## PCA --------
DefaultAssay(T_cell) <- "RNA"
T_cell <- FindVariableFeatures(T_cell)
T_cell <- ScaleData(T_cell, vars.to.regress = 'nCount_RNA')
T_cell <- RunPCA(T_cell)

T_cell_out_dir<-'result/T_cell_subtype/' 
dir.create(T_cell_out_dir)


### pca inspector ------
pdf(file.path(T_cell_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(T_cell, ndims = 50)
dev.off()

# pdf(file.path(T_cell_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
# DimHeatmap(T_cell, dims = 1:5)
# dev.off()
# pdf(file.path(T_cell_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
# DimHeatmap(T_cell, dims = 6:10)
# dev.off()
# pdf(file.path(T_cell_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
# DimHeatmap(T_cell, dims = 11:15)
# dev.off()
# pdf(file.path(T_cell_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
# DimHeatmap(T_cell, dims = 16:20)
# dev.off()


## dimension reduction --------
dims <- 15

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


T_cell <- RunUMAP(T_cell, dims = 1:dims,return.model = TRUE)

pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(T_cell, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

T_cell <- T_cell %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

T_cell <- RunUMAP(
  T_cell, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(T_cell, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'

### cluster -------
T_cell  <-
  FindNeighbors(
    T_cell ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


T_cell <-
  FindClusters(
    T_cell,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(T_cell, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(T_cell, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(T_cell, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(T_cell, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


### dotplot -------
T_cell_major_marker_ls <- c(
  "CD4",
  "CD8A",
  "CD8B",
  
  # Memory（记忆性T细胞）
  "Memory" = "CCR7",
  "Memory" = "TCF7",
  "Memory" = "CD69",
  "Memory" = "NR4A1",
  "Memory" = "MYADM",
  "Memory" = "GATA3",
  "Memory" = "TBX21",
  
  # Effector（效应T细胞）
  "Effector" = "BATF3",
  "Effector" = "IL2RA",
  "Effector" = "IFNG",
  "Effector" = "IRF4",
  "Effector" = "MYC",
  "Effector" = "SLC7A5",
  "Effector" = "SLC7A1",
  "Effector" = "XCL1",
  "Effector" = "GZMB",
  "Effector" = "CCL3",
  "Effector" = "CCL4",
  "Effector" = "IL2",
  "Effector" = "PRF1",
  "Effector" = "NKG7",
  "Effector" = "GNLY",
  
  # Tolerant（耐受性T细胞）
  "Tolerant" = "EGR2",
  "Tolerant" = "IKZF2",
  "Tolerant" = "IZUMO1R",
  "Tolerant" = "CD200",
  "Tolerant" = "DGKZ",
  "Tolerant" = "BTLA",
  "Tolerant" = "TOX",
  "Tolerant" = "CTLA4",
  
  # Exhausted（耗竭性T细胞）
  "Exhausted" = "PDCD1",
  "Exhausted" = "LAG3",
  "Exhausted" = "HAVCR2",
  "Exhausted" = "CD244",
  "Exhausted" = "CD160",
  "Exhausted" = "IL10R",
  "Exhausted" = "EOMES",
  "Exhausted" = "NR4A2",
  "Exhausted" = "PTGER4",
  "Exhausted" = "TOX",
  "Exhausted" = "TOX2",
  "Exhausted" = "TIGIT",
  "Exhausted" = "CTLA4",
  "Exhausted" = "ENTPD1",
  
  # Proliferation（增殖相关）
  "Proliferation" = "MKI67",
  "Proliferation" = "TK1",
  "Proliferation" = "STMN1",
  
  # Regulatory（调节性T细胞）
  "Regulatory" = "FOXP3",
  "Regulatory" = "IL2RA",
  "Regulatory" = "CTLA4",
  "Regulatory" = "TIGIT",
  "Regulatory" = "TNFRSF18",
  "Regulatory" = "IL10",
  "Regulatory" = "IKZF2",
  "Regulatory" = "CCR8",
  
  # Early_Memory（早期记忆T细胞）
  "Early_Memory" = "CCR2",
  "Early_Memory" = "CX3CR1",
  "Early_Memory" = "IL1RRAP",
  "Early_Memory" = "ZEB2",
  
  # Naive（初始T细胞）
  "Naive" = "IL7R",
  "Naive" = "CCR7",
  "Naive" = "SELL",
  "Naive" = "FOXO1",
  "Naive" = "KLF2",
  "Naive" = "KLF3",
  "Naive" = "LEF1",
  "Naive" = "TCF7",
  "Naive" = "ACTN1",
  "Naive" = "FOXP1",
  "Naive" = "GZMK",
  
  # Activation_Effector_function（活化效应功能）
  "Activation_Effector_function" = "FAS",
  "Activation_Effector_function" = "FASLG",
  "Activation_Effector_function" = "CD44",
  "Activation_Effector_function" = "CD69",
  "Activation_Effector_function" = "CD38",
  "Activation_Effector_function" = "KLRF1",
  "Activation_Effector_function" = "FOXP1",
  "Activation_Effector_function" = "KLRG1",
  "Activation_Effector_function" = "FCGR3A",
  "Activation_Effector_function" = "CX3CR1",
  "Activation_Effector_function" = "SELL",
  "Activation_Effector_function" = "CD40LG",
  "Activation_Effector_function" = "ANXA1",
  "Activation_Effector_function" = "IL2RA",
  
  # Exhaustion（耗竭相关，补充）
  "PD-1",
  "Exhaustion" = "PDCD1",
  "Exhaustion" = "LAYN",
  "Exhaustion" = "HAVCR2",
  "Exhaustion" = "LAG3",
  "Exhaustion" = "CTLA4",
  "Exhaustion" = "TOX",
  "Exhaustion" = "VSIR",
  "Exhaustion" = "BTLA",
  "Exhaustion" = "ENTPD1",
  "Exhaustion" = "CD160",
  
  # TCR_Signaling（TCR信号通路）
  "TCR_Signaling" = "CALM1",
  "TCR_Signaling" = "CALM2",
  "TCR_Signaling" = "CALM3",
  "TCR_Signaling" = "CAST",
  "TCR_Signaling" = "GNLY",
  "TCR_Signaling" = "PRF1",
  "TCR_Signaling" = "CD3D",
  "TCR_Signaling" = "CD3E",
  "TCR_Signaling" = "CD3G",
  "TCR_Signaling" = "CSK",
  
  # Cytotoxicity（细胞毒性）
  "Cytotoxicity" = "GZMA",
  "Cytotoxicity" = "GZMB",
  "Cytotoxicity" = "GZMH",
  "Cytotoxicity" = "GZMK",
  "Cytotoxicity" = "NKG7",
  "Cytotoxicity" = "PRF1",
  "Cytotoxicity" = "KLRF1",
  "Cytotoxicity" = "KLRK1",
  "Cytotoxicity" = "FCGR3A",
  "Cytotoxicity" = "IFNG",
  "Cytotoxicity" = "GNLY",
  
  # Cytokine_Cytokine_receptor（细胞因子-细胞因子受体）
  "Cytokine_Cytokine_receptor" = "CSF1",
  "Cytokine_Cytokine_receptor" = "IL10RA",
  "Cytokine_Cytokine_receptor" = "IL16",
  "Cytokine_Cytokine_receptor" = "IL17RA",
  "Cytokine_Cytokine_receptor" = "IL18RAP",
  "Cytokine_Cytokine_receptor" = "IL21R",
  "Cytokine_Cytokine_receptor" = "IL2RB",
  "Cytokine_Cytokine_receptor" = "IL2RG",
  "Cytokine_Cytokine_receptor" = "IL32",
  "Cytokine_Cytokine_receptor" = "IL9R",
  
  # Chemokine_Chemokine_receptor（趋化因子-趋化因子受体）
  "Chemokine_Chemokine_receptor" = "CCR4",
  "Chemokine_Chemokine_receptor" = "CCR5",
  "Chemokine_Chemokine_receptor" = "CCR7",
  "Chemokine_Chemokine_receptor" = "CXCR3",
  "Chemokine_Chemokine_receptor" = "CXCR4",
  "Chemokine_Chemokine_receptor" = "CXCR5",
  "Chemokine_Chemokine_receptor" = "CXCR6",
  "Chemokine_Chemokine_receptor" = "CXCL13",
  "Chemokine_Chemokine_receptor" = "CXCL8",
  "Chemokine_Chemokine_receptor" = "CCL3",
  
  # Senescence（衰老相关）
  "Senescence" = "NT5E",
  "Senescence" = "IZUMO1R",
  "Senescence" = "LAG3",
  "Senescence" = "CBLB",
  "Senescence" = "DOK1",
  "Senescence" = "DOK2",
  "Senescence" = "SERPINB1",
  "Senescence" = "SERPINB6",
  "Senescence" = "SERPINB9",
  "Senescence" = "ADAM10",
  
  # Anergy（无反应性）
  "Anergy" = "NFKB1",
  "Anergy" = "NFKB2",
  "Anergy" = "NFKBIA",
  "Anergy" = "NFKBIB",
  "Anergy" = "NFKBIZ",
  "Anergy" = "CHUK",
  "Anergy" = "IKBKB",
  "Anergy" = "IKBKG",
  "Anergy" = "RELA",
  "Anergy" = "RELB",
  
  # NFKB_Signaling（NFKB信号通路）
  "NFKB_Signaling" = "HIKESHI",
  "NFKB_Signaling" = "UBA52",
  "NFKB_Signaling" = "RPA2",
  "NFKB_Signaling" = "RPA3",
  "NFKB_Signaling" = "MAP2K1",
  "NFKB_Signaling" = "MAP2K2",
  "NFKB_Signaling" = "MAP2K3",
  "NFKB_Signaling" = "MAP3K4",
  "NFKB_Signaling" = "MAP3K5",
  "NFKB_Signaling" = "MAP3K8",
  
  # Stress_response（应激反应）
  "Stress_response" = "MAP2K1",
  "Stress_response" = "MAP2K2",
  "Stress_response" = "MAP2K3",
  "Stress_response" = "MAP3K4",
  "Stress_response" = "MAP3K5",
  "Stress_response" = "MAP3K8",
  "Stress_response" = "MAPK1",
  "Stress_response" = "MAPK11",
  "Stress_response" = "MAPK13",
  
  # MAPK_Signaling（MAPK信号通路）
  "MAPK_Signaling" = "ITGA1",
  "MAPK_Signaling" = "ITGA4",
  "MAPK_Signaling" = "ITGAE",
  "MAPK_Signaling" = "ITGAL",
  "MAPK_Signaling" = "ITGAM",
  "MAPK_Signaling" = "ITGB1",
  "MAPK_Signaling" = "ITGB2",
  "MAPK_Signaling" = "ITGB7",
  "MAPK_Signaling" = "SELL",
  "MAPK_Signaling" = "SELPLG",
  
  # Adhesion（黏附相关）
  "Adhesion" = "IFIT1",
  "Adhesion" = "IFIT2",
  "Adhesion" = "IFIT3",
  "Adhesion" = "IFIT5",
  "Adhesion" = "STAT1",
  "Adhesion" = "STAT2",
  "Adhesion" = "IRF1",
  "Adhesion" = "IRF4",
  "Adhesion" = "IRF7",
  "Adhesion" = "IRF8",
  
  # IFN_Response（IFN应答）
  "IFN_Response" = "ATP1B3",
  "IFN_Response" = "ATP2A3",
  "IFN_Response" = "ATP2B1",
  "IFN_Response" = "ATP2B4",
  "IFN_Response" = "ATP5A1",
  "IFN_Response" = "ATP5B",
  "IFN_Response" = "ATP5C1",
  "IFN_Response" = "ATP5D",
  "IFN_Response" = "ATP5E",
  "IFN_Response" = "ATP5EP2",
  
  # Oxidative_phosphorylation（氧化磷酸化）
  "Oxidative_phosphorylation" = "SLC2A3",
  "Oxidative_phosphorylation" = "SLC2A8",
  "Oxidative_phosphorylation" = "PFKFB3",
  "Oxidative_phosphorylation" = "ALDOA",
  "Oxidative_phosphorylation" = "PGAM1",
  "Oxidative_phosphorylation" = "PGK1",
  "Oxidative_phosphorylation" = "PKM",
  "Oxidative_phosphorylation" = "LDHA",
  "Oxidative_phosphorylation" = "LDHB",
  
  # Glycolysis（糖酵解）
  "Glycolysis" = "ACAA1",
  "Glycolysis" = "ACADM",
  "Glycolysis" = "ACADVL",
  "Glycolysis" = "ACSL5",
  "Glycolysis" = "ALDH16A1",
  "Glycolysis" = "ALDH9A1",
  "Glycolysis" = "BAK1",
  "Glycolysis" = "BAX",
  "Glycolysis" = "BID",
  
  # Fatty_acid_metabolism（脂肪酸代谢）
  "Fatty_acid_metabolism" = "BCL2L11",
  "Fatty_acid_metabolism" = "BCL11B",
  "Fatty_acid_metabolism" = "BCL2L1",
  "Fatty_acid_metabolism" = "BIN1",
  "Fatty_acid_metabolism" = "BIN2",
  "Fatty_acid_metabolism" = "BAG3",
  "Fatty_acid_metabolism" = "BAG4",
  "Fatty_acid_metabolism" = "BAK1",
  "Fatty_acid_metabolism" = "BAX",
  
  # Pro_apoptosis（促凋亡）
  "Pro_apoptosis" = "BCL2L11",
  "Pro_apoptosis" = "BCL11B",
  "Pro_apoptosis" = "BCL2L1",
  "Pro_apoptosis" = "BIN1",
  "Pro_apoptosis" = "BIN2",
  "Pro_apoptosis" = "BAG3",
  "Pro_apoptosis" = "BAG4",
  "Pro_apoptosis" = "BAK1",
  "Pro_apoptosis" = "BAX",
  
  # Anti_apoptosis（抗凋亡）
  "Anti_apoptosis" = "BCL2L1",
  "Anti_apoptosis" = "BCL11B",
  "Anti_apoptosis" = "BIN1",
  "Anti_apoptosis" = "BIN2",
  "Anti_apoptosis" = "BAG3",
  "Anti_apoptosis" = "BAG4",
  "Anti_apoptosis" = "BAK1",
  "Anti_apoptosis" = "BAX"
)


valid_T_cell_major_marker_ls <- T_cell_major_marker_ls[T_cell_major_marker_ls %in% rownames(T_cell)]
cat("有效基因:", length(valid_T_cell_major_marker_ls), "/", length(T_cell_major_marker_ls))

# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_gene <- setdiff(T_cell_major_marker_ls, valid_T_cell_major_marker_ls)
missing_gene

valid_T_cell_major_marker_lsunique <- valid_T_cell_major_marker_ls[!duplicated(valid_T_cell_major_marker_ls)]

# find markers
Idents(T_cell) <- "RNA_snn_res.0.4"
T_cell_markers <-
  FindAllMarkers(
    JoinLayers(T_cell),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

T_cell_top10_markers <- 
  T_cell_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10) %>%
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

names(T_cell_top10_markers) <- sub('.*\\.','',names(T_cell_top10_markers))


#提取deepseek
T_cell_top10_markers1 <-
  T_cell_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10)


# dotplot
pdf(file.path(T_cell_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 40, height = 6)
DotPlot(T_cell, group.by = "RNA_snn_res.0.4", features = major_marker)+RotatedAxis()
DotPlot(T_cell, group.by = "RNA_snn_res.0.4", features = valid_T_cell_major_marker_lsunique)+RotatedAxis()
# DotPlot(T_cell, group.by = "RNA_snn_res.0.2", features = immune_major_marker_ls)+RotatedAxis()
DotPlot(T_cell, group.by = "RNA_snn_res.0.4", features = T_cell_top10_markers)+RotatedAxis()
dev.off()

#valid_T_cell_major_marker_lsunique <- c("CD4","CD8A","CD8B")

lapply(valid_T_cell_major_marker_lsunique, function(gene) {
  p <- FeaturePlot(T_cell, features = gene, reduction = "harmonyumap") + RotatedAxis()
  ggsave(
    filename = file.path(T_cell_out_dir, "featureplot", paste0(gene, "_feature.png")),
    plot = p,
    width = 8, height = 6, dpi = 300
  )
})





##空间imagefeatureplot-------
gene_list <- c('CD3E',"CD4","CD8A")
valid_gene <- gene_list[gene_list %in% rownames(xenium)]
cat("有效基因:", length(valid_gene), "/", length(gene_list))

missing_gene <- setdiff(gene_list, valid_gene)

source("ImageFeature-plot.R") 

Imagefeature_dir <- 'result/T_cell_subtype/Imagefeature/'
dir.create(Imagefeature_dir)

all_colors <- c(
  "#A6CEE3", "#5B9EC9", "#2D82AF", "#7EBA98", "#98D277",
  "#52AF43", "#6F9E4C", "#DD9A88", "#FDBB69", "#FE982C",
  "#F78620", "#D9A295", "#feee7d", "#02c9c9", "#600473",
  "#B294C7", "#7D54A5", "#9E8099", "#F0EB99", "#DBB466",
  "#B15928", # 原始21色
  "#FF6B6B", "#4ECDC4", "#FF9F1C", "#7209B7", "#4361EE", 
  "#F15BB5", "#00BBF9", "#00F5D4", "#9B5DE5", "#00F5A7",
  "#FEE440", "#9C89B8", "#06D6A0"  # 新增13色
)

all_colors8 <- c(
  "#FF6B6B","#A6CEE3",
  "#F15BB5","#00F5A7", "#FDBB69","#feee7d", #"#9E8099",
  "#00BBF9","#7209B7"#,"#9B5DE5"
)
# all_colors10 <- c(
#   "gray90","#F15BB5","gray90","gray90","gray90","gray90","gray90","gray90","gray90","gray90"
# )

Idents(all_data2) <- "annotation"
names(all_colors8) <- levels(all_data2)

###所有细胞类型------
for (sample_id in unique(xenium@meta.data$sample)){
  for (gene in valid_gene) {
    
    P2 <- ImageFeaturePlot_wyt(xenium, 
                               features = gene,
                               size = 0.75, 
                               cols = c('#003366', '#FFD700'),
                               fov =sample_id,
                               nmols = 20000,
                               border.color = all_colors8, 
                               border.size = 0.2,
                               border.by = "celltype",
                               dark.background = F
    )
    
    ggsave(P2,filename = paste0(Imagefeature_dir,sample_id,"_",gene,"_Imagefeature.pdf"),width = 25,height = 25,dpi = 300)
    ggsave(P2,filename = paste0(Imagefeature_dir,sample_id,"_",gene,"_Imagefeature.png"),width = 25,height = 25,dpi = 300)
    
    message("成功保存: ", sample_id, "_", gene)
  }
}


###只画cells.use 细胞--------
for (sample_id in unique(xenium@meta.data$sample)){
  for (gene in valid_gene) {
    
    P4 <- ImageFeaturePlot_wyt(xenium, 
                               features = gene, 
                               size = 0.5, 
                               cols = c('#003366', '#FFD700'),
                               fov =sample_id,
                               nmols = 20000,
                               border.color = all_colors8, 
                               border.size = 0.2,
                               border.by = "celltype",
                               cells.use = c("Endothelial")
    )
    ggsave(P4,filename = paste0(Imagefeature_dir,sample_id,"_",gene,"_Imagefeature_ENDO.pdf"),width = 25,height = 25,dpi = 300)
    ggsave(P4,filename = paste0(Imagefeature_dir,sample_id,"_",gene,"_Imagefeature_ENDO.png"),width = 25,height = 25,dpi = 300)
    
  }
}




## annotation -----------
cell_cluster <- c(
  '0' = '0', 
  '1' = '1',
  '2' = '2',
  '3' = '3',
  '4' = '4',
  '5' = '5',
  '6' = '6',
  '7' = '7',
  '8' = '8',
  '9' = '9',
  '10' = '10',
  '11' = '11'
)

Idents(T_cell) <- "RNA_snn_res.0.4"

T_cell <- RenameIdents(T_cell, cell_cluster)
T_cell$subtype <- Idents(T_cell)

pdf(file.path(T_cell_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(T_cell, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(T_cell_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 7, width = 18)
DimPlot(T_cell, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()





## cell props barplot  --------
T_cell_meta <- T_cell[[]]

T_cell_meta$sample <- factor(T_cell_meta$sample, 
                               levels=c(T_cell_meta$sample %>% unique())
)
p_col_proj <-
  T_cell_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path("result/T_cell_subtype/cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  T_cell_meta %>%
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
  filename = file.path("result/T_cell_subtype/cell_freq_group.pdf"),
  height = 6,
  width  = 4
)


## volcano-------
T_cell_markers <-
  FindAllMarkers(
    JoinLayers(T_cell),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(T_cell_markers,'result/T_cell_subtype/T_cell_diff_gene.csv')

library(EnhancedVolcano)
type<-unique(T_cell_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- T_cell_markers[T_cell_markers$cluster == type, ]
    p<-EnhancedVolcano(markers,
                       lab = rownames(markers),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = str_glue('{name} volcano plot'),
                       pCutoff = 10e-2,
                       FCcutoff = 0,
                       pointSize = 3.0,
                       labSize = 6.0)
    p
  })
volcano_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path('result/T_cell_subtype', str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })

T_cell_markers_ls <- list(T_cell_markers)
names(T_cell_markers_ls) <- 'T_cell'

# scatter
p_volcanoplot_major_ls <- 
  T_cell_markers_ls %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 0.5
    dodge_width <- 0.1
    
    df_scatter <-
      dat %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      slice_sample(prop = 1) %>%
      group_by(cluster) %>%
      
      mutate(
        x_delta = 1:n(),
        x_delta = scales::rescale(x_delta, to = c(dodge_width, 1 - dodge_width)),
        xx = x_base - 1 + x_delta,
        y = ifelse(avg_log2FC > 0, avg_log2FC + rect_height, avg_log2FC - rect_height)
      ) %>%
      mutate(
        group = ifelse(p_val_adj > 0.01, "p > 0.01", "p < 0.01"),
        group = factor(group, levels = c("p > 0.01", "p < 0.01"))
      ) %>%
      dplyr::arrange(group)
    
    # df_label 
    df_label_up <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_max(avg_log2FC, n = 10)
    
    df_label_down <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_min(avg_log2FC, n = 10)
    
    
    
    # grp rect  
    df_grp_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise() %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1, xmax = x_base) %>%
      mutate(ymin = -rect_height, ymax = rect_height)
    
    
    # background rect
    df_bg_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise(ymax = max(avg_log2FC) + rect_height, ymin = min(avg_log2FC) - rect_height) %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1 + dodge_width, xmax = x_base - dodge_width)
    
    
    
    df_scatter %>%
      ggplot(aes(x = xx, y = y)) +
      geom_rect(
        data = df_bg_rect, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        inherit.aes = FALSE
      ) +
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_point(aes(color = group), size = 1) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text.x = element_blank()
      ) +
      labs(
        x = NULL,
        y = "Log2FC",
        title = name
      )
    
  })

p_volcanoplot_major_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path('result/T_cell_subtype', str_glue("{name} diffgene dot.pdf")),
      height = 10,
      width = 20
    )
  })


## enrichment ----------------------------
T_cell_de_genes_for_enrichmenet <- 
  T_cell_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.15) %>%
  dplyr::filter(abs(avg_log2FC) > 0.5) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")


enrich_ls <- 
  T_cell_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichGO(
      x,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
      ont='ALL',
      pAdjustMethod = 'BH',
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.05,
      keyType = 'SYMBOL'
    )
  })


er_output<-'result/T_cell_subtype/enrichment' 
dir.create(er_output)

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(er_output, name))
    
    write.csv(
      x@result,
      file = file.path(er_output, name, "Go enrichment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
# export dotplot
enrich_dot_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "go dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })

T_cell_de_genes_id_for_enrichmenet <-
  T_cell_de_genes_for_enrichmenet %>%
  map(function(x){ # browser()
    clusterProfiler::bitr(
      x,
      OrgDb = "org.Hs.eg.db",
      fromType = "SYMBOL",
      toType = "ENTREZID"
    ) %>%
      .$ENTREZID
  }) %>%
  purrr::discard(function(x) length(x) < 1)


enrich_ls <-
  T_cell_de_genes_id_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichKEGG(
      x,
      organism = "hsa", 
      keyType = "ncbi-geneid", 
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
    )
  })

# saveRDS(enrich_ls, file.path(er_output, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(er_output, name))
    
    write.csv(
      x@result,
      file = file.path(er_output, name, "KEGG enriment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
# export dotplot
enrich_dot_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "kegg dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })

saveRDS(T_cell, 'result/T_cell_subtype/T_cell.rds')




table(all_data2@meta.data$celltype, all_data2@meta.data$sample)
#MYELOID------
all_data2@meta.data$annotation %>% table()
MYELOID <- subset(all_data2, annotation %in% c('Macrophage'))
MYELOID@images <- list()

## PCA --------
DefaultAssay(MYELOID) <- "RNA"
MYELOID <- FindVariableFeatures(MYELOID)
MYELOID <- ScaleData(MYELOID, vars.to.regress = 'nCount_RNA')
MYELOID <- RunPCA(MYELOID)

MYELOID_out_dir<-'result/MYELOID_subtype/' 
dir.create(MYELOID_out_dir)


### pca------
pdf(file.path(MYELOID_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(MYELOID, ndims = 50)
dev.off()

pdf(file.path(MYELOID_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(MYELOID, dims = 1:5)
dev.off()
pdf(file.path(MYELOID_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(MYELOID, dims = 6:10)
dev.off()
pdf(file.path(MYELOID_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(MYELOID, dims = 11:15)
dev.off()
pdf(file.path(MYELOID_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(MYELOID, dims = 16:20)
dev.off()


## dimension reduction --------
dims <- 15

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


MYELOID <- RunUMAP(MYELOID, dims = 1:dims,return.model = TRUE)

pdf(file.path(MYELOID_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(MYELOID, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

MYELOID <- MYELOID %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

MYELOID <- RunUMAP(
  MYELOID, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(MYELOID_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(MYELOID, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'

### cluster -------
MYELOID  <-
  FindNeighbors(
    MYELOID ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


MYELOID <-
  FindClusters(
    MYELOID,
    resolution = c(0.2, 0.4, 0.6)
  )

# dimplot
pdf(file.path(MYELOID_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(MYELOID, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(MYELOID, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(MYELOID, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
#DimPlot(MYELOID, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


### dotplot -------
MYELOID_major_marker_ls <- c(
  # M1型巨噬细胞极化标记
  "M1 Macrophage" = "IL12",
  "M1 Macrophage" = "IL23",
  "M1 Macrophage" = "TNF",
  "M1 Macrophage" = "IL6",
  "M1 Macrophage" = "CD86",
  "M1 Macrophage" = "MHCII",
  "M1 Macrophage" = "IL1B",
  "M1 Macrophage" = "MARCO",
  "M1 Macrophage" = "NOS2",
  "M1 Macrophage" = "CD64",
  "M1 Macrophage" = "CD80",
  "M1 Macrophage" = "CXCR10",
  "M1 Macrophage" = "CXCL9",
  "M1 Macrophage" = "CXCL10",
  "M1 Macrophage" = "CXCL11",
  "M1 Macrophage" = "IL1A",
  "M1 Macrophage" = "CCL5",
  "M1 Macrophage" = "IRF5",
  "M1 Macrophage" = "IRF1",
  "M1 Macrophage" = "CD40",
  "M1 Macrophage" = "IDO1",
  "M1 Macrophage" = "KYNU",
  "M1 Macrophage" = "CCR7",
  "M1 Macrophage" = "PTPRC",
  "M1 Macrophage" = "CD68",
  "M1 Macrophage" = "CD115",
  "M1 Macrophage" = "HLA-DRA",
  "M1 Macrophage" = "CD205",
  "M1 Macrophage" = "CD14",
  
  # M2型巨噬细胞极化标记
  "M2 Macrophage" = "ARG1",
  "M2 Macrophage" = "ARG2",
  "M2 Macrophage" = "IL10",
  "M2 Macrophage" = "CD32",
  "M2 Macrophage" = "CD163",
  "M2 Macrophage" = "CD23",
  "M2 Macrophage" = "CD200R1",
  "M2 Macrophage" = "PDCD1LG2",
  "M2 Macrophage" = "CD274",
  "M2 Macrophage" = "CSF1R",
  "M2 Macrophage" = "CD206",
  "M2 Macrophage" = "IL1RN",
  "M2 Macrophage" = "IL1R2",
  "M2 Macrophage" = "IL4R",
  "M2 Macrophage" = "CCL4",
  "M2 Macrophage" = "CCL13",
  "M2 Macrophage" = "CCL20",
  "M2 Macrophage" = "CCL17",
  "M2 Macrophage" = "CCL18",
  "M2 Macrophage" = "CCL22",
  "M2 Macrophage" = "CCL24",
  "M2 Macrophage" = "LYVE1",
  "M2 Macrophage" = "VEGFA",
  "M2 Macrophage" = "VEGFB",
  "M2 Macrophage" = "VEGFC",
  "M2 Macrophage" = "VEGFD",
  "M2 Macrophage" = "EGF",
  "M2 Macrophage" = "CTSA",
  "M2 Macrophage" = "CTSB",
  "M2 Macrophage" = "CSTC",
  "M2 Macrophage" = "CTSD",
  "M2 Macrophage" = "TGFB1",
  "M2 Macrophage" = "TGFB2",
  "M2 Macrophage" = "TGFB3",
  "M2 Macrophage" = "MMP14",
  "M2 Macrophage" = "MMP19",
  "M2 Macrophage" = "MMP9",
  "M2 Macrophage" = "CLEC7A",
  "M2 Macrophage" = "WNT7B",
  "M2 Macrophage" = "FASL",
  "M2 Macrophage" = "TNFSF12",
  "M2 Macrophage" = "TNFSF8",
  "M2 Macrophage" = "CD276",
  "M2 Macrophage" = "VTCN1",
  "M2 Macrophage" = "MSR1",
  "M2 Macrophage" = "FN1",
  "M2 Macrophage" = "IRF4",
  # "M2 Macrophage" = "PTPRC",
  # "M2 Macrophage" = "CD68",
  # "M2 Macrophage" = "CD115",
  # "M2 Macrophage" = "HLA-DRA",
  # "M2 Macrophage" = "CD205",
  # "M2 Macrophage" = "CD14",
  
  # Fcγ受体
  "Fcγ Receptor" = "FCGR1A",
  "Fcγ Receptor" = "FCGR2A",
  "Fcγ Receptor" = "FCGR2B",
  "Fcγ Receptor" = "FCGR3A",
  "Fcγ Receptor" = "FCGR3B",
  
  # 免疫检查点分子
  # "Immune Checkpoint" = "TGFB1",
  # "Immune Checkpoint" = "CD274",
  # "Immune Checkpoint" = "IL10",
  "Immune Checkpoint" = "LAG3",
  "Immune Checkpoint" = "TIGIT",
  "Immune Checkpoint" = "HAVCR2"#,
  # "Immune Checkpoint" = "VEGFA"
)
valid_MYELOIDgenes <- MYELOID_major_marker_ls[MYELOID_major_marker_ls %in% rownames(MYELOID)]
cat("有效基因:", length(valid_MYELOIDgenes), "/", length(MYELOID_major_marker_ls))

# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_gene <- setdiff(MYELOID_major_marker_ls, valid_MYELOIDgenes)
missing_gene



# find markers
Idents(MYELOID) <- "RNA_snn_res.0.6"
MYELOID_markers <-
  FindAllMarkers(
    JoinLayers(MYELOID),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

MYELOID_top10_markers <- 
  MYELOID_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10) %>%
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

names(MYELOID_top10_markers) <- sub('.*\\.','',names(MYELOID_top10_markers))


#提取deepseek
MYELOID_top10_markers1 <-
  MYELOID_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10)


# dotplot
pdf(file.path(MYELOID_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 35, height = 6)
DotPlot(MYELOID, group.by = "RNA_snn_res.0.6", features = major_marker)+RotatedAxis()
DotPlot(MYELOID, group.by = "RNA_snn_res.0.6", features = MYELOID_major_marker_ls)+RotatedAxis()
#DotPlot(MYELOID, group.by = "RNA_snn_res.0.4", features = immune_major_marker_ls)+RotatedAxis()
DotPlot(MYELOID, group.by = "RNA_snn_res.0.6", features = MYELOID_top10_markers)+RotatedAxis()
dev.off()


for (gene in valid_MYELOIDgenes) {
  p <- FeaturePlot(MYELOID, features = gene, reduction = 'harmonyumap')+RotatedAxis()
  ggsave(
    filename = file.path(MYELOID_out_dir, paste0("featureplot/",gene, "_feature.png")),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}


## annotation -----------
cell_cluster <- c(
  '0' = "M2", 
  '1' = "MIX",
  '2' = "MIX",
  '3' = "M2", 
  '4' = "M2", 
  '5' = "MIX",
  '6' = "Neutrophil", 
  '7' = "proliferate", 
  '8' = "M0", 
  '9' = "M1", 
  '10' = "Mast", 
  '11' = "DC" 
)

Idents(MYELOID) <- "RNA_snn_res.0.6"

MYELOID <- RenameIdents(MYELOID, cell_cluster)
MYELOID$subtype <- Idents(MYELOID)

pdf(file.path(MYELOID_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(MYELOID, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(MYELOID_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 7, width = 18)
DimPlot(MYELOID, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()





## cell props barplot  --------
MYELOID_meta <- MYELOID[[]]

MYELOID_meta$sample <- factor(MYELOID_meta$sample, 
                              levels=c(MYELOID_meta$sample %>% unique())
)
p_col_proj <-
  MYELOID_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path("result/MYELOID_subtype/cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  MYELOID_meta %>%
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
  filename = file.path("result/MYELOID_subtype/cell_freq_group.pdf"),
  height = 6,
  width  = 4
)

## 箱线图------
MYELOID_meta <- MYELOID[[]]

boxplot_data <-  MYELOID_meta %>%
  group_by(sample, group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))  

types<- unique(boxplot_data$subtype)

lapply(types, function(type){
  p <- boxplot_data %>%
    filter(subtype == type) %>% 
    ggplot(aes(x = group, y = freq, fill = group)) +
    geom_boxplot(width = 0.6) +
    #facet_wrap(~annotation, scales = "free_y", ncol = 4) +
    stat_compare_means(
      comparisons = list(c("Infla","Noninfla")),  
      method = "t.test",  
      label = "p.signif",  # p.signif, "p.value"
      p.adjust.method = "none"  
    ) +
    
    theme_bw() +
    labs(
      title = "",
      x = "group",
      y = "freq"
    ) +
    scale_fill_manual(values = c("#4DBBD5", "#E64B35"))  #
  
  print(p)
  
  ggsave(file.path(MYELOID_out_dir,paste0("group_boxplot/",type,"_group_boxplot.pdf")), p, width =5, height = 5)
  ggsave(file.path(MYELOID_out_dir,paste0("group_boxplot/",type,"_group_boxplot.png")), p, width =5, height = 5)
  
})




#MYELOID <- readRDS("result/MYELOID_subtype/MYELOID.rds")
## volcano-------
MYELOID_markers <-
  FindAllMarkers(
    JoinLayers(MYELOID),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(MYELOID_markers,'result/MYELOID_subtype/MYELOID_diff_gene.csv')

library(EnhancedVolcano)
type<-unique(MYELOID_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- MYELOID_markers[MYELOID_markers$cluster == type, ]
    p<-EnhancedVolcano(markers,
                       lab = rownames(markers),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = str_glue('{name} volcano plot'),
                       pCutoff = 10e-2,
                       FCcutoff = 0,
                       pointSize = 3.0,
                       labSize = 6.0)
    p
  })
volcano_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path('result/MYELOID_subtype', str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })

MYELOID_markers_ls <- list(MYELOID_markers)
names(MYELOID_markers_ls) <- 'MYELOID'

# scatter
p_volcanoplot_major_ls <- 
  MYELOID_markers_ls %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 0.5
    dodge_width <- 0.1
    
    df_scatter <-
      dat %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      slice_sample(prop = 1) %>%
      group_by(cluster) %>%
      
      mutate(
        x_delta = 1:n(),
        x_delta = scales::rescale(x_delta, to = c(dodge_width, 1 - dodge_width)),
        xx = x_base - 1 + x_delta,
        y = ifelse(avg_log2FC > 0, avg_log2FC + rect_height, avg_log2FC - rect_height)
      ) %>%
      mutate(
        group = ifelse(p_val_adj > 0.01, "p > 0.01", "p < 0.01"),
        group = factor(group, levels = c("p > 0.01", "p < 0.01"))
      ) %>%
      dplyr::arrange(group)
    
    # df_label 
    df_label_up <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_max(avg_log2FC, n = 10)
    
    df_label_down <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_min(avg_log2FC, n = 10)
    
    
    
    # grp rect  
    df_grp_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise() %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1, xmax = x_base) %>%
      mutate(ymin = -rect_height, ymax = rect_height)
    
    
    # background rect
    df_bg_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise(ymax = max(avg_log2FC) + rect_height, ymin = min(avg_log2FC) - rect_height) %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1 + dodge_width, xmax = x_base - dodge_width)
    
    
    
    df_scatter %>%
      ggplot(aes(x = xx, y = y)) +
      geom_rect(
        data = df_bg_rect, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        inherit.aes = FALSE
      ) +
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_point(aes(color = group), size = 1) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text.x = element_blank()
      ) +
      labs(
        x = NULL,
        y = "Log2FC",
        title = name
      )
    
  })

p_volcanoplot_major_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path('result/MYELOID_subtype', str_glue("{name} diffgene dot.pdf")),
      height = 10,
      width = 20
    )
  })


## enrichment ----------------------------
MYELOID_de_genes_for_enrichmenet <- 
  MYELOID_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.15) %>%
  dplyr::filter(abs(avg_log2FC) > 0.5) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")


enrich_ls <- 
  MYELOID_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichGO(
      x,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
      ont='ALL',
      pAdjustMethod = 'BH',
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.05,
      keyType = 'SYMBOL'
    )
  })


er_output<-'result/MYELOID_subtype/enrichment' 
dir.create(er_output)

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(er_output, name))
    
    write.csv(
      x@result,
      file = file.path(er_output, name, "Go enrichment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
# export dotplot
enrich_dot_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "go dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })




MYELOID_de_genes_id_for_enrichmenet <-
  MYELOID_de_genes_for_enrichmenet %>%
  map(function(x){ # browser()
    clusterProfiler::bitr(
      x,
      OrgDb = "org.Hs.eg.db",
      fromType = "SYMBOL",
      toType = "ENTREZID"
    ) %>%
      .$ENTREZID
  }) %>%
  purrr::discard(function(x) length(x) < 1)


enrich_ls <-
  MYELOID_de_genes_id_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichKEGG(
      x,
      organism = "hsa", 
      keyType = "ncbi-geneid", 
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
    )
  })

# saveRDS(enrich_ls, file.path(er_output, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(er_output, name))
    
    write.csv(
      x@result,
      file = file.path(er_output, name, "KEGG enriment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
# export dotplot
enrich_dot_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "kegg dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })

saveRDS(MYELOID, 'result/MYELOID_subtype/MYELOID.rds')

#MESENCHYMAL------
xenium$annotation %>% table()
MESENCHYMAL <- subset(xenium, annotation %in% c('SMC', "Pericyte", "Fibroblast", "Pleura"))
MESENCHYMAL@images <- list()

## PCA --------
DefaultAssay(MESENCHYMAL) <- "Xenium"
MESENCHYMAL <- FindVariableFeatures(MESENCHYMAL)
MESENCHYMAL <- ScaleData(MESENCHYMAL, vars.to.regress = 'nCount_RNA')
MESENCHYMAL <- RunPCA(MESENCHYMAL)

MESENCHYMAL_out_dir<-'result/MESENCHYMAL_subtype/' 
dir.create(MESENCHYMAL_out_dir)


### pca------
pdf(file.path(MESENCHYMAL_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(MESENCHYMAL, ndims = 50)
dev.off()

pdf(file.path(MESENCHYMAL_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(MESENCHYMAL, dims = 1:5)
dev.off()
pdf(file.path(MESENCHYMAL_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(MESENCHYMAL, dims = 6:10)
dev.off()
pdf(file.path(MESENCHYMAL_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(MESENCHYMAL, dims = 11:15)
dev.off()
pdf(file.path(MESENCHYMAL_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(MESENCHYMAL, dims = 16:20)
dev.off()


## dimension reduction --------
dims <- 13

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


MESENCHYMAL <- RunUMAP(MESENCHYMAL, dims = 1:dims,return.model = TRUE)

pdf(file.path(MESENCHYMAL_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(MESENCHYMAL, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

MESENCHYMAL <- MESENCHYMAL %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

MESENCHYMAL <- RunUMAP(
  MESENCHYMAL, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(MESENCHYMAL_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(MESENCHYMAL, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'

### cluster -------
MESENCHYMAL  <-
  FindNeighbors(
    MESENCHYMAL ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


MESENCHYMAL <-
  FindClusters(
    MESENCHYMAL,
    resolution = c(0.2, 0.3, 0.4#, 0.6, 0.8
    )
  )

# dimplot
pdf(file.path(MESENCHYMAL_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(MESENCHYMAL, group.by = "Xenium_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(MESENCHYMAL, group.by = "Xenium_snn_res.0.3", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(MESENCHYMAL, group.by = "Xenium_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
#DimPlot(MESENCHYMAL, group.by = "Xenium_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


### dotplot -------
# 小鼠免疫细胞主要 marker 列表
MESENCHYMAL_major_marker_ls <- c(
  # Alveolar FBs（肺泡成纤维细胞）
  "Alveolar FBs" = "Wnt5a",
  "Alveolar FBs" = "Fgf10",
  "Alveolar FBs" = "Fgf7",
  "Alveolar FBs" = "Fgf2",
  "Alveolar FBs" = "Sfrp2",
  "Alveolar FBs" = "Sfrp4",
  "Alveolar FBs" = "Pi16",
  
  # SMCs/fibrocytes FBs（平滑肌/周细胞样成纤维细胞）
  "SMCs/fibrocytes FBs" = "Acta2",
  "SMCs/fibrocytes FBs" = "Cspg4",
  
  # # Activated fibroblasts FBs（活化成纤维细胞）
  # "Activated fibroblasts FBs" = "Dcn",
  # "Activated fibroblasts FBs" = "Lum",
  # "Activated fibroblasts FBs" = "Col1a1",
  # "Activated fibroblasts FBs" = "Col1a2",
  "Activated fibroblasts FBs" = "Col3a1",
  # "Activated fibroblasts FBs" = "Fn1",
  # "Activated fibroblasts FBs" = "Fap",
  "Activated fibroblasts FBs" = "Cthrc1",
  
  # Inflammatory FBs（炎症性成纤维细胞）
  "Inflammatory FBs" = "H2-Aa",
  "Inflammatory FBs" = "Lyz2",
  
  # Proliferating FBs（增殖性成纤维细胞）
  "Proliferating FBs" = "Mki67",
  "Proliferating FBs" = "Top2a",
  "Proliferating FBs" = "Cdk1",
  "Proliferating FBs" = "Cenpf",
  
  # Adventitial FBs（血管外膜成纤维细胞）
  "Adventitial FBs" = "Dcn",
  "Adventitial FBs" = "Lum",
  "Adventitial FBs" = "Col1a1",
  "Adventitial FBs" = "Col1a2",
  "Adventitial FBs" = "Fn1",
  
  # Myofibroblasts FBs（肌成纤维细胞）
  # "Myofibroblasts FBs" = "Acta2",
  # "Myofibroblasts FBs" = "Fap",
  # "Myofibroblasts FBs" = "Dcn",
  # "Myofibroblasts FBs" = "Lum",
  # "Myofibroblasts FBs" = "Col1a1",
  # "Myofibroblasts FBs" = "Col1a2",
  # "Myofibroblasts FBs" = "Fn1",
  
  # Subpleural mesothelial（胸膜下间皮细胞）
  "Subpleural mesothelial" = "Msln",
  "Subpleural mesothelial" = "Has1",
  "Subpleural mesothelial" = "Krt8",
  "Subpleural mesothelial" = "Krt18"
)

valid_MESENCHYMALgenes <- MESENCHYMAL_major_marker_ls[MESENCHYMAL_major_marker_ls %in% rownames(MESENCHYMAL)]
cat("有效基因:", length(valid_MESENCHYMALgenes), "/", length(MESENCHYMAL_major_marker_ls))

# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_gene <- setdiff(MESENCHYMAL_major_marker_ls, valid_MESENCHYMALgenes)
missing_gene

MESENCHYMAL_major_marker_ls1 <- c(
  "Col1a1",
  
  "Activated fibroblasts FBs" = "Cthrc1",
  "Activated fibroblasts FBs" = "Postn",
  
  "Adventitial FBs" = "Mfap5",
  
  "Inflammatory FBs" = "Ptgsd",
  
  "Mesothelial" = "Msln",
  
  "Myofibroblasts FBs" = "Wnt5a",
  
  "SMCs/fibrocytes FBs" = "Acta2",
  "SMCs/fibrocytes FBs" = "Cspg4",
  
  "Subpleural mesothelial" = "Plin2"
)
valid_MESENCHYMALgenes1 <- MESENCHYMAL_major_marker_ls1[MESENCHYMAL_major_marker_ls1 %in% rownames(MESENCHYMAL)]
cat("有效基因:", length(valid_MESENCHYMALgenes1), "/", length(MESENCHYMAL_major_marker_ls1))

# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_gene <- setdiff(MESENCHYMAL_major_marker_ls1, valid_MESENCHYMALgenes1)
missing_gene



# 小鼠肺/组织免疫细胞精细亚型 marker 列表（严格按指定分类）
immune_major_marker_ls <- c(
  # Interstitial macrophages
  "Interstitial macrophages" = "Mrc1",
  "Interstitial macrophages" = "Lyve1",
  "Interstitial macrophages" = "F13a1",
  "Interstitial macrophages" = "Cd163",
  
  # CD4 cells
  "CD4 cells" = "Cd4",
  "CD4 cells" = "Il7r",
  "CD4 cells" = "Trac",
  "CD4 cells" = "Cd3e",
  "CD4 cells" = "Lck",
  
  # CD8 cells
  "CD8 cells" = "Cd8a",
  "CD8 cells" = "Cd8b1",
  "CD8 cells" = "Gzma",
  "CD8 cells" = "Gzmb",
  #"CD8 cells" = "Nkg7",
  "CD8 cells" = "Ccl5",
  
  # Plasma cells
  "Plasma cells" = "Jchain",
  "Plasma cells" = "Xbp1",
  "Plasma cells" = "Sdc1",
  "Plasma cells" = "Tnfrsf17",
  "Plasma cells" = "Pim2",
  
  # Neutrophils
  "Neutrophils" = "S100a8",
  "Neutrophils" = "S100a9",
  "Neutrophils" = "Elane",
  "Neutrophils" = "Mpo",
  "Neutrophils" = "Ctsg",
  "Neutrophils" = "Camp",
  
  # Mast cells
  "Mast cells" = "Cpa3",
  "Mast cells" = "Kit",
  "Mast cells" = "Tpsab1",
  "Mast cells" = "Mcpt4",
  "Mast cells" = "Fcgr3",
  
  # Monocytes/MoDMs
  "Monocytes/MoDMs" = "Ly6c2",
  "Monocytes/MoDMs" = "Ccr2",
  "Monocytes/MoDMs" = "Fcgr1",
  "Monocytes/MoDMs" = "Itgam",
  "Monocytes/MoDMs" = "Cd14",
  
  # Spp1 macrophages
  "Spp1 macrophages" = "Spp1",
  "Spp1 macrophages" = "Apoe",
  "Spp1 macrophages" = "Trem2",
  "Spp1 macrophages" = "Lpl",
  "Spp1 macrophages" = "Ctsd",
  "Spp1 macrophages" = "Gpnmb",
  
  # Alveolar macrophages
  "Alveolar macrophages" = "Siglecf",
  #"Alveolar macrophages" = "Fcgr1",   # 替代 "Cd64"
  "Alveolar macrophages" = "Mafb",
  "Alveolar macrophages" = "Pparg",
  "Alveolar macrophages" = "Marco",
  
  # NK/NKT cells
  "NK/NKT cells" = "Nkg7",
  "NK/NKT cells" = "Klrb1c",
  #"NK/NKT cells" = "Gzmb",
  "NK/NKT cells" = "Cd247",
  "NK/NKT cells" = "Klrk1",
  #"NK/NKT cells" = "Trac",
  
  # Activated macrophages
  "Activated macrophages" = "Nos2",
  #"Activated macrophages" = "Il1b",
  "Activated macrophages" = "Tnf",
  "Activated macrophages" = "Cxcl9",
  "Activated macrophages" = "Cxcl10",
  
  # Inflammatory macrophages
  "Inflammatory macrophages" = "Il1b",
  "Inflammatory macrophages" = "Ccl2",
  "Inflammatory macrophages" = "Ccl3",
  "Inflammatory macrophages" = "Ccl4",
  # 注意：S100a8/S100a9 已归入 Neutrophils，此处可省略以减少跨类重复
  
  # cDCs
  "cDCs" = "Itgax",
  #"cDCs" = "H2-Ab1",
  "cDCs" = "Flt3",
  "cDCs" = "Zbtb46",
  
  # pDCs
  "pDCs" = "Siglech",
  "pDCs" = "Bst2",
  "pDCs" = "Tcf4",
  "pDCs" = "Irf7",
  "pDCs" = "Irf8",
  # 移除 "Lilra4"（小鼠中该基因复杂，且非核心 marker）
  
  # Migratory DCs
  "Migratory DCs" = "Ccr7",
  "Migratory DCs" = "Lamp3",
  "Migratory DCs" = "Ccl22",
  #"Migratory DCs" = "Cd207",
  "Migratory DCs" = "Xcr1",
  
  # Proliferating T cells
  "Proliferating T cells" = "Mki67",
  "Proliferating T cells" = "Top2a",
  "Proliferating T cells" = "Cdk1",
  "Proliferating T cells" = "Pcna",
  
  # B cells
  "B cells" = "Ms4a1",
  "B cells" = "Cd19",
  "B cells" = "Cd79a",
  "B cells" = "Bank1",
  "B cells" = "H2-Ab1",
  # 移除 "Ighm"（可变区，不稳定）
  
  # Basophils
  "Basophils" = "Mcpt8",
  "Basophils" = "Enpp3",      # 替代 "Cd203c"（小鼠官方基因名）
  "Basophils" = "Fcer1a",
  "Basophils" = "Il4",
  "Basophils" = "Il13",
  
  # Langerhans cells
  "Langerhans cells" = "Cd207",
  "Langerhans cells" = "Epcam"
  #"Langerhans cells" = "Itgax",
  #"Langerhans cells" = "Ccr7"
  #"Langerhans cells" = "H2-Ab1",
  
  # Proliferating NK/NKT
  # "Proliferating NK/NKT" = "Mki67",
  # "Proliferating NK/NKT" = "Top2a"
  # 其他如 Nkg7/Klrb1c/Trac 已在主类中定义，此处可省略以聚焦“增殖”特征
)
valid_immune_major_marker_ls <- immune_major_marker_ls[immune_major_marker_ls %in% rownames(MESENCHYMAL)]
cat("有效基因:", length(valid_immune_major_marker_ls), "/", length(immune_major_marker_ls))

# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_gene <- setdiff(immune_major_marker_ls, valid_immune_major_marker_ls)
missing_gene



# find markers
Idents(MESENCHYMAL) <- "Xenium_snn_res.0.3"
MESENCHYMAL_markers <-
  FindAllMarkers(
    JoinLayers(MESENCHYMAL),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

MESENCHYMAL_top10_markers <- 
  MESENCHYMAL_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10) %>%
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

names(MESENCHYMAL_top10_markers) <- sub('.*\\.','',names(MESENCHYMAL_top10_markers))


#提取deepseek
MESENCHYMAL_top10_markers1 <-
  MESENCHYMAL_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10)


# dotplot
pdf(file.path(MESENCHYMAL_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 35, height = 6)
DotPlot(MESENCHYMAL, group.by = "Xenium_snn_res.0.3", features = major_marker)+RotatedAxis()
DotPlot(MESENCHYMAL, group.by = "Xenium_snn_res.0.3", features = MESENCHYMAL_major_marker_ls)+RotatedAxis()
DotPlot(MESENCHYMAL, group.by = "Xenium_snn_res.0.3", features = MESENCHYMAL_major_marker_ls1)+RotatedAxis()
DotPlot(MESENCHYMAL, group.by = "Xenium_snn_res.0.3", features = immune_major_marker_ls)+RotatedAxis()
DotPlot(MESENCHYMAL, group.by = "Xenium_snn_res.0.3", features = MESENCHYMAL_top10_markers)+RotatedAxis()
dev.off()


for (gene in valid_MESENCHYMALgenes) {
  p <- FeaturePlot(MESENCHYMAL, features = gene, reduction = 'harmonyumap')+RotatedAxis()
  ggsave(
    filename = file.path(MESENCHYMAL_out_dir, paste0("featureplot/",gene, "_feature.png")),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}


## annotation -----------
cell_cluster <- c(
  '0' = "Activated fibroblasts FBs", # Dcn Lum Col1a1 Col1a2 Col3a1 Fn1 Fap
  '1' = "SMCs or fibrocytes FBs", # Acta2 Cspg4
  '2' = "Capillary", # PLVAP CLDN5 APLN
  '3' = "Alveolar FBs", # Wnt5a Fgf10 Fgf7 Pi16
  '4' = "Alveolar macrophage", # CD68 CD163 C1QA
  '5' = "Myofibroblasts FBs", # Acta2 Fap Dcn Lum Col1a1 Col1a2 Fn1
  '6' = "SMCs or fibrocytes FBs", # Acta2 Cspg4
  '7' = "delete", # no marker
  '8' = "Subpleural mesothelial", # Msln Has1 Krt8 Krt18
  '9' = "Alveolar FBs", # Wnt5a Fgf10 Fgf7 Pi16
  '10' = "delete", # no marker
  '11' = "delete" # no marker
)

Idents(MESENCHYMAL) <- "Xenium_snn_res.0.3"

MESENCHYMAL <- RenameIdents(MESENCHYMAL, cell_cluster)
MESENCHYMAL$subtype <- Idents(MESENCHYMAL)

pdf(file.path(MESENCHYMAL_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(MESENCHYMAL, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(MESENCHYMAL_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 7, width = 18)
DimPlot(MESENCHYMAL, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()





## cell props barplot  --------
MESENCHYMAL_meta <- MESENCHYMAL[[]]

MESENCHYMAL_meta$sample <- factor(MESENCHYMAL_meta$sample, 
                                  levels=c(MESENCHYMAL_meta$sample %>% unique())
)
p_col_proj <-
  MESENCHYMAL_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path("result/MESENCHYMAL_subtype/cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  MESENCHYMAL_meta %>%
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
  filename = file.path("result/MESENCHYMAL_subtype/cell_freq_group.pdf"),
  height = 6,
  width  = 4
)


## volcano-------
MESENCHYMAL_markers <-
  FindAllMarkers(
    JoinLayers(MESENCHYMAL),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(MESENCHYMAL_markers,'result/MESENCHYMAL_subtype/MESENCHYMAL_diff_gene.csv')

library(EnhancedVolcano)
type<-unique(MESENCHYMAL_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- MESENCHYMAL_markers[MESENCHYMAL_markers$cluster == type, ]
    p<-EnhancedVolcano(markers,
                       lab = rownames(markers),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = str_glue('{name} volcano plot'),
                       pCutoff = 10e-2,
                       FCcutoff = 0,
                       pointSize = 3.0,
                       labSize = 6.0)
    p
  })
volcano_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path('result/MESENCHYMAL_subtype', str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })

MESENCHYMAL_markers_ls <- list(MESENCHYMAL_markers)
names(MESENCHYMAL_markers_ls) <- 'MESENCHYMAL'

# scatter
p_volcanoplot_major_ls <- 
  MESENCHYMAL_markers_ls %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 0.5
    dodge_width <- 0.1
    
    df_scatter <-
      dat %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      slice_sample(prop = 1) %>%
      group_by(cluster) %>%
      
      mutate(
        x_delta = 1:n(),
        x_delta = scales::rescale(x_delta, to = c(dodge_width, 1 - dodge_width)),
        xx = x_base - 1 + x_delta,
        y = ifelse(avg_log2FC > 0, avg_log2FC + rect_height, avg_log2FC - rect_height)
      ) %>%
      mutate(
        group = ifelse(p_val_adj > 0.01, "p > 0.01", "p < 0.01"),
        group = factor(group, levels = c("p > 0.01", "p < 0.01"))
      ) %>%
      dplyr::arrange(group)
    
    # df_label 
    df_label_up <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_max(avg_log2FC, n = 10)
    
    df_label_down <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_min(avg_log2FC, n = 10)
    
    
    
    # grp rect  
    df_grp_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise() %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1, xmax = x_base) %>%
      mutate(ymin = -rect_height, ymax = rect_height)
    
    
    # background rect
    df_bg_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise(ymax = max(avg_log2FC) + rect_height, ymin = min(avg_log2FC) - rect_height) %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1 + dodge_width, xmax = x_base - dodge_width)
    
    
    
    df_scatter %>%
      ggplot(aes(x = xx, y = y)) +
      geom_rect(
        data = df_bg_rect, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        inherit.aes = FALSE
      ) +
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_point(aes(color = group), size = 1) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text.x = element_blank()
      ) +
      labs(
        x = NULL,
        y = "Log2FC",
        title = name
      )
    
  })

p_volcanoplot_major_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path('result/MESENCHYMAL_subtype', str_glue("{name} diffgene dot.pdf")),
      height = 10,
      width = 20
    )
  })


## enrichment ----------------------------
MESENCHYMAL_de_genes_for_enrichmenet <- 
  MESENCHYMAL_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.15) %>%
  dplyr::filter(abs(avg_log2FC) > 0.5) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")


enrich_ls <- 
  MESENCHYMAL_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichGO(
      x,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
      ont='ALL',
      pAdjustMethod = 'BH',
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.05,
      keyType = 'SYMBOL'
    )
  })


er_output<-'result/MESENCHYMAL_subtype/enrichment' 
dir.create(er_output)

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(er_output, name))
    
    write.csv(
      x@result,
      file = file.path(er_output, name, "Go enrichment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
# export dotplot
enrich_dot_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "go dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })

MESENCHYMAL_de_genes_id_for_enrichmenet <-
  MESENCHYMAL_de_genes_for_enrichmenet %>%
  map(function(x){ # browser()
    clusterProfiler::bitr(
      x,
      OrgDb = "org.Hs.eg.db",
      fromType = "SYMBOL",
      toType = "ENTREZID"
    ) %>%
      .$ENTREZID
  }) %>%
  purrr::discard(function(x) length(x) < 1)


enrich_ls <-
  MESENCHYMAL_de_genes_id_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichKEGG(
      x,
      organism = "hsa", 
      keyType = "ncbi-geneid", 
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
    )
  })

# saveRDS(enrich_ls, file.path(er_output, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(er_output, name))
    
    write.csv(
      x@result,
      file = file.path(er_output, name, "KEGG enriment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })
# export dotplot
enrich_dot_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "kegg dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })

saveRDS(MESENCHYMAL, 'result/MESENCHYMAL_subtype/MESENCHYMAL.rds')

