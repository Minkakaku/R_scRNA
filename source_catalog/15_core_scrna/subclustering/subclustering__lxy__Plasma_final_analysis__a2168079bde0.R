all_data <- readRDS("output/data/all_data.rds")


Plasma_out_dir <- "Plasma_subtype" 
dir.create(Plasma_out_dir)


library(Seurat)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Hs.eg.db)

#all_data <- readRDS(file.path("R_OUTPUT", "all_data.rds"))

# markers: https://cloud.tencent.com/developer/article/2115987
major_marker_ls <- c(
  "Adipocyte" = "PLIN1",
  "Adipocyte" = "PLIN4",
  "Adipocyte" = "ADIPOQ",
  "Adipocyte" = "FABP4",
  "Adipocyte" = "LPL",
  
  "cTEC" = "LY75", 
  "cTEC" = "CCL25", 
  "cTEC" = "CXCL12", 
  
  "mTEC" = "EPCAM", 
  "mTEC" = "CCL19", 
  "mTEC" = "CCL21", 
  
  "Endothelial" = "PECAM1",
  "Endothelial" = "VWF", 
  "Endothelial" = "AQP1", 
  "Endothelial" = "CDH5", 
  "Endothelial" = "PLVAP", 
  
  "Fibroblast" = "PDGFRA", 
  "Fibroblast" = "CCDC80", 
  "Fibroblast" = "COL1A1", 
  "Fibroblast" = "COL1A2", 
  "Fibroblast" = "DCN", 
  "Fibroblast" = "ACTA2", 
  
  "B cell" = "CD19", 
  "B cell" = "CD79A", 
  "B cell" = "MS4A1", 
  "Plasma" = "TNFRSF17", 
  "Plasma" = "JCHAIN", 
  "Plasma" = "IGJ", 
  "Plasma" = "MZB1", 
  "Plasma" = "XBP1", 
  
  "T cell" = "CD3D", 
  "T cell" = "CD3E", 
  "T cell" = "CD3G", 
  "T cell" = "TRAC", 
  "T cell" = "BCL11B",
  "T cell" = "CD4", 
  "T cell" = "FOXP3", 
  "T cell" = "CD8A", 
  "T cell" = "GZMB",
  
  "NK" = "GNLY", 
  "NK" = "NKG7", 
  "NK" = "NCAM1", 
  "NK" = "UNC13D", 
  
  "Proliferating T cell" = "MKI67", 
  "Proliferating T cell" = "TOP2A", 
  
  "Pericyte" = "PDGFRB", 
  "Pericyte" = "RGS5", 
  "Pericyte" = "MCAM", 
  #"Pericyte" = "PLVAP", 
  
  "Monocyte" = "CD14", 
  "Monocyte" = "FCGR3A", 
  "Monocyte" = "ITGAM", 
  
  "Macrophage" = "CD68", 
  "Macrophage" = "CD163", 
  "Macrophage" = "C1QA", 
  "Macrophage" = "C1QC", 
  "Macrophage" = "CSF1R", 
  "Macrophage" = "IL1B", 
  
  "Dendritic cell" = "CD83", 
  "Dendritic cell" = "LILRA4", 
  #"Dendritic cell" = "CLEC9A", 
  "Dendritic cell" = "IDO1", 
  "Dendritic cell" = "CD1C", 
  "Dendritic cell" = "FCER1A", 
  "Dendritic cell" = "CST3", 
  
  "pDC" = "IRF4", 
  "pDC" = "IRF7", 
  "pDC" = "IRF8", 
  
  "cDC1" = "XCR1",
  "cDC1" = "CLEC9A",
  
  "Lamp3_DC" = "LAMP3", 
  
  "Netrophil" = "CSF3R", 
  "Netrophil" = "FCGR3B", 
  "Netrophil" = "S100A8", 
  "Netrophil" = "S100A9", 
  "Netrophil" = "G0S2", 
  "Netrophil" = "RETNLG", 
  
  "Mast cell" = "MS4A2", 
  "Mast cell" = "KIT", 
  "Mast cell" = "HDC", 
  "Mast cell" = "TPSAB1", 
  "Mast cell" = "CTSG", 
  "Mast cell" = "CPA3", 
  
  "Smooth muscle cell" = "PLN", 
  "Smooth muscle cell" = "MYLK", 
  "Smooth muscle cell" = "MYH11", 
  "Smooth muscle cell" = "ACTG2", 
  "Smooth muscle cell" = "TAGLN", 
  
  "Schwann cell" = "MPZ", 
  "Schwann cell" = "SOX10", 
  
  "Neuron" = "TUBB3", 
  "Neuron" = "SNAP25", 
  
  "Mesothelial" = "MSLN", 
  
  #"Mesenchyme" = "PDGFRA", 
  
  "Erythrocytes" = "HBA-A1", 
  "Erythrocytes" = "HBA-A2"
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )



# Plasmas -------------------------------------------------------------------
#这次只提取一种细胞，不用替换
#all_data$annotation <- str_replace(all_data$annotation, 'Proliferating T cell', 'T cell')

Plasma <- subset(all_data, T_annotation == 'Plasma')

## PCA --------
Plasma <- FindVariableFeatures(Plasma)

options(future.globals.maxSize = 50 * 1024^3) # 把并行全局变量上限改成 50GB（足够用）
Plasma <- ScaleData(Plasma, vars.to.regress = 'nCount_RNA')
Plasma <- RunPCA(Plasma)

### pca inspector ------
pdf(file.path(Plasma_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Plasma, ndims = 50)
dev.off()


## dimension reduction --------
dims <- 12 # 18

library(future)
library(Seurat)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Hs.eg.db)
#options(future.globals.maxSize = 50 * 1024^3)
#options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
#plan(multisession, workers = 10)

Plasma <- RunUMAP(Plasma, dims = 1:dims,return.model = TRUE)

pdf(file.path(Plasma_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Plasma, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

Plasma <- Plasma %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Plasma <- RunUMAP(
  Plasma, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(Plasma_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Plasma, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'

#
### cluster -------
Plasma  <-
  FindNeighbors(
    Plasma ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Plasma <-
  FindClusters(
    Plasma,
    resolution = c(0.1,0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(Plasma_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Plasma, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Plasma, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Plasma, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Plasma, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Plasma, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Plasma, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Plasma, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Plasma, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

Idents(Plasma) <- "RNA_snn_res.0.6"

### dotplot -------
Plasma_major_marker_ls <- c(
  "PB" = "MKI67", #增殖质母细胞，表达MS4A1低
  "PB" = "UBE2C",
  "PB" = "BIRC5",
  "PB" = "SPN",
  "PB" = "TUBB",
  
  "early 1" = "CD52",
  "early 1" = "PEBP1",
  "early 1" = "HLA-DOB",
  
  "early 2" = "CD79A",
  "early 2" = "CD74",
  
  "trans" = "IGLV3-25",
  
  "late 1" = "SMOC1",
  "late 1" = "MOXD1",
  
  "late 2" = "CD9",
  "late 2" = "CST3",
  "late 2" = "CD63",
  "late 2" = "TIMP1",
  
  "late 3" = "RHOB",
  "late 3" = "CDKN1A",
  "late 3" = "RGCC",
  
  "late 4" = "JUNB",
  "late 4" = "FOS",
  "late 4" = "EGR1",
  "late 4" = "NR4A1",
  
  "pc9" = "ATF5",
  "pc9" = "PSAT1",
  "pc9" = "CEBPB",
  "pc9" = "ERN1",
  
  "IFN+" = "MX1",
  "IFN+" = "XAF1",
  "IFN+" = "IRF7",
  "IFN+" = "STAT1",
  
  "Mito-high" = "MALAT1",
  "Mito-high" = "IRF4",
  "Mito-high" = "PRDM1",
  "Mito-high" = "ZBTB20",
  "Mito-high" = "MTRNR2L12",
  "Mito-high" = "FCRL5",
  
  "pc12" = "JUND",
  "pc12" = "CXCR4",
  "pc12" = "SIK1B",
  
  "IgM1" = "CCDC88A",
  "IgM1" = "FOXP1",
  "IgM1" = "IGHM",
  "IgM1" = "KLHL14",
  "IgM1" = "MS4A1",
  
  "IgM2" = "JCHAIN",
  "IgM2" = "RGS2",
  "IgM2" = "TNFRSF4",
  
  "IgM3" = "CCL3",
  "IgM3" = "CCL4",
  
  "LLPC" = "CD19",
  "LLPC" = "CD38",
  "LLPC" = "CD138"
)



Idents(Plasma)<-'RNA_snn_res.0.6'
Plasma_markers <-
  FindAllMarkers(
    JoinLayers(Plasma),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(Plasma_markers, file = "Plasma_subtype/Plasma_all_markers.csv", row.names = TRUE)

Plasma_top5_markers <- 
  Plasma_markers %>% 
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

names(Plasma_top5_markers) <- sub('.*\\.','',names(Plasma_top5_markers))


# dotplot
# 如果已做其他细胞类型的亚群分析，这里就可以顺便画一下其他类型的细胞亚群基因，比如Fibro_major_marker_ls2，否则只用画T细胞
pdf(file.path(Plasma_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(Plasma, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Plasma, group.by = "RNA_snn_res.0.6", features = Plasma_major_marker_ls)+RotatedAxis()+theme_90angle
#DotPlot(Plasma, group.by = "RNA_snn_res.0.6", features = Plasma_markers)+RotatedAxis()
DotPlot(Plasma, group.by = "RNA_snn_res.0.6", features = Plasma_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(Plasma, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()

pdf(file.path(Plasma_out_dir, str_glue("dotplot dims{dims} T_major_marker.pdf")), width = 25, height = 6)
DotPlot(Plasma, group.by = "RNA_snn_res.0.6", features = T_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(Plasma, group.by = "RNA_snn_res.0.6", features = T_major_marker_ls_1)+RotatedAxis()+theme_90angle
dev.off()

saveRDS(Plasma ,file = "Plasma_subtype/Plasma_cluster.rds")
Plasma <- readRDS("Plasma_subtype/Plasma_cluster.rds")

some_genes <- c("CD19", "CD38", "CD138")

FeaturePlot(
  Plasma,
  features = some_genes,
  ncol = 3,          # 3个图并排
  pt.size = 0.2,     # 点大小
  label = TRUE,
  repel = TRUE
)

ggsave(
  "Plasma_UMAP_CD19_CD38_CD138.pdf",
  width = 15, height = 4
)


VlnPlot(
  Plasma,
  features = some_genes,
  ncol = 3,
  pt.size = 0      
)


ggsave(
  "Plasma_VlnPlot_CD19_CD38_CD138.pdf",
  width = 15, height = 6
)

## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "Naive_T", # 
  '1'  = "Mito-high", #
  '2'  = "Mito-high", #
  '3'  = "TEC", # 
  '4'  = "B cell", # 
  '5'  = "Fibroblast", # 
  '6'  = "Macrophage", # 
  '7'  = "Proliferative plasmablast", # 
  '8'  = "Adipocyte",
  '9'  = "B cell"
)

Idents(Plasma) <- "RNA_snn_res.0.6"

Plasma <- RenameIdents(Plasma, cell_cluster)
Plasma$subtype <- Idents(Plasma)


pdf(file.path(Plasma_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Plasma, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Plasma_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 15)
DimPlot(Plasma, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


# 更改原注释
# 1. 先把 all_data$T_annotation 复制到新列 Plasma_T_annotation
all_data$Plasma_T_annotation <- as.character(all_data$T_annotation)

# 2. 创建映射：细胞名 → 精细T亚型
cell_subtype <- Plasma$subtype
names(cell_subtype) <- colnames(Plasma)

# 3. 只把 "T cell" 替换成精细亚型，其他细胞保持原样
all_data$Plasma_T_annotation[colnames(Plasma)] <- as.character(cell_subtype)

# 查看结果（正确）
table(all_data$Plasma_T_annotation)
table(all_data$T_annotation)
table(Plasma$subtype)
saveRDS(all_data, "Plasma_subtype/all_data_Plasma_T_annotation.rds")


#只保留T细胞亚群
keep <- c(
  "Mito-high",
  "Proliferative plasmablast"
)

Plasma_clean <- subset(Plasma, subtype %in% keep)

Plasma_clean$subtype <- droplevels(Plasma_clean$subtype)
Idents(Plasma_clean) <- Plasma_clean$subtype

table(Plasma_clean$subtype)

Plasma <- Plasma_clean
#重新分布一下散点
dims <- 12
Plasma <- RunUMAP(Plasma, dims = 1:dims,return.model = TRUE)
pdf(file.path(Plasma_out_dir, str_glue("Dimplot dim{dims} Plasma annotation.pdf")), height = 6, width = 8)
DimPlot(Plasma, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Plasma_out_dir, str_glue("Dimplot dim{dims} Plasma annotation split.pdf")), height = 6, width = 15)
DimPlot(Plasma, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

saveRDS(Plasma_clean, "Plasma_subtype/Plasma_clean.rds")

#画一个注释完的点图
Plasma_major_marker_ls_1 <- c(
  "PB" = "MKI67", #增殖质母细胞，表达MS4A1低
  "PB" = "UBE2C",
  "PB" = "BIRC5",
  "PB" = "SPN",
  "PB" = "TUBB",
  
  "Mito-high" = "MALAT1",
  "Mito-high" = "IRF4",
  "Mito-high" = "PRDM1",
  "Mito-high" = "ZBTB20",
  "Mito-high" = "MTRNR2L12",
  "Mito-high" = "FCRL5"
)
pdf(file.path(Plasma_out_dir, str_glue("dotplot Plasma subtype major_marker.pdf")), width = 10, height = 6)
DotPlot(Plasma, group.by = "subtype", features = Plasma_major_marker_ls_1 )+RotatedAxis()+theme_90angle
dev.off()

## cell props barplot  ----------------------------------------------------------
Plasma_meta <- Plasma[[]]

Plasma_meta$sample <- factor(Plasma_meta$sample,levels=c(Plasma_meta$sample %>% unique()))

p_col_proj <-
  Plasma_meta %>%
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
  filename = file.path("Plasma_subtype", "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  Plasma_meta %>%
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
  filename = file.path(Plasma_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)


# 映射到空间切片
subtype_col <- 'subtype'

cells_position_ls <-
  cbind(
    Plasma[['sample']],
    Embeddings(Plasma[['spatial']])
  ) %>%
  rownames_to_column('cell_id') %>%
  group_by(sample) %>%
  nest() %>%
  deframe()

some_genes <- c("CD19", "CD38", "CD138")

plot_dat <-
  FetchData(
    Plasma, 
    vars = c(subtype_col, some_genes)
  ) %>%
  rownames_to_column('cell_id')


# plot celltype ------------------
all_cells <-  Plasma$subtype %>% unique %>% as.character()
cell_colors <- structure(
  names = all_cells,
  ggsci:: pal_igv()(length(all_cells))
)



p_sp_ct_ls <-
  cells_position_ls %>% # .[1] %>%
  imap(function(coords, sam){
    dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
    
    height <- 5
    width  <- 5
    p <-
      dat2 %>%
      ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
      geom_tile(
        aes(fill = .data[[subtype_col]]),
        height = height, 
        width = width,
        show.legend = c(fill = TRUE, alpha = FALSE)
      ) +
      theme_void() +
      coord_equal() +
      scale_y_reverse() +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values = cell_colors) +
      theme(panel.background = element_rect(fill = "black", color = NA)) +
      ggtitle(sam)
    p
    
  })


pdf(
  file.path("Plasma_subtype/spatial plot of celltype in each sample.pdf"),
  height = 10,
  width  = 10
)
p_sp_ct_ls %>% walk(print)
dev.off()


# plot each cell in all sample -----------
p_sp_ct_ls_ls <-
  cells_position_ls %>% # .[1] %>%
  imap(function(coords, sam){
    dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
    
    height <- 5
    width  <- 5
    
    p_ls <-
      all_cells %>%
      setNames(., .) %>%
      map(function(x){
        cell_colors[names(cell_colors) != x] <- NA
        
        p <-
          dat2 %>%
          ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
          geom_tile(
            aes(fill = .data[[subtype_col]]),
            height = height, 
            width = width,
            show.legend = c(fill = TRUE, alpha = FALSE)
          ) +
          theme_void() +
          coord_equal() +
          scale_y_reverse() +
          # ggsci::scale_fill_igv() +
          scale_fill_manual(values = cell_colors) +
          theme(panel.background = element_rect(fill = "black", color = NA)) +
          ggtitle(sam)
        p
      })
    p_ls
  })

# export by sample
p_sp_ct_ls_ls %>%
  iwalk(function(p_ls, sample){
    pdf(
      file.path(str_glue("Plasma_subtype/spatial plot of celltype in each sample - {sample}.pdf")),
      height = 10,
      width  = 10
    )
    p_ls %>% walk(print)
    dev.off()
  })






























