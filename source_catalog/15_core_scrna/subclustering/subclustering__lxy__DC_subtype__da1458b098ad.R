library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Hs.eg.db)


all_data <- readRDS('result_0622/T_cell_subtype_2/all_data_T_annotation.rds')

all_data$annotation1 %>% table()


DC_out_dir <- "result_0622/DC_subtype" 
dir.create(DC_out_dir)

all_data$annotation1 <- str_replace(all_data$annotation1, 'pDC', 'DC')
all_data$annotation1 <- str_replace(all_data$annotation1, 'cDC1', 'DC')
all_data$annotation1 <- str_replace(all_data$annotation1, 'LAMP3_DC', 'DC')

DC <- subset(all_data, annotation1 == 'DC')

## PCA --------
DC <- FindVariableFeatures(DC)

DC <- ScaleData(DC, vars.to.regress = 'nCount_RNA')
DC <- RunPCA(DC)

### pca inspector ------
pdf(file.path(DC_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(DC, ndims = 50)
dev.off()


## dimension reduction --------
dims <- 22 # 18

#library(future)
#options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
#plan(multisession, workers = 10)

DC <- RunUMAP(DC, dims = 1:dims,return.model = TRUE)

pdf(file.path(DC_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(DC, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

DC <- DC %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

DC <- RunUMAP(
  DC, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(DC_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(DC, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
DC  <-
  FindNeighbors(
    DC ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


DC <-
  FindClusters(
    DC,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(DC_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(DC, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(DC, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(DC, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(DC, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(DC, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(DC, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(DC, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

Idents(DC) <- "RNA_snn_res.0.8"

### dotplot -------
DC_major_marker_ls <- c(
  "DC" = "CD83",
  
  "pDC" = "LILRA4",
  "pDC" = "IRF4", 
  "pDC" = "IRF7", 
  "pDC" = "IRF8", 
  
  "cDC1" = "XCR1",
  "cDC1" = "CLEC9A",
  
  "cDC2" = "CD1C", 
  "cDC2" = "FCER1A", 
  "cDC2" = "CLEC10A",
  
  "mregDC" = "LAMP3",
  "mregDC" = "CCR7"
)

Idents(DC)<-'RNA_snn_res.0.8'
DC_markers <-
  FindAllMarkers(
    JoinLayers(DC),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(DC_markers, file = "result_0622/DC_subtype/DC_all_markers.csv", row.names = TRUE)

DC_top5_markers <- 
  DC_markers %>% 
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

names(DC_top5_markers) <- sub('.*\\.','',names(DC_top5_markers))


# dotplot
# 如果已做其他细胞类型的亚群分析，这里就可以顺便画一下其他类型的细胞亚群基因，比如Fibro_major_marker_ls2，否则只用画T细胞
pdf(file.path(DC_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(DC, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(DC, group.by = "RNA_snn_res.0.4", features = DC_major_marker_ls)+RotatedAxis()+theme_90angle
#DotPlot(DC, group.by = "RNA_snn_res.0.4", features = PT_major_marker_ls_1)+RotatedAxis()+theme_90angle
DotPlot(DC, group.by = "RNA_snn_res.0.4", features = DC_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(DC, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()


pdf(file.path(DC_out_dir, str_glue("dotplot dims{dims} major_marker_0.8.pdf")), width = 25, height = 6)
DotPlot(DC, group.by = "RNA_snn_res.0.8", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(DC, group.by = "RNA_snn_res.0.8", features = DC_major_marker_ls)+RotatedAxis()+theme_90angle
#DotPlot(DC, group.by = "RNA_snn_res.0.4", features = PT_major_marker_ls_1)+RotatedAxis()+theme_90angle
DotPlot(DC, group.by = "RNA_snn_res.0.8", features = DC_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(DC, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()




saveRDS(DC ,file = "result_0622/DC_subtype/DC_cluster.rds")
DC <- readRDS("result_0622/DC_subtype/DC_cluster.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "T_cell-like_DC", # 
  '1'  = "cDC2", #
  '2'  = "mregDC", # 
  '3'  = "cDC1", # 
  '4'  = "pDC", # 
  '5'  = "T_cell-like_DC",#  
  '6'  = "pDC", #  
  '7'  = "mix_DC", #  
  '8'  = "T_cell-like_DC", #  
  '9'  = "T_cell-like_DC", #  
  '10'  = "T_cell-like_DC", #  
  '11'  = "B cell", #  
  '12'  = "mregDC", #  
  '13'  = "mregDC", #  
  '14'  = "cDC1", #  
  '15'  = "Macrophage-like_DC", #  
  '16'  = "Macrophage-like_DC", #  
  '17'  = "Neutrophil-like_DC", #  
  '18'  = "Mast cell"
)

Idents(DC) <- "RNA_snn_res.0.8"
DC <- RenameIdents(DC, cell_cluster)
DC$subtype <- Idents(DC)


pdf(file.path(DC_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(DC, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(DC_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 15)
DimPlot(DC, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


keep <- c(
  "cDC1",
  "cDC2",
  "pDC",
  "mregDC"
)

DC_clean <- subset(DC, subtype %in% keep)

DC_clean$subtype <- droplevels(DC_clean$subtype)
Idents(DC_clean) <- DC_clean$subtype

table(DC_clean$subtype)
saveRDS(DC_clean, "result_0622/DC_subtype/DC_clean.rds")

DC <- DC_clean

DC <- RunUMAP(DC, dims = 1:22,n.neighbors = 50,return.model = TRUE)

pdf(file.path(DC_out_dir, str_glue("Dimplot dim22 DC annotation.pdf")), height = 6, width = 8)
DimPlot(DC, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(DC_out_dir, str_glue("Dimplot dim22 DC annotation split.pdf")), height = 6, width = 15)
DimPlot(DC, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()








DC_1 <- readRDS("result_0622/DC_subtype/DC_cluster.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "T_cell-like_DC", # 
  '1'  = "cDC2", #
  '2'  = "mregDC", # 
  '3'  = "cDC1", # 
  '4'  = "pDC", # 
  '5'  = "T_cell-like_DC",#  
  '6'  = "pDC", #  
  '7'  = "mix_DC", #  
  '8'  = "T_cell-like_DC", #  
  '9'  = "T_cell-like_DC", #  
  '10'  = "T_cell-like_DC", #  
  '11'  = "B cell", #  
  '12'  = "mregDC", #  
  '13'  = "mregDC", #  
  '14'  = "cDC1", #  
  '15'  = "Macrophage-like_DC", #  
  '16'  = "Macrophage-like_DC", #  
  '17'  = "Neutrophil-like_DC", #  
  '18'  = "Mast cell"
)

Idents(DC_1) <- "RNA_snn_res.0.8"
DC_1 <- RenameIdents(DC_1, cell_cluster)
DC_1$subtype <- Idents(DC_1)
saveRDS(DC_1, "result_0622/DC_subtype/DC_allannotation.rds")
######下面合并到alldata的注
all_data$annotation1 %>% table()
#all_data <- readRDS("output/data/all_data.rds")
# 1. 先把 all_data$annotation 复制到新列 T_annotation
all_data$annotation2 <- as.character(all_data$annotation1)

# 2. 创建映射：细胞名 → 精细T亚型
cell_subtype <- DC_1$subtype
names(cell_subtype) <- colnames(DC_1)

# 3. 只把 "T cell" 替换成精细亚型，其他细胞保持原样
all_data$annotation2[colnames(DC_1)] <- as.character(cell_subtype)

# 查看结果（正确）
table(all_data$annotation2)
table(all_data$annotation1)
table(DC_1$subtype)
saveRDS(all_data, "result_0622/DC_subtype/all_data_T_DC_annotation.rds")



DC <- readRDS("result_0622/DC_subtype/DC_clean.rds")
#画一个注释完的点图
DC_major_marker_ls_2 <- c(
  "DC" = "CD83",
  
  "pDC" = "LILRA4",
  #"pDC" = "IRF4", 
  "pDC" = "IRF7", 
  #"pDC" = "IRF8", 
  
  "cDC1" = "XCR1",
  "cDC1" = "CLEC9A",
  
  "cDC2" = "CD1C", 
  "cDC2" = "FCER1A", 
  "cDC2" = "CLEC10A",
  
  "mregDC" = "LAMP3",
  "mregDC" = "CCR7"
)

pdf(file.path(DC_out_dir, str_glue("dotplot DC subtype major_marker.pdf")), width = 12, height = 6)
DotPlot(DC, group.by = "subtype", features = DC_major_marker_ls_2)+RotatedAxis()+theme_90angle
dev.off()


DC$subtype %>% table()
DC$sample %>% table()
DC$group %>% table()
## cell props barplot  ----------------------------------------------------------
DC_meta <- DC[[]]

DC_meta$sample <- factor(DC_meta$sample,levels=c(DC_meta$sample %>% unique()))

p_col_proj <-
  DC_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  scale_fill_brewer(palette = "Set2", name = "DC Subtype") +
  #ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(DC_out_dir, "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  DC_meta %>%
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
  filename = file.path(DC_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)



# 映射到空间切片
subtype_col <- 'subtype'

cells_position_ls <-
  cbind(
    DC[['sample']],
    Embeddings(DC[['spatial']])
  ) %>%
  rownames_to_column('cell_id') %>%
  group_by(sample) %>%
  nest() %>%
  deframe()


plot_dat <-
  FetchData(
    DC, 
    vars = c(subtype_col)
  ) %>%
  rownames_to_column('cell_id')


# plot celltype ------------------
all_cells <-  DC$subtype %>% unique %>% as.character()
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
  file.path("result_0622/DC_subtype/spatial plot of celltype in each sample.pdf"),
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
      file.path(str_glue("result_0622/DC_subtype/spatial plot of celltype in each sample - {sample}.pdf")),
      height = 10,
      width  = 10
    )
    p_ls %>% walk(print)
    dev.off()
  })
