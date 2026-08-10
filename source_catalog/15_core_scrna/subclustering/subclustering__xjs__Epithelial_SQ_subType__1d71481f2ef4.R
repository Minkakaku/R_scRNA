library(Seurat)
library(tidyverse)

major_marker_ls <- c(
  'epithelials'  = 'EPCAM',
  'epithelials'  = 'KRT5',
  'epithelials'  = 'KRT14',
  
  "Foveolar cells" = "GKN1", #ACT
  "Foveolar cells" = "MUC5AC",
  "Foveolar cells" = "TFF1",
  
  'endothelials' = 'PECAM1', 
  'endothelials' = 'VWF', 
  
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "SFRP4",
  "Fibroblasts" = "COL1A1",
  "Fibroblasts" = "COL1A2",
  "Fibroblasts" = "DCN",
  "Fibroblasts" = "ACTA2",
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  "Plasma" = "TNFRSF17",
  "Plasma" = "JCHAIN",
  "Plasma" = "MZB1",
  "Plasma" = "XBP1",
  
  'T cells'           = 'CD3D', 
  'T cells'           = 'CD3E', 
  'T cells'           = 'CD3G', 
  'T cells'           = 'TRAC', 
  
  "NK" = "GNLY",
  "NK" = "NKG7",
  "NK" = "NCAM1",
  
  'Proliferating' = 'MKI67',
  'Proliferating' = 'TOP2A',
  
  "Pericytes" = "PDGFRB",
  "Pericytes" = "RGS5",
  "Pericytes" = "MCAM",
  
  "Monocytes" = "CD14", 
  "Monocytes" = "FCGR3A",
  
  "Macrophage" = "CD68",
  "Macrophage" = "CD163",
  "Macrophage" = "C1QA",
  "Macrophage" = "C1QC", 
  'Macrophage'  = 'CSF1R',
  
  "Dendritic cell" = "CD83",
  'Dendritic cell' = "LILRA4", 
  'pDCs' = "IRF4", 
  'pDCs' = "IRF7", 
  'pDCs' = "IRF8", 
  'Dendritic cell' = "CLEC9A",
  'Dendritic cell' = "IDO1",
  'Dendritic cell' = "CD1C",
  'Dendritic cell' = "FCER1A",
  'Lamp3_DCs' = "LAMP3",
  
  'netrophil' = "CSF3R",
  'netrophil' = "G0S2",
  'netrophil' = "RETNLG",
  
  "Myocytes" = "TNNC1", 
  "Myocytes" = "MYH7", 
  
  "Mast cell" = "CPA3",
  "Mast cell" = "CTSG",
  "Mast cell" = "TPSAB1",
  "Mast cell" = "KIT"
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

all_data <- readRDS(file.path("R_OUTPUT", "all_data.rds"))
all_data$sample %>% table()
all_data$annotation %>% table()
all_data$treat_result <- str_replace_all(all_data$sample, c(
  'CW_SH' = 'non-PCR', 
  'CW_SQ' = 'non-PCR', 
  'HJD_SH' = 'non-PCR', 
  'HJD_SQ' = 'non-PCR', 
  'LGRO_SH' = 'non-PCR', 
  'LGRO_SQ' = 'non-PCR', 
  'WBS_SH' = 'non-PCR', 
  'WBS_SQ' = 'non-PCR', 
  'WLJ_SH' = 'non-PCR', 
  'WLJ_SQ' = 'non-PCR', 
  'JJW_SH' = 'PCR', 
  'JJW_SQ' = 'PCR', 
  'YBK_SH' = 'PCR', 
  'YBK_SQ' = 'PCR'
))
all_data$sample <- factor(all_data$sample, levels=c(
  'CW_SQ', 
  'CW_SH', 
  'HJD_SQ',
  'HJD_SH',
  'LGRO_SQ',
  'LGRO_SH',
  'WBS_SQ',
  'WBS_SH',
  'WLJ_SQ',
  'WLJ_SH',
  'JJW_SQ',
  'JJW_SH',
  'YBK_SQ',
  'YBK_SH'
)
)
saveRDS(all_data, file.path("R_OUTPUT", "all_data.rds"))

Epithelial <- subset(all_data, annotation == 'Epithelial/Tumor cell' & group == 'SQ')

## PCA --------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3) # 80Gb

DefaultAssay(Epithelial) <- "RNA"
Epithelial <- FindVariableFeatures(Epithelial)
Epithelial <- ScaleData(Epithelial, vars.to.regress = 'nCount_RNA')
Epithelial <- RunPCA(Epithelial)

Epithelial_out_dir <-'Epithelial_SQ_subtype' 
dir.create(Epithelial_out_dir)


### pca inspector ------
pdf(file.path(Epithelial_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Epithelial, ndims = 50)
dev.off()


## dimension reduction --------
dims <- 20

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 20 * 1024 ^ 3) # 20Gb


Epithelial <- RunUMAP(Epithelial, dims = 1:dims,return.model = TRUE)

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} treat_result.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "treat_result", label = TRUE) + NoLegend()
dev.off()


Epithelial <- Epithelial %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
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

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} group - harmony.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "group", label = TRUE, reduction = "harmonyumap")
dev.off() 



reduction.use <- 'harmony' # 'harmony' 'pca'
### cluster -------
Epithelial <- 
  FindNeighbors(
    Epithelial ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Epithelial <-
  FindClusters(
    Epithelial,
    resolution = c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.5)
  )

# dimplot
pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()




Epithelial_major_marker_ls <- c(
  'AT1 cells' = "AGER", # 
  'AT1 cells' = "SEMA3B", # 
  
  'AT2 cells' = "ABCA3", #  
  'AT2 cells' = "SFTPA1", #  
  'AT2 cells' = "SFTPB", #  
  'AT2 cells' = "SFTPC", #  
  'AT2 cells' = "HHIP", #  
  
  'Basal cells' = "KRT5", # 
  'Basal cells' = "NGFR", # 
  'Basal cells' = "TP63", # 
  'Basal cells' = "KRT15", # 
  
  'Deuterosomal cells' = 'CDC208', 
  
  'Goblet cells' = 'MUC5AC', 
  'Goblet cells' = 'CEACAM5', 
  
  'Pulmonary neuroendocrine cells' = 'GRP', 
  'Pulmonary neuroendocrine cells' = 'CALCA', 
  
  'Tuft cells' = 'POU2F3', 
  'Tuft cells' = 'ASCL2', 
  
  'Lonocyes' = 'ASCL3', 
  'Lonocyes' = 'BSND', 
  'Lonocyes' = 'FOXI1', 
  
  'Hillock cells' = 'KRT14', 
  'Hillock cells' = 'KRT6A', 
  
  'Clumnar cells' = 'FOXJ1', 
  'Clumnar cells' = 'RSPH1', 
  'Clumnar cells' = 'CDHR3', 
  'Clumnar cells' = 'CDHR4', 
  
  'Mucous cells' = 'MUC5B', 
  'Mucous cells' = 'SPDEF', 
  'Mucous cells' = 'BPIFB2', 
  
  'Serous cells' = 'LYZ', 
  'Serous cells' = 'LTF', 
  
  'Secretory cells' = 'SCGB1A1', 
  'Secretory cells' = 'SCGB3A2', 
  
  'Tumor cell' = 'MIAT', 
  'Tumor cell' = 'DSG3', 
  'Tumor cell' = 'MUC4', 
  
  'Ciliated cells' = "AKAP14", # 
  'Ciliated cells' = "ALDH3B1", # 
  'Ciliated cells' = "ANKRD66", # 
  'Ciliated cells' = "C11orf97", # 
  'Ciliated cells' = "DNAI1", # 
  
  # 'Club cells' = "SCGB1A1", # 
  'Club cells' = "SCGB3A1", # 
  'Club cells' = "PIGR"
)

# find markers
Idents(Epithelial) <- "RNA_snn_res.1.5"
Epithelial_markers <-
  FindAllMarkers(
    JoinLayers(Epithelial),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

Epithelial_top5_markers <- 
  Epithelial_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>% 
  dplyr::filter(avg_log2FC > 0.75) %>% 
  dplyr::filter(pct.1 - pct.2 > 0.1) %>% 
  arrange(desc(avg_log2FC)) %>%
  # slice_head(n = 30) %>%
  # arrange(desc(pct.1)) %>%
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
pdf(file.path(Epithelial_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 35, height = 6)
DotPlot(Epithelial, group.by = "RNA_snn_res.1.5", features = major_marker_ls)+RotatedAxis()
DotPlot(Epithelial, group.by = "RNA_snn_res.1.5", features = Epithelial_major_marker_ls)+RotatedAxis()
DotPlot(Epithelial, group.by = "RNA_snn_res.1.5", features = Epithelial_top5_markers)+RotatedAxis()
# DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
# DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = B_major_marker_ls)+RotatedAxis()
# DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = macrophage_major_marker_ls)+RotatedAxis()
# DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = T_major_marker_ls)+RotatedAxis()
# DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = CD8_major_marker_ls)+RotatedAxis()
# DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = CD4_major_marker_ls)+RotatedAxis()
dev.off()

# DotPlot(Epithelial, group.by = "RNA_snn_res.0.8", features = Fibro_major_marker_ls2)+RotatedAxis()
VlnPlot(Epithelial, features = 'nCount_RNA', pt.size = 0)
pdf(file.path(Epithelial_out_dir, 'for REST tumor cell.pdf'), width = 8, height = 6)
VlnPlot(subset(Epithelial, REST > 0), features = 'REST', group.by = 'RNA_snn_res.1.5', pt.size = 0)
VlnPlot(Epithelial, features = 'nCount_RNA', pt.size = 0)
FeaturePlot(Epithelial, features = 'REST', reduction = 'harmonyumap')
DimPlot(Epithelial, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()




## annotation -----------
Epithelial$RNA_snn_res.1.5 %>% table()
table(Epithelial$RNA_snn_res.1.5, Epithelial$sample)
cell_cluster <- c(
  '0' = "Tumor cell", # MIAT
  '1' = "Tumor cell", # MIAT
  '2' = "Tumor cell", # MIAT
  '3' = "Tumor cell", # MIAT
  '4' = "Tumor cell", # MIAT
  '5' = "Tumor cell", # MIAT
  '6' = "Tumor cell", # MIAT
  '7' = "low-nCounts",
  '8' = "Tumor cell", # MIAT 
  '9' = "Tumor cell", # MIAT
  '10' = "Tumor cell", # MIAT
  '11' = "Tumor cell", # MIAT
  '12' = "Tumor cell", # MIAT
  '13' = "Tumor cell", # MIAT
  '14' = "REST_tumor cell", # DSG3
  '15' = "REST_tumor cell", # DSG3
  '16' = "low-nCounts",
  '17' = "REST_tumor cell",
  '18' = "REST_tumor cell", # MUC4, 肺鳞癌标志物
  '19' = "AT2 cell", # ABCA3
  '20' = "REST_tumor cell", # DSG3
  '21' = "REST_tumor cell", # DLEC1, DNAI1
  '22' = "AT1 cell", # AGER
  '23' = "Tumor cell", # MIAT
  '24' = "REST_tumor cell" # DSG3
)


Idents(Epithelial) <- "RNA_snn_res.1.5"

Epithelial <- RenameIdents(Epithelial, cell_cluster)
Epithelial$subtype <- Idents(Epithelial)

pdf(file.path(Epithelial_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 10)
DimPlot(Epithelial, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE)
dev.off()

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 7, width = 18)
DimPlot(Epithelial, group.by = "subtype",split.by = 'treat_result', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

# Epithelial$treat_result <- str_replace_all(Epithelial$sample, c(
#   'CW_SH' = 'non-PCR', 
#   'CW_SQ' = 'non-PCR', 
#   'HJD_SH' = 'non-PCR', 
#   'HJD_SQ' = 'non-PCR', 
#   'LGRO_SH' = 'non-PCR', 
#   'LGRO_SQ' = 'non-PCR', 
#   'WBS_SH' = 'non-PCR', 
#   'WBS_SQ' = 'non-PCR', 
#   'WLJ_SH' = 'non-PCR', 
#   'WLJ_SQ' = 'non-PCR', 
#   'JJW_SH' = 'PCR', 
#   'JJW_SQ' = 'PCR', 
#   'YBK_SH' = 'PCR', 
#   'YBK_SQ' = 'PCR'
# ))
# 
# pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} annotation split by treat_result.pdf")), height = 7, width = 18)
# DimPlot(Epithelial, group.by = "subtype",split.by = 'treat_result', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
# dev.off()
# 
# pdf(file.path(Epithelial_out_dir,str_glue("Dimplot dim{dims} annotation - PCR -SH.pdf")), height = 6, width = 8)
# DimPlot(subset(Epithelial, sample %in% c('JJW_SH', 'YBK_SH')), group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
# dev.off()




## cell props barplot  --------
Epithelial_meta <- Epithelial[[]]

Epithelial_meta$sample <- factor(Epithelial_meta$sample, levels=c(
  'CW_SQ', 
  'CW_SH', 
  'HJD_SQ',
  'HJD_SH',
  'LGRO_SQ',
  'LGRO_SH',
  'WBS_SQ',
  'WBS_SH',
  'WLJ_SQ',
  'WLJ_SH',
  'JJW_SQ',
  'JJW_SH',
  'YBK_SQ',
  'YBK_SH'
)
)

p_col_proj <-
  Epithelial_meta %>%
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
  filename = file.path(Epithelial_out_dir, "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  Epithelial_meta %>%
  group_by(treat_result, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = treat_result, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path(Epithelial_out_dir, "cell_freq_treat_result.pdf"),
  height = 6,
  width  = 4
)

cellType_marker_ls <- c(
  'Tumor cell' = 'MIAT', 
  'Tumor cell' = 'DSG3', 
  'Tumor cell' = 'MUC4', 
  
  'REST_tumor cell' = 'REST', 
  
  "AT2 cell" = 'ABCA3', 
  
  # 'Mix' = 'POSTN',
  # 'Mix' = 'EPCAM',
  # 
  # # "Ciliated cell" = 'DLEC1',
  # # "Ciliated cell" = 'DNAI1',
  
  "AT1 cell" = 'AGER'
)

pdf(file.path(Epithelial_out_dir, "dotplot celltype main_marker.pdf"), width = 10, height = 6)
DotPlot(Epithelial, group.by = "subtype", features = cellType_marker_ls, dot.scale = 10) + theme_90angle
dev.off()

pdf(file.path(Epithelial_out_dir, "dotplot celltype main_marker - pure.pdf"), width = 10, height = 5)
DotPlot(subset(Epithelial, subtype %in% c("Tumor cell", 'REST_tumor cell', "AT2 cell", "AT1 cell")), group.by = "subtype", features = cellType_marker_ls, dot.scale = 10) + theme_90angle
dev.off()

saveRDS(Epithelial, file.path(Epithelial_out_dir, 'Epithelial2.rds'))



# spatial plot -------------------
## 映射-------------------------------
library(tibble)
library(data.table)
all_data$annotation %>% table()
all_data$group %>% table()
subtype <- 
  Epithelial[['subtype']] %>% 
  rownames_to_column('cells')
all_data_SQ <- subset(all_data, group == 'SQ')

meta.data.t <- 
  subtype %>%
  right_join(
    all_data_SQ[['annotation']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation2 = ifelse((as.character(annotation) %in% c('Endothelial cell', 'Fibroblast', 'Macrophage', 'Plasma cell', 'T cell')),
                         as.character(annotation), 
                         as.character(subtype)))%>%
  dplyr::select(cells, annotation2) %>%
  column_to_rownames('cells')

all_data_SQ <- AddMetaData(all_data_SQ,meta.data.t['annotation2'])
all_data_SQ$annotation %>% table()
all_data_SQ$annotation2 %>% table()


## Spatial Dimplot -----------------
all_data_SQ_meta <- 
  all_data_SQ[[]] %>% 
  rownames_to_column('cell_id')
cells_boundary_ls <- readRDS('R_OUTPUT/data_for_analysis/cells_boundary_ls.rds')


ct_sp_out_dir <- file.path(Epithelial_out_dir, 'Spatial dimplot')
dir.create(ct_sp_out_dir)
ct_sp_p_ls <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){ # browser()
    if(sample %in% unique(all_data_SQ$sample)){
      anno_data   <- all_data_SQ_meta[c("sample", "treat_result")] %>% distinct()
      group_info  <- anno_data$treat_result[anno_data$sample == sample]
      
      dat <-
        cells_boundary %>%
        right_join(all_data_SQ_meta, by = "cell_id")
      
      p <- dat %>%
        ggplot(aes(x = vertex_x, y = vertex_y)) +
        geom_polygon(aes(fill = annotation2, group = cell_id))  +
        scale_y_reverse() +
        coord_equal() + 
        ggsci::scale_fill_igv() +
        theme_void() +
        labs(fill = NULL, title = paste(sample, group_info, sep = "_") )
    }
    else{
      return(NULL)
    }
    
    p
  }) %>%
  purrr::discard(function(x) length(x) < 1)

ct_sp_p_ls %>% 
  iwalk(function(p, name){
    if(all(is.na(p$data$vertex_x))) return(NULL)
    if(any(is.na(p$data$vertex_x))) warning("vertex_x has NAs ...", .immediate = T)
    pdf(
      file.path(ct_sp_out_dir, str_glue("{name} - spatial dimplot.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


ct_sp_out_dir2 <- file.path(Epithelial_out_dir, 'Spatial dimplot_split')
dir.create(ct_sp_out_dir2)
ct_sp_p_ls2 <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){
    if(sample %in% unique(all_data_SQ$sample)){
      anno_data   <- all_data_SQ_meta[c("sample", "treat_result")] %>% distinct()
      group_info  <- anno_data$treat_result[anno_data$sample == sample]
      
      dat <-
        cells_boundary %>%
        right_join(all_data_SQ_meta, by = "cell_id")
      
      ct_color <- structure(
        names = dat$annotation2 %>% unique(),
        ggsci::pal_igv()(20)
      ) %>%
        c(
          Others = "gray90"
        )
      
      p_ls <-
        dat$annotation2 %>% unique() %>% # .[1] %>%
        set_names(., .) %>%
        map(function(ct){
          dat2 <- dat
          
          dat2$annotation2 <- as.character(dat2$annotation2)
          dat2$annotation2[dat2$annotation2 != ct] <- "Others"
          
          p <- dat2 %>%
            ggplot(aes(x = vertex_x, y = vertex_y)) +
            geom_polygon(aes(fill = annotation2, group = cell_id))  +
            scale_y_reverse() +
            coord_equal() +
            # ggsci::scale_fill_igv() +
            scale_fill_manual(values = ct_color) +
            theme_void() +
            labs(fill = NULL, title = paste(sample, group_info, ct, sep = "_") )
          p
        })
    }
    else{
      return(NULL)
    }
    
    p_ls
  }) %>% 
  purrr::discard(function(x) length(x) < 1)

ct_sp_p_ls2 %>% 
  iwalk(function(p, name){ # browser()
    if(all(is.na(p[[1]]$data$vertex_x))) return(NULL)
    if(any(is.na(p[[1]]$data$vertex_x))) warning("vertex_x has NAs ...", .immediate = T)
    pdf(
      file.path(ct_sp_out_dir2, str_glue("{name} - spatial dimplot - split by annotation2.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })

all_data_SQ_meta$annotation2 %>% 
  unique() %>% 
  setNames(., .) %>% 
  # .[1] %>% 
  iwalk(function(ct, ct_name){
    pdf(
      file.path(ct_sp_out_dir2, str_glue("{ct_name} - spatial dimplot - split by sample.pdf")),
      width  = 20,
      height = 20
    )
    ct_sp_p_ls2 %>% 
      walk(function(p){ # browser()
        if(all(is.na(p[[1]]$data$vertex_x))) return(NULL)
        if(any(is.na(p[[1]]$data$vertex_x))) warning("vertex_x has NAs ...", .immediate = T)
        print(p[[ct_name]])
      })
    dev.off()
  })




# REST tumor annotation ----------
Epithelial_out_dir <-'Epithelial_SQ_subtype'
Epithelial <- readRDS("~/Epithelial_SQ_subtype/Epithelial.rds")
Epithelial$treat_result %>% table()

VlnPlot(subset(Epithelial, REST > 0 & RNA_snn_res.1.5 == '7'), features = 'REST', group.by = 'treat_result')
ggsave(file.path(Epithelial_out_dir, 'Vlnplot_REST_by_treat_result_cluster7.pdf'), width = 8, height = 6)

VlnPlot(subset(Epithelial, REST > 0 & RNA_snn_res.1.5 == '16'), features = 'REST', group.by = 'treat_result')
ggsave(file.path(Epithelial_out_dir, 'Vlnplot_REST_by_treat_result_cluster16.pdf'), width = 8, height = 6)

VlnPlot(subset(Epithelial, REST > 0 & RNA_snn_res.1.5 == '17'), features = 'REST', group.by = 'treat_result')
ggsave(file.path(Epithelial_out_dir, 'Vlnplot_REST_by_treat_result_cluster17.pdf'), width = 8, height = 6)
