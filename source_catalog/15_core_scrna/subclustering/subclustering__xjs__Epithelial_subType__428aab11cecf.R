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

Epithelial <- subset(all_data, annotation == 'Epithelial/Tumor cell')

## PCA --------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3) # 80Gb

DefaultAssay(Epithelial) <- "RNA"
Epithelial <- FindVariableFeatures(Epithelial)
Epithelial <- ScaleData(Epithelial, vars.to.regress = 'nCount_RNA')
Epithelial <- RunPCA(Epithelial)

Epithelial_out_dir <-'Epithelial_subtype' # dir.create(Epithelial_out_dir)


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

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} group.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "group", label = TRUE) + NoLegend()
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
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
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
Idents(Epithelial) <- "RNA_snn_res.0.6"
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
DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()
DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = Epithelial_major_marker_ls)+RotatedAxis()
DotPlot(Epithelial, group.by = "RNA_snn_res.0.6", features = Epithelial_top5_markers)+RotatedAxis()
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
VlnPlot(subset(Epithelial, REST > 0), features = 'REST', group.by = 'RNA_snn_res.0.8', pt.size = 0)
VlnPlot(Epithelial, features = 'nCount_RNA', pt.size = 0)
FeaturePlot(Epithelial, features = 'REST', reduction = 'harmonyumap')
DimPlot(Epithelial, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()




## inferCNV ---------------------------------------------------------------------
### set future -------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 50 * 1024 ^ 3) # 50Gb

library(infercnv)
infercnv_out_dir = "Epithelial_subtype/inferCNV_episub" # dir.create(infercnv_out_dir) 

### 映射-------------------------------
library(tibble)
library(data.table)
subtype <- 
  Epithelial[['RNA_snn_res.0.6']] %>% 
  rownames_to_column('cells') %>% 
  mutate(cluster = str_replace(.$RNA_snn_res.0.6, '^', 'cluster_'))

meta.data.t <- 
  subtype %>%
  right_join(
    all_data[['annotation']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation2 = ifelse((as.character(annotation) %in% c('Endothelial cell', 'Fibroblast', 'Macrophage', 'Plasma cell', 'T cell')),
                         as.character(annotation), 
                         as.character(cluster)))%>%
  dplyr::select(cells, annotation2) %>%
  column_to_rownames('cells')


all_data<-AddMetaData(all_data,meta.data.t['annotation2'])
all_data$annotation %>% table()
all_data$annotation2 %>% table()

meta_data_infer <- 
  all_data[[c('annotation2', 'sample')]] %>% 
  rownames_to_column('cell_id')

meta_data_infer_ls <- 
  meta_data_infer %>% 
  group_by(sample) %>%
  nest() %>% 
  deframe()


# cell is too much
# less than 50000
set.seed(1234)
meta_data_infer_subset_ls <-
  meta_data_infer_ls %>% 
  imap(function(meta_data_infer, sample){
    meta_data_infer2 <- 
      meta_data_infer %>% 
      group_by(annotation2) %>%
      slice_sample(prop = 50000/dim(meta_data_infer)[1])
  })

# export as file
meta_data_infer_subset_ls %>% 
  imap(function(meta_data_infer_subset, sample){
    if(!file.exists(file.path(infercnv_out_dir,sample))){
      dir.create(file.path(infercnv_out_dir,sample))
    }
    data.table::fwrite(
      meta_data_infer_subset, 
      file = file.path(infercnv_out_dir, sample, "meta_data_infer_subset.txt"),
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t"
    )
  })



infercnv_obj_ls <- 
  meta_data_infer_subset_ls %>% 
  imap(function(meta_data_infer_subset, sample){ # browser()
    infercnv_obj <- CreateInfercnvObject(
      raw_counts_matrix = JoinLayers(all_data)@assays$RNA$counts[VariableFeatures(Epithelial)[1:1000], meta_data_infer_subset$cell_id],
      #raw_counts_matrix = JoinLayers(all_data)@assays$RNA$counts[VariableFeatures(all_data)[1:500], meta_data_infer_subset$cell_id],
      annotations_file  = file.path(infercnv_out_dir, sample, "meta_data_infer_subset.txt"),
      gene_order_file   = file.path(infercnv_out_dir, "gene_pos47.txt"),
      delim             = "\t",
      min_max_counts_per_cell = c(30,Inf),
      ref_group_names = c('Macrophage', 'Plasma cell', 'T cell')
    )
    infercnv_obj
  })



infercnv_obj_ls[[2]]@count.data %>% rowMeans %>% { sum(. > 0.01) } # 788 <--- 
infercnv_obj_ls[[2]]@count.data %>% rowMeans %>% { sum(. < 0.01) } # 174

infercnv_obj_ls2 <- 
  infercnv_obj_ls %>% 
  imap(function(infercnv_obj, sample){
    infercnv_obj2 <- infercnv::run(
      infercnv_obj,
      cutoff              =    0.01,
      out_dir    =   file.path(infercnv_out_dir, sample),
      cluster_by_groups   =    TRUE,
      plot_steps          =    FALSE,
      denoise             =    TRUE,
      HMM                 =    FALSE,
      no_prelim_plot      =    TRUE,
      png_res             =    300,
      # useRaster   = FALSE,
      output_format = 'pdf',
      no_plot = T #  figure margin too large 
    )
    infercnv_obj2
  })

infercnv_obj_ls2 %>% 
  imap(function(infercnv_obj2, sample){
    infercnv::plot_cnv(
      infercnv_obj2,
      out_dir           =  file.path(infercnv_out_dir, sample),
      output_filename   =  'infercnv',
      output_format     =  "pdf",
      png_res = 300,
      cluster_by_groups =  TRUE,
      x.range           =  "auto",
      x.center          =  1,
      title             =  "infercnv",
      color_safe_pal    =  FALSE
    )
  })

infercnv_obj2_default_median_filtering_ls <-
  infercnv_obj_ls2 %>%
  map(function(infercnv_obj2){
    infercnv::apply_median_filtering(infercnv_obj2)
  })


infercnv_obj2_default_median_filtering_ls %>% 
  imap(function(infercnv_obj2_default_median_filtering, sample){
    infercnv::plot_cnv(
      infercnv_obj2_default_median_filtering, 
      out_dir           =  file.path(infercnv_out_dir, sample),
      output_filename   =  'infercnv_median_filtering',
      output_format="pdf",
      cluster_by_groups =  TRUE,
      x.range           =  "auto",
      x.center          =  1,
      title             =  "infercnv",
      color_safe_pal    =  FALSE
    )
  })

saveRDS(infercnv_obj2_default_median_filtering_ls, file.path(infercnv_out_dir, 'infercnv_obj2_default_median_filtering_ls.rds'))

saveRDS(infercnv_obj_ls2, file.path(infercnv_out_dir, 'infercnv_obj_ls2.rds'))
saveRDS(infercnv_obj_ls, file.path(infercnv_out_dir, 'infercnv_obj_ls.rds'))

### plot -------------------------------------------------------------------------
library(scales)
cnvScore <- function(data){
  sd_num <- apply(data, 1, sd)
  data <- data[sd_num > median(sd_num), ]
  
  rng <- quantile(data, c(0.01, 0.99))
  data <- pmax( pmin(data, rng[2]), rng[1] )
  
  data <- data %>% as.matrix() %>%
    t() %>% 
    scale() %>% 
    t()
  
  cnv_score <- as.data.frame(colSums(data * data))
  return(cnv_score)
}

cnv_score_ls <- 
  infercnv_obj_ls2 %>% 
  map(function(infercnv_obj2){
    cnv_score <- cnvScore(infercnv_obj2@expr.data)
    colnames(cnv_score)<-'CNV_score'
    cnv_score
  })


temp_ls <- 
  cnv_score_ls %>%
  imap(function(cnv_score, sample){
    cnv_score %>% 
      rownames_to_column('cell_id') %>% 
      left_join(.,as.data.frame(meta_data_infer_subset_ls[[sample]]),'cell_id') -> temp
    temp
  })


plot_ls <- 
  temp_ls %>% 
  map(function(temp){
    temp %>% 
      ggplot(aes(annotation2,log(CNV_score)))+
      geom_violin(aes(fill= annotation2),color="NA",drop = FALSE, trim = FALSE)+
      ggsci::scale_fill_igv()+
      # scale_fill_manual(values = color_v)+
      theme_bw()  +
      RotatedAxis()
  })

plot_ls %>% 
  iwalk(function(plot, sample){
    ggsave(
      plot,
      filename = file.path(infercnv_out_dir, sample, "CNV_score.pdf"),
      height = 8,
      width  = 15
    )
  })


plot_ls2 <-
  temp_ls %>%
  map(function(temp){
    temp %>% 
      group_by(annotation2) %>%
      
      ggplot(aes(x = annotation2, y = log(CNV_score), fill = annotation2)) +
      geom_boxplot(outlier.colour = NA) +
      ggsci::scale_fill_igv() +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      )  
  })


plot_ls2 %>% 
  iwalk(function(plot2, sample){
    ggsave(
      plot2,
      filename = file.path(infercnv_out_dir, sample, "CNV_score_barplot.pdf"),
      height = 6,
      width  = 8
    )
  })


cnvScore <- function(data){
  sd_num <- apply(data, 1, sd)
  data <- data[sd_num > median(sd_num), ]
  
  data <- data %>% as.matrix() %>%
    t() %>% 
    scale() %>% 
    t()
  
  cnv_score <- as.data.frame(colSums(data * data))
  return(cnv_score)
}

cnv_score_noQuantile_ls <- 
  infercnv_obj_ls2 %>% 
  map(function(infercnv_obj2){
    cnv_score <- cnvScore(infercnv_obj2@expr.data)
    colnames(cnv_score)<-'CNV_score'
    cnv_score
  })


temp_noQuantile_ls <- 
  cnv_score_noQuantile_ls %>%
  imap(function(cnv_score, sample){
    cnv_score %>% 
      rownames_to_column('cell_id') %>% 
      left_join(.,as.data.frame(meta_data_infer_subset_ls[[sample]]),'cell_id') -> temp
    temp
  })


plot_noQuantile_ls <- 
  temp_noQuantile_ls %>% 
  map(function(temp){
    temp %>% 
      ggplot(aes(annotation2,log(CNV_score)))+
      geom_violin(aes(fill= annotation2),color="NA",drop = FALSE, trim = FALSE)+
      ggsci::scale_fill_igv()+
      # scale_fill_manual(values = color_v)+
      theme_bw()  +
      RotatedAxis()
  })

plot_noQuantile_ls %>% 
  iwalk(function(plot, sample){
    ggsave(
      plot,
      filename = file.path(infercnv_out_dir, sample, "CNV_score_noQuantile.pdf"),
      height = 8,
      width  = 15
    )
  })


plot_noQuantile_ls2 <-
  temp_noQuantile_ls %>%
  map(function(temp){
    temp %>% 
      group_by(annotation2) %>%
      
      ggplot(aes(x = annotation2, y = log(CNV_score), fill = annotation2)) +
      geom_boxplot(outlier.colour = NA) +
      ggsci::scale_fill_igv() +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      )  
  })


plot_noQuantile_ls2 %>% 
  iwalk(function(plot2, sample){
    ggsave(
      plot2,
      filename = file.path(infercnv_out_dir, sample, "CNV_score_barplot_noQuantile.pdf"),
      height = 6,
      width  = 8
    )
  })




## annotation -----------
Epithelial$RNA_snn_res.0.6 %>% table()
table(Epithelial$RNA_snn_res.0.6, Epithelial$sample)
cell_cluster <- c(
  '0' = "Tumor cell", # MIAT https://pmc.ncbi.nlm.nih.gov/articles/PMC11910686/
  '1' = "Tumor cell", # STRIP2 https://news.bioon.com/article/6651e5992166.html
  '2' = "Tumor cell", # MIAT https://pmc.ncbi.nlm.nih.gov/articles/PMC11910686/
  '3' = "low-nCounts",
  '4' = "Tumor cell", # CALB2, https://www.bafang.org.cn/genes/99281.html#:~:text=%E4%BE%8B%E5%A6%82%EF%BC%8C%E5%9C%A8%E6%9F%90%E4%BA%9B%E7%B1%BB%E5%9E%8B%E7%9A%84%E8%82%BA%E7%99%8C%E5%92%8C%E7%BB%93%E7%9B%B4%E8%82%A0%E7%99%8C%E4%B8%AD%EF%BC%8CCALB2%E7%9A%84%E8%A1%A8%E8%BE%BE%E5%8F%AF%E8%83%BD%E4%B8%8E%E8%82%BF%E7%98%A4%E7%9A%84%E4%BE%B5%E8%A2%AD%E6%80%A7%E5%92%8C%E9%A2%84%E5%90%8E%E7%9B%B8%E5%85%B3%E3%80%82%20%E5%8F%91%E8%82%B2%E5%BC%82%E5%B8%B8%EF%BC%9ACALB2%E5%9F%BA%E5%9B%A0%E7%9A%84%E5%8F%98%E5%BC%82%E5%8F%AF%E8%83%BD%E4%B8%8E%E6%9F%90%E4%BA%9B%E5%8F%91%E8%82%B2%E5%BC%82%E5%B8%B8%E6%9C%89%E5%85%B3%EF%BC%8C%E5%A6%82%E5%85%88%E5%A4%A9%E6%80%A7%E8%84%91%E5%8F%91%E8%82%B2%E5%BC%82%E5%B8%B8%EF%BC%8C%E8%BF%99%E4%BA%9B%E5%8F%98%E5%BC%82%E5%8F%AF%E8%83%BD%E5%BD%B1%E5%93%8D%E7%A5%9E%E7%BB%8F%E5%85%83%E7%9A%84%E6%AD%A3%E5%B8%B8%E5%8F%91%E8%82%B2%E5%92%8C%E5%8A%9F%E8%83%BD%E3%80%82%20%E5%9C%A8,%E9%81%97%E4%BC%A0%E5%92%A8%E8%AF%A2%20%E4%B8%AD%EF%BC%8C%E4%BA%86%E8%A7%A3CALB2%E5%9F%BA%E5%9B%A0%E7%9A%84%E5%8F%98%E5%BC%82%E5%AF%B9%E4%BA%8E%E8%AF%84%E4%BC%B0%E7%A5%9E%E7%BB%8F%E7%B3%BB%E7%BB%9F%E7%96%BE%E7%97%85%E7%9A%84%E9%A3%8E%E9%99%A9%E5%85%B7%E6%9C%89%E9%87%8D%E8%A6%81%E6%84%8F%E4%B9%89%E3%80%82%20%E9%80%9A%E8%BF%87%20%E5%9F%BA%E5%9B%A0%E6%A3%80%E6%B5%8B%EF%BC%8C%E5%8F%AF%E4%BB%A5%E8%AF%86%E5%88%AB%E5%87%BA%E5%8F%AF%E8%83%BD%E5%AF%BC%E8%87%B4%E7%96%BE%E7%97%85%E7%9A%84%E7%89%B9%E5%AE%9A%E5%8F%98%E5%BC%82%EF%BC%8C%E4%BB%8E%E8%80%8C%E4%B8%BA%E6%82%A3%E8%80%85%E6%8F%90%E4%BE%9B%E4%B8%AA%E6%80%A7%E5%8C%96%E7%9A%84%E9%A2%84%E9%98%B2%E3%80%81%E8%AF%8A%E6%96%AD%E5%92%8C%E6%B2%BB%E7%96%97%E5%BB%BA%E8%AE%AE%E3%80%82
  '5' = "Tumor cell", # MAN1A1 https://www.sohu.com/a/903059746_99992284
  '6' = "Tumor cell", # EGLN3, https://pubmed.ncbi.nlm.nih.gov/39021988/
  '7' = "Tumor cell", # ASPM, https://www.ebiotrade.com/newsf/2025-9/20250906001230793.htm#:~:text=%E9%80%9A%E8%BF%87%E7%94%9F%E7%89%A9%E4%BF%A1%E6%81%AF%E5%AD%A6%E5%88%86%E6%9E%90%E5%92%8C%E5%8A%9F%E8%83%BD%E5%AE%9E%E9%AA%8C%E8%AF%81%E5%AE%9E%EF%BC%8CASPM%E9%80%9A%E8%BF%87%E8%B0%83%E6%8E%A7%E8%A1%A8%E7%9A%AE%E7%94%9F%E9%95%BF%E5%9B%A0%E5%AD%90%E5%8F%97%E4%BD%93%20%28EGFR%29%E8%9B%8B%E7%99%BD%E7%A8%B3%E5%AE%9A%E6%80%A7%EF%BC%8C%E5%A2%9E%E5%BC%BA%E8%82%BA%E8%85%BA%E7%99%8C%20%28LUAD%29%E7%BB%86%E8%83%9E%E7%9A%84%E5%A2%9E%E6%AE%96%E3%80%81%E4%BE%B5%E8%A2%AD%E5%92%8C%E8%80%90%E8%8D%AF%E6%80%A7%E3%80%82,ASPM%E5%9F%BA%E5%9B%A0%E6%B2%89%E9%BB%98%E5%8F%AF%E6%98%BE%E8%91%97%E6%81%A2%E5%A4%8D%E8%80%90%E8%8D%AF%E7%BB%86%E8%83%9E%E5%AF%B9%E7%AC%AC%E4%B8%89%E4%BB%A3EGFR-TKI%E7%9A%84%E6%95%8F%E6%84%9F%E6%80%A7%EF%BC%8C%E4%B8%BA%E5%85%8B%E6%9C%8D%E4%B8%B4%E5%BA%8A%E8%80%90%E8%8D%AF%E6%8F%90%E4%BE%9B%E4%BA%86%E6%96%B0%E9%9D%B6%E7%82%B9%E3%80%82%20ASPM%E5%9C%A8%E8%82%BA%E7%99%8C%E7%BB%84%E7%BB%87%E4%B8%AD%E7%9A%84%E9%AB%98%E8%A1%A8%E8%BE%BE%E4%B8%8E%E4%B8%8D%E8%89%AF%E9%A2%84%E5%90%8E%E6%98%BE%E8%91%97%E7%9B%B8%E5%85%B3%20%E7%94%9F%E7%89%A9%E4%BF%A1%E6%81%AF%E5%AD%A6%E5%88%86%E6%9E%90%E6%98%BE%E7%A4%BA%EF%BC%8CASPM%E5%9C%A833%E7%A7%8D%E5%AE%9E%E4%BD%93%E8%82%BF%E7%98%A4%E4%B8%AD%E5%BC%82%E5%B8%B8%E9%AB%98%E8%A1%A8%E8%BE%BE%EF%BC%8C%E5%B0%A4%E5%85%B6%E5%9C%A8%E8%82%BA%E8%85%BA%E7%99%8C%20%28LUAD%29%E7%BB%84%E7%BB%87%E4%B8%AD%E8%A1%A8%E8%BE%BE%E6%B0%B4%E5%B9%B3%E6%98%BE%E8%91%97%E9%AB%98%E4%BA%8E%E7%99%8C%E6%97%81%E7%BB%84%E7%BB%87%E3%80%82
  '8' = "Tumor cell", # MIAT https://pmc.ncbi.nlm.nih.gov/articles/PMC11910686/ 
  '9' = "Tumor cell", # SPHK1, https://www.sohu.com/a/888842370_121798269#:~:text=%E4%BF%83%E8%82%BF%E7%98%A4%E4%BD%9C%E7%94%A8%EF%BC%9ASPHK1%E5%9C%A8%E5%A4%9A%E7%A7%8D%E7%99%8C%E7%97%87%E4%B8%AD%E9%AB%98%E8%A1%A8%E8%BE%BE%EF%BC%88%E5%A6%82%E8%82%BA%E7%99%8C%E3%80%81%E7%BB%93%E8%82%A0%E7%99%8C%E3%80%81%E5%8D%B5%E5%B7%A2%E7%99%8C%EF%BC%89%EF%BC%8C%E9%80%9A%E8%BF%87S1P%E4%BF%83%E8%BF%9B%E8%82%BF%E7%98%A4%E7%BB%86%E8%83%9E%E5%A2%9E%E6%AE%96%E3%80%81%E4%BE%B5%E8%A2%AD%E3%80%81%E8%A1%80%E7%AE%A1%E7%94%9F%E6%88%90%E5%92%8C%E5%85%8D%E7%96%AB%E9%80%83%E9%80%B8%E3%80%82,%E8%82%BA%E7%99%8C%EF%BC%9ASPHK1%E7%BC%BA%E5%A4%B1%E5%8F%AF%E5%AF%BC%E8%87%B4%E8%84%82%E8%B4%A8%E4%BB%A3%E8%B0%A2%E7%B4%8A%E4%B9%B1%EF%BC%8C%E6%8B%9B%E5%8B%9FM2%E5%9E%8B%E5%B7%A8%E5%99%AC%E7%BB%86%E8%83%9E%EF%BC%8C%E4%BF%83%E8%BF%9B%E7%82%8E%E7%97%87%E5%BE%AE%E7%8E%AF%E5%A2%83%EF%BC%8C%E5%8A%A0%E9%80%9FKRAS%E7%AA%81%E5%8F%98%E8%82%BA%E8%85%BA%E7%99%8C%E8%BF%9B%E5%B1%95%E3%80%82%20%E7%BB%93%E8%82%A0%E7%99%8C%EF%BC%9ASPHK1%E9%80%9A%E8%BF%87%E4%B8%8A%E8%B0%83VEGF%E8%A1%A8%E8%BE%BE%EF%BC%8C%E4%BF%83%E8%BF%9B%E8%A1%80%E7%AE%A1%E7%94%9F%E6%88%90%E6%8B%9F%E6%80%81%EF%BC%88VM%EF%BC%89%E5%BD%A2%E6%88%90%EF%BC%8C%E5%A2%9E%E5%BC%BA%E8%82%BF%E7%98%A4%E4%BE%B5%E8%A2%AD%E6%80%A7%E3%80%82
  '10' = "Anti-tumor cell", # FAM107A,https://pmc.ncbi.nlm.nih.gov/articles/PMC9163664/  HSD17B2,https://pmc.ncbi.nlm.nih.gov/articles/PMC7170045/  SERPINA3,https://pmc.ncbi.nlm.nih.gov/articles/PMC10552721/
  '11' = "Mucous cell", # MUC5B
  '12' = "Basal cell", # TP63
  '13' = "AT1 cell", #AGER # ROS1 https://m.baidu.com/bh/m/detail/ar_10096883026808480312
  '14' = "Mix", # 
  '15' = "Ciliated cell", # DLEC1, DNAI1
  '16' = "Tumor cell", # MIAT https://pmc.ncbi.nlm.nih.gov/articles/PMC11910686/
  '17' = "low-nCounts",
  '18' = "low-nCounts", 
  '19' = "Tumor cell", # SLC2A1, https://www.sciencedirect.com/science/article/abs/pii/S0010482524002671
  '20' = "Tumor cell", # SLC7A7, https://www.xjishu.com/zhuanli/05/202411917427.html
  '21' = "low-nCounts", 
  '22' = "low-nCounts"
)


Idents(Epithelial) <- "RNA_snn_res.0.6"

Epithelial <- RenameIdents(Epithelial, cell_cluster)
Epithelial$subtype <- Idents(Epithelial)

pdf(file.path(Epithelial_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 10)
DimPlot(Epithelial, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE)
dev.off()

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 7, width = 18)
DimPlot(Epithelial, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

Epithelial$treat_result <- str_replace_all(Epithelial$sample, c(
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

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} annotation split by treat_result.pdf")), height = 7, width = 18)
DimPlot(Epithelial, group.by = "subtype",split.by = 'treat_result', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Epithelial_out_dir,str_glue("Dimplot dim{dims} annotation - PCR -SH.pdf")), height = 6, width = 8)
DimPlot(subset(Epithelial, sample %in% c('JJW_SH', 'YBK_SH')), group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()




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
  filename = file.path("Epithelial_subtype/cell_freq_sample.pdf"),
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
  filename = file.path("Epithelial_subtype/cell_freq_treat_result.pdf"),
  height = 6,
  width  = 4
)

cellType_marker_ls <- c(
  "Tumor cell" = 'MIAT', 
  # "Tumor cell" = 'STRIP2', 
  # "Tumor cell" = 'CALB2', 
  # "Tumor cell" = 'MAN1A1', 
  # "Tumor cell" = 'EGLN3', 
  # "Tumor cell" = 'ASPM', 
  # "Tumor cell" = 'SPHK1', 
  # "Tumor cell" = 'SLC2A1', 
  # "Tumor cell" = 'SLC7A7', 
   
  "Anti-tumor cell" = 'FAM107A', 
  "Anti-tumor cell" = 'HSD17B2', 
  "Anti-tumor cell" = 'SERPINA3', 
  
  "Mucous cell" = 'MUC5B', 
  "Basal cell" = 'TP63', 
  "AT1 cell" = 'AGER', 
  
  # 'Mix' = 'POSTN',
  # 'Mix' = 'EPCAM',
  
  "Ciliated cell" = 'DLEC1', 
  "Ciliated cell" = 'DNAI1'
)

pdf(file.path(Epithelial_out_dir, "dotplot celltype main_marker.pdf"), width = 20, height = 8)
DotPlot(Epithelial, group.by = "subtype", features = cellType_marker_ls, dot.scale = 10) + theme_90angle
dev.off()

pdf(file.path(Epithelial_out_dir, "dotplot celltype main_marker - pure.pdf"), width = 12, height = 7)
DotPlot(subset(Epithelial, subtype %in% c("Tumor cell", "Anti-tumor cell", "Mucous cell", "Basal cell", "AT1 cell", "Ciliated cell")), group.by = "subtype", features = cellType_marker_ls, dot.scale = 10) + theme_90angle
dev.off()

saveRDS(Epithelial, 'Epithelial_subtype/Epithelial.rds')

pdf(file.path(Epithelial_out_dir, "vlnplot REST-Tumor.pdf"), width = 10, height = 6)
VlnPlot(Epithelial, features = 'REST', group.by = 'subtype')
VlnPlot(Epithelial, features = 'REST', group.by = 'RNA_snn_res.0.6')
VlnPlot(subset(Epithelial, REST > 0), features = 'REST', group.by = 'RNA_snn_res.0.6', pt.size = 0)
dev.off()




# spatial plot -------------------
## 映射-------------------------------
library(tibble)
library(data.table)
all_data$annotation %>% table()
subtype <- 
  Epithelial[['subtype']] %>% 
  rownames_to_column('cells')

meta.data.t <- 
  subtype %>%
  right_join(
    all_data[['annotation']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation2 = ifelse((as.character(annotation) %in% c('Endothelial cell', 'Fibroblast', 'Macrophage', 'Plasma cell', 'T cell')),
                         as.character(annotation), 
                         as.character(subtype)))%>%
  dplyr::select(cells, annotation2) %>%
  column_to_rownames('cells')

all_data<-AddMetaData(all_data,meta.data.t['annotation2'])
all_data$annotation %>% table()
all_data$annotation2 %>% table()


## Spatial Dimplot -----------------
all_data_meta <- 
  all_data[[]] %>% 
  rownames_to_column('cell_id')
cells_boundary_ls <- readRDS('R_OUTPUT/data_for_analysis/cells_boundary_ls.rds')


ct_sp_out_dir <- file.path(Epithelial_out_dir, 'Spatial dimplot')
dir.create(ct_sp_out_dir)
ct_sp_p_ls <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){ # browser()
    anno_data   <- all_data_meta[c("sample", "treat_result")] %>% distinct()
    group_info  <- anno_data$treat_result[anno_data$sample == sample]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation2, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() + 
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL, title = paste(sample, group_info, sep = "_") )
    p
  })

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
    anno_data   <- all_data_meta[c("sample", "treat_result")] %>% distinct()
    group_info  <- anno_data$treat_result[anno_data$sample == sample]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
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
    
    p_ls
  })

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

all_data_meta$annotation2 %>% 
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




## Spatial Featureplot ----------------
col_val <- colorRampPalette(c("blue", "red")) (50)
alpha_val <- seq(0, 1, length.out = 50)
alpha(col_val, alpha_val) %>% scales::show_col()
legend_col <- alpha(col_val, alpha_val)

gene_sp_out_dir <- file.path(Epithelial_out_dir, 'Spatial featureplot')
dir.create(gene_sp_out_dir)

gene_cus <- c('REST', 'HTR1A')
dat_expr <- all_data %>% 
  FetchData(vars = unlist(gene_cus), layer = "data") %>% 
  rownames_to_column('cell_id')
stopifnot(all(all_data_meta$cell_id == dat_expr$cell_id))
dat_expr_coords_full <- left_join(all_data_meta, dat_expr, by = 'cell_id')

marker_p_ls <- 
  gene_cus %>%
  setNames(., .) %>% 
  map(function(x) {
    p_combined <-
      cells_boundary_ls %>%
      imap(function(cells_boundary, sample) { # browser()
        if (!x %in% colnames(dat_expr_coords_full)) return(NULL)
    
        dat <-
          cells_boundary %>%
          right_join(dat_expr_coords_full, by = "cell_id")
        
        min_max <- c(0, 1) # quantile(dat[[x]], c(0.01, 0.99))
        
        dat[[x]][ dat[[x]] < min_max[1] ] <- min_max[1]
        dat[[x]][ dat[[x]] > min_max[2] ] <- min_max[2]
        
        p <- dat %>%
          arrange(.data[[x]]) %>% 
          as.data.frame() %>%
          ggplot(aes(x = vertex_x, y = -vertex_y)) +
          geom_polygon(aes(fill = .data[[x]], group = cell_id, alpha = .data[[x]])) + 
          guides(alpha = "none") +
          theme_void() +
          coord_equal() +
          ggtitle(sample) 
        
        # all_values <- unlist(lapply(p_combined$patches$plots, function(p) p$data[[p$labels$colour]]))
        p <- p + scale_fill_gradientn(
          colours = legend_col,
          # limits = c(0, ceiling(max(dat_expr[-1]))),  # 统一颜色范围
          transform = "log1p", 
          # limits = c(0, 1), 
          name = "Expression"
        )
        
        p
      }) %>%
      purrr::discard(is.null) # %>% 
      # patchwork::wrap_plots(ncol = 1, guides = "collect") &
      # theme(legend.position = "right")  # 全局设置图例位置
    
    p_combined
  })


marker_p_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path(gene_sp_out_dir,  
                  str_glue("{grp} - spatial feature plot - sameLegend_red2blue.pdf")), 
        height = 8, 
        width = 8)
    print(p)
    dev.off()
  })    


gene_cus[1] %>%
  setNames(., .) %>% 
  map(function(x) {
    p_combined <-
      cells_boundary_ls[1] %>%
      imap(function(cells_boundary, sample) { # browser()
        if (!x %in% colnames(dat_expr_coords_full)) return(NULL)
        
        dat <-
          cells_boundary %>%
          right_join(dat_expr_coords_full, by = "cell_id")
        
        min_max <- c(0, 1) # quantile(dat[[x]], c(0.01, 0.99))
        
        dat[[x]][ dat[[x]] < min_max[1] ] <- min_max[1]
        dat[[x]][ dat[[x]] > min_max[2] ] <- min_max[2]
        
        p <- dat %>%
          arrange(.data[[x]]) %>% 
          as.data.frame() %>%
          ggplot(aes(x = vertex_x, y = -vertex_y)) +
          geom_polygon(aes(fill = .data[[x]], group = cell_id, alpha = .data[[x]])) + 
          guides(alpha = 'none') +
          theme_void() +
          coord_equal() +
          ggtitle(sample) 
        
        # all_values <- unlist(lapply(p_combined$patches$plots, function(p) p$data[[p$labels$colour]]))
        p <- p + scale_fill_gradientn(
          colours = legend_col,
          # limits = c(0, ceiling(max(dat_expr[-1]))),  # 统一颜色范围
          transform = "log1p", 
          # limits = c(0, 1), 
          name = "Expression"
        )
        
        p
      }) %>%
      purrr::discard(is.null) # %>% 
    # patchwork::wrap_plots(ncol = 1, guides = "collect") &
    # theme(legend.position = "right")  # 全局设置图例位置
    
    p_combined
  })
