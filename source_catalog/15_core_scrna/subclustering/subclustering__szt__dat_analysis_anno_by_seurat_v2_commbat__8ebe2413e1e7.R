library(Seurat)
library(tidyverse)

options(Seurat.memsafe = TRUE) # more GC frequently

options(future.globals.maxSize = 50 * 1024 ^ 3) # 10Gb
future::plan("multisession", workers = 10, gc = TRUE)
# future::plan("sequential")
# future:::ClusterRegistry("stop")


out_dir_cluster <- "dat_analysis_cluster_remove_DAPI_batch"
dir.create(out_dir_cluster)

low_quality_roi <- c(
  # 'A1', 'A2', 'A5', 'B5', 'A8', 'A14'
)

# NOT RUN. raw data from qupath -------------------
file_ls <- list.files("data/results2", pattern = "FWSP.*csv", full.names = TRUE)
names(file_ls) <- file_ls %>% basename() %>% str_remove("\\.csv$") %>% str_split('#1_') %>% map_chr(2)

all_dat_ls_full <- 
  file_ls %>% # [2:6] %>%
  imap(function(x, name){
    print(name)
    
    df <- data.table::fread(x) %>%
      as.data.frame() %>%
      rownames_to_column('cell_id') %>%
      mutate(
        ROI_Name = name,
        cell_id  = str_c(ROI_Name, cell_id, sep = "_")
      ) %>%
      column_to_rownames("cell_id")
    
    df
  })

all_dat_ls <- all_dat_ls_full[!names(all_dat_ls_full) %in% low_quality_roi]
all_dat_raw <- bind_rows(all_dat_ls )
all_dat_raw$cell_id <- rownames(all_dat_raw)

# measurement 
pattern_used <- "Cell: .*Mean$" # "Nucleus: .*Mean$" # "Cell: .*Mean$"
channel_columns <- all_dat_raw %>% colnames() %>% str_subset(pattern_used)

dat_meas_expr <- all_dat_raw[channel_columns]
colnames(dat_meas_expr) <- 
  dat_meas_expr %>% 
  colnames() %>%
  str_split(":") %>%
  map_chr(2) %>%
  str_trim()

dat_meas_expr[is.na(dat_meas_expr)] <- 0

colnames(dat_meas_expr) 


#  DAPI QC 
pdf(
  file.path(out_dir_cluster, "density plot of DAPI.pdf"),
  height = 6,
  width = 8
)
dat_meas_expr$DAPI %>% density() %>% plot()
abline(v = exp(7), col = "red")
dat_meas_expr$DAPI %>% log %>% density() %>% plot(main = "DAPI")
abline(v = 7, col = "red")
dev.off()

remove_dapi_cells <- rownames(dat_meas_expr)[dat_meas_expr$DAPI < exp(7)]

# remove DAPI
if ('DAPI' %in% colnames(dat_meas_expr)){
  dat_meas_expr <- dat_meas_expr[colnames(dat_meas_expr) != "DAPI"]
}


cols_sp <- c(
  'Centroid_X_µm','Centroid_Y_µm', 
  'ROI_Area_µm2', 'ROI_Name',
  'ObjectID', 'cell_id'
)

meta_df <-
  all_dat_raw %>%
  dplyr::select(all_of(cols_sp)) %>%
  dplyr::rename(
    spatial_X = Centroid_X_µm,
    spatial_Y = Centroid_Y_µm,
    ROI_Area  = ROI_Area_µm2
  )

dat_meas_expr2 <- t(dat_meas_expr)

stopifnot(all(colnames(dat_meas_expr2) == rownames(meta_df)))
sc <- CreateSeuratObject(
  CreateAssay5Object(
    counts = dat_meas_expr2[, !colnames(dat_meas_expr2) %in% remove_dapi_cells],
    data   = log1p(dat_meas_expr2)[, !colnames(dat_meas_expr2) %in% remove_dapi_cells],
  ),
  meta.data = meta_df[!rownames(meta_df) %in% remove_dapi_cells, ],
)

sc$Slide_ID <- 'Slide_1'
sc$ROI_ID <- sc$sample <- sc$project <- sc$ROI_Name


# raw data from HALO --------------
# DATA FROM HALO. USING CODE HERE.
# RUN code chunk:load_ini_data and annotation to acquire cd_final and meta_roi.
raw_dat <- GetAssayData(cd_final, layer = "counts")
stopifnot(all(colnames(raw_dat) == meta_roi$cell_id))


# # meta_roi$Slide_ID
# normalize_by_slide <- function(x, group_var, norm_func = c("q3", "median")){
#   norm_func <- match.arg(norm_func)
#   if(norm_func == "q3"){
#     norm_func  <- function(x) quantile(x, 0.75)
#   }else if(norm_func == "median"){
#     norm_func  <- median
#   }
#   
#   all_qc <- norm_func(x)
#   # group_qc <- tapply(x, group_var, norm_func)
#   
#   # norm_factor <- group_qc / all_qc
#   
#   xx <- vector(mode = class(x), length = length(x))
#   for(i in unique(group_var)){
#     idx <- x == i
#     temp_x <- x[idx]
#     xx[idx] <- temp_x / norm_func(temp_x) * all_qc
#   }
#   xx
# }
# 
# raw_dat[1,] %>% normalize_by_slide(group_var = meta_roi$Slide_ID)  

# combat
raw_dat_batch <- sva::ComBat(raw_dat, batch = cd_final$Slide_ID)
# raw_dat_batch2 <- limma::removeBatchEffect(raw_dat, batch = cd_final$Slide_ID)
sc <- CreateSeuratObject(
  CreateAssay5Object(
    raw_dat_batch %>% .[! rownames(.) %in% "DAPI", ] 
  ),
  meta.data = cd_final[[]]
)


# sc <-
#   cd_final %>%
#   .[! rownames(.) %in% "DAPI", ] %>% # remove DAPI
#   celldive2Seurat()


sc$sample <- sc$project <- sc$ROI_ID


# QC    ------
sc$log10GenesPerUMI <- log10(sc$nFeature_RNA + 1 + 1e5) / log10(sc$nCount_RNA + 1 + 1e5)

pdf(
  file.path(out_dir_cluster, "QC.pdf"),
  height = 6,
  width  = 8
)
VlnPlot(sc, features = "nCount_RNA", pt.size = 0, group.by = 'Slide_ID')
VlnPlot(sc, features = "nCount_RNA", pt.size = 0, log = TRUE, group.by = 'Slide_ID')
VlnPlot(sc, features = "log10GenesPerUMI", pt.size = 0, group.by = 'Slide_ID')
dev.off()



pdf(
  file.path(out_dir_cluster, "QC-feature.pdf"),
  height = 6,
  width  = 8
)
VlnPlot(
  sc, 
  features = rownames(sc), 
  pt.size = 0,
  group.by = 'Slide_ID',
  combine = FALSE
)  %>%
  walk(print)
dev.off()

pdf(
  file.path(out_dir_cluster, "QC-feature-log.pdf"),
  height = 6,
  width  = 8
)
VlnPlot(
  sc, 
  features = rownames(sc), 
  pt.size = 0,
  group.by = 'Slide_ID',
  log = TRUE,
  combine = FALSE
)  %>%
  walk(print)
dev.off()



# ncount_quantile_rng <- quantile(sc$nCount_RNA, c(0.001, 0.999))
# log10GeneperUMI_quantile_rng <- quantile(sc$log10GenesPerUMI, c(0.001, 0.999))
ncount_quantile_rng <- c(100, 1e6)
log10GeneperUMI_quantile_rng <- c(0.9, 0.999)


VlnPlot(sc, features = "log10GenesPerUMI", pt.size = 0, group.by = 'Slide_ID') +
  geom_hline(yintercept = log10GeneperUMI_quantile_rng)
VlnPlot(sc, features = "nCount_RNA", pt.size = 0, log = TRUE, group.by = 'Slide_ID')+
  geom_hline(yintercept = ncount_quantile_rng)


sc <- sc[ ,
          sc$nCount_RNA > ncount_quantile_rng[1] &
            sc$nCount_RNA < ncount_quantile_rng[2] &
            sc$log10GenesPerUMI > log10GeneperUMI_quantile_rng[1] &
            sc$log10GenesPerUMI < log10GeneperUMI_quantile_rng[2]
]


# hvg    ------
marker_celltype <- c(
  "CD11C",  "CD14",  "CD20",  "CD206",  "CD31",  "CD3e",  "CD4", "CD45",  "CD45RO",  "CD56",  "CD68",  "CD8",  
  "E-cadherin",  "FOXP3",  "GranzymeB",  "HLA-DR",  "Ki67",  "LAG3",  "Pan-cytokeratin",  "SMA",  "Vimentin"
  # "CTLA-4" 
  # "CD44",  
  # "EGFR"           
  # "HIF1a"
  # "Midkine"
  # "PD-1"
  # "PD-L1"
  # "TIM3"
  # "VISTA" 
)

sc <- sc %>% FindVariableFeatures(nfeatures = round(nrow(sc) * 0.75), selection.method = "vst")

p_hvg <- VariableFeaturePlot(sc)
p_hvg$data$gene <- rownames(p_hvg$data)
p_hvg2 <- p_hvg + ggrepel::geom_text_repel(aes(label = gene))
p_hvg2

# only used celltype marker
VariableFeatures(sc) <- VariableFeatures(sc) %>% intersect(marker_celltype)


pdf(
  file.path(out_dir_cluster, "hvg.pdf"),
  height = 6,
  width  = 8
)
p_hvg2 %>% print()
dev.off()


# norm   ------
normalize_q3 <- function(dat, log_transform = TRUE, min_val = 1e-5){
  stopifnot(length(dim(dat)) == 2)
  # stopifnot(is.numeric(dat))
  
  # dat <- dat + min_val
  
  q3 <- apply(dat, 2, quantile, 0.75)
  mean_q3  <- mean(q3)
  
  scale_fators <- mean_q3/q3
  dat2 <- sweep(dat, 2, scale_fators, '*')
  
  if(log_transform) dat2 <- log1p(dat2)
  dat2
}


# norm_assay <-
#   SetAssayData(
#     sc@assays$RNA,
#     slot = 'data',
#     new.data = normalize_q3(LayerData(sc, "counts"))
#   )
# 
# pdf(
#   file.path(out_dir_cluster, "normalize.pdf"),
#   height = 6,
#   width  = 8
# )
# 
# boxplot(as.matrix(norm_assay$counts[,1:50]))
# boxplot(as.matrix(norm_assay$data[,1:50]))
# boxplot(log1p(as.matrix(norm_assay$counts[,1:50])))
# dev.off()
# 
# # norm_assay$data[, "CWMDA_PCF_Slide_1...TS1_X2_Y1...TotalSegment...6233"]
# 
# sc[['RNA']] <- norm_assay
# # sc[['RNA']]@layers$counts <- log1p(sc[['RNA']]@layers$counts)


# # 20251107 only using typing marker
# typing_feature <- c(
#   "CD31",
#   "CD3e", "CD4", "CD8", "FoxP3",
#   "Ly6G", "S100A9",
#   "CD20", "CD45R/B220",
#   "F4/80", "CD206", "CD68", "CD11c", "Iba1",
#   "COL1A1", "Vimentin",
#   "PanCK",
#   "CD45", "Ki67", "Ter119"
# )
# sc[['type']] <- sc[typing_feature, ][['RNA']]
#
# DefaultAssay(sc) <- "type"


# # scaled by DAPI
# sc$DAPI <- FetchData(sc, vars = "DAPI")
# VariableFeatures(sc) <- VariableFeatures(sc) %>% setdiff("DAPI")
#
# sc <- sc %>% ScaleData( vars.to.regress = "DAPI" )



# pca ------
# sc <- sc %>% NormalizeData(normalization.method = "CLR", margin = 2)
sc <- sc %>% NormalizeData(normalization.method = NULL)

# na_idx <- sc@assays$RNA$data %>% is.na %>% apply(2, any)
# table(na_idx)
# sc <- sc[, !na_idx]

sc <- sc %>% ScaleData()
sc <- sc %>% RunPCA()
ElbowPlot(sc, ndims = 30)



pdf(
  file.path(out_dir_cluster, "pca.pdf"),
  height = 6,
  width  = 8
)
ElbowPlot(sc, ndims = 30)
dev.off()


pdf(
  file.path(out_dir_cluster, "DimHeatmap.pdf"),
  height = 6,
  width  = 8
)
DimHeatmap(sc, dims = 1:5, cells = 1000)
DimHeatmap(sc, dims = 6:10, cells = 1000)
DimHeatmap(sc, dims = 11:15, cells = 1000)
dev.off()



ndims <- 10


sc <- RunUMAP(sc, dims = 1:ndims, n.neighbors = 20, min.dist = 0.1, spread = 3)

# harmony    ------
sc <- sc %>%
  harmony::RunHarmony(
    dims.use = 1:ndims,
    group.by.vars = 'ROI_ID', # 'ROI_ID',
    reduction.save = "harmony",
    assay.use = "RNA"
  )

sc <- sc %>%
  RunUMAP(
    dims = 1:ndims,
    min.dist = 0.1,
    spread   = 3,
    n.neighbors = 20,
    # local.connectivity = 3,
    reduction = "harmony",
    reduction.name = "harmonyumap",
    reduction.key  = "harmonyUMAP_"
  )


pdf(
  file.path(out_dir_cluster, str_glue("umap - sample - dim{ndims}.pdf")),
  height = 6,
  width  = 8
)
DimPlot(sc, reduction = "umap", group.by = 'sample') # c("Group"))
DimPlot(sc, reduction = "harmonyumap", group.by = 'sample') # c("Group") )
dev.off()



# group
DimPlot(
  sc,
  reduction = "umap",
  group.by  = c("Group")
)

# group
DimPlot(
  sc,
  reduction = "harmonyumap",
  group.by  = c("Group")
)



FeaturePlot(
  sc,
  features = "CD3e",
  reduction = "harmonyumap"
) +
  scale_color_gradientn(
    colours = c("white", "white", RColorBrewer::brewer.pal(8, "Blues"))
  )

# cluster  ------------
sc <-
  FindNeighbors(
    sc,
    reduction = 'harmony', # 'harmony',
    dims = 1:ndims,
    k.param = 20, # 25, 20
    prune.SNN = 1/15, #1/15,
    # annoy.metric = "cosine"
  )

sc <-
  FindClusters(
    sc,
    resolution = c(0.2, 0.3, 0.4, 0.5, 0.6)
  )
sc <-
  FindClusters(
    sc,
    resolution = c(0.7, 0.8, 1.0)
  )

pdf(
  file.path(out_dir_cluster,  str_glue("umap - cluster - dim{ndims}.pdf")),
  height = 12,
  width  = 20
)
DimPlot(
  sc,
  group.by = c(str_glue('RNA_snn_res.0.{c(2:8)}'), 'RNA_snn_res.1'),
  reduction = 'umap', # "harmonyumap",
  label = TRUE,
  raster = TRUE,
  # raster.dpi = 100
) &
  NoLegend()
dev.off()


pdf(
  file.path(out_dir_cluster,  str_glue("umap - cluster - dim{ndims} - harmony.pdf")),
  height = 12,
  width  = 20
)
DimPlot(
  sc,
  group.by = c(str_glue('RNA_snn_res.0.{c(2:8)}'), 'RNA_snn_res.1'),
  reduction = 'harmonyumap', # "harmonyumap",
  label = TRUE,
  raster = TRUE,
  # raster.dpi = 100
) &
  NoLegend()
dev.off()


pdf(
  file.path(out_dir_cluster,  str_glue("umap - cluster - dim{ndims} - harmony - split.pdf")),
  height = 14,
  width  = 18
)
DimPlot(
  sc,
  group.by = 'RNA_snn_res.0.3',
  reduction = 'harmonyumap',  
  label = TRUE,
  raster = TRUE,
  ncol = 4,
  split.by = 'RNA_snn_res.0.3'
)
dev.off()

# sc$RNA_snn_res.0.3 %>% table()  %>%  .[str_sort(names(.), numeric = TRUE)]

# cell_gt_1000 <- sc$RNA_snn_res.0.3 %>% table() %>% .[. > 1000] %>% names
# sc <- sc[, sc$RNA_snn_res.0.3 %in% cell_gt_1000]
# sc$RNA_snn_res.0.3 <- sc$RNA_snn_res.0.3 %>% factor(., levels=str_sort(unique(.), numeric = T))

cl_lt_1000 <- sc$RNA_snn_res.0.3 %>% table %>% {names(.)[. < 1000]}
sc2 <- sc[, !sc$RNA_snn_res.0.3 %in% cl_lt_1000]


# plot cluster
pdf(
  file.path(out_dir_cluster,  str_glue("umap - cluster - dim{ndims} - v2.pdf")),
  height = 15,
  width  = 20
)
DimPlot(
  sc2,
  group.by = c(str_glue('RNA_snn_res.0.{c(2:8)}'), 'RNA_snn_res.1'),
  reduction = 'umap', # "harmonyumap",
  label = TRUE,
  raster = TRUE,
  # raster.dpi = 100
) &
  NoLegend()
dev.off()


pdf(
  file.path(out_dir_cluster,  str_glue("umap - cluster - dim{ndims} - harmony - v2.pdf")),
  height = 15,
  width  = 20
)
DimPlot(
  sc2,
  group.by = c(str_glue('RNA_snn_res.0.{c(2:8)}'), 'RNA_snn_res.1'),
  reduction = 'harmonyumap', # "harmonyumap",
  label = TRUE,
  raster = TRUE,
  # raster.dpi = 100
) &
  NoLegend()
dev.off()


# find markers -------
sc2  <- sc

select_res <- "RNA_snn_res.0.3"

Idents(sc2) <- select_res
# DefaultAssay(sc2) <- "RNA2"
all_markers <-
  FindAllMarkers(
    sc2,
    test.use = "wilcox",
    logfc.threshold = 0,
    return.thresh = 1.1,
    min.pct = 0.1,
    only.pos = TRUE
  )

# FindMarkers(
#   sc2,
#   ident.1 = '0',
#   ident.2 = NULL,
#   logfc.threshold = 0
# )


# # export marker gene
# all_markers %>%
#   dplyr::select(gene, everything()) %>%
#   write.csv(
#     file = file.path(mk_output, "marker.csv"),
#     row.names = FALSE
#   )

all_de_genes <-
  all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  dplyr::filter(abs(avg_log2FC) > 1) %>% ##
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  pull("gene")


all_top5_gene <-
  all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  # dplyr::filter(abs(avg_log2FC) > 1) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 5) %>%
  with(structure(
    gene,
    names = as.character(cluster)
  ))


# top5 dotplot
p_dot_top5 <-
  DotPlot(
    sc2,
    features = all_top5_gene %>% {structure(unique(.), names = names(.)[!duplicated(.)])}
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colours = c('white', 'white', 'blue', 'blue'))

ggsave(
  p_dot_top5,
  filename = file.path(out_dir_cluster, str_glue("top 5 gene dotplot - {select_res} - dim{ndims}.pdf")),
  height = 8,
  width  = 18
)


# DimPlot(
#   sc2[, !sc2$RNA_snn_res.0.3 %in% 12:14],
#   group.by = "RNA_snn_res.0.3",
#   reduction = "harmonyumap",
#   label = TRUE
# ) +
#   NoLegend()


# 上皮细胞（"Pan-Cytokeratin"  ）
# T cell（ "CD3e" "CD4" "CD8" "CD56" "Foxp3" "CD45RO" "GranzymeB"）
# B cell ( "CD20" "CD38")
# macrophage ( "CD68" "CD163")'
# neutrophil ( "MPO")
# monocyte ( "CD14" )
# mast cell ( "CD117" )
# DC ( "CD11C" )
# stroma ( "CD31" "SMA" "FAP" "CollagenIV" )
# FDC ( "CD21" )
# functional pr ( "Ki67" "PD-1" "PD-L1" "Histone H3-Lsy9" "Histone H3-Lys4"
#   "HLA-DR" "H2AX-Ser139" "HIF1a"
# )
# "IDO1" "LYVE-1" "SPP1"  


sc2@meta.data[[select_res]] %>% table
not_include_cluster <- c()

# not_include_cluster <- sc@meta.data[[select_res]] %>% table %>% {names(.)[. < 1000]} # less 1000


plot_data_dot <-
  sc2 %>%
  FetchData(
    vars = c(
      "RNA_snn_res.0.2",
      "CD11C",  "CD14",  "CD68",  "CD206", 
      "CD20",  
      "CD3e",  "CD4", "CD8", "FOXP3", "GranzymeB",  "CD56",  
      "CD45",  "CD45RO", 
      
      "E-cadherin",   "Pan-cytokeratin",
      "CD31", 
      
      "HLA-DR", 
      "Ki67", 
      "LAG3",
      "SMA",
      "Vimentin",
      
      "CTLA-4",
      "CD44",
      "EGFR",
      "HIF1a",
      "Midkine",
      "PD-1",
      "PD-L1",
      "TIM3",
      "VISTA"
    )
  )

plot_data_dot2 <-
  plot_data_dot %>%
  group_by(RNA_snn_res.0.2) %>% 
  summarise(across(everything(), mean))


p_dot_manual <-
  plot_data_dot2 %>%
  # dplyr::filter(!RNA_snn_res.0.2 %in% c(10:12)) %>%
  mutate(across(-RNA_snn_res.0.2, scale)) %>%
  pivot_longer(-RNA_snn_res.0.2, names_to = "gene", values_to = "value") %>%
  
  ggplot(aes(x = gene, y = RNA_snn_res.0.2, color = (value), alpha = (value))) +
  geom_point(size = 3, show.legend = c(color = T, alpha = F)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colours = c('white', 'white', 'blue'))  

ggsave(
  p_dot_manual,
  filename = file.path(out_dir_cluster, str_glue("dotplot - {select_res} - dim{ndims}.pdf")),
  height = 6,
  width  = 15
)




pdf(
  file.path(out_dir_cluster, str_glue("dotplot - {select_res} - dim{ndims}.pdf")),
  height = 8,
  width  = 18
)

DotPlot(
  sc2[, !sc2@meta.data[[select_res]] %in% not_include_cluster],
  group.by = select_res,
  features = c(
    
    "CD11C",  "CD14",  "CD68",  "CD206", 
    "CD20",  
    "CD3e",  "CD4", "CD8", "FOXP3", "GranzymeB",  "CD56",  
    "CD45",  "CD45RO", 
    
    "E-cadherin",   "Pan-cytokeratin",
    "CD31", 
    
    "HLA-DR", 
    "Ki67", 
    "LAG3",
    "SMA",
    "Vimentin",
    
    "CTLA-4",
    "CD44",
    "EGFR",
    "HIF1a",
    "Midkine",
    "PD-1",
    "PD-L1",
    "TIM3",
    "VISTA"
  )
) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colours = c('white', 'white', 'blue')) +
  ggtitle(select_res)


DotPlot(
  sc2[, !sc2@meta.data[[select_res]] %in% not_include_cluster],
  group.by = select_res,
  features = c(
    "CD11C",  "CD14",  "CD68",  "CD206", 
    "CD20",  
    "CD3e",  "CD4", "CD8", "FOXP3", "GranzymeB",  "CD56",  
    "CD45",  "CD45RO", 
    
    "E-cadherin",   "Pan-cytokeratin",
    "CD31", 
    
    "HLA-DR", 
    "Ki67", 
    "LAG3",
    "SMA",
    "Vimentin"
  )
) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colours = c(rep('white', 5), RColorBrewer::brewer.pal(9, "Blues"))) +
  ggtitle(select_res)
dev.off()


# # new assay
# sc2[['RNA2']] <- CreateAssayObject(
#   counts = sc2@assays$RNA$counts,
# )
# sc2[['RNA2']]$data <- log1p(sc2[['RNA2']]$counts)
# 
# 
# pdf(
#   file.path(out_dir_cluster, str_glue("dotplot - {select_res} - dim{ndims} - only log.pdf")),
#   height = 8,
#   width  = 18
# )
# DotPlot(
#   sc2[, !sc2@meta.data[[select_res]] %in% not_include_cluster],
#   group.by = select_res,
#   assay = "RNA2",
#   features = c(
#     "CD11C",  "CD14",  "CD68",  "CD206", 
#     "CD20",  
#     "CD3e",  "CD4", "CD8", "FOXP3", "GranzymeB",  "CD56",  
#     "CD45",  "CD45RO", 
#     
#     "E-cadherin",   "Pan-cytokeratin",
#     "CD31", 
#     
#     "HLA-DR", 
#     "Ki67", 
#     "LAG3",
#     "SMA",
#     "Vimentin",
#     
#     "CTLA-4",
#     "CD44",
#     "EGFR",
#     "HIF1a",
#     "Midkine",
#     "PD-1",
#     "PD-L1",
#     "TIM3",
#     "VISTA"
#   )
# ) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#   scale_color_gradientn(colours = c('white', 'white', 'blue')) +
#   ggtitle(select_res)
# 
# 
# DotPlot(
#   sc2[, !sc2@meta.data[[select_res]] %in% not_include_cluster],
#   group.by = select_res,
#   assay = "RNA2",
#   features = c(
#     "CD11C",  "CD14",  "CD68",  "CD206", 
#     "CD20",  
#     "CD3e",  "CD4", "CD8", "FOXP3", "GranzymeB",  "CD56",  
#     "CD45",  "CD45RO", 
#      
#     "E-cadherin",   "Pan-cytokeratin",
#     "CD31", 
# 
#     "HLA-DR", 
#     "Ki67", 
#     "LAG3",
#     "SMA",
#     "Vimentin"
#   )
# ) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#   scale_color_gradientn(colours = c(rep('white', 5), RColorBrewer::brewer.pal(9, "Blues"))) +
#   ggtitle(select_res)
# dev.off()




# DoHeatmap(
#   sc,
#   features = c(
#     "CD14" , "CD20"       , "CD31"    , "CD34" ,  "CD3e"        , "CD4"             ,
#     "CD45" , "CD45RO"     , "CD68"    , "CD8"  ,  "Collagen IV" , "E-cadherin"  ,
#     "FOXP3", "Ki67"       , "Pan-Cytokeratin" ,
#     "SMA"        , "Vimentin"
#   ),
#   group.by = "RNA_snn_res.0.3"
# )

# annotation --------
# "High endothelial venule endothelial cell"
anno_cluster <- c(
  '0'  = "Fibroblast",
  '1'  = "Myeloid cell",
  '2'  = "ECAD Epithelial",
  '3'  = "Memory CD8T",
  '4'  = "CD8T",
  '5'  = "Endothelial",
  '6'  = "PANCK Epithelial",
  '7'  = "Proliferating Epithelial",
  '8'  = "CD4T",
  '9'  = "B",
  '10' = "Others",
  '11' = "Others",
  '12' = "Others"
)


Idents(sc) <- select_res
sc <- sc %>% RenameIdents(anno_cluster)
sc$annotation <- Idents(sc)


# export umap of annotation
pdf(
  file.path(out_dir_cluster,  str_glue("umap - annotation - dim{ndims} - {select_res} - v2.pdf")),
  height = 6,
  width  = 8
)
DimPlot(
  sc,
  group.by = "annotation",
  reduction = "harmonyumap", # "harmonyumap",
  label = TRUE,
  raster = TRUE,
  # raster.dpi = 100
) +
  coord_cartesian(clip = "off") +
  NoLegend()
dev.off()

# export umap of group
pdf(
  file.path(out_dir_cluster,  str_glue("umap - group - dim{ndims} - v2.pdf")),
  height = 6,
  width  = 8
)
DimPlot(
  sc,
  group.by = "Group",
  reduction = "harmonyumap",
  label = TRUE,
  raster = TRUE,
  # raster.dpi = 100
) +
  coord_cartesian(clip = "off") +
  NoLegend()
dev.off()


# dotplot
pdf(
  file.path(out_dir_cluster, "dotplot - annotation.pdf"),
  height = 6,
  width  = 18
)
DotPlot(
  sc[, sc$annotation != "Others"],
  group.by = 'annotation',
  features = c(
    "CD11C",  "CD14",  "CD68",  "CD206", 
    "CD20",  
    "CD3e",  "CD4", "CD8", "FOXP3", "GranzymeB",  "CD56",  
    "CD45",  "CD45RO", 
    
    "E-cadherin",   "Pan-cytokeratin",
    "CD31", 
    
    "HLA-DR", 
    "Ki67", 
    "LAG3",
    "SMA",
    "Vimentin"
  )
) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colours = c('white', 'white', 'blue'))

dev.off()

# saveRDS(sc, file.path(out_dir_cluster, "sc.rds"))


# subcluster ------------
Idents(sc) <- select_res
sc <- sc %>%
  FindSubCluster(
    cluster = '4',
    graph.name = "RNA_snn",
    subcluster.name = "sub.cluster_4",
    resolution = 0.1,
    algorithm = 1
  )

Idents(sc) <- 'sub.cluster_4'
sc <- sc %>%
  FindSubCluster(
    cluster = '13',
    graph.name = "RNA_snn",
    subcluster.name = "sub.cluster_4_13",
    resolution = 0.1,
    algorithm = 1
  )

# Idents(sc) <- 'sub.cluster_4_11'
# sc <- sc %>%
#   FindSubCluster(
#     cluster = '1',
#     graph.name = "RNA_snn",
#     subcluster.name = "sub.cluster_1_4_11",
#     resolution = 0.2,
#     algorithm = 1
#   )

select_res_sub <- 'sub.cluster_4_13'

# subcluster umap 
pdf(
  file.path(out_dir_cluster,  str_glue("umap - {select_res_sub} - dim{ndims}.pdf")),
  height = 15,
  width  = 20
)
DimPlot(
  sc,
  group.by = select_res_sub,
  reduction = 'umap', # "harmonyumap",
  label = TRUE,
  raster = TRUE,
  # raster.dpi = 100
) &
  NoLegend()
dev.off()


pdf(
  file.path(out_dir_cluster,  str_glue("umap - {select_res_sub} - dim{ndims} - harmony.pdf")),
  height = 15,
  width  = 20
)
DimPlot(
  sc,
  group.by = select_res_sub,
  reduction = 'harmonyumap', # "harmonyumap",
  label = TRUE,
  raster = TRUE,
  # raster.dpi = 100
) &
  NoLegend()
dev.off()


# 无subcluster的必要

# subcluster dotplot 
pdf(
  file.path(out_dir_cluster, str_glue("dotplot - {select_res_sub} - dim{ndims}.pdf")),
  height = 8,
  width  = 18
)
DotPlot(
  sc[, !sc@meta.data[[select_res]] %in% not_include_cluster],
  group.by = select_res_sub,
  features = c(
    "CD31", "LYVE-1",
    "CD3e", "CD4", "CD8", "Foxp3", "CD56", "CD45RO", "GranzymeB",
    "CD20", "CD38",
    "CD14", "CD68", "CD163",
    "CD117", # Mast
    "CD11C",  # DC
    "CD21", # FDC
    "MPO", # "Neu"
    "Pan-Cytokeratin", # epi
    "FAP", "SMA", "CollagenIV", 
    
    "Histone H3-Lsy9", "Histone H3-Lys4", "H2AX-Ser139",
    "Ki67", 
    "HIF1a", 
    "HLA-DR", "IDO1", "PD-1", "PD-L1", "SPP1" 
  )
) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colours = c('white', 'white', 'blue')) +
  ggtitle(select_res)


DotPlot(
  sc[, !sc@meta.data[[select_res]] %in% not_include_cluster],
  group.by = select_res_sub,
  features = c(
    "CD31", "LYVE-1",
    "CD3e", "CD4", "CD8", "Foxp3", "CD56", "CD45RO", "GranzymeB",
    "CD20", "CD38",
    "CD14", "CD68", "CD163",
    "CD117", # Mast
    "CD11C", "HLA-DR", # DC
    "CD21", # FDC
    "MPO", # "Neu"
    "Pan-Cytokeratin", # epi
    "FAP", "SMA", "CollagenIV"
  )
) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colours = c(rep('white', 5), RColorBrewer::brewer.pal(9, "Blues"))) +
  ggtitle(select_res)
dev.off()

# anno_subcluster 
anno_subcluster <- c(
  '0'  = "Others",
  '1_0'  = "T_B", # 4-8
  '1_1'  = "CD8T", # 4-8
  '1_2'  = "fDC", # 4-8
  '2'  = "Epithelial",
  '3'  = "CD4T",
  '4_0'  = "DC", # CD11c
  '4_1'  = "High endothelial venule endothelial cell",  
  '4_2'  = "High endothelial venule endothelial cell", 
  '4_3'  = "Macrophage", # macro  
  '4_4'  = "Macrophage", # macro 
  '5'  = "B",
  '6'  = "Endothelial",
  '7'  = "CD4T",
  '8'  = "Epithelial",
  '9'  = "Plasma",
  '10' = "CD8T",
  '11_0' = "Treg",  
  '11_1' = "DC", # CD11c
  '12' = "Smooth muscle cell",
  '13' = "fDC",
  '14' = "Neutrophil",
  '15' = "CD8T"
)



Idents(sc) <- select_res_sub
sc <- sc %>% RenameIdents(anno_subcluster)
sc$annotation_sub <- Idents(sc)


# export umap of annotation
pdf(
  file.path(out_dir_cluster,  str_glue("umap - annotation - dim{ndims} - {select_res_sub} - v2.pdf")),
  height = 6,
  width  = 8
)
DimPlot(
  sc,
  group.by = "annotation_sub",
  reduction = "harmonyumap", # "harmonyumap",
  label = TRUE,
  raster = TRUE,
  # raster.dpi = 100
) +
  coord_cartesian(clip = "off") +
  NoLegend()
dev.off()




# dotplot
pdf(
  file.path(out_dir_cluster, "dotplot - annotation_sub.pdf"),
  height = 6,
  width  = 18
)
DotPlot(
  sc,
  group.by = 'annotation_sub',
  features = c(
    "CD31", "LYVE-1",
    "CD3e", "CD4", "CD8", "Foxp3", "CD56", "CD45RO", "GranzymeB",
    "CD20", "CD38",
    "CD14", "CD68", "CD163",
    "CD117", # Mast
    "CD11C", # DC
    "CD21", # FDC
    "MPO", # "Neu"
    "Pan-Cytokeratin", # epi
    "FAP", "SMA", "CollagenIV"
  )
) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colours = c('white', 'white', 'blue'))

dev.off()

# saveRDS(sc, file.path(out_dir_cluster, "sc.rds"))


# plot ------------------
## spatial plot of cluster ---------------
p_sp_cluster <-
  sc[[]] %>% #  head(10000) %>%
  ggplot() +
  geom_rect(aes(xmin = XMin, xmax = XMax, ymin = YMin, ymax = YMax, fill = .data[[select_res]])) +
  # geom_tile(
  #   aes(x = spatial_X, y = spatial_Y, fill = .data[[select_res]]),
  #   height = 15, width = 15
  # ) +
  # geom_text(
  #   data = ~ summarise(
  #     group_by(., ROI_ID),
  #     x = mean(spatial_X),
  #     y = max( spatial_Y)
  #   ),
  #   aes(x = x, y = y, label = ROI_ID),
#   show.legend = FALSE,
#   vjust = 1
# ) +
# annotate(
#   "text",
#   label = roi,
#   x = mean((meta$XMin + meta$XMax) * 0.5),
#   y = max((meta$YMin + meta$YMax) * 0.5),
#   vjust = 2
# ) +
ggsci::scale_fill_igv() +
  scale_y_reverse() +
  coord_equal(clip = 'off') +
  theme_void()

ggsave(
  p_sp_cluster,
  filename = file.path(out_dir_cluster, str_glue("spatial plot of annotation - {select_res}.pdf")),
  height = 22,
  width  = 10,
  limitsize = FALSE
)




# fator
cell_order <- c(
  'CD4T', 'CD8T', "Memory CD8T",
  'ECAD Epithelial', "PANCK Epithelial", "Proliferating Epithelial",
  
  'Fibroblast', 
  'Endothelial', 
  
  'B',
  'Myeloid cell',
  'Others'
)

protein_order <- c(
  'T' = "CD3e",
  'T' = "CD4",
  'T' = "CD8",
  'T' = "GranzymeB",
  'T' = "CD45RO",
  Treg = "FOXP3",
  
  Epithelial = "Pan-cytokeratin",
  Epithelial = "E-cadherin",
  Proliferating = "Ki67",
  
  Fibroblast = "Vimentin",
  Fibroblast = "SMA",
  
  Endothelial = "CD31",
  
  B = "CD20",
  
  
  'Myeloid cell' = "CD11C",
  'Myeloid cell' = "CD14",
  'Myeloid cell' = "CD68",
  'Myeloid cell' = "CD206"
)

sc$annotation  <- sc$annotation %>% factor(levels = cell_order)


## plot umap annotation ------
p_anno <-
  DimPlot(
    sc,
    label = TRUE,
    reduction = "harmonyumap",
    group.by = "annotation"
  ) +
  NoLegend() +
  ggsci::scale_color_igv()

ggsave(
  p_anno,
  filename = file.path(out_dir_cluster, "umap - annotation.pdf"),
  width  = 8,
  height = 6
)

## plot dotplot 1 ------
p_dot_annotation <-
  DotPlot(
    sc,
    group.by = "annotation",
    features = protein_order,
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.clip = 'off'
  ) +
  # scale_color_gradientn(colours = c(rep('white', 5), RColorBrewer::brewer.pal(9, "Blues"))) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu"))) +
  coord_cartesian(clip = "off")

p_dot_annotation2 <-
  p_dot_annotation +
  scale_color_gradientn(colours = c(rep('white', 5), RColorBrewer::brewer.pal(9, "Blues")))



pdf(
  file.path(out_dir_cluster, "dot - annotation.pdf"),
  height = 6,
  width  = 20
)
p_dot_annotation2 %>% print()
p_dot_annotation %>% print()
dev.off()


## plot dotplot 2 ------
p_dot_annotation3 <-
  DotPlot(
    sc,
    group.by = "annotation",
    features = protein_order %>% c(., setdiff(rownames(sc), c(., "DAPI")) ),
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.clip = 'off'
  ) +
  # scale_color_gradientn(colours = c(rep('white', 5), RColorBrewer::brewer.pal(9, "Blues"))) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu"))) +
  coord_cartesian(clip = "off")


p_dot_annotation4 <-
  p_dot_annotation3 +
  scale_color_gradientn(colours = c(rep('white', 5), RColorBrewer::brewer.pal(9, "Blues")))


pdf(
  file.path(out_dir_cluster, "dot - annotation - all markers.pdf"),
  height = 6,
  width  = 25
)
p_dot_annotation4 %>% print()
p_dot_annotation3 %>% print()
dev.off()


## plot avg_expr --------
avg_expr_rna <-
  sc %>%
  LayerData("data") %>%
  as.matrix() %>%
  apply(1, function(x) tapply(x, sc$annotation, function(y) log1p(mean(expm1(y))) ) )

pdf(
  file.path(out_dir_cluster, "heatmap - annotation.pdf"),
  height = 6,
  width  = 10
)
pheatmap::pheatmap(
  avg_expr_rna[cell_order, protein_order],
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(50),
  scale = "column",
  border_color = NA,
  cluster_rows = FALSE,
  cluster_cols = FALSE
)
dev.off()

# counts
avg_expr_rna_counts <-
  sc %>%
  LayerData("counts") %>%
  as.matrix() %>%
  # apply(1, function(x) tapply(x, sc$annotation, function(y) log1p(mean(expm1(y))) ) )
  apply(1, function(x) tapply(x, sc$annotation, function(y) log1p(mean(y)) ) )

pdf(
  file.path(out_dir_cluster, "heatmap - annotation - raw.pdf"),
  height = 6,
  width  = 10
)
pheatmap::pheatmap(
  avg_expr_rna_counts[cell_order, protein_order],
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(50),
  scale = "column",
  border_color = NA,
  cluster_rows = FALSE,
  cluster_cols = FALSE
)
dev.off()

## spatial plot of cluster ---------------
p_sp_ls <-
  sc$ROI_ID %>%
  str_c(sc$Slide_ID, sep = "-") %>%
  unique() %>%
  setNames(., .) %>%
  map(function(roi){ # browser()
    sc$roi_id2 <- str_c(sc$ROI_ID, sc$Slide_ID, sep = "-")
    meta <- sc[[]] %>% dplyr::filter(roi_id2 == roi)
    
    
    meta %>%
      ggplot() +
      geom_rect(aes(xmin = XMin, xmax = XMax, ymin = YMin, ymax = YMax, fill = annotation)) +
      # geom_tile(
      #   aes(x = spatial_X, y = spatial_Y, fill = annotation),
      #   height = 15, width = 15
      # ) +
      geom_text(
        data = ~ summarise(
          group_by(., ROI_ID),
          x = median( (XMin + XMax)/2 ), # mean(spatial_X),
          y = max(YMax), # max( spatial_Y)
        ),
        aes(x = x, y = y, label = ROI_ID),
        show.legend = FALSE,
        vjust = 1
      ) +
      ggtitle(roi) +
      ggsci::scale_fill_igv() +
      scale_y_reverse() +
      coord_equal(clip = 'off') +
      theme_void() +
      theme(plot.margin = margin(5, 5, 5, 5, 'mm'))
  })

pdf(
  file.path(out_dir_cluster, "spatial plot of annotation - split by sample.pdf"),
  height = 4,
  width  = 6
)
p_sp_ls %>% walk(print)
dev.off()


p_sp <-
  sc[[]] %>% #  head(10000) %>%
  ggplot() +
  geom_rect(aes(xmin = XMin, xmax = XMax, ymin = YMin, ymax = YMax, fill = annotation)) +
  geom_text(
    data = ~ summarise(
      group_by(., ROI_ID),
      x = mean((XMin + XMax) * 0.5),
      y = max((YMin + YMax) * 0.5)
    ),
    aes(x = x, y = y, label = ROI_ID),
    show.legend = FALSE,
    vjust = 1
  ) +
  # annotate(
  #   "text",
  #   label = roi,
  #   x = mean((meta$XMin + meta$XMax) * 0.5),
  #   y = max((meta$YMin + meta$YMax) * 0.5),
  #   vjust = 2
  # ) +
  ggsci::scale_fill_igv() +
  scale_y_reverse() +
  coord_equal(clip = 'off') +
  theme_void()

ggsave(
  p_sp,
  filename = file.path(out_dir_cluster, "spatial plot of annotation.pdf"),
  height = 22,
  width  = 10,
  limitsize = FALSE
)


# # export
# saveRDS(sc, str_glue("sc-dim{ndims}.rds"))

## spatial plot of intensity  ---------------
plot_dat <- 
  sc %>%
  FetchData(
    vars = c(rownames(sc), 'XMin', 'XMax', 'YMin', 'YMax', 'ROI_ID'),
    layer = "counts"
  ) %>%
  mutate(
    spatial_X = (XMin + XMin)/2,
    spatial_Y = (YMin + YMax)/2
  )
# FetchData(
#   vars = c(rownames(sc), 'spatial_X', 'spatial_Y', 'ROI_ID'),
#   layer = "counts"
# )

p_sp_intensity_ls <-
  rownames(sc) %>%
  setNames(., .) %>%
  map(function(x){
    plot_dat %>% #  head(10000) %>%
      ggplot() +
      # geom_rect(aes(xmin = XMin, xmax = XMax, ymin = YMin, ymax = YMax, fill = .data[[select_res]])) +
      geom_tile(
        aes(x = spatial_X, y = spatial_Y, fill = .data[[x]]),
        height = 15, width = 15
      ) +
      geom_text(
        data = ~ summarise(
          group_by(., ROI_ID),
          x = mean(spatial_X),
          y = max( spatial_Y)
        ),
        aes(x = x, y = y, label = ROI_ID),
        show.legend = FALSE,
        vjust = 1
      ) +
      scale_fill_gradientn(colours = Seurat:::SpatialColors(30), limits = c(0, 65536)) +
      scale_y_reverse() +
      coord_equal(clip = 'off') +
      theme_void() +
      ggtitle(x)
  })

intesity_dir <- file.path(out_dir_cluster, str_glue("spatial plot of intensity"))
dir.create(intesity_dir)


p_sp_intensity_ls %>%
  iwalk(function(p_ls, name){
    pdf(
      file.path(intesity_dir, str_glue("spatial plot of intensity - {name}.pdf")),
      height = 22,
      width  = 10,
    )
    p_ls %>% print()
    dev.off()
  })


# export cell porps table ----------------

sc[[]] %>%
  group_by(ID, Group, ROI_ID, annotation) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  dplyr::select(-count) %>%
  pivot_wider(
    names_from = "annotation",
    values_from = "freq"
  ) %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir_cluster, "cell props.xlsx"),
    asTable = TRUE,
    withFilter = FALSE
  )

sc[[]] %>%
  group_by(ID, Group, ROI_ID, annotation) %>%
  summarise(count = n()) %>%
  
  pivot_wider(
    names_from = "annotation",
    values_from = "count"
  ) %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir_cluster, "cell count.xlsx"),
    asTable = TRUE,
    withFilter = FALSE
  )



# 20260626 需求 -------
# 1 major annotation  
sc$annotation2 <-
  case_when(
    sc$annotation %in% c("ECAD Epithelial", "PANCK Epithelial", "Proliferating Epithelial", "Others") ~ "Epithelial",
    sc$annotation %in% c("CD4T", "CD8T", "Memory CD8T") ~ "T",
    .default = sc$annotation
  )

# export umap of group
pdf(
  file.path(out_dir_cluster,  str_glue("umap - annotation - major.pdf")),
  height = 6,
  width  = 8
)
DimPlot(
  sc,
  group.by = "annotation2",
  reduction = "harmonyumap",
  label = TRUE,
  raster = TRUE,
  # raster.dpi = 100
) +
  coord_cartesian(clip = "off") +
  NoLegend()
dev.off()


# dotplot
pdf(
  file.path(out_dir_cluster, "dotplot - annotation  - major.pdf"),
  height = 6,
  width  = 18
)
DotPlot(
  sc,
  group.by = 'annotation2',
  features = c(
    "CD11C",  "CD14",  "CD68",  "CD206", 
    "CD20",  
    "CD3e",  "CD4", "CD8", "FOXP3", "GranzymeB",  "CD56",  
    "CD45",  "CD45RO", 
    
    "E-cadherin",   "Pan-cytokeratin",
    "CD31", 
    
    "HLA-DR", 
    "Ki67", 
    "LAG3",
    "SMA",
    "Vimentin"
  )
) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu")))
dev.off()

# 2 cell proportion 
sc[['annotation2']] %>%
  group_by(annotation2) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  write.csv(
    row.names = FALSE,
    file = file.path(out_dir_cluster, 'cell count of each major annotation.csv')
  )

sc[['annotation']] %>%
  group_by(annotation) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  write.csv(
    row.names = FALSE,
    file = file.path(out_dir_cluster, 'cell count of each annotation.csv')
  )

# 3 cell befor and after QC 
data.frame(
  Before_QC = ncol(cd_final),
  After_QC = ncol(sc)
) %>%
  write.csv(
    row.names = FALSE,
    file = file.path(out_dir_cluster, 'cell count - before and after QC.csv')
  )



# test OSTU threshold -------------
otsu_manual <- function(x, levels = 256) {
  # 移除 NA、NaN、Inf 值
  x_clean <- x[is.finite(x)]
  if (length(x_clean) == 0) {
    warning("输入数据没有有效值")
    return(list(normalized = NA_real_, original = NA_real_))
  }
  
  # 保存原始数据的极值（用于反归一化）
  min_orig <- min(x_clean)
  max_orig <- max(x_clean)
  range_orig <- max_orig - min_orig
  
  # 如果所有值相同，阈值即为该值（归一化阈值固定为 0.5 或任意，但这里直接返回）
  if (range_orig == 0) {
    warning("所有数据值相同，阈值设置为该值")
    return(list(normalized = 0.5, original = min_orig))
  }
  
  # 将数据归一化到 [0,1] 区间
  x_norm <- (x_clean - min_orig) / range_orig
  x_norm <- pmax(0, pmin(1, x_norm))   # 防止微小越界
  
  # 计算直方图
  hist_info <- hist(x_norm, breaks = seq(0, 1, length.out = levels + 1), plot = FALSE)
  counts <- hist_info$counts
  
  N <- sum(counts)
  if (N == 0) {
    warning("直方图计数总和为零")
    return(list(normalized = NA_real_, original = NA_real_))
  }
  
  max_var <- 0
  best_threshold <- 1  # 默认值
  
  for (t in 1:levels) {
    w0 <- sum(counts[1:t], na.rm = TRUE) / N
    w1 <- sum(counts[(t+1):levels], na.rm = TRUE) / N
    
    if (w0 == 0 || w1 == 0) next
    
    mu0 <- sum((1:t) * counts[1:t], na.rm = TRUE) / sum(counts[1:t], na.rm = TRUE)
    mu1 <- sum(((t+1):levels) * counts[(t+1):levels], na.rm = TRUE) / sum(counts[(t+1):levels], na.rm = TRUE)
    
    between_var <- w0 * w1 * (mu0 - mu1)^2
    
    if (between_var > max_var) {
      max_var <- between_var
      best_threshold <- t
    }
  }
  
  # 归一化阈值 (在 [0,1] 范围内)
  thresh_norm <- best_threshold / levels
  
  # 反归一化得到原始尺度下的阈值
  thresh_orig <- min_orig + thresh_norm * range_orig
  
  # 返回两个值
  return(list(normalized = thresh_norm, original = thresh_orig))
}


expr_c1 <- sc@assays$RNA$counts[1, ]
otsu_manual(expr_c1, levels = 256)



expr_c1 %>% density() %>% plot
abline(v = 499)



expr_c2 <- sc@assays$RNA$counts['CD4', ]
otsu_manual(expr_c2, levels = 256)


expr_c2 %>% density() %>% plot
abline(v = 367)





