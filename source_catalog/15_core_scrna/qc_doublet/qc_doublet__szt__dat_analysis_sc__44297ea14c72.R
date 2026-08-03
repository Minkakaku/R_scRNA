library(Seurat)
library(tidyverse)

sc_out_dir <- "sc_output" 
dir.create(sc_out_dir)


col_plan_1 <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
                "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", "#de4307", 
                "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", "#fa6e57", "#f69e53", 
                "#248888", "#999999", "#E7475E", "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", 
                "#FF6D31", "#f1ac9d", "#f06966", "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", 
                "#FF7A5A", "#8ED2C9", "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", 
                "#ffc300", "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
                "#f17d80", "#737495")

col_plan_2 <- c("#1F77B4FF","#FF7F0EFF","#2CA02CFF","#D62728FF","#9467BDFF",
                "#8C564BFF","#E377C2FF","#7F7F7FFF","#BCBD22FF","#17BECFFF",
                "#AEC7E8FF","#FFBB78FF","#98DF8AFF","#FF9896FF","#C5B0D5FF",
                "#C49C94FF","#F7B6D2FF","#C7C7C7FF","#DBDB8DFF","#9EDAE5FF")

# file_dir <- "/data/outputdata/outputdata-sc/cellranger_out/FWSC20240319"
# all_files_dir <- list.files(file_dir, pattern = "^[XYZ]-.*$")
# 
# 
# all_files <- file.path(file_dir, all_files_dir, "step3/filtered_feature_bc_matrix")
# names(all_files) <- all_files_dir
# 
# sample_info <-
#   data.frame(
#     sample = all_files_dir,
#     group  = all_files_dir %>% str_extract("^[XYZ]"),
#     id = seq_along(all_files_dir)
#   )

# 
# # 6 samples
# stopifnot(length(all_files) == 6)

# all_data_ls <- 
#   all_files %>%
#   map(Read10X) %>%
#   imap(~CreateSeuratObject(
#     .x,
#     project = .y,
#     assay = "RNA",
#     names.field = 1,
#     names.delim = "_",
#     min.cells = 10,
#     min.features = 10
#   ))

# all_data_ls <-
#   all_data_ls %>%
#   map(function(x){
#     x$project <- x@project.name
#     x
#   })


meta_data <- read.csv("data/GSE151974_cell_metadata_postfilter.csv", row.names = 1) 
raw_count <- data.table::fread("data/GSE151974_raw_umi_matrix_postfilter.csv") %>% column_to_rownames("V1")


all_data_raw <- 
  CreateSeuratObject(
    raw_count,
    assay = "RNA",
    names.field = 1,
    names.delim = "_",
    min.cells = 10,
    min.features = 10
  )
Idents(all_data_raw) <- "s"
all_data_raw <- 
  all_data_raw %>%
  AddMetaData(meta_data[c('Barcode', 'Age', 'Oxygen', 'Sample', 'CellType')])


all_data_raw$project <- all_data_raw$Sample
# meta_data[c('Age', 'Oxygen')] %>% distinct()
# meta_data[c('Sample')] %>% distinct()


all_data_ls <- SplitObject(all_data_raw, split.by =  "project")

### QC ---------------------
all_data_ls <- 
  all_data_ls %>%
  map(function(x){
    x[["percent.mt"]]  <- PercentageFeatureSet(x, pattern = "^mt-")
    x[["percent.ribo"]]  <- PercentageFeatureSet(x, pattern = "^Rp[sl]")
    x$log10GenesPerUMI <- log10(x$nFeature_RNA) / log10(x$nCount_RNA)
    x
  })


p_qc_ls <-
  all_data_ls %>%
  imap(function(x, name){
    VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", 
                            "percent.mt", "percent.ribo", "log10GenesPerUMI"), 
            ncol = 5, log = FALSE) +
      patchwork::plot_annotation(title = name)
  })

p_qc_log_ls <-
  all_data_ls %>%
  imap(function(x, name){
    VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", 
                            "percent.mt", "percent.ribo", "log10GenesPerUMI"), 
            ncol = 5, log = TRUE) +
      patchwork::plot_annotation(title = str_glue('{name} + log'))
  })

pdf(
  file.path(sc_out_dir, "qc plot.pdf"),
  height = 6, 
  width = 20
)
p_qc_ls %>% walk(print)
p_qc_log_ls %>% walk(print)
dev.off()




### filter ---------------
qc_metric <- read.csv("qc_metric_sc.csv") %>% map_dfc(str_trim)
qc_metric_variable <- c(
  'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.ribo', 'log10GenesPerUMI'
)
qc_metric_ls <-
  qc_metric %>% 
  apply(1, name = names(.), simplify = FALSE, function(x, name = name){
    qc_ls <- map2(x, name, function(num, variable){
      if(!variable %in% qc_metric_variable) return()
      
      if(str_detect(num, ":")){
        num <- str_split(num, pattern = ":", n = 2) %>% .[[1]]
        qc <- str_glue("{variable} {num[1]} & {variable} {num[2]}")
      }else{
        qc <- str_glue("{variable} {num}")
      }
      qc
    })
    
    res <- qc_ls %>% discard(function(x) length(x) < 1) %>% str_c(collapse = " & ")
    res
  })

names(qc_metric_ls) <- qc_metric$sample

# filter
all_data_ls_pass <-
  names(all_data_ls) %>%
  set_names(.,.) %>%  
  map(function(sample){ # browser()
    all_data <- all_data_ls[[sample]]
    qc <- qc_metric_ls[[sample]]
    
    e <- str_glue("subset(
      all_data,
      {qc}
    )")
    
    eval(parse(text = e))
  })


# qc plot
qc_plot_filter_ls <-
  names(all_data_ls) %>%
  map(function(sample){
    all_data <- all_data_ls[[sample]]
    all_data_pass <- all_data_ls_pass[[sample]]
    
    # all_data
    plot1 <- FeatureScatter(all_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
      patchwork::plot_annotation(title = sample)
    plot2 <- FeatureScatter(all_data, feature1 = "nCount_RNA", feature2 = "percent.mt") +
      patchwork::plot_annotation(title = sample)
    plot3 <- FeatureScatter(all_data, feature1 = "nCount_RNA", feature2 = "log10GenesPerUMI") +
      patchwork::plot_annotation(title = sample)
    p1 <- plot1 + plot2 + plot3 + patchwork::plot_layout(ncol = 3)
    
    # alll_data_pass
    plot1_pass <- FeatureScatter(all_data_pass, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    plot2_pass <- FeatureScatter(all_data_pass, feature1 = "nCount_RNA", feature2 = "percent.mt")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    plot3_pass <- FeatureScatter(all_data_pass, feature1 = "nCount_RNA", feature2 = "log10GenesPerUMI")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    p2 <- plot1_pass + plot2_pass + plot3_pass  + patchwork::plot_layout(ncol = 3)
    
    p1 / p2
  })

# export 
pdf(
  file.path(sc_out_dir, "qc plot - pass.pdf"),
  height = 12, 
  width = 18
)
qc_plot_filter_ls %>% walk(print)
dev.off()


# cells filter
cell_counts_df <-
  names(all_data_ls) %>%
  set_names(., .) %>%
  map(function(sample){
    all_data <- all_data_ls[[sample]]
    all_data_pass <- all_data_ls_pass[[sample]]
    
    structure(
      c(ncol(all_data), ncol(all_data_pass)),
      names = c("Before filter", "After filter")
    )
  }) %>%
  as.data.frame(check.names = FALSE)

write.csv(cell_counts_df, file = file.path(sc_out_dir, "cell_counts_df.csv"))

# #### add cell_id prefix
# all_data_ls_pass <-
#   all_data_ls_pass %>%
#   names() %>%
#   set_names(., .) %>%
#   map(function(x){
#     all_data <- all_data_ls_pass[[x]]
#     ids <- sample_info$id[sample_info$sample == x]
#     
#     all_data <- RenameCells(all_data, add.cell.id = ids)
#     all_data
#   })

all_data <- merge(all_data_ls_pass[[1]], all_data_ls_pass[-1])


### normalize ------------
all_data <- NormalizeData(all_data)
all_data <- FindVariableFeatures(all_data)
all_data <- ScaleData(all_data, vars.to.regress = 'nCount_RNA')

# ### cellcycle score -------------
# all_data <- CellCycleScoring(
#   all_data,
#   s.features = cc.genes.updated.2019$s.genes,
#   g2m.features = cc.genes.updated.2019$g2m.genes
# )

# dscDblFinder ----------
# all_data_dbl <- 
#   all_data %>% 
#   JoinLayers() %>%
#   as.SingleCellExperiment() %>% 
#   scDblFinder::scDblFinder(
#     samples="project", 
#     dbr = 0.1,
#     BPPARAM=BiocParallel::MulticoreParam(5)
#   )
# 
# dbl_finder_df <- 
#   all_data_dbl@colData %>%
#   .[str_detect(colnames(.), "^scDblFinder")] %>% 
#   as.data.frame()
# 
# all_data <- 
#   all_data %>%
#   AddMetaData(dbl_finder_df)
# 
# all_data <- all_data[, all_data$scDblFinder.class == "singlet"]


### dimension reduction
all_data <- RunPCA(all_data)

pdf(file.path(sc_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_data, ndims = 30)
dev.off()

pdf(file.path(sc_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 1:5)
dev.off()
pdf(file.path(sc_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 6:10)
dev.off()
pdf(file.path(sc_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 11:15)
dev.off()
pdf(file.path(sc_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 16:20)
dev.off()

# saveRDS(all_data, "all_data_merge.rds")

# dr ---------------------
ndims <- 20
all_data <- RunUMAP(all_data, dims = 1:ndims)

pdf(file.path(sc_out_dir, str_glue("Dimplot dim{ndims} project.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE) + NoLegend()
dev.off()


all_data <- all_data %>% 
  harmony::RunHarmony(
    dims.use = 1:ndims,
    group.by.vars = "project",
    reduction.save = "harmony"
  )

all_data <- RunUMAP(
  all_data, 
  dims = 1:ndims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_"
)


pdf(file.path(sc_out_dir, str_glue("Dimplot dim{ndims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

### cluster  --------------------
all_data <-
  FindNeighbors(
    all_data,
    reduction = "pca",
    dims = 1:ndims,
    k.param = 15
  )
all_data <-
  FindClusters(
    all_data,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )


# dimplot
pdf(file.path(sc_out_dir, str_glue("Dimplot dim{ndims} cluster - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# dotplot
marker_ls <- c(
  "B cell" = "Cd19",
  "B cell" = "Cd79a", 
  "B cell" = "Ms4a1",
  
  "Plasma" = "Mzb1",
  "Plasma" = "Jchain",
  
  "T cell" = "Cd3e", 
  "T cell" = "Cd4",
  "T cell" = "Foxp3", 
  "T cell" = "Cd8a",
  "T cell" = "Gzmb",
  
  "NK" = "Ncam1",
  "NK" = "Gnly",
  "NK" = "Nkg7",
  
  "Endothelial" = "Pecam1",
  "Endothelial" = "Vwf",
  
  "Epithelial" = "Epcam",
  
  "Fibroblasts" = "Pdgfra",
  "Fibroblasts" = "Col1a1",
  "Fibroblasts" = "Dcn",
  "Fibroblasts" = "Acta2",
  
  "AT1" = "Ager",
  "AT2" = "Sftpa1",
  "AT2" = "Sftpd",
  
  "Ciliated cell" = "Cdhr3",
  "Ciliated cell" = "Cetn2",
  "Ciliated cell" = "Dnah11",
  
  "Club cell" = "Scgb1a1",
  "Club cell" = "Scgb3a1",
  "Club cell" = "Scgb3a2",
  
  
  "Smooth muscle cells" = "Tagln",
  "Smooth muscle cells" = "Myh11",
  
  "Pericyte" = "Rgs5",
  
  "Immune cells" = "Ptprc",
  
  "Monocyte" = "Cd14",
  "Monocyte" = "Fcgr3",
  "Monocyte" = "Vcan",
  "Monocyte" = "Fcna",
  "Monocyte" = "Csf1r",
  "Monocyte" = "Ly6c1", ### 
  
  "Macropahge" = "C1qa",
  "Macropahge" = "C1qb",
  "Macropahge" = "C1qc",
  "Macropahge" = "Cd68",
  
  "Alveolar Macrophage" = "Apoc1",
  "Alveolar Macrophage" = "Gpnmb",
  "Alveolar Macrophage" = "Inhba",
  "Alveolar Macrophage" = "Fcgr1", ### CD64
  "Alveolar Macrophage" = "Siglecf",
  # "Alveolar Macrophage" = "Itgax", # CD11c
  
  "Interstitial Macrophages" = "Adgre1", ### F4/80
  "Interstitial Macrophages" = "Mrc1", 
  
  "Neutrophil" = "Csf3r",
  "Neutrophil" = "S100a8",
  "Neutrophil" = "S100a9",
  "Neutrophil"  = "Ly6g", ### 
  "Neutrophil"  = "Itgam", # CD11b
  
  
  "Mast" = "Tpsab1",
  "Mast" = "Cpa3",
  
  "DC" = "Lyz1",
  "DC" = "Il3ra",
  "DC" = "Clec9a",
  "DC" = "Irf4",
  "cDC1" = "Xcr1",
  "cDC2" = "Itgax",
  "pDC"  = "Bst2",
  "pDC"  = "Siglech"
)



pdf(file.path(sc_out_dir, str_glue("dotplot dims{ndims}.pdf")), width = 35, height = 6)
DotPlot(all_data, group.by = "RNA_snn_res.0.2", features = marker_ls) + RotatedAxis()
DotPlot(all_data, group.by = "RNA_snn_res.0.4", features = marker_ls) + RotatedAxis()
DotPlot(all_data, group.by = "RNA_snn_res.0.6", features = marker_ls) + RotatedAxis()
DotPlot(all_data, group.by = "RNA_snn_res.0.8", features = marker_ls) + RotatedAxis()
dev.off()

# saveRDS(all_data, 'all_data_sc_dim15.rds')

# script sc_analysis.R --dims 15 --rds all_data_dim8.rds --save --k.param 15



library(Seurat)
library(tidyverse)

# all_data <- readRDS('all_data_sc_dim15.rds')

# major_cell_cluster <- c(
#   '0'  = "Fibroblast",
#   '1'  = "Endothelial",
#   '2'  = "Epithelial",
#   '3'  = "Smooth muscle cells",
#   '4'  = "Alveolar Macrophage",
#   '5'  = "NK/T",
#   '6'  = "Monocyte",
#   '7'  = "Monocyte",
#   '8'  = "Endothelial",
#   '9'  = "Fibroblast",
#   '10' = "B/Plasma",
#   '11' = "Smooth muscle cells",
#   '12' = "Epithelial",
#   '13' = "Macrophage",
#   '14' = "NK_Fibro",
#   '15' = "Epithelial",
#   '16' = "Endothelial",
#   '17' = "NK_Fibro",
#   '18' = "B/Plasma",
#   '19' = "NK_Monocyte",
#   '20' = "B_T"
# )


major_cell_cluster <- c(
  '0'  = "Fibroblasts",
  '1'  = "Endothelial",
  '2'  = "AT2",
  '3'  = "Fibroblasts",
  '4'  = "Interstitial Macrophage",
  '5'  = "NK/T",
  '6'  = "Neutrophil",
  '7'  = "Endothelial",
  '8'  = "Monocyte",
  '9'  = "Smooth muscle cells",
  '10' = "B/Plasma",
  '11' = "Fibroblasts",
  '12' = "AT1",
  '13' = "Alveolar Macrophage",
  '14' = "Fibroblasts",
  '15' = "Ciliated/Club cell",
  '16' = "Fibroblasts",
  '17' = "Endothelial",
  '18' = "DC",
  '19' = "B/Plasma" 
)

# doubele_cells <- c("NK_Fibro", "NK_Fibro", "NK_Monocyte", "B_T")

Idents(all_data) <- "RNA_snn_res.0.2"
all_data <- RenameIdents(all_data, major_cell_cluster)
all_data$annotation <- Idents(all_data)


pdf(file.path(sc_out_dir, str_glue("Dimplot dim{ndims} major annotation - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data[, !all_data$annotation %in% doubele_cells], group.by = "annotation", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data[, !all_data$annotation %in% doubele_cells], group.by = "annotation", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# find markers -------
Idents(all_data) <- "RNA_snn_res.0.2"
all_markers <-
  FindAllMarkers(
    JoinLayers(all_data),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )


all_de_genes <- 
  all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  pull("gene")

# DoHeatmap(sc, features = all_de_genes, group.by = "RNA_snn_res.0.8")
if(!exists('avg_expr')) avg_expr <- AverageExpression(JoinLayers(all_data))[[1]]


all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  nest() %>%
  deframe() %>%
  map(pull, 'gene') %>%
  imap(~{
    xx <- str_c(.x, collapse = ",")
    str_glue("cluster_{.y}:{xx}")
  }) %>%
  walk(cat, "\n")

avg_expr[all_de_genes, ] %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  setNames(., paste('cluster', (1:ncol(.) - 1) )) %>%  t %>%
  pheatmap::pheatmap(
    scale = "column",
    cluster_row = avg_expr[VariableFeatures(all_data), ] %>% t %>% dist %>% hclust,
    cluster_col = FALSE,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_col = NA
  )

### sub cluster ---------------

source("code base/sc_analysis_tools.R")

#### AT1 subcluster ----------------
sc_at1 <- all_data[, all_data$annotation == "AT1"]

sc_at1 <- sc_normalize(sc_at1)
sc_at1 <- sc_pca(
  sc_at1, 
  sc_out_dir = "subcluster/AT1/sc.1 pca",
  dir.create.force = TRUE
)

ndims_at1 <- 10

sc_at1 <- sc_dr(
  sc_at1, 
  ndims = ndims_at1,
  sc_out_dir = "subcluster/AT1/sc.2 dr",
  dir.create.force = TRUE
)


sc_at1 <- sc_cluster(
  sc_at1, 
  ndims = ndims_at1,
  reduction =  "pca", # "harmony",
  sc_out_dir = "subcluster/AT1/sc.3 cluster - pca",
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)


# sc_at1 <- sc_dotplot(
#   sc_at1, 
#   ndims = ndims_at1,
#   marker = at1_markers,
#   sc_out_dir = 'subcluster/AT1/sc.4 dotplot-pca', 
#   dir.create.force = TRUE
# )

all_markers_at1 <- 
  sc_findmarker(
    sc_at1,
    idents = "RNA_snn_res.0.4",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

# sc_findmarker(sc_at1, all_markers = all_markers_at1)

sc_mark_heatmap(
  sc_at1,
  markers = all_markers_at1,
  idents = "RNA_snn_res.0.4",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE
)

sc_mark_heatmap(
  sc_at1,
  markers = all_markers_at1,
  idents = "RNA_snn_res.0.4",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)

at1_markers_de <-
  all_markers_at1 %>% 
  group_by(cluster) %>%
  dplyr::filter(
    p_val_adj < 0.05,
    pct.1 - pct.2 > 0.1
  ) %>%
  dplyr::slice_max(avg_log2FC, n = 10) %>%
  dplyr::pull("gene") %>%
  unique()


DotPlot(
  sc_at1,
  features = at1_markers_de,
  group.by = "RNA_snn_res.0.4"
) +
  RotatedAxis()
  
sc_at1 <- sc_dotplot(
  sc_at1,
  ndims = ndims_at1,
  marker = at1_markers_de,
  sc_out_dir = 'subcluster/AT1/sc.4 dotplot-pca',
  dir.create.force = TRUE
)

dir.create('subcluster/AT1/sc.6 heatmap')
pdf(file.path('subcluster/AT1/sc.6 heatmap', str_glue("heatmap.pdf")), height = 4, width = 8)
sc_mark_heatmap(
  sc_at1,
  markers = all_markers_at1,
  idents = "RNA_snn_res.0.4",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE
)
dev.off()


ak1_cell_cluster <- c(
  '0'  = "Cldn3 AT1",
  '1'  = "Lars2 AT1",   
  '2'  = "Lars2 AT1",
  '3'  = "Mettl7a1 AT1",  
  '4'  = "Man1c1 AT1" 
)

FeaturePlot(sc_at1, features = "nCount_RNA")
FeaturePlot(sc_at1, features = "nFeature_RNA")
FeaturePlot(sc_at1, features = "percent.mt")

sc_at1_2 <- sc_at1[, !sc_at1$RNA_snn_res.0.4 %in% c('5', '6')]

Idents(sc_at1_2) <- "RNA_snn_res.0.4"
sc_at1_2 <- RenameIdents(sc_at1_2, ak1_cell_cluster)
sc_at1_2$annotation <- Idents(sc_at1_2)

dir.create('subcluster/AT1/sc.5 annotation')
pdf(file.path('subcluster/AT1/sc.5 annotation', str_glue("Dimplot dim{ndims_at1} at1 annotation.pdf")), height = 6, width = 8)
DimPlot(sc_at1_2, group.by = "annotation", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_at1_2, group.by = "annotation", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

saveRDS(sc_at1, file.path("subcluster", 'sc_at1.rds'))
saveRDS(sc_at1_2, file.path("subcluster", 'sc_at1_2.rds'))

##### cells props barplot -----
p_col_proj_at1 <-
  sc_at1_2[[]] %>%
  group_by(project, annotation) %>%
  summarise(count = n()) %>%
  group_by(project) %>% 
  mutate(freq = count / sum(count)) %>%
  mutate(project = factor(project, levels = c(
    str_glue("P{c(3, 7, 14)}_Normoxia"),
    str_glue("P{c(3, 7, 14)}_Hyperoxia")
  ))) %>% 
  
  ggplot(aes(x = project, y = freq, fill = annotation)) +
  geom_col() + 
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  RotatedAxis() +
  ggsci::scale_fill_aaas() +
  labs(
    x = NULL,
    y = NULL
  )

p_col_oxy_at1 <-
  sc_at1_2[[]] %>%
  group_by(Oxygen, annotation) %>%
  summarise(count = n()) %>%
  group_by(Oxygen) %>% 
  mutate(freq = count / sum(count)) %>%
  mutate(Oxygen = factor(Oxygen, levels = c("Normoxia", "Hyperoxia"))) %>% 
  
  ggplot(aes(x = Oxygen, y = freq, fill = annotation)) +
  geom_col() + 
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  RotatedAxis() +
  ggsci::scale_fill_aaas() +
  labs(
    x = NULL,
    y = NULL
  )

p_col_age_at1 <-
  sc_at1_2[[]] %>%
  group_by(Age, annotation) %>%
  summarise(count = n()) %>%
  group_by(Age) %>% 
  mutate(freq = count / sum(count)) %>%
  mutate(Age = factor(Age, levels = c("P3", "P7", "P14"))) %>% 
  
  ggplot(aes(x = Age, y = freq, fill = annotation)) +
  geom_col() + 
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  RotatedAxis() +
  ggsci::scale_fill_aaas() +
  labs(
    x = NULL,
    y = NULL
  )

# export
ggsave(
  p_col_proj_at1,
  file = file.path('subcluster/AT1', str_glue("cell proportion barplot - at1 annotation - proj.pdf") ),
  height = 6,
  width  = 8
)
ggsave(
  p_col_oxy_at1,
  file = file.path('subcluster/AT1', str_glue("cell proportion barplot - at1 annotation - oxygen.pdf") ),
  height = 6,
  width  = 3
)
ggsave(
  p_col_age_at1,
  file = file.path('subcluster/AT1', str_glue("cell proportion barplot - at1 annotation - age.pdf") ),
  height = 6,
  width  = 5
)




##### msigdb hallmark -----
h_geneset <- msigdbr::msigdbr(category = "H", species = "Mus musculus")

h_geneset$gs_name %>%
  unique() %>% 
  str_to_lower() %>%
  str_remove("^hallmark_")

h_geneset_ls <-
  h_geneset[c('gs_name', 'gene_symbol')] %>%
  with(split(gene_symbol, gs_name)) %>%
  {set_names(
    .,
    names(.) %>% 
      str_to_lower() %>%
      str_remove("^hallmark_")
  )}

sc_at1_score <-
  sc_at1_2 %>%
  JoinLayers() %>%
  AddModuleScore(
    features = h_geneset_ls,
    # name = names(h_geneset_ls)
  )

sc_at1_score %>% head

# msigdb anno
h_at1_score_anno <- 
  sc_at1_score[[]] %>%
  dplyr::select(annotation, Cluster1:Cluster50) %>%
  group_by(annotation) %>%
  summarise(across(everything(), mean))

p_at1_hallmark_1 <-
  pheatmap::pheatmap(
    h_at1_score_anno %>% dplyr::filter(annotation != "doublet") %>% column_to_rownames("annotation") %>% t,
    scale = 'row',
    labels_row = names(h_geneset_ls),
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(Inf, "RdBu")))(50),
    border_color = NA,
    silent = TRUE
  )

# export 
pdf(
  file.path("subcluster/AT1/", "Hallmark pathway score of 3 AT1 sub cluster cells.pdf"), 
  height = 12, 
  width  = 5
)
p_at1_hallmark_1$gtable %>% cowplot::plot_grid() %>% print()
dev.off()

##### de analysis of sub cluster of AT1 -------

all_markers_at1_df <-
  sc_at1_2 %>%
  .[, .$annotation != "doublet"] %>%
  JoinLayers() %>%
  {Idents(.) <- "annotation"; .} %>%
  FindAllMarkers(
    logfc.threshold = 0,
    min.pct = 0,
    return.thresh = 1.2
  )

all_markers_at1_df %>% head
# p_val avg_log2FC pct.1 pct.2     p_val_adj   cluster   gene
# Actg1  3.671699e-104  0.9636172 0.986 0.805 7.365061e-100 Cldn3 AT1  Actg1
# Malat1  1.330267e-93 -1.1808514 0.990 0.992  2.668383e-89 Cldn3 AT1 Malat1
# Cldn18  1.543076e-91  0.7951087 0.994 0.849  3.095256e-87 Cldn3 AT1 Cldn18
# H3f3b   2.011950e-82  0.7796797 0.994 0.859  4.035771e-78 Cldn3 AT1  H3f3b
# Ubb     5.249673e-78  0.8113941 0.980 0.734  1.053032e-73 Cldn3 AT1    Ubb
# Rplp0   4.980331e-77  0.8401779 0.975 0.736  9.990046e-73 Cldn3 AT1  Rplp0


# scatter
p_volcanoplot_AT1_ls <- 
  list("AT1" = all_markers_at1_df) %>%
  imap(function(dat, name){ # browser()
    
    set.seed(1234)
    rect_height <- 1
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
    
    
    # text label
    df_text_grp <-
      dat %>%
      group_by(cluster) %>%
      summarise() %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(x = x_base - 1) %>%
      mutate(x = (x_base + x)/2)
    
    df_scatter %>%
      ggplot(aes(x = xx, y = y)) +
      geom_rect(
        data = df_bg_rect, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        inherit.aes = FALSE
      ) +
      
      geom_point(aes(color = group), size = 1) +
      
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_text(
        data = df_text_grp,
        aes(x = x, y = 0, label = cluster)
      ) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      scale_fill_manual(values =  c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4")) +
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

pdf(
  file.path("subcluster/AT1/", "de scatter plot - AT1 annotation.pdf"), 
  height = 5, 
  width  = 16
)
p_volcanoplot_AT1_ls %>%  walk(print)
dev.off()



##### enrichment analysis of sub cluster of AT1 -------

all_geneid_AT1_ls <-
  list("AT1" = all_markers_at1_df) %>%
  map(function(dat){
    gene <-
      dat %>%
      dplyr::filter(avg_log2FC > 0.25, p_val_adj < 0.05, pct.1 > 0.1) %>%
      group_by(cluster) %>%
      nest %>%
      deframe() %>%
      map(~{
        print(nrow(.x))
        clusterProfiler::bitr(
          .x$gene,
          fromType = "SYMBOL",
          toType   = "ENTREZID",
          OrgDb    = "org.Mm.eg.db"
        ) %>%
          .$ENTREZID
      })
    
    gene  
  })

all_geneid_AT1_ls <- all_geneid_AT1_ls$AT1

all_enrich_AT1_ls <- 
  all_geneid_AT1_ls %>%
  map(function(x){
    clusterProfiler::enrichGO(
      gene = x,
      ont = "BP",
      OrgDb = "org.Mm.eg.db",
      pvalueCutoff  = 1, # 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff  = 1 # 0.2
    )
  })

# enrichment plot
p_all_enrich_AT1_ls <-
  all_enrich_AT1_ls %>%
  imap(function(enrich, name){ # browser()
    enrichplot::dotplot(
      enrich,
      x = "GeneRatio",
      color = "p.adjust",
      showCategory = 10,
      label_format = 40
    ) + 
      # scale_fill_gradient2(low = "gray", high = "red") +
      ggtitle(name)
  })



pdf(
  file.path("subcluster/AT1/", "de enrichment analysis dotplot - sub cluster of AT1.pdf"), 
  height = 6.5, 
  width  = 8
)
p_all_enrich_AT1_ls %>%  walk(print)
dev.off()


# export enrichemt table

all_enrich_AT1_df <-
  all_enrich_AT1_ls %>%
  map('result') %>%
  enframe(name = "cluster") %>%
  unnest(value) 

all_enrich_AT1_df$geneSymbol <-
  all_enrich_AT1_df$geneID %>%
  str_split("/") %>%
  map_chr(function(x){
    clusterProfiler::bitr(
      x,
      fromType = "ENTREZID",
      toType   = "SYMBOL",
      OrgDb    = "org.Mm.eg.db"
    ) %>%
      .$SYMBOL %>%
      str_c(collapse = "/")
  })

all_enrich_AT1_df %>%
  openxlsx::write.xlsx(
    file.path("subcluster/AT1/", "enrichment analysis table of sub cluster of AT1.xlsx"),
    asTable = TRUE
  )



#### AT2 cells subcluster ----------------
sc_at2 <- all_data[, all_data$annotation == "AT2"]

sc_at2 <- sc_normalize(sc_at2)
sc_at2 <- sc_pca(
  sc_at2, 
  sc_out_dir = "subcluster/AT2/sc.1 pca",
  dir.create.force = TRUE
)

ndims_at2 <- 20

sc_at2 <- sc_dr(
  sc_at2, 
  ndims = ndims_at2,
  sc_out_dir = "subcluster/AT2/sc.2 dr",
  dir.create.force = TRUE
)


sc_at2 <- sc_cluster(
  sc_at2, 
  ndims = ndims_at2,
  reduction =  "pca", # "harmony",
  sc_out_dir = "subcluster/AT2/sc.3 cluster",
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)


saveRDS(sc_at2, 'sc_at2.rds')

#### AT2_2 subcluster ----------------
sc_at2_2 <- sc_at2[, !sc_at2$RNA_snn_res.0.4 %in% c('4', '5', '6')]

sc_at2_2 <- sc_normalize(sc_at2_2)
sc_at2_2 <- sc_pca(
  sc_at2_2, 
  sc_out_dir = "subcluster/AT2_2/sc.1 pca",
  dir.create.force = TRUE
)

ndims_at2_2 <- 20

sc_at2_2 <- sc_dr(
  sc_at2_2, 
  ndims = ndims_at2_2,
  sc_out_dir = "subcluster/AT2_2/sc.2 dr",
  dir.create.force = TRUE
)


sc_at2_2 <- sc_cluster(
  sc_at2_2, 
  ndims = ndims_at2_2,
  reduction =  "pca", # "harmony",
  sc_out_dir = "subcluster/AT2_2/sc.3 cluster",
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.3, 0.4, 0.6, 0.8)
)


# sc_at2_2 <- sc_dotplot(
#   sc_at2_2, 
#   ndims = ndims_at2_2,
#   marker = at2_2_markers,
#   sc_out_dir = 'subcluster/AT2_2/sc.4 dotplot', 
#   dir.create.force = TRUE
# )

all_markers_at2_2 <- 
  sc_findmarker(
    sc_at2_2,
    idents = "RNA_snn_res.0.3",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

# sc_findmarker(sc_at2_2, all_markers = all_markers_at2_2)

at2_2_markers_de <-
  all_markers_at2_2 %>% 
  group_by(cluster) %>%
  dplyr::filter(
    p_val_adj < 0.05,
    pct.1 - pct.2 > 0.1
  ) %>%
  dplyr::slice_max(avg_log2FC, n = 10) %>%
  dplyr::pull("gene") %>%
  unique()


DotPlot(
  sc_at2_2,
  features = at2_2_markers_de,
  group.by = "RNA_snn_res.0.3"
) +
  RotatedAxis()

sc_at2_2 <- sc_dotplot(
  sc_at2_2,
  ndims = ndims_at2_2,
  marker = at2_2_markers_de,
  sc_out_dir = 'subcluster/AT2_2/sc.4 dotplot-pca',
  dir.create.force = TRUE
)

dir.create('subcluster/AT2_2/sc.6 heatmap')
pdf(file.path('subcluster/AT2_2/sc.6 heatmap', str_glue("heatmap.pdf")), height = 4, width = 8)
sc_mark_heatmap(
  sc_at2_2,
  markers = all_markers_at2_2,
  idents = "RNA_snn_res.0.3",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE
)
dev.off()



 
sc_mark_heatmap(
  sc_at2_2,
  markers = all_markers_at2_2,
  idents = "RNA_snn_res.0.3",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE
)

sc_mark_heatmap(
  sc_at2_2,
  markers = all_markers_at2_2,
  idents = "RNA_snn_res.0.3",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)


at2_2_cell_cluster <- c(
  '0'  = "Dlk1 AT2",
  '1'  = "Areg AT2",
  '2'  = "Acoxl AT2",
  '3'  = "Phgdh AT2"
)

VlnPlot(sc_at2_2, features = "Phgdh", group.by = "RNA_snn_res.0.3")


Idents(sc_at2_2) <- "RNA_snn_res.0.3"
sc_at2_2 <- RenameIdents(sc_at2_2, at2_2_cell_cluster)
sc_at2_2$annotation <- Idents(sc_at2_2)

dir.create('subcluster/AT2_2/sc.5 annotation')
pdf(file.path('subcluster/AT2_2/sc.5 annotation', str_glue("Dimplot dim{ndims_at2_2} harmony AT2 annotation.pdf")), height = 6, width = 8)
DimPlot(sc_at2_2, group.by = "annotation", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_at2_2, group.by = "annotation", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

saveRDS(sc_at2_2, file.path("subcluster", 'sc_at2_2.rds'))
 
##### cells props barplot -----
p_col_proj_at2 <-
  sc_at2_2[[]] %>%
  group_by(project, annotation) %>%
  summarise(count = n()) %>%
  group_by(project) %>% 
  mutate(freq = count / sum(count)) %>%
  mutate(project = factor(project, levels = c(
    str_glue("P{c(3, 7, 14)}_Normoxia"),
    str_glue("P{c(3, 7, 14)}_Hyperoxia")
  ))) %>% 
  
  ggplot(aes(x = project, y = freq, fill = annotation)) +
  geom_col() + 
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  RotatedAxis() +
  ggsci::scale_fill_aaas() +
  labs(
    x = NULL,
    y = NULL
  )

p_col_oxy_at2 <-
  sc_at2_2[[]] %>%
  group_by(Oxygen, annotation) %>%
  summarise(count = n()) %>%
  group_by(Oxygen) %>% 
  mutate(freq = count / sum(count)) %>%
  mutate(Oxygen = factor(Oxygen, levels = c("Normoxia", "Hyperoxia"))) %>% 
  
  ggplot(aes(x = Oxygen, y = freq, fill = annotation)) +
  geom_col() + 
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  RotatedAxis() +
  ggsci::scale_fill_aaas() +
  labs(
    x = NULL,
    y = NULL
  )

p_col_age_at2 <-
  sc_at2_2[[]] %>%
  group_by(Age, annotation) %>%
  summarise(count = n()) %>%
  group_by(Age) %>% 
  mutate(freq = count / sum(count)) %>%
  mutate(Age = factor(Age, levels = c("P3", "P7", "P14"))) %>% 
  
  ggplot(aes(x = Age, y = freq, fill = annotation)) +
  geom_col() + 
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  RotatedAxis() +
  ggsci::scale_fill_aaas() +
  labs(
    x = NULL,
    y = NULL
  )

# export
ggsave(
  p_col_proj_at2,
  file = file.path('subcluster/AT2_2/', str_glue("cell proportion barplot - at2 annotation - proj.pdf") ),
  height = 6,
  width  = 8
)
ggsave(
  p_col_oxy_at2,
  file = file.path('subcluster/AT2_2/', str_glue("cell proportion barplot - at2 annotation - oxygen.pdf") ),
  height = 6,
  width  = 3
)
ggsave(
  p_col_age_at2,
  file = file.path('subcluster/AT2_2/', str_glue("cell proportion barplot - at2 annotation - age.pdf") ),
  height = 6,
  width  = 5
)

##### msigdb hallmark -----
sc_at2_2_score <-
  sc_at2_2 %>%
  JoinLayers() %>%
  AddModuleScore(
    features = h_geneset_ls,
    # name = names(h_geneset_ls)
  )

sc_at2_2_score %>% head

# msigdb anno
h_at2_score_anno <- 
  sc_at2_2_score[[]] %>%
  dplyr::select(annotation, Cluster1:Cluster50) %>%
  group_by(annotation) %>%
  summarise(across(everything(), mean))

p_at2_hallmark_1 <-
  pheatmap::pheatmap(
    h_at2_score_anno %>% column_to_rownames("annotation") %>% t,
    scale = 'row',
    labels_row = names(h_geneset_ls),
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(Inf, "RdBu")))(50),
    border_color = NA,
    silent = TRUE
  )

# export 
pdf(
  file.path("subcluster/AT2_2/", "Hallmark pathway score of 3 AT2 sub cluster cells.pdf"), 
  height = 12, 
  width  = 5
)
p_at2_hallmark_1$gtable %>% cowplot::plot_grid() %>% print()
dev.off()

##### de analysis of sub cluster of AT2 -------

all_markers_at2_df <-
  sc_at2_2 %>%
  .[, .$annotation != "doublet"] %>%
  JoinLayers() %>%
  {Idents(.) <- "annotation"; .} %>%
  FindAllMarkers(
    logfc.threshold = 0,
    min.pct = 0,
    return.thresh = 1.2
  )

all_markers_at2_df %>% head
# p_val avg_log2FC pct.1 pct.2     p_val_adj        cluster     gene
# EBF1      0.000000e+00   3.184000 0.690 0.179  0.000000e+00 S3 fibroblasts     EBF1
# MIR99AHG 2.867768e-269   2.796024 0.619 0.153 7.660382e-265 S3 fibroblasts MIR99AHG
# KCNN3    5.304018e-180   4.171928 0.391 0.073 1.416809e-175 S3 fibroblasts    KCNN3
# ADGRL3   4.410816e-151  -4.038387 0.083 0.451 1.178217e-146 S3 fibroblasts   ADGRL3
# HPSE2    1.850722e-148  -3.198233 0.107 0.485 4.943647e-144 S3 fibroblasts    HPSE2
# C7       3.038981e-142   2.497229 0.484 0.165 8.117726e-138 S3 fibroblasts       C7


# scatter
p_volcanoplot_AT2_ls <- 
  list("AT2" = all_markers_at2_df) %>%
  imap(function(dat, name){ # browser()
    
    set.seed(1234)
    rect_height <- 1
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
    
    
    # text label
    df_text_grp <-
      dat %>%
      group_by(cluster) %>%
      summarise() %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(x = x_base - 1) %>%
      mutate(x = (x_base + x)/2)
    
    df_scatter %>%
      ggplot(aes(x = xx, y = y)) +
      geom_rect(
        data = df_bg_rect, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        inherit.aes = FALSE
      ) +
      
      geom_point(aes(color = group), size = 1) +
      
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_text(
        data = df_text_grp,
        aes(x = x, y = 0, label = cluster)
      ) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      scale_fill_manual(values =  c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4")) +
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

pdf(
  file.path("subcluster/AT2_2/", "de scatter plot - AT2 annotation.pdf"), 
  height = 5, 
  width  = 16
)
p_volcanoplot_AT2_ls %>%  walk(print)
dev.off()



##### enrichment analysis of sub cluster of AT1 -------

all_geneid_AT2_ls <-
  list("AT2" = all_markers_at2_df) %>%
  map(function(dat){
    gene <-
      dat %>%
      dplyr::filter(avg_log2FC > 0.25, p_val_adj < 0.05, pct.1 > 0.1) %>%
      group_by(cluster) %>%
      nest %>%
      deframe() %>%
      map(~{
        print(nrow(.x))
        clusterProfiler::bitr(
          .x$gene,
          fromType = "SYMBOL",
          toType   = "ENTREZID",
          OrgDb    = "org.Mm.eg.db"
        ) %>%
          .$ENTREZID
      })
    
    gene  
  })

all_geneid_AT2_ls <- all_geneid_AT2_ls$AT2

all_enrich_AT2_ls <- 
  all_geneid_AT2_ls %>%
  map(function(x){
    clusterProfiler::enrichGO(
      gene = x,
      ont = "BP",
      OrgDb = "org.Mm.eg.db",
      pvalueCutoff  = 1, # 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff  = 1 # 0.2
    )
  })

# enrichment plot
p_all_enrich_AT2_ls <-
  all_enrich_AT2_ls %>%
  imap(function(enrich, name){ # browser()
    enrichplot::dotplot(
      enrich,
      x = "GeneRatio",
      color = "p.adjust",
      showCategory = 10,
      label_format = 40
    ) + 
      # scale_fill_gradient2(low = "gray", high = "red") +
      ggtitle(name)
  })



pdf(
  file.path("subcluster/AT2_2/", "de enrichment analysis dotplot - sub cluster of AT2.pdf"), 
  height = 6.5, 
  width  = 8
)
p_all_enrich_AT2_ls %>%  walk(print)
dev.off()


# export enrichemt table

all_enrich_AT2_df <-
  all_enrich_AT2_ls %>%
  map('result') %>%
  enframe(name = "cluster") %>%
  unnest(value) 

all_enrich_AT2_df$geneSymbol <-
  all_enrich_AT2_df$geneID %>%
  str_split("/") %>%
  map_chr(function(x){
    clusterProfiler::bitr(
      x,
      fromType = "ENTREZID",
      toType   = "SYMBOL",
      OrgDb    = "org.Mm.eg.db"
    ) %>%
      .$SYMBOL %>%
      str_c(collapse = "/")
  })

all_enrich_AT2_df %>%
  openxlsx::write.xlsx(
    file.path("subcluster/AT2_2/", "enrichment analysis table of sub cluster of AT2.xlsx"),
    asTable = TRUE
  )



#### NK/T cells subcluster ----------------
sc_t <- all_data[, all_data$annotation == "NK/T"]

sc_t <- sc_normalize(sc_t)
sc_t <- sc_pca(
  sc_t, 
  sc_out_dir = "subcluster/t/sc.1 pca",
  dir.create.force = TRUE
)

ndims_t <- 20

sc_t <- sc_dr(
  sc_t, 
  ndims = ndims_t,
  sc_out_dir = "subcluster/t/sc.2 dr",
  dir.create.force = TRUE
)


sc_t <- sc_cluster(
  sc_t, 
  ndims = ndims_t,
  reduction =  "pca", # "harmony",
  sc_out_dir = "subcluster/t/sc.3 cluster",
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)

t_markers <- c(
  "Cd3d",
  "Cd3e",
  "Cd4",
  "Cd8a",
  "Foxp3",
  "Il7r",
  "Trac",
  "Trdc",
  "Ncam1",
  "Nkg7",
  
  'Mki67',
  'Top2a',
  
  # 'Malat1','H3f3b','Trbc2','Eif1','Satb1','H3f3a','Rbm39','Ubb','Sftpc',
  'Mier1',

  
  'Hbb-bt', 
  'Hbb-bs'
)

sc_t <- sc_dotplot(
  sc_t, 
  ndims = ndims_t,
  marker = t_markers,
  sc_out_dir = 'subcluster/t/sc.4 dotplot', 
  dir.create.force = TRUE,
  width = 25
)

all_markers_t <- 
  sc_findmarker(
    sc_t,
    idents = "RNA_snn_res.0.2",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

sc_findmarker(sc_t, all_markers = all_markers_t, sort_by = "pct")

sc_mark_heatmap(
  sc_t,
  markers = all_markers_t,
  idents = "RNA_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE
)

sc_mark_heatmap(
  sc_t,
  markers = all_markers_t,
  idents = "RNA_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)


t_cell_cluster <- c(
  '0'   = "CD4T cells",
  '1'   = "CD8T cells",
  '2'   = "NK cells",
  '3'   = "T proliferating cells",
  '4'   = "Treg cells",
  '5'   = "ILC2",
  '6'   = "gd T cells",
  '7'   = "Mier1 T cells",
  '8'   = "high-hbb cells"
)



Idents(sc_t) <- "RNA_snn_res.0.2"
sc_t <- RenameIdents(sc_t, t_cell_cluster)
sc_t$annotation <- Idents(sc_t)

dir.create('subcluster/t/sc.5 annotation')
pdf(file.path('subcluster/t/sc.5 annotation', str_glue("Dimplot dim{ndims} harmony t annotation.pdf")), height = 6, width = 8)
DimPlot(sc_t, group.by = "annotation", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_t, group.by = "annotation", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_t[, sc_t$annotation != "high-hbb cells"], group.by = "annotation", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_t[, sc_t$annotation != "high-hbb cells"], group.by = "annotation", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

saveRDS(sc_t, 'sc_t.rds')

#### B and Plasma cells subcluster ----------------
sc_b <- all_data[, all_data$annotation == "B/Plasma"]

sc_b <- sc_normalize(sc_b)
sc_b <- sc_pca(
  sc_b, 
  sc_out_dir = "subcluster/b/sc.1 pca",
  dir.create.force = TRUE
)

ndims_b <- 10

sc_b <- sc_dr(
  sc_b, 
  ndims = ndims_b,
  sc_out_dir = "subcluster/b/sc.2 dr",
  dir.create.force = TRUE
)


sc_b <- sc_cluster(
  sc_b, 
  ndims = ndims_b,
  reduction =  "pca", # "harmony",
  sc_out_dir = "subcluster/b/sc.3 cluster",
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)

b_markers <- c(
  "Cd19",
  "Ms4a1",
  "Cd79a",
  
  'Cd24a', # activate b
  'Ighd', # naive
  
  "Jchain", # plasma
  
  'Irf4', # plasmablast
  'Prdm1',
  
  'Cd38',
  
  "Igkc",
  "Igha",
  'Ighm',
  
  'Epcam', 'Pecam1', 'Col1a1', 'Cd3e'
)
sc_b <- sc_dotplot(
  sc_b, 
  ndims = ndims_b,
  marker = b_markers,
  sc_out_dir = 'subcluster/b/sc.4 dotplot', 
  dir.create.force = TRUE
)

all_markers_b <- 
  sc_findmarker(
    sc_b,
    idents = "RNA_snn_res.0.2",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

# sc_findmarker(sc_b, all_markers = all_markers_b)

sc_mark_heatmap(
  sc_b,
  markers = all_markers_b,
  idents = "RNA_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE
)

sc_mark_heatmap(
  sc_b,
  markers = all_markers_b,
  idents = "RNA_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)


b_cell_cluster <- c(
  '0'  = "B naive cells",
  '1'  = "B naive cells",
  '2'  = "B activated cells",
  '3'  = "Plasmablasts",
  '4'  = "doublet"
)



Idents(sc_b) <- "RNA_snn_res.0.2"
sc_b <- RenameIdents(sc_b, b_cell_cluster)
sc_b$annotation <- Idents(sc_b)

dir.create('subcluster/b/sc.5 annotation')
pdf(file.path('subcluster/b/sc.5 annotation', str_glue("Dimplot dim{ndims_b} harmony b annotation.pdf")), height = 6, width = 8)
DimPlot(sc_b, group.by = "annotation", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_b, group.by = "annotation", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_b[,sc_b$annotation!="doublet"], group.by = "annotation", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_b[,sc_b$annotation!="doublet"], group.by = "annotation", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

saveRDS(sc_b, 'sc_b.rds')

#### endoloid cells subcluster ----------------
sc_endo <- all_data[, all_data$annotation == "Endothelial"]

sc_endo <- sc_normalize(sc_endo)
sc_endo <- sc_pca(
  sc_endo, 
  sc_out_dir = "subcluster/endo/sc.1 pca",
  dir.create.force = TRUE
)

ndims_endo <- 20

sc_endo <- sc_dr(
  sc_endo, 
  ndims = ndims_endo,
  sc_out_dir = "subcluster/endo/sc.2 dr",
  dir.create.force = TRUE
)


sc_endo <- sc_cluster(
  sc_endo, 
  ndims = ndims_endo,
  reduction =  "pca", # "harmony",
  sc_out_dir = "subcluster/endo/sc.3 cluster",
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)

endo_markers <- c(
  'Epcam',
  'Pecam1',
  'Col1a1',
  'Cd3d',
  'Lyz',
  
  'Capillary 1 cells' = 'Plvap', # 'PLVAP',
  'Capillary 1 cells' = 'Gpihbp1', 
  'Capillary 2 cells' = 'Car4', # 'CAR4',
  'Capillary 2 cells' = 'Apln', 
  
  'Arterial endothelial cells' = "Dkk2", # DKK2,
  'Arterial endothelial cells' = "Gja5", # GJA5,
  'Arterial endothelial cells' = "Bmx", # BMX,
  
  'Venous endothelial cells' = "Ackr1", # ACKR1,
  'Venous endothelial cells' = "Ephb4", # EPHB4,
  'Venous endothelial cells' = "Hdac9", # HDAC9,
  
  'Lymphatic endothelial cells' = "Prox1", # PROX1
  'Lymphatic endothelial cells' = "Mmrn1", # MMRN1
  'Lymphatic endothelial cells' = "Ccl21", # CCL21
  'Lymphatic endothelial cells' = "Lyve1" # LYVE1
  
  # "ACKR1", # ACKR1
  # "Col15a1", # COL15A1
  # "Abcb1", # ABCB1
  # "Vwa1" # VWA1
  
  
)

sc_endo <- sc_dotplot(
  sc_endo, 
  ndims = ndims_endo,
  marker = endo_markers,
  sc_out_dir = 'subcluster/endo/sc.4 dotplot', 
  dir.create.force = TRUE,
  width = 15
)

all_markers_endo <- 
  sc_findmarker(
    sc_endo,
    idents = "RNA_snn_res.0.2",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

# sc_findmarker(sc_endo, all_markers = all_markers_endo)

sc_mark_heatmap(
  sc_endo,
  markers = all_markers_endo,
  idents = "RNA_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE
)

sc_mark_heatmap(
  sc_endo,
  markers = all_markers_endo,
  idents = "RNA_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)


endo_cell_cluster <- c(
  '0'   = "Capillary 1 cells",
  '1'   = "Capillary 2 cells",
  '2'   = "Capillary 1 cells",
  '3'   = "Arterial endothelial cells",
  '4'   = "Venous endothelial cells",
  '5'   = "Capillary 1 cells",
  '6'   = "Lymphatic endothelial cells",
  '7'   = "Capillary 2 cells",
  '8'   = "Capillary 2 cells"
)



Idents(sc_endo) <- "RNA_snn_res.0.2"
sc_endo <- RenameIdents(sc_endo, endo_cell_cluster)
sc_endo$annotation <- Idents(sc_endo)

dir.create('subcluster/endo/sc.5 annotation')
pdf(file.path('subcluster/endo/sc.5 annotation', str_glue("Dimplot dim{ndims_endo} harmony endo annotation.pdf")), height = 6, width = 8)
DimPlot(sc_endo, group.by = "annotation", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_endo, group.by = "annotation", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

saveRDS(sc_endo, 'sc_endo.rds')


#### stroma cells subcluster ----------------
sc_stma <- all_data[, all_data$annotation %in% c("Smooth muscle cells", "Fibroblasts")]

sc_stma <- sc_normalize(sc_stma)
sc_stma <- sc_pca(
  sc_stma, 
  sc_out_dir = "subcluster/stroma/sc.1 pca",
  dir.create.force = TRUE
)

ndims_stroma <- 20

sc_stma <- sc_dr(
  sc_stma, 
  ndims = ndims_stroma,
  sc_out_dir = "subcluster/stroma/sc.2 dr",
  dir.create.force = TRUE
)


sc_stma <- sc_cluster(
  sc_stma, 
  ndims = ndims_stroma,
  reduction =  "harmony", # "harmony",
  sc_out_dir = "subcluster/stroma/sc.3 cluster",
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)


stroma_markers <- c(
  'Pdgfra',
  'Tcf21',
  'Col1a1',
  'Dcn',
  # Alveolar Fibroblast 2 cells
  # https://www.lungmap.net/research/cell-cards/?cell_cards_id=LMCC0000000026
  'Mfap5',
  'Scara5',
  Alveolar = 'Scube2',
  Alveolar = 'Npnt',
  Alveolar = 'Ces1d',
  Alveolar = 'Wnt2',
  Infammatory = 'Saa3',
  Infammatory = 'Lcn2',
  'Stress-activated' = 'Cdkn1a',
  'Stress-activated' = 'Gdf15',
  Proliferating = 'Mki67',
  Fibrotic = 'Cthrc1',
  Adventitial = 'Pi16',
  Peribronchial = 'Hhip',
  Peribronchial = 'Wnt5a',
  'Smooth muscle' = 'Myh11',
  Pericyte = 'Cox4i2',
  Mesothelial = 'Msln'
)



sc_stma <- sc_dotplot(
  sc_stma,
  ndims = ndims_stroma,
  marker = stroma_markers,
  sc_out_dir = 'subcluster/stroma/sc.4 dotplot',
  dir.create.force = TRUE,
  width = 25
)

all_markers_stroma <- 
  sc_findmarker(
    sc_stma,
    idents = "RNA_snn_res.0.2",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

# sc_findmarker(sc_stma, all_markers = all_markers_stroma, sort_by = 'p')

sc_mark_heatmap(
  sc_stma,
  markers = all_markers_stroma,
  idents = "RNA_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE
)

sc_mark_heatmap(
  sc_stma,
  markers = all_markers_stroma,
  idents = "RNA_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)


stroma_cell_cluster <- c(
  '0'   = "Alveolar fibroblasts",
  '1'   = "Peribronchial fibroblasts",
  '2'   = "Adventitial fibroblasts",
  '3'   = "Smooth muscle cells",
  '4'   = "Alveolar fibroblasts",
  '5'   = "Proliferating fibroblasts",
  '6'   = "Pericyte",
  '7'   = "Mesothelial",
  '8'   = "Peribronchial fibroblasts"
)



Idents(sc_stma) <- "RNA_snn_res.0.2"
sc_stma <- RenameIdents(sc_stma, stroma_cell_cluster)
sc_stma$annotation <- Idents(sc_stma)

dir.create('subcluster/stroma/sc.5 annotation')
pdf(file.path('subcluster/stroma/sc.5 annotation', str_glue("Dimplot dim{ndims} harmony stroma annotation.pdf")), height = 6, width = 8)
DimPlot(sc_stma, group.by = "annotation", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_stma, group.by = "annotation", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

saveRDS(sc_stma, 'sc_stma.rds')


# all_ annotation ----------
meta_data_full <- rbind(
  sc_at1_2[['annotation']],
  sc_at2_2[['annotation']],
  sc_t[['annotation']],
  sc_b[['annotation']],
  sc_endo[['annotation']],
  sc_stma[['annotation']]
) %>%
  setNames("annotation2")

meta_data_full2 <-
  all_data[[c("annotation")]] %>%
  rownames_to_column("cell_id") %>%
  left_join(
    rownames_to_column(meta_data_full, "cell_id"),
    by = "cell_id"
  ) %>%
  mutate(
    annotation2 = case_when(
      annotation == "AT1" & is.na(annotation2) ~ "AT1 others cells",
      annotation == "AT2" & is.na(annotation2) ~ "AT2 others cells",
      is.na(annotation2) ~ annotation,
      TRUE ~ annotation2
    )
  )



all_data_sc <-
  all_data %>%
  AddMetaData(meta_data_full2)

other_cells <- c(
  "AT1 others cells",
  "AT2 others cells",
  "high-hbb cells",
  "doublet"
)
all_data_sc2 <- 
  all_data_sc %>%
  .[, ! .$annotation2 %in%other_cells]

all_data_sc %>% DimPlot(group.by = "annotation2", label = TRUE, repel = TRUE)
all_data_sc2 %>% DimPlot(group.by = "annotation2", label = TRUE, repel = TRUE)

saveRDS(all_data_sc,  'all_data_sc.rds')
saveRDS(all_data_sc2, 'all_data_sc2.rds')



# plot all annotation ---------
coors_mini_grob <-
  grid::gList(
    grid::linesGrob(
      x = unit(c(0, 0.3), "npc"),
      y = unit(c(0.03, 0.03), "npc"),
      arrow = grid::arrow(type = "closed"),
      gp = grid::gpar(fill = "black")
    ),
    grid::linesGrob(
      x = unit(c(0.03, 0.03), "npc"),
      y = unit(c(0, 0.3), "npc"),
      arrow = grid::arrow(type = "closed"),
      gp = grid::gpar(fill = "black")
    ),
    grid::textGrob(
      "UMAP_1",
      x = unit(0.15, "npc"),
      y = unit(0.0, "npc")
    ),
    grid::textGrob(
      "UMAP_2",
      x = unit(0.0, "npc"),
      y = unit(0.15, "npc"),
      rot = 90
    )
  )

p_umap_major <-
  all_data_sc2 %>%
  {
    .$annotation <-
    .$annotation %>%
    fct_recode(
      Fibroblasts = 'Smooth muscle cells',
      Fibroblasts = 'Fibroblasts'
    ) %>%
    as.character()
    .
  } %>%
  DimPlot(
    group.by = "annotation",
    cols = col_plan_1,
    label = TRUE
  ) +
  theme_void() +
  NoLegend() 

p_umap_final <-
  all_data_sc2 %>%
  DimPlot(
    group.by = "annotation2",
    cols = col_plan_1,
    label = TRUE,
    repel = TRUE
  ) +
  theme_void() +
  NoLegend() 

p_umap_final_combn <-
  cowplot::ggdraw() +
  cowplot::draw_plot(p_umap_final, x = 0.1, y = 0.1, height = 0.8, width = 0.8) +
  cowplot::draw_grob(coors_mini_grob, x = 0.05, y = 0.05, height = 0.95, width = 0.95, scale = 0.8, halign = 0, valign = 0)


ggsave(
  p_umap_final_combn,
  file = file.path(sc_out_dir, "Dimplot - annotation - final - minicoord.pdf"),
  height = 6, 
  width  = 8
)








