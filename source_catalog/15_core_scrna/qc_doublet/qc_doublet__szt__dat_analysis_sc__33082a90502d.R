library(Seurat)
library(tidyverse)

sc_out_dir <- "sc_output" # dir.create(sc_out_dir)

file_dir <- "/data/outputdata/outputdata-sc/cellranger_out/FWSC20240319"
all_files_dir <- list.files(file_dir, pattern = "^[XYZ]-.*$")


all_files <- file.path(file_dir, all_files_dir, "step3/filtered_feature_bc_matrix")
names(all_files) <- all_files_dir

sample_info <-
  data.frame(
    sample = all_files_dir,
    group  = all_files_dir %>% str_extract("^[XYZ]"),
    id = seq_along(all_files_dir)
  )

# 
# # 6 samples
# stopifnot(length(all_files) == 6)

all_data_ls <- 
  all_files %>%
  map(Read10X) %>%
  imap(~CreateSeuratObject(
    .x,
    project = .y,
    assay = "RNA",
    names.field = 1,
    names.delim = "_",
    min.cells = 10,
    min.features = 10
  ))



all_data_ls <-
  all_data_ls %>%
  map(function(x){
    x$project <- x@project.name
    x
  })

### QC ---------------------
all_data_ls <- 
  all_data_ls %>%
  map(function(x){
    x[["percent.mt"]]  <- PercentageFeatureSet(x, pattern = "^MT-")
    x[["percent.ribo"]]  <- PercentageFeatureSet(x, pattern = "^RP[SL]")
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

# #### add cell_id prefix
all_data_ls_pass <-
  all_data_ls_pass %>%
  names() %>%
  set_names(., .) %>%
  map(function(x){
    all_data <- all_data_ls_pass[[x]]
    ids <- sample_info$id[sample_info$sample == x]

    all_data <- RenameCells(all_data, add.cell.id = ids)
    all_data
  })

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
ndims <- 15
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
    reduction = "harmony",
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
    
    Endothelial = "PECAM1",
    Endothelial = "VWF",
    Endothelial = "AQP1",
    
    Epithelial  = "EPCAM", 
    Epithelial  = "KRT5",
    
    Enterocyte  = "SLC26A3",
    Enterocyte  = "CA1",
    Enterocyte  = "GUCA2A",
    Enterocyte  = "GUCA2B",
    
    Goblet = "MUC2",
    Goblet = "TFF3",
    
    Enteroendocrine = "CHGA",
    
    # "Fibroblasts" = "SFRP4", 
    "Fibroblasts" = "PDGFRA", 
    "Fibroblasts" = "CCDC80",
    "Fibroblasts" = "COL1A1",
    "Fibroblasts" = "DCN",
    "Fibroblasts" = "ACTA2",
    
    "Pericytes" = "PDGFRB",
    "Pericytes" = "RGS5",
    "Pericytes" = "MCAM",
    
    
    # "Smooth muscle cells" = "PDGFRB", 
    # "Smooth muscle cells" = "PLN", 
    # "Smooth muscle cells" = "MYLK", 
    # "Smooth muscle cells" = "MYH11", 
    # "Smooth muscle cells" = "ACTA2",
    
    
    # "Myefibroblast" = "PDGFRB", 
    "Myefibroblast" = "PLN", 
    "Myefibroblast" = "MYLK", 
    "Myefibroblast" = "MYH11", 
    # "Myefibroblast" = "ACTA2",
    
    # "Breast glandular cells" = "SERPINA3", 
    # "Breast glandular cells" = "DSC2",
    
    Immune = 'PTPRC',
    
    "B cells" = "CD19",
    "B cells" = "CD79A",
    "B cells" = "MS4A1",
    "Plasma" = "TNFRSF17",
    "Plasma" = "JCHAIN",
    "Plasma" = "MZB1",
    "Plasma" = "XBP1",
    
    # Plasmablast = "MZB1",
    # Plasmablast = "XBP1",
    
    "T cells" = "CD3E",
    "T cells" = "CD3G",
    "T cells" = "CD4",
    "T cells" = "CD8A",
    "T cells" = "TRAC",
    
    "NK" = "GNLY",
    "NK" = "NKG7",
    "NK" = "NCAM1",
    
    "Proliferating" = "MKI67",
    "Proliferating" = "TOP2A",
    
    "Myeloid" = "ITGAM",
    "Macrophage" = "CD68",
    "Macrophage" = "CD163",
    "Macrophage" = "C1QA",
    "Macrophage" = "C1QC",
    
    "Mast cell" = "CPA3",
    "Mast cell" = "CTSG",
    "Mast cell" = "TPSAB1",
    "Mast cell" = "MS4A2",
    "Mast cell" = "KIT",
    
    
    "Monocytes" = "CD14", 
    "Monocytes" = "FCGR3A",
    
    "Dendritic cell" = "CD83",
    'Dendritic cell' = "LILRA4",  
    
    'Dendritic cell' = "CLEC9A",
    'Dendritic cell' = "IDO1",
    
    'Dendritic cell' = "CD1C",
    'Dendritic cell' = "FCER1A",
    
    'Neutrophi' = "S100A8",
    'Neutrophi' = "S100A9",
    'Neutrophi' = "CSF3R",
    'Neutrophi' = "FCGR3B",
    
    Gial = "MPZ"
  )


pdf(file.path(sc_out_dir, str_glue("dotplot dims{ndims}.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "RNA_snn_res.0.2", features = marker_ls) + RotatedAxis()
DotPlot(all_data, group.by = "RNA_snn_res.0.4", features = marker_ls) + RotatedAxis()
DotPlot(all_data, group.by = "RNA_snn_res.0.6", features = marker_ls) + RotatedAxis()
DotPlot(all_data, group.by = "RNA_snn_res.0.8", features = marker_ls) + RotatedAxis()
dev.off()

# saveRDS(all_data, 'all_data_sc_dim15.rds')

# script sc_analysis.R --dims 15 --rds all_data_dim8.rds --save --k.param 15

# major cluster ---------------
major_marker <- c(
  'Epithelial cells' = 'EPCAM', 
  'Epithelial cells' = 'AQP8',
  'Epithelial cells' = 'BEST4',
  'Epithelial cells' = 'MUC2', 
  'Epithelial cells' = 'OLFM4', 
  # 'Epithelial cells' = 'PLCG2', 
  'Epithelial cells' = 'TRPM5', 
  'Epithelial cells' = 'ZG16',
  
  'B and Plasma cells' = 'CD79A',
  'B and Plasma cells' = 'BANK1',
  'B and Plasma cells' = 'CD19',
  'B and Plasma cells' = 'DERL3',
  'B and Plasma cells' = 'MS4A1',
  'B and Plasma cells' = 'MZB1',
  
  'T cells' = 'CD3D',
  'T cells' = 'CD3E',
  'T cells' = 'CD3G',
  'T cells' = 'CD8A',
  'T cells' = 'FOXP3',
  'T cells' = 'GZMA',
  'T cells' = 'GZMB',
  'T cells' = 'IL17A',
  'T cells' = 'NKG7',
  'T cells' = 'TRBC1',
  
  'Stromal cells' = 'ACTA2',
  'Stromal cells' = 'ADAMDEC1',
  'Stromal cells' = 'CHI3L1',
  'Stromal cells' = 'COL3A1',
  'Stromal cells' = 'NRXN1',
  'Stromal cells' = 'PVALP',
  'Stromal cells' = 'SOX6',
  'Stromal cells' = 'VWF',
  
  'Myeloid cells' = 'AIF1',
  'Myeloid cells' = 'C1QA',
  'Myeloid cells' = 'C1QB',
  'Myeloid cells' = 'CD14',
  'Myeloid cells' = 'CMTM2',
  'Myeloid cells' = 'FCGR3B',
  'Myeloid cells' = 'LYZ',
  'Myeloid cells' = 'MS4A2',
  'Myeloid cells' = 'TPSAB1',
  'Myeloid cells' = 'TPSAB2'
)



pdf(file.path(sc_out_dir, str_glue("dotplot dims{ndims} major marker.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "RNA_snn_res.0.2", features = major_marker) + RotatedAxis()
DotPlot(all_data, group.by = "RNA_snn_res.0.4", features = major_marker) + RotatedAxis()
DotPlot(all_data, group.by = "RNA_snn_res.0.6", features = major_marker) + RotatedAxis()
DotPlot(all_data, group.by = "RNA_snn_res.0.8", features = major_marker) + RotatedAxis()
dev.off()


library(Seurat)
library(tidyverse)

# all_data <- readRDS('all_data_sc_dim15.rds')


major_cell_cluster <- c(
  '0'  = "T cells",
  '1'  = "B and Plasma cells",
  '2'  = "Epithelial cells",
  '3'  = "Stromal cells",
  '4'  = "B and Plasma cells",
  '5'  = "Myeloid cells",
  '6'  = "Myeloid cells",
  '7'  = "Stromal cells",
  '8'  = "Stromal cells",
  '9'  = "Myeloid cells",
  '10' = "Epithelial cells",
  '11' = "Stromal cells",
  '12' = "T cells"
)



Idents(all_data) <- "RNA_snn_res.0.2"
all_data <- RenameIdents(all_data, major_cell_cluster)
all_data$annotation <- Idents(all_data)


pdf(file.path(sc_out_dir, str_glue("Dimplot dim{ndims} major annotation - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

### sub cluster ---------------

source("sc_analysis_tools.R")

#### Epi subcluster ----------------
sc_epi <- all_data[, all_data$annotation == "Epithelial cells"]

sc_epi <- sc_normalize(sc_epi)
sc_epi <- sc_pca(
  sc_epi, 
  sc_out_dir = "subcluster/epi/sc.1 pca",
  dir.create.force = TRUE
)

ndims_epi <- 10

sc_epi <- sc_dr(
  sc_epi, 
  ndims = ndims_epi,
  sc_out_dir = "subcluster/epi/sc.2 dr",
  dir.create.force = TRUE
)


sc_epi <- sc_cluster(
  sc_epi, 
  ndims = ndims_epi,
  reduction =  "pca", # "harmony",
  sc_out_dir = "subcluster/epi/sc.3 cluster - pca",
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)


sc_epi <- sc_dotplot(
  sc_epi, 
  ndims = ndims_epi,
  marker = epi_markers,
  sc_out_dir = 'subcluster/epi/sc.4 dotplot-pca', 
  dir.create.force = TRUE
)

all_markers_epi <- 
  sc_findmarker(
    sc_epi,
    idents = "RNA_snn_res.0.4",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

# sc_findmarker(sc_epi, all_markers = all_markers_epi)

sc_mark_heatmap(
  sc_epi,
  markers = all_markers_epi,
  idents = "RNA_snn_res.0.4",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE
)

sc_mark_heatmap(
  sc_epi,
  markers = all_markers_epi,
  idents = "RNA_snn_res.0.4",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)
epi_cell_cluster <- c(
  '0'  = "SLC10A2 Laminin colonocytes",
  '1'  = "Goblet cells",
  '2'  = "OLFM4 Goblet cells",
  '3'  = "Mature Goblet cells",
  '4'  = "PLCG colonocytes",
  '5'  = "Cycling transit-Amplifying cells",
  '6'  = "IL7R colonocytes",
  '7'  = "Laminin colonocytes",
  '8'  = "OLFM4 Goblet cells",
  '9'  = "Enteroendocrine cells"
)



Idents(sc_epi) <- "RNA_snn_res.0.4"
sc_epi <- RenameIdents(sc_epi, epi_cell_cluster)
sc_epi$annotation <- Idents(sc_epi)

dir.create('subcluster/epi/sc.5 annotation')
pdf(file.path('subcluster/epi/sc.5 annotation', str_glue("Dimplot dim{ndims} epi annotation.pdf")), height = 6, width = 8)
DimPlot(sc_epi, group.by = "annotation", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_epi, group.by = "annotation", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

saveRDS(sc_epi, 'sc_epi.rds')


#### B and Plasma cells subcluster ----------------
sc_b <- all_data[, all_data$annotation == "B and Plasma cells"]

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
  reduction =  "harmony", # "harmony",
  sc_out_dir = "subcluster/b/sc.3 cluster",
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
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
  '0'  = "Plama cells",
  '1'  = "B cells",
  '2'  = "Naive B cells",
  '3'  = "plasmablasts",
  '4'  = "plasmablasts",
  '5'  = "Cycling cells(PC and B cells)",
  '6'  = "Lambda cluster"
)



Idents(sc_b) <- "RNA_snn_res.0.2"
sc_b <- RenameIdents(sc_b, b_cell_cluster)
sc_b$annotation <- Idents(sc_b)

dir.create('subcluster/b/sc.5 annotation')
pdf(file.path('subcluster/b/sc.5 annotation', str_glue("Dimplot dim{ndims} harmony b annotation.pdf")), height = 6, width = 8)
DimPlot(sc_b, group.by = "annotation", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_b, group.by = "annotation", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

saveRDS(sc_b, 'sc_b.rds')

#### T cells subcluster ----------------
sc_t <- all_data[, all_data$annotation == "T cells"]

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
  reduction =  "harmony", # "harmony",
  sc_out_dir = "subcluster/t/sc.3 cluster",
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)


sc_t <- sc_dotplot(
  sc_t, 
  ndims = ndims_t,
  marker = t_markers,
  sc_out_dir = 'subcluster/t/sc.4 dotplot', 
  dir.create.force = TRUE,
  width = 35
)

all_markers_t <- 
  sc_findmarker(
    sc_t,
    idents = "RNA_snn_res.0.8",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

# sc_findmarker(sc_t, all_markers = all_markers_t)

sc_mark_heatmap(
  sc_t,
  markers = all_markers_t,
  idents = "RNA_snn_res.0.8",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE
)

sc_mark_heatmap(
  sc_t,
  markers = all_markers_t,
  idents = "RNA_snn_res.0.8",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)


t_cell_cluster <- c(
  '0'   = "CD8T.1",
  '1'   = "CD4T.1",
  '2'   = "CD4T naive",
  '3'   = "CD8T.2",
  '4'   = "Tregs",
  '5'   = "CD4T.2",
  '6'   = "CD8T TRM",
  '7'   = "CD8T effector cells",
  '8'   = "NK",
  '9'   = "CD4T naive",
  '10'  = "CD4T naive",
  '11'  = "Cycling cells (T cells)",
  '12'  = "ILC3",
  '13'  = "Others T"
)



Idents(sc_t) <- "RNA_snn_res.0.8"
sc_t <- RenameIdents(sc_t, t_cell_cluster)
sc_t$annotation <- Idents(sc_t)

dir.create('subcluster/t/sc.5 annotation')
pdf(file.path('subcluster/t/sc.5 annotation', str_glue("Dimplot dim{ndims} harmony t annotation.pdf")), height = 6, width = 8)
DimPlot(sc_t, group.by = "annotation", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_t, group.by = "annotation", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

saveRDS(sc_t, 'sc_t.rds')


#### myeloid cells subcluster ----------------
sc_mye <- all_data[, all_data$annotation == "Myeloid cells"]

sc_mye <- sc_normalize(sc_mye)
sc_mye <- sc_pca(
  sc_mye, 
  sc_out_dir = "subcluster/mye/sc.1 pca",
  dir.create.force = TRUE
)

ndims_mye <- 20

sc_mye <- sc_dr(
  sc_mye, 
  ndims = ndims_mye,
  sc_out_dir = "subcluster/mye/sc.2 dr",
  dir.create.force = TRUE
)


sc_mye <- sc_cluster(
  sc_mye, 
  ndims = ndims_mye,
  reduction =  "harmony", # "harmony",
  sc_out_dir = "subcluster/mye/sc.3 cluster",
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)


sc_mye <- sc_dotplot(
  sc_mye, 
  ndims = ndims_mye,
  marker = mye_markers,
  sc_out_dir = 'subcluster/mye/sc.4 dotplot', 
  dir.create.force = TRUE,
  width = 35
)

all_markers_mye <- 
  sc_findmarker(
    sc_mye,
    idents = "RNA_snn_res.0.4",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

# sc_findmarker(sc_mye, all_markers = all_markers_mye)

sc_mark_heatmap(
  sc_mye,
  markers = all_markers_mye,
  idents = "RNA_snn_res.0.4",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE
)

sc_mark_heatmap(
  sc_mye,
  markers = all_markers_mye,
  idents = "RNA_snn_res.0.8",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)


mye_cell_cluster <- c(
  '0'   = "M2 Macrophage",
  '1'   = "Inflammatory monocytes",
  '2'   = "Mast cells",
  '3'   = "Neutrophils",
  '4'   = "T cells", # <- 
  '5'   = "Fibroblast", # <-
  '6'   = "Basophil", # <-
  '7'   = "Plasma cell", # <-
  '8'   = "Memory B cell", # <-
  '9'   = "DCs",
  '10'  = "cDC1",
  '11'  = "Mast cells",
  '12'  = "Mast cells",
  '13'  = "Others",
  '14'  = "pDCs"
)



Idents(sc_mye) <- "RNA_snn_res.0.4"
sc_mye <- RenameIdents(sc_mye, mye_cell_cluster)
sc_mye$annotation <- Idents(sc_mye)

dir.create('subcluster/mye/sc.5 annotation')
pdf(file.path('subcluster/mye/sc.5 annotation', str_glue("Dimplot dim{ndims} harmony mye annotation.pdf")), height = 6, width = 8)
DimPlot(sc_mye, group.by = "annotation", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_mye, group.by = "annotation", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

saveRDS(sc_mye, 'sc_mye.rds')


#### stroma cells subcluster ----------------
sc_stma <- all_data[, all_data$annotation == "Stromal cells"]

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


sc_stma <- sc_dotplot(
  sc_stma, 
  ndims = ndims_stroma,
  marker = stroma_markers,
  sc_out_dir = 'subcluster/stroma/sc.4 dotplot', 
  dir.create.force = TRUE,
  width = 35
)

all_markers_stroma <- 
  sc_findmarker(
    sc_stma,
    idents = "RNA_snn_res.0.4",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

# sc_findmarker(sc_stma, all_markers = all_markers_stroma, sort_by = 'p')

sc_mark_heatmap(
  sc_stma,
  markers = all_markers_stroma,
  idents = "RNA_snn_res.0.4",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE
)

sc_mark_heatmap(
  sc_stma,
  markers = all_markers_stroma,
  idents = "RNA_snn_res.0.8",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)


stroma_cell_cluster <- c(
  '0'   = "S3 fibroblasts",
  '1'   = "S1 fibroblasts",
  '2'   = "Endothelial cells",
  '3'   = "Myofibroblasts",
  '4'   = "T cells",
  '5'   = "Immediate early-response (IER) fibroblasts",
  '6'   = "B cells",
  '7'   = "Glial cells",
  '8'   = "Pericytes",
  '9'   = "Macrophages",
  '10'  = "Memory B cell",
  '11'  = "Peri-cryptal S2 fibroblasts",
  '12'  = "Others",
  '13'  = "Lymphatic endothelium",
  '14'  = "Fat cell",
  '15'  = "Others"
)



Idents(sc_stma) <- "RNA_snn_res.0.4"
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
    sc_epi[['annotation']],
    sc_b[['annotation']],
    sc_t[['annotation']],
    sc_mye[['annotation']],
    sc_stma[['annotation']]
  ) %>%
  setNames("annotation2")

all_data_sc <-
  all_data %>%
  AddMetaData(meta_data_full)

all_data_sc %>% DimPlot(group.by = "annotation2")

saveRDS(all_data_sc, 'all_data_sc.rds')

