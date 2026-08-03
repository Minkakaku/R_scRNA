library(Seurat)
library(tidyverse)

sc_out_dir <- "sc_output" 
dir.create(sc_out_dir)

file_dir <- "GSE214611/"
all_files_dir <- list.files(file_dir, pattern = "sn")


all_files <- file.path(file_dir, all_files_dir)
names(all_files) <- all_files_dir

sample_info <-
  data.frame(
    sample = all_files_dir,
    group  = all_files_dir %>% str_remove("\\_.*$") %>%  str_remove("sn"),
    id = seq_along(all_files_dir)
  )

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


# remove snd1_2
all_data_ls_pass[['snd1_2']] <- NULL
all_data <- merge(all_data_ls_pass[[1]], all_data_ls_pass[-1])
all_data$group <- str_extract(all_data$project, "d[0-9]")

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

# saveRDS(all_data, "all_data_sc_merge.rds")

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
  Endothelial = "Pecam1",
  Endothelial = "Vwf",
  Endothelial = "Apoe",
  Endo_ECs = "Npr3", # endocardial ECs  # https://www.nature.com/articles/s41467-022-35691-7
  vas_ECs  = "Fabp4", # vascular ECs
  
  # Epi  = "Epcam", 
  Epi = 'Wt1', # epicardial cell 
  Epi = 'Tbx18',
  Epi = 'Aldh1a2',
  
  "Fibroblasts" = "Col1a1",
  "Fibroblasts" = "Col1a2",
  # "Fibroblasts" = "DCN",
  # "Fibroblasts" = "ACTA2",
  "Fibroblasts" = "Postn",
  "Mural cells" = "Pdgfrb",
  SMC = "Myh11",
  SMC = "Tagln",
  # Pdgfrb # Pericytes
  
  # "Pericytes" = "LFHP",
  "Pericytes" = "Rgs5",
  "Pericytes" = "Mcam",
  
  # "Cardiomyocytes" = "TNNT2",
  "CMs" = "Actn2",
  "CMs" = "Myh7",
  "CMs" = "Ttn",
  "Atrial_CMs" = "Sln",
  "Ven_CMs" = "Myl2",
  
  
  Immune = 'Ptprc',
  
  "B cells" = "Cd19",
  "B cells" = "Cd79a",
  "B cells" = "Ms4a1",
  
  Erythrocytes = "Hba-a1",
  Erythrocytes = "Hba-a2",
  
  "T cells" = "Cd3e",
  "T cells" = "Cd3g",
  "T cells" = "Trac",
  "T cells" = "Cd4",
  "T cells" = "Cd8A",
  "T cells" = "Bcl11b",
  
  "NK" = "Gnly",
  "NK" = "Nkg7",
  "NK" = "Ncam1",
  "NK" = "Unc13d",
  
  "Proliferating" = "Mki67",
  "Proliferating" = "Top2a",
  
  "Myeloid" = "Itgam",
  "Myeloid" = "Cd14", 
  "Myeloid" = "Fcgr3a",
  "Myeloid" = "Cd68",
  "Myeloid" = "Cd163",
  "Myeloid" = "C1qa",
  "Myeloid" = "Adgre1",
  
  # "Mast cell" = "CPA3",
  # "Mast cell" = "CTSG",
  # "Mast cell" = "TPSAB1",
  # "Mast cell" = "MS4A2",
  # "Mast cell" = "KIT",
  
  
  "Dendritic cell" = "H2-Oa",
  
  'Neutrophi' = "S100a8",
  'Neutrophi' = "S100a9",
  'Neutrophi' = "Csf3r",
  
  ribosome = "Rps23",
  ribosome = "Rpl11",
  ribosome = "Rplp1",
  ribosome = "Rps4x",
  ribosome = "Rps9",
  
  mitochondrial = "mt-Nd1",
  mitochondrial = "mt-Nd2",
  mitochondrial = "mt-Co1",
  mitochondrial = "mt-Co2",
  mitochondrial = "mt-Atp8",
  mitochondrial = "mt-Atp6",
  
  hsp = "Hspb1",
  hsp = "Hspb7",
  hsp = "Hspa1l",
  hsp = "Hspa8",
  hsp = "Hspb9",
  hsp = "Hsp90b1",
  
  "Nppa"
  # 'Actb',
  # 'Myl6',
  # 'Mgp'
  
)



pdf(file.path(sc_out_dir, str_glue("dotplot dims{ndims}.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "RNA_snn_res.0.2", features = marker_ls) + RotatedAxis()
DotPlot(all_data, group.by = "RNA_snn_res.0.4", features = marker_ls) + RotatedAxis()
DotPlot(all_data, group.by = "RNA_snn_res.0.6", features = marker_ls) + RotatedAxis()
DotPlot(all_data, group.by = "RNA_snn_res.0.8", features = marker_ls) + RotatedAxis()
dev.off()


# script sc_analysis.R --dims 15 --rds all_data_dim8.rds --save --k.param 15

# # major cluster ---------------
# major_marker <- c(
#   'Epithelial cells' = 'EPCAM', 
#   'Epithelial cells' = 'AQP8',
#   'Epithelial cells' = 'BEST4',
#   'Epithelial cells' = 'MUC2', 
#   'Epithelial cells' = 'OLFM4', 
#   # 'Epithelial cells' = 'PLCG2', 
#   'Epithelial cells' = 'TRPM5', 
#   'Epithelial cells' = 'ZG16',
#   
#   'B and Plasma cells' = 'CD79A',
#   'B and Plasma cells' = 'BANK1',
#   'B and Plasma cells' = 'CD19',
#   'B and Plasma cells' = 'DERL3',
#   'B and Plasma cells' = 'MS4A1',
#   'B and Plasma cells' = 'MZB1',
#   
#   'T cells' = 'CD3D',
#   'T cells' = 'CD3E',
#   'T cells' = 'CD3G',
#   'T cells' = 'CD8A',
#   'T cells' = 'FOXP3',
#   'T cells' = 'GZMA',
#   'T cells' = 'GZMB',
#   'T cells' = 'IL17A',
#   'T cells' = 'NKG7',
#   'T cells' = 'TRBC1',
#   
#   'Stromal cells' = 'ACTA2',
#   'Stromal cells' = 'ADAMDEC1',
#   'Stromal cells' = 'CHI3L1',
#   'Stromal cells' = 'COL3A1',
#   'Stromal cells' = 'NRXN1',
#   'Stromal cells' = 'PVALP',
#   'Stromal cells' = 'SOX6',
#   'Stromal cells' = 'VWF',
#   
#   'Myeloid cells' = 'AIF1',
#   'Myeloid cells' = 'C1QA',
#   'Myeloid cells' = 'C1QB',
#   'Myeloid cells' = 'CD14',
#   'Myeloid cells' = 'CMTM2',
#   'Myeloid cells' = 'FCGR3B',
#   'Myeloid cells' = 'LYZ',
#   'Myeloid cells' = 'MS4A2',
#   'Myeloid cells' = 'TPSAB1',
#   'Myeloid cells' = 'TPSAB2'
# )
# 
# 
# 
# pdf(file.path(sc_out_dir, str_glue("dotplot dims{ndims} major marker.pdf")), width = 25, height = 6)
# DotPlot(all_data, group.by = "RNA_snn_res.0.2", features = major_marker) + RotatedAxis()
# DotPlot(all_data, group.by = "RNA_snn_res.0.4", features = major_marker) + RotatedAxis()
# DotPlot(all_data, group.by = "RNA_snn_res.0.6", features = major_marker) + RotatedAxis()
# DotPlot(all_data, group.by = "RNA_snn_res.0.8", features = major_marker) + RotatedAxis()
# dev.off()


library(Seurat)
library(tidyverse)



major_cell_cluster <- c(
  '0'  = 'Fibroblast',
  '1'  = 'Vascular EC',
  '2'  = 'Immune cell',
  '3'  = 'CM',
  '4'  = 'CM',
  '5'  = 'Pericyte/SMC',
  '6'  = 'Fibroblast',
  '7'  = 'Endocardial EC',
  '8'  = 'Immune-like EC',
  '9'  = 'high-mt cell',
  '10' = 'CM',
  '11' = 'Pericyte/SMC',
  '12' = 'EC',
  '13' = 'Proliferating immune cell',
  '14' = 'Epicardial cell'
)



Idents(all_data) <- "RNA_snn_res.0.4"
all_data <- RenameIdents(all_data, major_cell_cluster)
all_data$annotation <- Idents(all_data)


pdf(file.path(sc_out_dir, str_glue("Dimplot dim{ndims} major annotation - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


saveRDS(all_data, str_glue('all_data_sc_dim{ndims}.rds'))
