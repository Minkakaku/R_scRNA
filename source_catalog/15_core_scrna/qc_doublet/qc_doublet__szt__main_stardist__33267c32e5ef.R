library(Seurat)
library(tidyverse)
library(SeuratDisk)

## set future -------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 20 * 1024 ^ 10) # 10Gb

# setting -------------------------------------------------------------
## set dir --------------
img_output <- "materials_for_rmd" # dir.create(img_output)
qc_outdir <- "1. QC" # dir.create(qc_outdir)
dr_outdir <- "2. dimension reduction and cluster"  # dir.create(dr_outdir)
sp_output <- "3. spatial plot" # dir.create(sp_output)
mk_output <- "4. marker gene" # dir.create(mk_output)
er_output <- "5. enriment analysis" # dir.create(er_output)

## set default ------------
orginism <- "hsa" # "hsa" # "mm"
orginism_kegg <- "hsa" # "hsa" # "mmu"

mt_pattern   <- if(orginism == 'hsa') "^MT-" else "^mt-"
ribo_pattern <- if(orginism == 'hsa') "^RP[SL]" else "^Rp[sl]"

org.db <- if(orginism == 'hsa') "org.Hs.eg.db" else "org.Mm.eg.db"


#h5ad文件转化成rds用conda env tools里的sceasy包convertFormat()
#library('sceasy')
#sc <- convertFormat('adata_celltype.h5ad', to = 'seurat')
#saveRDS(sc, 'adata_celltype.rds')

NEPC <- readRDS('stardist_output/NEPC/NEPC.rds')
NEPC <- RenameCells(NEPC, NEPC$id) 
NEPC$nCount_RNA <- NEPC$total_counts
NEPC$nFeature_RNA <- NEPC$nFeaturess_RNA_by_counts
ADPC <- readRDS('stardist_output/ADPC/ADPC.rds')
ADPC <- RenameCells(ADPC, ADPC$id)
ADPC$nCount_RNA <- ADPC$total_counts
ADPC$nFeature_RNA <- ADPC$nFeaturess_RNA_by_counts

visium_ls <- list(ADPC, NEPC)
names(visium_ls) <- c('ADPC', 'NEPC')

## add sample_id or sample_name -----------------------
visium_ls <-
  visium_ls %>%
  imap(function(visium, name){
    visium$project <- name
    Idents(visium) <- 'project'
    RenameCells(visium, add.cell.id = name)
  })


visium_ls <- 
  visium_ls %>%
  # map(~ {DefaultAssay(.x) <- assay_name; .x}) %>%
  map(NormalizeData) 




# qc ---------------------------------------------------------------------------
visium_ls <- 
  visium_ls %>%
  map(function(visium){
    subset(visium, nCount_RNA > 0)
  })

visium_ls <- 
  visium_ls %>%
  map(function(x){# browser()
    x[["percent.mt"]]    <- PercentageFeatureSet(x, pattern = mt_pattern)
    x[["percent.ribo"]]  <- PercentageFeatureSet(x, pattern = ribo_pattern)
    x$log10GenesPerUMI   <- log10(x[['nFeature_RNA']][[1]] + 1) / log10(x[['nCount_RNA']][[1]] + 1)
    x
  })

assay_name <- 'RNA'
nfeature_name <- 'nFeature_RNA'
ncount_name <- 'nCount_RNA'

p_qc_ls <-
  visium_ls %>%
  imap(function(x, name){ # browser()
    # only subset 10,000 cells
    x_sub <- subset(x, cells = Cells(x[[assay_name]]), downsample = 100000, seed = 1)
    p <- VlnPlot(
      x_sub, 
      features = c(nfeature_name, ncount_name, "percent.mt", "log10GenesPerUMI"), 
      ncol     = 4, 
      pt.size  = 0,
      log      = FALSE, 
      raster   = TRUE
    ) +
      patchwork::plot_annotation(title = name)
    
    # 精细划分y轴标度
    p & scale_y_continuous(
      breaks = function(x){
        xx <- round(x)
        
        # 5代表分为5个breaks
        # 10代表以10的倍数进行划分breaks
        steps <- if(xx[2] < 100){
          scaled_times <- ceiling(100/xx[2])
          
          floor((scaled_times * xx[2] / 5 ) / 10 ) * 10 / scaled_times
        }else{
          floor((xx[2] / 5 ) / 10 ) * 10
        }
        
        # 最小值都是0
        bk <- seq(0, xx[2], by = steps)
        
        # 第一间隔细分
        bk_minor <- seq(bk[1], bk[2], length.out = 5)
        bk <- c(bk_minor, bk[-(1:2)])
        
        # 最后间隔细分
        bk <- rev(bk)
        bk_minor <- seq(bk[1], bk[2], length.out = 5)
        bk <- c(bk_minor, bk[-(1:2)])
        bk <- rev(bk)
        
        bk
      }
    )
  })

p_qc_log_ls <-
  visium_ls %>%
  imap(function(x, name){
    # only subset 10,000 cells
    x_sub <- subset(x, cells = Cells(x[[assay_name]]), downsample = 100000, seed = 1)
    p <- VlnPlot(
      x_sub, 
      features = c(nfeature_name, ncount_name, "percent.mt", "log10GenesPerUMI"), 
      ncol     = 4, 
      pt.size  = 0,
      log      = TRUE,
      raster   = TRUE
    ) +
      patchwork::plot_annotation(title = name)
    
  })

p_qc_ls %>%
  iwalk(function(p, name){
    png(file.path(img_output, str_glue("1 qc plot - {name}_stardist.png")), width = 480 * 3, height = 480 )
    print(p)
    dev.off()
  })


pdf(
  file.path(qc_outdir, "qc plot_stardist.pdf"),
  height = 6, 
  width = 20
)
p_qc_ls %>% walk(print)
dev.off()

pdf(
  file.path(qc_outdir, "qc plot - log_stardist.pdf"),
  height = 6, 
  width = 20
)
p_qc_log_ls %>% walk(print)
dev.off()




# qc filter --------------------------------------------------------------------
qc_metric <- read.csv("qc_metric_stardist.csv") %>% map_dfc(str_trim)

qc_metric_variable <- c(
  nfeature_name, ncount_name,
  "percent.mt", "log10GenesPerUMI"
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
visim_ls_pass <-
  names(visium_ls) %>%
  set_names(.,.) %>%
  map(function(sample){ #browser()
    all_data <- visium_ls[[sample]]
    qc <- qc_metric_ls[[sample]]
    
    e <- str_glue("subset(
      all_data,
      {qc}
    )")
    
    eval(parse(text = e))
  })


# qc plot
qc_plot_filter_ls <-
  names(visium_ls) %>%
  setNames(., .) %>%
  map(function(sample){
    all_data <- visium_ls[[sample]] %>% subset(downsample = 100000, seed = 1)
    all_data_pass <- visim_ls_pass[[sample]] %>% subset(downsample = 100000, seed  = 1)
    
    # all_data
    plot1 <- FeatureScatter(all_data, feature1 = ncount_name, feature2 = nfeature_name) +
      patchwork::plot_annotation(title = sample)
    plot2 <- FeatureScatter(all_data, feature1 = ncount_name, feature2 = "percent.mt") +
      patchwork::plot_annotation(title = sample)
    plot3 <- FeatureScatter(all_data, feature1 = ncount_name, feature2 = "log10GenesPerUMI") +
      patchwork::plot_annotation(title = sample)
    p1 <- plot1 + plot2 + plot3 + patchwork::plot_layout(ncol = 3)
    
    # all_data_pass
    plot1_pass <- FeatureScatter(all_data_pass, feature1 = ncount_name, feature2 = nfeature_name)+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    plot2_pass <- FeatureScatter(all_data_pass, feature1 = ncount_name, feature2 = "percent.mt")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    plot3_pass <- FeatureScatter(all_data_pass, feature1 = ncount_name, feature2 = "log10GenesPerUMI")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    p2 <- plot1_pass + plot2_pass + plot3_pass  + patchwork::plot_layout(ncol = 3)
    
    p1 / p2
  })



# export for rmd
qc_plot_filter_ls %>%
  iwalk(function(p, name){
    png(file.path(img_output, str_glue("1 qc plot - pass - {name}_stardist.png")), width = 360 * 3, height = 360 * 2 )
    print(p)
    dev.off()
  })


# export 
pdf(
  file.path(qc_outdir, "qc plot - pass_stardist.pdf"),
  height = 12, 
  width  = 18
)
qc_plot_filter_ls %>% walk(print)
dev.off()

qc_metric_ls %>%
  enframe("sample", "qc") %>%
  unnest(qc) %>% 
  openxlsx::write.xlsx(
    file = file.path(qc_outdir, "qc_metric_stardist.xlsx"),
    asTable = TRUE,
    withFilter = FALSE
  )

# cells filter
cell_counts_df <-
  names(visium_ls) %>%
  set_names(., .) %>%
  map(function(sample){
    all_data <- visium_ls[[sample]]
    all_data_pass <- visim_ls_pass[[sample]]
    
    structure(
      c(ncol(all_data), ncol(all_data_pass)),
      names = c("Before filter", "After filter")
    )
  }) %>%
  as.data.frame(check.names = FALSE)


# export
cell_counts_df %>% 
  rownames_to_column(" ") %>%
  openxlsx::write.xlsx(
    file = file.path(qc_outdir, "cell count after filter_stardist.xlsx"),
    asTable = TRUE,
    withFilter = FALSE
  )
#saveRDS(visim_ls_pass, 'visim_ls_pass_stardist.rds')




# dimension reduction and cluster and annotation -------------------------------
all_data <- merge(visim_ls_pass[[1]], visim_ls_pass[-1])

### normalize for sketch data ---------------------------------------------------
all_data <- NormalizeData(all_data, scale.factor = median(all_data[[ncount_name]] %>% unlist))
all_data <- FindVariableFeatures(all_data, nfeatures = 2000)
all_data <- ScaleData(all_data,vars.to.regress = ncount_name)

all_data <- RunPCA(all_data, features = VariableFeatures(all_data))

pdf(file = file.path(dr_outdir, 'pca plot_stardist.pdf'), height = 6, width = 8)
ElbowPlot(all_data, ndims = 30)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 1-5_stardist.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 1:5)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 6-10_stardist.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 6:10)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 11-15_stardist.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 11:15)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 16-20_stardist.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 16:20)
# dev.off()

dims <- 10

all_data <- RunUMAP(all_data, dims = 1:dims, return.model = TRUE)

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} project_stardist.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "umap") + NoLegend()
dev.off()


all_data <- all_data %>%
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    reduction.save = "harmony"
  )

all_data <- RunUMAP(
  all_data,
  dims = 1:dims,
  return.model = TRUE,
  reduction = "harmony",
  reduction.name = "harmonyumap",
  reduction.key  = "harmonyUMAP_"
)

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} project - harmony_stardist.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster ------
all_data <-
  FindNeighbors(
    all_data,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )

all_data <-
  FindClusters(
    all_data,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )


# dimplot
pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_stardist.pdf")), height = 6, width = 8)
# DimPlot(all_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()




# reFiltering -------------------------------------------------------------------------------
all_data <-
  FindClusters(
    all_data,
    resolution = c(1.0, 1.2, 1.5)
  )
pdf(file.path(dr_outdir, str_glue("Dimplot dims{dims} {reduction.use} cluster_split_stardist.pdf")), width = 100, height = 5)
DimPlot(all_data, group.by = "RNA_snn_res.0.8", split.by = "RNA_snn_res.0.8", reduction = 'harmonyumap')
DimPlot(all_data, group.by = "RNA_snn_res.1", split.by = "RNA_snn_res.1", reduction = 'harmonyumap')
DimPlot(all_data, group.by = "RNA_snn_res.1.2", split.by = "RNA_snn_res.1.2", reduction = 'harmonyumap')
DimPlot(all_data, group.by = "RNA_snn_res.1.5", split.by = "RNA_snn_res.1.5", reduction = 'harmonyumap')
dev.off()

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_stardist.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.1", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

Idents(all_data) <- "RNA_snn_res.0.8"
all_data2 <- subset(all_data, idents = c(0:17))

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_stardist2.pdf")), height = 6, width = 8)
DimPlot(all_data2, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data2, group.by = "RNA_snn_res.1", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data2, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data2, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(dr_outdir, str_glue("Dimplot dims{dims} {reduction.use} cluster_split_stardist2.pdf")), width = 100, height = 5)
DimPlot(all_data2, group.by = "RNA_snn_res.0.8", split.by = "RNA_snn_res.0.8", reduction = 'harmonyumap')
DimPlot(all_data2, group.by = "RNA_snn_res.1", split.by = "RNA_snn_res.1", reduction = 'harmonyumap')
DimPlot(all_data2, group.by = "RNA_snn_res.1.2", split.by = "RNA_snn_res.1.2", reduction = 'harmonyumap')
DimPlot(all_data2, group.by = "RNA_snn_res.1.5", split.by = "RNA_snn_res.1.5", reduction = 'harmonyumap')
dev.off()

Idents(all_data2) <- "RNA_snn_res.1.5"
all_data3 <- subset(all_data2, idents = c(0:28))

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_stardist3.pdf")), height = 6, width = 8)
DimPlot(all_data3, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data3, group.by = "RNA_snn_res.1", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data3, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data3, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(dr_outdir, str_glue("Dimplot dims{dims} {reduction.use} cluster_split_stardist3.pdf")), width = 100, height = 5)
DimPlot(all_data3, group.by = "RNA_snn_res.0.8", split.by = "RNA_snn_res.0.8", reduction = 'harmonyumap')
DimPlot(all_data3, group.by = "RNA_snn_res.1", split.by = "RNA_snn_res.1", reduction = 'harmonyumap')
DimPlot(all_data3, group.by = "RNA_snn_res.1.2", split.by = "RNA_snn_res.1.2", reduction = 'harmonyumap')
DimPlot(all_data3, group.by = "RNA_snn_res.1.5", split.by = "RNA_snn_res.1.5", reduction = 'harmonyumap')
dev.off()




### normalize data --------------------
# saveRDS(all_data, 'all_data_stardist.rds')
# saveRDS(all_data3, 'all_data3_stardist.rds')
all_data <- readRDS('all_data3_stardist.rds')
all_data <- NormalizeData(all_data, scale.factor = median(all_data[[ncount_name]] %>% unlist))
all_data <- FindVariableFeatures(all_data, nfeatures = 2000)
all_data <- ScaleData(all_data,vars.to.regress = ncount_name)

all_data <- RunPCA(all_data, features = VariableFeatures(all_data))

pdf(file = file.path(dr_outdir, 'pca plot_stardist2.pdf'), height = 6, width = 8)
ElbowPlot(all_data, ndims = 30)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 1-5_stardist2.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 1:5)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 6-10_stardist2.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 6:10)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 11-15_stardist2.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 11:15)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 16-20_stardist2.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 16:20)
dev.off()

dims <- 10

all_data <- RunUMAP(all_data, dims = 1:dims, return.model = TRUE)

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} project_stardist2.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "umap") + NoLegend()
dev.off()


all_data <- all_data %>%
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    reduction.save = "harmony"
  )

all_data <- RunUMAP(
  all_data,
  dims = 1:dims,
  return.model = TRUE,
  reduction = "harmony",
  reduction.name = "harmonyumap",
  reduction.key  = "harmonyUMAP_"
)

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} project - harmony_stardist2.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster ------
all_data <-
  FindNeighbors(
    all_data,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )

all_data <-
  FindClusters(
    all_data,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )


# dimplot
pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_stardist4.pdf")), height = 6, width = 8)
# DimPlot(all_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()




# dotplot ------------------------------------------------------------------------
major_marker_ls <- c(
  'Mast cells' = 'TPSAB1',
  'Mast cells' = 'CPA3',
  'Mast cells' = 'HPGDS',
  'Mast cells' = 'RGS13',
  'Mast cells' = 'IL1RL1',
  'Mast cells' = 'KIT',
  'Mast cells' = 'GATA2',
  'Mast cells' = 'FOSB',
  'Mast cells' = 'RP11-354E11.2',
  'Mast cells' = 'VWA5A',
  'Mast cells' = 'SLC18A2',
  'Mast cells' = 'HDC',
  'Mast cells' = 'CAPG',
  
  'Endothelials' = 'RAMP2',
  'Endothelials' = 'TM4SF1',
  'Endothelials' = 'RNASE1',
  'Endothelials' = 'EGFL7',
  'Endothelials' = 'RAMP3',
  'Endothelials' = 'PLVAP',
  'Endothelials' = 'AQP1',
  'Endothelials' = 'ECSCR',
  'Endothelials' = 'FKBP1A',
  'Endothelials' = 'AC011526.1',
  'Endothelials' = 'EMP1',
  'Endothelials' = 'DARC',
  'Endothelials' = 'VWF',
  'Endothelials' = 'EMCN',
  
  'Pericytes' = 'RGS5',
  'Pericytes' = 'ACTA2',
  #  'Pericytes' = 'MYH11',
  'Pericytes' = 'MT1M',
  'Pericytes' = 'FRZB',
  'Pericytes' = 'MT1A',
  'Pericytes' = 'NDUFA4L2',
  'Pericytes' = 'PPP1R14A',
  'Pericytes' = 'MYLK',
  'Pericytes' = 'PHLDA1',
  
  'Fibroblasts' = 'DCN',
  'Fibroblasts' = 'LUM',
  'Fibroblasts' = 'PTN',
  'Fibroblasts' = 'IGF1',
  'Fibroblasts' = 'APOD',
  'Fibroblasts' = 'COL1A2',
  'Fibroblasts' = 'FBLN1',
  'Fibroblasts' = 'MEG3',
  'Fibroblasts' = 'CXCL12',
  
  'B cells' = 'MS4A1',
  'B cells' = 'CD79A',
  'B cells' = 'BTG1',
  'B cells' = 'VPREB3',
  'B cells' = 'BANK1',
  'B cells' = 'IGLL5',
  'B cells' = 'CD79B',
  
  'Plasma cells' = 'SEC11C',	
  'Plasma cells' = 'XBP1',
  'Plasma cells' = 'PRDX4',
  'Plasma cells' = 'SPCS2',
  'Plasma cells' = 'SSR3',
  'Plasma cells' = 'SDF2L1',
  'Plasma cells' = 'C19orf10',
  'Plasma cells' = 'MANF',
  'Plasma cells' = 'TMEM258',
  'Plasma cells' = 'DNAJB9',
  
  'Macrophages' = 'C1QA',	
  'Macrophages' = 'C1QC',
  'Macrophages' = 'C1QB',
  'Macrophages' = 'GPR34',
  'Macrophages' = 'CSF1R',
  'Macrophages' = 'MS4A4A',
  
  'Monocytes' = 'S100A9',
  'Monocytes' = 'FCN1',
  'Monocytes' = 'S100A8',
  'Monocytes' = 'CSTA',
  'Monocytes' = 'EREG',
  'Monocytes' = 'NEAT1',
  'Monocytes' = 'TKT',
  'Monocytes' = 'THBS1',
  'Monocytes' = 'VCAN',
  'Monocytes' = 'TSPO',
  'Monocytes' = 'LYZ',
  
  'T cells' = 'NKG7',
  'T cells' = 'KLRC1',
  'T cells' = 'CCR7',
  'T cells' = 'FOXP3',
  'T cells' = 'CLTA4',
  'T cells' = 'CD8B',
  'T cells' = 'XCR6',
  'T cells' = 'CD3D',
  
  #  'Smooth muscle cells' = 'ACTA2',
  'Smooth muscle cells' = 'MYH11',
  'Smooth muscle cells' = 'TAGLN',
  'Smooth muscle cells' = 'CNN1',
  'Smooth muscle cells' = 'DES',
  'Smooth muscle cells' = 'ACTG2',
  
  'Club cells' = 'SCGB1A3',
  'Club cells' = 'WFDC2',
  'Club cells' = 'LCN2',
  'Club cells' = 'MMP7',	
  'Club cells' = 'KRT4',	
  'Club cells' = 'TACSTD2',
  'Club cells' = 'SCGB3A1',
  'Club cells' = 'LTF',
  
  'Hillock cells' = 'KRT13',
  'Hillock cells' = 'S100A16',
  'Hillock cells' = 'S100A14',
  'Hillock cells' = 'KRT19',
  
  'Basal cells' = 'TP63',
  'Basal cells' = 'KRT14',
  'Basal cells' = 'KRT5',
  'Basal cells' = 'TRIM29',
  'Basal cells' = 'SLC14A1',
  
  'NE cells' = 'ASCL-1',
  'NE cells' = 'FoxA2',
  'NE cells' = 'CHGA',
  'NE cells' = 'CHGB',
  'NE cells' = 'SYB',
  'NE cells' = 'NKX2-1',
  'NE cells' = 'ENO2',
  'NE cells' = 'MYCN',
  'NE cells' = 'POUF3F2',
  'NE cells' = 'LMO3',
  'NE cells' = 'INSM1',
  'NE cells' = 'SOX2',
  'NE cells' = 'EZH2',
  'NE cells' = 'SIAH2',
  
  'Luminal cells' = 'KRT8',
  'Luminal cells' = 'KRT18',
  'Luminal cells' = 'NKX3-1',
  'Luminal cells' = 'KLK4',
  'Luminal cells' = 'KLK3',
  'Luminal cells' = 'KLK2'
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

major_marker_ls <- major_marker_ls[major_marker_ls %in% rownames(all_data@assays$RNA$data)]
pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} main_marker_stardist.pdf")), width = 35, height = 6)
DotPlot(all_data, group.by = "RNA_snn_res.0.2", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(all_data, group.by = "RNA_snn_res.0.4", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(all_data, group.by = "RNA_snn_res.0.6", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(all_data, group.by = "RNA_snn_res.0.8", features = major_marker_ls, dot.scale = 10) + theme_90angle
dev.off()

pdf(
  file.path(dr_outdir, str_glue("featureplot dims{dims} {reduction.use} main_marker_stardist.pdf")), 
  width = 25, 
  height = 70
)
FeaturePlot(
  all_data,
  features = major_marker_ls,
  ncol = 5,
  reduction = "harmonyumap",
  raster = TRUE
)  
dev.off()




# Annotaiton -------------------------------------------------------------------------------
Idents(all_data) <- "RNA_snn_res.0.8"
cell_cluster <- c(
  '0' = 'NE cells', # plasma
  '1' = 'NE cells', # plasma
  '2' = 'NE cells', # plasma
  '3' = 'Luminal cells',
  '4' = 'Luminal cells',
  '5' = 'Plasma cells', #
  '6' = 'Plasma cells', 
  '7' = 'Luminal cells',
  '8' = 'Smooth muscle cells',#
  '9' = 'NE cells',#
  '10' = 'Fibroblasts',
  '11' = 'NE cells', # plasma
  '12' = 'Club cells',
  '13' = 'Smooth muscle cells',
  '14' = 'NE cells',#
  '15' = 'NE cells',#
  '16' = 'Luminal cells', #
  '17' = 'Luminal cells',
  '18' = 'NE cells', #
  '19' = 'NE cells',#
  '20' = 'Luminal cells'
)

all_data <- RenameIdents(all_data, cell_cluster)
all_data$annotation <- Idents(all_data)


pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} annotation {reduction.use}_stardist.pdf")), height = 6, width = 8)
# DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()




# Luminal subype ---------------------------------------------------------------
luminal <- subset(all_data, annotation == 'Luminal cells')
luminal_out_dir <- "Luminal_subtype" # dir.create(luminal_out_dir)

## PCA --------
luminal <- FindVariableFeatures(luminal)
luminal <- ScaleData(luminal, vars.to.regress = 'nCount_RNA')
luminal <- RunPCA(luminal)

### pca inspector ------
pdf(file.path(luminal_out_dir, 'pca plot_stardist.pdf'), height = 6, width = 8)
ElbowPlot(luminal, ndims = 50)
dev.off()

pdf(file.path(luminal_out_dir, 'dimheatmap 1-5_stardist.pdf'), height = 6, width = 10)
DimHeatmap(luminal, dims = 1:5)
dev.off()
pdf(file.path(luminal_out_dir, 'dimheatmap 6-10_stardist.pdf'), height = 6, width = 10)
DimHeatmap(luminal, dims = 6:10)
dev.off()
pdf(file.path(luminal_out_dir, 'dimheatmap 11-15_stardist.pdf'), height = 6, width = 10)
DimHeatmap(luminal, dims = 11:15)
dev.off()
pdf(file.path(luminal_out_dir, 'dimheatmap 16-20_stardist.pdf'), height = 6, width = 10)
DimHeatmap(luminal, dims = 16:20)
dev.off()

## dimension reduction --------
dims <- 10

luminal <- RunUMAP(luminal, dims = 1:dims,return.model = TRUE)

pdf(file.path(luminal_out_dir, str_glue("Dimplot dim{dims} project_stardist.pdf")), height = 6, width = 8)
DimPlot(luminal, group.by = "project", label = TRUE) + NoLegend()
dev.off()

luminal <- luminal %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

luminal <- RunUMAP(
  luminal, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(luminal_out_dir, str_glue("Dimplot dim{dims} project - harmony_stardist.pdf")), height = 6, width = 8)
DimPlot(luminal, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
luminal  <-
  FindNeighbors(
    luminal ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


luminal <-
  FindClusters(
    luminal,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(luminal_out_dir, str_glue("Dimplot dim{dims} cluster_stardist.pdf")), height = 6, width = 8)
DimPlot(luminal, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(luminal, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(luminal, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(luminal, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

luminal_marker_ls <- c(
  'Adeno-malignant luminal' = 'AMACR',
  'Adeno-malignant luminal' = 'CACNA1D',
  'Adeno-malignant luminal' = 'PCA3',
  'Adeno-malignant luminal' = 'ERG',
  'Adeno-malignant luminal' = 'FABP5',
  'Adeno-malignant luminal' = 'COL9A2',
  'Adeno-malignant luminal' = 'GCNT1',
  'Adeno-malignant luminal' = 'PHGR1',
  'Adeno-malignant luminal' = 'KLK4',
  'Adeno-malignant luminal' = 'KLK3',
  'Adeno-malignant luminal' = 'KLK2',
  'Adeno-malignant luminal' = 'AR',
  
  'Normal/Normal-like luminal' = 'ACPP',
  'Normal/Normal-like luminal' = 'MSMB',
  'Normal/Normal-like luminal' = 'RDH11',
  'Normal/Normal-like luminal' = 'AZGP1'
)

luminal_markers <-
  FindAllMarkers(
    # JoinLayers(luminal),
    luminal, 
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )

luminal_top5_markers <- 
  luminal_markers %>% 
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

names(luminal_top5_markers) <- sub('.*\\.','',names(luminal_top5_markers))

# dotplot
pdf(file.path(luminal_out_dir, str_glue("dotplot dims{dims} major_marker_stardist.pdf")), width = 25, height = 6)
DotPlot(luminal, group.by = "RNA_snn_res.0.6", features = luminal_marker_ls)+RotatedAxis()
DotPlot(luminal, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()
DotPlot(luminal, group.by = "RNA_snn_res.0.6", features = luminal_top5_markers)+RotatedAxis()
dev.off()


## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = 'Adeno-malignant luminal',
  '1'  = 'Adeno-malignant luminal',
  '2'  = 'Adeno-malignant luminal',
  '3'  = 'Adeno-malignant luminal',
  '4'  = 'Adeno-malignant luminal',
  '5'  = 'Adeno-malignant luminal',
  '6'  = 'Adeno-malignant luminal',
  '7'  = 'Adeno-malignant luminal',
  '8'  = 'Adeno-malignant luminal',
  '9'  = 'Adeno-malignant luminal',
  '10' = 'Adeno-malignant luminal',
  '11'  = 'Adeno-malignant luminal',
  '12'  = 'Adeno-malignant luminal',
  '13'  = 'Adeno-malignant luminal',
  '14'  = 'Adeno-malignant luminal',
  '15'  = 'Adeno-malignant luminal',
  '16'  = 'NE cells',
  '17'  = 'Adeno-malignant luminal'
)

Idents(luminal) <- "RNA_snn_res.0.6"
luminal <- RenameIdents(luminal, cell_cluster)
luminal$subtype <- Idents(luminal)


pdf(file.path(luminal_out_dir, str_glue("Dimplot dim{dims} annotation_stardist.pdf")), height = 6, width = 8)
DimPlot(luminal, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(luminal_out_dir, str_glue("Dimplot dim{dims} annotation split_stardist.pdf")), height = 7, width = 15)
DimPlot(luminal, group.by = "subtype",split.by = 'project', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

# saveRDS(luminal, file.path(luminal_out_dir, 'luminal_stardist.rds'))




# 映射-------------------------------
library(tibble)
library(data.table)
subtype <- luminal[['subtype']] %>% rownames_to_column('cells')

meta.data.t <-
  subtype%>%
  right_join(
    all_data[['annotation']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation2 = ifelse((as.character(annotation) %in% c('Luminal cells')),
                         as.character(subtype), 
                         as.character(annotation)))%>%
  dplyr::select(cells, annotation2) %>%
  column_to_rownames('cells')


all_data<-AddMetaData(all_data,meta.data.t['annotation2'])




# cell type dotplot ------------------------------------------------------------
luminal_cellTye_markers <- c(
  'Fibroblasts' = 'COL1A2',
  'Fibroblasts' = 'APOD',
  'Fibroblasts' = 'LUM',
  
  'Plasma cells' = 'SEC11C',	
  'Plasma cells' = 'PRDX4',
  'Plasma cells' = 'XBP1',
  'Plasma cells' = 'MANF',
  
  'Smooth muscle cells' = 'MYH11',
  'Smooth muscle cells' = 'TAGLN',
  'Smooth muscle cells' = 'CNN1',
  'Smooth muscle cells' = 'DES',
  'Smooth muscle cells' = 'ACTG2',
  
  'Club cells' = 'WFDC2',
  'Club cells' = 'MMP7',	
  'Club cells' = 'TACSTD2',
  'Club cells' = 'LTF',
  
  'NE cells' = 'CHGA',
  'NE cells' = 'CHGB',
  
  'Luminal cells' = 'KRT8',
  'Luminal cells' = 'NKX3-1',
  # 'Luminal cells' = 'KLK3',
  # 'Luminal cells' = 'KLK2',
  
  'Adeno-malignant luminal' = 'AMACR',
  'Adeno-malignant luminal' = 'KLK3',
  'Adeno-malignant luminal' = 'KLK2',
  'Normal/Normal-like luminal' = 'AZGP1'
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

pdf(file.path(luminal_out_dir, "dotplot celltype main_marker_stardist.pdf"), width = 20, height = 6)
DotPlot(all_data, group.by = "annotation2", features = luminal_cellTye_markers, dot.scale = 10) + theme_90angle
dev.off()


# saveRDS(all_data, 'all_data_celltype_stardist.rds')




# cellType distribution ---------------------------------------------------------------
all_data_meta <-
  all_data[[]] %>%
  rownames_to_column("cell_id")


p_col_proj <-
  all_data_meta %>%
  group_by(project, annotation2) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = project, y = freq, fill = annotation2)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


p_box_group  <-
  all_data_meta %>%
  group_by(project, annotation2) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = annotation2, y = freq, fill = project)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.25)) +
  ggsci::scale_fill_igv() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(dr_outdir, "cell_freq_project_stardist.pdf"),
  height = 6,
  width  = 6
)

ggsave(
  p_box_group,
  filename = file.path(dr_outdir, "cell_freq_group_stardist.pdf"),
  height = 6,
  width  = 10
)

p_point_cellType <- 
  all_data_meta %>%
  group_by(annotation2,project) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = annotation2, y = freq, fill = project)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_point_cellType,
  filename = file.path(dr_outdir, "group_freq_cellType_stardist.pdf"),
  height = 6,
  width  = 10
)




## spatial plot ----------------------------------------------------------------
polygon_info <- 
  all_data[[]] %>%
  rownames_to_column('cell_id') %>% 
  mutate(label_id = id) %>% 
  select(cell_id, geometry, label_id, project) %>% 
  mutate(
    geometry = geometry %>% 
      str_replace('.*\\((.*)\\){2}', '\\1') 
  ) %>% 
  separate_rows(geometry, sep = ',') %>%
  map(str_trim) %>%
  as.data.frame() %>% 
  separate(geometry, c('vertex_x', 'vertex_y'), sep = ' ') %>%
  mutate(
    vertex_x = vertex_x %>% as.numeric(), 
    vertex_y = vertex_y %>% as.numeric()
    )


ADPC_boundary <- subset(polygon_info, project == 'ADPC')[1:4]
NEPC_boundary <- subset(polygon_info, project == 'NEPC')[1:4]
cells_boundary_ls <- list(ADPC_boundary, NEPC_boundary)
names(cells_boundary_ls) <- c('ADPC', 'NEPC')

# polygon_info %>% 
#   mutate(
#     geometry = sub('.*\\((.*)\\){2}', '\\1', geometry)
#     ) %>% 
#   separate_rows(geometry, sep = ',') %>%
#   mutate(geometry = sub('^ ', '', geometry)) %>%
#   separate(geometry, c('vertex_x', 'vertex_y'), sep = ' ') %>%
#   head()

all_data_meta <- 
  all_data@meta.data %>%
  rownames_to_column('cell_id')

sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation2, group = cell_id))  +
      scale_y_reverse() +
      #coord_fixed(ratio = 1/image_wh_r) +
      coord_equal() +
      ggsci::scale_fill_igv() +
      theme_void()
    p
  })

pdf(file.path(sp_output, str_glue("Spatial_cell_type_stardist.pdf")), width = 20, height = 10 )
sp_p_ls %>% 
  imap(function(p, name){
    p1 <- p + ggtitle(name)
   # p1
  }) %>% patchwork::wrap_plots(ncol = 2)
dev.off()

