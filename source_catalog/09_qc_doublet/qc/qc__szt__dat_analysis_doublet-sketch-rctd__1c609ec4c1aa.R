library(Seurat)
library(tidyverse)

# setting dir -------
img_output <- "materials_for_rmd_v2_16" # dir.create(img_output)
qc_outdir <- "1. QC_v2_16" # dir.create(qc_outdir)
dr_outdir <- "2. dimension reduction and cluster_v2_16"  # dir.create(dr_outdir)
sp_output <- "3. spatial plot_v2_16" # dir.create(sp_output)
mk_output <- "4. marker gene_v2_16" # dir.create(mk_output)
er_output <- "5. enriment analysis_v2_16" # dir.create(er_output)


# set default ------------
orginism <- "mm" # "hsa" # "mm"
orginism_kegg <- "mmu" # "hsa" # "mmu"

mt_pattern   <- if(orginism == 'hsa') "^MT-" else "^mt-"
ribo_pattern <- if(orginism == 'hsa') "^RP[SL]" else "^Rp[sl]"

org.db <- if(orginism == 'hsa') "org.Hs.eg.db" else "org.Mm.eg.db"

bin_size <- 16 # 16
assay_name <- case_when(
  bin_size == 8  ~ 'Spatial.008um',
  bin_size == 16 ~ 'Spatial.016um',
  TRUE ~ 'wrong bin size'
)

nfeature_name <- str_glue("nFeature_{assay_name}")
ncount_name   <- str_glue("nCount_{assay_name}")

# 4 11 20 hvg5000
# 12 hvg2000 regression out H cell


# set future -------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


# read data -------
sample_df <- # read.csv("data/Annoation.csv")
  data.frame(
    sample = c('DZZ_Saline34', 'SYZ_NH4_44'),
    group  = c('DZZ_Saline', 'SYZ_NH4')
  )

file_list <- 
  list.files(
    "/home/ztt/scRNA/spaceranger_out/FWSP20240508/spaceranger_output",
    pattern = paste(sample_df$sample, collapse = "|"), 
    full.names = TRUE
  ) 

names(file_list) <- sample_df$group # basename(file_list)


visim_ls <-
  file_list %>% 
  imap(function(path, name){
    path <- file.path(path, "outs")
    Load10X_Spatial(path, bin.size = bin_size, slice = name)
  })

# add sample_id or sample_name
visim_ls <-
  visim_ls %>%
  imap(function(visium, name){
    visium$project <- name
    Idents(visium) <- 'project'
    RenameCells(visium, add.cell.id = name)
  })


visim_ls <- 
  visim_ls %>%
  map(~ {DefaultAssay(.x) <- assay_name; .x}) %>%
  map(NormalizeData) 


# init qc, used for doublet test --------------
visim_ls <- 
  visim_ls %>%
  map(function(visium){
    if(assay_name == "Spatial.008um"){
      subset(visium, nCount_Spatial.008um > 10)
    }else{
      subset(visium, nCount_Spatial.016um > 10)
    }
    
  })

visim_ls <- 
  visim_ls %>%
  map(function(x){# browser()
    x[["percent.mt"]]    <- PercentageFeatureSet(x, pattern = mt_pattern)
    x[["percent.ribo"]]  <- PercentageFeatureSet(x, pattern = ribo_pattern)
    x$log10GenesPerUMI   <- log10(x[[nfeature_name]][[1]] + 1) / log10(x[[ncount_name]][[1]] + 1)
    x
  })

# doublet test -------------------
library(scDblFinder)
library(BiocParallel)

bp <- MulticoreParam(3, RNGseed=1234)
# bpstart(bp)
# bpstop(bp)

visim_ls <-
  visim_ls %>%
  map(function(visium){
    sce <- as.SingleCellExperiment(visium)
    sce <- 
      scDblFinder(
        sce, 
        dbr      = 0.8,
        clusters = NULL, 
        samples  = NULL,
        BPPARAM  = bp # SerialParam()
      )
    meta_data <-
      sce %>%
      colData() %>%
      as.data.frame() %>%
      rownames_to_column("cell_id") %>%
      dplyr::select(cell_id, starts_with('scDblFinder')) %>%
      column_to_rownames("cell_id")
    
    visium %>%
      AddMetaData(metadata = meta_data)
  })


# qc -------------
p_qc_ls <-
  visim_ls %>%
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
  visim_ls %>%
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


# export
# saveRDS(p_qc_ls, file = file.path(rds_output, "1. p_qc_ls.rds"))

p_qc_ls %>%
  iwalk(function(p, name){
    png(file.path(img_output, str_glue("1 qc plot - {name}.png")), width = 480 * 3, height = 480 )
    print(p)
    dev.off()
  })


pdf(
  file.path(qc_outdir, "qc plot.pdf"),
  height = 6, 
  width = 20
)
p_qc_ls %>% walk(print)
dev.off()

pdf(
  file.path(qc_outdir, "qc plot - log.pdf"),
  height = 6, 
  width = 20
)
p_qc_log_ls %>% walk(print)
dev.off()



# qc filter ----------
qc_metric <- read.csv("qc_metric.csv") %>% map_dfc(str_trim)

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
  names(visim_ls) %>%
  set_names(.,.) %>%
  map(function(sample){ #browser()
    all_data <- visim_ls[[sample]]
    qc <- qc_metric_ls[[sample]]
    
    e <- str_glue("subset(
      all_data,
      {qc}
    )")
    
    eval(parse(text = e))
  })


# qc plot
qc_plot_filter_ls <-
  names(visim_ls) %>%
  setNames(., .) %>%
  map(function(sample){
    all_data <- visim_ls[[sample]] %>% subset(downsample = 100000, seed = 1)
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
    png(file.path(img_output, str_glue("1 qc plot - pass - {name}.png")), width = 360 * 3, height = 360 * 2 )
    print(p)
    dev.off()
  })


# export 
pdf(
  file.path(qc_outdir, "qc plot - pass.pdf"),
  height = 12, 
  width  = 18
)
qc_plot_filter_ls %>% walk(print)
dev.off()

qc_metric_ls %>%
  enframe("sample", "qc") %>%
  unnest(qc) %>% 
  openxlsx::write.xlsx(
    file = file.path(qc_outdir, "qc_metric.xlsx"),
    asTable = TRUE,
    withFilter = FALSE
  )

# cells filter
cell_counts_df <-
  names(visim_ls) %>%
  set_names(., .) %>%
  map(function(sample){
    all_data <- visim_ls[[sample]]
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
    file = file.path(qc_outdir, "cell count after filter.xlsx"),
    asTable = TRUE,
    withFilter = FALSE
  )


# dimension reduction and cluster and annotation ------
all_data <- merge(visim_ls_pass[[1]], visim_ls_pass[-1])

## normalize for sketch data ------
all_data <- NormalizeData(all_data)
all_data <- FindVariableFeatures(all_data)


## sketch data --------
all_data <-
  SketchData(
    all_data,
    ncells         = 50000L,
    sketched.assay = "sketch",
    method         = "LeverageScore",
    var.name       = "leverage.score",
    seed = 123L
  )

singlet_cells <- 
  visim_ls_pass %>%
  map(function(visium){
    visium[[]] %>%
      rownames_to_column('cell_id') %>%
      dplyr::filter(scDblFinder.class == 'singlet') %>%
      dplyr::pull(cell_id)
  }) %>%
  unlist(use.names = FALSE)

all_data[['sketch']] <- all_data[['sketch']] %>% subset(cells = singlet_cells)

### dimension reduction
DefaultAssay(all_data) <- "sketch"
all_data <- FindVariableFeatures(all_data, nfeatures = 3000)
all_data <- ScaleData(
  all_data,
  vars.to.regress = c(
    ncount_name
  )
)

# all_data <- SCTransform(all_data, assay = "sketch")

all_data <- RunPCA(all_data, features = VariableFeatures(all_data))

pdf(file = file.path(dr_outdir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_data, ndims = 30)
dev.off()


set.seed(1234)
suset_cells <- 
  all_data[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(all_data) ) %>%
  dplyr::pull("cell_id")

pdf(file.path(dr_outdir, 'pca dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 1:5)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 6:10)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 11:15)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 16:20)
dev.off()


## dimension reduction --------
dims <- 20

all_data <- RunUMAP(all_data, dims = 1:dims, return.model = TRUE, assay  = 'sketch')


pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
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

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


reduction.use <- 'pca' # 'harmony' 'pca'
## cluster ------
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
pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()

# DimPlot(all_data, group.by = "SCT_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



## dotplot -------
major_marker_ls <- c(
  neron = 'Rbfox3', # 'NEUN',
  neron = 'Tubb3',
  neron = 'Map2',
  neron = 'Gap43',
  neron = 'Pcp4',

  'Glutamatergic neurons' = 'Slc17a7',
  'GABAergic neurons' = 'Gad2',
  'GABAergic neurons' = 'Gad1',
  'Dopaminergic neurons' = 'Th',
  'Dopaminergic neurons' = 'Aldh1a1',
  'Serotonergic neurons' = 'Tph1',
  'Serotonergic neurons' = 'Tph2',
  'Cholinergic neurons' = 'Chat',

  Astro = 'Gfap',
  Astro = 'Aldh1l1',
  Astro = 'S100b',

  Oligo = 'Mbp',
  Oligo = 'Plp1',
  Oligo = 'Mobp',
  Oligo = 'Cnp',

  'Microglial cells' = 'Dock8',
  'Microglial cells' = 'Sfmbt2',
  'Microglial cells' = 'Cx3cr1',
  'Microglial cells' = 'C1qc',

  'Endothelial cells' =	'Flt1',
  'Endothelial cells' =	'Vwf',
  'Endothelial cells' =	'Pecam1',
  'Endothelial cells' =	'Cdh5',
  'T cells' =	'Cd96',
  'T cells' =	'Skap1',
  'T cells' =	'Cd3g',
  'T cells' =	'Cd3e',
  'T cells' =	'Cd3d',
  Macrophages	= 'Plxdc2',
  Macrophages	= 'Cybb',
  Macrophages	= 'Apbb1ip'
)

major_marker_ls <- c(
  'Slc17a7',
  'Sst',
  'Pvalb',
  'Id2',
  'Vip',
  L4 = 'Rorb',
  L2_3 = 'Fezf2',
  L5_6 = 'Themis'
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} main_marker.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "sketch_snn_res.0.2", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(all_data, group.by = "sketch_snn_res.0.4", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(all_data, group.by = "sketch_snn_res.0.6", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(all_data, group.by = "sketch_snn_res.0.8", features = major_marker_ls, dot.scale = 10) + theme_90angle


# DotPlot(all_data, group.by = "SCT_snn_res.0.2", dot.scale = 10, features = major_marker_ls) + theme_90angle
# DotPlot(all_data, group.by = "SCT_snn_res.0.4", dot.scale = 10, features = major_marker_ls) + theme_90angle
# DotPlot(all_data, group.by = "SCT_snn_res.0.6", dot.scale = 10, features = major_marker_ls) + theme_90angle
# DotPlot(all_data, group.by = "SCT_snn_res.0.8", dot.scale = 10, features = major_marker_ls) + theme_90angle
dev.off()





pdf(
  file.path(dr_outdir, str_glue("featureplot dims{dims} {reduction.use} main_marker.pdf")), 
  width = 25, 
  height = 50
)
FeaturePlot(
  all_data,
  features = major_marker_ls,
  ncol = 5,
  raster = TRUE
)  
dev.off()


## annotation ----
Idents(all_data) <- "sketch_snn_res.0.2"



cell_cluster <- c(
  '0'   = "Oligo",
  '1'   = "Gap43+ Neuron",
  '2'   = "Endothelial",
  '3'   = "Cnp+ Oligo",
  '4'   = "Gap43+ Oligo",
  '5'   = "Oligo",
  '6'   = "Cholinergic Neuron",
  '7'   = "Glutamatergic Neuron",
  '8'   = "Unknown",
  '9'   = "Gfap+ Oligo",
  '10'  = "Cnp-dim Oligo",
  '11'  = "Map2+ Astro"
)


all_data <- RenameIdents(all_data, cell_cluster)
all_data$annotation <- Idents(all_data)


pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} annotation {reduction.use} .pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

# pdf(str_glue(str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2.pdf")), width = 25, height = 6)
# DotPlot(all_data, group.by = "annotation", features = major_cell_marker)
# dev.off()


## project data ------
full_umap <- ifelse(reduction.use == "pca", "umap", "harmonyumap")
all_data <- 
  all_data %>% 
  ProjectData(
    assay = "Spatial.008um",
    full.reduction = "pca.full",
    
    sketched.assay = "sketch",
    sketched.reduction = reduction.use, # "harmony",
    
    umap.model = full_umap, 
    dims = 1:dims, 
    
    refdata = list(
      annotation_full = "annotation",
      sketch_snn_res.0.2_full = "sketch_snn_res.0.2"
    )
  )

Idents(all_data) <- 'annotation_full'
DefaultAssay(all_data) <- "Spatial.008um"

pdf(
  file.path(dr_outdir, str_glue("Dimplot dim{dims} annotation - {reduction.use} - projected.pdf")),
  height = 6,
  width = 8
)
DimPlot(all_data, group.by = "annotation_full", label = TRUE, reduction = str_glue('full.{full_umap}')) + NoLegend()
dev.off()

pdf(
  file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use} - projected.pdf")),
  height = 6,
  width = 8
)
DimPlot(all_data, group.by = "sketch_snn_res.0.2_full", label = TRUE, reduction = str_glue('full.{full_umap}')) + NoLegend()
dev.off()

# saveRDS(all_data, "all_data_projected.rds")
# all_data <- readRDS("all_data_projected.rds")

pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2 projected.pdf")), width = 25, height = 6)
DotPlot(all_data, assay = "sketch", group.by = "annotation_full", features = major_marker_ls) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
dev.off()

# export for rmd
png(file.path(img_output, str_glue("2 dimplot - cluster.png")), width = 600, height = 480)
DimPlot(
  all_data, 
  group.by = "sketch_snn_res.0.2_full", 
  label = TRUE, 
  reduction = full_umap
) +
  NoLegend()
dev.off()

## find markers -------
Idents(all_data) <- "sketch_snn_res.0.2"
DefaultAssay(all_data) <- "sketch"
all_markers <-
  FindAllMarkers(
    JoinLayers(all_data),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

# export marker gene
all_markers %>%
  dplyr::select(gene, everything()) %>%
  write.csv(
    file = file.path(mk_output, "marker.csv"), 
    row.names = FALSE
  )

all_de_genes <-
  all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  dplyr::filter(abs(avg_log2FC) > 1) %>% ##
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  pull("gene")

all_de_genes_for_enrichmenet <- 
  all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  dplyr::filter(abs(avg_log2FC) > 1) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")


all_top5_gene <- 
  all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  dplyr::filter(abs(avg_log2FC) > 1) %>% ##
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 5) %>%
  with(structure(
    gene,
    names = as.character(cluster)
  ))


# top5 dotplot
p_dot_top5 <-
  DotPlot(
    all_data, 
    features = all_top5_gene
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(
  p_dot_top5,
  filename = file.path(mk_output, "top 5 gene dotplot.pdf"),
  height = 6,
  width  = 12
)

# export for rmd
png(file.path(img_output, str_glue("4 top5 gene dotplot.png")), width = 960, height = 480)
print(p_dot_top5)
dev.off()


# top1 gene
all_top1_gene <- 
  all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  dplyr::filter(abs(avg_log2FC) > 1) %>% ##
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 1)%>%
  with(structure(
    gene,
    names = as.character(cluster)
  ))

# top1 gene umap
pdf(
  file.path(mk_output, str_glue("top 1 gene featureplot.pdf")), 
  width = 25, 
  height = 5 * ceiling(length(all_top1_gene) / 5)
)
FeaturePlot(
  all_data,
  features = all_top1_gene,
  reduction = "harmonyumap",
  ncol = 5,
  raster = TRUE
)  
dev.off()

# export for rmd
png(file.path(img_output, str_glue("4 top1 gene umap.png")), width = 400 * 5, height = 320 * ceiling(length(all_top1_gene) / 5))
FeaturePlot(
  all_data,
  features = all_top1_gene,
  reduction = "harmonyumap",
  ncol = 5,
  raster = TRUE
) 
dev.off()


# spatial feature plot
DefaultAssay(all_data) <- "Spatial.008um"
Images(all_data) %>%
  setNames(., str_remove(., "\\.[0-9]{3}um$")) %>%
  iwalk(function(img, img_name){
    
    png(file.path(img_output, str_glue("4 top1 gene spatial - {img_name}.png")), width = 400 * 5, height = 320 * ceiling(length(all_top1_gene) / 5))
    print(SpatialFeaturePlot(
      all_data,
      features = all_top1_gene,
      images = img,
      ncol = 5
    ))
    dev.off()
  })




# DoHeatmap(sc, features = all_de_genes, group.by = "RNA_snn_res.0.8")
if(!exists('avg_expr')) avg_expr <- AverageExpression(JoinLayers(all_data))[['sketch']]


all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  nest() %>%
  deframe() %>%
  map(pull, 'gene') %>%
  imap(~{
    xx <- str_c(.x, collapse = ",")
    str_glue("cluster_{.y}:{xx}")
  }) %>%
  walk(cat, "\n")

# avg_expr[all_de_genes, ] %>% 
#   as.matrix() %>% 
#   as.data.frame() %>% 
#   setNames(., paste('cluster', (1:ncol(.) - 1) )) %>%  t %>%
#   pheatmap::pheatmap(
#     scale = "column",
#     cluster_row = F, # avg_expr[VariableFeatures(all_data), ] %>% t %>% dist %>% hclust,
#     cluster_col = FALSE,
#     color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
#     border_col = NA
#   )


# spatial plot --------------
p_sp_plot_anno <-
  Images(all_data) %>%
  setNames(., str_remove(., "\\.[0-9]{3}um$")) %>%
  imap(function(img, img_name){
    p <- SpatialDimPlot(
      all_data,
      images = img, 
      group.by = "annotation_full",
      pt.size.factor	= 4
    ) +
      ggsci::scale_fill_aaas() +
      # ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      ggtitle(img_name)
    p
  })

p_sp_plot_cluster <-
  Images(all_data) %>%
  setNames(., str_remove(., "\\.[0-9]{3}um$")) %>%
  imap(function(img, img_name){
    p <- SpatialDimPlot(
      all_data,
      images = img, 
      group.by = "sketch_snn_res.0.2_full",
      pt.size.factor	= 4
    ) +
      # ggsci::scale_fill_aaas() +
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      ggtitle(img_name)
    p
  })

# export for rmd
p_sp_plot_anno %>%
  iwalk(function(p, name){
    png(file.path(img_output, str_glue("3 spatial plot - annotation - {name}.png")), width = 720, height = 720)
    print(p)
    dev.off()
  })

p_sp_plot_cluster %>%
  iwalk(function(p, name){
    png(file.path(img_output, str_glue("3 spatial plot - cluster - {name}.png")), width = 720, height = 720)
    print(p)
    dev.off()
  })

# export
pdf(
  file.path(sp_output, "spatial plot annotation.pdf"),
  height = 8, 
  width  = 8
)
p_sp_plot_anno %>% walk(print)
dev.off()

pdf(
  file.path(sp_output, "spatial plot cluster.pdf"),
  height = 8, 
  width  = 8
)
p_sp_plot_cluster %>% walk(print)
dev.off()

# enrichemt analysis -----
library(clusterProfiler)

all_de_genes_id_for_enrichmenet <-
  all_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::bitr(
      x,
      OrgDb = org.db,
      fromType = "SYMBOL",
      toType = "ENTREZID"
    ) %>%
      .$ENTREZID
  }) %>%
  discard(function(x) length(x) < 1)


enrich_ls <-
  all_de_genes_id_for_enrichmenet %>%
  map(function(x){
    enrichKEGG(
      x,
      organism = orginism_kegg, 
      keyType = "ncbi-geneid", 
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
    )
  })

# saveRDS(enrich_ls, 'enrich_ls.rds')

# expot enrichment table
enrich_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    cluster_name <- paste0("cluster_", name)
    dir.create(file.path(er_output, cluster_name))
    
    write.csv(
      x@result,
      file = file.path(er_output, cluster_name, "enriment table.csv"),
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

# export dotplot
enrich_dot_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    cluster_name <- paste0("cluster_", name)
    dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, cluster_name, "dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

# export for rmd
enrich_dot_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    cluster_name <- paste0("cluster_", name)
    
    png(file.path(img_output, str_glue("5 enrichment dotplot - {cluster_name}.png")), width = 720, height = 720)
    print(x)
    dev.off()
  })



# patch code -----------------------
dr_outdir2 <- str_glue("{dr_outdir} - 2") 
dir.create(dr_outdir2)

all_files_dr_output <- list.files(dr_outdir, full.names = TRUE)

# dimplot
file.copy(
  all_files_dr_output %>% .[str_detect(basename(.), str_glue("Dimplot dim{dims}"))],
  dr_outdir2
)
# dotplot
file.copy(
  all_files_dr_output %>% .[str_detect(basename(.), str_glue("dotplot dims{dims}"))],
  dr_outdir2
)
# featureplot
file.copy(
  all_files_dr_output %>% .[str_detect(basename(.), str_glue("featureplot dims{dims}"))],
  dr_outdir2
)
# pca
file.copy(
  all_files_dr_output %>% .[str_detect(basename(.), "^pca")],
  dr_outdir2
)

