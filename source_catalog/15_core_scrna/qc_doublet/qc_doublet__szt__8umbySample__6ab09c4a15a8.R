library(Seurat)
library(tidyverse)

# setting -------------------------------------------------------------
## set dir --------------
img_output <- "result_8um/materials_for_rmd" # dir.create(img_output)
qc_outdir <- "result_8um/1. QC" # dir.create(qc_outdir)
dr_outdir <- "result_8um/2. dimension reduction and cluster"  # dir.create(dr_outdir)
sp_output <- "result_8um/3. spatial plot" # dir.create(sp_output)
mk_output <- "result_8um/4. marker gene" # dir.create(mk_output)
er_output <- "result_8um/5. enriment analysis" # dir.create(er_output)


## set default ------------
orginism <- "hsa" # "hsa" # "mm"
orginism_kegg <- "hsa" # "hsa" # "mmu"

mt_pattern   <- if(orginism == 'hsa') "^MT-" else "^mt-"
ribo_pattern <- if(orginism == 'hsa') "^RP[SL]" else "^Rp[sl]"
hb_genes <- if(orginism == 'hsa') c("HBB", "HBA2") else c("Hbb", "Hba2")

org.db <- if(orginism == 'hsa') "org.Hs.eg.db" else "org.Mm.eg.db"

bin_size <- 8
assay_name <- case_when(
  bin_size == 8  ~ 'Spatial.008um',
  bin_size == 16 ~ 'Spatial.016um',
  TRUE ~ 'wrong bin size'
)
nfeature_name <- str_glue("nFeature_{assay_name}")
ncount_name   <- str_glue("nCount_{assay_name}")

# 4 11 20 hvg5000
# 12 hvg2000 regression out H cell


## set future -------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 20 * 1024 ^ 3) # 20Gb


# read data --------------------------------------------------------------------
sample_df <- # read.csv("data/Annoation.csv")
  data.frame(
    sample = c('ADPC', 'NEPC'),
    group  = c('ADPC', 'NEPC')
  )

file_list <- 
  list.files(
    "/mnt/NexData_file/nas/outputdata/spaceranger_out/FWSP20202434_new", 
    pattern = paste(sample_df$sample, collapse = "|"), 
    full.names = TRUE
  ) 

names(file_list) <- sample_df$group # basename(file_list)

# visim_ls <-
#   file_list %>% 
#   imap(function(path, name){
#     path <- file.path(path, "outs")
#     Load10X_Spatial(path, bin.size = bin_size, slice = name)
#   })

visim_ls <-
  file_list %>% 
  imap(function(path, name){
    path <- file.path(path, "outs")
    
    Read10X_Image <-
      function (image.dir, spatial.dir = image.dir, image.name = "tissue_hires_image.png", 
                assay = "Spatial", slice = "slice1", filter.matrix = TRUE) {
        
        image <- png::readPNG(source = file.path(image.dir, image.name))
        scale.factors <- Read10X_ScaleFactors(filename = file.path(spatial.dir, 
                                                                   "scalefactors_json.json"))
        coordinates <- Read10X_Coordinates(filename = Sys.glob(file.path(spatial.dir, 
                                                                         "*tissue_positions*parquet")), filter.matrix)
        fov <- CreateFOV(coordinates[, c("imagerow", "imagecol")], 
                         type = "centroids", radius = scale.factors[["spot"]], 
                         assay = assay, key = Key(slice, quiet = TRUE))
        visium.fov <- new(Class = "VisiumV2", boundaries = fov@boundaries, 
                          molecules = fov@molecules, assay = fov@assay, key = fov@key, 
                          image = image, scale.factors = scale.factors)
        return(visium.fov)
      }
    
    Load10X_Spatial(
      path, 
      bin.size = bin_size, 
      slice = name, 
      image = Read10X_Image(
        file.path(path, "spatial"),
        spatial.dir = file.path(path, str_glue("binned_outputs/square_00{bin_size}um/spatial")),
        assay = nfeature_name
      )
    )
  })



## add sample_id or sample_name -----------------------
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


# qc ---------------------------------------------------------------------------
visim_ls <- 
  visim_ls %>%
  map(function(visium){
    if(assay_name == "Spatial.008um"){
      subset(visium, nCount_Spatial.008um > 0)
    }else{
      subset(visium, nCount_Spatial.016um > 0)
    }
  })

visim_ls <- 
  visim_ls %>%
  map(function(x){# browser()
    x[["percent.mt"]]    <- PercentageFeatureSet(x, pattern = mt_pattern)
    x[["percent.ribo"]]  <- PercentageFeatureSet(x, pattern = ribo_pattern)
    x[["percent.hb"]]  <- PercentageFeatureSet(x, features = hb_genes)
    x$log10GenesPerUMI   <- log10(x[[nfeature_name]][[1]] + 1) / log10(x[[ncount_name]][[1]] + 1)
    x
  })


p_qc_ls <-
  visim_ls %>%
  imap(function(x, name){ # browser()
    # only subset 10,000 cells
    x_sub <- subset(x, cells = Cells(x[[assay_name]]), downsample = 100000, seed = 1)
    p <- VlnPlot(
      x_sub, 
      features = c(nfeature_name, ncount_name, "percent.mt", "percent.hb", "log10GenesPerUMI"), 
      ncol     = 5, 
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
      features = c(nfeature_name, ncount_name, "percent.mt", "percent.hb", "log10GenesPerUMI"), 
      ncol     = 5, 
      pt.size  = 0,
      log      = TRUE,
      raster   = TRUE
    ) +
      patchwork::plot_annotation(title = name)
    
  })

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

saveRDS(visim_ls, 'result_8um/visim_ls.rds')




# NE cells qc -------------------------------------------------------------------
ADPC <- visim_ls$ADPC
NE_roi <- read.csv('result_8um/NE region_8um.csv')
NE_roi$Barcode <- str_replace(NE_roi$Barcode, '^', 'ADPC_')
ADPC_NE_region <- subset(ADPC, cells = NE_roi$Barcode)
p <- VlnPlot(
       ADPC_NE_region, 
       features = c(nfeature_name, ncount_name, "percent.mt", "percent.hb", "log10GenesPerUMI"), 
       ncol     = 5, 
       pt.size  = 0,
       log      = FALSE, 
       raster   = TRUE
   )
ggsave(file.path(qc_outdir, "qc plot_ADPC_NE_region.pdf"), height = 6, width = 20)

p1 <- VlnPlot(
  ADPC_NE_region, 
  features = c(nfeature_name, ncount_name, "percent.mt", "percent.hb", "log10GenesPerUMI"), 
  ncol     = 5, 
  pt.size  = 0,
  log      = TRUE, 
  raster   = TRUE
)
ggsave(file.path(qc_outdir, "qc plot - log_ADPC_NE_region.pdf"), height = 6, width = 20)




# ADPC --------------------------------------------------------------------------
## qc filter --------------------------------------------------------------------
qc_metric <- read.csv("result_8um/ADPC_qc_metric.csv") %>% map_dfc(str_trim)

qc_metric_variable <- c(
  nfeature_name, ncount_name,
  "percent.mt", "percent.hb", "log10GenesPerUMI"
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
  names(visim_ls[1]) %>%
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
  names(visim_ls[1]) %>%
  setNames(., .) %>%
  map(function(sample){
    all_data <- visim_ls[[sample]] %>% subset(downsample = 100000, seed = 1)
    all_data_pass <- visim_ls_pass[[sample]] %>% subset(downsample = 100000, seed  = 1)
    
    # all_data
    plot1 <- FeatureScatter(all_data, feature1 = ncount_name, feature2 = nfeature_name) +
      patchwork::plot_annotation(title = sample)
    plot2 <- FeatureScatter(all_data, feature1 = ncount_name, feature2 = "percent.mt") +
      patchwork::plot_annotation(title = sample)
    plot3 <- FeatureScatter(all_data, feature1 = ncount_name, feature2 = "percent.hb") +
      patchwork::plot_annotation(title = sample)
    plot4 <- FeatureScatter(all_data, feature1 = ncount_name, feature2 = "log10GenesPerUMI") +
      patchwork::plot_annotation(title = sample)
    p1 <- plot1 + plot2 + plot3 + plot4 + patchwork::plot_layout(ncol = 4)
    
    # all_data_pass
    plot1_pass <- FeatureScatter(all_data_pass, feature1 = ncount_name, feature2 = nfeature_name)+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    plot2_pass <- FeatureScatter(all_data_pass, feature1 = ncount_name, feature2 = "percent.mt")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    plot3_pass <- FeatureScatter(all_data_pass, feature1 = ncount_name, feature2 = "percent.hb")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    plot4_pass <- FeatureScatter(all_data_pass, feature1 = ncount_name, feature2 = "log10GenesPerUMI")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    p2 <- plot1_pass + plot2_pass + plot3_pass + plot4_pass + patchwork::plot_layout(ncol = 4)
    
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
  file.path(qc_outdir, "qc plot - pass_ADPC.pdf"),
  height = 12, 
  width  = 18
)
qc_plot_filter_ls %>% walk(print)
dev.off()

qc_metric_ls %>%
  enframe("sample", "qc") %>%
  unnest(qc) %>% 
  openxlsx::write.xlsx(
    file = file.path(qc_outdir, "qc_metric_ADPC.xlsx"),
    asTable = TRUE,
    withFilter = FALSE
  )

# cells filter
cell_counts_df <-
  names(visim_ls[1]) %>%
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
    file = file.path(qc_outdir, "cell count after filter_ADPC.xlsx"),
    asTable = TRUE,
    withFilter = FALSE
  )
# saveRDS(visim_ls_pass$ADPC, 'result_8um/ADPC_pass.rds')

SpatialDimPlot(visim_ls_pass$ADPC, group.by = 'project')




## dimension reduction and clustering -------------------------------------------
ADPC <- visim_ls_pass$ADPC
#SpatialFeaturePlot(ADPC, features = c('HBA2', 'HBB', ncount_name))

### normalize -------------------------------------------------------------------
ADPC <- NormalizeData(ADPC, scale.factor = median(ADPC[[ncount_name]] %>% unlist))
ADPC <- FindVariableFeatures(ADPC, nfeatures = 3000)
ADPC <- ScaleData(ADPC,vars.to.regress = ncount_name)

ADPC <- RunPCA(ADPC, features = VariableFeatures(ADPC))

pdf(file = file.path(dr_outdir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(ADPC, ndims = 30)
dev.off()


dims <- 10


#undebug(RunUMAP)


ADPC <- RunUMAP(ADPC, dims = 1:dims, return.model = TRUE)

DimPlot(ADPC, label = TRUE, reduction = "umap") + NoLegend()

reduction.use <- 'pca' # 'harmony' 'pca'




### cluster ---------------------------------------------------------------------
ADPC <-
  FindNeighbors(
    ADPC,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )

ADPC <-
  FindClusters(
    ADPC,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )


# dimplot
pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(ADPC, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(ADPC, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(ADPC, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(ADPC, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
dev.off()

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_split.pdf")), height = 6, width = 30)
DimPlot(ADPC, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.2') + NoLegend()
DimPlot(ADPC, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.4') + NoLegend()
DimPlot(ADPC, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.6') + NoLegend()
DimPlot(ADPC, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.8') + NoLegend()
dev.off()




## refiltering -----------------------------------------------------------------
Idents(ADPC) <- "Spatial.008um_snn_res.0.2"
ADPC2 <- subset(ADPC, idents = c(0:5))
pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_split2.pdf")), height = 6, width = 30)
DimPlot(ADPC2, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.2') + NoLegend()
DimPlot(ADPC2, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.4') + NoLegend()
DimPlot(ADPC2, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.6') + NoLegend()
DimPlot(ADPC2, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.8') + NoLegend()
dev.off()


Idents(ADPC2) <- "Spatial.008um_snn_res.0.6"
ADPC3 <- subset(ADPC2, idents = c(0:14))
pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_split3.pdf")), height = 6, width = 30)
DimPlot(ADPC3, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.2') + NoLegend()
DimPlot(ADPC3, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.4') + NoLegend()
DimPlot(ADPC3, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.6') + NoLegend()
DimPlot(ADPC3, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.8') + NoLegend()
dev.off()


Idents(ADPC3) <- "Spatial.008um_snn_res.0.8"
ADPC4 <- subset(ADPC3, idents = c(0:19))
pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_split4.pdf")), height = 6, width = 30)
DimPlot(ADPC4, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.2') + NoLegend()
DimPlot(ADPC4, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.4') + NoLegend()
DimPlot(ADPC4, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.6') + NoLegend()
DimPlot(ADPC4, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.8') + NoLegend()
dev.off()

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use} 4.pdf")), height = 6, width = 8)
DimPlot(ADPC4, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(ADPC4, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(ADPC4, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(ADPC4, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
dev.off()


## dimension reduction and clustering -------------------------------------------
ADPC4 <- RunPCA(ADPC4, features = VariableFeatures(ADPC4))

pdf(file = file.path(dr_outdir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(ADPC4, ndims = 30)
dev.off()

dims <- 10


ADPC4 <- RunUMAP(ADPC4, dims = 1:dims, return.model = TRUE)

DimPlot(ADPC4, label = TRUE, reduction = "umap") + NoLegend()

reduction.use <- 'pca' # 'harmony' 'pca'




### cluster ---------------------------------------------------------------------
ADPC4 <-
  FindNeighbors(
    ADPC4,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )

ADPC4 <-
  FindClusters(
    ADPC4,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )


# dimplot
pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use} 4_2.pdf")), height = 6, width = 8)
DimPlot(ADPC4, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(ADPC4, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(ADPC4, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(ADPC4, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
dev.off()

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_split 4_2.pdf")), height = 6, width = 30)
DimPlot(ADPC4, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.2') + NoLegend()
DimPlot(ADPC4, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.4') + NoLegend()
DimPlot(ADPC4, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.6') + NoLegend()
DimPlot(ADPC4, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.8') + NoLegend()
dev.off()




## dotplot ------------------------------------------------------------------------
major_marker_ls <- c(
  'Epithelial cells'  = 'EPCAM',
  'Epithelial cells'  = 'KRT5',
  'Epithelial cells'  = 'KRT14',
  
  'Club cells' = 'SCGB1A3',
  #  'Club cells' = 'WFDC2',
  'Club cells' = 'LCN2',
  'Club cells' = 'MMP7',	
  'Club cells' = 'KRT4',	
  'Club cells' = 'TACSTD2',
  'Club cells' = 'SCGB3A1',
  'Club cells' = 'LTF',
  'Club cells' = 'RARRES1',
  'Club cells' = 'KRT7',
  'Club cells' = 'AGR2',
  
  'Hillock cells' = 'KRT13',
  'Hillock cells' = 'S100A16',
  'Hillock cells' = 'S100A14',
  'Hillock cells' = 'KRT19',
  'Hillock cells' = 'WFDC2',
  
  'Basal cells' = 'TP63',
  #  'Basal cells' = 'KRT14',
  #  'Basal cells' = 'KRT5',
  'Basal cells' = 'TRIM29',
  'Basal cells' = 'SLC14A1',
  'Basal cells' = 'KRT17',
  'Basal cells' = 'KRT15',
  
  'Luminal cells' = 'KRT8',
  'Luminal cells' = 'KRT18',
  'Luminal cells' = 'NKX3-1',
  'Luminal cells' = 'KLK4',
  'Luminal cells' = 'KLK3',
  'Luminal cells' = 'KLK2',
  'Luminal cells' = 'MSMB',
  'Luminal cells' = 'ACPP',
  
  'NE cells' = 'ASCL-1',
  'NE cells' = 'FoxA2',
  'NE cells' = 'CHGA',
  'NE cells' = 'CHGB',
  'NE cells' = 'SYP',
  'NE cells' = 'NKX2-1',
  'NE cells' = 'ENO2',
  'NE cells' = 'MYCN',
  'NE cells' = 'POUF3F2',
  'NE cells' = 'LMO3',
  'NE cells' = 'INSM1',
  'NE cells' = 'SOX2',
  'NE cells' = 'EZH2',
  'NE cells' = 'SIAH2',
  
  'Fibroblasts' = 'DCN',
  'Fibroblasts' = 'LUM',
  'Fibroblasts' = 'VIM',
  'Fibroblasts' = 'IGFBP7',
  'Fibroblasts' = 'COL1A1',
  'Fibroblasts' = 'COL1A2',
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "SFRP4",
  "Fibroblasts" = "IGF1",
  "Fibroblasts" = "APOD",
  "Fibroblasts" = "FBLN1",
  "Fibroblasts" = "MEG3",
  "Fibroblasts" = "CXCL12",
  
  'B cells' = 'MS4A1',
  'B cells' = 'CD79A',
  'B cells' = 'CD19',
  'B cells' = 'BTG1',
  'B cells' = 'VPREB3',
  'B cells' = 'BANK1',
  'B cells' = 'IGLL5',
  'B cells' = 'CD79B',
  'B cells' = 'IGHM',
  'B cells' = 'CD22',
  
  "Plasma cells" = "IGHG1",
  "Plasma cells" = "IGHA1",
  "Plasma cells" = "IGHG3",
  "Plasma cells" = "IGKC",
  "Plasma cells" = "TNFRSF17",
  "Plasma cells" = "JCHAIN",
  "Plasma cells" = "MZB1",
  "Plasma cells" = "XBP1",
  'Plasma cells' = 'SEC11C',	
  'Plasma cells' = 'PRDX4',
  'Plasma cells' = 'SPCS2',
  'Plasma cells' = 'SSR3',
  'Plasma cells' = 'SDF2L1',
  'Plasma cells' = 'C19orf10',
  'Plasma cells' = 'MANF',
  'Plasma cells' = 'TMEM258',
  'Plasma cells' = 'DNAJB9',
  'Plasma cells' = 'CD138',
  'Plasma cells' = 'CD38',
  'Plasma cells' = 'BCMA',
  'Plasma cells' = 'IRF4',
  'Plasma cells' = 'BLIMP1',
  
  'T cells' = 'CD3D', 
  'T cells' = 'CD3E', 
  'T cells' = 'CD3G', 
  'T cells' = 'TRAC', 
  'T cells' = 'NKG7',
  'T cells' = 'KLRC1',
  'T cells' = 'CCR7',
  'T cells' = 'FOXP3',
  'T cells' = 'CLTA4',
  'T cells' = 'CD8B',
  'T cells' = 'XCR6',
  'T cells' = 'IL32',
  'T cells' = 'CXCR6',
  
  "Macrophages" = "CD68",
  "Macrophages" = "CD163",
  "Macrophages" = "C1QA",
  "Macrophages" = "C1QC", 
  'Macrophages'  = 'CSF1R',
  'Macrophages' = 'C1QB',
  'Macrophages' = 'GPR34',
  'Macrophages' = 'MS4A4A',
  
  "NK cells" = "GNLY",
  #  "NK cells" = "NKG7",
  "NK cells" = "NCAM1",
  "NK cells" = "GZMB",
  "NK cells" = "KLRD1",
  
  "Monocytes" = "CD14", 
  "Monocytes" = "FCGR3A",
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
  
  "Dendritic cell" = "CD83",
  'Dendritic cell' = "LILRA4",  
  'Dendritic cell' = "CLEC9A",
  'Dendritic cell' = "IDO1",
  'Dendritic cell' = "CD1C",
  'Dendritic cell' = "FCER1A",
  
  'Netrophils' = "CSF3R",
  'Netrophils' = "G0S2",
  'Netrophils' = "RETNLG",
  
  'Mast cells' = 'KIT',
  'Mast cells'= 'TPSAB1',
  'Mast cells' = 'CPA3',
  'Mast cells' = 'HPGDS',
  'Mast cells' = 'RGS13',
  'Mast cells' = 'IL1RL1',
  'Mast cells' = 'GATA2',
  'Mast cells' = 'FOSB',
  'Mast cells' = 'RP11-354E11.2',
  'Mast cells' = 'VWA5A',
  'Mast cells' = 'SLC18A2',
  'Mast cells' = 'HDC',
  'Mast cells' = 'CAPG',
  'Mast cells' = 'FCER1G',
  
  "Pericytes" = "PDGFRB",
  "Pericytes" = "RGS5",
  "Pericytes" = "MCAM",
  'Pericytes' = 'ACTA2',
  #  'Pericytes' = 'MYH11',
  'Pericytes' = 'MT1M',
  'Pericytes' = 'FRZB',
  'Pericytes' = 'MT1A',
  'Pericytes' = 'NDUFA4L2',
  'Pericytes' = 'PPP1R14A',
  'Pericytes' = 'MYLK',
  'Pericytes' = 'PHLDA1',
  
  #  'Smooth muscle cells' = 'ACTA2',
  'Smooth muscle cells' = 'MYH11',
  'Smooth muscle cells' = 'TAGLN',
  'Smooth muscle cells' = 'CNN1',
  'Smooth muscle cells' = 'DES',
  'Smooth muscle cells' = 'ACTG2',
  
  'Endothelials' = 'PECAM1', 
  'Endothelials' = 'ACKR1', 
  'Endothelials' = 'CLDN5',
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
  'Endothelials' = 'EMCN'
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

major_marker_ls <- major_marker_ls[major_marker_ls %in% rownames(LayerData(JoinLayers(ADPC4), layer = 'data'))]

pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} multi res_ADPC4.pdf")), width = 55, height = 8)
DotPlot(ADPC4, group.by = "Spatial.008um_snn_res.0.2", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(ADPC4, group.by = "Spatial.008um_snn_res.0.4", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(ADPC4, group.by = "Spatial.008um_snn_res.0.6", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(ADPC4, group.by = "Spatial.008um_snn_res.0.8", features = major_marker_ls, dot.scale = 10) + theme_90angle
dev.off()

## refiltering -----------------------------------------------------------------
Idents(ADPC4) <- "Spatial.008um_snn_res.0.6"
ADPC5 <- subset(ADPC4, idents = c(0:11, 13, 14))
pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_split5.pdf")), height = 6, width = 30)
DimPlot(ADPC5, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.2') + NoLegend()
DimPlot(ADPC5, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.4') + NoLegend()
DimPlot(ADPC5, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.6') + NoLegend()
DimPlot(ADPC5, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.8') + NoLegend()
dev.off()


Idents(ADPC5) <- "Spatial.008um_snn_res.0.8"
ADPC6 <- subset(ADPC5, idents = c(0:13, 16))
pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_split6.pdf")), height = 6, width = 30)
DimPlot(ADPC6, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.2') + NoLegend()
DimPlot(ADPC6, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.4') + NoLegend()
DimPlot(ADPC6, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.6') + NoLegend()
DimPlot(ADPC6, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.8') + NoLegend()
dev.off()


## dimension reduction and clustering -------------------------------------------
ADPC6 <- RunPCA(ADPC6, features = VariableFeatures(ADPC6))

pdf(file = file.path(dr_outdir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(ADPC6, ndims = 30)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(ADPC6, dims = 1:5)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(ADPC6, dims = 6:10)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(ADPC6, dims = 11:15)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(ADPC6, dims = 16:20)
dev.off()

dims <- 10


ADPC6 <- RunUMAP(ADPC6, dims = 1:dims, return.model = TRUE)

DimPlot(ADPC6, label = TRUE, reduction = "umap") + NoLegend()

reduction.use <- 'pca' # 'harmony' 'pca'




### cluster ---------------------------------------------------------------------
ADPC6 <-
  FindNeighbors(
    ADPC6,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )

ADPC6 <-
  FindClusters(
    ADPC6,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )


# dimplot
pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use} 6.pdf")), height = 6, width = 8)
DimPlot(ADPC6, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(ADPC6, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(ADPC6, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(ADPC6, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
dev.off()

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_split 6.pdf")), height = 6, width = 30)
DimPlot(ADPC6, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.2') + NoLegend()
DimPlot(ADPC6, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.4') + NoLegend()
DimPlot(ADPC6, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.6') + NoLegend()
DimPlot(ADPC6, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.8') + NoLegend()
dev.off()

pdf(
  file.path(dr_outdir, str_glue("featureplot dims{dims} {reduction.use} major_marker.pdf")), 
  width = 25, 
  height = 100
)
FeaturePlot(
  ADPC6,
  features = major_marker_ls,
  ncol = 5,
  reduction = 'umap',
  raster = TRUE
)  
dev.off()


Idents(ADPC6) <- "Spatial.008um_snn_res.0.8"
all_markers <-
  FindAllMarkers(
    JoinLayers(ADPC6),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )

all_top5_markers <- 
  all_markers %>% 
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

names(all_top5_markers) <- sub('.*\\.','',names(all_top5_markers))

# dotplot
pdf(file.path(dr_outdir, str_glue("dotplot dims {dims} {reduction.use} majorMarker_topDEGs.pdf")), width = 40, height = 6)
DotPlot(ADPC6, group.by = "Spatial.008um_snn_res.0.8", features = major_marker_ls)+RotatedAxis()
DotPlot(ADPC6, group.by = "Spatial.008um_snn_res.0.8", features = all_top5_markers)+RotatedAxis()
dev.off()

## Annotaiton -------------------------------------------------------------------------------
Idents(ADPC6) <- "Spatial.008um_snn_res.0.8"
cell_cluster <- c(
  '0' = 'Luminal cells', #Crisp3 genecard
  '1' = 'Luminal cells', 
  '2' = 'Luminal cells', 
  '3' = 'Luminal cells', 
  '4' = 'Fibroblasts', 
  '5' = 'Luminal cells', 
  '6' = 'Luminal cells', #Crisp3 genecard
  '7' = 'Luminal cells', 
  '8' = 'Luminal cells', #Crisp3 genecard
  '9' = 'Fibroblasts', 
  '10' = 'Luminal cells', 
  '11' = 'Luminal cells', 
  '12' = 'Luminal cells', 
  '13' = 'Fibroblasts', 
  '14' = 'Luminal cells', 
  '15' = 'Luminal cells', 
  '16' = 'Luminal cells', 
  '17' = 'Luminal cells', 
  '18' = 'Luminal cells', 
  '19' = 'Fibroblasts' 
)

ADPC6 <- RenameIdents(ADPC6, cell_cluster)
ADPC6$annotation <- Idents(ADPC6)

pdf(file.path(dr_outdir, str_glue('Dimplot dim{dims} annotation {reduction.use}.pdf')), width = 8, height = 6)
DimPlot(ADPC6, group.by = 'annotation', reduction = 'umap', label = T) + NoLegend()
dev.off()

cellType_markers <-
  FindAllMarkers(
    JoinLayers(ADPC6),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )

cellType_top5_markers <- 
  cellType_markers %>% 
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

names(cellType_top5_markers) <- sub('.*\\.','',names(cellType_top5_markers))


major_marker <- c(
  'Luminal cells' = 'KRT8',
  # 'Luminal cells' = 'KRT18',
  'Luminal cells' = 'NKX3-1',
  'Luminal cells' = 'KLK4',
  'Luminal cells' = 'KLK3',
  'Luminal cells' = 'KLK2',
  'Luminal cells' = 'MSMB',
  'Luminal cells' = 'ACPP',

  'Fibroblasts' = 'DCN',
  'Fibroblasts' = 'LUM',
  'Fibroblasts' = 'VIM',
  'Fibroblasts' = 'IGFBP7',
  'Fibroblasts' = 'COL1A1',
  'Fibroblasts' = 'COL1A2'
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} annotation main_marker.pdf")), width = 25, height = 10)
DotPlot(ADPC6, group.by = "annotation", features = major_marker)+ RotatedAxis() + theme_90angle
DotPlot(ADPC6, group.by = "annotation", features = cellType_top5_markers)+ RotatedAxis() + theme_90angle
dev.off()

# saveRDS(ADPC6, 'result_8um/ADPC6_mainCellType.rds')




## Luminal subype ---------------------------------------------------------------
luminal <- subset(ADPC6, annotation == 'Luminal cells')
luminal_out_dir <- "result_8um/Luminal_subtype" # dir.create(luminal_out_dir)

### PCA --------
luminal <- FindVariableFeatures(luminal)
luminal <- ScaleData(luminal, vars.to.regress = ncount_name)
luminal <- RunPCA(luminal)

#### pca inspector ------
pdf(file.path(luminal_out_dir, 'pca plot_ADPC.pdf'), height = 6, width = 8)
ElbowPlot(luminal, ndims = 50)
dev.off()

### dimension reduction --------
dims <- 10

library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

luminal <- RunUMAP(luminal, dims = 1:dims,return.model = TRUE)

reduction.use <- 'pca' # 'harmony' 'pca'


#### cluster -------
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
pdf(file.path(luminal_out_dir, str_glue("Dimplot dim{dims} cluster_ADPC.pdf")), height = 6, width = 8)
DimPlot(luminal, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(luminal, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(luminal, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(luminal, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = 'umap') + NoLegend()
dev.off()

Idents(luminal) <- "Spatial.008um_snn_res.0.6"

### refiltering -----------------------------------------------------------------
Idents(luminal) <- "Spatial.008um_snn_res.0.2"
luminal2 <- subset(luminal, idents = c(0:3, 6, 15, 16))
pdf(file.path(luminal_out_dir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_ADPC_split2.pdf")), height = 6, width = 30)
DimPlot(luminal2, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.2') + NoLegend()
DimPlot(luminal2, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.4') + NoLegend()
DimPlot(luminal2, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.6') + NoLegend()
DimPlot(luminal2, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.8') + NoLegend()
dev.off()

### PCA --------
luminal2 <- FindVariableFeatures(luminal2)
luminal2 <- ScaleData(luminal2, vars.to.regress = ncount_name)
luminal2 <- RunPCA(luminal2)

#### pca inspector ------
pdf(file.path(luminal2_out_dir, 'pca plot_ADPC.pdf'), height = 6, width = 8)
ElbowPlot(luminal2, ndims = 50)
dev.off()

### dimension reduction --------
dims <- 10

library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

luminal2 <- RunUMAP(luminal2, dims = 1:dims,return.model = TRUE)

reduction.use <- 'pca' # 'harmony' 'pca'


#### cluster -------
luminal2  <-
  FindNeighbors(
    luminal2 ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


luminal2 <-
  FindClusters(
    luminal2,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(luminal_out_dir, str_glue("Dimplot dim{dims} cluster_ADPC_2.pdf")), height = 6, width = 8)
DimPlot(luminal2, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(luminal2, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(luminal2, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(luminal2, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = 'umap') + NoLegend()
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

pdf(file.path(luminal_out_dir, str_glue("dotplot dims{dims} major_marker_ADPC_2_multi_res.pdf")), width = 45, height = 10)
DotPlot(luminal2, group.by = "Spatial.008um_snn_res.0.2", features = luminal_marker_ls)+RotatedAxis()
DotPlot(luminal2, group.by = "Spatial.008um_snn_res.0.4", features = luminal_marker_ls)+RotatedAxis()
DotPlot(luminal2, group.by = "Spatial.008um_snn_res.0.6", features = luminal_marker_ls)+RotatedAxis()
DotPlot(luminal2, group.by = "Spatial.008um_snn_res.0.8", features = luminal_marker_ls)+RotatedAxis()
dev.off()




### refiltering -----------------------------------------------------------------
Idents(luminal2) <- "Spatial.008um_snn_res.0.6"
luminal3 <- subset(luminal2, idents = c(0:10, 12:15))
pdf(file.path(luminal_out_dir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_ADPC_split3.pdf")), height = 6, width = 30)
DimPlot(luminal3, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.2') + NoLegend()
DimPlot(luminal3, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.4') + NoLegend()
DimPlot(luminal3, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.6') + NoLegend()
DimPlot(luminal3, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.8') + NoLegend()
dev.off()

### PCA --------
luminal3 <- FindVariableFeatures(luminal3)
luminal3 <- ScaleData(luminal3, vars.to.regress = ncount_name)
luminal3 <- RunPCA(luminal3)

#### pca inspector ------
pdf(file.path(luminal_out_dir, 'pca plot_ADPC.pdf'), height = 6, width = 8)
ElbowPlot(luminal3, ndims = 50)
dev.off()

pdf(file.path(luminal_out_dir, 'dimheatmap 1-5_ADPC.pdf'), height = 6, width = 10)
DimHeatmap(luminal3, dims = 1:5)
dev.off()
pdf(file.path(luminal_out_dir, 'dimheatmap 6-10_ADPC.pdf'), height = 6, width = 10)
DimHeatmap(luminal3, dims = 6:10)
dev.off()
pdf(file.path(luminal_out_dir, 'dimheatmap 11-15_ADPC.pdf'), height = 6, width = 10)
DimHeatmap(luminal3, dims = 11:15)
dev.off()
pdf(file.path(luminal_out_dir, 'dimheatmap 16-20_ADPC.pdf'), height = 6, width = 10)
DimHeatmap(luminal3, dims = 16:20)
dev.off()

### dimension reduction --------
dims <- 10

library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

luminal3 <- RunUMAP(luminal3, dims = 1:dims,return.model = TRUE)

reduction.use <- 'pca' # 'harmony' 'pca'


#### cluster -------
luminal3  <-
  FindNeighbors(
    luminal3 ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


luminal3 <-
  FindClusters(
    luminal3,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(luminal_out_dir, str_glue("Dimplot dim{dims} cluster_ADPC_3.pdf")), height = 6, width = 8)
DimPlot(luminal3, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(luminal3, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(luminal3, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(luminal3, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = 'umap') + NoLegend()
dev.off()

Idents(luminal3) <- "Spatial.008um_snn_res.0.6"
luminal3_markers <-
  FindAllMarkers(
    JoinLayers(luminal3),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )

luminal3_top5_markers <- 
  luminal3_markers %>% 
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

names(luminal3_top5_markers) <- sub('.*\\.','',names(luminal3_top5_markers))

# dotplot
pdf(file.path(luminal_out_dir, str_glue("dotplot dims{dims} major_marker_ADPC_3.pdf")), width = 45, height = 10)
DotPlot(luminal3, group.by = "Spatial.008um_snn_res.0.6", features = major_marker_ls)+RotatedAxis()
DotPlot(luminal3, group.by = "Spatial.008um_snn_res.0.6", features = luminal_marker_ls)+RotatedAxis()
DotPlot(luminal3, group.by = "Spatial.008um_snn_res.0.6", features = luminal3_top5_markers)+RotatedAxis()
dev.off()





### annotation ----------------------------------------------------------------
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
  '13'  = 'Adeno-malignant luminal'
)

# Idents(luminal3) <- "Spatial.008um_snn_res.0.6"

luminal3 <- RenameIdents(luminal3, cell_cluster)
luminal3$subtype <- Idents(luminal3)

pdf(file.path(luminal_out_dir, str_glue("Dimplot dim{dims} annotation ADPC.pdf")), height = 6, width = 8)
DimPlot(luminal3, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

# saveRDS(luminal3, file.path(luminal_out_dir, 'luminal3_ADPC.rds'))

ADPC7 <- subset(ADPC6, cells = setdiff(Cells(ADPC), setdiff(Cells(luminal), Cells(luminal3))))




## 映射-------------------------------
library(tibble)
library(data.table)
subtype <- luminal3[['subtype']] %>% rownames_to_column('cells')

meta.data.t <-
  subtype%>%
  right_join(
    ADPC7[['annotation']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation2 = ifelse((as.character(annotation) %in% c('Luminal cells')),
                         as.character(subtype), 
                         as.character(annotation)))%>%
  dplyr::select(cells, annotation2) %>%
  column_to_rownames('cells')


ADPC7 <- AddMetaData(ADPC7,meta.data.t['annotation2'])




## cell type dotplot ------------------------------------------------------------
luminal3_cellTye_markers <- c(
  'Adeno-malignant luminal3' = 'AMACR',
  'Adeno-malignant luminal3' = 'COL9A2',
  'Adeno-malignant luminal3' = 'KLK4',
  'Adeno-malignant luminal3' = 'KLK3',
  'Adeno-malignant luminal3' = 'KLK2',

  'Fibroblasts' = 'DCN',
  'Fibroblasts' = 'LUM',
  'Fibroblasts' = 'VIM',
  'Fibroblasts' = 'IGFBP7',
  'Fibroblasts' = 'COL1A1',
  'Fibroblasts' = 'COL1A2'
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

pdf(file.path(luminal_out_dir, "dotplot celltype main_marker_ADPC.pdf"), width = 25, height = 8)
DotPlot(ADPC7, group.by = "annotation2", features = luminal3_cellTye_markers, dot.scale = 10) + theme_90angle
dev.off()


# saveRDS(ADPC7, 'result_8um/ADPC7_subCellType.rds')




## cellType distribution ---------------------------------------------------------------
ADPC7_meta <-
  ADPC7[[]] %>%
  rownames_to_column("cell_id")


p_col_proj <-
  ADPC7_meta %>%
  group_by(project, annotation2) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = project, y = freq, fill = annotation2)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(dr_outdir, "cell_freq_project_ADPC.pdf"),
  height = 6,
  width  = 6
)




### spatial plot ------
ggsave(
  SpatialDimPlot(ADPC7, group.by = 'annotation2', pt.size.factor = 6, image.scale = 'hires'),
  filename = file.path(sp_output, "Saptial_cell_type_ADPC.pdf"),
  height = 10,
  width  = 20
)




# NEPC --------------------------------------------------------------------------
## qc filter --------------------------------------------------------------------
qc_metric <- read.csv("result_8um/NEPC_qc_metric.csv") %>% map_dfc(str_trim)

qc_metric_variable <- c(
  nfeature_name, ncount_name,
  "percent.mt", "percent.hb", "log10GenesPerUMI"
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
  names(visim_ls[2]) %>%
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
  names(visim_ls[2]) %>%
  setNames(., .) %>%
  map(function(sample){
    all_data <- visim_ls[[sample]] %>% subset(downsample = 100000, seed = 1)
    all_data_pass <- visim_ls_pass[[sample]] %>% subset(downsample = 100000, seed  = 1)
    
    # all_data
    plot1 <- FeatureScatter(all_data, feature1 = ncount_name, feature2 = nfeature_name) +
      patchwork::plot_annotation(title = sample)
    plot2 <- FeatureScatter(all_data, feature1 = ncount_name, feature2 = "percent.mt") +
      patchwork::plot_annotation(title = sample)
    plot3 <- FeatureScatter(all_data, feature1 = ncount_name, feature2 = "percent.hb") +
      patchwork::plot_annotation(title = sample)
    plot4 <- FeatureScatter(all_data, feature1 = ncount_name, feature2 = "log10GenesPerUMI") +
      patchwork::plot_annotation(title = sample)
    p1 <- plot1 + plot2 + plot3 + plot4 + patchwork::plot_layout(ncol = 4)
    
    # all_data_pass
    plot1_pass <- FeatureScatter(all_data_pass, feature1 = ncount_name, feature2 = nfeature_name)+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    plot2_pass <- FeatureScatter(all_data_pass, feature1 = ncount_name, feature2 = "percent.mt")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    plot3_pass <- FeatureScatter(all_data_pass, feature1 = ncount_name, feature2 = "percent.hb")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    plot4_pass <- FeatureScatter(all_data_pass, feature1 = ncount_name, feature2 = "log10GenesPerUMI")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    p2 <- plot1_pass + plot2_pass + plot3_pass + plot4_pass + patchwork::plot_layout(ncol = 4)
    
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
  file.path(qc_outdir, "qc plot - pass_NEPC.pdf"),
  height = 12, 
  width  = 18
)
qc_plot_filter_ls %>% walk(print)
dev.off()

qc_metric_ls %>%
  enframe("sample", "qc") %>%
  unnest(qc) %>% 
  openxlsx::write.xlsx(
    file = file.path(qc_outdir, "qc_metric_NEPC.xlsx"),
    asTable = TRUE,
    withFilter = FALSE
  )

# cells filter
cell_counts_df <-
  names(visim_ls[2]) %>%
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
    file = file.path(qc_outdir, "cell count after filter_NEPC.xlsx"),
    asTable = TRUE,
    withFilter = FALSE
  )
# saveRDS(visim_ls_pass$NEPC, 'result_8um/NEPC_pass.rds')

SpatialDimPlot(visim_ls_pass$NEPC, group.by = 'project')




## dimension reduction and clustering -------------------------------------------
NEPC <- visim_ls_pass$NEPC
#SpatialFeaturePlot(NEPC, features = c('HBA2', 'HBB', ncount_name))

### normalize -------------------------------------------------------------------
NEPC <- NormalizeData(NEPC, scale.factor = median(NEPC[[ncount_name]] %>% unlist))
NEPC <- FindVariableFeatures(NEPC, nfeatures = 3000)
NEPC <- ScaleData(NEPC,vars.to.regress = ncount_name)

NEPC <- RunPCA(NEPC, features = VariableFeatures(NEPC))

pdf(file = file.path(dr_outdir, 'pca plot_NEPC.pdf'), height = 6, width = 8)
ElbowPlot(NEPC, ndims = 30)
dev.off()


dims <- 10


#undebug(RunUMAP)


NEPC <- RunUMAP(NEPC, dims = 1:dims, return.model = TRUE)

DimPlot(NEPC, label = TRUE, reduction = "umap") + NoLegend()

reduction.use <- 'pca' # 'harmony' 'pca'




### cluster ---------------------------------------------------------------------
NEPC <-
  FindNeighbors(
    NEPC,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )

NEPC <-
  FindClusters(
    NEPC,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )


# dimplot
pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_NEPC.pdf")), height = 6, width = 8)
DimPlot(NEPC, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(NEPC, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(NEPC, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(NEPC, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
dev.off()

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_NEPC_split.pdf")), height = 6, width = 30)
DimPlot(NEPC, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.2') + NoLegend()
DimPlot(NEPC, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.4') + NoLegend()
DimPlot(NEPC, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.6') + NoLegend()
DimPlot(NEPC, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.8') + NoLegend()
dev.off()




## refiltering -----------------------------------------------------------------
Idents(NEPC) <- "Spatial.008um_snn_res.0.2"
NEPC2 <- subset(NEPC, idents = c(0:6))
pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_NEPC_split2.pdf")), height = 6, width = 30)
DimPlot(NEPC2, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.2') + NoLegend()
DimPlot(NEPC2, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.4') + NoLegend()
DimPlot(NEPC2, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.6') + NoLegend()
DimPlot(NEPC2, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.8') + NoLegend()
dev.off()

Idents(NEPC2) <- "Spatial.008um_snn_res.0.8"
NEPC3 <- subset(NEPC2, idents = c(0:14))
pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_NEPC_split3.pdf")), height = 6, width = 30)
DimPlot(NEPC3, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.2') + NoLegend()
DimPlot(NEPC3, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.4') + NoLegend()
DimPlot(NEPC3, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.6') + NoLegend()
DimPlot(NEPC3, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.8') + NoLegend()
dev.off()

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_NEPC3.pdf")), height = 6, width = 8)
DimPlot(NEPC3, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(NEPC3, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(NEPC3, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(NEPC3, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
dev.off()




## dimension reduction and clustering -------------------------------------------
NEPC3 <- RunPCA(NEPC3, features = VariableFeatures(NEPC3))

pdf(file = file.path(dr_outdir, 'pca plot_NEPC3.pdf'), height = 6, width = 8)
ElbowPlot(NEPC3, ndims = 30)
dev.off()

dims <- 10


NEPC3 <- RunUMAP(NEPC3, dims = 1:dims, return.model = TRUE)

DimPlot(NEPC3, label = TRUE, reduction = "umap") + NoLegend()

reduction.use <- 'pca' # 'harmony' 'pca'




#### cluster ---------------------------------------------------------------------
NEPC3 <-
  FindNeighbors(
    NEPC3,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )

NEPC3 <-
  FindClusters(
    NEPC3,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )


# dimplot
pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_NEPC3_2.pdf")), height = 6, width = 8)
DimPlot(NEPC3, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(NEPC3, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(NEPC3, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(NEPC3, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
dev.off()

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_NEPC_split3_2.pdf")), height = 6, width = 30)
DimPlot(NEPC3, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.2') + NoLegend()
DimPlot(NEPC3, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.4') + NoLegend()
DimPlot(NEPC3, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.6') + NoLegend()
DimPlot(NEPC3, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.8') + NoLegend()
dev.off()




## refiltering -----------------------------------------------------------------
Idents(NEPC3) <- "Spatial.008um_snn_res.0.6"
NEPC4 <- subset(NEPC3, idents = c(0:11, 13))
pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_NEPC_split4.pdf")), height = 6, width = 30)
DimPlot(NEPC4, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.2') + NoLegend()
DimPlot(NEPC4, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.4') + NoLegend()
DimPlot(NEPC4, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.6') + NoLegend()
DimPlot(NEPC4, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.8') + NoLegend()
dev.off()

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_NEPC4.pdf")), height = 6, width = 8)
DimPlot(NEPC4, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(NEPC4, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(NEPC4, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(NEPC4, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
dev.off()




## dimension reduction and clustering -------------------------------------------
NEPC4 <- RunPCA(NEPC4, features = VariableFeatures(NEPC4))

pdf(file = file.path(dr_outdir, 'pca plot_NEPC4.pdf'), height = 6, width = 8)
ElbowPlot(NEPC4, ndims = 30)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 1-5_NEPC.pdf'), height = 6, width = 10)
DimHeatmap(NEPC4, dims = 1:5)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 6-10_NEPC.pdf'), height = 6, width = 10)
DimHeatmap(NEPC4, dims = 6:10)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 11-15_NEPC.pdf'), height = 6, width = 10)
DimHeatmap(NEPC4, dims = 11:15)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 16-20_NEPC.pdf'), height = 6, width = 10)
DimHeatmap(NEPC4, dims = 16:20)
dev.off()

dims <- 10


NEPC4 <- RunUMAP(NEPC4, dims = 1:dims, return.model = TRUE)

DimPlot(NEPC4, label = TRUE, reduction = "umap") + NoLegend()

reduction.use <- 'pca' # 'harmony' 'pca'




#### cluster ---------------------------------------------------------------------
NEPC4 <-
  FindNeighbors(
    NEPC4,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )

NEPC4 <-
  FindClusters(
    NEPC4,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )


# dimplot
pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_NEPC4_2.pdf")), height = 6, width = 8)
DimPlot(NEPC4, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(NEPC4, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(NEPC4, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(NEPC4, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
dev.off()

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_NEPC_split4_2.pdf")), height = 6, width = 30)
DimPlot(NEPC4, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.2') + NoLegend()
DimPlot(NEPC4, group.by = "Spatial.008um_snn_res.0.4", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.4') + NoLegend()
DimPlot(NEPC4, group.by = "Spatial.008um_snn_res.0.6", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.6') + NoLegend()
DimPlot(NEPC4, group.by = "Spatial.008um_snn_res.0.8", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.8') + NoLegend()
dev.off()

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

major_marker_ls <- major_marker_ls[major_marker_ls %in% rownames(LayerData(JoinLayers(NEPC4), layer = 'data'))]

pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} multi res_NEPC.pdf")), width = 55, height = 8)
DotPlot(NEPC4, group.by = "Spatial.008um_snn_res.0.2", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(NEPC4, group.by = "Spatial.008um_snn_res.0.4", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(NEPC4, group.by = "Spatial.008um_snn_res.0.6", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(NEPC4, group.by = "Spatial.008um_snn_res.0.8", features = major_marker_ls, dot.scale = 10) + theme_90angle
dev.off()




### dotplot ---------------------------------------------------------------------
pdf(
  file.path(dr_outdir, str_glue("featureplot dims{dims} {reduction.use} major_marker_NEPC.pdf")), 
  width = 25, 
  height = 100
)
FeaturePlot(
  NEPC4,
  features = major_marker_ls,
  ncol = 5,
  reduction = 'umap',
  raster = TRUE
)  
dev.off()


Idents(NEPC4) <- "Spatial.008um_snn_res.0.6"
NEPC4_markers <-
  FindAllMarkers(
    JoinLayers(NEPC4),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )

NEPC4_top5_markers <- 
  NEPC4_markers %>% 
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

names(NEPC4_top5_markers) <- sub('.*\\.','',names(NEPC4_top5_markers))

# dotplot
pdf(file.path(dr_outdir, str_glue("dotplot dims {dims} {reduction.use} majorMarker_topDEGs_NEPC.pdf")), width = 40, height = 6)
DotPlot(NEPC4, group.by = "Spatial.008um_snn_res.0.6", features = major_marker_ls)+RotatedAxis()
DotPlot(NEPC4, group.by = "Spatial.008um_snn_res.0.6", features = NEPC4_top5_markers)+RotatedAxis()
dev.off()




## Annotaiton -------------------------------------------------------------------------------
Idents(NEPC4) <- "Spatial.008um_snn_res.0.6"
cell_cluster <- c(
  '0' = 'NE cells', 
  '1' = 'NE cells', 
  '2' = 'NE cells',  
  '3' = 'NE cells', 
  '4' = 'Club cells', 
  '5' = 'Plasma cells', 
  '6' = 'Smooth muscle cells',
  '7' = 'NE cells', 
  '8' = 'Fibroblasts', #Crisp3 genecard
  '9' = 'Macrophages', 
  '10' = 'Pericytes', 
  '11' = 'NE cells', 
  '12' = 'Plasma cells'
)

NEPC4 <- RenameIdents(NEPC4, cell_cluster)
NEPC4$annotation <- Idents(NEPC4)

pdf(file.path(dr_outdir, str_glue('Dimplot dim{dims} annotation {reduction.use}_NEPC.pdf')), width = 8, height = 6)
DimPlot(NEPC4, group.by = 'annotation', reduction = 'umap', label = T) + NoLegend()
dev.off()

cellType_markers <-
  FindAllMarkers(
    JoinLayers(NEPC4),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )

cellType_top5_markers <- 
  cellType_markers %>% 
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

names(cellType_top5_markers) <- sub('.*\\.','',names(cellType_top5_markers))


major_marker <- c(
  'NE cells' = 'ASCL-1',
  'NE cells' = 'FoxA2',
  'NE cells' = 'CHGA',
  'NE cells' = 'CHGB',
  'NE cells' = 'SYP',
  'NE cells' = 'NKX2-1',
  'NE cells' = 'ENO2',
  'NE cells' = 'MYCN',
  'NE cells' = 'POUF3F2',
  'NE cells' = 'LMO3',
  'NE cells' = 'INSM1',
  'NE cells' = 'SOX2',
  'NE cells' = 'EZH2',
  'NE cells' = 'SIAH2',
  
  'Club cells' = 'SCGB1A3',
  #  'Club cells' = 'WFDC2',
  'Club cells' = 'LCN2',
  'Club cells' = 'MMP7',	
  'Club cells' = 'KRT4',	
  'Club cells' = 'TACSTD2',
  'Club cells' = 'SCGB3A1',
  'Club cells' = 'LTF',
  'Club cells' = 'RARRES1',
  'Club cells' = 'KRT7',
  'Club cells' = 'AGR2',
  
  "Plasma cells" = "IGHG1",
  "Plasma cells" = "IGHA1",
  "Plasma cells" = "IGHG3",
  "Plasma cells" = "IGKC",
  "Plasma cells" = "TNFRSF17",
  "Plasma cells" = "JCHAIN",
  "Plasma cells" = "MZB1",
  "Plasma cells" = "XBP1",
  'Plasma cells' = 'SEC11C',	
  'Plasma cells' = 'PRDX4',
  'Plasma cells' = 'SPCS2',
  'Plasma cells' = 'SSR3',
  'Plasma cells' = 'SDF2L1',
  'Plasma cells' = 'C19orf10',
  'Plasma cells' = 'MANF',
  'Plasma cells' = 'TMEM258',
  'Plasma cells' = 'DNAJB9',
  'Plasma cells' = 'CD138',
  'Plasma cells' = 'CD38',
  'Plasma cells' = 'BCMA',
  'Plasma cells' = 'IRF4',
  'Plasma cells' = 'BLIMP1',
  
  #  'Smooth muscle cells' = 'ACTA2',
  'Smooth muscle cells' = 'MYH11',
  'Smooth muscle cells' = 'TAGLN',
  'Smooth muscle cells' = 'CNN1',
  'Smooth muscle cells' = 'DES',
  'Smooth muscle cells' = 'ACTG2',
  
  'Fibroblasts' = 'DCN',
  'Fibroblasts' = 'LUM',
  'Fibroblasts' = 'VIM',
  'Fibroblasts' = 'IGFBP7',
  'Fibroblasts' = 'COL1A1',
  'Fibroblasts' = 'COL1A2',
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "SFRP4",
  "Fibroblasts" = "IGF1",
  "Fibroblasts" = "APOD",
  "Fibroblasts" = "FBLN1",
  "Fibroblasts" = "MEG3",
  "Fibroblasts" = "CXCL12",
  
  "Macrophages" = "CD68",
  "Macrophages" = "CD163",
  "Macrophages" = "C1QA",
  "Macrophages" = "C1QC", 
  'Macrophages'  = 'CSF1R',
  'Macrophages' = 'C1QB',
  'Macrophages' = 'GPR34',
  'Macrophages' = 'MS4A4A',
  
  "Pericytes" = "PDGFRB",
  "Pericytes" = "RGS5",
  "Pericytes" = "MCAM",
  'Pericytes' = 'ACTA2',
  #  'Pericytes' = 'MYH11',
  'Pericytes' = 'MT1M',
  'Pericytes' = 'FRZB',
  'Pericytes' = 'MT1A',
  'Pericytes' = 'NDUFA4L2',
  'Pericytes' = 'PPP1R14A',
  'Pericytes' = 'MYLK',
  'Pericytes' = 'PHLDA1'
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} annotation main_marker_NEPC.pdf")), width = 25, height = 10)
DotPlot(NEPC4, group.by = "annotation", features = major_marker)+ RotatedAxis() + theme_90angle
DotPlot(NEPC4, group.by = "annotation", features = cellType_top5_markers)+ RotatedAxis() + theme_90angle
dev.off()

# saveRDS(NEPC4, 'result_8um/NEPC4_mainCellType.rds')




13## cellType distribution ---------------------------------------------------------------
NEPC4_meta <-
  NEPC4[[]] %>%
  rownames_to_column("cell_id")


p_col_proj <-
  NEPC4_meta %>%
  group_by(project, annotation) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = project, y = freq, fill = annotation)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(dr_outdir, "cell_freq_project_NEPC.pdf"),
  height = 6,
  width  = 6
)




### spatial plot ------
ggsave(
  SpatialDimPlot(NEPC4, group.by = 'annotation', pt.size.factor = 6, image.scale = 'hires'),
  filename = file.path(sp_output, "Saptial_cell_type_NEPC.pdf"),
  height = 10,
  width  = 20
)




# Annotation info output -------------------------------------------------------
ADPC7$annotation2 %>% table()
cellType_anno_ADPC <- 
  ADPC7@meta.data[c('project', 'annotation2')] %>% 
  rownames_to_column('Barcode') %>% 
  mutate(
    cell_id = str_replace(Barcode, '^.*?_', ''), 
    group = annotation2
  ) %>% 
  dplyr::select(cell_id, group)

cellType_anno_ADPC %>% 
  write.csv(str_glue('result_8um/cell type ADPC.csv'), row.names = F)


NEPC4$annotation %>% table()
cellType_anno_NEPC <- 
  NEPC4@meta.data[c('project', 'annotation')] %>% 
  rownames_to_column('Barcode') %>% 
  mutate(
    cell_id = str_replace(Barcode, '^.*?_', ''), 
    group = annotation
  ) %>% 
  dplyr::select(cell_id, group)

cellType_anno_NEPC %>% 
  write.csv(str_glue('result_8um/cell type NEPC.csv'), row.names = F)




# RCTD --------------------------------------------------------------------------
library(spacexr)
rctd_out_dir <- 'result_8um/RCTD/'
sc_NE <- readRDS(file.path(rctd_out_dir, 'sc_NE.rds'))
sc_primary <- readRDS(file.path(rctd_out_dir, 'sc_primary.rds'))


## NEPC -------------------------------------------------------------------------
sc_NE@meta.data %>% head()
sc_NE$tissue_type %>% table()
sc_NE$cell_type_1 %>% table()
sc_NE$predicted_labels %>% table()

sc_NE_anno <- as.factor(sc_NE$cell_type_1)
names(sc_NE_anno) <- rownames(sc_NE@meta.data)

### make sc ref -----------------------------------------------------------------
sc_NE_ref <- Reference(
  GetAssayData(sc_NE, layer = 'counts'), 
  sc_NE_anno, 
  sc_NE$nCount_RNA
)

### prepare spatial data --------------------------------------------------------
spatial_rna_NE <- 
  SpatialRNA(
    GetTissueCoordinates(NEPC4)[1:2], 
    GetAssayData(NEPC4, layer = 'counts', assay = 'Spatial.008um'), 
    colSums(GetAssayData(NEPC4, layer = 'counts', assay = 'Spatial.008um'))
  )

### run RCTD --------------------------------------------------------------------
rctd_obj_NE <- 
  create.RCTD(
    spatial_rna_NE, 
    sc_NE_ref, 
    max_cores = 10
  )

rctd_result_NE <- 
  run.RCTD(
    rctd_obj_NE,
    doublet_mode = 'doublet'
  )




## ADPC -------------------------------------------------------------------------
sc_primary@meta.data %>% head()
sc_primary$cell_type_1 %>% table()

sc_primary_anno <- as.factor(sc_primary$cell_type_1)
names(sc_primary_anno) <- rownames(sc_primary@meta.data)

### make sc ref -----------------------------------------------------------------
sc_primary_ref <- Reference(
  GetAssayData(sc_primary, layer = 'counts'), 
  sc_primary_anno, 
  sc_primary$total_counts
)

### prepare spatial data --------------------------------------------------------
spatial_rna_AD <- 
  SpatialRNA(
    GetTissueCoordinates(ADPC7)[1:2], 
    GetAssayData(ADPC7, layer = 'counts', assay = 'Spatial.008um'), 
    colSums(GetAssayData(ADPC7, layer = 'counts', assay = 'Spatial.008um'))
  )

### run RCTD --------------------------------------------------------------------
rctd_obj_AD <- 
  create.RCTD(
    spatial_rna_AD, 
    sc_primary_ref, 
    max_cores = 10
  )

rctd_result_AD <- 
  run.RCTD(
    rctd_obj_AD,
    doublet_mode = 'doublet'
  )




## visualization ---------------------------------------------------------------
allcolour=c('B cell' = "#DC143C", 
            'DC' = "#0000FF", 
            'Endothelial' = "#20B2AA",
            'Fibroblast' = "#FFA500",
            'Macrophage' = "#9370DB",
            'Monocyte' = "#98FB98",
            'NE cell' = "#8B4513",
            'T cell' = "#F08080",
            'Unknown' = "#808000",
            'Basal cel' = "#FF00FF",
            'Club cel' = "#FA8072",
            'Luminal cel' = "#87CEEB",
            'Mast' = "#40E0D0")
NEPC4 <- AddMetaData(NEPC4, metadata = rctd_result_NE@results$results_df)
NEPC4$first_type <- as.character(NEPC4$first_type)
NEPC4$first_type[is.na(NEPC4$first_type)] <- "Unknown"

ADPC7 <- AddMetaData(ADPC7, metadata = rctd_result_AD@results$results_df)
ADPC7$first_type <- as.character(ADPC7$first_type)
ADPC7$first_type[is.na(ADPC7$first_type)] <- "Unknown"
ADPC7$first_type %>% table()

visium_rctd_ls <- list(NEPC4, ADPC7)
names(visium_rctd_ls) <- c('NEPC', 'ADPC')


image_scale_ls <-
  visium_rctd_ls %>%
  # visium_roi_ls %>%
  map(GetTissueCoordinates) %>%
  map('x') %>%
  map(range) %>%
  map(diff) %>%
  map(log) %>%
  {map(., function(x, max_val)x/max_val, max_val = max(unlist(.)))}

pdf(
  file.path(rctd_out_dir, 'spatial plot of RCTD.pdf'), 
  height = 15, 
  width = 15
)
visium_rctd_ls %>%
  imap(function(visium, name){
    if(name == "HC") name_fix <- "NML" else name_fix <- name
    
    # visium$first_type[visium$first_type == "Unknown"] <- NA
    # visium$first_type[visium$first_type == "Other"] <- NA
    
    SpatialDimPlot(
      visium,
      group.by       = 'first_type',
      # cols = allcolour, 
      image.scale    = 'hires',
      pt.size.factor = 5 * image_scale_ls[[name]],
      shape          = 22,
      image.alpha = 0.5
      # alpha          = 0.5
    ) +
      ggtitle(name_fix) +
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5)))
  })
dev.off()


## save --------------------------------------------------------------------------
visium_rctd <- merge(visium_rctd_ls[[1]], visium_rctd_ls[-1])
saveRDS(visium_rctd, file.path(rctd_out_dir, "visium_rctd.rds"))

# export csv
visium_rctd %>%
  GetTissueCoordinates() %>%
  write.csv(file = file.path(rctd_out_dir, "TissueCoordinates.csv"), row.names = FALSE)


visium_rctd[[]] %>%
  write.csv(file = file.path(rctd_out_dir, "rctd predicted table.csv"), row.names = FALSE)




# Annotation info output -------------------------------------------------------
ADPC7$first_type %>% table()
cellType_anno_ADPC <- 
  ADPC7@meta.data[c('project', 'first_type')] %>% 
  rownames_to_column('Barcode') %>% 
  mutate(
    cell_id = str_replace(Barcode, '^.*?_', ''), 
    group = first_type
  ) %>% 
  dplyr::select(cell_id, group)

cellType_anno_ADPC %>% 
  write.csv(str_glue('result_8um/RCTD/cell type ADPC_RCTD.csv'), row.names = F)


NEPC4$first_type %>% table()
cellType_anno_NEPC <- 
  NEPC4@meta.data[c('project', 'first_type')] %>% 
  rownames_to_column('Barcode') %>% 
  mutate(
    cell_id = str_replace(Barcode, '^.*?_', ''), 
    group = first_type
  ) %>% 
  dplyr::select(cell_id, group)

cellType_anno_NEPC %>% 
  write.csv(str_glue('result_8um/RCTD/cell type NEPC_RCTD.csv'), row.names = F)





## cellType distribution ---------------------------------------------------------------
visium_rctd_meta <-
  visium_rctd[[]] %>%
  rownames_to_column("cell_id")


p_col_proj <-
  visium_rctd_meta %>%
  group_by(project, first_type) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = project, y = freq, fill = first_type)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(rctd_out_dir, "cell_freq_project_RCTD.pdf"),
  height = 6,
  width  = 6
)


p_bar_group  <-
  visium_rctd_meta %>%
  group_by( first_type,project) %>%
  summarise(count = n()) %>% 
  ungroup() %>% 
  tidyr::complete(first_type, project, fill = list(count = 0)) %>%
  mutate(freq = count/sum(count)) %>% 
  
  
  ggplot(aes(x = first_type, y = freq, fill = project)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
  ggsci::scale_fill_igv() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_bar_group,
  filename = file.path(rctd_out_dir, "box_cell_freq_RCTD.pdf"),
  height = 8,
  width  = 15
)




# reAnnotation ------------------------------------------------------------------
visium_rctd <- readRDS("result_8um/RCTD/visium_rctd.rds")
ADPC_firstType <- 
  visium_rctd@meta.data %>% 
  rownames_to_column('cell_id') %>% 
  dplyr::filter(project == 'ADPC') %>% 
  dplyr::select(cell_id, first_type)

ADPC7 <- readRDS('result_8um/ADPC7_subCellType.rds')
ADPC7_meta <- 
  ADPC7@meta.data %>% 
  rownames_to_column('cell_id') %>% 
  left_join(ADPC_firstType, by = 'cell_id') %>%
  column_to_rownames('cell_id')

ADPC7 <- AddMetaData(ADPC7, metadata = ADPC7_meta['first_type'])
ADPC7$first_type %>% table()

ADPC7$annotation_full <- 
  ifelse(
    (as.character(ADPC7$first_type) %in% c('Unknown')), 
    as.character(ADPC7$annotation), 
    as.character(ADPC7$first_type)
) %>%
  str_replace_all(c(
    'Luminal cells' = 'Luminal cell', 
    'Fibroblasts' = 'Fibroblast'
  ))


pdf(file.path(dr_outdir, str_glue('Dimplot dim{dims} annotation {reduction.use}_ADPC7.pdf')), width = 9, height = 6)
DimPlot(ADPC7, group.by = 'annotation_full', reduction = 'umap', label = T)
dev.off()


ggsave(
  SpatialDimPlot(ADPC7, group.by = 'annotation_full', pt.size.factor = 6, image.scale = 'hires'),
  filename = file.path(sp_output, "Saptial_cell_type_ADPC7.pdf"),
  height = 10,
  width  = 20
)

ADPC7$annotation_full %>% table()
cellType_anno_ADPC <- 
  ADPC7@meta.data[c('project', 'annotation_full')] %>% 
  rownames_to_column('Barcode') %>% 
  mutate(
    Barcode = str_replace(Barcode, '^.*?_', ''), 
    group = annotation_full
  ) %>% 
  dplyr::select(Barcode, group)

cellType_anno_ADPC %>% 
  write.csv(str_glue('result_8um/cell type ADPC_annotation_full.csv'), row.names = F)

saveRDS(ADPC7, 'result_8um/ADPC7_annotation.rds')


NEPC_firstType <- 
  visium_rctd@meta.data %>% 
  rownames_to_column('cell_id') %>% 
  dplyr::filter(project == 'NEPC') %>% 
  dplyr::select(cell_id, first_type)

NEPC4 <- readRDS('result_8um/NEPC4_mainCellType.rds')
NEPC4_meta <- 
  NEPC4@meta.data %>% 
  rownames_to_column('cell_id') %>% 
  left_join(NEPC_firstType, by = 'cell_id') %>%
  column_to_rownames('cell_id')

NEPC4 <- AddMetaData(NEPC4, metadata = NEPC4_meta['first_type'])

NEPC4$first_type %>% table()
NEPC4$annotation %>% table()

NEPC4$annotation2 <- 
  ifelse(
    (as.character(NEPC4$annotation %in% c('Club cells', 'Plasma cells', 'Smooth muscle cells'))), 
    as.character(NEPC4$annotation), 
    as.character(NEPC4$first_type)
  )
NEPC4$annotation_full <- 
  ifelse(
    (as.character(NEPC4$annotation2) %in% c('Unknown')), 
     as.character(NEPC4$annotation), 
     as.character(NEPC4$annotation2)
  ) %>%
  str_replace_all(c(
    'Club cells' = 'Club cell', 
    'Fibroblasts' = 'Fibroblast', 
    'Macrophages' = 'Macrophage',
    'NE cells' = 'NE cell', 
    'Pericytes' = 'Pericyte', 
    'Plasma cells' = 'Plasma cell',
    'Smooth muscle cells' = 'Smooth muscle cell'
  ))


pdf(file.path(dr_outdir, str_glue('Dimplot dim{dims} annotation {reduction.use}_NEPC4.pdf')), width = 10, height = 6)
DimPlot(NEPC4, group.by = 'annotation_full', reduction = 'umap', label = T)
dev.off()


ggsave(
  SpatialDimPlot(NEPC4, group.by = 'annotation_full', pt.size.factor = 6, image.scale = 'hires'),
  filename = file.path(sp_output, "Saptial_cell_type_NEPC4.pdf"),
  height = 10,
  width  = 20
)

NEPC4$annotation_full %>% table()
cellType_anno_NEPC <- 
  NEPC4@meta.data[c('project', 'annotation_full')] %>% 
  rownames_to_column('Barcode') %>% 
  mutate(
    Barcode = str_replace(Barcode, '^.*?_', ''), 
    group = annotation_full
  ) %>% 
  dplyr::select(Barcode, group)

cellType_anno_NEPC %>% 
  write.csv(str_glue('result_8um/cell type NEPC_annotation_full.csv'), row.names = F)
saveRDS(NEPC4, 'result_8um/NEPC4_annotation.rds')




# ADPC luminal reannotation -----------------------------------------------------
ADPC7$annotation_full %>% table()
Luminal <- subset(ADPC7, annotation_full == 'Luminal cell')
Luminal_out_dir <- "result_8um/Luminal_reannotation" # dir.create(Luminal_out_dir)

### PCA --------
Luminal <- FindVariableFeatures(Luminal)
Luminal <- ScaleData(Luminal, vars.to.regress = ncount_name)
Luminal <- RunPCA(Luminal)

#### pca inspector ------
pdf(file.path(Luminal_out_dir, 'pca plot_ADPC.pdf'), height = 6, width = 8)
ElbowPlot(Luminal, ndims = 50)
dev.off()

### dimension reduction --------
dims <- 10

library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

Luminal <- RunUMAP(Luminal, dims = 1:dims,return.model = TRUE)

reduction.use <- 'pca' # 'harmony' 'pca'


#### cluster -------
Luminal  <-
  FindNeighbors(
    Luminal ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Luminal <-
  FindClusters(
    Luminal,
    resolution = c(0.1, 0.12, 0.15, 0.18, 0.2)
  )



# dimplot
pdf(file.path(Luminal_out_dir, str_glue("Dimplot dim{dims} cluster_ADPC.pdf")), height = 6, width = 8)
DimPlot(Luminal, group.by = "Spatial.008um_snn_res.0.1", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(Luminal, group.by = "Spatial.008um_snn_res.0.12", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(Luminal, group.by = "Spatial.008um_snn_res.0.15", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(Luminal, group.by = "Spatial.008um_snn_res.0.18", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(Luminal, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = 'umap') + NoLegend()
dev.off()


Luminal_markers_ls <- 
  list(0.1, 0.12, 0.15, 0.18, 0.2) %>% 
  map(function(res){ # browser()
    resolution <- print(str_glue("Spatial.008um_snn_res.{res}"))
    Idents(Luminal) <- resolution
    Luminal_markers <-
      FindAllMarkers(
        JoinLayers(Luminal),
        test.use = "wilcox",
        logfc.threshold = 0,
        min.pct = 0,
        only.pos = TRUE
      )
    
    Luminal_markers
  })
names(Luminal_markers_ls) <- c(0.1, 0.12, 0.15, 0.18, 0.2)
Luminal_markers_ls$`0.18` %>% openxlsx::write.xlsx(file.path(Luminal_out_dir, "temp.xlsx"), asTable = T)

Luminal_top5_markers_ls <- 
  Luminal_markers_ls %>% 
  map(function(Luminal_markers){ # browser()
    Luminal_top5_markers <-
      Luminal_markers %>% 
      group_by(cluster) %>%
      dplyr::filter(
        p_val_adj < 0.05,
        avg_log2FC > 1,
        pct.1-pct.2 > 0.1
      ) %>%
      # arrange(desc(avg_log2FC)) %>%
      # slice_head(n = 30) %>% 
      arrange(desc((pct.1-pct.2))) %>%
      # slice_head(n = 5) %>% 
      # arrange(desc(pct.1)) %>%
      # slice_head(n = 5) %>%
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
    
    names(Luminal_top5_markers) <- sub('.*\\.','',names(Luminal_top5_markers))
    Luminal_top5_markers  
  })


p_dotplot_ls <-
  Luminal_top5_markers_ls %>% 
  imap(function(Luminal_top5_markers, res){
    p <- DotPlot(Luminal, group.by = str_glue("Spatial.008um_snn_res.{res}"), features = Luminal_top5_markers) + RotatedAxis()
    p
  })

pdf(file.path(Luminal_out_dir, str_glue("dotplot dims{dims} multiRes_ADPC.pdf")), width = 20, height = 6)
p_dotplot_ls %>% 
  iwalk(function(p_dotplot, res){
    print(p_dotplot + ggtitle(res))
  })
dev.off()




Luminal_marker_ls <- c(
  'Adeno-malignant Luminal' = 'AMACR',
  'Adeno-malignant Luminal' = 'CACNA1D',
  'Adeno-malignant Luminal' = 'PCA3',
  'Adeno-malignant Luminal' = 'ERG',
  'Adeno-malignant Luminal' = 'FABP5',
  'Adeno-malignant Luminal' = 'COL9A2',
  'Adeno-malignant Luminal' = 'GCNT1',
  'Adeno-malignant Luminal' = 'PHGR1',
  'Adeno-malignant Luminal' = 'KLK4',
  'Adeno-malignant Luminal' = 'KLK3',
  'Adeno-malignant Luminal' = 'KLK2',
  'Adeno-malignant Luminal' = 'AR',
  
  'Normal/Normal-like Luminal' = 'ACPP',
  'Normal/Normal-like Luminal' = 'MSMB',
  'Normal/Normal-like Luminal' = 'RDH11',
  'Normal/Normal-like Luminal' = 'AZGP1'
)

# dotplot
pdf(file.path(Luminal_out_dir, str_glue("dotplot dims{dims} major_marker_ADPC.pdf")), width = 40, height = 6)
DotPlot(Luminal, group.by = "Spatial.008um_snn_res.0.18", features = major_marker_ls)+RotatedAxis()
DotPlot(Luminal, group.by = "Spatial.008um_snn_res.0.18", features = Luminal_marker_ls)+RotatedAxis()
DotPlot(Luminal, group.by = "Spatial.008um_snn_res.0.18", features = Luminal_top5_markers_ls$`0.18`)+RotatedAxis()
dev.off()




## volcano -------------------------------------------------------------------
Luminal_markers <-
  FindAllMarkers(
    JoinLayers(Luminal),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(Luminal_markers,file.path(Luminal_out_dir, 'Luminal_diff_gene_1.csv'))


library(EnhancedVolcano)
type<-unique(Luminal_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- Luminal_markers[Luminal_markers$cluster == type, ]
    p<-EnhancedVolcano(markers,
                       lab = rownames(markers),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = str_glue('{name} volcano plot'),
                       pCutoff = 10e-2,
                       FCcutoff = 0,
                       pointSize = 3.0,
                       labSize = 6.0)
    p
  })
volcano_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(Luminal_out_dir, 'volcano1', str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })


all_markers_major <- list(Luminal = Luminal_markers)
# scatter
p_volcanoplot_major_ls <- 
  all_markers_major %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 1
    dodge_width <- 0.2
    
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
    
    
    
    df_scatter %>%
      ggplot(aes(x = xx, y = y)) +
      geom_rect(
        data = df_bg_rect, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        inherit.aes = FALSE
      ) +
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_point(aes(color = group), size = 1) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
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

p_volcanoplot_major_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(Luminal_out_dir, 'volcano1', str_glue("{name} diffgene dot.pdf")),
      height = 8,
      width = 15
    )
  })




## enrichment -------------------------------------------------------------------
Luminal_de_genes_for_enrichmenet <- 
  Luminal_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.15) %>%
  dplyr::filter(abs(avg_log2FC) > 0.25) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")


enrich_ls <-
  Luminal_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichGO(
      x,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
      ont='ALL',
      pAdjustMethod = 'BH',
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.05,
      keyType = 'SYMBOL'
    )
  })


Luminal_er_output <- file.path(Luminal_out_dir, 'enrichment1') # dir.create(Luminal_er_output)
# expot enrichment table
enrich_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Luminal_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Luminal_er_output, name, "Go enriment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

# export dotplot
enrich_dot_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Luminal_er_output, name, "go dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Luminal_er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })

Luminal_de_genes_id_for_enrichmenet <-
  Luminal_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::bitr(
      x,
      OrgDb = "org.Hs.eg.db",
      fromType = "SYMBOL",
      toType = "ENTREZID"
    ) %>%
      .$ENTREZID
  }) %>%
  purrr::discard(function(x) length(x) < 1)


enrich_ls <-
  Luminal_de_genes_id_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichKEGG(
      x,
      organism = "hsa", 
      keyType = "ncbi-geneid", 
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
    )
  })

# saveRDS(enrich_ls,file.path(Luminal_er_output, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Luminal_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Luminal_er_output, name, "KEGG enriment table.csv"),
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

enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

# export dotplot
enrich_dot_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Luminal_er_output, name, "kegg dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Luminal_er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })
# saveRDS(Luminal, file.path(Luminal_out_dir, "Luminal.rds"))




### refiltering -----------------------------------------------------------------
Idents(Luminal) <- "Spatial.008um_snn_res.0.18"
Luminal2 <- subset(Luminal, idents = c(0:4))
pdf(file.path(Luminal_out_dir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_ADPC_split2.pdf")), height = 6, width = 30)
DimPlot(Luminal2, group.by = "Spatial.008um_snn_res.0.1", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.1') + NoLegend()
DimPlot(Luminal2, group.by = "Spatial.008um_snn_res.0.12", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.12') + NoLegend()
DimPlot(Luminal2, group.by = "Spatial.008um_snn_res.0.15", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.15') + NoLegend()
DimPlot(Luminal2, group.by = "Spatial.008um_snn_res.0.18", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.18') + NoLegend()
DimPlot(Luminal2, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = "umap", split.by = 'Spatial.008um_snn_res.0.2') + NoLegend()
dev.off()

### PCA --------
Luminal2 <- FindVariableFeatures(Luminal2)
Luminal2 <- ScaleData(Luminal2, vars.to.regress = ncount_name)
Luminal2 <- RunPCA(Luminal2)

#### pca inspector ------
pdf(file.path(Luminal_out_dir, 'pca plot_ADPC2.pdf'), height = 6, width = 8)
ElbowPlot(Luminal2, ndims = 50)
dev.off()

### dimension reduction --------
dims <- 10

library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

Luminal2 <- RunUMAP(Luminal2, dims = 1:dims,return.model = TRUE)

reduction.use <- 'pca' # 'harmony' 'pca'


#### cluster -------
Luminal2  <-
  FindNeighbors(
    Luminal2 ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Luminal2 <-
  FindClusters(
    Luminal2,
    resolution = c(0.1, 0.12, 0.15, 0.18, 0.2)
  )


# dimplot
pdf(file.path(Luminal_out_dir, str_glue("Dimplot dim{dims} cluster_ADPC_2.pdf")), height = 6, width = 8)
DimPlot(Luminal2, group.by = "Spatial.008um_snn_res.0.1", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(Luminal2, group.by = "Spatial.008um_snn_res.0.12", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(Luminal2, group.by = "Spatial.008um_snn_res.0.15", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(Luminal2, group.by = "Spatial.008um_snn_res.0.18", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(Luminal2, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = 'umap') + NoLegend()
dev.off()

Luminal2_markers_ls <- 
  list(0.1, 0.12, 0.15, 0.18, 0.2) %>% 
  map(function(res){ # browser()
    resolution <- print(str_glue("Spatial.008um_snn_res.{res}"))
    Idents(Luminal2) <- resolution
    Luminal2_markers <-
      FindAllMarkers(
        JoinLayers(Luminal2),
        test.use = "wilcox",
        logfc.threshold = 0,
        min.pct = 0,
        only.pos = TRUE
      )
    
    Luminal2_markers
  })
names(Luminal2_markers_ls) <- c(0.1, 0.12, 0.15, 0.18, 0.2)
Luminal2_markers_ls$`0.18` %>% openxlsx::write.xlsx(file.path(Luminal2_out_dir, "temp.xlsx"), asTable = T)

Luminal2_top5_markers_ls <- 
  Luminal2_markers_ls %>% 
  map(function(Luminal2_markers){ # browser()
    Luminal2_top5_markers <-
      Luminal2_markers %>% 
      group_by(cluster) %>%
      dplyr::filter(
        p_val_adj < 0.05,
        avg_log2FC > 1,
        pct.1-pct.2 > 0.1
      ) %>%
      # arrange(desc(avg_log2FC)) %>%
      # slice_head(n = 30) %>% 
      arrange(desc((pct.1-pct.2))) %>%
      # slice_head(n = 5) %>% 
      # arrange(desc(pct.1)) %>%
      # slice_head(n = 5) %>%
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
    
    names(Luminal2_top5_markers) <- sub('.*\\.','',names(Luminal2_top5_markers))
    Luminal2_top5_markers  
  })


p_dotplot_ls <-
  Luminal2_top5_markers_ls %>% 
  imap(function(Luminal2_top5_markers, res){
    p <- DotPlot(Luminal2, group.by = str_glue("Spatial.008um_snn_res.{res}"), features = Luminal2_top5_markers) + RotatedAxis()
    p
  })

pdf(file.path(Luminal_out_dir, str_glue("dotplot dims{dims} multiRes_ADPC_2.pdf")), width = 20, height = 6)
p_dotplot_ls %>% 
  iwalk(function(p_dotplot, res){
    print(p_dotplot + ggtitle(res))
  })
dev.off()

DotPlot(Luminal2, group.by = "Spatial.008um_snn_res.0.18", features = major_marker_ls)+RotatedAxis()
DotPlot(Luminal2, group.by = "Spatial.008um_snn_res.0.2", features = major_marker_ls)+RotatedAxis()




## sub cluster ------------------------------------------------------------------
Idents(Luminal2) <- "Spatial.008um_snn_res.0.18"
Luminal2 <- FindSubCluster(Luminal2, cluster = 3, graph.name = "Spatial.008um_snn", algorithm = 1,resolution = 0.6)
Idents(Luminal2) <- "sub.cluster"
#find markers
Luminal2_markers_sub <-
  FindAllMarkers(
    JoinLayers(Luminal2),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
Luminal2_markers_sub %>% 
  openxlsx::write.xlsx(file.path(Luminal_out_dir, "Luminal_sub.xlsx"), 
                       asTable = T)

Luminal2_top5_markers_sub <- 
  Luminal2_markers_sub %>% 
  group_by(cluster) %>%
  dplyr::filter(
                p_val_adj < 0.05,
                avg_log2FC > 1,
                pct.1-pct.2 > 0.1
                ) %>%
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

names(Luminal2_top5_markers_sub) <- sub('.*\\.','',names(Luminal2_top5_markers_sub))

DimPlot(Luminal2, label = TRUE) + NoLegend()
DimPlot(Luminal2, label = TRUE, split.by = "sub.cluster") + NoLegend()
DimPlot(Luminal2, label = TRUE, split.by = "Spatial.008um_snn_res.0.18") + NoLegend()

dims <- 10
pdf(file.path(Luminal_out_dir, str_glue("dotplot dims{dims} major_marker_sub.pdf")), width = 35, height = 6)
DotPlot(Luminal2, group.by = "sub.cluster", features = major_marker_ls)+RotatedAxis()
DotPlot(Luminal2, group.by = "sub.cluster", features = Luminal_marker_ls)+RotatedAxis()
DotPlot(Luminal2, group.by = "sub.cluster", features = Luminal2_top5_markers_sub)+theme_90angle
dev.off()


Idents(Luminal2) <- "Spatial.008um_snn_res.0.18"
#find markers
Luminal2_markers_sub <-
  FindAllMarkers(
    JoinLayers(Luminal2),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
Luminal2_markers_sub %>% 
  openxlsx::write.xlsx(file.path(Luminal_out_dir, "Luminal.xlsx"), 
                       asTable = T)

Luminal2_top5_markers_sub <- 
  Luminal2_markers_sub %>% 
  group_by(cluster) %>%
  dplyr::filter(
    p_val_adj < 0.05,
    # avg_log2FC > 1,
    pct.1-pct.2 > 0.1
  ) %>%
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

names(Luminal2_top5_markers_sub) <- sub('.*\\.','',names(Luminal2_top5_markers_sub))

DimPlot(Luminal2, label = TRUE) + NoLegend()

dims <- 10
pdf(file.path(Luminal_out_dir, str_glue("dotplot dims{dims} major_marker_ADPC_2.pdf")), width = 35, height = 6)
DotPlot(Luminal2, group.by = "Spatial.008um_snn_res.0.18", features = major_marker_ls)+RotatedAxis()
DotPlot(Luminal2, group.by = "Spatial.008um_snn_res.0.18", features = luminal_marker_ls)+RotatedAxis()
DotPlot(Luminal2, group.by = "Spatial.008um_snn_res.0.18", features = Luminal2_top5_markers_sub)+theme_90angle
dev.off()




### annotation ------------------------------------------------------------------
Idents(Luminal2) <- "Spatial.008um_snn_res.0.18"
cell_cluster <- c(
  '0'  = 'cluster0_luminal',
  '1'  = 'cluster1_luminal',
  '2'  = 'CAMK2N1-positive Luminal',
  '3'  = 'Fibroblast-like Luminal',
  '4'  = 'AMACR-positive Luminal'
)

Luminal2 <- RenameIdents(Luminal2, cell_cluster)
Luminal2$subtype <- Idents(Luminal2)

pdf(file.path(Luminal_out_dir, str_glue("Dimplot dim{dims} annotation ADPC.pdf")), height = 6, width = 8)
DimPlot(Luminal2, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

# saveRDS(Luminal2, file.path(Luminal_out_dir, 'Luminal2_ADPC.rds'))




# Annotation info output -------------------------------------------------------
Luminal2$subtype %>% table()
cellType_anno_Luminal2 <- 
  Luminal2@meta.data[c('project', 'subtype')] %>% 
  rownames_to_column('cell_id') %>% 
  mutate(
    Barcode = str_replace(cell_id, '^.*?_', ''), 
    group = subtype
  ) %>% 
  dplyr::select(Barcode, group)

cellType_anno_Luminal2 %>% 
  write.csv(str_glue('result_8um/Luminal_reannotation/cell type Luminal2.csv'), row.names = F)


Luminal2$sub.cluster %>% table()
sub.cluster_anno_Luminal2 <- 
  Luminal2@meta.data[c('project', 'sub.cluster')] %>% 
  rownames_to_column('cell_id') %>% 
  mutate(
    Barcode = str_replace(cell_id, '^.*?_', ''), 
    group = sub.cluster
  ) %>% 
  dplyr::select(Barcode, group)

sub.cluster_anno_Luminal2 %>% 
  write.csv(str_glue('result_8um/Luminal_reannotation/fibroblast sub cluster Luminal2.csv'), row.names = F)




## volcano -------------------------------------------------------------------
Luminal2_markers <-
  FindAllMarkers(
    JoinLayers(Luminal2),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(Luminal2_markers,file.path(Luminal_out_dir, 'Luminal2_diff_gene.csv'))


library(EnhancedVolcano)
type<-unique(Luminal2_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- Luminal2_markers[Luminal2_markers$cluster == type, ]
    p<-EnhancedVolcano(markers,
                       lab = rownames(markers),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = str_glue('{name} volcano plot'),
                       pCutoff = 10e-2,
                       FCcutoff = 0,
                       pointSize = 3.0,
                       labSize = 6.0)
    p
  })
volcano_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(Luminal_out_dir, 'volcano2', str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })


all_markers_major <- list(Luminal2 = Luminal2_markers)
# scatter
p_volcanoplot_major_ls <- 
  all_markers_major %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 1
    dodge_width <- 0.2
    
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
    
    
    
    df_scatter %>%
      ggplot(aes(x = xx, y = y)) +
      geom_rect(
        data = df_bg_rect, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        inherit.aes = FALSE
      ) +
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_point(aes(color = group), size = 1) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
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

p_volcanoplot_major_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(Luminal_out_dir, 'volcano2', str_glue("{name} diffgene dot.pdf")),
      height = 8,
      width = 15
    )
  })




## enrichment -------------------------------------------------------------------
Luminal2_de_genes_for_enrichmenet <- 
  Luminal2_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.1) %>%
  dplyr::filter(abs(avg_log2FC) > 0.25) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")


enrich_ls <-
  Luminal2_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichGO(
      x,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
      ont='ALL',
      pAdjustMethod = 'BH',
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.05,
      keyType = 'SYMBOL'
    )
  })


Luminal2_er_output <- file.path(Luminal_out_dir, 'enrichment') # dir.create(Luminal2_er_output)
# expot enrichment table
enrich_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Luminal2_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Luminal2_er_output, name, "Go enriment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

# export dotplot
enrich_dot_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Luminal2_er_output, name, "go dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Luminal2_er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })

Luminal2_de_genes_id_for_enrichmenet <-
  Luminal2_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::bitr(
      x,
      OrgDb = "org.Hs.eg.db",
      fromType = "SYMBOL",
      toType = "ENTREZID"
    ) %>%
      .$ENTREZID
  }) %>%
  purrr::discard(function(x) length(x) < 1)


enrich_ls <-
  Luminal2_de_genes_id_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichKEGG(
      x,
      organism = "hsa", 
      keyType = "ncbi-geneid", 
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
    )
  })

# saveRDS(enrich_ls,file.path(Luminal2_er_output, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Luminal2_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Luminal2_er_output, name, "KEGG enriment table.csv"),
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

enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

# export dotplot
enrich_dot_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Luminal2_er_output, name, "kegg dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Luminal2_er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })
# saveRDS(Luminal2, file.path(Luminal_out_dir, "Luminal2.rds"))




ADPC8 <- subset(ADPC7, cells = setdiff(Cells(ADPC7), setdiff(Cells(Luminal), Cells(Luminal2))))


## 映射-------------------------------
library(tibble)
library(data.table)
subtype <- Luminal2[['subtype']] %>% rownames_to_column('cells')

meta.data.t <-
  subtype%>%
  right_join(
    ADPC8[['annotation_full']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation4 = ifelse((as.character(annotation_full) %in% c('Luminal cell')),
                         as.character(subtype), 
                         as.character(annotation_full)))%>%
  dplyr::select(cells, annotation4) %>%
  column_to_rownames('cells')


ADPC8 <- AddMetaData(ADPC8,meta.data.t['annotation4'])




# subcluster filtering ----------------------------------------------------------
Luminal3 <- subset(Luminal2, sub.cluster %in% c('0', '1', '2', '4', '3_0', '3_1', '3_2', '3_3', '3_4', '3_5', '3_6'))
Luminal_out_dir <- 'result_8um/Luminal_reannotation3'  # dir.create(Luminal_out_dir)

### PCA --------
Luminal3 <- FindVariableFeatures(Luminal3)
Luminal3 <- ScaleData(Luminal3, vars.to.regress = ncount_name)
Luminal3 <- RunPCA(Luminal3)

#### pca inspector ------
pdf(file.path(Luminal_out_dir, 'pca plot_ADPC.pdf'), height = 6, width = 8)
ElbowPlot(Luminal3, ndims = 50)
dev.off()

### dimension reduction --------
dims <- 10

library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

Luminal3 <- RunUMAP(Luminal3, dims = 1:dims,return.model = TRUE)

reduction.use <- 'pca' # 'harmony' 'pca'


#### cluster -------
Luminal3  <-
  FindNeighbors(
    Luminal3 ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Luminal3 <-
  FindClusters(
    Luminal3,
    resolution = c(0.1, 0.18, 0.2, 0.25, 0.3)
  )


# dimplot
pdf(file.path(Luminal_out_dir, str_glue("Dimplot dim{dims} cluster_ADPC_old.pdf")), height = 6, width = 8)
DimPlot(Luminal3, group.by = "Spatial.008um_snn_res.0.1", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(Luminal3, group.by = "Spatial.008um_snn_res.0.18", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(Luminal3, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(Luminal3, group.by = "Spatial.008um_snn_res.0.25", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(Luminal3, group.by = "Spatial.008um_snn_res.0.3", label = TRUE, reduction = 'umap') + NoLegend()
dev.off()

saveRDS(Luminal3, 'result_8um/Luminal_reannotation2/Luminal3.rds')




## REFILTERING ------------------------------------------------------------------
Luminal4 <- subset(Luminal3, Spatial.008um_snn_res.0.25 %in% c(0:5))

pdf(file.path(Luminal_out_dir, str_glue("Dimplot dim{dims} cluster_ADPC.pdf")), height = 6, width = 8)
DimPlot(Luminal4, group.by = "Spatial.008um_snn_res.0.25", label = TRUE, reduction = 'umap') + NoLegend()
dev.off()

Luminal4_markers_ls <- 
  list(0.25) %>% 
  map(function(res){ # browser()
    resolution <- print(str_glue("Spatial.008um_snn_res.{res}"))
    Idents(Luminal4) <- resolution
    Luminal4_markers <-
      FindAllMarkers(
        JoinLayers(Luminal4),
        test.use = "wilcox",
        logfc.threshold = 0,
        min.pct = 0,
        only.pos = TRUE
      )
    
    Luminal4_markers
  })
names(Luminal4_markers_ls) <- c(0.25)
# Luminal4_markers_ls$`0.2` %>% openxlsx::write.xlsx(file.path(Luminal_out_dir, "temp0.2.xlsx"), asTable = T)
Luminal4_markers_ls$`0.25` %>% openxlsx::write.xlsx(file.path(Luminal_out_dir, "temp0.25.xlsx"), asTable = T)

Luminal4_top5_markers_ls <- 
  Luminal4_markers_ls %>% 
  map(function(Luminal4_markers){ # browser()
    Luminal4_top5_markers <-
      Luminal4_markers %>% 
      group_by(cluster) %>%
      dplyr::filter(
        p_val_adj < 0.05,
        avg_log2FC > 1,
        pct.1-pct.2 > 0.05
      ) %>%
      # arrange(desc(avg_log2FC)) %>%
      # slice_head(n = 30) %>% 
      arrange(desc((pct.1-pct.2))) %>%
      # slice_head(n = 5) %>% 
      # arrange(desc(pct.1)) %>%
      # slice_head(n = 5) %>%
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
    
    names(Luminal4_top5_markers) <- sub('.*\\.','',names(Luminal4_top5_markers))
    Luminal4_top5_markers  
  })


p_dotplot_ls <-
  Luminal4_top5_markers_ls %>% 
  imap(function(Luminal4_top5_markers, res){
    p <- DotPlot(Luminal4, group.by = str_glue("Spatial.008um_snn_res.{res}"), features = Luminal4_top5_markers) + RotatedAxis()
    p
  })

pdf(file.path(Luminal_out_dir, str_glue("dotplot dims{dims} multiRes_ADPC.pdf")), width = 30, height = 6)
p_dotplot_ls %>% 
  iwalk(function(p_dotplot, res){
    print(p_dotplot + ggtitle(res))
  })
dev.off()

#DotPlot(Luminal4, group.by = "Spatial.008um_snn_res.0.2", features = major_marker_ls)+RotatedAxis()
DotPlot(Luminal4, group.by = "Spatial.008um_snn_res.0.25", features = major_marker_ls)+RotatedAxis()


# dotplot
pdf(file.path(Luminal_out_dir, str_glue("dotplot dims{dims} major_marker_ADPC.pdf")), width = 45, height = 10)
DotPlot(Luminal4, group.by = "Spatial.008um_snn_res.0.25", features = major_marker_ls)+RotatedAxis()
DotPlot(Luminal4, group.by = "Spatial.008um_snn_res.0.25", features = luminal_marker_ls)+RotatedAxis()
DotPlot(Luminal4, group.by = "Spatial.008um_snn_res.0.25", features = Luminal4_top5_markers_ls$`0.25`)+RotatedAxis()
dev.off()




## volcano -------------------------------------------------------------------
Idents(Luminal4) <- 'Spatial.008um_snn_res.0.25'
Luminal4_markers <-
  FindAllMarkers(
    JoinLayers(Luminal4),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(Luminal4_markers,file.path(Luminal_out_dir, 'Luminal4_diff_gene.csv'))


library(EnhancedVolcano)
type<-unique(Luminal4_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- Luminal4_markers[Luminal4_markers$cluster == type, ]
    p<-EnhancedVolcano(markers,
                       lab = rownames(markers),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = str_glue('{name} volcano plot'),
                       pCutoff = 10e-2,
                       FCcutoff = 0,
                       pointSize = 3.0,
                       labSize = 6.0)
    p
  })
volcano_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(Luminal_out_dir, 'volcano', str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })


all_markers_major <- list(Luminal4 = Luminal4_markers)
# scatter
p_volcanoplot_major_ls <- 
  all_markers_major %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 1
    dodge_width <- 0.2
    
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
    
    
    
    df_scatter %>%
      ggplot(aes(x = xx, y = y)) +
      geom_rect(
        data = df_bg_rect, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        inherit.aes = FALSE
      ) +
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_point(aes(color = group), size = 1) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
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

p_volcanoplot_major_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(Luminal_out_dir, 'volcano', str_glue("{name} diffgene dot.pdf")),
      height = 8,
      width = 15
    )
  })




## enrichment -------------------------------------------------------------------
Luminal4_de_genes_for_enrichmenet <- 
  Luminal4_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.1) %>%
  dplyr::filter(abs(avg_log2FC) > 0.25) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")


enrich_ls <-
  Luminal4_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichGO(
      x,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
      ont='ALL',
      pAdjustMethod = 'BH',
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.05,
      keyType = 'SYMBOL'
    )
  })


Luminal4_er_output <- file.path(Luminal_out_dir, 'enrichment') # dir.create(Luminal4_er_output)




# expot enrichment table
enrich_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Luminal4_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Luminal4_er_output, name, "Go enriment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

# export dotplot
enrich_dot_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Luminal4_er_output, name, "go dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Luminal4_er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })

Luminal4_de_genes_id_for_enrichmenet <-
  Luminal4_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::bitr(
      x,
      OrgDb = "org.Hs.eg.db",
      fromType = "SYMBOL",
      toType = "ENTREZID"
    ) %>%
      .$ENTREZID
  }) %>%
  purrr::discard(function(x) length(x) < 1)


enrich_ls <-
  Luminal4_de_genes_id_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichKEGG(
      x,
      organism = "hsa", 
      keyType = "ncbi-geneid", 
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
    )
  })

# saveRDS(enrich_ls,file.path(Luminal4_er_output, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Luminal4_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Luminal4_er_output, name, "KEGG enriment table.csv"),
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

enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

# export dotplot
enrich_dot_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Luminal4_er_output, name, "kegg dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Luminal4_er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })
# saveRDS(Luminal4, file.path(Luminal_out_dir, "Luminal4.rds"))




# Annotation info output -------------------------------------------------------
Luminal4$Spatial.008um_snn_res.0.25 %>% table()
cellType_anno_Luminal2 <- 
  Luminal4@meta.data[c('project', 'Spatial.008um_snn_res.0.25')] %>% 
  rownames_to_column('cell_id') %>% 
  mutate(
    Barcode = str_replace(cell_id, '^.*?_', ''), 
    group = Spatial.008um_snn_res.0.25
  ) %>% 
  dplyr::select(Barcode, group)

cellType_anno_Luminal2 %>% 
  write.csv(str_glue('result_8um/Luminal_reannotation3/cell type Luminal4.csv'), row.names = F)




Luminal4$clusters <- 
  Luminal4$Spatial.008um_snn_res.0.25 %>% 
  as.character() 
Luminal4$clusters <- factor(Luminal4$clusters, c(0:5))

saveRDS(Luminal4, 'result_8um/Luminal_reannotation3/Luminal4_final.rds')


  

# Gene scoring ------------------------------------------------------------------
gene_set <- openxlsx::read.xlsx('result_8um/Luminal_reannotation/GSVA Gene Sets for Prostate Cancer Transdifferentiation (40 Sets).xlsx')

pathway <- list()
for (i in 1:nrow(gene_set)) { # browser()
  variable_name <- gene_set[i, 2] %>% str_replace(' ', '_')
  values <- gene_set[i, 4:ncol(gene_set)][!is.na(gene_set[i, 4:ncol(gene_set)])]
  pathway[[variable_name]] <- values
  assign(variable_name, values)
}

pathway_names <- gene_set[[2]] %>% str_replace(' ', '_')
feature_names <- str_replace_all(pathway_names, '$', '1')

n_na_idx <- pathway %>% map_lgl(function(x) all(is.na(x))) %>% {! .}%>%  which
pathway <- pathway[n_na_idx] %>% map(function(x) x[!is.na(x)])

## Luminal4 ---------------------------------------------------------------------
for ( i in seq_along(names(pathway))){
  Luminal4 <- AddModuleScore(JoinLayers(Luminal4), pathway[i], name =  names(pathway)[i],nbin=5,ctrl=10)
}

p <- VlnPlot(Luminal4,
             features = feature_names,
             group.by='Spatial.008um_snn_res.0.25',
             pt.size = 0)

ggsave(
  p,
  filename = file.path(Luminal_out_dir, 'gene set score.pdf'),
  height = 70,
  limitsize = FALSE,
  width  = 30
) 


library(ComplexHeatmap)
plot_data <- Luminal4@meta.data[feature_names]
# 配色方案
col <- c("pink", "lightblue", "green", "yellow", "orange", "#99CCFF")
# col <- c(
#   "#FFCCCC", # 浅红
#   "#CCFFFF", # 浅蓝
#   "#CCFFCC", # 浅绿
#   "#FFFFCC", # 浅黄
#   "#FFCCFF", # 浅紫
#   "#CCCCFF", # 浅粉
#   "#CCFFCC", # 浅青
#   "#FFCC99", # 浅橙
#   "#CCFF99", # 浅绿黄
#   "#99CCFF"  # 浅天蓝
# )
names(col) <- Luminal4$Spatial.008um_snn_res.0.25 %>% as.character() %>% unique() %>% sort()

# 列注释
top_anno <- HeatmapAnnotation(
  subtype = anno_block(
    gp = gpar(fill = col), 
    labels = names(col)
  )
)
# 绘图
pdf(file.path(Luminal_out_dir, 'Heatmap gene set score.pdf'), height = 8, width  = 20)
Heatmap(
  matrix = t(as.matrix(plot_data)), 
  cluster_rows = T, 
  cluster_columns = FALSE, 
  show_column_names = FALSE, 
  show_row_names = T, 
  show_heatmap_legend = F, 
  column_split = Luminal4$Spatial.008um_snn_res.0.25, 
  column_title = NULL, 
  top_annotation = top_anno,
  row_names_max_width = unit(20, "cm")
)
dev.off()




## Luminal ---------------------------------------------------------------------
for ( i in seq_along(names(pathway))){
  Luminal <- AddModuleScore(JoinLayers(Luminal), pathway[i], name =  names(pathway)[i],nbin=5,ctrl=10)
}

p <- VlnPlot(Luminal,
             features = feature_names,
             group.by='Spatial.008um_snn_res.0.18',
             pt.size = 0)

ggsave(
  p,
  filename = file.path(Luminal_out_dir, 'gene set score1.pdf'),
  height = 70,
  limitsize = FALSE,
  width  = 30
) 


library(ComplexHeatmap)
plot_data <- Luminal@meta.data[feature_names]
# 配色方案
col <- c(
  "#FFCCCC", # 浅红
  "#CCFFFF", # 浅蓝
  "#CCFFCC", # 浅绿
  "#FFFFCC", # 浅黄
  "#FFCCFF", # 浅紫
  "#CCCCFF", # 浅粉
  "#CCFFCC", # 浅青
  "#FFCC99", # 浅橙
  # "#CCFF99", # 浅绿黄
  "#99CCFF"  # 浅天蓝
)
names(col) <- unique(Luminal$Spatial.008um_snn_res.0.18)
# 列注释
top_anno <- HeatmapAnnotation(
  Spatial.008um_snn_res.0.18 = anno_block(
    gp = gpar(fill = col), 
    labels = levels(Luminal$Spatial.008um_snn_res.0.18)
  )
)
# 绘图
pdf(file.path(Luminal_out_dir, 'Heatmap gene set score1.pdf'), height = 10, width  = 30)
Heatmap(
  matrix = t(as.matrix(plot_data)), 
  cluster_rows = T, 
  cluster_columns = FALSE, 
  show_column_names = FALSE, 
  show_row_names = T, 
  show_heatmap_legend = F, 
  column_split = Luminal$Spatial.008um_snn_res.0.18, 
  column_title = NULL, 
  top_annotation = top_anno
)
dev.off()




## NE ---------------------------------------------------------------------
NE_out_dir <- "result_8um/NE_subtype/filtered"
for ( i in seq_along(names(pathway))){
  NE <- AddModuleScore(JoinLayers(NE), pathway[i], name =  names(pathway)[i],nbin=5,ctrl=10)
}

p <- VlnPlot(NE,
             features = feature_names,
             group.by='Spatial.008um_snn_res.0.25',
             pt.size = 0)

ggsave(
  p,
  filename = file.path(NE_out_dir, 'gene set score.pdf'),
  height = 70,
  limitsize = FALSE,
  width  = 30
) 


library(ComplexHeatmap)
plot_data <- NE@meta.data[feature_names]
# 配色方案
col <- c("pink", "lightblue", "green", "yellow", "orange")
# col <- c(
#   "#FFCCCC", # 浅红
#   "#CCFFFF", # 浅蓝
#   "#CCFFCC", # 浅绿
#   "#FFFFCC", # 浅黄
#   "#FFCCFF", # 浅紫
#   "#CCCCFF", # 浅粉
#   "#CCFFCC", # 浅青
#   "#FFCC99", # 浅橙
#   "#CCFF99", # 浅绿黄
#   "#99CCFF"  # 浅天蓝
# )
names(col) <- NE$Spatial.008um_snn_res.0.25 %>% as.character() %>% unique() %>% sort()

# 列注释
top_anno <- HeatmapAnnotation(
  subtype = anno_block(
    gp = gpar(fill = col), 
    labels = names(col)
  )
)
# 绘图
pdf(file.path(NE_out_dir, 'Heatmap gene set score.pdf'), height = 8, width  = 20)
Heatmap(
  matrix = t(as.matrix(plot_data)), 
  cluster_rows = T, 
  cluster_columns = FALSE, 
  show_column_names = FALSE, 
  show_row_names = T, 
  show_heatmap_legend = F, 
  column_split = NE$Spatial.008um_snn_res.0.25, 
  column_title = NULL, 
  top_annotation = top_anno,
  row_names_max_width = unit(20, "cm")
)
dev.off()




# NEPC NE subclustering ---------------------------------------------------------
NEPC4 <- readRDS('result_8um/NEPC4_annotation.rds')
NEPC4$annotation_full %>% table()
Idents(NEPC4) <- 'annotation_full'
Idents(NEPC4) %>% table()

NE <- subset(NEPC4, annotation_full == 'NE cell')
NE$annotation_full %>% table()

NE_out_dir <- "result_8um/NE_subtype" # dir.create(NE_out_dir)

### PCA --------
NE <- FindVariableFeatures(NE)
NE <- ScaleData(NE, vars.to.regress = ncount_name)
NE <- RunPCA(NE)

#### pca inspector ------
pdf(file.path(NE_out_dir, 'pca plot_NE.pdf'), height = 6, width = 8)
ElbowPlot(NE, ndims = 50)
dev.off()

### dimension reduction --------
dims <- 10

library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

NE <- RunUMAP(NE, dims = 1:dims,return.model = TRUE)

reduction.use <- 'pca' # 'harmony' 'pca'


#### cluster -------
NE  <-
  FindNeighbors(
    NE ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


NE <-
  FindClusters(
    NE,
    resolution = c(0.1, 0.15, 0.2, 0.25, 0.3)
  )



# dimplot
pdf(file.path(NE_out_dir, str_glue("Dimplot dim{dims} cluster_NEPC.pdf")), height = 6, width = 8)
DimPlot(NE, group.by = "Spatial.008um_snn_res.0.1", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(NE, group.by = "Spatial.008um_snn_res.0.15", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(NE, group.by = "Spatial.008um_snn_res.0.2", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(NE, group.by = "Spatial.008um_snn_res.0.25", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(NE, group.by = "Spatial.008um_snn_res.0.3", label = TRUE, reduction = 'umap') + NoLegend()
dev.off()

pdf(file.path(NE_out_dir, str_glue("Dimplot dim{dims} cluster_NEPC_split.pdf")), height = 6, width = 15)
DimPlot(NE, group.by = "Spatial.008um_snn_res.0.25", label = TRUE, reduction = 'umap', split.by = "Spatial.008um_snn_res.0.25") + NoLegend()
DimPlot(NE, group.by = "Spatial.008um_snn_res.0.3", label = TRUE, reduction = 'umap', split.by = "Spatial.008um_snn_res.0.3") + NoLegend()
dev.off()

NE_markers_ls <- 
  list(0.25) %>% 
  map(function(res){ # browser()
    resolution <- print(str_glue("Spatial.008um_snn_res.{res}"))
    Idents(NE) <- resolution
    NE_markers <-
      FindAllMarkers(
        JoinLayers(NE),
        test.use = "wilcox",
        logfc.threshold = 0,
        min.pct = 0,
        only.pos = TRUE
      )
    
    NE_markers
  })
names(NE_markers_ls) <- c(0.25)
NE_markers_ls$`0.25` %>% openxlsx::write.xlsx(file.path(NE_out_dir, "NE_DEGs.xlsx"), asTable = T)

NE_top5_markers_ls <- 
  NE_markers_ls %>% 
  map(function(NE_markers){ # browser()
    NE_top5_markers <-
      NE_markers %>% 
      group_by(cluster) %>%
      dplyr::filter(
        p_val_adj < 0.05,
        avg_log2FC > 1,
        pct.1-pct.2 > 0.1
      ) %>%
      # arrange(desc(avg_log2FC)) %>%
      # slice_head(n = 30) %>% 
      arrange(desc((pct.1-pct.2))) %>%
      # slice_head(n = 5) %>% 
      # arrange(desc(pct.1)) %>%
      # slice_head(n = 5) %>%
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
    
    names(NE_top5_markers) <- sub('.*\\.','',names(NE_top5_markers))
    NE_top5_markers  
  })


p_dotplot_ls <-
  NE_top5_markers_ls %>% 
  imap(function(NE_top5_markers, res){
    p <- DotPlot(NE, group.by = str_glue("Spatial.008um_snn_res.{res}"), features = NE_top5_markers) + RotatedAxis()
    p
  })

pdf(file.path(NE_out_dir, str_glue("dotplot dims{dims} multiRes_NEPC.pdf")), width = 20, height = 6)
p_dotplot_ls %>% 
  iwalk(function(p_dotplot, res){
    print(p_dotplot + ggtitle(res))
  })
dev.off()




# dotplot
pdf(file.path(NE_out_dir, str_glue("dotplot dims{dims} major_marker_NEPC.pdf")), width = 40, height = 6)
DotPlot(NE, group.by = "Spatial.008um_snn_res.0.25", features = major_marker_ls)+RotatedAxis()
DotPlot(NE, group.by = "Spatial.008um_snn_res.0.25", features = NE_top5_markers_ls$`0.25`)+RotatedAxis()
dev.off()




## sub cluster ------------------------------------------------------------------
Idents(NE) <- "Spatial.008um_snn_res.0.25"
NE <- FindSubCluster(NE, cluster = 3, graph.name = "Spatial.008um_snn", algorithm = 1, resolution = 0.6)
pdf(file.path(NE_out_dir, str_glue('Dimplot dim{dims} sub.cluster_NEPC_res0.6.pdf')), width = 8, height = 6)
DimPlot(NE, group.by = 'sub.cluster', label = T) + NoLegend()
dev.off()

Idents(NE) <- "sub.cluster"
#find markers
NE_markers_sub <-
  FindAllMarkers(
    JoinLayers(NE),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
NE_markers_sub %>% 
  openxlsx::write.xlsx(file.path(NE_out_dir, "NE_sub.xlsx"), 
                       asTable = T)

NE_top5_markers_sub <- 
  NE_markers_sub %>% 
  group_by(cluster) %>%
  dplyr::filter(
    p_val_adj < 0.05,
    avg_log2FC > 1,
    pct.1-pct.2 > 0.1
  ) %>%
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

names(NE_top5_markers_sub) <- sub('.*\\.','',names(NE_top5_markers_sub))

DimPlot(NE, label = TRUE) + NoLegend()
DimPlot(NE, label = TRUE, split.by = "sub.cluster") + NoLegend()
DimPlot(NE, label = TRUE, split.by = "Spatial.008um_snn_res.0.25") + NoLegend()

pdf(file.path(NE_out_dir, str_glue("dotplot dim{dims} major_marker_sub.pdf")), width = 35, height = 6)
DotPlot(NE, group.by = "sub.cluster", features = major_marker_ls)+theme_90angle
DotPlot(NE, group.by = "sub.cluster", features = NE_top5_markers_sub)+theme_90angle
dev.off()




# Annotation info output -------------------------------------------------------
NE$sub.cluster %>% table()
cellType_anno_NE <- 
  NE@meta.data[c('project', 'sub.cluster')] %>% 
  rownames_to_column('cell_id') %>% 
  mutate(
    Barcode = str_replace(cell_id, '^.*?_', ''), 
    group = sub.cluster
  ) %>% 
  dplyr::select(Barcode, group)

cellType_anno_NE %>% 
  write.csv(file.path(NE_out_dir, str_glue('cell type NE.csv')), row.names = F)
saveRDS(NE, 'result_8um/NE_subtype/NE.rds')




## refiltering ------------------------------------------------------------------
NE$sub.cluster %>% table()
NE2 <- subset(NE, sub.cluster %in% c('0', '1', '2', '3_0', '3_1', '3_2', '3_3', '3_4', '3_5', '3_6', '3_7', '4'))
NE2$sub.cluster %>% table()
NE_out_dir <- file.path(NE_out_dir, 'filtered')
dir.create(NE_out_dir)
pdf(file.path(NE_out_dir, str_glue('Dimplot dim{dims} NEPC_refiltered.pdf')), width = 8, height = 6)
DimPlot(NE2, group.by = 'Spatial.008um_snn_res.0.25', label = T) + NoLegend()
DimPlot(NE2, group.by = 'sub.cluster', label = T) + NoLegend()
dev.off()




NE <- NE2
## volcano -------------------------------------------------------------------
Idents(NE) <- 'Spatial.008um_snn_res.0.25'
NE_markers <- 
  FindAllMarkers(
    JoinLayers(NE),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(NE_markers,file.path(NE_out_dir, 'NE_diff_gene.csv'))


library(EnhancedVolcano)
type<-unique(NE_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- NE_markers[NE_markers$cluster == type, ]
    p<-EnhancedVolcano(markers,
                       lab = rownames(markers),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = str_glue('{name} volcano plot'),
                       pCutoff = 10e-2,
                       FCcutoff = 0,
                       pointSize = 3.0,
                       labSize = 6.0)
    p
  })
volcano_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(NE_out_dir, 'volcano', str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })


all_markers_major <- list(NE = NE_markers)
# scatter
p_volcanoplot_major_ls <- 
  all_markers_major %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 1
    dodge_width <- 0.2
    
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
    
    
    
    df_scatter %>%
      ggplot(aes(x = xx, y = y)) +
      geom_rect(
        data = df_bg_rect, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        inherit.aes = FALSE
      ) +
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_point(aes(color = group), size = 1) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
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

p_volcanoplot_major_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(NE_out_dir, 'volcano', str_glue("{name} diffgene dot.pdf")),
      height = 8,
      width = 15
    )
  })




## enrichment -------------------------------------------------------------------
NE_de_genes_for_enrichmenet <- 
  NE_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.15) %>%
  dplyr::filter(abs(avg_log2FC) > 0.25) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")


enrich_ls <-
  NE_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichGO(
      x,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
      ont='ALL',
      pAdjustMethod = 'BH',
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.05,
      keyType = 'SYMBOL'
    )
  })


NE_er_output <- file.path(NE_out_dir, 'enrichment') # dir.create(NE_er_output)
# expot enrichment table
enrich_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(NE_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(NE_er_output, name, "Go enriment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

# export dotplot
enrich_dot_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(NE_er_output, name, "go dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(NE_er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })

NE_de_genes_id_for_enrichmenet <-
  NE_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::bitr(
      x,
      OrgDb = "org.Hs.eg.db",
      fromType = "SYMBOL",
      toType = "ENTREZID"
    ) %>%
      .$ENTREZID
  }) %>%
  purrr::discard(function(x) length(x) < 1)


enrich_ls <-
  NE_de_genes_id_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichKEGG(
      x,
      organism = "hsa", 
      keyType = "ncbi-geneid", 
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
    )
  })

# saveRDS(enrich_ls,file.path(NE_er_output, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(NE_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(NE_er_output, name, "KEGG enriment table.csv"),
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

enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

# export dotplot
enrich_dot_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(NE_er_output, name, "kegg dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(NE_er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })
# saveRDS(NE, file.path(NE_out_dir, "NE_filterd.rds"))




# cnv score cluster cell prop ---------------------------------------------------
## luminal ----------------------------------------------------------------------
cnv_meta <- read.csv('result_8um/Luminal_reannotation3/adata_manual_obs_res0.4.csv')
cnv_meta <- 
  cnv_meta %>% 
  rename(cell_id = X) %>% 
  select(cell_id, cnv_leiden)

Luminal4 <- readRDS('result_8um/Luminal_reannotation3/Luminal4_final.rds')
Luminal_meta<- 
  Luminal4@meta.data %>% 
  rownames_to_column('cell_id') %>% 
  left_join(cnv_meta, by = 'cell_id') %>% 
  column_to_rownames('cell_id')
Luminal4 <- AddMetaData(Luminal4, Luminal_meta['cnv_leiden'])

Luminal4$clusters <- 
  as.character(Luminal4$clusters) %>% 
  str_replace_all(c(
    '0' = 'Seurat_cluster_0', 
    '1' = 'Seurat_cluster_1', 
    '2' = 'Seurat_cluster_2', 
    '3' = 'Seurat_cluster_3', 
    '4' = 'Seurat_cluster_4', 
    '5' = 'Seurat_cluster_5'
  )) %>% 
  factor(levels = c(
    'Seurat_cluster_0', 
    'Seurat_cluster_1', 
    'Seurat_cluster_2', 
    'Seurat_cluster_3', 
    'Seurat_cluster_4', 
    'Seurat_cluster_5'
  ))

Luminal4$cnv_leiden %>% table()
Luminal4$cnv_leiden <- 
  Luminal4$cnv_leiden %>%
  str_replace_all(c(
    '0' = 'cnv_cluster_0', 
    '1' = 'cnv_cluster_1', 
    '2' = 'cnv_cluster_2', 
    '3' = 'cnv_cluster_3', 
    '4' = 'cnv_cluster_4', 
    '5' = 'cnv_cluster_5', 
    '6' = 'cnv_cluster_6', 
    '7' = 'cnv_cluster_7'
  )) %>% 
  factor(levels = c(
    'cnv_cluster_0', 
    'cnv_cluster_1', 
    'cnv_cluster_2', 
    'cnv_cluster_3', 
    'cnv_cluster_4', 
    'cnv_cluster_5', 
    'cnv_cluster_6', 
    'cnv_cluster_7'
  ))

# saveRDS(Luminal4, 'result_8um/Luminal_reannotation3/Luminal4_final_cnv.rds')

Luminal4_meta <- 
  Luminal4@meta.data %>%
  rownames_to_column('cell_id')

p_cnv_col <- 
  Luminal4_meta %>% 
  group_by(cnv_leiden, clusters) %>% 
  summarize(count = n()) %>% 
  mutate(freq = count/sum(count)) %>% 
  
  ggplot(aes(x = cnv_leiden, y = freq, fill = clusters)) + 
  geom_col() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  ggtitle('Luminal seurat cluster propotion by cnv cluster')

ggsave(
  p_cnv_col, 
  filename = 'result_8um/Luminal_reannotation3/Luminal seurat cluster propotion by cnv cluster_res0.4.pdf', 
  height = 6, 
  width = 6
  )


## NE cell ----------------------------------------------------------------------
cnv_meta <- read.csv('result_8um/NE_subtype/filtered/adata_manual_obs_NEPC_res0.5.csv')
cnv_meta <- 
  cnv_meta %>% 
  rename(cell_id = X) %>% 
  select(cell_id, cnv_leiden)
NE <- readRDS('result_8um/NE_subtype/filtered/NE_filterd.rds')
NE_meta <- 
  NE@meta.data %>%
  rownames_to_column('cell_id') %>% 
  left_join(cnv_meta, by = 'cell_id') %>% 
  column_to_rownames('cell_id')
NE <- AddMetaData(NE, NE_meta['cnv_leiden'])
NE$clusters <- 
  as.character(NE$Spatial.008um_snn_res.0.25) %>% 
  str_replace_all(c(
    '0' = 'Seurat_cluster_0', 
    '1' = 'Seurat_cluster_1', 
    '2' = 'Seurat_cluster_2', 
    '3' = 'Seurat_cluster_3', 
    '4' = 'Seurat_cluster_4'
    )) %>% 
  factor(levels = c(
    'Seurat_cluster_0', 
    'Seurat_cluster_1', 
    'Seurat_cluster_2', 
    'Seurat_cluster_3', 
    'Seurat_cluster_4'
  ))

NE$cnv_leiden <- 
  as.character(NE$cnv_leiden) %>% 
  str_replace_all(c(
    '0' = 'cnv_cluster_0', 
    '1' = 'cnv_cluster_1', 
    '2' = 'cnv_cluster_2', 
    '3' = 'cnv_cluster_3', 
    '4' = 'cnv_cluster_4', 
    '5' = 'cnv_cluster_5'
  )) %>% 
  factor(levels = c(
    'cnv_cluster_0', 
    'cnv_cluster_1', 
    'cnv_cluster_2', 
    'cnv_cluster_3', 
    'cnv_cluster_4', 
    'cnv_cluster_5'
  ))
saveRDS(NE, 'result_8um/NE_subtype/filtered/NE_filterd_cnv.rds')

NE_meta <- 
  NE@meta.data %>%
  rownames_to_column('cell_id')

p_col_cnv <- 
  NE_meta %>% 
  group_by(cnv_leiden, clusters) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count)) %>% 
  
  ggplot(aes(x = cnv_leiden, y = freq, fill = clusters)) + 
  geom_col() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  ggtitle('NE cell seurat cluster propotion by cnv cluster')

ggsave(
  p_col_cnv, 
  filename = 'result_8um/NE_subtype/filtered/NE cell seurat cluster propotion by cnv cluster.pdf', 
  height = 6, 
  width = 6
  )
