library(Seurat)
library(tidyverse)

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
hb_genes <- if(orginism == 'hsa') c("HBB", "HBA2") else c("Hbb", "Hba2")

org.db <- if(orginism == 'hsa') "org.Hs.eg.db" else "org.Mm.eg.db"

bin_size <- 16
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
options(future.globals.maxSize = 20 * 1024 ^ 3) # 10Gb


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
        spatial.dir = file.path(path, str_glue("binned_outputs/square_0{bin_size}um/spatial")),
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


# qc filter --------------------------------------------------------------------
qc_metric <- read.csv("qc_metric.csv") %>% map_dfc(str_trim)

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
#saveRDS(visim_ls_pass, 'visim_ls_pass.rds')





# dimension reduction and cluster and annotation -------------------------------
all_data <- merge(visim_ls_pass[[1]], visim_ls_pass[-1])
#SpatialFeaturePlot(all_data, features = c('HBA2', 'HBB', ncount_name))

## normalize for sketch data --------------------
all_data <- NormalizeData(all_data, scale.factor = median(all_data[[ncount_name]] %>% unlist))
all_data <- FindVariableFeatures(all_data, nfeatures = 3000)
all_data <- ScaleData(all_data,vars.to.regress = ncount_name)

all_data <- RunPCA(all_data, features = VariableFeatures(all_data))

pdf(file = file.path(dr_outdir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_data, ndims = 30)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 1:5)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 6:10)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 11:15)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 16:20)
dev.off()

dims <- 20

all_data <- RunUMAP(all_data, dims = 1:dims, return.model = TRUE)

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
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'
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
# DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_split.pdf")), height = 6, width = 30)
DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.2", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.2') + NoLegend()
DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.4", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.4') + NoLegend()
DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.6", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.6') + NoLegend()
DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.8", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.8') + NoLegend()
dev.off()




# reFiltering -------------------------------------------------------------------------------
Idents(all_data) <- "Spatial.016um_snn_res.0.2"
all_data2 <- subset(all_data, idents = c(0:6))

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data2, group.by = "Spatial.016um_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data2, group.by = "Spatial.016um_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data2, group.by = "Spatial.016um_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data2, group.by = "Spatial.016um_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_split.pdf")), height = 6, width = 30)
DimPlot(all_data2, group.by = "Spatial.016um_snn_res.0.2", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.2') + NoLegend()
DimPlot(all_data2, group.by = "Spatial.016um_snn_res.0.4", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.4') + NoLegend()
DimPlot(all_data2, group.by = "Spatial.016um_snn_res.0.6", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.6') + NoLegend()
DimPlot(all_data2, group.by = "Spatial.016um_snn_res.0.8", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.8') + NoLegend()
dev.off()

Idents(all_data2) <- "Spatial.016um_snn_res.0.4"
all_data3 <- subset(all_data2, idents = c(0:9, 11, 12))

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data3, group.by = "Spatial.016um_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data3, group.by = "Spatial.016um_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data3, group.by = "Spatial.016um_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data3, group.by = "Spatial.016um_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_split.pdf")), height = 6, width = 30)
DimPlot(all_data3, group.by = "Spatial.016um_snn_res.0.2", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.2') + NoLegend()
DimPlot(all_data3, group.by = "Spatial.016um_snn_res.0.4", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.4') + NoLegend()
DimPlot(all_data3, group.by = "Spatial.016um_snn_res.0.6", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.6') + NoLegend()
DimPlot(all_data3, group.by = "Spatial.016um_snn_res.0.8", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.8') + NoLegend()
dev.off()

Idents(all_data3) <- "Spatial.016um_snn_res.0.6"
all_data4 <- subset(all_data3, idents = c(0:13))

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data4, group.by = "Spatial.016um_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data4, group.by = "Spatial.016um_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data4, group.by = "Spatial.016um_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data4, group.by = "Spatial.016um_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_split.pdf")), height = 6, width = 30)
DimPlot(all_data4, group.by = "Spatial.016um_snn_res.0.2", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.2') + NoLegend()
DimPlot(all_data4, group.by = "Spatial.016um_snn_res.0.4", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.4') + NoLegend()
DimPlot(all_data4, group.by = "Spatial.016um_snn_res.0.6", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.6') + NoLegend()
DimPlot(all_data4, group.by = "Spatial.016um_snn_res.0.8", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.8') + NoLegend()
dev.off()

Idents(all_data4) <- "Spatial.016um_snn_res.0.8"
all_data5 <- subset(all_data4, idents = c(0:15, 17, 19:21))

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data5, group.by = "Spatial.016um_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data5, group.by = "Spatial.016um_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data5, group.by = "Spatial.016um_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data5, group.by = "Spatial.016um_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_split.pdf")), height = 6, width = 30)
DimPlot(all_data5, group.by = "Spatial.016um_snn_res.0.2", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.2') + NoLegend()
DimPlot(all_data5, group.by = "Spatial.016um_snn_res.0.4", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.4') + NoLegend()
DimPlot(all_data5, group.by = "Spatial.016um_snn_res.0.6", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.6') + NoLegend()
DimPlot(all_data5, group.by = "Spatial.016um_snn_res.0.8", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.8') + NoLegend()
dev.off()

saveRDS(all_data5, 'all_data5.rds')


### normalize data --------------------
all_data <- readRDS('all_data5.rds')
all_data <- NormalizeData(all_data, scale.factor = median(all_data[[ncount_name]] %>% unlist))
all_data <- FindVariableFeatures(all_data, nfeatures = 2000)
all_data <- ScaleData(all_data,vars.to.regress = ncount_name)

all_data <- RunPCA(all_data, features = VariableFeatures(all_data))

pdf(file = file.path(dr_outdir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_data, ndims = 30)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 1:5)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 6:10)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 11:15)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 16:20)
dev.off()

dims <- 20

all_data <- RunUMAP(all_data, dims = 1:dims, return.model = TRUE)

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
pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_split.pdf")), height = 6, width = 30)
DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.2", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.2') + NoLegend()
DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.4", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.4') + NoLegend()
DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.6", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.6') + NoLegend()
DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.8", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.8') + NoLegend()
dev.off()

Idents(all_data) <- "Spatial.016um_snn_res.0.2"
all_data6 <- subset(all_data, idents = c(0:7))

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data6, group.by = "Spatial.016um_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data6, group.by = "Spatial.016um_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data6, group.by = "Spatial.016um_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data6, group.by = "Spatial.016um_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_split.pdf")), height = 6, width = 30)
DimPlot(all_data6, group.by = "Spatial.016um_snn_res.0.2", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.2') + NoLegend()
DimPlot(all_data6, group.by = "Spatial.016um_snn_res.0.4", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.4') + NoLegend()
DimPlot(all_data6, group.by = "Spatial.016um_snn_res.0.6", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.6') + NoLegend()
DimPlot(all_data6, group.by = "Spatial.016um_snn_res.0.8", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.8') + NoLegend()
dev.off()

Idents(all_data6) <- "Spatial.016um_snn_res.0.4"
all_data7 <- subset(all_data6, idents = c(0:10))

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data7, group.by = "Spatial.016um_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data7, group.by = "Spatial.016um_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data7, group.by = "Spatial.016um_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data7, group.by = "Spatial.016um_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_split.pdf")), height = 6, width = 30)
DimPlot(all_data7, group.by = "Spatial.016um_snn_res.0.2", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.2') + NoLegend()
DimPlot(all_data7, group.by = "Spatial.016um_snn_res.0.4", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.4') + NoLegend()
DimPlot(all_data7, group.by = "Spatial.016um_snn_res.0.6", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.6') + NoLegend()
DimPlot(all_data7, group.by = "Spatial.016um_snn_res.0.8", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.8') + NoLegend()
dev.off()


saveRDS(all_data7, 'all_data7.rds')


### normalize data --------------------
all_data <- readRDS('all_data7.rds')
all_data <- NormalizeData(all_data, scale.factor = median(all_data[[ncount_name]] %>% unlist))
all_data <- FindVariableFeatures(all_data, nfeatures = 2000)
all_data <- ScaleData(all_data,vars.to.regress = ncount_name)

all_data <- RunPCA(all_data, features = VariableFeatures(all_data))

pdf(file = file.path(dr_outdir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_data, ndims = 30)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 1:5)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 6:10)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 11:15)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 16:20)
dev.off()

dims <- 20

all_data <- RunUMAP(all_data, dims = 1:dims, return.model = TRUE)

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
pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}_split.pdf")), height = 6, width = 30)
DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.2", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.2') + NoLegend()
DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.4", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.4') + NoLegend()
DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.6", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.6') + NoLegend()
DimPlot(all_data, group.by = "Spatial.016um_snn_res.0.8", label = TRUE, reduction = "harmonyumap", split.by = 'Spatial.016um_snn_res.0.8') + NoLegend()
dev.off()




# dotplot ------------------------------------------------------------------------
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

major_marker_ls <- major_marker_ls[major_marker_ls %in% rownames(LayerData(JoinLayers(all_data), layer = 'data'))]
pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} main_marker.pdf")), width = 55, height = 8)
DotPlot(all_data, group.by = "Spatial.016um_snn_res.0.2", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(all_data, group.by = "Spatial.016um_snn_res.0.4", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(all_data, group.by = "Spatial.016um_snn_res.0.6", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(all_data, group.by = "Spatial.016um_snn_res.0.8", features = major_marker_ls, dot.scale = 10) + theme_90angle
dev.off()

pdf(
  file.path(dr_outdir, str_glue("featureplot dims{dims} {reduction.use} main_marker.pdf")), 
  width = 25, 
  height = 100
)
FeaturePlot(
  all_data,
  features = major_marker_ls,
  ncol = 5,
  reduction = "harmonyumap",
  raster = TRUE
)  
dev.off()

Idents(all_data) <- "Spatial.016um_snn_res.0.8"
all_markers <-
  FindAllMarkers(
    JoinLayers(all_data),
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
DotPlot(all_data, group.by = "Spatial.016um_snn_res.0.8", features = major_marker_ls)+RotatedAxis()
DotPlot(all_data, group.by = "Spatial.016um_snn_res.0.8", features = all_top5_markers)+RotatedAxis()
dev.off()

# Annotaiton -------------------------------------------------------------------------------
Idents(all_data) <- "Spatial.016um_snn_res.0.8"
cell_cluster <- c(
  '0' = 'Plasma cells', #IGKC
  '1' = 'NE cells', 
  '2' = 'Luminal cells', # KLK2
  '3' = 'NE cells',
  '4' = 'NE cells',
  '5' = 'Club cells', # MMP7
  '6' = 'Luminal cells', 
  '7' = 'Luminal cells',
  '8' = 'Smooth muscle cells', # TAGLN DES
  '9' = 'Endothelials', # VMF PECAM1
  '10' = 'Fibroblasts', # COL1A1
  '11' = 'Luminal cells', #NPY https://advanced.onlinelibrary.wiley.com/doi/full/10.1002/advs.202305724 tumor
  '12' = 'Macrophages', # CD68
  '13' = 'NE cells', # CHGB
  '14' = 'Luminal cells', #CRISP3 https://erc.bioscientifica.com/view/journals/erc/27/7/ERC-20-0092.xml tumor
  '15' = 'Luminal cells', #IGFBP2 https://www.spandidos-publications.com/10.3892/or.2020.7892 tumor
  '16' = 'Luminal cells', # KLK3
  '17' = 'Luminal cells', 
  '18' = 'Luminal cells', # AMACR
  '19' = 'Luminal cells' # REG4 https://pmc.ncbi.nlm.nih.gov/articles/PMC9511136/ tumor
)

all_data <- RenameIdents(all_data, cell_cluster)
all_data$annotation <- Idents(all_data)

pdf(file.path(dr_outdir, str_glue('Dimplot dim{dims} annotation {reduction.use}.pdf')), width = 8, height = 6)
DimPlot(all_data, group.by = 'annotation', reduction = 'harmonyumap', label = T) + NoLegend()
dev.off()

cellType_markers <-
  FindAllMarkers(
    JoinLayers(all_data),
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


major_marker2 <- c(
  "Plasma cells" = "IGHG1",
  "Plasma cells" = "IGHA1",
  "Plasma cells" = "IGHG3",
  "Plasma cells" = "IGKC",
  
  # 'NE cells' = 'ASCL-1',
  # 'NE cells' = 'FoxA2',
  'NE cells' = 'CHGA',
  'NE cells' = 'CHGB',
  'NE cells' = 'SYP',
  'NE cells' = 'NKX2-1',
  'NE cells' = 'ENO2',
  # 'NE cells' = 'MYCN',
  # 'NE cells' = 'POUF3F2',
  # 'NE cells' = 'LMO3',
  'NE cells' = 'INSM1',
  'NE cells' = 'SOX2',
  'NE cells' = 'EZH2',
  # 'NE cells' = 'SIAH2',
  
  'Luminal cells' = 'KRT8',
  # 'Luminal cells' = 'KRT18',
  'Luminal cells' = 'NKX3-1',
  'Luminal cells' = 'KLK4',
  'Luminal cells' = 'KLK3',
  'Luminal cells' = 'KLK2',
  'Luminal cells' = 'MSMB',
  'Luminal cells' = 'ACPP',
  
  # 'Club cells' = 'SCGB1A3',
  # 'Club cells' = 'WFDC2',
  'Club cells' = 'LCN2',
  'Club cells' = 'MMP7',	
  'Club cells' = 'KRT4',	
  'Club cells' = 'TACSTD2',
  'Club cells' = 'SCGB3A1',
  'Club cells' = 'LTF',
  'Club cells' = 'RARRES1',
  'Club cells' = 'KRT7',
  'Club cells' = 'AGR2',
  
  # 'Smooth muscle cells' = 'ACTA2',
  'Smooth muscle cells' = 'MYH11',
  'Smooth muscle cells' = 'TAGLN',
  'Smooth muscle cells' = 'CNN1',
  'Smooth muscle cells' = 'DES',
  'Smooth muscle cells' = 'ACTG2',
  
  'Endothelials' = 'PECAM1', 
  # 'Endothelials' = 'ACKR1', 
  # 'Endothelials' = 'CLDN5',
  # 'Endothelials' = 'RAMP2',
  # 'Endothelials' = 'TM4SF1',
  'Endothelials' = 'RNASE1',
  # 'Endothelials' = 'EGFL7',
  'Endothelials' = 'RAMP3',
  'Endothelials' = 'PLVAP',
  # 'Endothelials' = 'AQP1',
  # 'Endothelials' = 'ECSCR',
  # 'Endothelials' = 'FKBP1A',
  # 'Endothelials' = 'AC011526.1',
  # 'Endothelials' = 'EMP1',
  # 'Endothelials' = 'DARC',
  'Endothelials' = 'VWF',
  'Endothelials' = 'EMCN',
  
  'Fibroblasts' = 'DCN',
  'Fibroblasts' = 'LUM',
  'Fibroblasts' = 'VIM',
  'Fibroblasts' = 'IGFBP7',
  'Fibroblasts' = 'COL1A1',
  'Fibroblasts' = 'COL1A2',
  
  "Macrophages" = "CD68",
  "Macrophages" = "CD163",
  "Macrophages" = "C1QA",
  "Macrophages" = "C1QC", 
  'Macrophages'  = 'CSF1R',
  'Macrophages' = 'C1QB'
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2.pdf")), width = 25, height = 10)
DotPlot(all_data, group.by = "annotation", features = major_marker2)+ RotatedAxis() + theme_90angle
DotPlot(all_data, group.by = "annotation", features = cellType_top5_markers)+ RotatedAxis() + theme_90angle
dev.off()

# saveRDS(all_data, 'all_data_mainCellType.rds')




# Luminal subype ---------------------------------------------------------------
luminal <- subset(all_data, annotation == 'Luminal cells')
luminal_out_dir <- "Luminal_subtype" # dir.create(luminal_out_dir)

## PCA --------
luminal <- FindVariableFeatures(luminal)
luminal <- ScaleData(luminal, vars.to.regress = ncount_name)
luminal <- RunPCA(luminal)

### pca inspector ------
pdf(file.path(luminal_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(luminal, ndims = 50)
dev.off()

pdf(file.path(luminal_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(luminal, dims = 1:5)
dev.off()
pdf(file.path(luminal_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(luminal, dims = 6:10)
dev.off()
pdf(file.path(luminal_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(luminal, dims = 11:15)
dev.off()
pdf(file.path(luminal_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(luminal, dims = 16:20)
dev.off()

## dimension reduction --------
dims <- 20

library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

luminal <- RunUMAP(luminal, dims = 1:dims,return.model = TRUE)

pdf(file.path(luminal_out_dir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
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

pdf(file.path(luminal_out_dir, str_glue("Dimplot dim{dims} project - harmony.pdf")), height = 6, width = 8)
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
pdf(file.path(luminal_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(luminal, group.by = "Spatial.016um_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(luminal, group.by = "Spatial.016um_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(luminal, group.by = "Spatial.016um_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(luminal, group.by = "Spatial.016um_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

Idents(luminal) <- "Spatial.016um_snn_res.0.6"

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
    JoinLayers(luminal),
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
pdf(file.path(luminal_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 45, height = 10)
DotPlot(luminal, group.by = "Spatial.016um_snn_res.0.6", features = luminal_marker_ls)+RotatedAxis()
DotPlot(luminal, group.by = "Spatial.016um_snn_res.0.6", features = major_marker_ls)+RotatedAxis()
DotPlot(luminal, group.by = "Spatial.016um_snn_res.0.6", features = luminal_top5_markers)+RotatedAxis()
dev.off()




## reFiltering ----------------------------------------------------------------
luminal2 <- subset(luminal, idents = c(0:22, 26))
luminal_raw <- luminal
luminal <- luminal2

## PCA --------
luminal <- FindVariableFeatures(luminal)
luminal <- ScaleData(luminal, vars.to.regress = ncount_name)
luminal <- RunPCA(luminal)

### pca inspector ------
pdf(file.path(luminal_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(luminal, ndims = 50)
dev.off()

pdf(file.path(luminal_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(luminal, dims = 1:5)
dev.off()
pdf(file.path(luminal_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(luminal, dims = 6:10)
dev.off()
pdf(file.path(luminal_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(luminal, dims = 11:15)
dev.off()
pdf(file.path(luminal_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(luminal, dims = 16:20)
dev.off()

## dimension reduction --------
dims <- 20

luminal <- RunUMAP(luminal, dims = 1:dims,return.model = TRUE)

pdf(file.path(luminal_out_dir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
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

pdf(file.path(luminal_out_dir, str_glue("Dimplot dim{dims} project - harmony.pdf")), height = 6, width = 8)
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
pdf(file.path(luminal_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(luminal, group.by = "Spatial.016um_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(luminal, group.by = "Spatial.016um_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(luminal, group.by = "Spatial.016um_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(luminal, group.by = "Spatial.016um_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

Idents(luminal) <- "Spatial.016um_snn_res.0.6"

luminal_markers <-
  FindAllMarkers(
    JoinLayers(luminal),
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
pdf(file.path(luminal_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 45, height = 10)
DotPlot(luminal, group.by = "Spatial.016um_snn_res.0.6", features = luminal_marker_ls)+RotatedAxis()
DotPlot(luminal, group.by = "Spatial.016um_snn_res.0.6", features = major_marker_ls)+RotatedAxis()
DotPlot(luminal, group.by = "Spatial.016um_snn_res.0.6", features = luminal_top5_markers)+RotatedAxis()
dev.off()




## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = 'Normal/Normal-like luminal',
  '1'  = 'Normal/Normal-like luminal', # AZGP1
  '2'  = 'Adeno-malignant luminal',
  '3'  = 'Adeno-malignant luminal',
  '4'  = 'Normal/Normal-like luminal',
  '5'  = 'Adeno-malignant luminal',
  '6'  = 'Normal/Normal-like luminal',
  '7'  = 'Adeno-malignant luminal',
  '8'  = 'Adeno-malignant luminal',
  '9'  = 'Adeno-malignant luminal',
  '10' = 'Normal/Normal-like luminal',
  '11'  = 'Plasma cells',
  '12'  = 'Adeno-malignant luminal',
  '13'  = 'Normal/Normal-like luminal', # ACPP
  '14'  = 'NE cells',
  '15'  = 'Adeno-malignant luminal',
  '16'  = 'Adeno-malignant luminal',
  '17'  = 'Adeno-malignant luminal',
  '18'  = 'Adeno-malignant luminal',
  '19'  = 'Normal/Normal-like luminal',
  '20'  = 'Adeno-malignant luminal',
  '21'  = 'Normal/Normal-like luminal',
  '22'  = 'Adeno-malignant luminal',
  '23'  = 'Normal/Normal-like luminal'
)

# Idents(luminal) <- "Spatial.016um_snn_res.0.6"

luminal <- RenameIdents(luminal, cell_cluster)
luminal$subtype <- Idents(luminal)

pdf(file.path(luminal_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(luminal, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(luminal_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 7, width = 15)
DimPlot(luminal, group.by = "subtype",split.by = 'project', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

#saveRDS(luminal, file.path(luminal_out_dir, 'luminal.rds'))

all_data2 <- subset(all_data, cells = setdiff(Cells(all_data), Cells(subset(luminal_raw, idents = c(23:25)))))
all_data_raw <- all_data
all_data <- all_data2




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
  'Adeno-malignant luminal' = 'AMACR',
  'Adeno-malignant luminal' = 'COL9A2',
  'Adeno-malignant luminal' = 'KLK4',
  'Adeno-malignant luminal' = 'KLK3',
  'Adeno-malignant luminal' = 'KLK2',
  
  'Club cells' = 'MMP7',	
  'Club cells' = 'KRT4',	
  'Club cells' = 'TACSTD2',
  'Club cells' = 'SCGB3A1',
  'Club cells' = 'LTF',
  'Club cells' = 'RARRES1',
  'Club cells' = 'KRT7',
  'Club cells' = 'AGR2',
  'Endothelials' = 'PECAM1', 
  # 'Endothelials' = 'ACKR1', 
  # 'Endothelials' = 'CLDN5',
  # 'Endothelials' = 'RAMP2',
  # 'Endothelials' = 'TM4SF1',
  'Endothelials' = 'RNASE1',
  # 'Endothelials' = 'EGFL7',
  'Endothelials' = 'RAMP3',
  'Endothelials' = 'PLVAP',
  # 'Endothelials' = 'AQP1',
  # 'Endothelials' = 'ECSCR',
  # 'Endothelials' = 'FKBP1A',
  # 'Endothelials' = 'AC011526.1',
  # 'Endothelials' = 'EMP1',
  # 'Endothelials' = 'DARC',
  'Endothelials' = 'VWF',
  'Endothelials' = 'EMCN',
  
  'Fibroblasts' = 'DCN',
  'Fibroblasts' = 'LUM',
  'Fibroblasts' = 'VIM',
  'Fibroblasts' = 'IGFBP7',
  'Fibroblasts' = 'COL1A1',
  'Fibroblasts' = 'COL1A2',
  
  "Macrophages" = "CD68",
  "Macrophages" = "CD163",
  "Macrophages" = "C1QA",
  "Macrophages" = "C1QC", 
  'Macrophages'  = 'CSF1R',
  'Macrophages' = 'C1QB',
  
  # 'NE cells' = 'ASCL-1',
  # 'NE cells' = 'FoxA2',
  'NE cells' = 'CHGA',
  'NE cells' = 'CHGB',
  'NE cells' = 'SYP',
  'NE cells' = 'NKX2-1',
  'NE cells' = 'ENO2',
  # 'NE cells' = 'MYCN',
  # 'NE cells' = 'POUF3F2',
  # 'NE cells' = 'LMO3',
  'NE cells' = 'INSM1',
  'NE cells' = 'SOX2',
  'NE cells' = 'EZH2',
  # 'NE cells' = 'SIAH2',
  
  'Normal/Normal-like luminal' = 'ACPP',
  'Normal/Normal-like luminal' = 'MSMB',
  'Normal/Normal-like luminal' = 'RDH11',
  'Normal/Normal-like luminal' = 'AZGP1', 
  
  "Plasma cells" = "IGHG1",
  "Plasma cells" = "IGHA1",
  "Plasma cells" = "IGHG3",
  "Plasma cells" = "IGKC",
  
  # 'Smooth muscle cells' = 'ACTA2',
  'Smooth muscle cells' = 'MYH11',
  'Smooth muscle cells' = 'TAGLN',
  'Smooth muscle cells' = 'CNN1',
  'Smooth muscle cells' = 'DES',
  'Smooth muscle cells' = 'ACTG2'
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

pdf(file.path(luminal_out_dir, "dotplot celltype main_marker.pdf"), width = 25, height = 8)
DotPlot(all_data, group.by = "annotation2", features = luminal_cellTye_markers, dot.scale = 10) + theme_90angle
dev.off()


# saveRDS(all_data, 'all_data_celltype.rds')




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
  filename = file.path(dr_outdir, "cell_freq_project.pdf"),
  height = 6,
  width  = 6
)

ggsave(
  p_box_group,
  filename = file.path(dr_outdir, "cell_freq_group.pdf"),
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
  filename = file.path(dr_outdir, "group_freq_cellType.pdf"),
  height = 6,
  width  = 10
)




## spatial plot ------
ggsave(
  SpatialDimPlot(all_data, group.by = 'annotation2', pt.size.factor = 6, image.scale = 'hires'),
  filename = file.path(sp_output, "Saptial_cell_type2.pdf"),
  height = 10,
  width  = 20
)


# Annotation info output -------------------------------------------------------
all_data@meta.data['annotation2'] %>% 
  rownames_to_column('Barcode') %>% 
  write.csv('annotation.csv', row.names = F)

