library(Seurat)
library(tidyverse)

# setting dir -------
img_output <- "materials_for_rmd" # dir.create(img_output)
qc_outdir <- "1. QC" # dir.create(qc_outdir)
dr_outdir <- "2. dimension reduction and cluster"  # dir.create(dr_outdir)
sp_output <- "3. spatial plot" # dir.create(sp_output)
mk_output <- "4. marker gene" # dir.create(mk_output)
er_output <- "5. enriment analysis" # dir.create(er_output)

dir.create(img_output)
dir.create(qc_outdir)
dir.create(dr_outdir)
dir.create(sp_output)
dir.create(mk_output)
dir.create(er_output)

# 4 11 20 hvg5000
# 12 hvg2000 regression out H cell


# set future -------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 20 * 1024 ^ 3) # 10Gb


# read data -------
sample_df <- 
  # read.csv("data/Annoation.csv")
  data.frame(
    sample = c('4M-1','4M-2','24M-1','24M-2'),
    group  = c('4M-7','4M-7','24M-16','24M-16')
  )

file_list <- 
  list.files(
    "/mnt/NexData_file/nas/outputdata/spaceranger_out/FWSP20250661_new",
    pattern = paste(sample_df$group, collapse = "|"), 
    full.names = TRUE
  ) %>%
  str_subset("\\..*$", negate = TRUE)

stopifnot(nrow(file_list) == nrow(sample_df))

# TODO: fix bugs
# Samples order in file_list may not identical to sample_df
# names(file_list) <- sample_df$group # basename(file_list)
names(file_list) <- sample_df$group[match(basename(file_list), sample_df$group)]


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
      bin.size = 8, 
      slice = name, 
      image = Read10X_Image(
        file.path(path, "spatial"),
        spatial.dir = file.path(path, "binned_outputs/square_008um/spatial"),
        assay = 'Spatial.008um'
      )
    )
  })




# parse ROI----------------
roi_ls <- list.files(
  'ROI_2',
  full.names = TRUE
) 

names(roi_ls) <- 
  roi_ls %>%
  basename() %>%
  str_remove("\\..*")


roi_coords_ls <-
  roi_ls %>%
  map(function(path){
    dat <- read.csv(path)
    colnames(dat) <- dat %>% colnames() %>% str_to_lower()
    dat<-dat %>% mutate(cell=barcode)%>%
      column_to_rownames('barcode')
    dat
  }) 


visium_roi_ls <-
  roi_coords_ls %>%
  imap(function(roi, name){ # browser()
    # get Visium Object
    visium_name <- sample_df$group[sample_df$sample == name]
    visium <- visim_ls[[visium_name]]
    
    # subset visium
    cells <- roi$cell
    visium_roi <- subset(visium, cells = cells)
    
    visium_roi
  })


# add sample_id or sample_name
visium_roi_ls <-
  visium_roi_ls %>%
  imap(function(visium, name){
    visium$project <- name
    Idents(visium) <- 'project'
    RenameCells(visium, add.cell.id = name)
  })

# visium_roi override visium
visim_ls <- visium_roi_ls


visim_ls <- 
  visim_ls %>%
  map(~ {DefaultAssay(.x) <- 'Spatial.008um'; .x}) %>%
  map(NormalizeData) 







# # read coordinate ------------
# visim_coords_ls <-
#   file_list %>% 
#   imap(function(path, name){ 
#     proj <- basename(path)
#     coord_file_name <- paste0("coords_", proj, ".csv")
#     
#     df <- read.csv(file.path("data", coord_file_name), header = TRUE)
#     df$Barcode
#   })



# qc -----
visim_ls <- 
  visim_ls %>%
  map(function(visium){
    subset(visium, nCount_Spatial.008um > 0)
  })

visim_ls <- 
  visim_ls %>%
  map(function(x){# browser()
    x[["percent.mt"]]    <- PercentageFeatureSet(x, pattern = "^mt-")
    x[["percent.ribo"]]  <- PercentageFeatureSet(x, pattern = "^Rp[SL]")
    x$log10GenesPerUMI   <- log10(x$nFeature_Spatial.008um) / log10(x$nCount_Spatial.008um)
    x
  })


p_qc_ls <-
  visim_ls %>%
  imap(function(x, name){ # browser()
    # only subset 10,000 cells
    x_sub <- subset(x, cells = Cells(x[['Spatial.008um']]), downsample = 100000, seed = 1)
    p <- VlnPlot(
      x_sub, 
      features = c("nFeature_Spatial.008um", "nCount_Spatial.008um", "percent.mt", "log10GenesPerUMI"), 
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
    x_sub <- subset(x, cells = Cells(x[['Spatial.008um']]), downsample = 100000, seed = 1)
    p <- VlnPlot(
      x_sub, 
      features = c("nFeature_Spatial.008um", "nCount_Spatial.008um", "percent.mt", "log10GenesPerUMI"), 
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
  "nFeature_Spatial.008um", "nCount_Spatial.008um",
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
    
    res <- qc_ls %>% purrr::discard(function(x) length(x) < 1) %>% str_c(collapse = " & ")
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
    plot1 <- FeatureScatter(all_data, feature1 = "nCount_Spatial.008um", feature2 = "nFeature_Spatial.008um") +
      patchwork::plot_annotation(title = sample)
    plot2 <- FeatureScatter(all_data, feature1 = "nCount_Spatial.008um", feature2 = "percent.mt") +
      patchwork::plot_annotation(title = sample)
    plot3 <- FeatureScatter(all_data, feature1 = "nCount_Spatial.008um", feature2 = "log10GenesPerUMI") +
      patchwork::plot_annotation(title = sample)
    p1 <- plot1 + plot2 + plot3 + patchwork::plot_layout(ncol = 3)
    
    # all_data_pass
    plot1_pass <- FeatureScatter(all_data_pass, feature1 = "nCount_Spatial.008um", feature2 = "nFeature_Spatial.008um")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    plot2_pass <- FeatureScatter(all_data_pass, feature1 = "nCount_Spatial.008um", feature2 = "percent.mt")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    plot3_pass <- FeatureScatter(all_data_pass, feature1 = "nCount_Spatial.008um", feature2 = "log10GenesPerUMI")+
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


# # ROI分析 ---------------------
# roi_out_dir <- "roi_analysis" # dir.create(roi_out_dir)
# 
# all_data <- readRDS("all_data_projected.rds")
# 
# roi_df <-
#   list.files('ROI', full.names = TRUE, pattern = "csv") %>%
#   # setNames(., str_remove_all(basename(.), "^.*?-|-.*?$")) %>%
#   setNames(., str_extract(., "(?<=ROI/)[^.]+")) %>% 
#   map(~ read.csv(.x, col.names = c('Barcode', 'ROI'))) %>%
#   enframe("group") %>%
#   unnest(value) %>%
#   mutate(Barcode = str_c(group, Barcode, sep = "_")) %>%
#   # mutate(roi_group = str_c(group, ROI, sep = "_")) %>%
#   dplyr::select(group, Barcode) %>%
#   column_to_rownames("Barcode")
# 
# # roi_ls <- 
# #   split(
# #     rownames(roi_df),
# #     roi_df$roi_group
# #   )
# 
# 
# roi_df <- 
#   read.csv('ROI/roi.csv') %>%
#   dplyr::filter(group != "")
# 
# roi_group <- matrix(
#   c(
#     'KB2209859', 'Infection',
#     'KB2108287', 'Drug',
#     'KB2319088', 'NSAIDs',
#     'KB2305877', 'Diesel oil'
#   ),
#   ncol  = 2,
#   byrow = TRUE
# ) %>% 
#   as.data.frame() %>%
#   setNames(c("sample_id", "sample_name"))
# 
# meta_roi_df <- 
#   roi_df %>%
#   left_join(
#     roi_group,
#     by = c(group = 'sample_id')
#   ) %>%
#   mutate(Barcode = str_c("sample", Barcode, sep = "_")) %>%
#   column_to_rownames("Barcode")
# 
# all_data <-
#   all_data %>%
#   AddMetaData(meta_roi_df)
# 
# all_data <- all_data %>% subset(cells = rownames(meta_roi_df))
# 




# dimension reduction and cluster and annotation ------
# all_data <- merge(visim_ls_pass[[1]], visim_ls_pass[-1])
all_data <- merge(visim_ls[[1]], visim_ls[-1])

## normalize for sketch data ------
all_data <- NormalizeData(all_data, scale.factor = floor(mean(all_data$nCount_Spatial.008um)))
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

sketch_cells_pass_qc <- all_data[['sketch']] %>% colnames() %>% .[. %in%  unlist(map(visim_ls, colnames))]
all_data[['sketch']]  <- all_data[['sketch']] %>% subset(cells = sketch_cells_pass_qc)

### dimension reduction
DefaultAssay(all_data) <- "sketch"
all_data <- FindVariableFeatures(all_data, nfeatures = 3000)
all_data <- ScaleData(
  all_data,
  vars.to.regress = c(
    'nCount_Spatial.008um'
  )
)


# all_data <- SCTransform(all_data, assay = "sketch")

all_data <- RunPCA(all_data, features = VariableFeatures(all_data))

pdf(file = file.path(dr_outdir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_data, ndims = 30)
dev.off()


# 加入分组信息
group_map<-read.table(file="sample_map_group.txt",header=T,sep="\t",comment.char = "")
sample.type<-vector("logical",length=ncol(all_data))
for (k in 1:ncol(all_data)){
  sample.type[k]<-as.vector(group_map[which(group_map[,1]==as.vector(all_data@meta.data$project[k]),arr.ind = TRUE),2])
}
all_data[["group"]] <- sample.type
all_data[["sample"]] <- all_data[["project"]]

set.seed(1234)
suset_cells <- 
  all_data[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(group) %>%
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
dims <- 15

all_data <- RunUMAP(all_data, dims = 1:dims, return.model = TRUE)


pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "group", label = TRUE, reduction = "umap") + NoLegend()
dev.off()


all_data <- all_data %>%
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "group",
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
DimPlot(all_data, group.by = "group", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

# validate QC metric
FeaturePlot(
  all_data,
  features  = c('nCount_Spatial.008um', 'nFeature_Spatial.008um', 
                'percent.mt', 'log10GenesPerUMI'),
  reduction = "harmonyumap"
)


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
    algorithm = 1, 
    resolution = c(0.05, 0.1, 0.2, 0.4, 0.6)
  )


# dimplot
pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "sketch_snn_res.0.05", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.1", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.05", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.1", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()

# DimPlot(all_data, group.by = "SCT_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


## project data ------
all_data <- 
  all_data %>% 
  ProjectData(
    assay = "Spatial.008um",
    full.reduction = "harmony.full",
    
    sketched.assay = "sketch",
    sketched.reduction = "harmony",
    
    umap.model = "harmonyumap", 
    dims = 1:dims, 
    
    refdata = list(
      sketch_snn_res.0.05_full = "sketch_snn_res.0.05",
      sketch_snn_res.0.1_full = "sketch_snn_res.0.1",
      sketch_snn_res.0.2_full = "sketch_snn_res.0.2",
      sketch_snn_res.0.4_full = "sketch_snn_res.0.4"
    )
  )


saveRDS(all_data, file = "data_for_analysis/all_data_projected.rds")

## dotplot ---------
major_marker_ls <- c(
  
  Endothelial = "PECAM1",
  Endothelial = "VWF",
  Endothelial = "AQP1",
  
  # Epithelial  = "EPCAM", 
  # Epithelial  = "CEACAM5", 
  # Epithelial  = "CEACAM6", 
  # Cholangiocyte  = "KRT7",
  # Cholangiocyte  = "KRT19",
  # Hepatocyte = 'ALB',
  # Hepatocyte = 'APOC3',
  
  'Follicular epithelial cell' =  'KRT17',
  'Follicular epithelial cell' =  'S100A2',
  
  'Sebaceous gland cell' = 'DCD',
  'Sebaceous gland cell' = 'FST',
  'Sebaceous gland cell' = 'MUCL1',
  
  'Keratinocytes undiff' = 'KRT5',
  'Keratinocytes undiff' = 'KRT14',
  'Keratinocytes undiff' = 'TP63',
  'Keratinocytes undiff' = 'ITGB1',
  'Keratinocytes undiff' = 'ITGA6',
  
  
  'Keratinocytes diff' = 'KRT1',
  'Keratinocytes diff' = 'KRT10',
  'Keratinocytes diff' = 'SBSN',
  'Keratinocytes diff' = 'KRTDAP',
  
  # "Fibroblasts" = "SFRP4", 
  # "Fibroblasts" = "CCDC80",
  # "Fibroblasts" = "COL1A1",
  
  "Fibroblasts" = "PDGFRA",
  "Fibroblasts" = "LUM",
  "Fibroblasts" = "DCN",
  "Fibroblasts" = "VIM",
  "Fibroblasts" = "COL1A2",
  
  # CAF = "ACTA2",
  # CAF = "COL1A2",
  # iCAF = "CXCL12",
  # iCAF = 'PDGFRA',
  # mCAF = 'RGS5',
  # mCAF = 'TAGLN',
  
  "Pericytes" = "PDGFRB",
  Pericytes = "ACTA2",
  Pericytes = 'RGS5',
  
  Erythrocytes = 'HBA1',
  Erythrocytes = 'HBA2',
  Erythrocytes = 'HBB',
  
  Melanocytes = 'PMEL',
  Melanocytes = 'MLANA',
  Melanocytes = 'TYRP1',
  Melanocytes = 'DCT',
  
  
  # "Smooth muscle cells" = "PDGFRB", 
  # "Smooth muscle cells" = "MYH11", 
  # "Smooth muscle cells" = "ACTA2",
  
  # "Breast glandular cells" = "SERPINA3", 
  # "Breast glandular cells" = "DSC2",
  
  
  
  Immune = 'PTPRC',
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  "Plasma" = "IGHG1",
  "Plasma" = "IGHA1",
  "Plasma" = "IGHG3",
  "Plasma" = "IGKC",
  "Plasma" = "TNFRSF17",
  "Plasma" = "JCHAIN",
  
  Plasmablast = "MZB1",
  Plasmablast = "XBP1",
  
  "T cells" = "CD3D",
  "T cells" = "CD3E",
  "T cells" = "CD3G",
  "T cells" = "LCK",
  "T cells" = "TRAC",
  "T cells" = "CD4",
  "T cells" = "CD8A",
  "T cells" = "GZMK",
  
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
  
  "Monocytes" = "CD14", 
  "Monocytes" = "FCGR3A",
  
  "Dendritic cell" = "CD83",
  'Dendritic cell' = "LILRA4",  
  
  'Dendritic cell' = "CLEC9A",
  'Dendritic cell' = "IDO1",
  
  'Dendritic cell' = "CD1C",
  'Dendritic cell' = "FCER1A",
  
  'Fat cell' = 'ADIPOQ',
  'Fat cell' = 'FABP4',
  'Fat cell' = 'G0S2'
)


skin_marker_ls <- c(
  Keratinocyte = 'KRT14',
  Keratinocyte = 'KRT1',
  Keratinocyte = 'DMKN',
  Keratinocyte = 'KRT10',
  Keratinocyte = 'KRT5',
  Keratinocyte = 'KRTDAP',
  
  Melanocyte = 'DCT',
  Melanocyte = 'TYRP1',
  Melanocyte = 'PMEL',
  Melanocyte = 'MLANA',
  Melanocyte = 'QPCT',
  Melanocyte = 'MIFT',
  
  'Eccrine gland' = "SCGB1D2",
  'Eccrine gland' = "PIP",
  'Eccrine gland' = "SCGB1B2P",
  'Eccrine gland' = "SCGB2A2",
  'Eccrine gland' = "MUCL1",
  'Eccrine gland' = "DCD",
  
  Endothelial = "CLDN5",
  Endothelial = "FABP4",
  Endothelial = "PECAM1",
  Endothelial = "CDH5",
  Endothelial = "TM4SF1",
  Endothelial = "CCL21",
  
  Fibroblast = "DCN",
  Fibroblast = "COL1A1",
  Fibroblast = "CFD",
  Fibroblast = "COL1A2",
  Fibroblast = "COL3A1",
  Fibroblast = "APOD",
  
  'Smooth muscle' = "TAGLN",
  'Smooth muscle' = "ACTA2",
  'Smooth muscle' = "MYL9",
  'Smooth muscle' = "RGS55",
  'Smooth muscle' = "TPM2",
  'Smooth muscle' = "CALD1",
  
  'Nerve cell' = 'MPZ',
  'Nerve cell' = 'PLP1',
  'Nerve cell' = 'S100B',
  'Nerve cell' = 'PMP22',
  'Nerve cell' = 'CRYAB',
  'Nerve cell' = 'MBP',
  
  'T cell' = 'IL32',
  'T cell' = 'CD52',
  'T cell' = 'CXCR4',
  'T cell' = 'CD3E',
  'T cell' = 'CD3D',
  'T cell' = 'TRAC',
  
  'Myeloid cell' = 'HLA-DRA',
  'Myeloid cell' = 'HLA-DPB1',
  'Myeloid cell' = 'CD74',
  'Myeloid cell' = 'LYZ',
  'Myeloid cell' = 'HLA-DPA1',
  'Myeloid cell' = 'HLA-DRB1',
  
  'Mast cell' = "TPSAB1",
  'Mast cell' = "CPA3",
  'Mast cell' = "CTSG",
  'Mast cell' = "HPGD",
  'Mast cell' = "HPGDS"
)

major_marker_ls <- c(
  
  Endothelial = "PECAM1",
  Endothelial = "CLDN5",
  
  Epithelial  = "EPCAM",
  
  "Fibroblasts" = "DCN",
  "Fibroblasts" = "COL1A2",
  
  
  "Pericytes" = "PDGFRB",
  Pericytes = "ACTA2",
  Pericytes = 'RGS5',
  
  
  Immune = 'PTPRC',
  Immune = "CD3D", 
  Immune = "TRAC",
  Immune = "LYZ", 
  Immune = "CD68"
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} main_marker.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "Spatial.008um_snn_res.0.2", features = major_marker_ls, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "sketch_snn_res.0.4", features = major_marker_ls, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "sketch_snn_res.0.6", features = major_marker_ls, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "sketch_snn_res.0.8", features = major_marker_ls, dot.scale = 10) + theme_90angle


# DotPlot(all_data, group.by = "SCT_snn_res.0.2", dot.scale = 10, features = major_marker_ls) + theme_90angle
# DotPlot(all_data, group.by = "SCT_snn_res.0.4", dot.scale = 10, features = major_marker_ls) + theme_90angle
# DotPlot(all_data, group.by = "SCT_snn_res.0.6", dot.scale = 10, features = major_marker_ls) + theme_90angle
# DotPlot(all_data, group.by = "SCT_snn_res.0.8", dot.scale = 10, features = major_marker_ls) + theme_90angle
dev.off()


pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} skin_marker.pdf")), width = 25, height = 6)

DotPlot(all_data, group.by = "sketch_snn_res.0.2", features = skin_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(all_data, group.by = "sketch_snn_res.0.4", features = skin_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(all_data, group.by = "sketch_snn_res.0.6", features = skin_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(all_data, group.by = "sketch_snn_res.0.8", features = skin_marker_ls, dot.scale = 10) + theme_90angle

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


## esophagus
eosinophils_markers <- c(
  "ITGA4",
  "IL5RA", 
  'CEACAM8',
  'FUT4',
  'CCR3',
  'CD125',
  "PTGDR2",
  "SIGLEC8",
  "ITGAM",
  
  'KRT5'
)

DotPlot(
  all_data, 
  features = eosinophils_markers, 
  assay    = 'sketch',
  group.by = "sketch_snn_res.0.2"
)


## annotation ----
Idents(all_data) <- "Spatial.008um_snn_res.0.2"


# cell_cluster <- c(
#   '0'   = "Follicular epithelial cells",
#   '1'   = "Macrophage/T",
#   '2'   = "Endothelial cells",
#   '3'   = "Fibroblasts",
#   '4'   = "Keratinocytes undiff",
#   '5'   = "Sebaceous gland cell",
#   '6'   = "Keratinocytes diff",
#   '7'   = "TMSB4X supporting cells",
#   '8'   = "Fat cell"
# )

cell_cluster <- c(
  '0'    = "Non_Immune",
  '1'    = "Non_Immune",
  '2'    = "Immune",
  '3'    = "Non_Immune",
  '4'    = "Non_Immune",
  '5'    = "Non_Immune",
  '6'    = "Non_Immune",
  '7'    = "Non_Immune",
  '8'    = "Non_Immune"
)


all_data <- RenameIdents(all_data, cell_cluster)
all_data$annotation <- Idents(all_data)


pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} annotation {reduction.use} .pdf")), height = 6, width = 8)
# DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

# pdf(str_glue(str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2.pdf")), width = 25, height = 6)
# DotPlot(all_data, group.by = "annotation", features = major_cell_marker)
# dev.off()

## sub-cluster -----------
DefaultAssay(all_data) <- "sketch"
Idents(all_data) <- "sketch_snn_res.0.2"

all_data <-
  FindSubCluster(
    all_data,
    cluster = '1',
    graph.name = 'sketch_snn',
    subcluster.name = "immune.sub.cluster",
    resolution = 0.4,
    algorithm = 1
  )

DimPlot(
  all_data,
  reduction = "harmonyumap",
  group.by = "immune.sub.cluster"
)


pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} immune_sub_marker.pdf")), width = 25, height = 6)
DotPlot(
  all_data,
  features = major_marker_ls %>% .[
    names(.) %in% c("Immune", "B cells", "Plasma",
                    "T cells", # "NK",
                    "Macrophage",
                    "Monocytes", "Dendritic cell")
  ] %>% 
    c("EPCAM", "PECAM1", "DCN"),
  group.by = "immune.sub.cluster"
) +
  theme_90angle

dev.off()

### annotation ----
Idents(all_data) <- "immune.sub.cluster"

immune_cell_sub_cluster <- c(
  '0'   = "Keratinocytes undiff", #"Follicular epithelial cells",
  '1_0'   = "T",
  '1_1'   = "Macrophage",
  '1_2'   = "Undefined",
  '1_3'   = "Undefined",
  '1_4'   = "Macrophage",
  '1_5'   = "Plasma",
  '1_6'   = "Macrophage",
  '2'   = "Endothelial cells",
  '3'   = "Fibroblasts",
  '4'   = "Keratinocytes undiff",
  '5'   = "Sebaceous gland cell",
  '6'   = "Keratinocytes diff",
  '7'   = "TMSB4X supporting cells",
  '8'   = "Fat cell"
)


all_data <- RenameIdents(all_data, immune_cell_sub_cluster)
all_data$annotation2 <- Idents(all_data)


pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} annotation2 {reduction.use} .pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation2", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(all_data, group.by = "annotation2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

# pdf(str_glue(str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2.pdf")), width = 25, height = 6)
# DotPlot(all_data, group.by = "annotation", features = major_cell_marker)
# dev.off()



## project data ------
# 2023.10.30 不做immune subcluster
all_data <- 
  all_data %>% 
  ProjectData(
    assay = "Spatial.008um",
    full.reduction = "pca.full",
    
    sketched.assay = "sketch",
    sketched.reduction = "harmony",
    
    umap.model = "harmonyumap", 
    dims = 1:dims, 
    
    refdata = list(
      annotation_full = "annotation",
      annotation2_full = "annotation2",
      sketch_snn_res.0.2_full = "sketch_snn_res.0.2",
      immune.sub.cluster_full = "immune.sub.cluster"
    )
  )

Idents(all_data) <- 'annotation2_full'
DefaultAssay(all_data) <- "Spatial.008um"

# annotation2
all_data$annotation2_full2 <- all_data$annotation2_full
all_data$annotation2_full2[all_data$annotation2_full.score < 0.9] <- NA

pdf(
  file.path(dr_outdir, str_glue("Dimplot dim{dims} annotation2 - {reduction.use} - projected.pdf")),
  height = 6,
  width = 8
)
DimPlot(all_data, group.by = "annotation2_full2", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
dev.off()

pdf(
  file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use} - projected.pdf")),
  height = 6,
  width = 8
)
DimPlot(all_data, group.by = "sketch_snn_res.0.2_full", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
dev.off()

# saveRDS(all_data, "all_data_projected.rds")
# all_data <- readRDS("all_data_projected.rds")

pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} annotation2 major_marker2 projected.pdf")), width = 25, height = 6)
DotPlot(all_data, assay = "sketch", group.by = "annotation2_full", features = major_marker_ls) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
dev.off()

## 2023.10.30 不做immune subcluster ------------ 
all_data <- 
  all_data %>% 
  ProjectData(
    assay = "Spatial.008um",
    full.reduction = "pca.full",
    
    sketched.assay = "sketch",
    sketched.reduction = "harmony",
    
    umap.model = "harmonyumap", 
    dims = 1:dims, 
    
    refdata = list(
      annotation_full = "annotation",
      sketch_snn_res.0.2_full = "sketch_snn_res.0.2"
    )
  )

Idents(all_data) <- 'annotation_full'
DefaultAssay(all_data) <- "Spatial.008um"

# annotation
all_data$annotation_full2 <- all_data$annotation_full
all_data$annotation_full2[all_data$annotation_full.score < 0.9] <- NA

pdf(
  file.path(dr_outdir, str_glue("Dimplot dim{dims} annotation - {reduction.use} - projected.pdf")),
  height = 6,
  width = 8
)
DimPlot(all_data, group.by = "annotation_full", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
DimPlot(all_data, group.by = "annotation_full2", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
dev.off()

pdf(
  file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use} - projected.pdf")),
  height = 6,
  width = 8
)
DimPlot(all_data, group.by = "sketch_snn_res.0.2_full", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
dev.off()

# saveRDS(all_data, "all_data_projected_v2.rds")




# 20240930 ---
# dotplot, remove cell name of cell markers

pdf(str_glue("dotplot dims{dims} {reduction.use} annotation2 unname-major_marker2 projected.pdf"), width = 25, height = 6)
DotPlot(all_data, 
        group.by = "annotation2_full",
        features = unname(major_marker_ls),
        dot.scale = 6) +
  theme_90angle
dev.off()


# export for rmd
png(file.path(img_output, str_glue("2 dimplot - cluster.png")), width = 600, height = 480)
DimPlot(
  all_data, 
  group.by = "sketch_snn_res.0.2_full", 
  label = TRUE, 
  reduction = 'full.harmonyumap'
) +
  NoLegend()
dev.off()

## QC pass assay -----------
all_data[['Spatial.qc.pass']] <- 
  all_data[['Spatial.008um']] %>%
  subset(cells = visim_ls_pass %>% map(colnames) %>% unlist(use.names = FALSE))

saveRDS(all_data, 'all_data.full.annotation.rds')


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
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
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
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
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
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
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
  # arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  arrange(desc(pct.1), .by_group = TRUE) %>%
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

# RCTD -----------------------
library(spacexr)
rctd_out_dir <- 'RCTD' # dir.create(rctd_out_dir)

# all_data <- readRDS("all_data.full.annotation.rds")

#1. sketch data  
# NOT RUN. ROI data, only 50k cells
# visium_roi_ls %>%
#   map(function(visium){
#     # sketch the cortical subset of the Visium HD dataset
#     DefaultAssay(visium) <- "Spatial.008um"
#     visium <- FindVariableFeatures(visium)
#     visium <- SketchData(
#       object = cortex,
#       ncells = 50000,
#       method = "LeverageScore",
#       sketched.assay = "sketch"
#     )
#     
#     DefaultAssay(visium) <- "sketch"
#     visium <- ScaleData(visium)
#     visium <- RunPCA(visium, assay = "sketch", reduction.name = "pca.sketch", verbose = T)
#     visium <- FindNeighbors(visium, reduction = "pca.sketch", dims = 1:50)
#     visium <- RunUMAP(visium, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50, verbose = T)
#     
#     visium
#   })

#2. single data
# all_data_sc <- readRDS('BP_AD_Pso-v2.rds')
all_data_sc[['RNA']] <- all_data_sc[['RNA']] %>% as("Assay5")
all_data_sc$annotation <- all_data_sc$celltype 

Idents(all_data_sc) <- "annotation"

count_sc_data <- readRDS('adata_counts.rds')
count_sc_data[['RNA']] <- count_sc_data[['RNA']] %>% as("Assay5")
counts_sc  <- count_sc_data %>% JoinLayers() %>% LayerData(layer = "counts")
cluster_sc <- as.factor(all_data_sc$annotation)
nUMI_sc    <- all_data_sc$total_counts

sc_cell_filter_idx <-  cluster_sc != "Unknown"

# create the RCTD reference object
reference <- Reference(
  counts_sc[, sc_cell_filter_idx], 
  droplevels(cluster_sc[sc_cell_filter_idx]), 
  nUMI_sc[sc_cell_filter_idx]
)

# 3. create the RCTD query object
rctd_query_ls <-
  visim_ls %>%
  map(function(visium){# browser()
    counts_hd <- visium %>% LayerData(layer = "counts")
    coords    <- visium %>% GetTissueCoordinates() %>% .[, 1:2]
    
    # create the RCTD query object
    query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))
    query
  })

#4. run RCTD
rctd_result_ls<- 
  rctd_query_ls %>%
  imap(function(query, name){
    # run RCTD
    RCTD <- create.RCTD(
      query, reference, max_cores = 10,
      UMI_min = 50
      # UMI_max
    )
    RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
    
    RCTD
  })

# saveRDS(rctd_result_ls, 'rctd_result_ls.rds')

#5. addMetaData
visium_rctd_ls <-
  visim_ls %>%
  imap(function(visium, name){ # browser()
    rctd <- rctd_result_ls[[name]]
    visium <- AddMetaData(visium, metadata = rctd@results$results_df)
    
    visium$first_type <- as.character(visium$first_type)
    visium$first_type[is.na(visium$first_type)] <- "Other" # "Unknown"
    visium
  })

#6. spatial plot 
# SpatialFeaturePlot(
#   visium_roi_ls[[1]], 
#   image.scale = 'hires', 
#   pt.size.factor = 5,
#   shape = 22,
#   alpha = c(0.5, 0.8),
#   features = "nFeature_Spatial.008um"
# )

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
    visium$first_type[visium$first_type == "Other"] <- NA
    
    SpatialDimPlot(
      visium,
      group.by       = 'first_type',
      image.scale    = 'lowres',
      pt.size.factor = 5 * image_scale_ls[[name]],
      shape          = 22,
      # alpha          = 0.5
    ) +
      ggtitle(name_fix) +
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5)))
  })
dev.off()

## RCTD remove non-tissue area -------
pdf(
  file.path(rctd_out_dir, 'spatial plot of RCTD - only tissue area.pdf'), 
  height = 15, 
  width = 15
)
visium_rctd_ls %>%
  imap(function(visium, name){ # browser()
    if(name == "HC") name_fix <- "NML" else name_fix <- name
    visium$first_type[visium$first_type == "Other"] <- NA
    
    
    image  <- visim_image_he_ls[[name]]
    image  <- RenameCells(image, new.names = paste0(name, "_", Cells(image)))
    
    coords <- visim_coords_ls[[name]]
    coords_fixed <- paste(name, coords, sep = "_")
    
    image_sub  <- subset(image,  cells = coords_fixed)
    visium_sub <- subset(visium, cells = coords_fixed)
    visium_sub@images[[1]] <- image_sub
    
    
    # all_celltype <- visium_sub$first_type %>% unique
    # cell_colors  <- structure(
    #   names = all_celltype,
    #   ggsci::pal_igv()(length(all_celltype))
    # )
    # cell_colors <- c(cell_colors, Other = NA)
    
    # remove others 
    visium_sub <- visium_sub[, !is.na(visium_sub$first_type)]
    
    SpatialDimPlot(
      visium_sub,
      group.by       = 'first_type',
      image.scale    = 'hires',
      pt.size.factor = 5 * image_scale_ls[[name]],
      shape          = 22,
      # alpha          = 0.5
    ) +
      ggtitle(name_fix) +
      ggsci::scale_fill_igv() +
      # scale_fill_manual(values = cell_colors) +
      guides(fill = guide_legend(override.aes = list(size = 5)))
  })
dev.off()

# spatia plot of some genes ---------------

gene_to_plot <- c(
  'CD3D',
  'CD69',
  'ITGAE',
  'CXCR6',
  
  'CD4',
  'CD8A',
  "IL13",
  "IL17A",
  "IL9",
  "IL22"
)

table(gene_to_plot %in% rownames(all_data))

pdf(
  file.path(sp_output, 'spatial plot of genes.pdf'), 
  height = 15, 
  width = 15
)



# all_data2 <- all_data
# change Image in all_data to hires.
for(img_idx in seq_along(visim_image_he_ls)){
  
  image_new <- visim_image_he_ls[[img_idx]]
  name <- names(visim_image_he_ls)[[img_idx]]
  
  image_name <- paste0(name, '.008um')
  all_image_name <- Images(all_data)
  
  if(!image_name %in% all_image_name && name == "NML"){
    name <- "HC"
    image_name <- paste0(name, '.008um')
  }
  
  image <- all_data@images[[image_name]]
  
  image_new  <- RenameCells(image_new, new.names = paste0(name, "_", Cells(image_new)))
  Key(image_new) <- Key(image)
  
  all_data@images[[image_name]] <- image_new
}

# coords subset
all_coords <-
  visim_coords_ls %>% 
  imap(function(coords, name){
    if(name == "NML") name <- "HC"
    paste(name, coords, sep = "_")
  }) %>%
  unlist(use.names = FALSE)

all_data <- all_data %>% subset(cells = all_coords)
names(all_data@images) <-
  Images(all_data) %>% 
  str_remove("\\.[0-9]{3}um$") %>% 
  str_replace("^HC$", "NML")



p_all_spatial_ls <-
  SpatialFeaturePlot(
    all_data,
    images         = NULL,
    features       = gene_to_plot,
    image.scale    = 'hires',
    pt.size.factor = 5 * 2, # * image_scale_ls[[name]],
    shape          = 22,
    keep.scale     = 'all',
    combine        = FALSE,
    alpha          = c(0, 1) # 0.5
  )


idx_row2col <-
  matrix(seq_along(p_all_spatial_ls), byrow = TRUE, ncol = 4) %>% 
  c()


p_spatial_combined_all <- 
  patchwork::wrap_plots(
    p_all_spatial_ls[idx_row2col],
    nrow = 4,
    # guides = "collect"
  ) &
  scale_fill_gradient(low = "#FFFFFF11", high = "red", trans = 'log1p') 


pdf(
  file.path(sp_output, 'spatial plot of genes_all_combined.pdf'), 
  height = 15, 
  width = 25
)
print(p_spatial_combined_all)
dev.off()


## spatial plot of some genes combined ------
pdf(
  file.path(sp_output, 'spatial plot of genes.pdf'), 
  height = 15, 
  width = 15
)
all_data %>%
  Images() %>%
  set_names(., str_remove(., "\\..*$")) %>%
  imap(function(image, name){
    if(name == "HC") name_fix <- "NML" else name_fix <- name
    
    image_scale_ls <-
      all_data %>%
      Images() %>%
      set_names(., str_remove(., "\\..*$")) %>%
      map(~ all_data@images[[.x]]) %>%
      map(GetTissueCoordinates) %>%
      map('x') %>%
      map(range) %>%
      map(diff) %>%
      map(log) %>%
      {map(., function(x, max_val)x/max_val, max_val = max(unlist(.)))}
    
    SpatialFeaturePlot(
      all_data,
      images         = image,
      features       = gene_to_plot,
      image.scale    = 'lowres',
      pt.size.factor = 5 * image_scale_ls[[name]],
      shape          = 22,
      alpha          = c(0, 1) # 0.5
    ) +
      patchwork::plot_annotation(title = name) &
      scale_fill_gradient(low = "#FFFFFF11", high = "red", trans = 'log1p') 
    # guides(fill = guide_legend(override.aes = list(size = 5)))
  })
dev.off()




# spatial plot --------------
# Idents(all_data) <- 'Spatial.008um'
# Idents(all_data) <- 'Spatial.qc.pass'


p_sp_plot_anno <-
  Images(all_data) %>%
  setNames(., c('DZZ_Saline_0um','DZZ_Saline_100um','SYZ_NH4_0um','SYZ_NH4_100um')) %>%
  imap(function(img, img_name){
    p <- SpatialDimPlot(
      all_data,
      images = img, 
      image.scale = "hires", 
      group.by = "annotation_full",
      pt.size.factor	= 10
    ) +
      #ggsci::scale_fill_aaas() +
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      ggtitle(img_name)
    p
  })

p_sp_plot_anno <-
  Images(all_data2) %>%
  setNames(., str_remove(., "\\.[0-9]{3}um$")) %>%
  imap(function(img, img_name){
    p <- SpatialDimPlot(
      all_data2,
      images = img, 
      group.by = "annotation2_full",
      pt.size.factor	= 4
    ) +
      # ggsci::scale_fill_aaas() +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#7A65A5FF","#CE3D32FF","#749B58FF","#F0E685FF","#466983FF","#BA6338FF",
                                 "#5DB1DDFF","#802268FF","#924822FF","#D595A7FF","#6BD76BFF","#5050FFFF"))+
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      ggtitle(img_name)
    p
  })


p_sp_plot_anno <-
  Images(all_data) %>%
  setNames(., str_remove(., "\\.[0-9]{3}um$")) %>%
  imap(function(img, img_name){
    p <- SpatialDimPlot(
      all_data,
      images = img, 
      image.scale = "hires",
      group.by = "first_type",
      pt.size.factor	= 3,
      shape = 22
    ) +
      # ggsci::scale_fill_aaas() +
      # scale_shape_manual(values =22) +
      ggsci::scale_fill_igv() +
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

# spatial plot2 ----------
p_sp_plot_anno2 <-
  Images(all_data) %>%
  setNames(., str_remove(., "\\.[0-9]{3}um$")) %>%
  imap(function(img, img_name){
    if(img_name == "HC") img_name <- "NML"
    dat <- 
      all_data[['annotation2_full2']] %>%
      rownames_to_column('cell_id') %>%
      inner_join(
        all_data@images[[img]] %>% GetTissueCoordinates() %>% rownames_to_column('cell_id'),
        by = 'cell_id'
      )
    dat$annotation2_full2[dat$annotation2_full2 == 'Undefined'] <- NA
    p <-
      dat %>%
      ggplot(aes(x = y, y = -x, color = annotation2_full2)) +
      geom_point(size = 0.1, shape = 15) +
      theme_void() +
      ggsci::scale_color_igv() +
      guides(color = guide_legend(override.aes = list(size = 3))) +
      ggtitle(img_name)
    p
  })

# export
pdf(
  file.path(sp_output, "spatial plot annotation - point.pdf"),
  height = 8, 
  width  = 8
)
p_sp_plot_anno2 %>% walk(print)
dev.off()


# enrichemt analysis -----
library(clusterProfiler)

all_de_genes_id_for_enrichmenet <-
  all_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::bitr(
      x,
      OrgDb = "org.Hs.eg.db",
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
dr_outdir2 <- str_glue("{dr_outdir} - 2") # dir.create(dr_outdir2)
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


