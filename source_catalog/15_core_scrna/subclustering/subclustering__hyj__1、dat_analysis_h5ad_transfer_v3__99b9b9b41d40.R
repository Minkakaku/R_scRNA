library(tidyverse)
library(Seurat)
library(patchwork)


# setting variable   --------------------
data_dir  <- "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2512177577-01/rawdata"
list.files(data_dir)
anno_file <- "data/Annotation.csv"

dir.create("R_OUTPUT")
dir.create("R_OUTPUT/data_for_analysis")



# library --------------
options(java.parameters = '-Xmx10g') # options('java.parameters')
library(rJava)
#library(RBioFormats)
# source("xenium.R")

# multi workers  -----
# library(future)
# plan("multisession", workers = 10, gc = TRUE)
# 
# options(future.globals.maxSize = 10 * 1024 ^ 3)
# 
# nbrOfFreeWorkers()


# define function --------------
parseProjName_xenium <- 
  function(x,
           pattern = "^output.*",
           pattern2 = ".*?__.*?__.*?__.*?__.*?",
           split_pattern = "__", 
           extract_postion =  c(3),
           check.unique = TRUE) { 
    
    xenium_out_idx <- check_xenium_pattern(
      x, 
      pattern = pattern, 
      pattern2 = pattern2,
      use_basename = TRUE
    )
    
    if(sum(xenium_out_idx) < 1){
      stop(str_glue("No output diretory match the pattern: {pattern} and pattern: {pattern2}"))
    }
    
    x <- x[xenium_out_idx]
    
    
    stopifnot(length(x) > 0) 
    
    df <- tibble::tibble(.path = x)
    df <- df %>% 
      mutate(
        .dir = dirname(.path),
        .base = basename(.path),
        .proj = str_remove_all(.base, pattern)
      ) %>% 
      mutate(
        .proj = str_remove_all(.proj, "[._-]?$"),
        # .proj = .proj %>% str_split("__") %>% map_chr(3)
        .proj = .proj %>% extract_pattern(split_pattern = split_pattern, extract_postion = extract_postion)
      ) %>%
      dplyr::filter(
        str_detect(.base, pattern = pattern)
      )
    
    # if proj duplicated, send a message
    is_duplicated <- duplicated(df$.proj)
    if(sum(is_duplicated) > 0){
      message("Not all project names be uniquely identified...")
    }
    
    if(check.unique && sum(is_duplicated) > 0){
      
      message("Duplicated project name(s) set as unique label.")
      df$.proj <- make.names(df$.proj, unique = TRUE)
    }
    
    # return
    df <- as.data.frame(df)
    res_tryCatch <- 
      tryCatch({
        rownames(df) <- df$.proj
        "work"
      },
      error = function(x) return("error")
      )
    if(res_tryCatch == "error"){
      rownames(df) <- make.names(df$.path, unique = TRUE)
    }
    return(df)
  }


check_xenium_pattern <-
  function(x, 
           pattern = "^output", 
           pattern2 = ".*?__.*?__.*?__.*?__.*?",
           use_basename = TRUE){
    
    stopifnot(length(x) > 0)
    
    if(use_basename){
      x <- basename(x)
    }
    
    stopifnot(length(pattern)  < 2)
    stopifnot(length(pattern2) < 2)
    
    if(length(pattern) == 1){
      idx <- x %>% str_detect(pattern = pattern)
    }else{
      idx <- TRUE
    }
    
    if(length(pattern2) == 1){
      idx2 <- x %>% str_detect(pattern = pattern2)
    }else{
      idx2 <- TRUE
    }
    
    # return
    return(
      idx & idx2
    )
    
  }


extract_pattern <- function(
    x, 
    split_pattern = "__", 
    extract_postion =  c(3)){
  
  stopifnot(is.character(x))
  stopifnot(length(extract_postion) > 0)
  stopifnot(is.numeric(extract_postion))
  
  
  if(any(!str_detect(x, split_pattern))) stop(str_glue("Not found pattern {split_pattern} in one or more of {x}"))
  
  xx <- x %>% str_split(split_pattern)
  xx_max_len <- xx %>% map_dbl(length) %>% max
  
  
  if(max(extract_postion) > xx_max_len){
    stop("extract_postion must be wrong, causing maximum value of extract_postion bigger than split charater.")
  }
  
  res <- xx %>% map(~ .[extract_postion]) %>% map_chr(str_c, collapse = split_pattern)
  
  res
  
}




# reading data -------------
## data path and project name  -----
parse_data_path <-
  dir(data_dir, full.names = TRUE) %>% 
  parseProjName_xenium(pattern = "output[-]?") %>%
  mutate(.counts_path = file.path(.path, "cell_feature_matrix")) %>%
  identity()

## image w/h ratio ------
# used to fix width and heigth ratio in In Situ Spatial plot.
# image_wh_r_ls <-
#   with(parse_data_path, split(.path, .proj)) %>%
#   map(function(x){
#     imag_wh <-
#       RBioFormats::read.metadata(
#         file.path(x, "morphology.ome.tif")
#       ) %>%
#       .[[1]] %>%
#       .$coreMetadata %>%
#       .[c("sizeX", "sizeY")]
#     
#     image_wh_r <- imag_wh$sizeX/imag_wh$sizeY
#     image_wh_r
#   })


## read into annotation -------
if(file.exists(anno_file)){
  anno_data <- readr::read_csv(anno_file)
  User_columns <- colnames(anno_data)
}else{
  anno_data <- NULL
  User_columns <- "project"
}


## seurat object ------
all_data_ls <-
  parse_data_path %>% 
  dplyr::select(.counts_path, .proj) %>%
  pmap(function(.counts_path, .proj){
    # merge anno_data
    if(!is.null(anno_data)){
      row_idx <- which(anno_data$project == .proj)
      col_idx <- anno_data %>% colnames() 
      
      # 没有找到当前样本的注释数据
      #<<<<<< 所以Annotation.csv也有筛选文件的作用 >>>>>#
      if(length(row_idx) < 1){
        message(str_glue("Not found sample {.proj} in column: project of Annotation table.") )
        return(NULL)
      }
      # 没有任何数据，报错
      if(length(col_idx) < 1) stop("Not found any data in Annotation table.")
      
      meta_dat <- anno_data %>% .[row_idx, col_idx] %>% as.list()
    }else{
      meta_dat <- list(project = .proj)
    }
    
    # 构建Seurat
    sc <-
      Read10X(.counts_path) %>% 
      .[["Gene Expression"]] %>%
      CreateSeuratObject(
        project = .proj,
        min.cells    = 5,
        min.features = 5
      ) %>%
      AddMetaData(meta_dat) %>%
      # RenameCells(add.cell.id = .id)
      RenameCells(add.cell.id = .proj)
    sc
  }) 
all_data_ls <- all_data_ls %>% purrr::discard(is.null)
names(all_data_ls) <- map_chr(all_data_ls, "project.name")

# project number, used as condition judges
proj_name  <- names(all_data_ls)


# merged object
all_data <- merge(all_data_ls[[1]], all_data_ls[-1])


saveRDS(all_data, 'all_data.rds')


## read coords -------------------
# cell coordinate and cell boundary
cells_position_ls <- 
  # with(parse_data_path, split(.path, .proj)) %>%
  parse_data_path %>%
  group_by(.proj) %>%
  nest %>%
  deframe() %>%
  imap(function(x, proj){
    xx <-
      readr::read_csv(
        file.path(x$.path, "cells.csv.gz"),
        col_types = list(cell_id = "c")
      ) %>%
      mutate(cell_id = str_c(proj, cell_id, sep = "_"))
    xx
  })
cells_position_ls <- cells_position_ls[proj_name]

cells_boundary_ls <- 
  # with(parse_data_path, split(.path, .proj)) %>%
  parse_data_path %>%
  group_by(.proj) %>%
  nest %>%
  deframe() %>%
  imap(function(x, proj){
    xx <-
      readr::read_csv(
        file.path(x$.path, "cell_boundaries.csv.gz"),
        col_types = list(cell_id = "c")
      ) %>%
      mutate(cell_id = str_c(proj, cell_id, sep = "_"))
    xx
  })
cells_boundary_ls <- cells_boundary_ls[proj_name]

# export data
saveRDS(cells_boundary_ls, file.path("R_OUTPUT", "data_for_analysis", "cells_boundary_ls_raw.rds"))
saveRDS(cells_position_ls, file.path("R_OUTPUT", "data_for_analysis", "cells_position_ls_raw.rds"))




# FORMAT CONVERSION ---------
# config python
library(reticulate)
use_condaenv(condaenv = "/home/szt/miniforge3/envs/rpy/bin/python")
py_config()
py_available()
py_module_available("pandas")


sc <- import("scanpy")
adata <- sc$read_h5ad("/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/AE2512177577-01/output/data/adata_merged.h5ad")

adata_meta <- adata$obs
adata_umap <- adata$obsm['X_umap']


# sample match
if("sample_raw" %in% colnames(adata_meta)){
  adata_meta_new <-
    adata_meta %>% 
    # left_join(sample_match_df2, by = c(sample = "to"))
    mutate(
      from = sample_raw
    )
}else{
  # find original sample_id
  sample_match_df <- 
    data.frame(
      from = c("FWSP20250913", "FWSP20250913",	"FWSP20250913",	"FWSP20250913",	"FWSP20250913",	"FWSP20250913",	"FWSP20250913",	"FWSP20250913",	"20_19345_C3",	"20_37222_B2",	"22_010650_A1",	"23_010125_A1_23_099429_C1",	"23_010125_A1_23_099429_C1",	"23_020064_A2_20_26379_C1",	"23_020064_A2_20_26379_C1",	"19_19186_B3",	"19_49709_B2",	"19_55625_E2",	"21_58073_C4",	"21_76852_C1",	"22_013510_C2",	"22_014073_B1",	"21_36100_H4",	"22_079225_C1",	"23_011238_C1",	"23_03024_B3",	"23_036910_B1",	"23_037718_C3",	"23_093096_B2",	"23_29818_D2",	"18_11211_C1",	"19_60381_B3",	"21_04911_C3",	"22_033641_F1",	"22_094337_C2",	"23_020665_C2",	"23_036986_C2",	"201501994_2019120979",	"201501994_2019120979",	"201639643_Y202272139",	"201639643_Y202272139",	"201750656",	"2018103363",	"2021005471",	"2022136838",	"2019115971",	"2020008909",	"2021047701",	"2021050076",	"2022036398",	"2023082036",	"Y202116587",	"Z202209248",	"201710508",	"2019110264",	"2020036918",	"2020107233",	"2020108622",	"2020111518",	"2021134652",	"2022016790",	"201539405",	"201747058",	"2019089885",	"2022101285",	"2022124807-2022095049",	"2022124807-2022095049",	"2022137338",	"2023121844"
) ,
      to <- c("202157344",	"2024012910",	"20240365013",	"2024036703",	"20-72372",	"240724373",	"2436291-3",	"Z2021409303",	"20_19345_C3",	"20_37222_B2",	"22_010650_A1",	"23_010125_A1",	"23_099429_C1",	"23_020064_A2",	"20_26379_C1",	"19_19186_B3",	"19_49709_B2",	"19_55625_E2",	"21_58073_C4",	"21_76852_C1",	"22_013510_C2",	"22_014073_B1",	"21_36100_H4",	"22_079225_C1",	"23_011238_C1",	"23_03024_B3",	"23_036910_B1",	"23_037718_C3",	"23_093096_B2",	"23_29818_D2",	"18_11211_C1",	"19_60381_B3",	"21_04911_C3",	"22_033641_F1",	"22_094337_C2",	"23_020665_C2",	"23_036986_C2",	"201501994",	"2019120979",	"201639643",	"Y202272139",	"201750656",	"2018103363",	"2021005471",	"2022136838",	"2019115971",	"2020008909",	"2021047701",	"2021050076",	"2022036398",	"2023082036",	"Y202116587",	"Z202209248",	"201710508",	"2019110264",	"2020036918",	"2020107233",	"2020108622",	"2020111518",	"2021134652",	"2022016790",	"201539405",	"201747058",	"2019089885",	"2022101285",	"2022124807",	"2022095049",	"2022137338",	"2023121844"
) 
    ) 
  
  if(nrow(sample_match_df) > 0){
    sample_match_df2 <-
      sample_match_df %>%
      right_join(parse_data_path['.proj'], by = c(from = '.proj')) %>%
      mutate(to = ifelse(is.na(to), from, to))
  }else{
    sample_match_df2 <- 
      parse_data_path %>%
      dplyr::select(
        from = .proj, 
        to = .proj,
      )
  }
  
  
  adata_meta_new <-
    adata_meta %>% 
    left_join(sample_match_df2, by = c(sample = "to"))
}


T1 <- data.frame(adata_meta_new$from %>% unique())

stopifnot(all(!is.na(adata_meta_new$from)))


if("cell_id_full" %in% colnames(adata_meta_new)){
  adata_meta_new$cell_id2 <- adata_meta_new$cell_id_full
}else{
  adata_meta_new$cell_id2 <- str_c(adata_meta_new$from, adata_meta_new$cell_id, sep = "_")
}

# some cells duplicated 
dup_cell_ids <- adata_meta_new$cell_id2 %>% .[duplicated(.)]
dup_cell_idx <- adata_meta_new$cell_id2 %in% dup_cell_ids

if(sum(dup_cell_idx) > 0){
  warning(str_glue("{sum(dup_cell_idx)} cell_ids are duplicated and have been removed"))
  adata_meta_new <- adata_meta_new[!dup_cell_idx, ]
  adata_umap <- adata_umap[!dup_cell_idx, ]
}


# add METADATA
all_data2 <- all_data %>% subset(cells = adata_meta_new$cell_id2)
all_data2@meta.data <-
  all_data2@meta.data %>%
  rownames_to_column("cell_id2") %>%
  left_join(adata_meta_new, by = "cell_id2") %>%
  column_to_rownames("cell_id2")


# add UMAP
adata_umap2 <- adata_umap
colnames(adata_umap2) <- str_c("UMAP", 1:ncol(adata_umap2), sep = "_")
rownames(adata_umap2) <- adata_meta_new$cell_id2
adata_umap2 <- adata_umap2[colnames(all_data2), ]

all_data2[['umap']] <- 
  CreateDimReducObject(
    adata_umap2,
    key    = "UMAP_",
    assay  = DefaultAssay(all_data2),
    global = TRUE
  )

# add useful columns
cluster_res_cols <- all_data2[[]] %>% colnames() %>% str_subset("^res_")
for(col_name in cluster_res_cols){
  col_name_new <- col_name %>% str_replace("res_", str_glue("{DefaultAssay(all_data2)}_snn_res.")) 
  all_data2[[col_name_new]] <- all_data2[[col_name]][[1]]
}
all_data2$annotation <- all_data2$celltype


# plot
all_data2 %>% DimPlot(reduction = "umap", group.by = "celltype", label = TRUE) + NoLegend()
all_data2 %>% DimPlot(reduction = "umap", group.by = "annotation", label = TRUE) + NoLegend()
all_data2 %>% DimPlot(reduction = "umap", group.by = "res_0.8", label = TRUE) + NoLegend()
all_data2 %>% DimPlot(reduction = "umap", group.by = "RNA_snn_res.0.8", label = TRUE) + NoLegend()




#sketch data --------
all_data2 <-
  SketchData(
    all_data2,
    ncells         = 10000L,
    sketched.assay = "sketch",
    method         = "LeverageScore",
    var.name       = "leverage.score",
    seed = 123L
  )


## dimension reduction - pca ------------
dr_dir <- "dimension reduction" # dir.create(dr_dir)

DefaultAssay(all_data2) <- "sketch"
all_data2 <- FindVariableFeatures(all_data2, nfeatures = 3000)
all_data2 <- ScaleData(all_data2, vars.to.regress = 'nCount_RNA')

all_data2 <- RunPCA(all_data2)

pdf(file.path(dr_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_data2, ndims = 30)
dev.off()


set.seed(1234)
suset_cells <- 
  all_data2[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(all_data2) ) %>%
  dplyr::pull("cell_id")

pdf(file.path(dr_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(all_data2[, suset_cells], dims = 1:5)
dev.off()

pdf(file.path(dr_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(all_data2[, suset_cells], dims = 6:10)
dev.off()

pdf(file.path(dr_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(all_data2[, suset_cells], dims = 11:15)
dev.off()

pdf(file.path(dr_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(all_data2[, suset_cells], dims = 16:20)
dev.off()


## dimension reduction --------
dims <- 20


library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


all_data2 <- RunUMAP(all_data2, dims = 1:dims, return.model = TRUE)


pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
DimPlot(all_data2, group.by = "project", label = TRUE, reduction = "umap") + NoLegend()
dev.off()


all_data2 <- all_data2 %>%
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    reduction.save = "harmony"
  )

all_data2 <- RunUMAP(
  all_data2,
  dims = 1:dims,
  reduction = "harmony",
  reduction.name = "harmonyumap",
  reduction.key  = "harmonyUMAP_", 
  return.model = TRUE
)

pdf(dr_dir, file.path(str_glue("Dimplot dim{dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data2, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


reduction.use <- 'harmony' # 'harmony' 'pca'

## cluster ------
all_data2 <-
  FindNeighbors(
    all_data2,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )

all_data2 <-
  FindClusters(
    all_data2,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data2, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data2, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data2, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data2, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data2, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data2, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data2, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data2, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


## # cell marker ----
major_marker_ls <- c(
  
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
  
  # "Fibroblasts" = 'COL1A2',
  # "Fibroblasts" = 'COL3A1',
  # "Fibroblasts" = 'LUM',
  # "Fibroblasts" = 'FN1',
  # "Fibroblasts" = 'COL6A1',
  # "Fibroblasts" = 'FBLN1',
  # 
  # "Fibroblasts" = 'VIM',
  # "Fibroblasts" = 'THY',
  # "Fibroblasts" = 'S100A4',
  # "Fibroblasts" = 'PDPN',
  # "Fibroblasts" = 'FAP',
  
  
  # 
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


fibro_markers <- c(
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "COL1A1",
  "Fibroblasts" = "DCN",
  "Fibroblasts" = "ACTA2",
  
  "Fibroblasts" = 'COL1A2',
  "Fibroblasts" = 'COL3A1',
  "Fibroblasts" = 'LUM',
  "Fibroblasts" = 'FN1',
  "Fibroblasts" = 'COL6A1',
  "Fibroblasts" = 'FBLN1',
  
  "Fibroblasts" = 'VIM',
  "Fibroblasts" = 'THY',
  "Fibroblasts" = 'S100A4',
  "Fibroblasts" = 'PDPN',
  "Fibroblasts" = 'FAP'
)

fibro_markers %in% rownames(all_data)

# dotplot
pdf(file.path(dr_dir, str_glue("dotplot dims{dims} {reduction.use} major_marker2.pdf")), width = 25, height = 6)
DotPlot(all_data2, group.by = "sketch_snn_res.0.2", features = major_marker_ls)
DotPlot(all_data2, group.by = "sketch_snn_res.0.4", features = major_marker_ls)
DotPlot(all_data2, group.by = "sketch_snn_res.0.6", features = major_marker_ls)
DotPlot(all_data2, group.by = "sketch_snn_res.0.8", features = major_marker_ls)
dev.off()

## annotation ----
Idents(all_data2) <- "sketch_snn_res.0.2"



cell_cluster <- c(
  '0'  = "Epithelial",
  '1'  = "Fibroblast",
  '2'  = "Epithelial",
  '3'  = "Macrophage",
  '4'  = "Plasma",
  '5'  = "T",
  '6'  = "Endothelial",
  '7'  = "Myefibroblast",
  '8'  = "Pericyte",
  '9'  = "B",
  '10' = "Gial",
  '11' = "Mast cell",
  '12' = "NK",
  '13' = "Others",
  '14' = "Others"
)


all_data2 <- RenameIdents(all_data2, cell_cluster)
all_data2$annotation <- Idents(all_data2)


pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} annotation {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data2, group.by = "annotation", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(all_data2, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

pdf(file.path(dr_dir, str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2.pdf")), width = 25, height = 6)
DotPlot(all_data2, group.by = "annotation", features = major_marker_ls)
dev.off()


## find markers -------
Idents(all_data2) <- "sketch_snn_res.0.2"
all_markers <-
  FindAllMarkers(
    JoinLayers(all_data2),
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


# project data ------
DefaultAssay(all_data2) <- "sketch"
all_data2 <- 
  all_data2 %>% 
  ProjectData(
    assay = "RNA",
    full.reduction = "harmony.full",
    
    sketched.assay = "sketch",
    sketched.reduction = "harmony",# "pca",
    
    umap.model = "harmonyumap", 
    dims = 1:dims, 
    
    refdata = list(
      annotation_full = "annotation",
      RNA_snn_res.0.2 = 'sketch_snn_res.0.2',
      RNA_snn_res.0.4 = 'sketch_snn_res.0.4',
      RNA_snn_res.0.6 = 'sketch_snn_res.0.6',
      RNA_snn_res.0.8 = 'sketch_snn_res.0.8'
    )
  )

Idents(all_data2) <- 'annotation_full'
DefaultAssay(all_data2) <- "RNA"

pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} annotation - {reduction.use} - projected.pdf")), height = 6, width = 8)
DimPlot(all_data2, group.by = "annotation_full", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
dev.off()

# saveRDS(all_data2, "all_data_projected.rds")
# all_data2 <- readRDS("all_data_projected.rds")

pdf(file.path(dr_dir, str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2 projected.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "annotation_full", features = major_cell_marker)
dev.off()

# # 发现T与plasma标错顺序了
# all_data$annotation_full2 <- all_data$annotation_full
# all_data$annotation_full2[all_data$annotation_full == "T"]      <- "Plasma"
# all_data$annotation_full2[all_data$annotation_full == "Plasma"] <- "T"
# all_data$annotation_full  <- all_data$annotation_full2

saveRDS(all_data2, 'all_data.rds')

## cell props barplot  --------

all_data_meta <-
  all_data2[[]] %>%
  rownames_to_column("cell_id")


p_col_proj <-
  all_data_meta %>%
  group_by(project, annotation_full) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = project, y = freq, fill = annotation_full)) +
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


p_box_group  <-
  all_data_meta %>%
  group_by(group, project, annotation_full) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = annotation_full, y = freq, fill = group)) +
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
  filename = file.path(dr_dir, "cell_freq_project.pdf"),
  height = 6,
  width  = 10
)

ggsave(
  p_box_group,
  filename = file.path(dr_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 10
)

## spatial plot --------------

# all_data_meta <-
#   all_data[[]] %>%
#   rownames_to_column("cell_id")

sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation_full, group = cell_id))  +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      ggsci::scale_fill_igv() +
      theme_void()
    p
  })


# export image
sp_dir <- "spatial plot" # dir.create(sp_dir)

sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(sp_dir, str_glue("spatial plot {grp}.pdf")),
      height = 12,
      width  = 12
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })
































# split ROI ------------
# 20250421 add

#1 split layers
all_data2 <- all_data2 %>% JoinLayers() %>% split(all_data2$sample)


#2 cells_position_ls
#3 cells_boundary_ls
roi_match_ls <-
  all_data2[[c('from', 'sample')]] %>%
  distinct() %>%
  group_by(from) %>%
  nest() %>%
  deframe() %>%
  map(unlist, use.names = FALSE)

roi_cell_id_ls <-
  all_data2[['sample']] %>% 
  rownames_to_column("cell_id") %>%
  group_by(sample) %>%
  nest() %>%
  deframe() %>%
  map(unlist, use.names = FALSE)



cells_position_ls2  <- list()
cells_boundary_ls2  <- list()

for(idx in seq_along(roi_match_ls)){
  from_name <- names(roi_match_ls)[idx]
  to_name   <- roi_match_ls[[idx]]
  
  print(str_glue("###### {from_name} ->"))
  for(sub in to_name){
    print(str_glue("  -> {sub} ######"))
    cell_id_sub <- roi_cell_id_ls[[sub]]
    
    cells_position_new <- cells_position_ls[[from_name]] %>% dplyr:::filter(cell_id %in% cell_id_sub)
    cells_boundary_new <- cells_boundary_ls[[from_name]] %>% dplyr:::filter(cell_id %in% cell_id_sub)
    
    cells_position_ls2[[sub]] <- cells_position_new
    cells_boundary_ls2[[sub]] <- cells_boundary_new
  }
}

# export data
saveRDS(cells_boundary_ls2, file.path("R_OUTPUT", "data_for_analysis", "cells_boundary_ls.rds"))
saveRDS(cells_position_ls2, file.path("R_OUTPUT", "data_for_analysis", "cells_position_ls.rds"))


p1 <- cells_boundary_ls2[['S9']] %>%
  slice_sample(prop = 0.5) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()

p2 <- cells_boundary_ls2[['S3']] %>%
  slice_sample(prop = 0.5) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()


p <- cells_boundary_ls[['S3']] %>%
  slice_sample(prop = 0.25) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()

library(patchwork)
(p2 | p1)  / p




# normalize ----------------------------------------------------------------------
all_data2 <- NormalizeData(all_data2, scale.factor = median(all_data2$nCount_RNA))
all_data2 <- FindVariableFeatures(all_data2)

# PCA ----------------------------------------------------------------------------
## PCA --------
#DefaultAssay(all_data2) <- "sketch"
#all_data2 <- FindVariableFeatures(all_data2)
VariableFeaturePlot(all_data2)

all_data2 <- ScaleData(all_data2, vars.to.regress = 'nCount_RNA')
all_data2 <- RunPCA(all_data2)

dr_dir <- '2. dimension reduction and cluster/'
dir.create(dr_dir)

### pca inspector ------
pdf(file.path(dr_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_data2, ndims = 50)
dev.off()

all_data2@meta.data %>% head()

all_data2$cell_id <- NULL

set.seed(1234)
suset_cells <-
  all_data2[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(sample) %>%
  slice_sample(prop = 50000/ncol(all_data2) ) %>%
  dplyr::pull("cell_id")

pdf(file.path(dr_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(all_data2[, suset_cells], dims = 1:5)
dev.off()

pdf(file.path(dr_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(all_data2[, suset_cells], dims = 6:10)
dev.off()

pdf(file.path(dr_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(all_data2[, suset_cells], dims = 11:15)
dev.off()

pdf(file.path(dr_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(all_data2[, suset_cells], dims = 16:20)
dev.off()




## dimension reduction -----------------------------------------------------------------
dims <- 30

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


#all_data2 <- RunUMAP(all_data2, dims = 1:dims, return.model = TRUE, n.neighbors = 30)

all_data2 <- all_data2 %>%
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
  )

all_data2 <- RunUMAP(
  all_data2,
  dims = 1:dims,
  reduction = "harmony",
  reduction.name = "harmonyumap",
  reduction.key  = "harmonyUMAP_", 
  return.model = TRUE
)

pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data2, group.by = "sample", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data2, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} group - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data2, group.by = "group", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data2, group.by = "group", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

# export 
saveRDS(all_data2, file = file.path("R_OUTPUT", "all_data.rds"))


# 生成列联表：统计group1下各sample的细胞数量
table(all_data2@meta.data$annotation, all_data2@meta.data$"celltype") %>% 
  as.data.frame() %>%
  filter(Freq > 0) %>%  # 过滤0值
  dplyr::rename(annotation = Var1, "celltype" = Var2) %>%
  arrange(annotation, desc(Freq))


