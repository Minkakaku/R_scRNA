library(tidyverse)
library(Seurat)
library(patchwork)


# setting variable   --------------------
dr_outdir <- "2. dimension reduction and cluster"  
dir.create(dr_outdir)
sp_output <- "3. spatial plot" 
dir.create(sp_output)


data_dir  <- "/data/xenium/FWSP20254331-03"
list.files(data_dir)
anno_file <- "data/Annotation.csv" # 非必须

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
  parseProjName_xenium(pattern = '^20251108', pattern2 = ".*?__.*?__.*?__.*?", extract_postion =  c(2)) %>%
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
adata <- sc$read_h5ad("/home/xjs/case/FWSP20254331-03/output/data/adata_celltype.h5ad")

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
      from = c( # project name
        'csy2_ln_zyj', 'csy2_ln_zyj', 'csy2_ln_zyj', 
        'lzf', 'sq3', 'wjy', 
        'wyt', 'ysg', 'ysp', 
        'zbl_ljf_gxh_ozr', 'zbl_ljf_gxh_ozr', 'zbl_ljf_gxh_ozr', 'zbl_ljf_gxh_ozr'),
      to   = c( # sample name
        "A3", "B2", "C1", 
        "C3", "B3", "B1", 
        "A2", "D3", "C4", 
        "A1", "C2", "D1", "D2")
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
all_data2 %>% DimPlot(reduction = "umap", group.by = "res_1", label = TRUE) + NoLegend()
all_data2 %>% DimPlot(reduction = "umap", group.by = "RNA_snn_res.0.8", label = TRUE) + NoLegend()




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


p1 <- cells_boundary_ls2[['A1']] %>%
  slice_sample(prop = 0.5) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()

p2 <- cells_boundary_ls2[['C2']] %>%
  slice_sample(prop = 0.5) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()

p3 <- cells_boundary_ls2[['D1']] %>%
  slice_sample(prop = 0.5) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()

p4 <- cells_boundary_ls2[['D2']] %>%
  slice_sample(prop = 0.5) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()

p <- cells_boundary_ls[['zbl_ljf_gxh_ozr']] %>%
  slice_sample(prop = 0.25) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()

library(patchwork)
((p2 | p4) / (p1 | p3))  / p




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
### pca inspector ------
pdf(file.path(dr_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_data2, ndims = 50)
dev.off()

all_data2@meta.data %>% head()

all_data2$cell_id <- NULL

DimPlot(all_data2, group.by = "sample", label = TRUE, reduction = "umap") + NoLegend()




# ## dimension reduction -----------------------------------------------------------------
# dims <- 20
# 
# library(future)
# plan(multisession, workers = 10)
# options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
# 
# 
# #all_data2 <- RunUMAP(all_data2, dims = 1:dims, return.model = TRUE, n.neighbors = 30)
# 
# all_data2 <- all_data2 %>%
#   harmony::RunHarmony(
#     dims.use = 1:dims,
#     group.by.vars = "sample",
#     reduction.save = "harmony"
#   )
# 
# all_data2 <- RunUMAP(
#   all_data2,
#   dims = 1:dims,
#   reduction = "harmony",
#   reduction.name = "harmonyumap",
#   reduction.key  = "harmonyUMAP_", 
#   return.model = TRUE
# )
# 
# pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
# DimPlot(all_data2, group.by = "sample", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data2, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
# dev.off()
# 
# pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} group - harmony.pdf")), height = 6, width = 8)
# DimPlot(all_data2, group.by = "group", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data2, group.by = "group", label = TRUE, reduction = "harmonyumap") + NoLegend()
# dev.off()

# export 
saveRDS(all_data2, file = file.path("R_OUTPUT", "all_data.rds"))

DimPlot(all_data2, group.by = "sample", label = TRUE, reduction = "umap")
ggsave(file.path(dr_dir, 'Dimplot_sample.pdf'), width = 8, height = 6)

