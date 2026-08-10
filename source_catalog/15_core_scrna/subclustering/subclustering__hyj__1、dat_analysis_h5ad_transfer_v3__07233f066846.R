library(tidyverse)
library(Seurat)
library(patchwork)


# setting variable   --------------------
data_dir  <- "/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/XE2511036121-01/XE2511036121-01"
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
#use_condaenv(condaenv = "/home/szt/miniforge3/envs/rpy/bin/python")
use_condaenv(condaenv = "/home/liuran/software/miniconda3/envs/rpy/bin/python")
py_config()
py_available()
py_module_available("pandas")


sc <- import("scanpy")
adata <- sc$read_h5ad("/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/XE2511036121-01/output/data/adata_celltype.h5ad")

adata_meta <- adata$obs
adata_umap <- adata$obsm['X_umap']

#from = c(rep("TMA2", 24),rep("TM1", 24)) 

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
      from = c(rep("TMA2", 24),rep("TM1", 24)) ,
      to <- c("B1_1","B1_2","B1_3","B1_4","B1_5","F1_1","F1_2","F1_3","F1_4","F1_5","F2_1","F2_2","F2_3","F2_4","F2_5","C6_1","C6_2","C6_3","C6_4","C6_5","A1_1","A1_2","A1_3","A1_4","A1_5","D2_1","D2_2","D2_3","D2_4","D2_5","A3_1","A3_2","A3_3","A3_4","A3_5","C6_1R","C6_4R","C6_5R","F2_1R","F2_2R","F2_4R","F1_3R","F1_4R","F1_5R","B1_5R","A1_2R","D2_3R","A3_2R") 
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
all_data2 %>% DimPlot(reduction = "umap", group.by = "res_0.8", label = TRUE) + NoLegend()
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


p1 <- cells_boundary_ls2[['F2_2R']] %>%
  slice_sample(prop = 0.5) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()

p2 <- cells_boundary_ls2[['A3_2']] %>%
  slice_sample(prop = 0.5) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()


p <- cells_boundary_ls[['TM1']] %>%
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
dims <- 20

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



#spatial plot -----
rownames(all_data2@meta.data) %>% head()
# all_data2@meta.data$cell_id %>% head()
# all_data2@meta.data$cell_id_raw <- all_data2@meta.data$cell_id
all_data2@meta.data$cell_id <- NULL

all_data_meta <- all_data2@meta.data %>%
  rownames_to_column("cell_id")

sp_p_ls <- 
  cells_boundary_ls2 %>% 
  #.[1] %>%
  .[names(all_data_sample_ls)] %>%
  imap(function(cells_boundary, project){ #browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    dat <- 
      cells_boundary %>% 
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>% 
      ggplot(aes(x = vertex_x, y = vertex_y)) + 
      geom_polygon(aes(fill = annotation, group = cell_id)) + 
      scale_y_reverse() + 
      coord_equal() + 
      #scale_fill_manual(values = plan_rctd) + 
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5, color = "black"))) + 
      theme_void()
    p
  })

spatial_dir <- "result/spatial_plot"
dir.create(spatial_dir)

# export
sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path("result/spatial_plot", str_glue("spatial plot annotation {grp}.pdf")), height = 18, width = 18)
    print(p)
    dev.off()
  })




# all_data_sample_ls <- SplitObject(all_data2, split.by = "sample")
# all_data_sample_ls$F2_2R@meta.data$cell_id %>% head()
# saveRDS(all_data_sample_ls,"R_OUTPUT/all_data_sample_ls.rds")

# sp_p_ls <-
#   cells_boundary_ls2 %>%
#   .[names(all_data_sample_ls)] %>%
#   imap(function(cells_boundary, project){ #browser()
#     #image_wh_r <- image_wh_r_ls[[project]]
#     #all_data_meta <- all_data_sample_ls[[project]][[]] %>% rownames_to_column('cell_id')
#     
#     all_data_meta <- all_data2[ ,all_data2$sample == project][[]] %>%
#       rownames_to_column("cell_id")
# 
#     dat <-
#       cells_boundary %>%
#       right_join(all_data_meta, by = "cell_id")
#     
#     p <- dat %>%
#       ggplot(aes(x = vertex_x, y = vertex_y)) +
#       geom_polygon(aes(fill = annotation, group = cell_id))  +
#       scale_y_reverse() +
#       #coord_fixed(ratio = 1/image_wh_r) +
#       coord_equal()+
#       ggsci::scale_fill_igv() +
#       # scale_fill_gradient(low = "gray90", high = "red") +
#       theme_void() +
#       labs(
#         title = project
#       )
#     p
#   })
# 
# 
# spatial_dir <- "result/spatial_plot"
# dir.create(spatial_dir)
# 
# # export
# sp_p_ls %>% 
#   iwalk(function(p, grp){
#     pdf(file.path("result/spatial_plot", str_glue("spatial plot annotation {grp}.pdf")), height = 18, width = 18)
#     print(p)
#     dev.off()
#   })




