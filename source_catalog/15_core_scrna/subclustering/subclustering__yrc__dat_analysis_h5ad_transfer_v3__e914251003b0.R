library(tidyverse)
library(Seurat)
library(patchwork)


# setting variable   --------------------
dr_outdir <- "2. dimension reduction and cluster"  # dir.create(dr_outdir)
sp_output <- "3. spatial plot" # dir.create(sp_output)


data_dir  <- "/data/xenium/XE2511046228-04"
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
use_condaenv(condaenv = "/home/szt/miniforge3/envs/rpy/bin/python")
py_config()
py_available()
py_module_available("pandas")


sc <- import("scanpy")
adata <- sc$read_h5ad("/mnt/NexData_file/home/yrc/XE2511046228/output/data/adata_celltype.h5ad")

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
        'TMA2','TMA2','TMA2','TMA2','TMA2','TMA2','TMA2','TMA2','TMA2','TMA2','TMA2','TMA2',
        'TMA1', 'TMA1', 'TMA1', 'TMA1', 'TMA1', 'TMA1', 'TMA1', 'TMA1', 'TMA1', 'TMA1', 'TMA1', 'TMA1'
        ),
      to   = c( # sample name
        "Con1", "Con2", "Con3", "Con4","Con5","Con6","Con7","Con8","Con9","Con10","Con11","Con12", 
        "CR1", "CR2", "NearCR3_1", "NearCR3_2","NearCR4","NR1","NR2","NR3","NR4","PR1","PR2","PR3"
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


p1 <- cells_boundary_ls2[['CR1']] %>%
  slice_sample(prop = 0.5) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()

p2 <- cells_boundary_ls2[['CR2']] %>%
  slice_sample(prop = 0.5) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()


p <- cells_boundary_ls[['TMA1']] %>%
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
ggsave(file.path(dr_dir, 'Dimplot_sample.pdf'), width = 6, height = 8)


all_data2@meta.data$group1 <- "other"
all_data2@meta.data$group2 <- "other"
all_data2@meta.data$group3 <- all_data2$group

all_data2@meta.data[all_data2$group%in%c("Con"),]$group1 <- "Con"
all_data2@meta.data[all_data2$group%in%c("NR","PR","CR"),]$group1 <- "Tre"

all_data2@meta.data[all_data2$group%in%c("Con"),]$group2 <- "Con"
all_data2@meta.data[all_data2$group%in%c("NR"),]$group2 <- "Non-Responder"
all_data2@meta.data[all_data2$group%in%c("PR","CR"),]$group2 <- "Responder"

all_data2$group1 <- factor(all_data2$group1,levels = c("Con","Tre"))
all_data2$group2 <- factor(all_data2$group2,levels = c("Con","Non-Responder","Responder"))
all_data2$group3 <- factor(all_data2$group3,levels = c("Con","NR","PR","CR"))
all_data2$sample <- factor(all_data2$sample,levels = c("Con1","Con2","Con3","Con4","Con5","Con6","Con7","Con8","Con9","Con10","Con11","Con12",
                                                       "NR1","NR2","NR3","NR4","PR1","PR2","PR3","CR1","CR2","NearCR3_1","NearCR3_2","NearCR4"))
saveRDS(all_data2, file = file.path("R_OUTPUT", "ALL.RDS"))

DimPlot(all_data2,group.by = "celltype")
ggsave("/mnt/NexData_file/home/yrc/XE2511046228/R_OUTPUT/umap_celltype.png",width = 6,height = 5)
ggsave("/mnt/NexData_file/home/yrc/XE2511046228/R_OUTPUT/umap_celltype.pdf",width = 6,height = 5)

DimPlot(all_data2,group.by = "celltype",label = T)+NoLegend()
ggsave("/mnt/NexData_file/home/yrc/XE2511046228/R_OUTPUT/umap_celltype_label.png",width = 8,height = 8)
ggsave("/mnt/NexData_file/home/yrc/XE2511046228/R_OUTPUT/umap_celltype_label.pdf",width = 8,height = 8)

DimPlot(all_data2,group.by = "sample")
ggsave("/mnt/NexData_file/home/yrc/XE2511046228/R_OUTPUT/umap_sample.png",width = 7,height = 5)
ggsave("/mnt/NexData_file/home/yrc/XE2511046228/R_OUTPUT/umap_sample.pdf",width = 7,height = 5)

DimPlot(all_data2,group.by = "group1")
ggsave("/mnt/NexData_file/home/yrc/XE2511046228/R_OUTPUT/umap_group1.png",width = 6,height = 5)
ggsave("/mnt/NexData_file/home/yrc/XE2511046228/R_OUTPUT/umap_group1.pdf",width = 6,height = 5)
DimPlot(all_data2,split.by = "group1")
ggsave("/mnt/NexData_file/home/yrc/XE2511046228/R_OUTPUT/umap_group1_split.png",width = 10,height = 5)
ggsave("/mnt/NexData_file/home/yrc/XE2511046228/R_OUTPUT/umap_group1_split.pdf",width = 10,height = 5)

DimPlot(all_data2,group.by = "group2")
ggsave("/mnt/NexData_file/home/yrc/XE2511046228/R_OUTPUT/umap_group2.png",width = 6,height = 5)
ggsave("/mnt/NexData_file/home/yrc/XE2511046228/R_OUTPUT/umap_group2.pdf",width = 6,height = 5)
DimPlot(all_data2,split.by =  "group2")
ggsave("/mnt/NexData_file/home/yrc/XE2511046228/R_OUTPUT/umap_group2_split.png",width = 15,height = 5)
ggsave("/mnt/NexData_file/home/yrc/XE2511046228/R_OUTPUT/umap_group2_split.pdf",width = 15,height = 5)

DimPlot(all_data2,group.by = "group3")
ggsave("/mnt/NexData_file/home/yrc/XE2511046228/R_OUTPUT/umap_group3.png",width = 6,height = 5)
ggsave("/mnt/NexData_file/home/yrc/XE2511046228/R_OUTPUT/umap_group3.pdf",width = 6,height = 5)
DimPlot(all_data2,split.by =  "group3")
ggsave("/mnt/NexData_file/home/yrc/XE2511046228/R_OUTPUT/umap_group3_split.png",width = 20,height = 5)
ggsave("/mnt/NexData_file/home/yrc/XE2511046228/R_OUTPUT/umap_group3_split.pdf",width = 20,height = 5)


####堆叠图####
p_col_proj <-
  all_data2@meta.data %>%
  group_by(group3, celltype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group3, y = freq, fill = celltype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave(
  p_col_proj,
  filename = file.path("R_OUTPUT","prop_group3.png"),
  height = 6,
  width  = 4
)
ggsave(
  p_col_proj,
  filename = file.path("R_OUTPUT","prop_group3.pdf"),
  height = 6,
  width  = 4
)

p_col_proj <-
  all_data2@meta.data %>%
  group_by(sample, celltype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = celltype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave(
  p_col_proj,
  filename = file.path("R_OUTPUT","prop_sample.png"),
  height = 6,
  width  = 18
)
ggsave(
  p_col_proj,
  filename = file.path("R_OUTPUT","prop_sample.pdf"),
  height = 6,
  width  = 18
)

####箱线图####
all_data2$group <- all_data2$group3

Ratio <- all_data2@meta.data %>%
  group_by(sample,group,celltype,) %>%
  count() %>%
  group_by(sample) %>%
  mutate(Freq = n/sum(n)*100)

my_comparisons <- list(c("Con","NR"),c("Con","PR"),c("Con","CR"),c("NR","PR"),c("NR","CR"),c("PR","CR"))#差异比较组

library(ggpubr)
ggplot(Ratio,aes(group,Freq))+
  geom_boxplot(aes(fill = group))+
  geom_point(size=1)+
  scale_x_discrete(limits=levels(Ratio$group))+ #很重要！保证X轴的label和实际level顺序对应
  scale_fill_discrete()+
  stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, aes(label = ..p.signif..)) +  
  scale_y_continuous("Proportion(%)")+
  facet_grid(~celltype)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  xlab("")+
  theme(panel.grid = element_blank())

ggsave(filename = file.path("R_OUTPUT","prop_boxplot_group3.png"),height = 6,width  = 18)
ggsave(filename = file.path("R_OUTPUT","prop_boxplot_group3.pdf"),height = 6,width  = 18)


####ROE####
library(tidyverse)

# AT1 
# roe_out_dir <- "subcluster/AT1"
# 
# anno_column <- "annotation"
# grp_column  <- "project"
# 
# all_data_meta_roe <- 
#   sc_at1[[]] %>%
#   dplyr::filter(.data[[anno_column]]  != "doublet") %>%
#   mutate({{anno_column}} := fct_drop(  .data[[anno_column]] ))

roe_out_dir <- "./R_OUTPUT/roe_analysis/group3"
dir.create(roe_out_dir)

anno_column <- "celltype"
grp_column  <- "group3"

all_data_meta_roe <-
  all_data2[[]] %>%
  dplyr::filter(.data[[grp_column]] != "Con") %>%
  droplevels()


# Cell Preference
library(ComplexHeatmap)
# 生成细胞数目分布频数统计表observe.data
observe.data <- unclass(table(all_data_meta_roe[[anno_column]], all_data_meta_roe[[grp_column]]))
stopifnot(is.matrix(observe.data))

expected <-
  function(x){
    if(!is.matrix(x)){stop("Must be a matrix")}
    rtot <- margin.table(x, 1)
    ctot <- margin.table(x, 2)
    tot  <- margin.table(x)
    outer(rtot, ctot, "*")/tot
  }
expected.data <- expected(observe.data) # 生成期望频数表
plot.data <- observe.data/expected.data # 计算实际细胞数目与期望细胞数目的比值
# 绘图数据的列顺序
col.order <- ggplot2:::`%||%`(
  levels(all_data_meta_roe[[grp_column]]),
  sort(unique(all_data_meta_roe[[grp_column]]))
) 
plot.data <- plot.data[, col.order] # 对绘图数据列进行排序
plot.data <- plot.data[order(apply(plot.data, 1, which.max), decreasing = T), ] # 对绘图数据行进行排序
# 如果想深色色块在主对角线使用which.min，在副对角线则使用which.max

# 热图颜色设定
roe_cols <- structure(
  c("#FEE6CE", "#FDC08C", "#F5904B", "#E6550D"),
  names = c("1", "1.5", "2", ">3")
)

cell_fun = function(j, i, x, y, width, height, fill) { # browser()
  
  # print(c(j, i, x, y, width, height, fill))
  
  if(plot.data[i,j] < as.numeric(names(roe_cols)[1])) {
    grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = roe_cols[1], fill = roe_cols[1]))
    grid::grid.text("\u00B1", x = x, y = y)
  } else if(plot.data[i,j] < as.numeric(names(roe_cols)[2])) {
    grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = roe_cols[2], fill = roe_cols[2]))
    grid::grid.text("+", x = x, y = y)
  } else if(plot.data[i,j] < as.numeric(names(roe_cols)[3])) {
    grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = roe_cols[3], fill = roe_cols[3]))
    grid::grid.text("++", x = x, y = y)
  } else {
    grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = roe_cols[4], fill = roe_cols[4]))
    grid::grid.text("+++", x = x, y = y)
  }
  
}

# 生成热图
roi_hm <-
  ComplexHeatmap::Heatmap(
    plot.data, 
    cell_fun = cell_fun,
    cluster_rows = F, 
    cluster_columns = F,
    width = unit(5, "cm"), 
    height = unit(17, "cm"), 
    show_heatmap_legend = F
  )

# stat comparasion
chi2.T <- sum((observe.data - expected.data)^2/expected.data)
p_chisq <- pchisq(q = chi2.T, df = (nrow(observe.data)-1)*(ncol(observe.data)-1), lower.tail = FALSE)
print(p_chisq)

# 生成图例
lgd <-
  ComplexHeatmap::Legend(
    labels = names(roe_cols), 
    title = expression(R[o/e]), 
    legend_gp = grid::gpar(fill = roe_cols)
  )


png(file.path(roe_out_dir, "cellprop_Roe.png"), 
    width = 8, height = 8, units = "in", res = 300)
ComplexHeatmap::draw(roi_hm, annotation_legend_list = lgd)
dev.off()

pdf(file.path(roe_out_dir, "cellprop_Roe.pdf"), width = 8, height = 8)
draw(roi_hm, annotation_legend_list = lgd)
dev.off()


####空间定位####
library(tidyverse)
library(Seurat)

out_dir <- 'R_OUTPUT/spatial'
dir.create(out_dir)

all_data <- readRDS(file.path("R_OUTPUT", "ALL.RDS"))
all_data@meta.data %>% head()

position <- readRDS("./R_OUTPUT/data_for_analysis/cells_position_ls.rds")


all_coords <- bind_rows(position, .id = "sample")
all_coords <- all_coords %>% select(cell_id, x_centroid, y_centroid)

# 将 all_coords 转换为以 cell_id 为行名的矩阵
coords_matrix <- all_coords %>%
  column_to_rownames("cell_id") %>%
  select(x_centroid, y_centroid)

# 只保留 Seurat 对象中存在的细胞
coords_matrix <- coords_matrix[colnames(all_data), , drop = FALSE]

# 添加到 meta.data
all_data <- AddMetaData(all_data, coords_matrix)

plot_data <- all_data@meta.data %>%
  select(x_centroid, y_centroid, celltype, sample) %>%
  filter(!is.na(x_centroid))  # 去除可能缺失坐标的细胞


for (sm in unique(plot_data$sample)) {
  p <- ggplot(plot_data %>% filter(sample == sm), aes(x = x_centroid, y = y_centroid, color = celltype)) +
    geom_point(size = 0.2, alpha = 0.8) +
    theme_classic() +   # 白色背景，无网格线
    theme(
      axis.text.x = element_blank(),   # 去掉x轴数值
      axis.text.y = element_blank(),   # 去掉y轴数值
      axis.ticks = element_blank(),    # 去掉刻度线（可选）
      axis.line = element_blank(),     # 去掉默认的轴线（因为我们要加方框）
      panel.border = element_rect(color = "black", fill = NA, size = 1)) +
    labs(title = sm, x = "spatial_1", y = "spatial_2") +
    coord_fixed()+
    guides(color = guide_legend(override.aes = list(size = 3)))
  ggsave(file.path(out_dir, paste0("spatial_", sm, ".png")), p, width = 8, height = 6)
  ggsave(file.path(out_dir, paste0("spatial_", sm, ".pdf")), p, width = 8, height = 6)
}
