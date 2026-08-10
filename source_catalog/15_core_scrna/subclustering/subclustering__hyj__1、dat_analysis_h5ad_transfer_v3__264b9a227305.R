library(tidyverse)
library(Seurat)
library(patchwork)


# setting variable   --------------------
data_dir  <- "/data/xenium/FWSP20254307"
list.files(data_dir)
anno_file <- "data/Annotation.csv"

dir.create("R_OUTPUT")
dir.create("R_OUTPUT/data_for_analysis")



# library --------------
options(java.parameters = '-Xmx10g') # options('java.parameters')
library(rJava)
#library(RBioFormats)
# source("xenium.R")

#multi workers  -----
library(future)
plan("multisession", workers = 10, gc = TRUE)

options(future.globals.maxSize = 10 * 1024 ^ 3)

nbrOfFreeWorkers()


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
use_condaenv(condaenv = "/home/liuran/software/miniconda3/envs/rpy/bin/python")
py_config()
py_available()
py_module_available("pandas")


sc <- import("scanpy")
adata <- sc$read_h5ad("/mnt/NexData_file/home/hyj/Spatial_project/xenium_5k/FWSP20254307_yurongbing_20260414/output/data/adata_celltype.h5ad")

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
      from = c(rep("TMA-1", 8)) ,
      to <- c("Ctrl1","Ctrl2","CtrMSC1","CtrMSC2","IR1","IR2","IRMSC1","IRMSC2") 
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


p1 <- cells_boundary_ls2[['CtrMSC1']] %>%
  slice_sample(prop = 0.5) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()

p2 <- cells_boundary_ls2[['IRMSC1']] %>%
  slice_sample(prop = 0.5) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()


p <- cells_boundary_ls[['TMA-1']] %>%
  slice_sample(prop = 0.25) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()

library(patchwork)
(p2 | p1)  / p


saveRDS(all_data2, file = file.path("R_OUTPUT", "all_data.rds"))


# 
# # normalize ----------------------------------------------------------------------
# all_data2 <- NormalizeData(all_data2, scale.factor = median(all_data2$nCount_RNA))
# all_data2 <- FindVariableFeatures(all_data2)
# 
# # PCA ----------------------------------------------------------------------------
# ## PCA --------
# #DefaultAssay(all_data2) <- "sketch"
# #all_data2 <- FindVariableFeatures(all_data2)
# VariableFeaturePlot(all_data2)
# 
# all_data2 <- ScaleData(all_data2, vars.to.regress = 'nCount_RNA')
# all_data2 <- RunPCA(all_data2)
# 
# dr_dir <- '2. dimension reduction and cluster/'
# dir.create(dr_dir)
# 
# ### pca inspector ------
# pdf(file.path(dr_dir, 'pca plot.pdf'), height = 6, width = 8)
# ElbowPlot(all_data2, ndims = 50)
# dev.off()
# 
# all_data2@meta.data %>% head()
# 
# all_data2$cell_id <- NULL
# 
# set.seed(1234)
# suset_cells <-
#   all_data2[[]] %>%
#   rownames_to_column("cell_id") %>%
#   group_by(sample) %>%
#   slice_sample(prop = 50000/ncol(all_data2) ) %>%
#   dplyr::pull("cell_id")
# 
# pdf(file.path(dr_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
# DimHeatmap(all_data2[, suset_cells], dims = 1:5)
# dev.off()
# 
# pdf(file.path(dr_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
# DimHeatmap(all_data2[, suset_cells], dims = 6:10)
# dev.off()
# 
# pdf(file.path(dr_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
# DimHeatmap(all_data2[, suset_cells], dims = 11:15)
# dev.off()
# 
# pdf(file.path(dr_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
# DimHeatmap(all_data2[, suset_cells], dims = 16:20)
# dev.off()
# 
# 
# 
# 
# ### dimension reduction -----------------------------------------------------------------
# dims <- 30
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
# 
# # export 
# saveRDS(all_data2, file = file.path("R_OUTPUT", "all_data.rds"))
# 
# 
# # 生成列联表：统计group1下各sample的细胞数量
# table(all_data2@meta.data$annotation, all_data2@meta.data$"celltype") %>% 
#   as.data.frame() %>%
#   filter(Freq > 0) %>%  # 过滤0值
#   dplyr::rename(annotation = Var1, "celltype" = Var2) %>%
#   arrange(annotation, desc(Freq))

all_data2$celltype %>% unique()

#1、 找Lgr5+ 小肠干细胞-------------------------------------------------------------------
Epithelial <- subset(all_data2, annotation %in% c('Epithelial progenitor', 'Paneth', 'Intestinal Stem Cells'))

### PCA --------
DefaultAssay(Epithelial) <- "RNA"
Epithelial <- FindVariableFeatures(Epithelial)
Epithelial <- NormalizeData(Epithelial)
Epithelial <- ScaleData(Epithelial, vars.to.regress = 'nCount_RNA')
Epithelial <- RunPCA(Epithelial)

Epithelial_out_dir <-'result/Find_lgr5_Epithelial_subtype' 
dir.create(Epithelial_out_dir)


#### pca inspector ------
pdf(file.path(Epithelial_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Epithelial, ndims = 50)
dev.off()

# pdf(file.path(Epithelial_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
# DimHeatmap(Epithelial, dims = 1:5)
# dev.off()
# pdf(file.path(Epithelial_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
# DimHeatmap(Epithelial, dims = 6:10)
# dev.off()
# pdf(file.path(Epithelial_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
# DimHeatmap(Epithelial, dims = 11:15)
# dev.off()
# pdf(file.path(Epithelial_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
# DimHeatmap(Epithelial, dims = 16:20)
# dev.off()


### dimension reduction --------
dims <- 20

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


Epithelial <- RunUMAP(Epithelial, dims = 1:dims,return.model = TRUE)



pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "sample", label = TRUE) + NoLegend()
dev.off()


Epithelial <- Epithelial %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

Epithelial <- RunUMAP(
  Epithelial, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


#### cluster -------
Epithelial  <-
  FindNeighbors(
    Epithelial ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Epithelial <-
  FindClusters(
    Epithelial,
    resolution = c(0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
#DimPlot(Epithelial, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

pdf(file.path(Epithelial_out_dir, str_glue("FeaturePlot dim{dims} Lgr5.pdf")), height = 6, width = 8)
FeaturePlot(Epithelial,feature = c("Lgr5"),reduction = 'harmonyumap')
dev.off()


Epithelial$res_1.2 %>% table()
pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} res_1.2.pdf")), height = 6, width = 8)
#DimPlot(Epithelial, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "res_1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()

dev.off()

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} celltype.pdf")), height = 6, width = 8)
#DimPlot(Epithelial, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Epithelial, group.by = "celltype", label = TRUE, reduction = 'harmonyumap') + NoLegend()

dev.off()



#### dotplot -------
# 1. Epi2: Xiao S, Zhou L. Gastric Stem Cells: Physiological and Pathological Perspectives. Front Cell Dev Biol. 2020
# 2. Epi3: Kim TH et al., Stomach development, stem cells and disease. Development. 2016
# 3. Yaser Greish et al., Biodegradable Scaffolds for Gastric Tissue Regeneration. 2017
# 4. Epi4: Busslinger GA et al., Human gastrointestinal epithelial of the esophagus, stomach, and duodenum resolved at single-cell resolution. Cell Rep. 2021
# 5. Epi5: Engevik AC et al., The Physiology of the Gastric Parietal Cell. Physiol Rev. 2020
# 6. Epi6: Zhang P et al.,Dissecting the Single-Cell Transcriptome Network Underlying Gastric Premalignant Lesions and Early Gastric Cancer. Cell Rep. 2019
# 7. Epi7: Sathe A et al., Single-Cell Genomic Characterization Reveals the Cellular Reprogramming of the Gastric Tumor Microenvironment. Clin Cancer Res. 2020
# 8. Epi8: Zhang M et al., Dissecting transcriptional heterogeneity in primary gastric adenocarcinoma by single cell RNA sequencing. Gut. 2021
Epithelial_major_marker_ls <- c(
  'Surface mucous cells' = "Muc5ac", # 1
  'Surface mucous cells' = "Tff1", # 1 
  'Surface mucous cells' = "Pas", # 1 
  'Surface mucous cells' = "Gkn1", # 9 
  'Surface mucous cells' = "Gkn2", # 9 
  
  'Mucous gland cells' = "Muc6", # 1 # 5:Tff2 
  
  'Mucous neck cells' = "Tff2", # 1 
  
  'Parietal cells' = "Gif", # 1 
  'Parietal cells' = "Esrrg", # 2
  'Parietal cells' = "Atp4b", # 2
  
  'Enteroendocrine cells' = "Chga", # 3,6
  'Enteroendocrine cells' = "Chgb", # 6
  'Enteroendocrine cells' = "Gast", # 9
  'Enteroendocrine cells' = "Prox1", # 9
  
  'Acid-secreting parietal cells' = "Atp4a", # 4 
  #'Acid-secreting parietal cells' = "Atp4b", # 4 
  
  # 'Proliferating isthmus cells' = "Mki67", # 4 
  'Proliferating isthmus cells' = "Stmn1", # 4
  
  'Proliferative cells' = "Mki67", # 6
  'Proliferative cells' = "Bric5", # 6
  'Proliferative cells' = "Cdkn2a", # 9
  'Proliferative cells' = "Rbp4", # 9
  
  'Chief cells' = "Lipf", # 4
  'Chief cells' = "Pga5", # 4
  'Chief cells' = "Mist1", # 5
  'Chief cells' = "Pgc", # 5
  'Chief cells' = "Pga3", # 6
  'Chief cells' = "Pga4", # 6
  
  # 'Endocrine cells' = "Gast", # 5
  'Endocrine cells' = "Ghrl", # 5
  'Endocrine cells' = "Sst", # 5
  'Endocrine cells' = "Hdc", # 5
  
  'Tuft cells' = "Trpm5", # 5
  'Tuft cells' = "Sh2d6", # 5
  
  'Tumor cells' = "Ceacam5", # 6
  'Tumor cells' = "Ceacam6", # 6
  'Tumor cells' = "Epcam", # 9 → 
  'Tumor cells' = "Cdh17", # 9
  'Tumor cells' = "Mdm2", # P53 pathway
  
  'Enterocytes' = "Fabp1", # 6
  'Enterocytes' = "Apoa1", # 6
  'Enterocytes' = "Fabp2", # 8
  
  'Primary cancer cells' = "Notch2", # 7
  'Primary cancer cells' = "Notch2nl", # 7
  'Primary cancer cells' = "Kif5b", # 7
  'Primary cancer cells' = "Erbb4", # 7
  
  'Metastatic cancer cells' = "Cdk12", # 7
  'Metastatic cancer cells' = "Erbb2", # 7
  'Metastatic cancer cells' = "Cldn11", # 7 
  
  'Driver gene' = "Fos", # 7
  'Driver gene' = "Jun", # 7
  
  'Metaplasia genes' = "Cd44", # 8 癌前
  # 'Metaplasia genes' = "Tff2", # 8
  # 'Metaplasia genes' = "Sox9", # 8
  'Metaplasia genes' = "Hla-Dpa1", # 8 → 小鼠：H2-Dpa1
  'Metaplasia genes' = "Hla-Dra", # 8 → 小鼠：H2-Dra
  
  'Intestinal metaplasia cells' = "Tff3", # 9
  'Intestinal metaplasia cells' = "Cdx1", # 9
  'Intestinal metaplasia cells' = "Cdx2", # 9
  
  'Goblet cells' = "Muc2", # 6, 9
  'Goblet cells' = "Itln1", # 6, 9
  'Goblet cells' = "Hes6", # 9
  
  'Metaplastic-stem cells' = "Ephb2", # 6
  'Metaplastic-stem cells' = "Sox9", # 6  
  'Metaplastic-stem cells' = "Olfm4", # 6, 9
  'Metaplastic-stem cells' = "Reg1a", # 9
  'Metaplastic-stem cells' = "Cldn3" # 9
)

intestinal_marker_ls <- c(
  # 上皮细胞群
  'Epithelial Stem Cell' = "Lgr5",
  'Epithelial Stem Cell' = "Ascl2",
  'Epithelial Stem Cell' = "Sox9",
  
  "Epcam",
  
  'Paneth Cell' = "Mmp7",
  
  'Enterocyte' = "Slc26a3",
  'Enterocyte' = "Vil1",
  
  'Goblet Cell' = "Muc2",
  'Goblet Cell' = "Spdef",
  
  'Enteroendocrine Cell (EEC)' = "Neurod1",
  'Enteroendocrine Cell (EEC)' = "Tph1",
  'Enteroendocrine Cell (EEC)' = "Gcg",
  'Enteroendocrine Cell (EEC)' = "Sst",
  
  'Tuft Cell' = "Il25",
  'Tuft Cell' = "Trpm5",
  'Tuft Cell' = "Avil",
  
  'M Cell' = "Tnfrsf11a",
  
  # 免疫细胞群
  'CD4+ T Cell' = "Cd3d",
  'CD4+ T Cell' = "Cd3e",
  'CD4+ T Cell' = "Cd3g",
  'CD4+ T Cell' = "Cd4",
  'CD4+ T Cell' = "Foxp3",
  
  # 'CD8+ T Cell' = "Cd3d",
  # 'CD8+ T Cell' = "Cd3e",
  # 'CD8+ T Cell' = "Cd3g",
  'CD8+ T Cell' = "Cd8a",
  'CD8+ T Cell' = "Gzmb",
  
  'B Cell' = "Cd19",
  'B Cell' = "Cd79a",
  'B Cell' = "Ighm",
  
  'Plasma Cell' = "Sdc1",
  'Plasma Cell' = "Jchain",
  
  'Dendritic Cell (DC)' = "Itgax",
  'Dendritic Cell (DC)' = "Zbtb46",
  'Dendritic Cell (DC)' = "Clec9a",
  
  'Macrophage' = "Adgre1",
  'Macrophage' = "Csf1r",
  'Macrophage' = "Cd68",
  'Macrophage' = "Cx3cr1",
  
  'Innate Lymphoid Cell (ILC)' = "Il7r",
  'Innate Lymphoid Cell (ILC)' = "Id2",
  'Innate Lymphoid Cell (ILC)' = "Rorc",
  'Innate Lymphoid Cell (ILC)' = "Gata3",
  
  # 基质细胞群
  'Fibroblast' = "Col1a1",
  'Fibroblast' = "Pdgfra",
  'Fibroblast' = "Pdgfrb",
  
  'Endothelial Cell' = "Pecam1",
  'Endothelial Cell' = "Cdh5",
  'Endothelial Cell' = "Lyve1",
  
  'Pericyte' = "Cspg4",
  #'Pericyte' = "Pdgfrb",
  
  'Enteric Glial Cell (EGC)' = "Gfap",
  'Enteric Glial Cell (EGC)' = "S100b"
)

intestinal_marker_ls1 <- c(
  # Intestinal Stem Cells (小肠干细胞)
  'Intestinal Stem Cell Core Markers' = "Ascl2",
  'Intestinal Stem Cell Core Markers' = "Axin2",
  'Intestinal Stem Cell Core Markers' = "Cdk6",
  'Intestinal Stem Cell Core Markers' = "Edn1",
  'Intestinal Stem Cell Core Markers' = "Ephb2",
  'Intestinal Stem Cell Core Markers' = "Ephb3",
  'Intestinal Stem Cell Core Markers' = "Lgr5",
  'Intestinal Stem Cell Core Markers' = "Lrig1",
  'Intestinal Stem Cell Core Markers' = "Mecom",
  'Intestinal Stem Cell Core Markers' = "Msi1",
  'Intestinal Stem Cell Core Markers' = "Rgmb",
  'Intestinal Stem Cell Core Markers' = "Smoc2",
  'Intestinal Stem Cell Core Markers' = "Sox9",
  
  # ISC Proliferation Pathway Genes (小肠干细胞增殖通路相关基因)
  'ISC Proliferation Pathway' = "Acvr1",
  'ISC Proliferation Pathway' = "Acvr2a",
  'ISC Proliferation Pathway' = "Acvr2b",
  'ISC Proliferation Pathway' = "Akt1",
  'ISC Proliferation Pathway' = "Akt2",
  'ISC Proliferation Pathway' = "Apc",
  'ISC Proliferation Pathway' = "Axin1",
  #'ISC Proliferation Pathway' = "Axin2",
  'ISC Proliferation Pathway' = "Bmp2",
  'ISC Proliferation Pathway' = "Bmp4",
  'ISC Proliferation Pathway' = "Bmp5",
  'ISC Proliferation Pathway' = "Bmp7",
  'ISC Proliferation Pathway' = "Bmpr1a",
  'ISC Proliferation Pathway' = "Bmpr1b",
  'ISC Proliferation Pathway' = "Bmpr2",
  'ISC Proliferation Pathway' = "Ccnd1",
  'ISC Proliferation Pathway' = "Ccnd2",
  'ISC Proliferation Pathway' = "Cdk4",
  #'ISC Proliferation Pathway' = "Cdk6",
  'ISC Proliferation Pathway' = "Cdkn1a",
  'ISC Proliferation Pathway' = "Cdkn1b",
  'ISC Proliferation Pathway' = "Cdkn2a",
  'ISC Proliferation Pathway' = "Ctnnb1",
  'ISC Proliferation Pathway' = "Dhh",
  'ISC Proliferation Pathway' = "Dll1",
  'ISC Proliferation Pathway' = "Dll4",
  'ISC Proliferation Pathway' = "Egfr",
  'ISC Proliferation Pathway' = "Erbb2",
  'ISC Proliferation Pathway' = "Erbb3",
  'ISC Proliferation Pathway' = "Erbb4",
  'ISC Proliferation Pathway' = "Ereg",
  'ISC Proliferation Pathway' = "Fgf1",
  'ISC Proliferation Pathway' = "Fgf10",
  'ISC Proliferation Pathway' = "Fgf2",
  'ISC Proliferation Pathway' = "Fgf4",
  'ISC Proliferation Pathway' = "Fgf7",
  'ISC Proliferation Pathway' = "Fgfr1",
  'ISC Proliferation Pathway' = "Fgfr2",
  'ISC Proliferation Pathway' = "Fgfr3",
  'ISC Proliferation Pathway' = "Fgfr4",
  'ISC Proliferation Pathway' = "Fzd1",
  'ISC Proliferation Pathway' = "Fzd10",
  'ISC Proliferation Pathway' = "Fzd2",
  'ISC Proliferation Pathway' = "Fzd3",
  'ISC Proliferation Pathway' = "Fzd4",
  'ISC Proliferation Pathway' = "Fzd5",
  'ISC Proliferation Pathway' = "Fzd6",
  'ISC Proliferation Pathway' = "Fzd7",
  'ISC Proliferation Pathway' = "Fzd8",
  'ISC Proliferation Pathway' = "Gli1",
  'ISC Proliferation Pathway' = "Gli2",
  'ISC Proliferation Pathway' = "Gli3",
  'ISC Proliferation Pathway' = "Gsk3b",
  'ISC Proliferation Pathway' = "Hes1",
  'ISC Proliferation Pathway' = "Hes5",
  'ISC Proliferation Pathway' = "Id1",
  'ISC Proliferation Pathway' = "Id2",
  'ISC Proliferation Pathway' = "Ihh",
  'ISC Proliferation Pathway' = "Jag1",
  'ISC Proliferation Pathway' = "Jag2",
  'ISC Proliferation Pathway' = "Lats1",
  'ISC Proliferation Pathway' = "Lats2",
  'ISC Proliferation Pathway' = "Lrp5",
  'ISC Proliferation Pathway' = "Lrp6",
  'ISC Proliferation Pathway' = "Mapk1",
  'ISC Proliferation Pathway' = "Mapk14",
  'ISC Proliferation Pathway' = "Mapk3",
  'ISC Proliferation Pathway' = "Mapk8",
  'ISC Proliferation Pathway' = "Mst1",
  'ISC Proliferation Pathway' = "Mtor",
  'ISC Proliferation Pathway' = "Myc",
  'ISC Proliferation Pathway' = "Notch1",
  'ISC Proliferation Pathway' = "Notch2",
  'ISC Proliferation Pathway' = "Notch3",
  'ISC Proliferation Pathway' = "Notch4",
  'ISC Proliferation Pathway' = "Pik3ca",
  'ISC Proliferation Pathway' = "Ptch1",
  'ISC Proliferation Pathway' = "Rbpj",
  'ISC Proliferation Pathway' = "Rspo1",
  'ISC Proliferation Pathway' = "Rspo2",
  'ISC Proliferation Pathway' = "Rspo3",
  'ISC Proliferation Pathway' = "Shh",
  'ISC Proliferation Pathway' = "Smad1",
  'ISC Proliferation Pathway' = "Smad2",
  'ISC Proliferation Pathway' = "Smad3",
  'ISC Proliferation Pathway' = "Smad4",
  'ISC Proliferation Pathway' = "Smad5",
  'ISC Proliferation Pathway' = "Smo",
  'ISC Proliferation Pathway' = "Tgfa",
  'ISC Proliferation Pathway' = "Tgfb1",
  'ISC Proliferation Pathway' = "Tgfb2",
  'ISC Proliferation Pathway' = "Tgfb3",
  'ISC Proliferation Pathway' = "Tgfbr1",
  'ISC Proliferation Pathway' = "Tgfbr2",
  'ISC Proliferation Pathway' = "Tgfbr3",
  'ISC Proliferation Pathway' = "Wnt3",
  'ISC Proliferation Pathway' = "Wnt3a",
  'ISC Proliferation Pathway' = "Wnt5a",
  'ISC Proliferation Pathway' = "Wnt5b",
  'ISC Proliferation Pathway' = "Yap1"
)


# find markers
Idents(Epithelial) <- "RNA_snn_res.0.4"
Epithelial_markers <-
  FindAllMarkers(
    JoinLayers(Epithelial),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

Epithelial_top5_markers <- 
  Epithelial_markers %>% 
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

names(Epithelial_top5_markers) <- sub('.*\\.','',names(Epithelial_top5_markers))

Epithelial_top10_markers <- 
  Epithelial_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10) 



# dotplot
pdf(file.path(Epithelial_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 35, height = 6)
#DotPlot(Epithelial, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()
DotPlot(Epithelial, group.by = "RNA_snn_res.0.4", features = Epithelial_major_marker_ls)+RotatedAxis()
DotPlot(Epithelial, group.by = "RNA_snn_res.0.4", features = intestinal_marker_ls)+RotatedAxis()
DotPlot(Epithelial, group.by = "RNA_snn_res.0.4", features = intestinal_marker_ls1)+RotatedAxis()
DotPlot(Epithelial, group.by = "RNA_snn_res.0.4", features = Epithelial_top5_markers)+RotatedAxis()
dev.off()

# DotPlot(Epithelial, group.by = "RNA_snn_res.0.8", features = Fibro_major_marker_ls2)+RotatedAxis()
VlnPlot(Epithelial, features = 'nCount_RNA', pt.size = 0)


### annotation -----------
cell_cluster <- c(
  '0' = "Proliferate Epithlial",
  '1' = "delete",
  '2' = "Enterocytes",
  '3' = "Lgr5 Intestinal Stem Cell",
  '4' = "Enterocytes",
  '5' = "Panth",
  '6' = "Goblet",
  '7' = "delete",
  '8' = "Proliferate Epithlial",
  '9' = "delete",
  '10' = "Panth",
  '11' = "delete",
  '12' = "delete"
)


Idents(Epithelial) <- "RNA_snn_res.0.4"

Epithelial <- RenameIdents(Epithelial, cell_cluster)
Epithelial$subtype <- Idents(Epithelial)

pdf(file.path(Epithelial_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Epithelial, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Epithelial_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 7, width = 18)
DimPlot(Epithelial, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


### cell props barplot  --------
Epithelial_meta <- Epithelial[[]]

Epithelial_meta$sample <- factor(Epithelial_meta$sample, 
                                 levels=c(Epithelial_meta$sample %>% unique())
)
p_col_proj <-
  Epithelial_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(Epithelial_out_dir,"cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  Epithelial_meta %>%
  group_by(group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path(Epithelial_out_dir,"cell_freq_group.pdf"),
  height = 6,
  width  = 4
)

saveRDS(Epithelial, file.path(Epithelial_out_dir, 'Epithelial.rds'))


##加入 lgr5 Item 到 all_data2----
library(tibble)
library(data.table)
sub <- Epithelial[['subtype']] %>% rownames_to_column('cells')

all_data2$annotation %>% table()
meta.data.t <-
  sub%>%
  right_join(
    all_data2[['annotation']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation_Lgr5 = ifelse((as.character(annotation) %in% c('Epithelial progenitor', 'Paneth', 'Intestinal Stem Cells')),
                         as.character(subtype), 
                         as.character(annotation)))%>%
  dplyr::select(cells, annotation_Lgr5) %>%
  column_to_rownames('cells')

meta.data.t$annotation_Lgr5 %>% table()

all_data2 <- AddMetaData(all_data2,meta.data.t['annotation_Lgr5'])

saveRDS(all_data2, 'result/all_data_Lgr5.rds')



cellType_marker_ls <- c(
  'Intestinal metaplasia cells' = "CDX1", # 9
  'Intestinal metaplasia cells' = "CDX2", # 9
  
  'Surface mucous cells' = "MUC5AC", # 1
  'Surface mucous cells' = "GKN1", # 9 
  
  'Proliferating genes' = "STMN1", # 4
  'Proliferating genes' = "MKI67", # 6
  
  'Metastatic cancer cells' = "ERBB2", # 7
  
  'Tumor cells' = "EPCAM", # 9
  'Tumor cells' = "MDM2", # P53 pathway
  
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "SFRP4",
  
  'myCAF' = 'MMP11',
  'myCAF' = 'CTHRC1',
  'myCAF' = 'COL11A1',
  'myCAF' = 'POSTN',
  'myCAF' = 'COL10A1',
  'myCAF' = 'TNC',
  
  'Chief cells' = "LIPF", # 4
  'Chief cells' = "PGA3", # 6
  
  'T cells' = 'CD3E', 
  'T cells' = 'CD3G', 
  'CD8_Tex'='HAVCR2',
  'CD8_Tex'='CTLA4',
  'CD8_Tex'='LAG3',
  'CD8_Tex'='TIGIT',
  
  'Goblet cells' = "MUC2", # 6, 9
  'Goblet cells' = "ITLN1", # 6, 9
  'Goblet cells' = "HES6", # 9
  
  'Enteroendocrine cells' = "PROX1", # 9
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "Plasma" = "MZB1",
  "Plasma" = "XBP1",
  
  'Primary cancer cells' = "NOTCH2", # 7
  'Primary cancer cells' = "KIF5B", # 7
  
  'Metaplasia genes' = "CD44", # 8 癌前
  
  'Metaplastic-stem cells' = "EPHB2", # 6
  'Metaplastic-stem cells' = "SOX9" # 6  
)

pdf(file.path(Epithelial_out_dir, "dotplot celltype main_marker.pdf"), width = 20, height = 8)
DotPlot(Epithelial, group.by = "subtype", features = cellType_marker_ls, dot.scale = 10) + theme_90angle
dev.off()




## volcano-------
Epithelial_markers <-
  FindAllMarkers(
    JoinLayers(Epithelial),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(Epithelial_markers,file.path(Epithelial_out_dir,"Epithelial_diff_gene.csv"))

library(EnhancedVolcano)
type<-unique(Epithelial_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- Epithelial_markers[Epithelial_markers$cluster == type, ]
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
      filename = file.path(Epithelial_out_dir, str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })

Epithelial_markers_ls <- list(Epithelial_markers)
names(Epithelial_markers_ls) <- 'anno cells'

# scatter
p_volcanoplot_major_ls <- 
  Epithelial_markers_ls %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 0.5
    dodge_width <- 0.1
    
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
      filename = file.path(Epithelial_out_dir, str_glue("{name} diffgene dot.pdf")),
      height = 10,
      width = 20
    )
  })




## enrichment ----------------------------
Epithelial_de_genes_for_enrichmenet <- 
  Epithelial_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.15) %>%
  dplyr::filter(abs(avg_log2FC) > 0.5) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")


enrich_ls <- 
  Epithelial_de_genes_for_enrichmenet %>%
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


er_output<- file.path(Epithelial_out_dir,'enrichment') 
dir.create(er_output)

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(er_output, name))
    
    write.csv(
      x@result,
      file = file.path(er_output, name, "Go enrichment table.csv"),
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
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, name, "go dotplot.pdf"),
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
      filename = file.path(er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })

Epithelial_de_genes_id_for_enrichmenet <-
  Epithelial_de_genes_for_enrichmenet %>%
  map(function(x){ # browser()
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
  Epithelial_de_genes_id_for_enrichmenet %>%
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

# saveRDS(enrich_ls, file.path(er_output, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(er_output, name))
    
    write.csv(
      x@result,
      file = file.path(er_output, name, "KEGG enriment table.csv"),
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
      filename = file.path(er_output, name, "kegg dotplot.pdf"),
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
      filename = file.path(er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })


all_data2$annotation %>% table()
#2、找骨髓来源间充质细胞-------
Fibroblast3_out_dir<-'result/Find_BM-MSC_Fibroblast_subtype/' 
dir.create(Fibroblast3_out_dir)
Fibroblast3<-subset(all_data2, subset = annotation %in% c("BM_MSCs","Fibroblast","Myofibroblast"))
Fibroblast3$annotation %>% table()

## PCA --------
DefaultAssay(Fibroblast3) <- "RNA"
Fibroblast3 <- FindVariableFeatures(Fibroblast3)
Fibroblast3 <- NormalizeData(Fibroblast3)
Fibroblast3 <- ScaleData(Fibroblast3, vars.to.regress = 'nCount_RNA')
Fibroblast3 <- RunPCA(Fibroblast3)


### pca inspector ------
pdf(file.path(Fibroblast3_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Fibroblast3, ndims = 50)
dev.off()

# pdf(file.path(Fibroblast3_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
# DimHeatmap(Fibroblast3, dims = 1:5)
# dev.off()
# pdf(file.path(Fibroblast3_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
# DimHeatmap(Fibroblast3, dims = 6:10)
# dev.off()
# pdf(file.path(Fibroblast3_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
# DimHeatmap(Fibroblast3, dims = 11:15)
# dev.off()
# pdf(file.path(Fibroblast3_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
# DimHeatmap(Fibroblast3, dims = 16:20)
# dev.off()


## dimension reduction --------
dims <- 20

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb

Fibroblast3 <- RunUMAP(Fibroblast3, dims = 1:dims,return.model = TRUE)

pdf(file.path(Fibroblast3_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Fibroblast3, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

table(all_data2@meta.data$celltype, all_data2@meta.data$project)

Fibroblast3 <- Fibroblast3 %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

Fibroblast3 <- RunUMAP(
  Fibroblast3, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(Fibroblast3_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Fibroblast3, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'

### cluster -------
Fibroblast3  <-
  FindNeighbors(
    Fibroblast3 ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )

Fibroblast3 <-
  FindClusters(
    Fibroblast3,
    resolution = c(0.2, 0.3, 0.4)
  )

# dimplot
pdf(file.path(Fibroblast3_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Fibroblast3, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast3, group.by = "RNA_snn_res.0.3", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast3, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


pdf(file.path(Fibroblast3_out_dir, str_glue("FeaturePlot dim{dims} Pdgfra.pdf")), height = 6, width = 8)
FeaturePlot(Fibroblast3,feature = c("Pdgfra"),reduction = 'harmonyumap')
dev.off()


Fibroblast3$res_1.2 %>% table()
pdf(file.path(Fibroblast3_out_dir, str_glue("Dimplot dim{dims} res_1.2.pdf")), height = 6, width = 8)

DimPlot(Fibroblast3, group.by = "res_1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()

dev.off()

pdf(file.path(Fibroblast3_out_dir, str_glue("Dimplot dim{dims} celltype.pdf")), height = 6, width = 8)

DimPlot(Fibroblast3, group.by = "celltype", label = TRUE, reduction = 'harmonyumap') + NoLegend()

dev.off()




### dotplot -------
intestinal_marker_ls <- c(
  # 上皮细胞群
  'Epithelial Stem Cell' = "Lgr5",
  'Epithelial Stem Cell' = "Ascl2",
  'Epithelial Stem Cell' = "Sox9",
  
  "Epcam",
  
  'Paneth Cell' = "Mmp7",
  
  'Enterocyte' = "Slc26a3",
  'Enterocyte' = "Vil1",
  
  'Goblet Cell' = "Muc2",
  'Goblet Cell' = "Spdef",
  
  'Enteroendocrine Cell (EEC)' = "Neurod1",
  'Enteroendocrine Cell (EEC)' = "Tph1",
  'Enteroendocrine Cell (EEC)' = "Gcg",
  'Enteroendocrine Cell (EEC)' = "Sst",
  
  'Tuft Cell' = "Il25",
  'Tuft Cell' = "Trpm5",
  'Tuft Cell' = "Avil",
  
  'M Cell' = "Tnfrsf11a",
  
  # 免疫细胞群
  'CD4+ T Cell' = "Cd3d",
  'CD4+ T Cell' = "Cd3e",
  'CD4+ T Cell' = "Cd3g",
  'CD4+ T Cell' = "Cd4",
  'CD4+ T Cell' = "Foxp3",
  
  # 'CD8+ T Cell' = "Cd3d",
  # 'CD8+ T Cell' = "Cd3e",
  # 'CD8+ T Cell' = "Cd3g",
  'CD8+ T Cell' = "Cd8a",
  'CD8+ T Cell' = "Gzmb",
  
  'B Cell' = "Cd19",
  'B Cell' = "Cd79a",
  'B Cell' = "Ighm",
  
  'Plasma Cell' = "Sdc1",
  'Plasma Cell' = "Jchain",
  
  'Dendritic Cell (DC)' = "Itgax",
  'Dendritic Cell (DC)' = "Zbtb46",
  'Dendritic Cell (DC)' = "Clec9a",
  
  'Macrophage' = "Adgre1",
  'Macrophage' = "Csf1r",
  'Macrophage' = "Cd68",
  'Macrophage' = "Cx3cr1",
  
  'Innate Lymphoid Cell (ILC)' = "Il7r",
  'Innate Lymphoid Cell (ILC)' = "Id2",
  'Innate Lymphoid Cell (ILC)' = "Rorc",
  'Innate Lymphoid Cell (ILC)' = "Gata3",
  
  # 基质细胞群
  'Fibroblast' = "Col1a1",
  'Fibroblast' = "Pdgfra",
  'Fibroblast' = "Pdgfrb",
  
  'Endothelial Cell' = "Pecam1",
  'Endothelial Cell' = "Cdh5",
  'Endothelial Cell' = "Lyve1",
  
  'Pericyte' = "Cspg4",
  #'Pericyte' = "Pdgfrb",
  
  'Enteric Glial Cell (EGC)' = "Gfap",
  'Enteric Glial Cell (EGC)' = "S100b"
)

msc_marker_ls <- c(
  # MSC 经典标记基因
  'MSC Classic Markers' = "Cd248",
  'MSC Classic Markers' = "Cd34",
  'MSC Classic Markers' = "Cd44",
  'MSC Classic Markers' = "Eng",
  'MSC Classic Markers' = "Itgb1",
  'MSC Classic Markers' = "Mcam",
  'MSC Classic Markers' = "Nt5e",
  'MSC Classic Markers' = "Pdgfra",
  'MSC Classic Markers' = "Sdc2",
  'MSC Classic Markers' = "Vcam1",
  
  # MSC 分泌细胞因子
  'MSC Secreted Cytokines' = "Angpt1",
  'MSC Secreted Cytokines' = "Angpt2",
  'MSC Secreted Cytokines' = "Bmp2",
  'MSC Secreted Cytokines' = "Bmp4",
  # 'MSC Secreted Cytokines' = "Ccl11",
  # 'MSC Secreted Cytokines' = "Ccl2",
  # 'MSC Secreted Cytokines' = "Ccl7",
  # 'MSC Secreted Cytokines' = "Cxcl1",
  'MSC Secreted Cytokines' = "Cxcl10",
  'MSC Secreted Cytokines' = "Cxcl12",
  'MSC Secreted Cytokines' = "Cxcl9",
  'MSC Secreted Cytokines' = "Fgf2",
  #'MSC Secreted Cytokines' = "Gdf10",
  'MSC Secreted Cytokines' = "Gdnf",
  'MSC Secreted Cytokines' = "Hgf",
  'MSC Secreted Cytokines' = "Ido1",
  'MSC Secreted Cytokines' = "Igf1",
  'MSC Secreted Cytokines' = "Il10",
  'MSC Secreted Cytokines' = "Il1rn",
  'MSC Secreted Cytokines' = "Lif",
  'MSC Secreted Cytokines' = "Pdgfa",
  'MSC Secreted Cytokines' = "Pdgfb",
  'MSC Secreted Cytokines' = "Ptgs2",
  #'MSC Secreted Cytokines' = "Serpine1",
  'MSC Secreted Cytokines' = "Tgfb1",
  'MSC Secreted Cytokines' = "Tgfb2",
  'MSC Secreted Cytokines' = "Tgfb3",
  'MSC Secreted Cytokines' = "Timp1",
  'MSC Secreted Cytokines' = "Vegfa",
  #'MSC Secreted Cytokines' = "Vegfd",
  
  # MSC 中实际表达的细胞因子
  'MSC Expressed Cytokines' = "Ccl11",
  'MSC Expressed Cytokines' = "Ccl2",
  'MSC Expressed Cytokines' = "Ccl7",
  'MSC Expressed Cytokines' = "Cxcl1",
  'MSC Expressed Cytokines' = "Gdf10",
  'MSC Expressed Cytokines' = "Serpine1",
  'MSC Expressed Cytokines' = "Vegfd"
)



Idents(Fibroblast3)<-'RNA_snn_res.0.4'
Fibroblast3_markers <-
  FindAllMarkers(
    JoinLayers(Fibroblast3),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

Fibroblast3_top10_markers <- 
  Fibroblast3_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10) %>%
  nest() %>%
  deframe() %>% 
  map(pull, 'gene') %>%
  imap(function(genes, name){
    re <- genes %>%
      map(function(gene){
        # Fibroblastrowser()
        names(gene) = str_glue("cluster_{name}")
        gene
      })
    re
  }) %>%
  unlist() %>% 
  .[!duplicated(.)]

names(Fibroblast3_top10_markers) <- sub('.*\\.','',names(Fibroblast3_top10_markers))

Fibroblast3_top10_markers1 <-
  Fibroblast3_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10)

# dotplot
pdf(file.path(Fibroblast3_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(Fibroblast3, group.by = "RNA_snn_res.0.4", features = intestinal_marker_ls)+RotatedAxis()
DotPlot(Fibroblast3, group.by = "RNA_snn_res.0.4", features = msc_marker_ls)+RotatedAxis()
DotPlot(Fibroblast3, group.by = "RNA_snn_res.0.4", features = Fibroblast3_top10_markers)+RotatedAxis()
dev.off()

# Idents(Fibroblast3) <- "RNA_snn_res.0.6"
# Fibroblast3 <- FindSubCluster(Fibroblast3, cluster = 1, graph.name = "RNA_snn", algorithm = 1,resolution = 0.3)
# Idents(Fibroblast3) <- "sub.cluster"
# DimPlot(Fibroblast3, label = TRUE, reduction = 'harmonyumap') + NoLegend()
# # DimPlot(Fibroblast3, label = TRUE, reduction = 'harmonyumap', split.by = "sub.cluster") + NoLegend()
# # DimPlot(Fibroblast3, label = TRUE, reduction = 'harmonyumap', split.by = "RNA_snn_res.0.6") + NoLegend()
# DotPlot(Fibroblast3, group.by = "sub.cluster", features = Fibro_major_marker_ls2)+RotatedAxis()

###score--------------------
MSC_Classic <- c("Cd248","Cd34","Cd44","Eng","Itgb1","Mcam","Nt5e","Pdgfra","Sdc2","Vcam1")
MSC_Secreted <- c("Angpt1","Angpt2","Bmp2","Bmp4","Cxcl10","Cxcl12","Cxcl9","Fgf2","Gdnf","Hgf","Ido1","Igf1","Il10","Il1rn","Lif","Pdgfa","Pdgfb","Ptgs2","Tgfb1","Tgfb2","Tgfb3","Timp1","Vegfa")
MSC_Expressed <- c("Ccl11","Ccl2","Ccl7","Cxcl1","Gdf10","Serpine1","Vegfd")

pathway<-list(MSC_Classic, MSC_Secreted, MSC_Expressed)
names(pathway)<-c('MSC_Classic','MSC_Secreted','MSC_Expressed')

# M1_M2<-read.table("/home/lr/ATAC_thirteen/Immune/MG/M1_M2.txt",header=TRUE,fill=TRUE)
n_na_idx <- pathway %>% map_lgl(function(x) all(is.na(x))) %>% {! .}%>%  which
pathway <- pathway[n_na_idx] %>% map(function(x) x[!is.na(x)])

for ( i in seq_along(names(pathway))){
  Fibroblast3 <- AddModuleScore(JoinLayers(Fibroblast3), pathway[i], name =  names(pathway)[i],nbin=5,ctrl=100)
}

features_to_plot <- paste0(c(
  'MSC_Classic','MSC_Secreted','MSC_Expressed'), '1')

p1 <- VlnPlot(Fibroblast3, features = features_to_plot[1], group.by = "RNA_snn_res.0.3", pt.size = 0)
p2 <- VlnPlot(Fibroblast3, features = features_to_plot[2], group.by = "RNA_snn_res.0.3", pt.size = 0)
p3 <- VlnPlot(Fibroblast3, features = features_to_plot[3], group.by = "RNA_snn_res.0.3", pt.size = 0)

p <- (p1 | p2 | p3) & theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p,file = file.path(Fibroblast3_out_dir,"MSC_score_cluster_Vlnplot.pdf"),width = 12, height = 5)


p1 <- VlnPlot(Fibroblast3, features = features_to_plot[1], group.by = "group", pt.size = 0)
p2 <- VlnPlot(Fibroblast3, features = features_to_plot[2], group.by = "group", pt.size = 0)
p3 <- VlnPlot(Fibroblast3, features = features_to_plot[3], group.by = "group", pt.size = 0)

p <- (p1 | p2 | p3) & theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p,file = file.path(Fibroblast3_out_dir,"MSC_score_group_Vlnplot.pdf"),width = 12, height = 5)

p1 <- DotPlot(Fibroblast3, features = features_to_plot[1], group.by = "group")
p2 <- DotPlot(Fibroblast3, features = features_to_plot[2], group.by = "group")
p3 <- DotPlot(Fibroblast3, features = features_to_plot[3], group.by = "group")

p <- (p1 | p2 | p3) & theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p,file = file.path(Fibroblast3_out_dir,"MSC_score_group_Dotplot.pdf"),width = 12, height = 5)




DimPlot(Fibroblast3, group.by = "celltype")  # 先看聚类

FeaturePlot(Fibroblast3, features = features_to_plot, ncol = 3, order = TRUE)


library(dplyr)

avg_scores <- Fibroblast3@meta.data %>%
  group_by(group) %>%
  summarise(
    MSC_Classic = mean(MSC_Classic1, na.rm = TRUE),
    MSC_Secreted = mean(MSC_Secreted1, na.rm = TRUE),
    MSC_Expressed = mean(MSC_Expressed1, na.rm = TRUE)
  ) %>%
  column_to_rownames("group")

pheatmap::pheatmap(avg_scores, 
                   scale = "row",          # 按行（亚群）标准化
                   cluster_rows = TRUE,
                   cluster_cols = FALSE,
                   display_numbers = TRUE)

ggplot(Fibroblast3@meta.data, aes(x = MSC_Classic1, y = MSC_Secreted1)) +
  geom_point(aes(color = celltype), alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_minimal()





## annotation ----
cell_cluster <- c(
  '0' = "iCAF", #CXCL12 
  '1' = "mCAF", #POSTN COL10A1 COL11A1 CTHRC1 MMP11
  '2'  = "mCAF", #POSTN COL10A1 COL11A1 CTHRC1 MMP11
  '3'  = "mCAF", #POSTN CTHRC1 MMMP11
  '4'  = "Epi_like CAF", #肿瘤上皮细胞
  '5'  = "iCAF", # CXCL9 CXCL10 CXCL11 
  '6'  = "Fib", #成纤维CCDC80 SFRP4
  '7'  = "Mac_like CAF", # 巨噬细胞
  '8'  = "mCAF", #POSTN CTHRC1 MMMP11 
  '9'  = "Plasma_like CAF",#浆细胞
  '10' = "meCAF" # BNIP3 SLC2A1 ADM ANGPTL4
)

Idents(Fibroblast3) <- "RNA_snn_res.0.6"
#Idents(Fibroblast3) <- "sub.cluster"

Fibroblast3 <- RenameIdents(Fibroblast3, cell_cluster)
Fibroblast3$subtype <- Idents(Fibroblast3)

Fib_marker_ls <- c(
  "iCAF" = 'CXCL12',
  "iCAF" = 'CXCL9',
  "iCAF" = 'CXCL10',
  "iCAF" = 'CXCL11',
  
  "mCAF" = 'POSTN',
  "mCAF" = 'COL10A1',
  "mCAF" = 'COL11A1',
  "mCAF" = 'CTHRC1',
  "mCAF" = 'MMP11',
  
  "Epi_like CAF" = 'EPCAM',
  "Epi_like CAF" = 'CDH1',
  "Epi_like CAF" = 'STMN1',
  
  "Fib" = 'CCDC80',
  "Fib" = 'SFRP4',
  
  "Mac_like CAF" = 'CD68',
  
  "Plasma_like CAF" = 'MZB1',
  "Plasma_like CAF" = 'XBP1',
  
  "meCAF" = "BNIP3" ,
  "meCAF" ='SLC2A1' ,
  "meCAF" ='ADM',
  "meCAF" ='ANGPTL4'
)

pdf(file.path(Fibroblast3_out_dir, "dotplot celltype marker.pdf"), width = 20, height = 8)
DotPlot(Fibroblast3, group.by = "subtype", features = Fib_marker_ls, dot.scale = 10) + theme_90angle
dev.off()

pdf(file.path(Fibroblast3_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
# DimPlot(Fibroblast3, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
DimPlot(Fibroblast3, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Fibroblast3_out_dir, str_glue("Dimplot dim{dims} annotation split3.pdf")), height = 6, width = 14)
DimPlot(Fibroblast3, group.by = "subtype",split.by = 'group3', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

saveRDS(Fibroblast3, file.path(Fibroblast3_out_dir,"Fibroblast3.rds"))

## cell props barplot  --------
Fibroblast3_meta <- Fibroblast3[[]]

Fibroblast3_meta$sample <- factor(Fibroblast3_meta$sample, 
                                 levels=c(Fibroblast3_meta$sample %>% unique())
)
p_col_proj <-
  Fibroblast3_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(Fibroblast3_out_dir,"cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  Fibroblast3_meta %>%
  group_by(group3, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group3, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path(Fibroblast3_out_dir,"cell_freq_group3.pdf"),
  height = 6,
  width  = 4
)


## 映射-------------------------------
library(tibble)
library(data.table)
sub <- all_CAFs[['subtype']] %>% rownames_to_column('cells')

meta.data.t <-
  sub%>%
  right_join(
    Fibroblast3[['subtype']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    subtype2 = ifelse((as.character(subtype.x) %in% c('myCAFs', 'iCAFs', 'CAFs', 'CAFs2', 'vCAFs', 'tCAFs', 'rCAFs', 'dCAFs')),
                      as.character(subtype.x), 
                      as.character(subtype.y)))%>%
  dplyr::select(cells, subtype2) %>%
  column_to_rownames('cells')


Fibroblast3 <- AddMetaData(Fibroblast3,meta.data.t['subtype2'])
Fibroblast3$subtype_raw <- Fibroblast3$subtype
Fibroblast3$subtype <- Fibroblast3$subtype2
saveRDS(Fibroblast3, file.path(Fibroblast3_out_dir,'all_Fibroblast.rds'))



#3、BM-MSC再分群、打分------
BMMSC_out_dir<-'result/BMMSC_subtype/' 
dir.create(BMMSC_out_dir)
BMMSC<-subset(all_data2, subset = annotation %in% c("BM_MSCs"))
BMMSC$annotation %>% table()

## PCA --------
DefaultAssay(BMMSC) <- "RNA"
BMMSC <- FindVariableFeatures(BMMSC)
BMMSC <- NormalizeData(BMMSC)
BMMSC <- ScaleData(BMMSC, vars.to.regress = 'nCount_RNA')
BMMSC <- RunPCA(BMMSC)


### pca inspector ------
pdf(file.path(BMMSC_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(BMMSC, ndims = 50)
dev.off()

# pdf(file.path(BMMSC_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
# DimHeatmap(BMMSC, dims = 1:5)
# dev.off()
# pdf(file.path(BMMSC_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
# DimHeatmap(BMMSC, dims = 6:10)
# dev.off()
# pdf(file.path(BMMSC_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
# DimHeatmap(BMMSC, dims = 11:15)
# dev.off()
# pdf(file.path(BMMSC_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
# DimHeatmap(BMMSC, dims = 16:20)
# dev.off()


## dimension reduction --------
dims <- 20

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb

BMMSC <- RunUMAP(BMMSC, dims = 1:dims,return.model = TRUE)

pdf(file.path(BMMSC_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(BMMSC, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

table(all_data2@meta.data$celltype, all_data2@meta.data$project)

BMMSC <- BMMSC %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

BMMSC <- RunUMAP(
  BMMSC, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(BMMSC_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(BMMSC, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'

### cluster -------
BMMSC  <-
  FindNeighbors(
    BMMSC ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )

BMMSC <-
  FindClusters(
    BMMSC,
    resolution = c(0.2, 0.3, 0.4)
  )

# dimplot
pdf(file.path(BMMSC_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(BMMSC, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(BMMSC, group.by = "RNA_snn_res.0.3", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(BMMSC, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


pdf(file.path(BMMSC_out_dir, str_glue("FeaturePlot dim{dims} Pdgfra.pdf")), height = 6, width = 8)
FeaturePlot(BMMSC,feature = c("Pdgfra"),reduction = 'harmonyumap')
dev.off()


BMMSC$res_1.2 %>% table()
pdf(file.path(BMMSC_out_dir, str_glue("Dimplot dim{dims} res_1.2.pdf")), height = 6, width = 8)

DimPlot(BMMSC, group.by = "res_1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()

dev.off()

pdf(file.path(BMMSC_out_dir, str_glue("Dimplot dim{dims} celltype.pdf")), height = 6, width = 8)

DimPlot(BMMSC, group.by = "celltype", label = TRUE, reduction = 'harmonyumap') + NoLegend()

dev.off()




### dotplot -------
intestinal_marker_ls <- c(
  # 上皮细胞群
  'Epithelial Stem Cell' = "Lgr5",
  'Epithelial Stem Cell' = "Ascl2",
  'Epithelial Stem Cell' = "Sox9",
  
  "Epcam",
  
  'Paneth Cell' = "Mmp7",
  
  'Enterocyte' = "Slc26a3",
  'Enterocyte' = "Vil1",
  
  'Goblet Cell' = "Muc2",
  'Goblet Cell' = "Spdef",
  
  'Enteroendocrine Cell (EEC)' = "Neurod1",
  'Enteroendocrine Cell (EEC)' = "Tph1",
  'Enteroendocrine Cell (EEC)' = "Gcg",
  'Enteroendocrine Cell (EEC)' = "Sst",
  
  'Tuft Cell' = "Il25",
  'Tuft Cell' = "Trpm5",
  'Tuft Cell' = "Avil",
  
  'M Cell' = "Tnfrsf11a",
  
  # 免疫细胞群
  'CD4+ T Cell' = "Cd3d",
  'CD4+ T Cell' = "Cd3e",
  'CD4+ T Cell' = "Cd3g",
  'CD4+ T Cell' = "Cd4",
  'CD4+ T Cell' = "Foxp3",
  
  # 'CD8+ T Cell' = "Cd3d",
  # 'CD8+ T Cell' = "Cd3e",
  # 'CD8+ T Cell' = "Cd3g",
  'CD8+ T Cell' = "Cd8a",
  'CD8+ T Cell' = "Gzmb",
  
  'B Cell' = "Cd19",
  'B Cell' = "Cd79a",
  'B Cell' = "Ighm",
  
  'Plasma Cell' = "Sdc1",
  'Plasma Cell' = "Jchain",
  
  'Dendritic Cell (DC)' = "Itgax",
  'Dendritic Cell (DC)' = "Zbtb46",
  'Dendritic Cell (DC)' = "Clec9a",
  
  'Macrophage' = "Adgre1",
  'Macrophage' = "Csf1r",
  'Macrophage' = "Cd68",
  'Macrophage' = "Cx3cr1",
  
  'Innate Lymphoid Cell (ILC)' = "Il7r",
  'Innate Lymphoid Cell (ILC)' = "Id2",
  'Innate Lymphoid Cell (ILC)' = "Rorc",
  'Innate Lymphoid Cell (ILC)' = "Gata3",
  
  # 基质细胞群
  'Fibroblast' = "Col1a1",
  'Fibroblast' = "Pdgfra",
  'Fibroblast' = "Pdgfrb",
  
  'Endothelial Cell' = "Pecam1",
  'Endothelial Cell' = "Cdh5",
  'Endothelial Cell' = "Lyve1",
  
  'Pericyte' = "Cspg4",
  #'Pericyte' = "Pdgfrb",
  
  'Enteric Glial Cell (EGC)' = "Gfap",
  'Enteric Glial Cell (EGC)' = "S100b"
)

msc_marker_ls <- c(
  # MSC 经典标记基因
  'MSC Classic Markers' = "Cd248",
  'MSC Classic Markers' = "Cd34",
  'MSC Classic Markers' = "Cd44",
  'MSC Classic Markers' = "Eng",
  'MSC Classic Markers' = "Itgb1",
  'MSC Classic Markers' = "Mcam",
  'MSC Classic Markers' = "Nt5e",
  'MSC Classic Markers' = "Pdgfra",
  'MSC Classic Markers' = "Sdc2",
  'MSC Classic Markers' = "Vcam1",
  
  # MSC 分泌细胞因子
  'MSC Secreted Cytokines' = "Angpt1",
  'MSC Secreted Cytokines' = "Angpt2",
  'MSC Secreted Cytokines' = "Bmp2",
  'MSC Secreted Cytokines' = "Bmp4",
  # 'MSC Secreted Cytokines' = "Ccl11",
  # 'MSC Secreted Cytokines' = "Ccl2",
  # 'MSC Secreted Cytokines' = "Ccl7",
  # 'MSC Secreted Cytokines' = "Cxcl1",
  'MSC Secreted Cytokines' = "Cxcl10",
  'MSC Secreted Cytokines' = "Cxcl12",
  'MSC Secreted Cytokines' = "Cxcl9",
  'MSC Secreted Cytokines' = "Fgf2",
  #'MSC Secreted Cytokines' = "Gdf10",
  'MSC Secreted Cytokines' = "Gdnf",
  'MSC Secreted Cytokines' = "Hgf",
  'MSC Secreted Cytokines' = "Ido1",
  'MSC Secreted Cytokines' = "Igf1",
  'MSC Secreted Cytokines' = "Il10",
  'MSC Secreted Cytokines' = "Il1rn",
  'MSC Secreted Cytokines' = "Lif",
  'MSC Secreted Cytokines' = "Pdgfa",
  'MSC Secreted Cytokines' = "Pdgfb",
  'MSC Secreted Cytokines' = "Ptgs2",
  #'MSC Secreted Cytokines' = "Serpine1",
  'MSC Secreted Cytokines' = "Tgfb1",
  'MSC Secreted Cytokines' = "Tgfb2",
  'MSC Secreted Cytokines' = "Tgfb3",
  'MSC Secreted Cytokines' = "Timp1",
  'MSC Secreted Cytokines' = "Vegfa",
  #'MSC Secreted Cytokines' = "Vegfd",
  
  # MSC 中实际表达的细胞因子
  'MSC Expressed Cytokines' = "Ccl11",
  'MSC Expressed Cytokines' = "Ccl2",
  'MSC Expressed Cytokines' = "Ccl7",
  'MSC Expressed Cytokines' = "Cxcl1",
  'MSC Expressed Cytokines' = "Gdf10",
  'MSC Expressed Cytokines' = "Serpine1",
  'MSC Expressed Cytokines' = "Vegfd"
)



Idents(BMMSC)<-'RNA_snn_res.0.2'
BMMSC_markers <-
  FindAllMarkers(
    JoinLayers(BMMSC),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

BMMSC_top10_markers <- 
  BMMSC_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10) %>%
  nest() %>%
  deframe() %>% 
  map(pull, 'gene') %>%
  imap(function(genes, name){
    re <- genes %>%
      map(function(gene){
        # BMMSCrowser()
        names(gene) = str_glue("cluster_{name}")
        gene
      })
    re
  }) %>%
  unlist() %>% 
  .[!duplicated(.)]

names(BMMSC_top10_markers) <- sub('.*\\.','',names(BMMSC_top10_markers))

BMMSC_top10_markers1 <-
  BMMSC_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10)

# dotplot
pdf(file.path(BMMSC_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(BMMSC, group.by = "RNA_snn_res.0.2", features = intestinal_marker_ls)+RotatedAxis()
DotPlot(BMMSC, group.by = "RNA_snn_res.0.2", features = msc_marker_ls)+RotatedAxis()
DotPlot(BMMSC, group.by = "RNA_snn_res.0.2", features = BMMSC_top10_markers)+RotatedAxis()
dev.off()

# Idents(BMMSC) <- "RNA_snn_res.0.6"
# BMMSC <- FindSubCluster(BMMSC, cluster = 1, graph.name = "RNA_snn", algorithm = 1,resolution = 0.3)
# Idents(BMMSC) <- "sub.cluster"
# DimPlot(BMMSC, label = TRUE, reduction = 'harmonyumap') + NoLegend()
# # DimPlot(BMMSC, label = TRUE, reduction = 'harmonyumap', split.by = "sub.cluster") + NoLegend()
# # DimPlot(BMMSC, label = TRUE, reduction = 'harmonyumap', split.by = "RNA_snn_res.0.6") + NoLegend()
# DotPlot(BMMSC, group.by = "sub.cluster", features = Fibro_major_marker_ls2)+RotatedAxis()

###AddModuleScore--------------------
MSC_Classic <- c("Cd248","Cd34","Cd44","Eng","Itgb1","Mcam","Nt5e","Pdgfra","Sdc2","Vcam1")
MSC_Secreted <- c("Angpt1","Angpt2","Bmp2","Bmp4","Cxcl10","Cxcl12","Cxcl9","Fgf2","Gdnf","Hgf","Ido1","Igf1","Il10","Il1rn","Lif","Pdgfa","Pdgfb","Ptgs2","Tgfb1","Tgfb2","Tgfb3","Timp1","Vegfa")
MSC_Expressed <- c("Ccl11","Ccl2","Ccl7","Cxcl1","Gdf10","Serpine1","Vegfd")

pathway<-list(MSC_Classic, MSC_Secreted, MSC_Expressed)
names(pathway)<-c('MSC_Classic','MSC_Secreted','MSC_Expressed')

# M1_M2<-read.table("/home/lr/ATAC_thirteen/Immune/MG/M1_M2.txt",header=TRUE,fill=TRUE)
n_na_idx <- pathway %>% map_lgl(function(x) all(is.na(x))) %>% {! .}%>%  which
pathway <- pathway[n_na_idx] %>% map(function(x) x[!is.na(x)])

for ( i in seq_along(names(pathway))){
  BMMSC <- AddModuleScore(JoinLayers(BMMSC), pathway[i], name =  names(pathway)[i],nbin=5,ctrl=100)
}

features_to_plot <- paste0(c(
  'MSC_Classic','MSC_Secreted','MSC_Expressed'), '1')

p1 <- VlnPlot(BMMSC, features = features_to_plot[1], group.by = "RNA_snn_res.0.2", pt.size = 0)
p2 <- VlnPlot(BMMSC, features = features_to_plot[2], group.by = "RNA_snn_res.0.2", pt.size = 0)
p3 <- VlnPlot(BMMSC, features = features_to_plot[3], group.by = "RNA_snn_res.0.2", pt.size = 0)

p <- (p1 | p2 | p3) & theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p,file = file.path(BMMSC_out_dir,"MSC_score_cluster_Vlnplot.pdf"),width = 12, height = 5)


p1 <- VlnPlot(BMMSC, features = features_to_plot[1], group.by = "group", pt.size = 0)
p2 <- VlnPlot(BMMSC, features = features_to_plot[2], group.by = "group", pt.size = 0)
p3 <- VlnPlot(BMMSC, features = features_to_plot[3], group.by = "group", pt.size = 0)

p <- (p1 | p2 | p3) & theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p,file = file.path(BMMSC_out_dir,"MSC_score_group_Vlnplot.pdf"),width = 12, height = 5)

p1 <- DotPlot(BMMSC, features = features_to_plot[1], group.by = "group")
p2 <- DotPlot(BMMSC, features = features_to_plot[2], group.by = "group")
p3 <- DotPlot(BMMSC, features = features_to_plot[3], group.by = "group")

p <- (p1 | p2 | p3) & theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p,file = file.path(BMMSC_out_dir,"MSC_score_group_Dotplot.pdf"),width = 12, height = 5)


DimPlot(BMMSC, group.by = "RNA_snn_res.0.2")  # 先看聚类

FeaturePlot(BMMSC, features = features_to_plot, ncol = 3, order = TRUE)


library(dplyr)

avg_scores <- BMMSC@meta.data %>%
  group_by(group) %>%
  summarise(
    MSC_Classic = mean(MSC_Classic1, na.rm = TRUE),
    MSC_Secreted = mean(MSC_Secreted1, na.rm = TRUE),
    MSC_Expressed = mean(MSC_Expressed1, na.rm = TRUE)
  ) %>%
  column_to_rownames("group")

pheatmap::pheatmap(avg_scores, 
                   scale = "row",          # 按行（亚群）标准化
                   cluster_rows = TRUE,
                   cluster_cols = FALSE,
                   display_numbers = TRUE)

ggplot(BMMSC@meta.data, aes(x = MSC_Classic1, y = MSC_Secreted1)) +
  geom_point(aes(color = RNA_snn_res.0.2), alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_minimal()


###GSVA-------
# library(GSVA)
# library(Seurat)
# 
# MSC_Classic <- c("Cd248","Cd34","Cd44","Eng","Itgb1","Mcam","Nt5e","Pdgfra","Sdc2","Vcam1")
# MSC_Secreted <- c("Angpt1","Angpt2","Bmp2","Bmp4","Cxcl10","Cxcl12","Cxcl9","Fgf2","Gdnf","Hgf","Ido1","Igf1","Il10","Il1rn","Lif","Pdgfa","Pdgfb","Ptgs2","Tgfb1","Tgfb2","Tgfb3","Timp1","Vegfa")
# MSC_Expressed <- c("Ccl11","Ccl2","Ccl7","Cxcl1","Gdf10","Serpine1","Vegfd")
# 
# pathway <- list(
#   MSC_Classic = MSC_Classic,
#   MSC_Secreted = MSC_Secreted,
#   MSC_Expressed = MSC_Expressed
# )
# 
# #expr_mat <- GetAssayData(BMMSC, slot = "data") 
# counts_mat <- GetAssayData(BMMSC, slot = "counts")
# 
# all_genes <- unique(unlist(pathway))
# common_genes <- intersect(all_genes, rownames(counts_mat))
# counts_sub <- counts_mat[common_genes, ]
# 
# # 运行 GSVA
# gsva_scores <- gsva(
#   expr = counts_sub,
#   gset.idx.list = pathway,
#   method = "gsva",          # or "ssgsea"
#   kcdf = "Poisson",         # for count data
#   mx.diff = TRUE,
#   verbose = FALSE
# )
# 
# gsva_df <- as.data.frame(t(gsva_scores))
# colnames(gsva_df) <- paste0(colnames(gsva_df), "_gsva")  # 避免重名
# 
# BMMSC@meta.data <- cbind(BMMSC@meta.data, gsva_df)
# 
# # 可视化
# features_gsva <- paste0(names(pathway), "_gsva")
# p1 <- VlnPlot(BMMSC, features = features_gsva[1], group.by = "RNA_snn_res.0.2", pt.size = 0)
# p2 <- VlnPlot(BMMSC, features = features_gsva[2], group.by = "RNA_snn_res.0.2", pt.size = 0)
# p3 <- VlnPlot(BMMSC, features = features_gsva[3], group.by = "RNA_snn_res.0.2", pt.size = 0)
# 
# (p1 | p2 | p3) & theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# 
# p1 <- VlnPlot(BMMSC, features = features_gsva[1], group.by = "group", pt.size = 0.1)
# p2 <- VlnPlot(BMMSC, features = features_gsva[2], group.by = "group", pt.size = 0.1)
# p3 <- VlnPlot(BMMSC, features = features_gsva[3], group.by = "group", pt.size = 0.1)
# 
# (p1 | p2 | p3) & theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# p1 <- DotPlot(BMMSC, features = features_gsva[1], group.by = "group")
# p2 <- DotPlot(BMMSC, features = features_gsva[2], group.by = "group")
# p3 <- DotPlot(BMMSC, features = features_gsva[3], group.by = "group")
# 
# (p1 | p2 | p3) & theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# 
# p1 <-FeaturePlot(BMMSC[,BMMSC$group == "CTRL"], features = features_gsva, ncol = 3, order = TRUE)
# p2 <-FeaturePlot(BMMSC[,BMMSC$group == "CtrlMSC"], features = features_gsva, ncol = 3, order = TRUE)
# p3 <-FeaturePlot(BMMSC[,BMMSC$group == "IR"], features = features_gsva, ncol = 3, order = TRUE)
# p4 <-FeaturePlot(BMMSC[,BMMSC$group == "IRMSC"], features = features_gsva, ncol = 3, order = TRUE)
# 
# (p1 | p2 | p3 | p4) & theme(axis.text.x = element_text(angle = 45, hjust = 1))


saveRDS(BMMSC, file.path(BMMSC_out_dir,'BMMSC.rds'))



#4、髓系亚类注释------
all_data2$annotation %>% table()
#all_data$annotation <- str_replace(all_data$annotation, 'Monocytes', 'Macrophages')

Myeloids_out_dir <- "result/Myeloids_subtype" 
dir.create(Myeloids_out_dir)
Myeloids <- subset(all_data2, annotation == 'Myeloid')

## PCA --------
DefaultAssay(Myeloids) <- "RNA"
Myeloids <- FindVariableFeatures(Myeloids)
Myeloids <- NormalizeData(Myeloids)
Myeloids <- ScaleData(Myeloids, vars.to.regress = 'nCount_RNA')
Myeloids <- RunPCA(Myeloids)



### pca inspector ------
pdf(file.path(Myeloids_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Myeloids, ndims = 50)
dev.off()

pdf(file.path(Myeloids_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(Myeloids, dims = 1:5)
dev.off()
pdf(file.path(Myeloids_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(Myeloids, dims = 6:10)
dev.off()
pdf(file.path(Myeloids_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(Myeloids, dims = 11:15)
dev.off()
pdf(file.path(Myeloids_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(Myeloids, dims = 16:20)
dev.off()

## dimension reduction --------
dims <- 20 # 18

library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

Myeloids <- RunUMAP(Myeloids, dims = 1:dims,return.model = TRUE)

pdf(file.path(Myeloids_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Myeloids, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

Myeloids <- Myeloids %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Myeloids <- RunUMAP(
  Myeloids, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(Myeloids_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Myeloids, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Myeloids  <-
  FindNeighbors(
    Myeloids ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Myeloids <-
  FindClusters(
    Myeloids,
    resolution = c(0.2, 0.4)
  )

# dimplot
pdf(file.path(Myeloids_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Myeloids, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Myeloids, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Myeloids, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Myeloids, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

Idents(Myeloids) <- "RNA_snn_res.0.4"


### dotplot -------
Myeloids_marker_ls <- c(
  M1 = 'Nos2',
  M1 = 'Il12a', # IL12 (p35) 在小鼠中通常指 Il12a
  M1 = 'Fcgr1',
  M1 = 'Cd80',
  M1 = 'Il23a', # IL23 (p19) 在小鼠中通常指 Il23a
  M1 = 'Cxcl9',
  M1 = 'Cxcl10',
  M1 = 'Cxcl11',
  M1 = 'Cd86',
  M1 = 'Il1a',
  M1 = 'Il1b',
  M1 = 'Il6',
  M1 = 'Tnf',
  M1 = 'H2-Aa', # MHCII 在小鼠中主要指 H2 类分子，如 H2-Aa 或 H2-Ab1
  M1 = 'Ccl5',
  M1 = 'Irf5',
  M1 = 'Irf1',
  M1 = 'Cd40',
  M1 = 'Kynu',
  M1 = 'Ccr7',
  M1 = 'Vcan',
  M1 = 'Cd300e',
  M1 = 'Inhba',
  M1 = 'Tnip3',
  M2 = 'Arg1',
  M2 = 'Il10',
  M2 = 'Fcgr2b', # CD32 在小鼠中主要对应 Fcgr2b
  M2 = 'Cd163',
  M2 = 'Fcer2a', # FCER2 在小鼠中为 Fcer2a
  M2 = 'Cd200r1',
  M2 = 'Pdcd1lg2',
  M2 = 'Cd274',
  M2 = 'Marco',
  M2 = 'Csf1r',
  M2 = 'Mrc1',
  M2 = 'Il1rn',
  M2 = 'Il1r2',
  M2 = 'Il4ra', # IL4R 在小鼠中为 Il4ra
  M2 = 'Ccl4',
  M2 = 'Ccl13', # 小鼠无直接 Ccl13 同源物，功能上可能由 Ccl17/22 等替代，此处保留原名
  M2 = 'Ccl20',
  M2 = 'Ccl17',
  M2 = 'Ccl18', # 小鼠无直接 Ccl18 同源物，此处保留原名
  M2 = 'Ccl22',
  M2 = 'Ccl24',
  M2 = 'Lyve1',
  M2 = 'Vegfa',
  M2 = 'Vegfb',
  M2 = 'Vegfc',
  M2 = 'Vegfd',
  M2 = 'Egf',
  M2 = 'Ctsa',
  M2 = 'Ctsb',
  M2 = 'Cst3', # CSTC 应为 Cystatin C，在小鼠中为 Cst3
  M2 = 'Ctsd',
  M2 = 'Tgfb1',
  M2 = 'Cd209a', # CD209 在小鼠中有多个同源物，如 Cd209a
  M2 = 'Cd163l1',
  M2 = 'Folr2',
  M0 = 'Cd68',
  M0 = 'C1qa',
  M0 = 'C1qb',
  M0 = 'Selenop',
  M0 = 'Lyz1', # LYZ 在小鼠中主要为 Lyz1
  M0 = 'H2-Dpb1', # HLA-DPB1 在小鼠中对应 H2-Dpb1
  M0 = 'H2-Dpa1', # HLA-DPA1 在小鼠中对应 H2-Dpa1
  M0 = 'Aif1', 
  TAM = 'Spp1', 
  TAM = 'C1qc' # CQC 应为 C1QC 的笔误，其小鼠同源基因为 C1qc
)


Myeloids_major_marker_ls <- c(
  # 单核细胞亚群
  'Classical Monocyte' = "Cd14",
  'Classical Monocyte' = "S100a8",
  'Non-classical Monocyte' = "Fcgr3",  # CD16
  
  # 巨噬细胞亚群
  'M1-lMacrophage' = "Cd86",
  'M1-Macrophage' = "Cxcl9",
  'M2-Macrophage' = "Cd163",
  'M2-Macrophage' = "Msr1",
  'Pro-angiogenic TAM' = "Spp1",       # 促血管生成型[5](@ref)
  'Anti-inflammatory TAM' = "C1qc",    # 抗炎型[5](@ref)
  'Anti-inflammatory TAM' = "Apoe",
  
  # 树突状细胞亚群
  "Dendritic" = "Fcer1a",
  "Dendritic" = "Jchain",
  'cDC1' = "Clec9a",                   # 抗原呈递型[2,5](@ref)
  'cDC1' = "Xcr1",
  'cDC2' = "Cd1c",
  'LAMP3+ DC' = "Lamp3",               # 迁移型[5](@ref)
  
  # 粒细胞亚群
  'Neutrophil' = "Csf3r",
  'Neutrophil' = "Cxcr2",
  'Eosinophil' = "Siglec8", # 小鼠 Siglec8 功能与人不同，有时用 Siglecf 替代，此处保留
  'Basophil' = "Il3ra", # CD123 即 IL3RA
  
  # 肥大细胞
  'Mast' = "Kit",                 # 干细胞因子受体[2](@ref)
  'Mast' = "Tpsab1",              # 类胰蛋白酶[2](@ref)
  "Mast" = "Ms4a2",
  "Mast" = "Cpa3",
  
  # 髓系祖细胞
  'CMP' = "Cd34",                      # 共同髓系祖细胞[1](@ref)
  'GMP' = "Cd33"                       # 粒-单核祖细胞[4](@ref)
)

myeloids_cell_markers <- list(
  "Mast_KIT" = c("Kit", "Tpsab1", "Cpa3"),
  "pDC_LILRA4" = c("Lilra4", "Gzmb", "Il3ra"), # LILRA4 在小鼠中无直接同源物，此处按规则转换
  "cDC1_CLEC9A" = c("Clec9a", "Flt3", "Ido1"),
  "cDC2_CD1C" = c("Cd1c", "Fcer1a", "H2-Dqa1"), # HLA-DQA1 在小鼠中对应 H2-Dqa1
  "cDC3_LAMP3" = c("Lamp3", "Ccr7", "Fscn1"), # CCFR7 应为 CCR7 的笔误
  "Mono_CD14" = c("Fcn1", "S100a8", "S100a9"),
  "Mono_CD16" = c("Fcgr3", "Lst1", "Lilrb2"), # LILRB2 在小鼠中有多个同源物，如 Lilrb4
  "Macro_INHBA" = c("Inhba","Il1rn", "Ccl4"),
  "Macro_NLRP3" = c("Nlrp3", "Ereg", "Il1b"),
  "Macro_LYVE1" = c("Lyve1", "Pltp", "Sepp1"),
  "Macro_C1QC" = c("C1qc", "C1qa", "Apoe")
)


Idents(Myeloids)<-'RNA_snn_res.0.4'
Myeloids_markers <-
  FindAllMarkers(
    JoinLayers(Myeloids),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )

Myeloids_top10_markers <- 
  Myeloids_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10) %>%
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

names(Myeloids_top10_markers) <- sub('.*\\.','',names(Myeloids_top10_markers))


#提取deepseek
Myeloids_top10_markers1 <- 
  Myeloids_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10) 


# dotplot
pdf(file.path(Myeloids_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 30, height = 6)
DotPlot(Myeloids, group.by = "RNA_snn_res.0.4", features = intestinal_marker_ls)+RotatedAxis()
DotPlot(Myeloids, group.by = "RNA_snn_res.0.4", features = Myeloids_marker_ls)+RotatedAxis()
DotPlot(Myeloids, group.by = "RNA_snn_res.0.4", features = Myeloids_major_marker_ls)+RotatedAxis()
DotPlot(Myeloids, group.by = "RNA_snn_res.0.4", features = myeloids_cell_markers)+RotatedAxis()
DotPlot(Myeloids, group.by = "RNA_snn_res.0.4", features = Myeloids_top10_markers)+RotatedAxis()
dev.off()



###score--------------------
M1 <- c(
  'Nos2',
  'Il12a',      # IL12 (p35 subunit)
  'Fcgr1',      # FCGR1A/B/C all map to Fcgr1 in mouse
  'Cd80',
  'Il23a',      # IL23 (p19 subunit)
  'Cxcl9',
  'Cxcl10',
  'Cxcl11',
  'Cd86',
  'Il1a',
  'Il1b',
  'Il6',
  'Tnf',
  'H2-Ab1',     # MHC Class II (e.g., I-A in C57BL/6)
  'Ccl5',
  'Irf5',
  'Irf1',
  'Cd40',
  'Kynu',
  'Ccr7',
  'Vcan',
  'Cd300e',
  'Inhba',
  'Tnip3'
)

M2 <- c(
  'Arg1',
  'Il10',
  'Fcgr2b',     # CD32
  'Cd163',
  'Fcer2a',     # FCER2
  'Cd200r1',
  'Pdcd1lg2',
  'Cd274',
  'Marco',
  'Csf1r',
  'Mrc1',
  'Il1rn',
  'Il1r2',
  'Il4ra',      # IL4R
  'Ccl4',
  'Ccl13',      # Note: No direct ortholog, functional analogs exist
  'Ccl20',
  'Ccl17',
  'Ccl18',      # Note: No direct ortholog in mouse
  'Ccl22',
  'Ccl24',
  'Lyve1',
  'Vegfa',
  'Vegfb',
  'Vegfc',
  'Vegfd',
  'Egf',
  'Ctsa',
  'Ctsb',
  'Cst3',       # CSTC (Cystatin C)
  'Ctsd',
  'Tgfb1',
  'Cd209a',     # CD209 (multiple paralogs in mouse)
  'Cd163l1',
  'Folr2'
)

Neutrophils <- c("Fcgr4", "Ifitm2", "Csf3r") # FCGR3B maps to Fcgr4 in mouse

Plasma <- c("Tnfrsf17", "Jchain", "Mzb1", "Xbp1", "Ighg1", "Igha1", "Ighg3", "Igkc")

Tcells <- c("Cd3d", "Cd3e", "Cd3g", "Lck", "Trac", "Cd4", "Cd8a", "Gzmk")

Fibroblasts <- c("Sfrp4", "Pdgfra", "Ccdc80", "Col1a1", "Acta2", 'Dcn', "Fap", "Tnc", "Podoplanin", "Pdgfrb", "Col4a2", "Col5a2", "Stc2", "Loxl2", "Lum", "Vim", "Col1a2")

Bcells <- c("Cd19", "Cd79a", "Ms4a1")

DC <- c("Cd83", "Lilra4", "Clec9a", "Ido1", "Cd1c", "Fcer1a")

pathway<-list(M1,M2,Neutrophils,Plasma,Tcells,Fibroblasts,Bcells,DC)
names(pathway)<-c('M1','M2','Neutrophils','Plasma','Tcells','Fibroblasts','Bcells','DC')

# M1_M2<-read.table("/home/lr/ATAC_thirteen/Immune/MG/M1_M2.txt",header=TRUE,fill=TRUE)
n_na_idx <- pathway %>% map_lgl(function(x) all(is.na(x))) %>% {! .}%>%  which
pathway <- pathway[n_na_idx] %>% map(function(x) x[!is.na(x)])

for ( i in seq_along(names(pathway))){
  Myeloids <- AddModuleScore(JoinLayers(Myeloids), pathway[i], name =  names(pathway)[i],nbin=5,ctrl=10)
}

p1<-VlnPlot(Myeloids,
           features = c('M11','M21','Neutrophils1','Plasma1','Tcells1','Fibroblasts1','Bcells1','DC1'),
           group.by='RNA_snn_res.0.4',
           pt.size = 0)

p2 <- FeaturePlot(Myeloids,
            features = c('M11','M21','Neutrophils1','Plasma1','Tcells1','Fibroblasts1','Bcells1','DC1'),
            reduction = 'harmonyumap')

p <- p1/p2

ggsave(
  p,
  filename = 'result/Myeloids_subtype/M1_M2 score.pdf',
  height = 30,
  width  = 20
) 

FeaturePlot(Myeloids,features = "Tpsab1",reduction = "harmonyumap")


## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0' = 'M2',          # Lum, Cd163, Ccl7, F13a1
  '1' = 'M2',          # Wnt4, Ocstamp, Lama5, Smoc2 (Stromal-like)
  '2' = 'Enterocytes',   # Krt20, Vil1, Prap1
  '3' = 'M1',          # Il1rn, Cxcl10, Chil3, Arg1, Acod1
  '4' = 'BEC',   # Kdr, Podxl, Cdh5, Vwf
  '5' = 'Plasma',# Igkc, Igha, Jchain, Mzb1
  '6' = 'Goblet',# Cps1, Itln1, Prom1, Cldn2
  '7' = 'Cycling_Macro',    # Mki67, Top2a, Ccna2, Birc5
  '8' = 'T_cell',# Lat, Cd3d, Cd3e, Gzmb
  '9' = 'Neutrophi', # Alox15, Cxcr2, Mmp25, Ffar2
  '10' = 'DC',           # Cd209a, Flt3, Ccl17, Zbtb46
  '11' = 'DC',           # Xcr1, Clec9a, Flt3, Ly75
  '12' = 'DC',          # Ccl19, Ccl21a, Ccl21b, Il7
  '13' = 'LEC'   # Lyve1, Flt4, Prox1, Reln
)


Idents(Myeloids) <- "RNA_snn_res.0.4"

Myeloids <- RenameIdents(Myeloids, cell_cluster)
Myeloids$subtype <- Idents(Myeloids)


pdf(file.path(Myeloids_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Myeloids, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Myeloids_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 20)
DimPlot(Myeloids, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

Myeloids_cellTye_markers <- c( 
  "M0-Mac" = "Cd68",
  
  "M2-Mac" = "Cd163",
  "M2-Mac" = "Mrc1",
  
  "Caf-like_MoMa" = "Postn",
  "Caf-like_MoMa" = "Col5a1",
  
  # "cDC2" = "Cd1c", # Note: No direct mouse ortholog for human CD1C
  "cDC2" = "Fcer1a",
  
  'Epi-like_MoMa' = "Epcam",
  
  'M1-Mac' = "Cxcl9",
  'M1-Mac' = "Cxcl10",
  
  'Pro_MoMa' = 'Mki67',
  'Pro_MoMa' = 'Top2a',
  
  'Mono' = "Fcgbp",
  'Mono' = "Cdc42",
  
  "MAST" = "Kit",
  "MAST" = "Ms4a2",
  
  'Plasma-like_MoMa' = 'Mzb1',
  'Plasma-like_MoMa' = 'Xbp1',
  
  'LAMP3_DC' = "Lamp3",
  'LAMP3_DC' = 'Fscn1', 
  
  # 'pDC' = 'Clec4c', # Note: No direct mouse ortholog for human CLEC4C (BDCA-2)
  # 'pDC' = 'Lilra4', # Note: No direct mouse ortholog for human LILRA4 (ILT7)
  
  'cDC1s' = 'Clec9a',
  'cDC1s' = 'Ido1'
)

pdf(file.path(Myeloids_out_dir, "dotplot celltype main_marker.pdf"), width = 20, height = 6)
DotPlot(Myeloids, group.by = "subtype", features = Myeloids_cellTye_markers, dot.scale = 10) 
#+ theme_90angle
dev.off()


## cell props barplot  ----------------------------------------------------------
Myeloids_meta <- Myeloids[[]]

#Myeloids_meta$sample <- factor(Myeloids_meta$sample,levels=c(Myeloids_meta$sample %>% unique()))

p_col_proj <-
  Myeloids_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path("result/Myeloids_subtype", "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  Myeloids_meta %>%
  group_by(group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path(Myeloids_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)


## vocanol -------------------
Myeloids_markers <-
  FindAllMarkers(
    JoinLayers(Myeloids),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(Myeloids_markers,'Myeloids_subtype/Myeloids_diff_gene.csv')

library(EnhancedVolcano)
type<-unique(Myeloids_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- Myeloids_markers[Myeloids_markers$cluster == type, ]
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
      filename = file.path(Myeloids_out_dir, str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })


Myeloids_markers_major<-list(Myeloids = Myeloids_markers)
# scatter
p_volcanoplot_major_ls <- 
  Myeloids_markers_major %>%
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
      filename = file.path(Myeloids_out_dir, str_glue("{name} diffgene dot.pdf")),
      height = 10,
      width = 25
    )
  })


## enrrichment ---------------------------------------
Myeloids_de_genes_for_enrichmenet <- 
  Myeloids_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.15) %>%
  dplyr::filter(abs(avg_log2FC) > 0.75) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")

Myeloids_er_output<-'Myeloids_subtype/enrichment' # dir.create(Myeloids_er_output)
enrich_ls <-
  Myeloids_de_genes_for_enrichmenet %>%
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


# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Myeloids_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Myeloids_er_output, name, "Go enriment table.csv"),
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
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(Myeloids_er_output, name, "go dotplot.pdf"),
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
      filename = file.path(Myeloids_er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })


Myeloids_de_genes_id_for_enrichmenet <-
  Myeloids_de_genes_for_enrichmenet %>%
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
  Myeloids_de_genes_id_for_enrichmenet %>%
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

# saveRDS(Myeloids_er_output, 'enrich_ls.rds')

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(Myeloids_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(Myeloids_er_output, name, "KEGG enriment table.csv"),
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
      filename = file.path(Myeloids_er_output, name, "kegg dotplot.pdf"),
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
      filename = file.path(Myeloids_er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })

saveRDS(Myeloids, 'result/Myeloids_subtype/Myeloids.rds')




#5、成纤维、肌成纤维亚类注释------
Fibroblast2_out_dir<-'result/Fibroblast_subtype/' 
dir.create(Fibroblast2_out_dir)
Fibroblast2<-subset(all_data2, subset = annotation %in% c("Fibroblast","Myofibroblast"))
Fibroblast2$annotation %>% table()

## PCA --------
DefaultAssay(Fibroblast2) <- "RNA"
Fibroblast2 <- FindVariableFeatures(Fibroblast2)
Fibroblast2 <- NormalizeData(Fibroblast2)
Fibroblast2 <- ScaleData(Fibroblast2, vars.to.regress = 'nCount_RNA')
Fibroblast2 <- RunPCA(Fibroblast2)


### pca inspector ------
pdf(file.path(Fibroblast2_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Fibroblast2, ndims = 50)
dev.off()

# pdf(file.path(Fibroblast2_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
# DimHeatmap(Fibroblast2, dims = 1:5)
# dev.off()
# pdf(file.path(Fibroblast2_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
# DimHeatmap(Fibroblast2, dims = 6:10)
# dev.off()
# pdf(file.path(Fibroblast2_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
# DimHeatmap(Fibroblast2, dims = 11:15)
# dev.off()
# pdf(file.path(Fibroblast2_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
# DimHeatmap(Fibroblast2, dims = 16:20)
# dev.off()


## dimension reduction --------
dims <- 20

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb

Fibroblast2 <- RunUMAP(Fibroblast2, dims = 1:dims,return.model = TRUE)

pdf(file.path(Fibroblast2_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Fibroblast2, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

table(all_data2@meta.data$celltype, all_data2@meta.data$project)

Fibroblast2 <- Fibroblast2 %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

Fibroblast2 <- RunUMAP(
  Fibroblast2, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(Fibroblast2_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Fibroblast2, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'

### cluster -------
Fibroblast2  <-
  FindNeighbors(
    Fibroblast2 ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )

Fibroblast2 <-
  FindClusters(
    Fibroblast2,
    resolution = c(0.2, 0.3, 0.4)
  )

# dimplot
pdf(file.path(Fibroblast2_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Fibroblast2, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast2, group.by = "RNA_snn_res.0.3", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast2, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


pdf(file.path(Fibroblast2_out_dir, str_glue("FeaturePlot dim{dims} Pdgfra.pdf")), height = 6, width = 8)
FeaturePlot(Fibroblast2,feature = c("Pdgfra"),reduction = 'harmonyumap')
dev.off()


Fibroblast2$res_1.2 %>% table()
pdf(file.path(Fibroblast2_out_dir, str_glue("Dimplot dim{dims} res_1.2.pdf")), height = 6, width = 8)

DimPlot(Fibroblast2, group.by = "res_1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()

dev.off()

pdf(file.path(Fibroblast2_out_dir, str_glue("Dimplot dim{dims} celltype.pdf")), height = 6, width = 8)

DimPlot(Fibroblast2, group.by = "celltype", label = TRUE, reduction = 'harmonyumap') + NoLegend()

dev.off()




### dotplot -------
intestinal_marker_ls <- c(
  # 上皮细胞群
  'Epithelial Stem Cell' = "Lgr5",
  'Epithelial Stem Cell' = "Ascl2",
  'Epithelial Stem Cell' = "Sox9",
  
  "Epcam",
  
  'Paneth Cell' = "Mmp7",
  
  'Enterocyte' = "Slc26a3",
  'Enterocyte' = "Vil1",
  
  'Goblet Cell' = "Muc2",
  'Goblet Cell' = "Spdef",
  
  'Enteroendocrine Cell (EEC)' = "Neurod1",
  'Enteroendocrine Cell (EEC)' = "Tph1",
  'Enteroendocrine Cell (EEC)' = "Gcg",
  'Enteroendocrine Cell (EEC)' = "Sst",
  
  'Tuft Cell' = "Il25",
  'Tuft Cell' = "Trpm5",
  'Tuft Cell' = "Avil",
  
  'M Cell' = "Tnfrsf11a",
  
  # 免疫细胞群
  'CD4+ T Cell' = "Cd3d",
  'CD4+ T Cell' = "Cd3e",
  'CD4+ T Cell' = "Cd3g",
  'CD4+ T Cell' = "Cd4",
  'CD4+ T Cell' = "Foxp3",
  
  # 'CD8+ T Cell' = "Cd3d",
  # 'CD8+ T Cell' = "Cd3e",
  # 'CD8+ T Cell' = "Cd3g",
  'CD8+ T Cell' = "Cd8a",
  'CD8+ T Cell' = "Gzmb",
  
  'B Cell' = "Cd19",
  'B Cell' = "Cd79a",
  'B Cell' = "Ighm",
  
  'Plasma Cell' = "Sdc1",
  'Plasma Cell' = "Jchain",
  
  'Dendritic Cell (DC)' = "Itgax",
  'Dendritic Cell (DC)' = "Zbtb46",
  'Dendritic Cell (DC)' = "Clec9a",
  
  'Macrophage' = "Adgre1",
  'Macrophage' = "Csf1r",
  'Macrophage' = "Cd68",
  'Macrophage' = "Cx3cr1",
  
  'Innate Lymphoid Cell (ILC)' = "Il7r",
  'Innate Lymphoid Cell (ILC)' = "Id2",
  'Innate Lymphoid Cell (ILC)' = "Rorc",
  'Innate Lymphoid Cell (ILC)' = "Gata3",
  
  # 基质细胞群
  'Fibroblast' = "Col1a1",
  'Fibroblast' = "Pdgfra",
  'Fibroblast' = "Pdgfrb",
  
  'Endothelial Cell' = "Pecam1",
  'Endothelial Cell' = "Cdh5",
  'Endothelial Cell' = "Lyve1",
  
  'Pericyte' = "Cspg4",
  #'Pericyte' = "Pdgfrb",
  
  'Enteric Glial Cell (EGC)' = "Gfap",
  'Enteric Glial Cell (EGC)' = "S100b"
)

msc_marker_ls <- c(
  # MSC 经典标记基因
  'MSC Classic Markers' = "Cd248",
  'MSC Classic Markers' = "Cd34",
  'MSC Classic Markers' = "Cd44",
  'MSC Classic Markers' = "Eng",
  'MSC Classic Markers' = "Itgb1",
  'MSC Classic Markers' = "Mcam",
  'MSC Classic Markers' = "Nt5e",
  'MSC Classic Markers' = "Pdgfra",
  'MSC Classic Markers' = "Sdc2",
  'MSC Classic Markers' = "Vcam1",
  
  # MSC 分泌细胞因子
  'MSC Secreted Cytokines' = "Angpt1",
  'MSC Secreted Cytokines' = "Angpt2",
  'MSC Secreted Cytokines' = "Bmp2",
  'MSC Secreted Cytokines' = "Bmp4",
  # 'MSC Secreted Cytokines' = "Ccl11",
  # 'MSC Secreted Cytokines' = "Ccl2",
  # 'MSC Secreted Cytokines' = "Ccl7",
  # 'MSC Secreted Cytokines' = "Cxcl1",
  'MSC Secreted Cytokines' = "Cxcl10",
  'MSC Secreted Cytokines' = "Cxcl12",
  'MSC Secreted Cytokines' = "Cxcl9",
  'MSC Secreted Cytokines' = "Fgf2",
  #'MSC Secreted Cytokines' = "Gdf10",
  'MSC Secreted Cytokines' = "Gdnf",
  'MSC Secreted Cytokines' = "Hgf",
  'MSC Secreted Cytokines' = "Ido1",
  'MSC Secreted Cytokines' = "Igf1",
  'MSC Secreted Cytokines' = "Il10",
  'MSC Secreted Cytokines' = "Il1rn",
  'MSC Secreted Cytokines' = "Lif",
  'MSC Secreted Cytokines' = "Pdgfa",
  'MSC Secreted Cytokines' = "Pdgfb",
  'MSC Secreted Cytokines' = "Ptgs2",
  #'MSC Secreted Cytokines' = "Serpine1",
  'MSC Secreted Cytokines' = "Tgfb1",
  'MSC Secreted Cytokines' = "Tgfb2",
  'MSC Secreted Cytokines' = "Tgfb3",
  'MSC Secreted Cytokines' = "Timp1",
  'MSC Secreted Cytokines' = "Vegfa",
  #'MSC Secreted Cytokines' = "Vegfd",
  
  # MSC 中实际表达的细胞因子
  'MSC Expressed Cytokines' = "Ccl11",
  'MSC Expressed Cytokines' = "Ccl2",
  'MSC Expressed Cytokines' = "Ccl7",
  'MSC Expressed Cytokines' = "Cxcl1",
  'MSC Expressed Cytokines' = "Gdf10",
  'MSC Expressed Cytokines' = "Serpine1",
  'MSC Expressed Cytokines' = "Vegfd"
)



Idents(Fibroblast2)<-'RNA_snn_res.0.3'
Fibroblast2_markers <-
  FindAllMarkers(
    JoinLayers(Fibroblast2),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

Fibroblast2_top10_markers <- 
  Fibroblast2_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10) %>%
  nest() %>%
  deframe() %>% 
  map(pull, 'gene') %>%
  imap(function(genes, name){
    re <- genes %>%
      map(function(gene){
        # Fibroblastrowser()
        names(gene) = str_glue("cluster_{name}")
        gene
      })
    re
  }) %>%
  unlist() %>% 
  .[!duplicated(.)]

names(Fibroblast2_top10_markers) <- sub('.*\\.','',names(Fibroblast2_top10_markers))

Fibroblast2_top10_markers1 <-
  Fibroblast2_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10)

# dotplot
pdf(file.path(Fibroblast2_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(Fibroblast2, group.by = "RNA_snn_res.0.4", features = intestinal_marker_ls)+RotatedAxis()
DotPlot(Fibroblast2, group.by = "RNA_snn_res.0.4", features = msc_marker_ls)+RotatedAxis()
DotPlot(Fibroblast2, group.by = "RNA_snn_res.0.4", features = Fibroblast2_top10_markers)+RotatedAxis()
dev.off()




## annotation ----
cell_cluster <- c(
  '0'  = "Inflammatory Fibro",  # 炎症型成纤维细胞
  '1'  = "Reparative Fibro",   # 修复型成纤维细胞
  '2'  = "Stromal Fibro",      # 基质成纤维细胞
  '3'  = "Pericyte",           # 周细胞
  '4'  = "Plasma",             # 浆细胞
  '5'  = "Myofibroblasts",     # 肌成纤维细胞
  '6'  = "Enterocytes",        # 肠上皮细胞
  '7'  = "BEC",                # 血管内皮细胞
  '8'  = "Myofibroblasts",     # 肌成纤维细胞
  '9'  = "LEC",                # 淋巴管内皮细胞
  '10' = "BEC",                # 血管内皮细胞
  '11' = "Cycling_Fibro",      # 增殖型成纤维细胞
  '12' = "Reparative Fibro",   # 修复型成纤维细胞
  '13' = "delete"              # 剔除
)

Idents(Fibroblast2) <- "RNA_snn_res.0.4"
#Idents(Fibroblast2) <- "sub.cluster"

Fibroblast2 <- RenameIdents(Fibroblast2, cell_cluster)
Fibroblast2$subtype <- Idents(Fibroblast2)

Fib_marker_ls <- c(
  "iCAF" = 'CXCL12',
  "iCAF" = 'CXCL9',
  "iCAF" = 'CXCL10',
  "iCAF" = 'CXCL11',
  
  "mCAF" = 'POSTN',
  "mCAF" = 'COL10A1',
  "mCAF" = 'COL11A1',
  "mCAF" = 'CTHRC1',
  "mCAF" = 'MMP11',
  
  "Epi_like CAF" = 'EPCAM',
  "Epi_like CAF" = 'CDH1',
  "Epi_like CAF" = 'STMN1',
  
  "Fib" = 'CCDC80',
  "Fib" = 'SFRP4',
  
  "Mac_like CAF" = 'CD68',
  
  "Plasma_like CAF" = 'MZB1',
  "Plasma_like CAF" = 'XBP1',
  
  "meCAF" = "BNIP3" ,
  "meCAF" ='SLC2A1' ,
  "meCAF" ='ADM',
  "meCAF" ='ANGPTL4'
)

pdf(file.path(Fibroblast2_out_dir, "dotplot celltype marker.pdf"), width = 20, height = 8)
DotPlot(Fibroblast2, group.by = "subtype", features = Fib_marker_ls, dot.scale = 10) + theme_90angle
dev.off()

pdf(file.path(Fibroblast2_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
# DimPlot(Fibroblast2, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
DimPlot(Fibroblast2, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Fibroblast2_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 14)
DimPlot(Fibroblast2, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


## cell props barplot  --------
Fibroblast2_meta <- Fibroblast2[[]]

p_col_proj <-
  Fibroblast2_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(Fibroblast2_out_dir,"cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  Fibroblast2_meta %>%
  group_by(group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path(Fibroblast2_out_dir,"cell_freq_group.pdf"),
  height = 6,
  width  = 4
)

saveRDS(Fibroblast2, file.path(Fibroblast2_out_dir,'Fibroblast2.rds'))

## 映射-------------------------------
library(tibble)
library(data.table)
sub <- all_CAFs[['subtype']] %>% rownames_to_column('cells')

meta.data.t <-
  sub%>%
  right_join(
    Fibroblast2[['subtype']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    subtype2 = ifelse((as.character(subtype.x) %in% c('myCAFs', 'iCAFs', 'CAFs', 'CAFs2', 'vCAFs', 'tCAFs', 'rCAFs', 'dCAFs')),
                      as.character(subtype.x), 
                      as.character(subtype.y)))%>%
  dplyr::select(cells, subtype2) %>%
  column_to_rownames('cells')


Fibroblast2 <- AddMetaData(Fibroblast2,meta.data.t['subtype2'])
Fibroblast2$subtype_raw <- Fibroblast2$subtype
Fibroblast2$subtype <- Fibroblast2$subtype2

