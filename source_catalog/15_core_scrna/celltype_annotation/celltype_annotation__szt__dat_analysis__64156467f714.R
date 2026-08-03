library(tidyverse)
library(Seurat)
library(patchwork)

source("xenium.R")

# multi workers  -----
# library(future)
# plan("multisession", workers = 10, gc = TRUE)
# 
# options(future.globals.maxSize = 10 * 1024 ^ 3)
# 
# nbrOfFreeWorkers()


# setting variable   --------------------
data_dir <- "data"
anno_file <- "data/Annotation.csv"

# create directory ----------------
dir.create('data_for_analysis')

qc_dir <- "QC" 
dir.create(qc_dir)


dr_dir <- "dimension reduction" 
dir.create(dr_dir)

# reading data -------------
## data path and project name  -----
parse_data_path <-
  dir(data_dir, full.names = TRUE) %>% 
  parseProjName_xenium(pattern = "^output[-]?") %>%
  mutate(.counts_path = file.path(.path, "cell_feature_matrix")) %>%
  identity()

parse_data_path$.proj <- parse_data_path$.proj %>%  str_remove("^.*_")

## add id for parse_data_path --------
parse_data_path <- 
  parse_data_path %>%
  mutate(.id = seq_along(.proj))


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
  dplyr::select(.counts_path, .proj, .id) %>%
  pmap(function(.counts_path, .proj, .id){
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
      RenameCells(add.cell.id = .id)
    sc
  }) 
all_data_ls <- all_data_ls %>% purrr::discard(is.null)
names(all_data_ls) <- map_chr(all_data_ls, "project.name")

# project number, used as condition judges
proj_count <- length(all_data_ls)
proj_name  <- names(all_data_ls)


if(proj_count > 1){
  #merge
  all_data_ls[["Merge"]] <- mergeSeurat(all_data_ls)
}


## read coords -------------------
# cell coordinate and cell bundary
cells_position_ls <- 
  # with(parse_data_path, split(.path, .proj)) %>%
  parse_data_path %>%
  group_by(.proj) %>%
  nest %>%
  deframe() %>%
  map(function(x){
    xx <-
      readr::read_csv(
        file.path(x$.path, "cells.csv.gz"),
        col_types = list(cell_id = "c")
      ) %>%
      mutate(cell_id = str_c(x$.id, cell_id, sep = "_"))
    xx
  })
cells_position_ls <- cells_position_ls[proj_name]

cells_boundary_ls <- 
  # with(parse_data_path, split(.path, .proj)) %>%
  parse_data_path %>%
  group_by(.proj) %>%
  nest %>%
  deframe() %>%
  map(function(x){
    xx <-
      readr::read_csv(
        file.path(x$.path, "cell_boundaries.csv.gz"),
        col_types = list(cell_id = "c")
      ) %>%
      mutate(cell_id = str_c(x$.id, cell_id, sep = "_"))
    xx
  })
cells_boundary_ls <- cells_boundary_ls[proj_name]

# export data
saveRDS(cells_boundary_ls, file.path("data_for_analysis", "cells_boundary_ls.rds"))
saveRDS(cells_position_ls, file.path("data_for_analysis", "cells_position_ls.rds"))


# QC ------------------
## QC_feature ----------
qc_metric <- c(
  "nFeature_RNA",
  "nCount_RNA", 
  "log10GenesPerUMI"
)


## add qc metric ----------
all_data_ls <- all_data_ls %>%
  map(function(sc){
    sc[["log10GenesPerUMI"]] <- log10(sc$nFeature_RNA) / log10(sc$nCount_RNA)
    sc
  })

## nRNA nFeature percent.mt ----------
mito_qc_ls <-
  all_data_ls %>% 
  map(function(sc){
    
    # use to determine the filtering criteria
    p1 <- VlnPlot(object = sc,
                  features = qc_metric,
                  group.by  =  "project",
                  ncol = 3,
                  # log = T,
                  pt.size = 0)
    # geom_hline(yintercept = 0.985)
    p1 <-
      p1 & scale_y_continuous(
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
    p1
  })



pdf(file.path(qc_dir, "Vlnplot.pdf"), height = 5, width = 12)
mito_qc_ls %>%
  walk2(names(.), function(p, grp){
    print(p + patchwork::plot_annotation(title = grp))
  })
dev.off()

pdf(file.path(qc_dir, "Vlnplot-log.pdf"), height = 5, width = 12)
mito_qc_ls %>%
  walk2(names(.), function(p, grp){
    print(p + patchwork::plot_annotation(title = grp) & scale_y_log10())
  })
dev.off()


# export mito_qc.csv
mito_qc_path <- file.path('data', "mito_qc.csv")
if(!file.exists(mito_qc_path)){

  if(proj_count < 2){
    nFeature_cutoff <- all_data_ls[[1]][["nFeature_RNA"]] %>% unlist %>% quantile(0.9) %>% round()
    nCount_cutoff <- all_data_ls[[1]][["nCount_RNA"]] %>% unlist %>% quantile(0.9) %>% round()
  }else{
    nFeature_cutoff <- all_data_ls[["Merge"]][["nFeature_RNA"]] %>% unlist %>% quantile(0.9) %>% round()
    nCount_cutoff <- all_data_ls[["Merge"]][["nCount_RNA"]] %>% unlist %>% quantile(0.9) %>% round()
  }
  
  
  origin_mito_qc <-
    data.frame(
      proj               = proj_name,
      nFeature_RNA       = paste0(">10"), # paste0(">10; <", nFeature_cutoff),
      # nCount_RNA         = paste0(">20"), # paste0(">20; <", nCount_cutoff),
      log10GenesPerUMI   = ">0.6"
    )
  
  readr::write_csv(origin_mito_qc, mito_qc_path)
}


#############################################
######## check `mito_qc.csv`,       #########
######## and adjust QC metric value #########
######## if needed.                 #########
#############################################


## cell filter -----------
# read into qc criteria
origin_mito_qc <- 
  readr::read_csv(mito_qc_path) %>%
  group_by(proj) %>%
  nest %>%
  deframe()

parsed_mito_qc <-
  origin_mito_qc%>%
  map(function(x){
    x %>% map(
      function(col){
        col %>% 
          str_split(pattern = "[ ;]+") %>%
          .[[1]] %>%
          map(parse_qc_criteria)
      }
    )
  })


if(proj_count > 1){
  all_data_ls$Merge <- NULL
}

# befor mito filter
before_filter_qc <- all_data_ls %>%
  map2(names(.), function(sc, proj){
    # befor mito filter
    p1 <- FeatureScatter(object = sc, 
                         feature1 = "nCount_RNA", 
                         feature2 = "nFeature_RNA",
                         group.by = "project")
    
    
    p2 <- FeatureScatter(object = sc, 
                         feature1 = "nCount_RNA", 
                         feature2 = "log10GenesPerUMI",
                         group.by = "project")
    
    p <- p1 + p2 +
      plot_annotation(title = str_glue("{proj}: Before.qc.filter")) + 
      plot_layout(nrow=1, guides = "collect") &
      theme(legend.position = "top")
    p
  })

# filter out qc failed cells
all_data_ls_pass <- 
  all_data_ls %>%
  imap(qc_ls = parsed_mito_qc, function(sc, name, qc_ls){
    
    # qc参数
    if(!name %in% names(qc_ls)){
      stop(str_glue("Please check QC table, {name} not found qc criteria..."))
    }
    
    qc <- qc_ls[[name]]
    # browser()
    sc_have_name <- is.null(names(sc))
    qc_have_name <- is.null(names(qc))
    
    if(all(sc_have_name, qc_have_name)){
      stopifnot(all(names(sc) == names(qc)))
    }
    
    combined_idx <-
      qc %>%
      map2(names(.), function(conds, var){#  browser()
        # var: nFeture_RNA percent.mt
        combined_idx_per_var <-
          conds %>% 
          map(~eval_parsed_criteria(.x, others = sc[[var]], verbose = FALSE)) %>%
          purrr::reduce(`&`)
        combined_idx_per_var
      }) %>%
      purrr::reduce(`&`)
    
    if(sum(combined_idx) < 1){
      stop(str_glue("No cell pass the qc criteria, is qc criteria of {name} too much strict?"))
    }
    # return
    sc[, which(combined_idx)]
  })

# after mito filter
after_filter_qc <- all_data_ls_pass %>%
  map2(names(.), function(sc, proj){
    # befor mito filter
    p1 <- FeatureScatter(object = sc, 
                         feature1 = "nCount_RNA", 
                         feature2 = "nFeature_RNA",
                         group.by = "project") 
    
    p2 <- FeatureScatter(object = sc, 
                         feature1 = "nCount_RNA", 
                         feature2 = "log10GenesPerUMI",
                         group.by = "project")
    
    p <- p1 + p2 +
      plot_annotation(title = str_glue("{proj}: After.qc.filter")) + 
      plot_layout(nrow = 1, guides = "collect") &
      theme(legend.position = "top")
    p
  })


# print the image
qc_filter <- 
  map2(before_filter_qc, 
       after_filter_qc, 
       function(before, after){
         cowplot::plot_grid(before, after, ncol = 1)
       })



pdf(file.path(qc_dir, "qc plot - pass.pdf"), height = 12, width = 18)
qc_filter %>%
  walk2(names(.), 
        function(p, name){
          print(p)
        })
dev.off()


## cell filter table ------
cell_count_before <- all_data_ls      %>% map_dbl(ncol) %>% enframe(name = "Sample", value = "Before filter")
cell_count_after  <- all_data_ls_pass %>% map_dbl(ncol) %>% enframe(name = "Sample", value = "After filter")

# print 
cell_count_df <- cell_count_before %>% left_join(cell_count_after, by = "Sample")
write.csv(cell_count_df, file = file.path(qc_dir, 'filter table.csv'), row.names = FALSE)


# merge
all_data <- merge(all_data_ls_pass[[1]], all_data_ls_pass[-1])


# normalize ------
all_data <- NormalizeData(all_data)
all_data <- FindVariableFeatures(all_data)

# sketch data --------
all_data <-
  SketchData(
    all_data,
    ncells         = 20000L,
    sketched.assay = "sketch",
    method         = "LeverageScore",
    var.name       = "leverage.score",
    seed = 123L
  )


## dimension reduction - pca ------------
DefaultAssay(all_data) <- "sketch"
all_data <- FindVariableFeatures(all_data, nfeatures = 3000)
all_data <- ScaleData(all_data, vars.to.regress = 'nCount_RNA')

all_data <- RunPCA(all_data)

pdf(file.path(dr_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_data, ndims = 50)
dev.off()


set.seed(1234)
suset_cells <- 
  all_data[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(all_data) ) %>%
  dplyr::pull("cell_id")

pdf(file.path(dr_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 1:5)
dev.off()

pdf(file.path(dr_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 6:10)
dev.off()

pdf(file.path(dr_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 11:15)
dev.off()

pdf(file.path(dr_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 16:20)
dev.off()


# dimension reduction --------
dims <- 10


library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 50 * 1024 ^ 3) # 10Gb


all_data <- RunUMAP(all_data, dims = 1:dims, return.model = TRUE)


pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
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
  reduction = "harmony",
  reduction.name = "harmonyumap",
  reduction.key  = "harmonyUMAP_", 
  return.model = TRUE
)

pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


reduction.use <- 'pca' # 'harmony' 'pca'

# cluster ------
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
pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(all_data, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(all_data, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(all_data, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# # cell marker ----
major_marker_ls <- c(
  
  Endothelial = "PECAM1",
  Endothelial = "VWF",
  Endothelial = "AQP1",
  
  Epithelial  = "EPCAM", 
  Epithelial  = "KRT5",
  
  
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
  
  
  "AT1" = "AGER",
  "AT1" = "PDPN",
  "AT1" = "CLIC5",
  
  "AT2" = "SFTPB",
  "AT2" = "SFTPC",
  "AT2" = "SFTPD",
  "AT2" = "MUC1",
  "AT2" = "ETV5",
  
  
  # "Myefibroblast" = "PDGFRB", 
  "Myofibroblast" = "PLN", 
  "Myofibroblast" = "MYLK", 
  "Myofibroblast" = "MYH11", 
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


major_marker <- c(
  "Airway epithelial cells" = "ADH7",
  "Airway epithelial cells" = "AQP1",
  "Airway epithelial cells" = "CDH1",
  "Airway epithelial cells" = "SEC14L3",
  "Airway goblet cells" = "AGR2",
  "Airway goblet cells" = "AQP5",
  "Airway goblet cells" = "CEACAM1",
  "Airway goblet cells" = "DMBT1",
  "Airway goblet cells" = "DUSP4",
  "Mesothelial cells" = "C2",
  "Mesothelial cells" = "CALB2",
  "Mesothelial cells" = "CD44",
  # "Mesothelial cells" = "CDH1",
  "Mesothelial cells" = "DES",
  "Fibroblasts" = "COL3A1",
  "Fibroblasts" = "COL5A2",
  "Fibroblasts" = "DPT",
  "Fibroblasts" = "FN1",
  "Fibroblasts" = "GSN",
  "Basal cells (Airway progenitor cells)" = "ABI3BP",
  "Basal cells (Airway progenitor cells)" = "AQP3",
  "Basal cells (Airway progenitor cells)" = "DAPL1",
  "Basal cells (Airway progenitor cells)" = "GSTM2",
  "Basal cells (Airway progenitor cells)" = "HPGD",
  "Alveolar macrophages" = "ABCG1",
  "Alveolar macrophages" = "CCL3",
  "Alveolar macrophages" = "CD36",
  "Alveolar macrophages" = "CLEC7A",
  "Alveolar macrophages" = "CSF2RB",
  "Ciliated cells" = "APPL2",
  "Ciliated cells" = "ATP5MD",
  "Ciliated cells" = "CCDC153",
  "Ciliated cells" = "CCDC17",
  "Ciliated cells" = "CCDC181",
  "Clara cells" = "AHR",
  "Clara cells" = "ALDH1A1",
  "Clara cells" = "BPIFA1",
  "Clara cells" = "CTSE",
  "Clara cells" = "CYP2E1",
  "Immune system cells" = "ADGRE1",
  "Immune system cells" = "ARG1",
  "Immune system cells" = "BIRC5",
  "Immune system cells" = "CCL17",
  "Immune system cells" = "CCL18",
  "Endothelial cell" = "CD34",
  "Endothelial cell" = "EGFL7",
  "Endothelial cell" = "EMCN",
  "Endothelial cell" = "ESAM",
  "Endothelial cell" = "FLT1",
  "Epithelial cells" = "ANPEP",
  "Epithelial cells" = "EPCAM",
  "Epithelial cells" = "IL10",
  "Epithelial cells" = "IL6R",
  "Ionocytes" = "CFTR",
  "Ionocytes" = "CLCNKB",
  "Ionocytes" = "FOXI1",
  "Ionocytes" = "KCNMA1",
  "Ionocytes" = "SCGB1A1",
  "Pulmonary alveolar type I cells" = "AGER",
  "Pulmonary alveolar type I cells" = "AKAP5",
  # "Pulmonary alveolar type I cells" = "AQP3",
  # "Pulmonary alveolar type I cells" = "AQP5",
  "Pulmonary alveolar type I cells" = "CCN2",
  "Pulmonary alveolar type II cells" = "ABCA3",
  "Pulmonary alveolar type II cells" = "ADGRF5",
  # "Pulmonary alveolar type II cells" = "AGER",
  # "Pulmonary alveolar type II cells" = "CD36",
  "Pulmonary alveolar type II cells" = "CD3G",
  "Secretory cell" = "MUC5B",
  "Secretory cell" = "PIGR",
  # "Secretory cell" = "SCGB1A1",
  "Cancer stem cells" = "ABCG2",
  "Cancer stem cells" = "ALCAM",
  # "Cancer stem cells" = "ALDH1A1",
  "Cancer stem cells" = "BMI1",
  "Cancer stem cells" = "CD24"
)

# dotplot
pdf(file.path(dr_dir, str_glue("dotplot dims{dims} {reduction.use} major_marker.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "sketch_snn_res.0.2", features = major_marker_ls)
DotPlot(all_data, group.by = "sketch_snn_res.0.4", features = major_marker_ls)
DotPlot(all_data, group.by = "sketch_snn_res.0.6", features = major_marker_ls)
DotPlot(all_data, group.by = "sketch_snn_res.0.8", features = major_marker_ls)
dev.off()

pdf(file.path(dr_dir, str_glue("dotplot dims{dims} {reduction.use} major_marker2.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "sketch_snn_res.0.2", features = major_marker)
DotPlot(all_data, group.by = "sketch_snn_res.0.4", features = major_marker)
DotPlot(all_data, group.by = "sketch_snn_res.0.6", features = major_marker)
DotPlot(all_data, group.by = "sketch_snn_res.0.8", features = major_marker)
dev.off()

# find markers -------
Idents(all_data) <- "sketch_snn_res.0.8"
all_markers <-
  FindAllMarkers(
    JoinLayers(all_data),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

write.csv(all_markers, file = file.path(dr_dir, 'all_markers.csv'), row.names = FALSE)

all_de_genes <-
  all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  pull("gene")


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


pdf(file.path(dr_dir, str_glue("dotplot dims{dims} {reduction.use} - de gene.pdf")), width = 25, height = 6)
DotPlot(
  all_data,
  features = all_de_genes %>% unique,
  group.by = "sketch_snn_res.0.8"
) +
  RotatedAxis()
dev.off()


# annotation ----

cell_cluster <- c(
  '0'  = "Endothelial",
  '1'  = "Myofibroblast1",
  '2'  = "Myofibroblast2",
  '3'  = "Immune cell",
  '4'  = "AT1",
  '5'  = "AT2",
  '6'  = "AT2",
  '7'  = "Endothelial-Fibroblast", ##
  '8'  = "Immune cell",
  '9'  = "Pericyte",
  '10' = "Fibroblast",
  '11' = "Airway epithelial cells",
  '12' = "Immune cell",
  '13' = "Others"
)

Idents(all_data) <- "sketch_snn_res.0.8"

all_data <- RenameIdents(all_data, cell_cluster)
all_data$annotation <- Idents(all_data)


pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} annotation {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

pdf(file.path(dr_dir, str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "annotation", features = major_marker_ls)
dev.off()



# project data ------
DefaultAssay(all_data) <- "sketch"
all_data <- 
  all_data %>% 
  ProjectData(
    assay = "RNA",
    full.reduction = "harmony.full",
    
    sketched.assay = "sketch",
    sketched.reduction = "pca",# "pca", "harmony"
    
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

Idents(all_data) <- 'annotation_full'
DefaultAssay(all_data) <- "RNA"

pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} annotation - {reduction.use} - projected.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation_full", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
dev.off()

saveRDS(all_data, "all_data_projected.rds")
all_data <- readRDS("all_data_projected.rds")




