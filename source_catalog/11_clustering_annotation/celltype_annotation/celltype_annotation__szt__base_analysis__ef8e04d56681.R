library(tidyverse)
library(Seurat)
library(patchwork)

# library --------------
options(java.parameters = '-Xmx10g') # options('java.parameters')
library(rJava)
library(RBioFormats)
# source("xenium.R")

# define function ----------------
parseProjName <- 
  function(x,
           pattern = "matrix.*",
           check.unique = TRUE) {
    
    stopifnot(length(x) > 0) 
    
    df <- tibble::tibble(.path = x)
    df <- df %>% 
      mutate(
        .dir = dirname(.path),
        .base = basename(.path),
        .proj = str_remove_all(.base, pattern)
      ) %>% 
      mutate(
        .proj = str_remove_all(.proj, "[._-]?$")
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

# tool class for result of  parsed mito_qc criteria
# class: parsed_criteria
print.parsed_criteria <- function(x, ...){
  cat("parsed_criteria: ")
  x <- unclass(x)
  cat(str_c(x, collapse = " "))
  cat("\n")
}

# parse mito_qc criteria
# only support >, <, =
# x: string use as input
# values: a list start with a function,  then argument
# example:
# '> 200' would be parsed as: list('>', 200)


parse_qc_criteria <- function(x){
  
  # x first by trimed, also will be transformed a string if not
  orign_x <- x
  x <- str_trim(x)
  
  if(str_length(x) < 1){
    stop("x cannot be a  null character...")
  }
  
  # string split
  operator <- str_sub(x, end = 1)
  number   <- str_sub(x, start = 2) 
  if(str_trim(str_to_lower(number)) == "inf") {
    number <- Inf
  }else if(str_trim(str_to_lower(number)) %in% c("-inf", "- inf", "-  inf")) {
    number <- -Inf
  }else{
    number <- parse_number(number)
  }
  
  if(is.na(number)){
    stop(str_glue("Error number style after operator, your input: '{orign_x}'"))
  }
  
  # parse 
  res <- switch (operator,
                 ">" = list(">", number),
                 "<" = list("<", number),
                 "=" = list("==", number)
  )
  
  return(
    structure(res, class = "parsed_criteria")
  )
}


# tools function to eval parsed_criteria 
# parsed_criteria record two things: (1) fucntion name; (2) the argument
# this fucntion also can add one more argument to eval the function

eval_parsed_criteria  <-
  function(x, 
           others = NULL,
           position = c("before", "after"),
           verbose = TRUE
  ){
    stopifnot("parsed_criteria" %in% class(x) )
    
    position <- match.arg(position)
    
    # extract data
    func <- x[[1]]
    argument <- x[-1]
    
    
    # add any argument?
    if(!is.null(others)){
      argument <- switch(
        position,
        before = c(list(others), argument),
        after = c(argument, list(others))
      )
    }
    # do.call
    if(verbose){
      cat("Eval expression: \n  ")
      print(rlang::call2(func, !!!argument))
    }
    
    do.call(func, argument)
    
  }





# setting variable   --------------------
data_dir <- "/data/xenium/20241010__051140__FWSP20241029"
anno_file <- "data/Annotation.csv"
ROI_anno_file <- "../data/ROI_Annotation.xlsx"
panel_file <- "SAM000626 xenium-human-gene-expression-panel-gene-list.xlsx"


# reading data -------------
## data path and project name  -----
parse_data_path <-
  dir(data_dir, full.names = TRUE) %>% 
  parseProjName_xenium(pattern = "output[-]?") %>%
  mutate(.counts_path = file.path(.path, "cell_feature_matrix")) %>%
  identity()

## add id for parse_data_path --------
parse_data_path <- 
  parse_data_path %>%
  mutate(.id = seq_along(.proj))



## image w/h ratio ------
# used to fix width and heigth ratio in In Situ Spatial plot.
image_wh_r_ls <-
  with(parse_data_path, split(.path, .proj)) %>%
  map(function(x){
    imag_wh <-
      RBioFormats::read.metadata(
        file.path(x, "morphology.ome.tif")
      ) %>%
      .[[1]] %>%
      .$coreMetadata %>%
      .[c("sizeX", "sizeY")]
    
    image_wh_r <- imag_wh$sizeX/imag_wh$sizeY
    image_wh_r
  })


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
dir.create('data_for_analysis')
saveRDS(cells_boundary_ls, file.path("data_for_analysis", "cells_boundary_ls.rds"))
saveRDS(cells_position_ls, file.path("data_for_analysis", "cells_position_ls.rds"))
saveRDS(image_wh_r_ls, file.path("data_for_analysis", "image_wh_r_ls.rds"))



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


pdf(file.path("QC", "Vlnplot.pdf"), height = 5, width = 12)
mito_qc_ls %>%
  walk2(names(.), function(p, grp){
    # if(proj_count > 1){
    #   cat(str_glue("####  {grp}\n\n"))
    # }
    print(p + patchwork::plot_annotation(title = grp))
    # cat("  \n\n")
  })
dev.off()

pdf(file.path("QC", "Vlnplot-log.pdf"), height = 5, width = 12)
mito_qc_ls %>%
  walk2(names(.), function(p, grp){
    # if(proj_count > 1){
    #   cat(str_glue("####  {grp}\n\n"))
    # }
    print(p + patchwork::plot_annotation(title = grp) & scale_y_log10())
    # cat("  \n\n")
  })
dev.off()


# export mito_qc.csv
mito_qc_path <- file.path('data', "mito_qc.csv")
if(!file.exists(mito_qc_path)){
  # 上限
  # nFeature_cutoff <- all_data_ls[["Merge"]][["nFeature_RNA"]] %>% unlist %>% {mean(.)+3*sd(.)} %>% round()
  # nCount_cutoff <- all_data_ls[["Merge"]][["nCount_RNA"]] %>% unlist %>% {mean(.)+3*sd(.)} %>% round()
  
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

## cell filter -----------

# read into qc criteria
origin_mito_qc <- 
  # readxl::read_excel(mito_qc_path) %>%
  readr::read_csv(mito_qc_path) %>%
  group_by(proj) %>%
  nest %>%
  deframe()

parsed_mito_qc <-
  origin_mito_qc%>%
  map(function(x){
    # browser()
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

# 去除基因总量异常的cell
all_data_ls_pass <- 
  all_data_ls %>%
  imap(qc_ls = parsed_mito_qc, function(sc, name, qc_ls){ # browser()
    
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
      map2(names(.), function(conds, var){ #  browser()
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





pdf(file.path("QC", "qc plot - pass.pdf"), height = 12, width = 18)
qc_filter %>%
  walk2(names(.), 
        function(p, name){
          # if(proj_count > 1){
          #   cat(str_glue("####  {name}\n\n"))
          # }
          print(p)
          # cat("  \n\n")
        })
dev.off()


## cell filter table ------
cell_count_before <- all_data_ls      %>% map_dbl(ncol) %>% enframe(name = "Sample", value = "Before filter")
cell_count_after  <- all_data_ls_pass %>% map_dbl(ncol) %>% enframe(name = "Sample", value = "After filter")

# print 
cell_count_df <- cell_count_before %>% left_join(cell_count_after, by = "Sample")



# merge
all_data <- merge(all_data_ls_pass[[1]], all_data_ls_pass[-1])


# normalize ------
all_data <- NormalizeData(all_data, scale.factor = median(all_data$nCount_RNA))
all_data <- FindVariableFeatures(all_data)
# all_data <- ScaleData(all_data, vars.to.regress = 'nCount_Spatial.008um')

# sketch data --------
all_data <-
  SketchData(
    all_data,
    ncells         = 10000L,
    sketched.assay = "sketch",
    method         = "LeverageScore",
    var.name       = "leverage.score",
    seed = 123L
  )


## dimension reduction - pca ------------
dr_dir <- "dimension reduction" # dir.create(dr_dir)

DefaultAssay(all_data) <- "sketch"

# all_data <- NormalizeData(all_data, scale.factor = median(all_data$nCount_RNA))

all_data <- FindVariableFeatures(all_data, nfeatures = 1500)
VariableFeaturePlot(all_data)

all_data <- ScaleData(all_data, vars.to.regress = c('nCount_RNA'))

all_data <- RunPCA(all_data)

pdf(file.path(dr_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_data, ndims = 30)
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
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


all_data <- RunUMAP(all_data, dims = 1:dims, return.model = TRUE, n.neighbors = 30)


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


reduction.use <- 'harmony' # 'harmony' 'pca'

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
DimPlot(all_data, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# # cell marker ----

panel_gene_ls <- read.csv("XeniumPrimeHuman5Kpan_tissue_pathways_metadata.csv") %>%
  dplyr::select(gene_name, cell_type) %>%
  mutate(cell_type = cell_type %>% str_split(";") %>% map(str_trim)) %>%
  unnest(cell_type) %>%
  dplyr::filter(cell_type != "")  %>%
  group_by(cell_type) %>%
  nest() %>%
  deframe() %>%
  map(unlist, use.names = FALSE)

hepa_panel_gene_ls <-
  panel_gene_ls[c(
    "hepatocyte",
    "fibroblast",
    "hepatic stellate cell" ,
    "midzonal region hepatocyte",
    "centrilobular region hepatocyte",
    "hepatoblast",
    "endothelial cell of hepatic sinusoid",
    "endothelial cell of pericentral hepatic sinusoid"
  )]




major_marker <- c(
  'epithelial cells'  = 'EPCAM',
  'epithelial cells'  = 'KRT5',
  'epithelial cells'  = 'KRT14',
  
  Hepatocyte          = 'ALB',
  Hepatocyte          = 'APOB',
  Hepatocyte          = 'CYP3A5',
  Hepatocyte          = 'CPB2',
  
  'endothelial cells' = 'PECAM1', 
  'endothelial cells' = 'VWF', 
  
  # 'fibroblasts'       = 'COL1A1', 
  # 'fibroblasts'       = 'DCN', 
  
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "COL1A1",
  "Fibroblasts" = "DCN",
  "Fibroblasts" = "ACTA2",
  
  'B/Plasma cells'    = 'CD79', 
  'B/Plasma cells'    = 'CD19', 
  'B/Plasma cells'    = 'JCHAIN', 
  'B/Plasma cells'    = 'TNFRSF17', 
  
  'myeloid cells'       = 'CD14', 
  'myeloid cells'       = 'LYZ', 
  'myeloid cells'       = 'CSF1R',  # macrophage
  
  'T cells'           = 'CD3D', 
  'T cells'           = 'CD3E', 
  'T cells'           = 'CD3G', 
  'T cells'           = 'TRAC', 
  
  'Proliferating' = 'MKI67',
  'Proliferating' = 'TOP2A',
  
  'myocyte cells'     = 'MYH7', 
  'myocyte cells'     = 'TNNC1', 
  
  'pericytes/SMCs'    = 'RGS5',  # pericytes
  'pericytes/SMCs'    = 'TAGLN', # SMCs
  
  'mast cells'        = 'KIT',
  'mast cells'        = 'TPSAB1',
  
  monocyte = "VCAN",
  monocyte = "FCN1",
  
  macrophage = "C1QA",
  macrophage = "C1QB",
  
  netrophil = "CSF3R",
  netrophil = "G0S2",
  netrophil = "RETNLG"
)

major_marker_ls <- c(
  
  Endothelial = "PECAM1",
  Endothelial = "VWF",
  Endothelial = "AQP1",
  
  Epithelial  = "EPCAM", 
  Epithelial  = "KRT5",
  
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "COL1A1",
  "Fibroblasts" = "DCN",
  "Fibroblasts" = "ACTA2",

  "Pericytes" = "PDGFRB",
  "Pericytes" = "RGS5",
  "Pericytes" = "MCAM",


  "Myefibroblast" = "PLN", 
  "Myefibroblast" = "MYLK", 
  "Myefibroblast" = "MYH11", 

  
  
  Immune = 'PTPRC',
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  "Plasma" = "TNFRSF17",
  "Plasma" = "JCHAIN",
  "Plasma" = "MZB1",
  "Plasma" = "XBP1",
  
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
  'Neutrophi' = "FCGR3B"
)


# dotplot
pdf(file.path(dr_dir, str_glue("dotplot dims{dims} {reduction.use} major_marker.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "sketch_snn_res.0.2", features = major_marker)
DotPlot(all_data, group.by = "sketch_snn_res.0.4", features = major_marker)
DotPlot(all_data, group.by = "sketch_snn_res.0.6", features = major_marker)
DotPlot(all_data, group.by = "sketch_snn_res.0.8", features = major_marker)
dev.off()

# dotplot

pdf(file.path(dr_dir, str_glue("dotplot dims{dims} {reduction.use} major_marker_ls.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "sketch_snn_res.0.2", features = major_marker_ls)
DotPlot(all_data, group.by = "sketch_snn_res.0.4", features = major_marker_ls)
DotPlot(all_data, group.by = "sketch_snn_res.0.6", features = major_marker_ls)
DotPlot(all_data, group.by = "sketch_snn_res.0.8", features = major_marker_ls)
dev.off()

# dotplot

pdf(file.path(dr_dir, str_glue("dotplot dims{dims} {reduction.use} hepa_panel_marker.pdf")), width = 20, height = 6)
hepa_panel_gene_ls %>%
  iwalk(function(x, y){
    # (DotPlot(all_data, group.by = "sketch_snn_res.0.2", features = x) + ggtitle(y) + RotatedAxis()) %>% print()
    (DotPlot(all_data, group.by = "sketch_snn_res.0.4", features = x) + ggtitle(y) + RotatedAxis()) %>% print()
    # (DotPlot(all_data, group.by = "sketch_snn_res.0.6", features = x) + ggtitle(y) + RotatedAxis()) %>% print()
    # (DotPlot(all_data, group.by = "sketch_snn_res.0.8", features = x) + ggtitle(y) + RotatedAxis()) %>% print()
  })  
dev.off() 


# annotation ----
Idents(all_data) <- "sketch_snn_res.0.4"

cell_cluster <- c(
  '0' = "Hepatocyte",
  '1' = "T",
  '2' = "Myeloid cells",
  '3' = "Fibroblasts",
  '4' = "Endothelial cells",
  '5' = "B/Plasma cells",
  '6' = "Proliferating cells",
  '7' = "Pericytes/SMCs",
  '8' = "Epithelial cells",
  '9' = "Epithelial cells"
)

all_data <- RenameIdents(all_data, cell_cluster)
all_data$annotation <- Idents(all_data)


pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} annotation {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

pdf(file.path(dr_dir, str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "annotation", features = major_marker_ls)
dev.off()


# find markers -------
Idents(all_data) <- "sketch_snn_res.0.4"
all_markers <-
  FindAllMarkers(
    JoinLayers(all_data),
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
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  pull("gene")

# DoHeatmap(sc, features = all_de_genes, group.by = "RNA_snn_res.0.8")
if(!exists('avg_expr')) avg_expr <- AverageExpression(JoinLayers(all_data))[[1]]


all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  # arrange(desc(pct.1)) %>%
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
DefaultAssay(all_data) <- "sketch"
all_data <- 
  all_data %>% 
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

Idents(all_data) <- 'annotation_full'
DefaultAssay(all_data) <- "RNA"

pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} annotation - {reduction.use} - projected.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation_full", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
dev.off()

# saveRDS(all_data, "all_data_projected.rds")
# all_data <- readRDS("all_data_projected.rds")

pdf(file.path(dr_dir, str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2 projected.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "annotation_full", features = major_marker_ls) + RotatedAxis()
dev.off()



# cell props barplot  --------

all_data_meta <-
  all_data[[]] %>%
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


# p_box_group  <-
#   all_data_meta %>%
#   group_by(group, project, annotation_full) %>%
#   summarise(count = n()) %>%
#   mutate(freq = count/sum(count)) %>%
#   
#   ggplot(aes(x = annotation_full, y = freq, fill = group)) +
#   geom_boxplot(outlier.colour = NA) +
#   geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.25)) +
#   ggsci::scale_fill_igv() +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
#   )


# export
ggsave(
  p_col_proj,
  filename = file.path(dr_dir, "cell_freq_project.pdf"),
  height = 6,
  width  = 10
)

# ggsave(
#   p_box_group,
#   filename = file.path(dr_dir, "cell_freq_group.pdf"),
#   height = 6,
#   width  = 10
# )

# spatial plot --------------

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


