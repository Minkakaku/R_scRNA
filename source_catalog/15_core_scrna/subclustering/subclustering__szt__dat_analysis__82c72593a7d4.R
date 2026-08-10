
library(tidyverse)
library(Seurat)
library(patchwork)

# setting dir -------
# img_output <- "materials_for_rmd" # dir.create(img_output)
qc_outdir <- "1. QC" # dir.create(qc_outdir)
dr_outdir <- "2. dimension reduction and cluster"  # dir.create(dr_outdir)
sp_output <- "3. spatial plot" # dir.create(sp_output)
sp_dir <- "spatial plot" # dir.create(sp_dir)
mk_output <- "4. marker gene" # dir.create(mk_output)
er_output <- "5. enriment analysis" # dir.create(er_output)

# dir.create(img_output)
dir.create(qc_outdir)
dir.create(dr_outdir)
dir.create(sp_output)
dir.create(mk_output)
dir.create(er_output)


list.files("/data/xenium/FWSP20253952-03")

# library --------------
# copy jar file.
bf_cache_dir <- tools::R_user_dir("RBioFormats", which = "cache")
dir.create(bf_cache_dir)
bf_jar_path <- list.files(system.file("java", package = "RBioFormats"), pattern = "jar", full.names = TRUE)
bf_ver <- read.dcf(file.path(system.file(package = "RBioFormats"), "DESCRIPTION"), "BioFormats")
bf_cache_name <- sprintf("bioformats_package_%s.jar", bf_ver)
file.copy(bf_jar_path, file.path(bf_cache_dir, bf_cache_name), overwrite = TRUE)

options(timeout = 10 ^ 4)
options(java.parameters = '-Xmx10g') # options('java.parameters')
library(rJava)
library(RBioFormats)
source("xenium.R")
options(timeout = 60)


# multi workers  -----
# library(future)
# plan("multisession", workers = 10, gc = TRUE)
# options(future.globals.maxSize = 10 * 1024 ^ 3)
# 
# nbrOfFreeWorkers()

# options(timeout = 10000)

# setting variable   --------------------
data_dir <- "/data/xenium/FWSP20253952-03"
anno_file <- "data/Annotation.csv"
# ROI_anno_file <- "../data/ROI_Annotation.xlsx"
# panel_file <- "SAM000626 xenium-human-gene-expression-panel-gene-list.xlsx"


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
# 这段代码的主要目的是计算每个图像的长宽比 (width-to-height ratio)，并将其存储在列表 image_wh_r_ls 中
# 现在可以用 coord_equal 调整宽高比，这段代码可以不跑
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


# parse ROI----------------
# sample_df <-  
#   data.frame(
#     sample = c("kxz1","kxz6","kxz7","kcz2","kcz3","kcz8","gxz10","gxz11","gxz12","gcz13","gcz14","gcz15"),
#     group  = c("AirMale","AirMale","AirMale","AirFemale","AirFemale","AirFemale","HyperoxiaMale","HyperoxiaMale","HyperoxiaMale","HyperoxiaFemale","HyperoxiaFemale","HyperoxiaFemale")
#   )

### sample_df中的内容：sample为样本名，group为对应FOV名
sample_df <-  
  data.frame(
    group = c("TMA-5","TMA-5","TMA-5","TMA-5","TMA-5","TMA-5","TMA-5","TMA-5","TMA-5","TMA-5","TMA-5","TMA-5","TMA-5","TMA-5","TMA_4","TMA_4","TMA_4","TMA_4","TMA_4","TMA_4","TMA_4","TMA_4","TMA_4","TMA_4","TMA_4","TMA_4","TMA_4","TMA_4","TMA_4","TMA_3","TMA_3","TMA_3","TMA_3","TMA_3","TMA_3","TMA_3","TMA_3","TMA_3","TMA_3","TMA_3","TMA_3","TMA_3","TMA_3","TMA_3","TMA_3","TMA_3","TMA_3","TMA_3","TMA_3","TMA_3","TMA_3","TMA_1","TMA_1","TMA_1","TMA_1","TMA_1","TMA_1","TMA_1","TMA_1","TMA_1","TMA_1","TMA_1","TMA_1","TMA_1","TMA_1","TMA_1","TMA_1","TMA_1","TMA_1","TMA_1","TMA_1","TMA","TMA","TMA","TMA","TMA","TMA","TMA","TMA","TMA","TMA","TMA","TMA","TMA","TMA","TMA"),
    sample  = c("F1803468-2","Z1635989-12","Z1802235-10","Z1814675-2","Z1826868-6","Z1952876-9","Z2057500-35","Z2164107-11","Z2171766-3","Z2323818-8","Z2336500-6","Z2345459-7","Z2353385-10","ZF2304103-8","Z2047816-12","Z2118307-9","Z2123354-1-1","Z2123354-1-2","Z2123354-1-3","Z2262137-5","Z2262884-12","Z2269045-13","Z2271577-12","Z2332034-5","Z2417954-8","Z2489020-13","ZF2202916-8","ZF2207170-8","ZF2208790-1","H2510098-1","H2510098-2","Q2233699-2-1","Q2233699-2-2","Q2239242-1-1","Q2239242-1-2","Q2239242-1-3","Q2243285-1-1","Q2243285-1-2","Q2319383-1-1","Q2319383-1-2","Q2322393-1-1","Q2322393-1-2","Q2337205-2-1","Q2337205-2-2","Q2337350-1-1","Q2337350-1-2","Z1915172-7","Z1954693-9","Z2279427-2-1","Z2279427-2-2","Z2370717-2","F1706967-11","F1805049-2","Q2220872-1-1","Q2220872-1-2","Q2220872-1-3","Q2236500-1-1","Q2236500-1-2","Q2236500-1-3","Q2324341-1-1","Q2324341-1-2","Q2324341-1-3","Z1530765-11","Z1642635-14","Z2038176-15","Z2101519-1-1","Z2101519-1-2","Z2101519-1-3","Z2158905-4","Z2315540-8","Z2459257-20","F1706834-1","Q1810471-1-1","Q1810471-1-2","Q2406275-2-1","Q2406275-2-2","Q2406275-2-3","Z1430967-9","Z1736838-15","Z1944219-14","Z2007432-9","Z2224197-10","Z2314730-5","Z2354531-1","Z2489020-3","Z2489129-3")
  )


roi_ls <- list.files(
  'ROI',
  full.names = TRUE
) 

names(roi_ls) <- 
  roi_ls %>%
  basename() %>%
  str_remove("\\..*")


roi_coords_ls <-
  roi_ls %>%
  map(function(path){
    # 读取文件内容，跳过前两行
    dat <- read.csv(path,skip = 2)
    # 将列名转换为小写
    colnames(dat) <- dat %>% colnames() %>% str_to_lower()
    # 添加一个新的列 'cell'，其值与 'cell.id' 列相同
    dat<-dat %>% mutate(cell=cell.id)%>%
      # 将 'cell.id' 列转换为行名
      column_to_rownames('cell.id')
    dat
  }) 

### 根据每个样本所在芯片的顺序编号，来对应all_data_ls中的cell id
.id <- c("5_","4_","1_","4_","3_","3_","5_","5_","4_","4_","4_","3_","3_","4_","4_","4_","3_","3_","3_","3_","3_","3_","3_","3_","3_","4_","4_","4_","3_","3_","3_","3_","5_","5_","5_","5_","4_","1_","4_","5_","1_","1_","1_","3_","5_","1_","3_","5_","4_","2_","1_","4_","4_","4_","2_","2_","2_","2_","4_","1_","1_","5_","2_","2_","2_","2_","3_","3_","5_","4_","1_","2_","1_","1_","1_","5_","3_","2_","4_","2_","5_","5_","2_","2_","2_","1_")
# .id<- str_c(names(roi_coords_ls), "_")
# .id
names(.id) <- names(roi_coords_ls)
roi_coords_ls <- roi_coords_ls %>%
  map2(.id, ~ {
    .x$cell <- str_c(.y, .x$cell)
    return(.x)
  })


xenium_roi_ls <-
  roi_coords_ls %>%
  imap(function(roi, name){
    # print(roi)
    # return()
    # browser()
    # get Visium Object
    xenium_name <- sample_df$group[sample_df$sample == name]
    xenium <- all_data_ls[[xenium_name]]
    
    # subset visium
    cells <- roi$cell
    xenium_roi <- subset(xenium, cells = cells)
    xenium_roi
  })


# # 计算 roi_coords_ls 中每个样本的细胞数
# cell_counts <- sapply(roi_coords_ls, function(x) length(x$cell))
# 
# # 输出每个样本的细胞数
# cell_counts


# xenium_roi_ls$`F21-6989A`<-all_data_ls$`F21-6989A`
# xenium_roi_ls$`F22-2468A`<-all_data_ls$`F22-2468A`

# add sample_id or sample_name -------------------
xenium_roi_ls <-
  xenium_roi_ls %>%
  imap(function(xenium, name){ #browser()
    xenium$project <- name
    Idents(xenium) <- 'project'
    # RenameCells(visium, add.cell.id = name)
    xenium
  })

# visium_roi override visium
all_data_ls <- xenium_roi_ls





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

roi_position_ls <-
  roi_coords_ls %>%
  imap(function(roi, name){  #browser()
    # get Visium Object
    xenium_name <- sample_df$group[sample_df$sample == name]
    xenium <- cells_position_ls[[xenium_name]]
    
    # subset visium
    cells <- roi$cell
    position_roi <-  xenium[which(xenium$cell_id %in% cells),]
    position_roi
  })

cells_position_ls <- cells_position_ls[proj_name]
cells_position_ls <- cells_position_ls %>% purrr::discard(is.null)
cells_position_ls <- c(roi_position_ls,cells_position_ls)


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

roi_boundary_ls <-
  roi_coords_ls %>%
  imap(function(roi, name){  #browser()
    # get Visium Object
    xenium_name <- sample_df$group[sample_df$sample == name]
    xenium <- cells_boundary_ls[[xenium_name]]
    
    # subset visium
    cells <- roi$cell
    boundary_roi <- xenium[which(xenium$cell_id %in% cells),]
    boundary_roi
  })

cells_boundary_ls <- cells_boundary_ls[proj_name]
cells_boundary_ls <- cells_boundary_ls %>% purrr::discard(is.null)
cells_boundary_ls <- c(roi_boundary_ls,cells_boundary_ls)


# export data
dir.create("data_for_analysis") 
saveRDS(cells_boundary_ls, file.path("data_for_analysis", "cells_boundary_ls.rds"))
saveRDS(cells_position_ls, file.path("data_for_analysis", "cells_position_ls.rds"))
saveRDS(all_data_ls, file.path("data_for_analysis", "all_data_ls.rds"))

# all_data <- merge(all_data_ls[[1]], all_data_ls[-1])
# saveRDS(all_data, "data_for_analysis/all_data.rds")




























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

# if(proj_count > 1){
#   mito_qc_ls[["Merge_vlnplot"]] <- mito_qc_ls[["Merge"]]
#   mito_qc_ls[["Merge"]] <- NULL
#   mito_qc_ls[["Merge_density"]] <-
#     all_data_ls$Merge[[]] %>%
#     dplyr::select(project, all_of(qc_metric)) %>%
#     pivot_longer(
#       -project,
#       names_to = "QC_metric",
#       values_to = "Value"
#     ) %>%
#     mutate(QC_metric = factor(QC_metric, levels = qc_metric)) %>%
#     ggplot(aes(x = Value, fill = project)) +
#     geom_density(alpha = 0.4)+
#     facet_wrap(~QC_metric, scales = "free", ncol = 3)+
#     ggpubr::theme_pubclean() +
#     labs(x = NULL,
#          y = NULL)
# }


pdf(file.path("1.QC", "Vlnplot.pdf"), height = 5, width = 12)
mito_qc_ls %>%
  walk2(names(.), function(p, grp){
    # if(proj_count > 1){
    #   cat(str_glue("####  {grp}\n\n"))
    # }
    print(p + patchwork::plot_annotation(title = grp))
    # cat("  \n\n")
  })
dev.off()

pdf(file.path("1.QC", "Vlnplot-log.pdf"), height = 5, width = 12)
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





pdf(file.path("1.QC", "qc plot - pass.pdf"), height = 12, width = 18)
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

# export
cell_count_df %>% 
  rownames_to_column(" ") %>%
  openxlsx::write.xlsx(
    file = file.path("1.QC", "cell count after filter.xlsx"),
    asTable = TRUE,
    withFilter = FALSE
  )

# export data
saveRDS(all_data_ls_pass, file.path("data_for_analysis", "all_data_ls_pass.rds"))


# merge
all_data <- merge(all_data_ls_pass[[1]], all_data_ls_pass[-1])
all_data_sp <- merge(all_data_ls[[1]], all_data_ls[-1])

# normalize ------
all_data <- NormalizeData(all_data)
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
DefaultAssay(all_data) <- "sketch"
all_data <- FindVariableFeatures(all_data, nfeatures = 3000)
all_data <- ScaleData(all_data, vars.to.regress = 'nCount_RNA')

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
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 50 * 1024 ^ 3) # 10Gb


dims <- 15


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
    resolution = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
  )

# dimplot
pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.3", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.5", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.7", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.3", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.5", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.7", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# # cell marker ----
major_marker_ls <- c(
  "Epithelial"  = "Epcam",
  AT1 = "Ager",
  AT1 = "Pdpn",
  AT1 = "Csf1r",
  AT1 = "Mafb",
  AT1 = "Batf3",
  AT1 = "Krt79",
  AT1 = "Msr1",
  AT2 = "Abca3",
  AT2 = "Wif1",
  AT2 = "Etv5",
  "Basel cell" = "Ngfr",
  "B cell" = "Cd79a",
  "B cell" = "Cd79b",
  "B cell" = "Cd19",
  "B cell" = "Ms4a1",
  "T cells" = "Cd3d",
  "T cells" = "Cd3e",
  "T cells" = "Cd3g",
  "CD4T cell" = "Cd4",
  "CD4T cell" = "Cxcl13",
  "CD4T cell" = "Il21",
  "CD8T cell" = "Cd2",
  "CD8T cell" = "Cd8a",
  "CD8T cell" = "Cd8b1",
  "CD8T cell" = "Il2ra",
  "CD8T cell" = "Ptprc",
  Cilia = "Gsta1",
  Cilia = "Pifo",
  Cilia = "Foxj1",
  # "Club cell" = "PECAM1",
  Cilliated = "Muc5ac",
  "Dendritic cell" = "Ccl17",
  "Dendritic cell" = "Cd83",
  "Dendritic cell" = "Clec9a",
  "Dendritic cell" = "Ido1",
  Endothelial = "Pecam1",
  Endothelial = "Cdh5",
  Endothelial = "Vwf",
  Endothelial = "Plvap",
  Fibroblast = "Col1a1", 
  Fibroblast = "Col1a2",
  Fibroblast = "Pdgfra",
  Fibroblast = "Lum",
  Granulocyte = "Csf3r",
  "Mono/Mac" = "Cd68",
  "Mono/Mac" = "Fcgr1",
  "Mono/Mac" = "Cd163",
  "Mono/Mac" = "Cd14",
  "Mono/Mac" = "Fcgr3",
  "Mono/Mac" = "Vcan",
  "Mono/Mac" = "Mrc1",
  "Mast cell" = "Gata2",
  "Mast cell" = "Kit",
  "Mast cell" = "Ms4a2",
  "Mast cell" = "Tpsab1",
  "Mast cell" = "Tpsb2",
  "NK cell" = "Ncam1",
  "NK cell" = "Klrd1",
  Pericytes = "Pdgfrb",
  Pericytes = "Mcam",
  "Smooth muscle cells" = "Pln",
  "Smooth muscle cells" = "Myh11",
  "Smooth muscle cells" = "Tagln"
)

fibro_markers %in% rownames(all_data)

# dotplot
pdf(file.path(dr_dir, str_glue("dotplot dims{dims} {reduction.use} major_marker2.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "sketch_snn_res.0.2", features = major_marker_ls)
DotPlot(all_data, group.by = "sketch_snn_res.0.3", features = major_marker_ls)
DotPlot(all_data, group.by = "sketch_snn_res.0.4", features = major_marker_ls)
DotPlot(all_data, group.by = "sketch_snn_res.0.5", features = major_marker_ls) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
DotPlot(all_data, group.by = "sketch_snn_res.0.6", features = major_marker_ls)
DotPlot(all_data, group.by = "sketch_snn_res.0.7", features = major_marker_ls)
DotPlot(all_data, group.by = "sketch_snn_res.0.8", features = major_marker_ls)
dev.off()

# annotation ----
Idents(all_data) <- "sketch_snn_res.0.5"



cell_cluster <- c(
  '0'  = "Endothelial",
  '1'  = "SMC",
  '2'  = "AT2",
  '3'  = "Fibroblast",
  '4'  = "Mono_Mac",
  '5'  = "AT1",
  '6'  = "SMC",
  '7'  = "Cilia",
  '8'  = "T cell",
  '9'  = "AT1",
  '10' = "Endothelial",
  '11' = "B cell",
  '12' = "Fibroblast",
  '13' = "SMC",
  '14' = "Endothelial",
  '15' = "Endothelial"
)


all_data <- RenameIdents(all_data, cell_cluster)
all_data$annotation <- Idents(all_data)


pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} annotation {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

DefaultAssay(all_data) <- "sketch"
Idents(all_data) <- "annotation"
pdf(file.path(dr_dir, str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "annotation", features = major_marker_ls, dot.scale = 10) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
dev.off()


# find markers -------
Idents(all_data) <- "sketch_snn_res.0.5"
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
  width  = 20
)

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
  setNames(., paste('cluster', (1:ncol(.) - 1) )) %>%  t() %>%
  pheatmap::pheatmap(
    filename = "4. marker gene/all_markers_heatmap.pdf", width = 20, height = 5,
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

saveRDS(all_data, "data_for_analysis/all_data_projected.rds")
# all_data <- readRDS("all_data_projected.rds")

pdf(file.path(dr_dir, str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2 projected.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "annotation_full", features = major_marker_ls, dot.scale = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
dev.off()

# # 发现T与plasma标错顺序了
# all_data$annotation_full2 <- all_data$annotation_full
# all_data$annotation_full2[all_data$annotation_full == "T"]      <- "Plasma"
# all_data$annotation_full2[all_data$annotation_full == "Plasma"] <- "T"
# all_data$annotation_full  <- all_data$annotation_full2

saveRDS(all_data, 'data_for_analysis/all_data.rds')

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

# spatial plot --------------

all_data_meta <-
  all_data[[]] %>%
  rownames_to_column("cell_id")

plan <- c("#00B4F0","#7CAE00","#FF64B0","#619CFF","#F8766D","#00BA38","#00C08B","#DE8C00","#C77CFF")
plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
          "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", "#de4307", 
          "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", "#fa6e57", "#f69e53", 
          "#248888", "#999999", "#E7475E", "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", 
          "#FF6D31", "#f1ac9d", "#f06966", "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", 
          "#FF7A5A", "#8ED2C9", "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", 
          "#ffc300", "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
          "#f17d80", "#737495")

sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation_full, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() +
      scale_fill_manual(values = plan_rctd)+
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5, color = "black"))) +
      theme_void()
    p
  })


# export image
sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(sp_output, str_glue("spatial plot {grp}.pdf")),
      height = 12,
      width  = 12
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })



# Fibroblast subcluster ---------
fb_out_dir <- "fibroblast_subcluster" # dir.create(fb_out_dir)

sc_fibro <- all_data[, all_data$annotation_full == "Fibroblast"]
DefaultAssay(sc_fibro) <- "RNA"
sc_fibro[['sketch']] <- NULL


sc_fibro <-
  sc_fibro %>%
  # SCTransform()
  FindVariableFeatures(nfeatures = 1000) %>%
  NormalizeData() %>%
  ScaleData(vars.to.regress = 'nCount_RNA')



## pca -----
sc_fibro <- sc_fibro %>% RunPCA()

pdf(file.path(fb_out_dir, 'pca plot-fibro.pdf'), height = 6, width = 8)
ElbowPlot(sc_fibro, ndims = 30)
dev.off()

set.seed(1234)
suset_cells_fibro <- 
  sc_fibro[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(sc_fibro) ) %>%
  dplyr::pull("cell_id")


pdf(file.path(fb_out_dir, 'fibro dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(sc_fibro[, suset_cells_fibro], dims = 1:5)
dev.off()
pdf(file.path(fb_out_dir, 'fibro dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(sc_fibro[, suset_cells_fibro], dims = 6:10)
dev.off()
pdf(file.path(fb_out_dir, 'fibro dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(sc_fibro[, suset_cells_fibro], dims = 11:15)
dev.off()
pdf(file.path(fb_out_dir, 'fibro dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(sc_fibro[, suset_cells_fibro], dims = 16:20)
dev.off()


## dimension reduction --------
fibro_dims <- 10 #  20 # 20


# library(future)
# plan(multisession, workers = 10)
# options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


sc_fibro <- RunUMAP(sc_fibro, dims = 1:fibro_dims)


pdf(file.path(fb_out_dir, str_glue("fibro Dimplot dim{fibro_dims} project.pdf")), height = 6, width = 8)
DimPlot(sc_fibro, group.by = "project", label = TRUE) + NoLegend()
dev.off()


sc_fibro <- sc_fibro %>% 
  harmony::RunHarmony(
    dims.use = 1:fibro_dims,
    group.by.vars = "project",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

sc_fibro <- RunUMAP(
  sc_fibro, 
  dims = 1:fibro_dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(fb_out_dir, str_glue("fibro Dimplot dim{fibro_dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(sc_fibro, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'pca' # 'harmony' 'pca'

## cluster ------
sc_fibro <-
  FindNeighbors(
    sc_fibro,
    reduction = reduction.use,
    dims = 1:fibro_dims,
    k.param = 10
  )

sc_fibro <-
  FindClusters(
    sc_fibro,
    resolution = c(0.05, 0.1, 0.2) #, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(fb_out_dir, str_glue("fibro Dimplot dim{fibro_dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(sc_fibro, group.by = "RNA_snn_res.0.05", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_fibro, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_fibro, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_fibro, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_fibro, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_fibro, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(sc_fibro, group.by = "RNA_snn_res.0.05", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_fibro, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_fibro, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(sc_fibro, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(sc_fibro, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(sc_fibro, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# pdf(file.path(fb_out_dir, str_glue("fibro Dimplot dim{fibro_dims} group2 - {reduction.use}.pdf")), height = 6, width = 8)
# DimPlot(sc_fibro, group.by = "group2")
# dev.off()

## find markers -------
Idents(sc_fibro) <- "RNA_snn_res.0.2"
DefaultAssay(sc_fibro) <- "RNA"

# all_markers <-
#   FindAllMarkers(
#     JoinLayers(sc_fibro),
#     test.use = "wilcox",
#     logfc.threshold = 0.25,
#     min.pct = 0.1,
#     only.pos = TRUE
#   )
all_markers <-
  FindAllConservedMarkers(
    JoinLayers(sc_fibro),
    grouping.var = "project",
    test.use = "wilcox",
    # logfc.threshold = 0.25,
    # min.pct = 0.1,
    only.pos = TRUE
  ) %>%
  map(rownames_to_column, "gene") %>%
  enframe("cluster") %>%
  unnest(value) %>%
  rowwise() %>%
  mutate(
    avg_log2FC = mean(c_across(ends_with("avg_log2FC")), na.rm = TRUE),
    pct.1 = mean(c_across(ends_with("pct.1")), na.rm = TRUE)
  ) %>%
  mutate(p_val_adj  = minimump_p_val)

all_de_genes <-
  all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  # dplyr::filter(pct.1 > 0.25) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  pull("gene")

# DoHeatmap(sc, features = all_de_genes, group.by = "RNA_snn_res.0.8")
if(!exists('avg_expr')) avg_expr <- AverageExpression(JoinLayers(sc_fibro))[[1]]

# export all_markers
all_markers %>%
  as.data.frame() %>%
  openxlsx::write.xlsx(
    file = file.path(fb_out_dir, "all_marker_table.xlsx"),
    asTable = TRUE
  )

all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  # arrange(p_val_adj, .by_group = TRUE) %>%
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

avg_expr[all_de_genes %>% unique, -(7:9)] %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  setNames(., paste('cluster', (1:ncol(.) - 1) )) %>%  t %>%
  pheatmap::pheatmap(
    scale = "column",
    cluster_row = FALSE, # avg_expr[VariableFeatures(all_data), ] %>% t %>% dist %>% hclust,
    cluster_col = FALSE,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_col = NA
  )

fibro_genes <- c(
  'VIM', 'THY', 'S100A4', 'DCN',
  'COL1A1', 'ACTA2', 'PDPN', 'FAP'
)

avg_expr[rownames(avg_expr) %in% fibro_genes, ] %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  setNames(., paste('cluster', (1:ncol(.) - 1) )) %>%  t %>%
  pheatmap::pheatmap(
    scale = "column",
    cluster_row = FALSE, # avg_expr[VariableFeatures(all_data), ] %>% t %>% dist %>% hclust,
    cluster_col = FALSE,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_col = NA
  )

# saveRDS(sc_fibro, file = file.path(fb_out_dir, "sc_fibro.rds"))
## remove 6 7 8--------
sc_fibro2 <- sc_fibro[, !sc_fibro$RNA_snn_res.0.2 %in% 6:8]

pdf(file.path(fb_out_dir, str_glue("fibro Dimplot dim{fibro_dims} cluster - {reduction.use} - trimmed.pdf")), height = 6, width = 8)
DimPlot(sc_fibro2, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
dev.off()


# heatmap 7-9是小群，剔除
pdf(file.path(fb_out_dir, str_glue("fibro marker heatmapdim dim{fibro_dims} - {reduction.use} - trimmed.pdf")), height = 4, width = 8)
avg_expr[all_de_genes %>% unique %>% setdiff(c('CCDC80', 'TCF21')), -(7:9)] %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  setNames(., paste('cluster', (1:ncol(.) - 1) )) %>%  t %>%
  pheatmap::pheatmap(
    scale = "column",
    cluster_row = FALSE, # avg_expr[VariableFeatures(all_data), ] %>% t %>% dist %>% hclust,
    cluster_col = FALSE,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_col = NA
  )
dev.off()

# cell props 
sc_fibro_cell_props_df <-
  sc_fibro2[[]] %>% 
  dplyr::select(project, group, RNA_snn_res.0.2) %>%
  group_by(project, group, RNA_snn_res.0.2) %>%
  summarise(count = n()) %>%
  group_by(project) %>%
  mutate(freq_project = count/sum(count)) %>%
  group_by(group) %>%
  mutate(freq_group = count/sum(count))

p_sc_fibro_props_project <-
  sc_fibro_cell_props_df %>%
  ggplot(aes(x = project, y = freq_project, fill = RNA_snn_res.0.2)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggsci::scale_fill_bmj()

p_sc_fibro_props_group <-
  sc_fibro_cell_props_df %>%
  ggplot(aes(x = group, y = freq_group, fill = RNA_snn_res.0.2)) +
  geom_col() +
  theme_bw() +
  ggsci::scale_fill_bmj()

pdf(file.path(fb_out_dir, str_glue("fibro - cell props baplot - trimmed.pdf")), height = 4, width = 8)
p_sc_fibro_props_project %>% print()
p_sc_fibro_props_group   %>% print()
dev.off()


## dotplot ----------
FeaturePlot(
  sc_fibro,
  features = c(
    "ACTA2", "S100A4", "COL1A1", "VIM", 
    "PDGFRA", "PDGFRB", "THY1",
    "DCN", "FAP", "PDPN"
  ),
  reduction = "harmonyumap"
)


DotPlot(
  sc_fibro,
  features = c("ACTA2", "S100A4", "COL1A1", "VIM", "PDGFRA", "PDGFRB", "CD3D", "CD8A"),
  group.by = "RNA_snn_res.0.6"
)


DimPlot(sc_fibro, reduction = "harmonyumap", label = TRUE) + NoLegend() 



# testing analysis below -------------

# spatial plot of some gene ------
sp_genes <- c('EGFR', 'EPCAM', 'ASCL1', 'SYP', 'TP63')
sp_p_ls <-
  cells_boundary_ls %>%
  # .[1] %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    all_data <- all_data_ls[[project]]
    dat_expr <- FetchData(all_data, vars = sp_genes) %>% rownames_to_column("cell_id")
    
    
    dat <-
      cells_boundary %>%
      right_join(dat_expr, by = "cell_id")
    
    p_ls <-
      sp_genes %>%
      set_names(., .) %>%
      map(function(gene){
        q_val_rng <- quantile(dat[[gene]], c(0.01, 0.99))
        dat[[gene]] <- pmax(pmin(dat[[gene]], q_val_rng[2]), q_val_rng[1])
        
        p <- dat %>%
          ggplot(aes(x = vertex_x, y = vertex_y)) +
          geom_polygon(aes(fill = .data[[gene]], group = cell_id))  +
          scale_y_reverse() +
          coord_equal() +
          # ggsci::scale_fill_igv() +
          scale_fill_gradientn(colours = Seurat:::SpatialColors(30), transform = "log1p") +
          theme_void() +
          labs(
            title = project
          )
        p
      })
    p_ls
  })

# export
pdf(file.path(str_glue("spatial plot - gene.pdf")), height = 18, width = 18)
print(sp_p_ls)
dev.off()

sp_p_ls$FWSP20250913 %>% 
  iwalk(function(p, gene){
    pdf(file.path(str_glue("spatial plot - {gene}.pdf")), height = 18, width = 18)
    print(p)
    dev.off()
  })


# some nFeature_RNA metric -------
all_data_ls %>%
  map(function(x){
    (x$nFeature_RNA > 1000) %>% table
  }) %>%
  bind_rows()


all_data_ls %>%
  map(function(x){
    (x$nFeature_RNA > 500) %>% table
  }) %>%
  bind_rows() %>%
  dplyr::rename(
    "< 500" = 1,
    "> 500" = 2
  )


all_data_ls %>%
  map(function(x){
    summary(x$nFeature_RNA)
  })  %>%
  enframe("sample") %>%
  mutate(value = bind_rows(value)) %>%
  unnest(value)


# Spatial plot by nFeature_RNA --------

sp_p_ls <-
  cells_boundary_ls %>%
  # .[1] %>%
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    all_data_meta <- all_data_ls[[project]][[]] %>% rownames_to_column('cell_id')
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = nFeature_RNA, group = cell_id))  +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_gradient(low = "gray90", high = "red") +
      theme_void() +
      labs(
        title = project
      )
    p
  })

# export
sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path("QC", str_glue("spatial plot nFeature {grp}.pdf")), height = 18, width = 18)
    print(p)
    dev.off()
  })


# two samples analysis ---------
dir.create("two samples analysis")
all_data_ls %>% names

all_data_ls_susbet <- all_data_ls[c("X_202413076_32", "Y_2241905_12")]


all_data_ls_susbet <- 
  all_data_ls_susbet %>%
  map(function(x){
    subset(
      x,
      log10GenesPerUMI < 0.99 & log10GenesPerUMI > 0.9 
    )
  })

## normalize ------
all_data_ls_susbet <-
  all_data_ls_susbet %>%
  map(NormalizeData) %>%
  map(FindVariableFeatures) %>%
  map(ScaleData, vars.to.regress = 'nCount_RNA')



all_data_ls_susbet <-
  all_data_ls_susbet %>%
  map(RunPCA)


all_data_ls_susbet %>%
  iwalk(function(sc, name){# browser()
    
    pdf(file.path("two samples analysis", str_glue('pca plot {name}.pdf')), height = 6, width = 8)
    ElbowPlot(sc, ndims = 30) %>% print()
    dev.off()
    
    
    set.seed(1234)
    suset_cells <-
      sc[[]] %>%
      rownames_to_column("cell_id") %>%
      group_by(project) %>%
      slice_sample(prop = 50000/ncol(sc) ) %>%
      dplyr::pull("cell_id")
    
    pdf(file.path("two samples analysis", str_glue('dimheatmap 1-5 {name}.pdf')), height = 6, width = 10)
    DimHeatmap(sc[, suset_cells], dims = 1:5)
    dev.off()
    
    pdf(file.path("two samples analysis", str_glue('dimheatmap 6-10 {name}.pdf')), height = 6, width = 10)
    DimHeatmap(sc[, suset_cells], dims = 6:10)
    dev.off()
    
    pdf(file.path("two samples analysis", str_glue('dimheatmap 11-15 {name}.pdf')), height = 6, width = 10)
    DimHeatmap(sc[, suset_cells], dims = 11:15)
    dev.off()
    
    pdf(file.path("two samples analysis", str_glue('dimheatmap 16-20 {name}.pdf')), height = 6, width = 10)
    DimHeatmap(sc[, suset_cells], dims = 16:20)
    dev.off()
  })

## dimension reduction -------
two_samples_dims <- c(
  X_202413076_32 = 10,
  Y_2241905_12   = 9
)


all_data_ls_susbet <- 
  all_data_ls_susbet %>%
  imap(function(sc, name){
    dims <- two_samples_dims[name]
    # sc <- RunUMAP(sc, dims = 1:dims)
    
    pdf(file.path("two samples analysis", str_glue("Dimplot dim{dims} project {name}.pdf")), height = 6, width = 8)
    print(DimPlot(sc, group.by = "project", label = TRUE, reduction = "umap") + NoLegend())
    dev.off()
    
    sc
  })


## cluster -------
all_data_ls_susbet <- 
  all_data_ls_susbet %>%
  imap(function(sc, name){
    reduction.use <- 'pca' # 'harmony' 'pca'
    
    dims <- two_samples_dims[name]
    
    # cluster  
    sc <-
      FindNeighbors(
        sc,
        reduction = reduction.use,
        dims = 1:dims,
        k.param = 20
      )
    
    sc <-
      FindClusters(
        sc,
        resolution = c(0.2, 0.4, 0.6, 0.8)
      )
    
    # dimplot
    pdf(file.path("two samples analysis", str_glue("Dimplot dim{dims} cluster - {reduction.use} {name}.pdf")), height = 6, width = 8)
    print(DimPlot(sc, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend())
    print(DimPlot(sc, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend())
    print(DimPlot(sc, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend())
    print(DimPlot(sc, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend())
    dev.off()
    
    sc
  })


## dotplot -------
major_marker_ls <- c(
  
  Endothelial = "PECAM1",
  Endothelial = "VWF",
  Endothelial = "AQP1",
  
  Epithelial  = "EPCAM", 
  Epithelial  = "KRT5",
  
  "Fibroblasts" = "SFRP4", 
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "COL1A1",
  # "Fibroblasts" = "ACTA2",
  
  # "Pericytes" = "PDGFRB",
  
  
  "Smooth muscle cells" = "PDGFRB", 
  "Smooth muscle cells" = "MYH11", 
  "Smooth muscle cells" = "ACTA2",
  
  "Breast glandular cells" = "SERPINA3", 
  "Breast glandular cells" = "DSC2",
  
  Immune = 'PTPRC',
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  "Plasma" = "TNFRSF17",
  "Plasma" = "JCHAIN",
  
  Plasmablast = "MZB1",
  Plasmablast = "XBP1",
  
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
  
  "Monocytes" = "CD14", 
  "Monocytes" = "FCGR3A",
  
  "Dendritic cell" = "CD83",
  'Dendritic cell' = "LILRA4",  
  
  'Dendritic cell' = "CLEC9A",
  'Dendritic cell' = "IDO1",
  
  'Dendritic cell' = "CD1C",
  'Dendritic cell' = "FCER1A",
  
  'Fat cell' = 'ADIPOQ'
)


# dotplot
all_data_ls_susbet %>%
  iwalk(function(sc, name){
    reduction.use <- 'pca' # 'harmony' 'pca'
    dims <- two_samples_dims[name]
    
    pdf(file.path('two samples analysis', str_glue("dotplot dims{dims} major_marker - {reduction.use} {name}.pdf")), width = 25, height = 6)
    print(DotPlot(sc, group.by = "RNA_snn_res.0.2", features = major_marker_ls))
    print(DotPlot(sc, group.by = "RNA_snn_res.0.4", features = major_marker_ls))
    print(DotPlot(sc, group.by = "RNA_snn_res.0.6", features = major_marker_ls))
    print(DotPlot(sc, group.by = "RNA_snn_res.0.8", features = major_marker_ls))
    dev.off()
  })


## find markers ---------

all_markers_ls_susbet <-
  all_data_ls_susbet %>%
  map(function(sc){
    Idents(sc) <- "RNA_snn_res.0.2"
    all_markers <-
      FindAllMarkers(
        sc,
        test.use = "wilcox",
        logfc.threshold = 0.25,
        min.pct = 0.1,
        only.pos = TRUE
      )
    all_markers
  })



all_de_genes_ls_susbet <-
  all_markers_ls_susbet %>%
  map(function(all_markers){
    all_markers %>% 
      group_by(cluster) %>%
      dplyr::filter(p_val_adj < 0.05) %>%
      dplyr::filter(pct.1 > 0.25) %>%
      arrange(desc(avg_log2FC), .by_group = TRUE) %>%
      slice_head(n = 10) %>%
      pull("gene")
  })

all_markers_ls_susbet %>%
  walk(function(all_markers){
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
  })


de_heatmap_ls <-
  all_data_ls_susbet %>%
  imap(function(sc, name){
    Idents(sc) <- "RNA_snn_res.0.2"
    
    de_genes <- all_de_genes_ls_susbet[[name]]
    
    avg_expr <- AverageExpression(sc)[['RNA']]
    
    p <-
      avg_expr[de_genes, ] %>% 
      as.matrix() %>% 
      as.data.frame() %>% 
      setNames(., paste('cluster', (1:ncol(.) - 1) )) %>%  t %>%
      pheatmap::pheatmap(
        scale = "column",
        cluster_row = avg_expr[VariableFeatures(sc), ] %>% t %>% dist %>% hclust,
        cluster_col = FALSE,
        color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
        border_col = NA,
        silent = TRUE
      )
    cowplot::plot_grid(p$gtable)
  })


## annotation -------

annotation_ls <- list(
  "X_202413076_32" = c(
    '0'  = "Fibroblast",
    '1'  = "CD14 monocyte",
    '2'  = "Plasmablast",
    '3'  = "Endothelial",
    '4'  = "CD8T",
    '5'  = "Mast cell",
    '6'  = "Fat cell",
    '7'  = "Others"
  ),
  "Y_2241905_12"   = c(
    '0'  = "CD8T",
    '1'  = "Others",
    '2'  = "Plasmablast",
    '3'  = "Macrophage",
    '4'  = "B",
    '5'  = "CD14-dim monocyte",
    '6'  = "Epithelial",
    '7'  = "Endothelial",
    '8'  = "Fibroblast",
    '9'  = "Proliferating lymphocyte"
  )
)


# anno
all_data_ls_susbet <-
  all_data_ls_susbet %>%
  imap(function(sc, name){
    reduction.use <- "pca"
    cell_cluster <- annotation_ls[[name]]
    dims <- two_samples_dims[name]
    
    Idents(sc) <- "RNA_snn_res.0.2"
    sc <- RenameIdents(sc, cell_cluster)
    sc$annotation <- Idents(sc)
    
    
    
    pdf(file.path("two samples analysis", str_glue("Dimplot dim{dims} annotation {reduction.use} .pdf")), height = 6, width = 8)
    print(DimPlot(sc, group.by = "annotation", label = TRUE, reduction = 'umap') + NoLegend())
    dev.off()
    
    pdf(file.path("two samples analysis", str_glue("dotplot dims{dims} {reduction.use} annotation major_marker.pdf")) , width = 25, height = 6)
    print(DotPlot(sc, group.by = "annotation", features = major_marker_ls))
    dev.off()
    
    sc
  })

## spatial plot -----
sp_p_ls_subset <-
  cells_boundary_ls %>%
  .[names(all_data_ls_susbet)] %>%
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    all_data_meta <- all_data_ls_susbet[[project]][[]] %>% rownames_to_column('cell_id')
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation, group = cell_id))  +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      ggsci::scale_fill_igv() +
      # scale_fill_gradient(low = "gray90", high = "red") +
      theme_void() +
      labs(
        title = project
      )
    p
  })

# export
sp_p_ls_subset %>% 
  iwalk(function(p, grp){
    pdf(file.path("two samples analysis", str_glue("spatial plot annotation {grp}.pdf")), height = 18, width = 18)
    print(p)
    dev.off()
  })


