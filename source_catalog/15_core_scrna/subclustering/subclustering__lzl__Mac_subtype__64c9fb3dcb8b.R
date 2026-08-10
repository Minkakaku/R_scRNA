dir.create("rds_data")
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(purrr)
library(tidyverse)
out_dir <- "/mnt/NexData_file/home/lzl/20260625_XE2510095462-01/"
all_data2_TLS <- readRDS("/mnt/NexData_file/home/lzl/XE2510095462-01data/all_data2_TLS.rds")
ALL_T <- readRDS("/mnt/NexData_file/home/lzl/XE2510095462-01data/All_T.rds")
cells_boundary_ls_4slides_filtered <- readRDS("/mnt/NexData_file/home/lzl/XE2510095462-01data/cells_boundary_ls_4slides_filtered.rds")
cells_position_ls_4slides_filtered <- readRDS("/mnt/NexData_file/home/lzl/XE2510095462-01data/cells_position_ls_4slides_filtered.rds")
table(ALL_T$group)
table(all_data2_TLS$sample)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

####Mac
ALL_Mac <- subset(all_data2_TLS, annotation == 'Macrophage')


## PCA --------
DefaultAssay(ALL_Mac) <- "RNA"
ALL_Mac <- FindVariableFeatures(ALL_Mac)
ALL_Mac <- ScaleData(ALL_Mac, vars.to.regress = 'nCount_RNA')
ALL_Mac <- RunPCA(ALL_Mac)

Mac_out_dir <- file.path(out_dir, 'Mac_subtype') 
dir.create(Mac_out_dir)


### pca inspector ------
pdf(file.path(Mac_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(ALL_Mac, ndims = 50)
dev.off()

## dimension reduction --------
dims <- 20

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 20 * 1024 ^ 3) # 20Gb


ALL_Mac <- RunUMAP(ALL_Mac, dims = 1:dims,return.model = TRUE)

pdf(file.path(Mac_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(ALL_Mac, group.by = "sample", label = TRUE) + NoLegend()
dev.off()


ALL_Mac <- ALL_Mac %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

ALL_Mac <- RunUMAP(
  ALL_Mac, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(Mac_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(ALL_Mac, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
ALL_Mac  <-
  FindNeighbors(
    ALL_Mac ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


ALL_Mac <-
  FindClusters(
    ALL_Mac,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(Mac_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(ALL_Mac, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(ALL_Mac, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(ALL_Mac, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(ALL_Mac, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

Macrophage_major_marker_ls <- c(
  "M1"= "NOS2",
  "M1"= "IL12",
  "M1"= "FCGR1A",
  "M1"= "FCGR1B",
  "M1"= "FCGR1C",
  "M1"= "CD80",
  "M1"= "CXCR10",
  "M1"= "IL23",
  "M1"= "CXCL9",
  "M1"= "CXCL10",
  "M1"= "CXCL11",
  "M1"= "CD86",
  "M1"= "IL1A",
  "M1"= "IL1B",
  "M1"= "IL6",
  "M1"= "TNF",
  "M1"= "MHCII",
  "M1"= "CCL5",
  "M1"= "IRF5",
  "M1"= "IRF1",
  "M1"= "CD40",
  "M1"= "KYNU",
  "M1"= "CCR7",
  "M1"= "VCAN",
  "M1"= "CD300E",
  "M1"= "INHBA",
  "M1"= "TNIP3",　
  
  "M2"= "ARG1",
  "M2"= "IL10",
  "M2"= "CD32",
  "M2"= "CD163",
  "M2"= "FCER2",
  "M2"= "CD200R1",
  "M2"= "PDCD1LG2",
  "M2"= "CD274",
  "M2"= "MARCO",
  "M2"= "CSF1R",
  "M2"= "MRC1",
  "M2"= "IL1RN",
  "M2"= "Il1R2",
  "M2"= "IL4R",
  "M2"= "CCL4",
  "M2"= "CCL13",
  "M2"= "CCL20",
  "M2"= "CCL17",
  "M2"= "CCL18",
  "M2"= "CCL22",
  "M2"= "CCL24",
  "M2"= "LYVE1",
  "M2"= "VEGFA",
  "M2"= "VEGFB",
  "M2"= "VEGFC",
  "M2"= "VEGFD",
  "M2"= "EGF",
  "M2"= "CTSA",
  "M2"= "CTSB",
  "M2"= "CSTC",
  "M2"= "CTSD",
  "M2"= "TGFB1",
  "M2"= "CD209",
  "M2"= "CD163L1",
  "M2"= "FOLR2",　
  
  "M0"= "CD68",
  "M0"= "C1QA",
  "M0"= "C1QB",
  "M0"= "SELENOP",
  "M0"= "LYZ",
  "M0"= "HLA-DPB1",
  "M0"= "HLA-DPA1",
  "M0"= "AIF1", 
  
  "TAM" = "SPP1", 
  "TAM" = "CQC"
)

##
# find markers
Idents(ALL_Mac) <- "RNA_snn_res.0.4"
Mac_markers <-
  FindAllMarkers(
    JoinLayers(ALL_Mac),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(Mac_markers, file.path(Mac_out_dir, 'Mac_markers.csv'))

Mac_top5_markers <- 
  Mac_markers %>% 
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

names(Mac_top5_markers) <- sub('.*\\.','',names(Mac_top5_markers))



# dotplot
pdf(file.path(Mac_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 35, height = 6)
DotPlot(ALL_Mac, group.by = "RNA_snn_res.0.4", features = Major_marker_ls)+RotatedAxis()
DotPlot(ALL_Mac, group.by = "RNA_snn_res.0.4", features = Macrophage_major_marker_ls)+RotatedAxis()
DotPlot(ALL_Mac, group.by = "RNA_snn_res.0.4", features = Mac_top5_markers)+RotatedAxis()
dev.off()

## annotation -----------
cell_cluster <- c(
  '0' = "M1.Macro", # T, Epi
  '1' = "M2.STAB1",
  '2' = "Proliferating Macro", # 
  '3' = "Fibroblast-like Macro", # 
  '4' = "EMT", # 
  '5' = "Neutrophil", # 
  '6' = "Pericyte", # 
  '7' = "M0.Macro", # 
  '8' = "M2.MMP9", # 
  '9' = "cDC1", # 
  '10' = "DC.LAMP3", # 
  '11' = "pDC", # 
  '12' = "Unknown" ,#
  '13' = "Mast.cell"
)


Idents(ALL_Mac) <- "RNA_snn_res.0.4"

ALL_Mac <- RenameIdents(ALL_Mac, cell_cluster)
ALL_Mac$subtype <- Idents(ALL_Mac)

pdf(file.path(Mac_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(ALL_Mac, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Mac_out_dir, str_glue("Dimplot dim{dims} annotation split by group.pdf")), height = 7, width = 18)
DimPlot(ALL_Mac, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()



saveRDS(ALL_Mac, file.path(Mac_out_dir, 'ALL_Mac.rds'))

table(ALL_Mac$subtype)
####
## pure -----------
ALL_Mac_pure <- ALL_Mac[, ALL_Mac$subtype %in% c('M1.Macro', 'M2.STAB1', 'Neutrophil', 'M0.Macro','Proliferating Macro','M2.MMP9','cDC1',"DC.LAMP3","pDC","Mast.cell")]
ALL_Mac_pure <- RunUMAP(ALL_Mac_pure, dims = 1:dims, return.model = TRUE)


pdf(file.path(Mac_out_dir,str_glue("Dimplot dim{dims} annotation - pure.pdf")), height = 6, width = 8)
DimPlot(ALL_Mac_pure, group.by = "subtype", label = TRUE, repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Mac_out_dir, str_glue("Dimplot dim{dims} annotation split by group - pure.pdf")), height = 7, width = 18)
DimPlot(ALL_Mac_pure, group.by = "subtype",split.by = 'group', label = TRUE, repel = TRUE) + NoLegend()
dev.off()


cellType_marker_ls <- c(
  "M1.Macro"='CXCL9', 
  "M1.Macro"='CXCL10',
  "M1.Macro"='IRF1',
  
  "M2.STAB1"="STAB1",
  "M2.STAB1"="CD163",
  "M2.STAB1"="MRC1",
  
  "Proliferating Macro"="MKI67", # 
  "Proliferating Macro"="TOP2A",
  
  "Neutrophil"="CSF3R", # 
  
  "M0.Macro"="CD68", # 
  
  "M2.MMP9"="CSF1R", # 
  "M2.MMP9"="MMP9",
  
  "cDC1"="CLEC9A", #
  "cDC1"="IDO1",
  
  "DC.LAMP3"='CD83', #
  "DC.LAMP3"='LAMP3',
  
  "pDC"="IRF4",
  "pDC"="IRF7",
  "pDC"="IRF8",

  "Mast.cell"="MS4A2",
  "Mast.cell"="KIT",
  "Mast.cell"="HDC",
  "Mast.cell"="CTSG"
)


pdf(file.path(Mac_out_dir, "dotplot celltype main_marker - pure.pdf"), width = 20, height = 8)
DotPlot(ALL_Mac_pure, group.by = "subtype", features = cellType_marker_ls, dot.scale = 10) + theme_90angle
dev.off()

saveRDS(ALL_Mac_pure, file.path(Mac_out_dir, 'ALL_Mac_pure.rds'))
table(ALL_Mac_pure$sample)
#######
## sample -----------
meta <- ALL_Mac_pure@meta.data
meta$subtype <- as.character(meta$subtype)
meta$subtype <- as.factor(meta$subtype)
meta$sample <- factor(meta$sample, levels = c(
  'A31', 'A32', 
  'P54', 'A131', 'A71')
)
meta_data_ls <- list(meta)
names(meta_data_ls) <- 'by_sample'

meta_data_ls %>%
  iwalk(function(meta, name){
    
    Mac_out_dir_cn <- str_glue("{Mac_out_dir}/{name}/") 
    if(!file.exists(Mac_out_dir_cn)){
      dir.create(Mac_out_dir_cn)
    }
    
    anno_column <- "subtype" # <---------------------
    grp_column  <- "sample"    # <---------------------
    
    all_data_meta_roe <- meta
    
    
    ## Ro/e -------
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
    
    row_max_idx <- apply(plot.data, 1, function(x) {
      wm <- which.max(x)
      if(length(wm) == 0) return(1)
      return(as.numeric(wm[1]))
    })
    if(any(sapply(row_max_idx, is.list))) {
      row_max_idx <- sapply(row_max_idx, function(x) as.numeric(x)[1])
    }
    
    plot.data <- plot.data[order(row_max_idx, decreasing = TRUE), , drop = FALSE]# 对绘图数据行进行排序
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
    
    pdf(
      file.path(Mac_out_dir_cn, "cell props boxplot - Roe plot.pdf"), 
      width = 8, 
      height = 8
    )
    ComplexHeatmap::draw(roi_hm, annotation_legend_list = lgd)
    dev.off()
    
  })


## group -----------
meta <- ALL_Mac_pure@meta.data
meta$subtype <- as.character(meta$subtype)
meta$subtype <- as.factor(meta$subtype)
meta_data_ls <- list(meta)
names(meta_data_ls) <- 'by_group'

meta_data_ls %>%
  iwalk(function(meta, name){
    
    Mac_out_dir_cn <- str_glue("{Mac_out_dir}/{name}/") 
    if(!file.exists(Mac_out_dir_cn)){
      dir.create(Mac_out_dir_cn)
    }
    
    anno_column <- "subtype" # <---------------------
    grp_column  <- "group"    # <---------------------
    
    all_data_meta_roe <- meta
    
    
    ## Ro/e -------
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
   row_max_idx <- apply(plot.data, 1, function(x) {
      wm <- which.max(x)
      if(length(wm) == 0) return(1)
      return(as.numeric(wm[1]))
    })
    if(any(sapply(row_max_idx, is.list))) {
      row_max_idx <- sapply(row_max_idx, function(x) as.numeric(x)[1])
    }
    
    plot.data <- plot.data[order(row_max_idx, decreasing = TRUE), , drop = FALSE]# 对绘图数据行进行排序
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
    
    pdf(
      file.path(Mac_out_dir_cn, "cell props boxplot - Roe plot.pdf"), 
      width = 8, 
      height = 8
    )
    ComplexHeatmap::draw(roi_hm, annotation_legend_list = lgd)
    dev.off()
    
  })




# cell prop by annotation ----------
# 单样本进行细胞比例统计后，按分组取均值画堆叠图，百分比标在条的右侧
plot_data <- 
  ALL_Mac_pure[[]] %>% 
  ungroup() %>% 
  group_by(sample, subtype, group) %>% 
  summarise(count = n(), .groups = "drop") %>% 
  group_by(sample) %>% 
  mutate(freq = count / sum(count)) %>% 
  ungroup() %>% 
  group_by(group, subtype) %>% 
  summarise(aver_freq = mean(freq), .groups = "drop") %>% 
  group_by(subtype) %>% 
  mutate(global_mean = mean(aver_freq),
         total = sum(aver_freq)) %>% 
  ungroup() %>% 
  mutate(subtype = forcats::fct_reorder(subtype, global_mean, .desc = TRUE))

# 准备标签数据：每个 subtype 下两个组的百分比，用 '/' 拼接
label_data <- 
  plot_data %>%
  group_by(subtype) %>%
  summarise(label = paste(sprintf("%.1f%%", rev(aver_freq) * 100), collapse = "/"),
            total = sum(aver_freq))  # 总长度用于定位标签x坐标

p_col_proj <- 
  ggplot(plot_data, aes(x = aver_freq, y = subtype, fill = group)) + 
  geom_col(position = "stack", width = 0.7) +  # 堆叠条形
  geom_text(data = label_data, 
            aes(x = total + 0.02, y = subtype, label = label), 
            inherit.aes = FALSE,   # 关键：不继承主图的 aes
            hjust = 0, size = 3.5) +
  ggsci::scale_fill_igv() +
  labs(x = "cell proportion (%)", y = NULL) +
  scale_x_continuous(limits = c(0, 0.9)) +    
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# export
ggsave(
  p_col_proj,
  filename = file.path(Mac_out_dir, "cell_proportion_subtype.pdf"),
  height = 5,
  width  = 8
)
table(ALL_T$subtype)
table(all_data2_TLS$annotation)
table(ALL_Mac_pure$subtype)
###提取DC
ALL_Mac_pure$subtype <- as.character(ALL_Mac_pure$subtype)
ALL_Mac_pure$subtype <- as.factor(ALL_Mac_pure$subtype)
table(ALL_Mac_pure$subtype)

DC_data_from_Mac <- subset(ALL_Mac_pure,subtype%in%c("cDC1","DC.LAMP3","pDC"))
dim(DC_data_from_Mac)
DC_cell <- DC_data_from_Mac$cell_id_sample
all_data2_TLS$annotation[all_data2_TLS$cell_id_sample %in% DC_cell] <- "DC"
table(all_data2_TLS$annotation)
#
saveRDS(DC_data_from_Mac,"/mnt/NexData_file/home/lzl/20260625_XE2510095462-01/rds_data/DC_from_mac.rds")
saveRDS(all_data2_TLS,"/mnt/NexData_file/home/lzl/20260625_XE2510095462-01/rds_data/all_data2_TLS_add_DC.rds")

table(all_data2_TLS$annotation)
#
ALL_Mac_pure_filter <- subset(ALL_Mac_pure,subtype%in%c("M0.Macro","M1.Macro","M2.MMP9","M2.STAB1","Proliferating Macro"))
ALL_Mac_pure_filter <- RunUMAP(ALL_Mac_pure_filter, dims = 1:dims, return.model = TRUE)


pdf(file.path(Mac_out_dir,str_glue("Dimplot dim{dims} annotation - pure fil.pdf")), height = 6, width = 8)
DimPlot(ALL_Mac_pure_filter, group.by = "subtype", label = TRUE, repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Mac_out_dir, str_glue("Dimplot dim{dims} annotation split by group - pure fil.pdf")), height = 7, width = 18)
DimPlot(ALL_Mac_pure_filter, group.by = "subtype",split.by = 'group', label = TRUE, repel = TRUE) + NoLegend()
dev.off()

####DC UMAP 未聚类
DC_data_from_Mac

DC_data_from_Mac <- RunUMAP(DC_data_from_Mac, dims = 1:dims, return.model = TRUE)

pdf(file.path(Mac_out_dir,str_glue("Dimplot dim{dims} annotation - DC.pdf")), height = 6, width = 8)
DimPlot(DC_data_from_Mac, group.by = "subtype", label = TRUE, repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Mac_out_dir, str_glue("Dimplot dim{dims} annotation split by group DC.pdf")), height = 7, width = 18)
DimPlot(DC_data_from_Mac, group.by = "subtype",split.by = 'group', label = TRUE, repel = TRUE) + NoLegend()
dev.off()
