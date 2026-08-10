####B cell
ALL_B <- subset(all_data2_TLS, annotation %in%c('B cell',"Plasma") )


## PCA --------
DefaultAssay(ALL_B) <- "RNA"
ALL_B <- FindVariableFeatures(ALL_B)
ALL_B <- ScaleData(ALL_B, vars.to.regress = 'nCount_RNA')
ALL_B <- RunPCA(ALL_B)

B_out_dir <- file.path(out_dir, 'B_subtype') 
dir.create(B_out_dir)


### pca inspector ------
pdf(file.path(B_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(ALL_B, ndims = 50)
dev.off()

## dimension reduction --------
dims <- 15

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 20 * 1024 ^ 3) # 20Gb


ALL_B <- RunUMAP(ALL_B, dims = 1:dims,return.model = TRUE)

pdf(file.path(B_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(ALL_B, group.by = "sample", label = TRUE) + NoLegend()
dev.off()


ALL_B <- ALL_B %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

ALL_B <- RunUMAP(
  ALL_B, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(B_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(ALL_B, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
ALL_B  <-
  FindNeighbors(
    ALL_B ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


ALL_B <-
  FindClusters(
    ALL_B,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(B_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(ALL_B, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(ALL_B, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(ALL_B, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(ALL_B, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()



# find markers
Idents(ALL_B) <- "RNA_snn_res.0.6"
##
B_major_marker_ls <- c(
  "B cell" = "CD19",
  "B cell" = "CD79A",
  "B cell" = "MS4A1", 
  
  "naive B" = "CD20",
  "naive B" = "CD24",
  "naive B" = "IGHM",
  "naive B" = "BANK1",
  "naive B" = "FCRL1",
  "naive B" = "HLA-DQB1",
  "naive B" = "HLA-DQA1",
  "naive B" = "HVCN1",
  "naive B" = "IGHD",
  "naive B" = "FCER2",
  "naive B" = "IL4R", 
  
  "Follicular B" = "CD40", 
  "Follicular B" = "BCL6", 
  
  "memory B"= "CD27",
  "memory B"= "FCRL4",
  "memory B"= "TACI",
  "memory B"= "AIM2",
  "memory B"= "TXNIP",
  "memory B"= "GPR183",
  "memory B"= "TNFRSF13B",
  "memory B"= "CD22",
  "memory B"= "FCRL2",
  "memory B"= "ADAM28",
  "memory B"= "IGLL5",
  "memory B"= "P2RX5",
  "memory B"= "BLK", 
  
  # "Marginal zone B" = "CD27", 
  "Marginal zone B" = "CR2",
  "Marginal zone B" = "CD1C", 
  
  "cycling B"="MKI67",
  "cycling B"="STMN1",
  "cycling B"="HMGB2",
  "cycling B"="BIRC5", 
  
  "Plasma" = "CD138",
  "Plasma" = "PRDM1",
  "Plasma" = "BLIMP1",
  "Plasma" = "IGHA1",
  "Plasma" = "IGHA2",
  "Plasma" = "IGHG1",
  "Plasma" = "IGHG2",
  "Plasma" = "IGHG3",
  "Plasma" = "IGHG4",
  "Plasma" = "IGKC",
  "Plasma" = "TNFRSF17",
  "Plasma" = "JCHAIN",
  "Plasma" = "MZB1",
  "Plasma" = "XBP1",
  "Plasma" = "SDC1", 
  
  #"Plasmablast" = "MKI67", 
  "Plasmablast" = "CD38", 
  "Plasmablast" = "IRF4", 
  
  "Breg" = "IL10", 
  "Breg" = "CD1D", 
  "Breg" = "TGFB1", 
  
  "GC B"= "S1PI2",
  # "GC B"="BCL6",
  "GC B"="RGS13",
  "GC B"="TCL1A",
  "GC B"= "LRMP",
  "GC B"= "SUGCT",
  "GC B"= "MME",
  "GC B"="AICDA"
)
B_markers <-
  FindAllMarkers(
    JoinLayers(ALL_B),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(B_markers, file.path(B_out_dir, 'B_markers.csv'))

B_top5_markers <- 
  B_markers %>% 
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

names(B_top5_markers) <- sub('.*\\.','',names(B_top5_markers))



# dotplot
pdf(file.path(B_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 35, height = 6)
DotPlot(ALL_B, group.by = "RNA_snn_res.0.6", features = Major_marker_ls)+RotatedAxis()
DotPlot(ALL_B, group.by = "RNA_snn_res.0.6", features = B_major_marker_ls)+RotatedAxis()
DotPlot(ALL_B, group.by = "RNA_snn_res.0.6", features = B_top5_markers)+RotatedAxis()
dev.off()

## annotation -----------
cell_cluster <- c(
  '0' = "Plasma", # T, Epi
  '1' = "naive B",
  '2' = "Fibro-like plasma", # 
  '3' = "Plasmablast", # 
  '4' = "GC B/cycling B", # 
  '5' = "CCND1+ B", # 
  '6' = "T-like B", # 
  '7' = "Marginal zone B", # 
  '8' = "Fibro-like plasma", # 
  '9' = "Fibro-like plasma", # 
  '10' = "Macro-like B",
  '11' = "T-like B",
  '12' = "Endothelial",
  '13' = "DC",
  '14' = "Epithelial",
  '15' = "B-plasma mix"
)


Idents(ALL_B) <- "RNA_snn_res.0.6"

ALL_B <- RenameIdents(ALL_B, cell_cluster)
ALL_B$subtype <- Idents(ALL_B)

pdf(file.path(B_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(ALL_B, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(B_out_dir, str_glue("Dimplot dim{dims} annotation split by group.pdf")), height = 7, width = 18)
DimPlot(ALL_B, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()



saveRDS(ALL_B, file.path(B_out_dir, 'ALL_B.rds'))

###
## pure -----------
ALL_B_pure <- ALL_B[, ALL_B$subtype %in% c('Plasma', 'naive B', 'Plasmablast', 'GC B/cycling B','CCND1+ B','Marginal zone B')]
ALL_B_pure <- RunUMAP(ALL_B_pure, dims = 1:dims, return.model = TRUE)


pdf(file.path(B_out_dir,str_glue("Dimplot dim{dims} annotation - pure.pdf")), height = 6, width = 8)
DimPlot(ALL_B_pure, group.by = "subtype", label = TRUE, repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(B_out_dir, str_glue("Dimplot dim{dims} annotation split by group - pure.pdf")), height = 7, width = 18)
DimPlot(ALL_B_pure, group.by = "subtype",split.by = 'group', label = TRUE, repel = TRUE) + NoLegend()
dev.off()

####re_cluster
## PCA --------
DefaultAssay(ALL_B_pure_pure) <- "RNA"
ALL_B_pure <- FindVariableFeatures(ALL_B_pure)
ALL_B_pure <- ScaleData(ALL_B_pure, vars.to.regress = 'nCount_RNA')
ALL_B_pure <- RunPCA(ALL_B_pure)

B_out_dir <- file.path(out_dir, 'B_subtype') 
dir.create(B_out_dir)


### pca inspector ------
pdf(file.path(B_out_dir, 'pca plot_pure.pdf'), height = 6, width = 8)
ElbowPlot(ALL_B_pure, ndims = 50)
dev.off()

## dimension reduction --------
dims <- 10

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 20 * 1024 ^ 3) # 20Gb


ALL_B_pure <- RunUMAP(ALL_B_pure, dims = 1:dims,return.model = TRUE)

pdf(file.path(B_out_dir, str_glue("Dimplot dim{dims} sample_pure.pdf")), height = 6, width = 8)
DimPlot(ALL_B_pure, group.by = "sample", label = TRUE) + NoLegend()
dev.off()


ALL_B_pure <- ALL_B_pure %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

ALL_B_pure <- RunUMAP(
  ALL_B_pure, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(B_out_dir, str_glue("Dimplot dim{dims} sample_pure - harmony.pdf")), height = 6, width = 8)
DimPlot(ALL_B_pure, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
ALL_B_pure  <-
  FindNeighbors(
    ALL_B_pure ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


ALL_B_pure <-
  FindClusters(
    ALL_B_pure,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(B_out_dir, str_glue("Dimplot dim{dims} cluster_pure.pdf")), height = 6, width = 8)
DimPlot(ALL_B_pure, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(ALL_B_pure, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(ALL_B_pure, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(ALL_B_pure, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

###
# find markers
Idents(ALL_B_pure) <- "RNA_snn_res.0.6"
##
B_markers_pure <-
  FindAllMarkers(
    JoinLayers(ALL_B_pure),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(B_markers_pure, file.path(B_out_dir, 'B_markers_pure.csv'))

B_top5_markers_pure <- 
  B_markers_pure %>% 
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

names(B_top5_markers_pure) <- sub('.*\\.','',names(B_top5_markers_pure))



# dotplot
pdf(file.path(B_out_dir, str_glue("dotplot dims{dims} major_marker_pure.pdf")), width = 35, height = 6)
DotPlot(ALL_B_pure, group.by = "RNA_snn_res.0.6", features = Major_marker_ls)+RotatedAxis()
DotPlot(ALL_B_pure, group.by = "RNA_snn_res.0.6", features = B_major_marker_ls)+RotatedAxis()
DotPlot(ALL_B_pure, group.by = "RNA_snn_res.0.6", features = B_top5_markers_pure)+RotatedAxis()
dev.off()

###
cell_cluster <- c(
  '0' = "Plasma", # T, Epi
  '1' = "Naive B/Memory B",
  '2' = "Plasma", # 
  '3' = "Macro-like B", # 
  '4' = "GC B", # 
  '5' = "Marginal zone B", # 
  '6' = "Marginal zone B",# 
  '7' = "Macro-like B",
  '8' = "Plasmablast",
  '9' = "Follicular B",
  '10' = "Cycling B"
)


Idents(ALL_B_pure) <- "RNA_snn_res.0.6"

ALL_B_pure <- RenameIdents(ALL_B_pure, cell_cluster)
ALL_B_pure$subtype <- Idents(ALL_B_pure)

pdf(file.path(B_out_dir,str_glue("Dimplot dim{dims} annotation_pure.pdf")), height = 6, width = 8)
DimPlot(ALL_B_pure, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(B_out_dir, str_glue("Dimplot dim{dims} annotation split by group_pure.pdf")), height = 7, width = 18)
DimPlot(ALL_B_pure, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

saveRDS(ALL_B_pure, file.path(B_out_dir, 'ALL_B_pure.rds'))
####
ALL_B_pure_fil <- subset(ALL_B_pure,subtype!="Macro-like B")
ALL_B_pure_fil <- RunUMAP(ALL_B_pure_fil, dims = 1:dims, return.model = TRUE)


pdf(file.path(B_out_dir,str_glue("Dimplot dim{dims} annotation - pure filt.pdf")), height = 6, width = 8)
DimPlot(ALL_B_pure_fil, group.by = "subtype", label = TRUE, repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(B_out_dir, str_glue("Dimplot dim{dims} annotation split by group - pure filt.pdf")), height = 7, width = 18)
DimPlot(ALL_B_pure_fil, group.by = "subtype",split.by = 'group', label = TRUE, repel = TRUE) + NoLegend()
dev.off()

###
saveRDS(ALL_B_pure_fil, file.path(B_out_dir, 'ALL_B_pure_fil.rds'))

###########
## sample -----------
meta <- ALL_B_pure_fil@meta.data
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
    
    B_out_dir_cn <- str_glue("{B_out_dir}/{name}/") 
    if(!file.exists(B_out_dir_cn)){
      dir.create(B_out_dir_cn)
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
      file.path(B_out_dir_cn, "cell props boxplot - Roe plot.pdf"), 
      width = 8, 
      height = 8
    )
    ComplexHeatmap::draw(roi_hm, annotation_legend_list = lgd)
    dev.off()
    
  })


## group -----------
meta <- B_out_dir@meta.data
meta$subtype <- as.character(meta$subtype)
meta$subtype <- as.factor(meta$subtype)
meta_data_ls <- list(meta)
names(meta_data_ls) <- 'by_group'

meta_data_ls %>%
  iwalk(function(meta, name){
    
    B_out_dir_cn <- str_glue("{B_out_dir}/{name}/") 
    if(!file.exists(B_out_dir_cn)){
      dir.create(B_out_dir_cn)
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
      file.path(B_out_dir_cn, "cell props boxplot - Roe plot.pdf"), 
      width = 8, 
      height = 8
    )
    ComplexHeatmap::draw(roi_hm, annotation_legend_list = lgd)
    dev.off()
    
  })




# cell prop by annotation ----------
# 单样本进行细胞比例统计后，按分组取均值画堆叠图，百分比标在条的右侧
plot_data <- 
  ALL_B_pure_fil[[]] %>% 
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
  filename = file.path(B_out_dir, "cell_proportion_subtype.pdf"),
  height = 5,
  width  = 8
)
