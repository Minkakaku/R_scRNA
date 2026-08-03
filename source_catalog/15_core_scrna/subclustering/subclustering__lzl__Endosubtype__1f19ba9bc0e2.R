all_data2_TLS <- readRDS("/mnt/NexData_file/home/lzl/20260625_XE2510095462-01/rds_data/all_data2_TLS_add_DC.rds")
table(all_data2_TLS$annotation)
####B cell
ALL_ENDO <- subset(all_data2_TLS, annotation %in%c('Endothelial') )


## PCA --------
DefaultAssay(ALL_ENDO) <- "RNA"
ALL_ENDO <- FindVariableFeatures(ALL_ENDO)
ALL_ENDO <- ScaleData(ALL_ENDO, vars.to.regress = 'nCount_RNA')
ALL_ENDO <- RunPCA(ALL_ENDO)

Endo_out_dir <- file.path(out_dir, 'Endo_subtype') 
dir.create(Endo_out_dir)


### pca inspector ------
pdf(file.path(Endo_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(ALL_ENDO, ndims = 50)
dev.off()

## dimension reduction --------
dims <- 15

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 20 * 1024 ^ 3) # 20Gb


ALL_ENDO <- RunUMAP(ALL_ENDO, dims = 1:dims,return.model = TRUE)

pdf(file.path(Endo_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(ALL_ENDO, group.by = "sample", label = TRUE) + NoLegend()
dev.off()


ALL_ENDO <- ALL_ENDO %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

ALL_ENDO <- RunUMAP(
  ALL_ENDO, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(Endo_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(ALL_ENDO, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
ALL_ENDO  <-
  FindNeighbors(
    ALL_ENDO ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


ALL_ENDO <-
  FindClusters(
    ALL_ENDO,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(Endo_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(ALL_ENDO, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(ALL_ENDO, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(ALL_ENDO, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(ALL_ENDO, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()



# find markers
Idents(ALL_ENDO) <- "RNA_snn_res.0.4"

major_marker_ls <- c(
  # 第0行：B cell
  "B cell" = "CD19", 
  "B cell" = "CD79A", 
  "B cell" = "MS4A1", 
  "Plasma" = "TNFRSF17", 
  "Plasma" = "JCHAIN", 
  "Plasma" = "IGJ", 
  "Plasma" = "MZB1", 
  "Plasma" = "XBP1", 
  
  "T cell" = "CD3D", 
  "T cell" = "CD3E", 
  "T cell" = "CD3G", 
  "T cell" = "TRAC", 
  "T cell" = "BCL11B",
  "T cell" = "CD4", 
  "T cell" = "FOXP3", 
  "T cell" = "CD8A", 
  "T cell" = "GZMB",
  
  "NK" = "GNLY", 
  "NK" = "NKG7", 
  "NK" = "NCAM1", 
  "NK" = "UNC13D", 
  
  # 第11行：Plasma
  
  "Adipocyte" = "PLIN1",
  "Adipocyte" = "PLIN4",
  "Adipocyte" = "ADIPOQ",
  "Adipocyte" = "FABP4",
  "Adipocyte" = "LPL",
  
  
  "Endothelial" = "PECAM1",
  "Endothelial" = "VWF", 
  "Endothelial" = "AQP1", 
  "Endothelial" = "CDH5", 
  "Endothelial" = "PLVAP", 
  
  # 第4行：Fibroblast
  "Fibroblast" = "PDGFRA", 
  "Fibroblast" = "CCDC80", 
  "Fibroblast" = "COL1A1", 
  "Fibroblast" = "COL1A2", 
  "Fibroblast" = "DCN", 
  "Fibroblast" = "ACTA2",
  
  "Pericyte" = "PDGFRB", 
  "Pericyte" = "RGS5", 
  "Pericyte" = "MCAM", 
  
  # 第3行：Macrophage
  "Macrophage" = "CD68", 
  "Macrophage" = "CD163", 
  "Macrophage" = "C1QA", 
  "Macrophage" = "C1QC", 
  "Macrophage" = "CSF1R", 
  #"Macrophage" = "IL1B", 
  
  # 第10行：Lamp3_DC
  "Lamp3_DC" = "LAMP3",
  
  # 第9行：Smooth muscle cell
  "Smooth muscle cell" = "PLN", 
  "Smooth muscle cell" = "MYLK", 
  "Smooth muscle cell" = "MYH11", 
  "Smooth muscle cell" = "ACTG2", 
  "Smooth muscle cell" = "TAGLN", 
  
  # 其他已在图中出现的细胞类型
  "Dendritic cell" = "CD83", 
  "Dendritic cell" = "LILRA4", 
  #"Dendritic cell" = "CLEC9A", 
  "Dendritic cell" = "IDO1", 
  "Dendritic cell" = "CD1C", 
  "Dendritic cell" = "FCER1A", 
  "Dendritic cell" = "CST3", 
  
  "Epithelial" = "EPCAM",
  "Epithelial" = "KRT20",
  
  "Monocyte" = "CD14",
  "Monocyte" = "FCGR3A",
  "Monocyte" = "ITGAM",
  
  "Neutrophil" = "CSF3R", 
  "Neutrophil" = "FCGR3B", 
  "Neutrophil" = "S100A8", 
  "Neutrophil" = "S100A9", 
  "Neutrophil" = "G0S2", 
  "Neutrophil" = "RETNLG", 
  
  
  "Neuron" = "TUBB3",
  "Neuron" = "SNAP25",
  
  
  "pDC" = "IRF4",
  "pDC" = "IRF7",
  "pDC" = "IRF8",
  
  "cDC1" = "XCR1",
  "cDC1" = "CLEC9A",
  
  "Proliferating" = "MKI67",
  "Proliferating" = "TOP2A",
  
  "Schwann cell" = "MPZ",
  "Schwann cell" = "SOX10",
  "Schwann cell" = "NGFR",
  
  "Mesothelial" = "MSLN",
  
  "Mast cell" = "MS4A2", 
  "Mast cell" = "KIT", 
  "Mast cell" = "HDC", 
  "Mast cell" = "TPSAB1", 
  "Mast cell" = "CTSG", 
  "Mast cell" = "CPA3",
  
  "Erythrocytes" = "HBA-A1", 
  "Erythrocytes" = "HBA-A2",
  
  'AM'="ITGAX",
  'AM'="FCGR1A",
  'AM'='SIGLEC1'
  
)
##
EC_major_marker_ls <- c(
  "LEC"= "PROM1",
  "LEC"= "PDPN",
  "LEC"= "LYVE1",
  
  "VEC"="PECAM1", 
  "VEC"="CDH5", 
  "VEC"="ENG"
)
##
Endo_markers <-
  FindAllMarkers(
    JoinLayers(ALL_ENDO),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(Endo_markers, file.path(Endo_out_dir, 'Endo_markers.csv'))

Endo_top5_markers <- 
  Endo_markers %>% 
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

names(Endo_top5_markers) <- sub('.*\\.','',names(Endo_top5_markers))



# dotplot
pdf(file.path(Endo_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 35, height = 6)
DotPlot(ALL_ENDO, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()
DotPlot(ALL_ENDO, group.by = "RNA_snn_res.0.4", features = EC_major_marker_ls)+RotatedAxis()
DotPlot(ALL_ENDO, group.by = "RNA_snn_res.0.4", features = Endo_top5_markers)+RotatedAxis()
dev.off()

## annotation -----------
cell_cluster <- c(
  '0' = "Epi-like EC", # 
  '1' = "PRCP+ VEC",
  '2' = "ESM1+ VEC", # 
  '3' = "Epi-like EC", # 
  '4' = "Pericyte-ilke EC", # 
  '5' = "Macro-like EC", # 
  '6' = "Porliferating VEC", # 
  '7' = "VEGFC+ VEC", # 
  '8' = "LEC", # 
  '9' = "Fibro-like EC"
)


Idents(ALL_ENDO) <- "RNA_snn_res.0.4"

ALL_ENDO <- RenameIdents(ALL_ENDO, cell_cluster)
ALL_ENDO$subtype <- Idents(ALL_ENDO)

pdf(file.path(Endo_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(ALL_ENDO, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Endo_out_dir, str_glue("Dimplot dim{dims} annotation split by group.pdf")), height = 7, width = 18)
DimPlot(ALL_ENDO, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()



saveRDS(ALL_ENDO, file.path(Endo_out_dir, 'ALL_ENDO.rds'))

###
## pure -----------
ALL_ENDO_pure <- ALL_ENDO[, ALL_ENDO$subtype %in% c('LEC',"PRCP+ VEC","ESM1+ VEC","Porliferating VEC","VEGFC+ VEC")]
ALL_ENDO_pure <- RunUMAP(ALL_ENDO_pure, dims = 1:dims, return.model = TRUE)


pdf(file.path(Endo_out_dir,str_glue("Dimplot dim{dims} annotation - pure.pdf")), height = 6, width = 8)
DimPlot(ALL_ENDO_pure, group.by = "subtype", label = TRUE, repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Endo_out_dir, str_glue("Dimplot dim{dims} annotation split by group - pure.pdf")), height = 7, width = 18)
DimPlot(ALL_ENDO_pure, group.by = "subtype",split.by = 'group', label = TRUE, repel = TRUE) + NoLegend()
dev.off()

####re_cluster
## PCA --------
DefaultAssay(ALL_ENDO_pure) <- "RNA"
ALL_ENDO_pure <- FindVariableFeatures(ALL_ENDO_pure)
ALL_ENDO_pure <- ScaleData(ALL_ENDO_pure, vars.to.regress = 'nCount_RNA')
ALL_ENDO_pure <- RunPCA(ALL_ENDO_pure)

#dir.create(Endo_out_dir)


### pca inspector ------
pdf(file.path(Endo_out_dir, 'pca plot_pure.pdf'), height = 6, width = 8)
ElbowPlot(ALL_ENDO_pure, ndims = 50)
dev.off()

## dimension reduction --------
dims <- 10

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 20 * 1024 ^ 3) # 20Gb


ALL_ENDO_pure <- RunUMAP(ALL_ENDO_pure, dims = 1:dims,return.model = TRUE)

pdf(file.path(Endo_out_dir, str_glue("Dimplot dim{dims} sample_pure.pdf")), height = 6, width = 8)
DimPlot(ALL_ENDO_pure, group.by = "sample", label = TRUE) + NoLegend()
dev.off()


ALL_ENDO_pure <- ALL_ENDO_pure %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

ALL_ENDO_pure <- RunUMAP(
  ALL_ENDO_pure, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(Endo_out_dir, str_glue("Dimplot dim{dims} sample_pure - harmony.pdf")), height = 6, width = 8)
DimPlot(ALL_ENDO_pure, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
ALL_ENDO_pure  <-
  FindNeighbors(
    ALL_ENDO_pure ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


ALL_ENDO_pure <-
  FindClusters(
    ALL_ENDO_pure,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(Endo_out_dir, str_glue("Dimplot dim{dims} cluster_pure.pdf")), height = 6, width = 8)
DimPlot(ALL_ENDO_pure, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(ALL_ENDO_pure, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(ALL_ENDO_pure, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(ALL_ENDO_pure, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

###
# find markers
Idents(ALL_ENDO_pure) <- "RNA_snn_res.0.2"
##
Endo_markers_pure <-
  FindAllMarkers(
    JoinLayers(ALL_ENDO_pure),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(Endo_markers_pure, file.path(Endo_out_dir, 'Endo_markers_pure.csv'))

Endo_top5_markers_pure <- 
  Endo_markers_pure %>% 
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

names(Endo_top5_markers_pure) <- sub('.*\\.','',names(Endo_top5_markers_pure))



# dotplot
pdf(file.path(Endo_out_dir, str_glue("dotplot dims{dims} major_marker_pure.pdf")), width = 35, height = 6)
DotPlot(ALL_ENDO_pure, group.by = "RNA_snn_res.0.2", features = major_marker_ls)+RotatedAxis()
DotPlot(ALL_ENDO_pure, group.by = "RNA_snn_res.0.2", features = EC_major_marker_ls)+RotatedAxis()
DotPlot(ALL_ENDO_pure, group.by = "RNA_snn_res.0.2", features = Endo_top5_markers_pure)+RotatedAxis()
dev.off()

###

cell_cluster <- c(
  '0' = "ESM1+ VEC", # T, Epi
  '1' = "PRCP+ VEC",
  '2' = "EFNB2+ VEC", # 
  '3' = "Porliferating VEC", # 
  '4' = "LEC", # 
  '5' = "Fibro-like EC"
)


Idents(ALL_ENDO_pure) <- "RNA_snn_res.0.2"

ALL_ENDO_pure <- RenameIdents(ALL_ENDO_pure, cell_cluster)
ALL_ENDO_pure$subtype <- Idents(ALL_ENDO_pure)

pdf(file.path(Endo_out_dir,str_glue("Dimplot dim{dims} annotation_pure.pdf")), height = 6, width = 8)
DimPlot(ALL_ENDO_pure, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Endo_out_dir, str_glue("Dimplot dim{dims} annotation split by group_pure.pdf")), height = 7, width = 18)
DimPlot(ALL_ENDO_pure, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

saveRDS(ALL_ENDO_pure, file.path(Endo_out_dir, 'ALL_ENDO_pure.rds'))
####
ALL_ENDO_pure_fil <- subset(ALL_ENDO_pure,subtype!="Fibro-like EC")
ALL_ENDO_pure_fil <- RunUMAP(ALL_ENDO_pure_fil, dims = 1:dims, return.model = TRUE)


pdf(file.path(Endo_out_dir,str_glue("Dimplot dim{dims} annotation - pure filt.pdf")), height = 6, width = 8)
DimPlot(ALL_ENDO_pure_fil, group.by = "subtype", label = TRUE, repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(Endo_out_dir, str_glue("Dimplot dim{dims} annotation split by group - pure filt.pdf")), height = 7, width = 18)
DimPlot(ALL_ENDO_pure_fil, group.by = "subtype",split.by = 'group', label = TRUE, repel = TRUE) + NoLegend()
dev.off()

###
saveRDS(ALL_ENDO_pure_fil, file.path(Endo_out_dir, 'ALL_ENDO_pure_fil.rds'))
###########
## sample -----------
meta <- ALL_ENDO_pure_fil@meta.data
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
    
    Endo_out_dir_cn <- str_glue("{Endo_out_dir}/{name}/") 
    if(!file.exists(Endo_out_dir_cn)){
      dir.create(Endo_out_dir_cn)
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
      file.path(Endo_out_dir_cn, "cell props boxplot - Roe plot.pdf"), 
      width = 8, 
      height = 8
    )
    ComplexHeatmap::draw(roi_hm, annotation_legend_list = lgd)
    dev.off()
    
  })


## group -----------
meta <- Endo_out_dir@meta.data
meta$subtype <- as.character(meta$subtype)
meta$subtype <- as.factor(meta$subtype)
meta_data_ls <- list(meta)
names(meta_data_ls) <- 'by_group'

meta_data_ls %>%
  iwalk(function(meta, name){
    
    Endo_out_dir_cn <- str_glue("{Endo_out_dir}/{name}/") 
    if(!file.exists(Endo_out_dir_cn)){
      dir.create(Endo_out_dir_cn)
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
      file.path(Endo_out_dir_cn, "cell props boxplot - Roe plot.pdf"), 
      width = 8, 
      height = 8
    )
    ComplexHeatmap::draw(roi_hm, annotation_legend_list = lgd)
    dev.off()
    
  })




# cell prop by annotation ----------
# 单样本进行细胞比例统计后，按分组取均值画堆叠图，百分比标在条的右侧
plot_data <- 
  ALL_ENDO_pure_fil[[]] %>% 
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
  filename = file.path(Endo_out_dir, "cell_proportion_subtype.pdf"),
  height = 5,
  width  = 8
)
