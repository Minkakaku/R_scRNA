setwd("/mnt/NexData_file/home/lzl/20260508_XE2511106427-01/")
library(tidyverse)
library(Seurat)
library(ragg)
# set default ---------
#out_dir <- 'KC_subtype'
dir.create(out_dir)

# data loading ---------
all_data <- readRDS('/mnt/NexData_file/home/lzl/20260508_XE2511106427-01/new_celltype/alldata_celltype.rds')
out_dir <- "KC_subtype_1"  
dir.create(out_dir)
KC_data <- subset(all_data, annotation %in% c("Bas KC","Grn KC","Spn KC","Ecc Duct","HFE","Seb","Melano","Ecc Gland"))

KC_data <- JoinLayers(KC_data)
## PCA --------
KC_data <- FindVariableFeatures(KC_data)
length(VariableFeatures(KC_data)) #2000
KC_data <- ScaleData(KC_data)
KC_data <- RunPCA(KC_data)

### pca inspector ------
pdf(file.path(out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(KC_data, ndims = 50)
dev.off()

pdf(file.path(out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(KC_data, dims = 1:5)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(KC_data, dims = 6:10)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(KC_data, dims = 11:15)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(KC_data, dims = 16:20)
dev.off()

## dimension reduction --------
dims <- 20 # 18

library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

KC_data <- RunUMAP(KC_data, dims = 1:dims,return.model = TRUE)

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 8, width = 8)
DimPlot(KC_data, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

KC_data <- KC_data %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

KC_data <- RunUMAP(
  KC_data, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 8, width = 8)
DimPlot(KC_data, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
KC_data  <-
  FindNeighbors(
    KC_data ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


KC_data <-
  FindClusters(
    KC_data,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 6)
DimPlot(KC_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(KC_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(KC_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(KC_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

Idents(KC_data) <- "RNA_snn_res.0.6"

### dotplot -------
KC_major_marker_ls <- c(
  # 角质形成细胞亚型
  ###基底层1基底干细胞/增殖态
  KRT_B1 = 'KRT5',
  KRT_B1 = 'KRT14',
  KRT_B1 = 'CYR61',
  KRT_B1 = 'ADRB2',
  
  ###基底层2基底细胞
  KRT_B2 = 'KRT5',
  KRT_B2 = 'KRT14',
  
  ###棘层分化细胞
  KRT_SP = 'KRT1',
  KRT_SP = 'KRT10',
  
  ###颗粒层
  KRT_GR = 'KRT2',
  KRT_GR = 'KRT11',
  
  ###汗腺
  KRT_ECR = 'CFTR',
  KRT_ECR = 'MUCL1'
)


KC_major_marker_ls <- c(
  # 角质形成细胞
  Keratinocyte = "KRT14",
  Keratinocyte = "KRT1",
  Keratinocyte = "DMKN",
  Keratinocyte = "KRT10",
  Keratinocyte = "KRT5",
  Keratinocyte = "KDRTAP",
  
  # 黑素细胞
  Melanocyte = "DCT",
  Melanocyte = "TYRP1",
  Melanocyte = "PMEL",
  Melanocyte = "MLANA",
  Melanocyte = "QPCT",
  Melanocyte = "MITF",
  Melanocyte = "MART1",
  Melanocyte = "SOX10",
  Melanocyte = "PAX3",
  Melanocyte = "TYR",
  Melanocyte = "TYRP2",
  
  # 外泌汗腺细胞
  Eccrine_gland_cell = "SCGB1D2",
  Eccrine_gland_cell = "PIP",
  Eccrine_gland_cell = "SCGB1B2P",
  Eccrine_gland_cell = "SCGB2A2",
  Eccrine_gland_cell = "MUCL1",
  Eccrine_gland_cell = "DCD",
  
  # 神经细胞
  Nerve_cell = "MPZ",
  Nerve_cell = "PLP1",
  Nerve_cell = "S100B",
  Nerve_cell = "PMP22",
  Nerve_cell = "CRYAB",
  Nerve_cell = "MBP",
  
  # 毛囊干细胞
  Hair_follicle_stem_cell = "LGR5",
  Hair_follicle_stem_cell = "LGR6",
  Hair_follicle_stem_cell = "CD200",
  Hair_follicle_stem_cell = "SCD1",
  Hair_follicle_stem_cell = "ELOVL4",
  Hair_follicle_stem_cell = "FA2H",
  
  # 内根鞘细胞
  Inner_root_sheath_cell = "KRT71",
  Inner_root_sheath_cell = "KRT25",
  
  # 外根鞘细胞
  Out_root_sheath_cell = "KRT17",
  
  # 伴随层细胞
  Companion_layer_cell = "KRT6A",
  Companion_layer_cell = "KRT14",
  
  # 脂肪细胞
  Adipocyte = "PLIN1",
  Adipocyte = "PLIN4",
  
  # 基底层细胞
  Stratum_basal_cell = "COL17A1",
  Stratum_basal_cell = "ITGB4",
  
  # 疾病相关细胞
  disease_associated = "DSG1",
  disease_associated = "LAD1",
  
  # 循环过渡扩增细胞
  Cycling_transit_Amplifying_cell = "MKI67",
  Cycling_transit_Amplifying_cell = "TOP2A",
  Cycling_transit_Amplifying_cell = "PCNA",
  
  # 皮脂腺细胞
  Sebocyte = "FASN",
  Sebocyte = "PPARG"
)



KC_major_marker_ls <- c(
  # 角质形成细胞
  `Cyc KC` = "MKI67",
  `Cyc KC` = "TOP2A",
  
  `Bas KC` = "COL17A1",
  `Bas KC` = "COL7A1",
  
  `Spn KC` = "RHOV",
  `Spn KC` = "KLF5",
  
  `Grn KC` = "DSC1",
  `Grn KC` = "FLG",
  
  HFE = "GJB6",
  HFE = "SOX9",
  
  `Ecc Duct` = "KRT77",
  `Ecc Duct` = "SEMA3C",
  
  `Ecc Gland` = "SFRP1",
  `Ecc Gland` = "DNER"
)

KC_major_marker_ls <- list(
  `Cyc KC`          = c("TOP2A", "MKI67", "CCNB1", "CDC20", "BUB1"),
  `Bas KC`          = c("COL17A1", "WNT4", "ITGA2", "COL7A1", "ITGB4"),
  `Spn KC I`        = c("DSC1", "GATA3", "WNT7B", "FGFR3", "IL1RN"),
  `Spn KC II`       = c("SERPINB4", "S100A8", "GJB2", "C1orf54", "WNT5A"),
  `Grn KC`          = c("FLG", "IL37", "EREG", "KRT77", "CD27"),
  `Ecc Duct`        = c("LAMB3", "KRT77", "RSPO1", "SEMA3C", "INHBA"),
  `Ecc Gland`       = c("DNER", "CLDN10", "SCGB1B2P", "LTF", "EPCAM"),
  `Ecc Myoepi`      = c("SCGB1B2P", "SFRP1", "DNER", "CCN2", "SOX10"),
  `Bas Inf`         = c("ALDH3A1", "PTN", "PTHLH", "MIR205HG"),
  Melano            = c("MLANA", "MITF", "KIT", "IRF4", "EDNRB"),
  `Diff Inf`        = c("IGFL2", "KLK7", "PTN", "IL34"),
  `Diff Seb`        = c("KRT79", "MGST1", "APOC1", "CCL17", "KRT4"),
  `Bas Seb`         = c("IL1R2", "WFDC2", "SCARB1", "FERMT1", "TKT"),
  Bulge             = c("IL31RA", "WIF1", "DKK3", "DIO2"),
  `ORS Basal`       = c("LGR5", "IL31RA", "COL11A1", "EDN2", "GLI1"),
  `ORS Suprabasal`  = c("GJB6", "KRT75", "SOX9", "FZD7", "LGR4"),
  `HS`          = c( "S100A3")
)
KC_major_marker_ls <- lapply(KC_major_marker_ls, function(genes) {
  genes[genes %in% rownames(KC_data)]
})


KC_major_marker_ls <- KC_major_marker_ls[KC_major_marker_ls%in%rownames(KC_data)]

#KC_major_marker_ls <- c(
#  KC_basal = 'COL17A1',
#  KC_basal = 'POSTN',
#  KC_basal = 'ITGB4',
#  KC_spinous = 'KRT1',
#  KC_granular = 'KRT2',
#  KC_activated = 'KRT17',
#  KC_eccrine = 'CFTR'
#)



Idents(KC_data)<-'RNA_snn_res.0.6'
KC_data_markers <-
  FindAllMarkers(
    JoinLayers(KC_data),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )

KC_data_top5_markers <- 
  KC_data_markers %>% 
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

names(KC_data_top5_markers) <- sub('.*\\.','',names(KC_data_top5_markers))

# dotplot
pdf(file.path(out_dir, str_glue("dotplot dims{dims} major_marker_new.pdf")), width = 36, height = 6)
#table(unlist(KC_major_marker_ls))[table(unlist(KC_major_marker_ls)) > 1]
DotPlot(KC_data, group.by = "RNA_snn_res.0.6", features = KC_major_marker_ls)+RotatedAxis()
dev.off()


pdf(file.path(out_dir, str_glue("dotplot dims{dims} top5_marker.pdf")), width = 30, height = 6)

DotPlot(KC_data, group.by = "RNA_snn_res.0.6", features = KC_data_top5_markers)+RotatedAxis()
dev.off()

"POSTN⁻"
KC_data$subtype <- case_when(
  KC_data$RNA_snn_res.0.6 == 0 ~ "Spn KC I",
  KC_data$RNA_snn_res.0.6 == 1 ~ "Bas KC",
  KC_data$RNA_snn_res.0.6 == 2 ~ "Seb",
  KC_data$RNA_snn_res.0.6 == 3 ~ "Seb",
  KC_data$RNA_snn_res.0.6 == 4 ~ "Ecc duct",
  KC_data$RNA_snn_res.0.6 == 5 ~ "ORS basal",
  KC_data$RNA_snn_res.0.6 == 6 ~ "ORS suprabasal",
  KC_data$RNA_snn_res.0.6 == 7 ~ "Grn KC",
  KC_data$RNA_snn_res.0.6 == 8 ~ "Cyc KC",
  KC_data$RNA_snn_res.0.6 == 9 ~ "Ecc gland",
  KC_data$RNA_snn_res.0.6 == 10 ~ "Ecc Myoepi",
  KC_data$RNA_snn_res.0.6 == 11 ~ "Spn KC II",
  KC_data$RNA_snn_res.0.6 == 12 ~ "Others",
  KC_data$RNA_snn_res.0.6 == 13 ~ "Seb",
  KC_data$RNA_snn_res.0.6 == 14 ~ "Spn KC I",
  KC_data$RNA_snn_res.0.6 == 15 ~ "Melano",
  KC_data$RNA_snn_res.0.6 == 16 ~ "Spn KC II",
  KC_data$RNA_snn_res.0.6 == 17 ~ "IRS/HS",
  KC_data$RNA_snn_res.0.6 == 18 ~ "Inf",
  TRUE ~ "Unknown"  # 其他情况
)
epi_colors <- c(
  "Spn KC I"         = "#8B4513",  # SaddleBrown
  "Bas KC"           = "#CD5C5C",  # IndianRed
  "Grn KC"           = "#D2B48C",  # Tan
  "Spn KC II"        = "#2E8B57",  # SeaGreen
  "Cyc KC"           = "#F08080",  # LightCoral
  #"Bas Inf"          = "#FF7F50",  # Coral
  "Inf"         = "#FFD700",  # Gold
  #"Diff Seb"         = "#3CB371",  # MediumSeaGreen
  "Seb"          = "#90EE90",  # LightGreen
  #"Bulge"            = "#9ACD32", # YellowGreen
  "Melano"           = "#7CFC00",  # LawnGreen
  "Ecc duct"         = "#00FA9A",  # MediumSpringGreen
  "Ecc gland"        = "#00FF7F",  # SpringGreen
  "Ecc Myoepi"       = "#00CED1",  # DarkTurquoise
  "ORS basal"        = "#4682B4",  # SteelBlue
  "ORS suprabasal"   = "#5F9EA0",  # CadetBlue
  "IRS/HS"               = "#48D1CC",  # MediumTurquoise
  "Others"           = "#808080"
)
out_dir <- "KC_subtype_1"
pdf(file.path(out_dir, str_glue("Dimplot dim{dims} subtype.pdf")), height = 6, width = 6)
DimPlot(KC_data, group.by = "subtype", label = TRUE, reduction = 'harmonyumap') +
  scale_color_manual(values = epi_colors) +
  NoLegend()
dev.off()

out_dir <- "KC_subtype_1"
DotPlot(KC_data, group.by = "RNA_snn_res.0.6", features = KC_major_marker_ls) +
  RotatedAxis() +
  theme(
    strip.text = element_text(
      angle = 90,      # 垂直旋转
      hjust = 0,       # 左对齐（0=左，0.5=中，1=右）
      vjust = 0.5,     # 垂直居中
      size = 14
    ),
    strip.placement = "outside",
    panel.spacing.x = unit(0.2, "lines"),  # 减小面板间距，让标签更紧凑
    plot.margin = margin(t = 50)           # 顶部留空间
  ) +
  coord_cartesian(clip = "off") 
pdf(file.path(out_dir, str_glue("nCount_RNA.pdf")), height = 6, width = 6)
VlnPlot(KC_data, features = "nCount_RNA", pt.size = 0) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

FeaturePlot(KC_data, features = "KRT1",reduction = 'harmonyumap')


########ALL_CELL
all_data$subtype <- as.character(all_data$annotation)  # 或 NA

# 提取 KC_data 的细胞名和对应注释
kc_cells <- colnames(KC_data)
kc_subtype <- KC_data$subtype  # 如 KC_data$cell_type 或 KC_data$RNA_snn_res.0.4

unique(KC_data$subtype)
# 创建命名向量
subtype_map <- setNames(as.character(kc_subtype), kc_cells)

# 赋值到 all_data
all_data$subtype[match(kc_cells, colnames(all_data))] <- subtype_map[kc_cells]

# 非 KC 细胞保留原注释
non_kc <- setdiff(colnames(all_data), kc_cells)
all_data$subtype[non_kc] <- as.character(all_data$annotation[non_kc])
unique(all_data$subtype)

rename_vector <- c(
  "Fibroblast" = "Fibro",        
  "Mast cell" = "Mast",          
  "Dendritic cell" = "DC",       
  "T cell" = "T cell",           # 保持不变或按需改为 Tc
  "Macrophage" = "Mac",          
  "HFE" = "HFE",                 # 暂未提供简称，保持原样
  "Out root sheath cell" = "ORS",
  "Melanocyte" = "Melano",   
  "Endothelial" = "EC",          
  "Sebocyte" = "Seb",       
  "Pericyte" = "Peri",           
  "Schwann cell" = "Schwann",    
  "Adipocyte" = "Adipo",         
  "Lymphatic endothelial cell" = "LEC"
)
rename_vector <- c(
  "Ecc" = "Ecc Duct"        
)
all_data$subtype <- dplyr::recode(all_data$subtype, !!!rename_vector)

out_dir <- "Celltype_spatial"
pdf(file.path(out_dir, str_glue("Dimplot dim{dims} alltype.pdf")), height = 6, width = 6)
DimPlot(all_data_1, group.by = "annotation", label = TRUE) +
  scale_color_manual(values = cell_colors) +
  NoLegend()
dev.off()
###HRAS突变基因
out_dir <- "HRAS_mut"
dir.create(out_dir)
grep("HRAS", rownames(all_data), value = TRUE)
c("HRAS","HRAS-ALT-37:C","HRAS-WT-37:C")


KC_data$HRAS_status <- ifelse(
  FetchData(KC_data, vars = "HRAS-ALT-37:C")[, 1] > 0,
  "HRAS_ALT",
  "other"
)

FeaturePlot(KC_data, 
            features = "HRAS-ALT-37:C",
            reduction = 'harmonyumap',
            order = TRUE,                       # 高表达点画在最上层
            min.cutoff = 0,                     # 最小值截断
            max.cutoff = "q95")

hras_colors <- c("other" = "grey80", "HRAS_ALT" = "red")

plot_data <- FetchData(KC_data, vars = c("harmonyumap_1", "harmonyumap_2", "HRAS_status"))




pdf(file.path(out_dir, str_glue("dimplot dims{dims} HRAS.pdf")), width = 6, height = 6)
FeaturePlot(KC_data, 
            features = "HRAS",
            reduction = 'harmonyumap',
            order = TRUE,                       # 高表达点画在最上层
            min.cutoff = 0,                     # 最小值截断
            max.cutoff = "q95")
dev.off()

hras_alt_expr <- FetchData(KC_data, vars = "HRAS")[, 1]

# 判定表达（> 0 即认为表达）
is_expressed <- hras_alt_expr > 0

# 计算比例
total_cells <- ncol(KC_data)
expressed_cells <- sum(is_expressed)
proportion_expressed <- expressed_cells / total_cells

cat("表达 HRAS 的细胞数:", expressed_cells, "\n")
cat("总细胞数:", total_cells, "\n")
cat("表达比例:", round(proportion_expressed * 100, 2), "%\n")


# 提取表达细胞的亚群信息
expressed_subtype <- KC_data$subtype[is_expressed]

# 计算各亚群比例
subtype_counts <- table(expressed_subtype)
subtype_proportions <- prop.table(subtype_counts)

# 整理为数据框
result_df <- data.frame(
  Subtype = names(subtype_counts),
  Count = as.numeric(subtype_counts),
  Proportion = as.numeric(subtype_proportions),
  Percentage = round(as.numeric(subtype_proportions) * 100, 2)
)

# 按数量排序
result_df <- result_df[order(result_df$Count, decreasing = TRUE), ]
rownames(result_df) <- NULL

print(result_df)
####Immnue_subtype
out_dir <- "Immnue_subtype"  
dir.create(out_dir)
unique(all_data$annotation)
Immnue_data <- subset(all_data, annotation %in% c("Lym","Mac","DC","Mast"))

Immnue_data <- JoinLayers(Immnue_data)
## PCA --------
Immnue_data <- FindVariableFeatures(Immnue_data)
length(VariableFeatures(Immnue_data)) #2000
Immnue_data <- ScaleData(Immnue_data)
Immnue_data <- RunPCA(Immnue_data)

### pca inspector ------
pdf(file.path(out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Immnue_data, ndims = 50)
dev.off()

pdf(file.path(out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(Immnue_data, dims = 1:5)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(Immnue_data, dims = 6:10)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(Immnue_data, dims = 11:15)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(Immnue_data, dims = 16:20)
dev.off()

## dimension reduction --------
dims <- 25 # 18

library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

Immnue_data <- RunUMAP(Immnue_data, dims = 1:dims,return.model = TRUE)

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 8, width = 8)
DimPlot(Immnue_data, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

Immnue_data <- Immnue_data %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Immnue_data <- RunUMAP(
  Immnue_data, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 8, width = 8)
DimPlot(Immnue_data, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Immnue_data  <-
  FindNeighbors(
    Immnue_data ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Immnue_data <-
  FindClusters(
    Immnue_data,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 6)
DimPlot(Immnue_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Immnue_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Immnue_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Immnue_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

Idents(Immnue_data) <- "RNA_snn_res.0.6"

### dotplot -------

Immnue_major_marker_ls <- Immnue_major_marker_ls[Immnue_major_marker_ls%in%rownames(Immnue_data)]
Immnue_major_marker_ls <- list(
  `CCR7+ DC`       = c("CCR7", "CD83", "IL7R", "CXCR4", "BIRC3"),
  `CD1C+ DC`       = c("CLEC10A", "ITGAX", "IL1B", "HAVCR2", "MMP9"),
  `CD4+ Th`        = c("IL7R", "CCR4", "CD40LG", "CD3E", "CD28"),
  `CD4+ Treg`      = c("FOXP3", "CTLA4", "IL2RA", "TIGIT", "CD27"),
  `CD8+ Tc`        = c("GZMK", "CCL5", "CD8A", "NKG7", "GZMA"),
  `CLEC9A+ DC`     = c("IRF8", "CLEC9A", "FLT3", "IDO1", "BATF3"),
  `Cyc Imm`        = c("FGF18", "MKI67", "CRABP1", "JCHAIN", "JAK3"),
  `LC`             = c("CD207", "CD1A", "CSF1R", "TNF", "CD83"),
  `Mac`            = c("CD163", "C1QB", "C1QA", "F13A1", "MRC1"),
  `Mast`           = c("HDC", "CPA3", "CMA1", "IL1RL1", "KIT"),
  `Mono`           = c("OSM", "ITGAX", "CXCR1", "CXCR2", "IL1B"),
  `Naive T`        = c("CCR4", "IL7R", "CD3E", "CD3D", "CD28"),
  `Plasma`         = c("JCHAIN", "KLRG1", "CCL22", "IDO1", "TCL1A")
)
Immnue_major_marker_ls <- list(
  `CCR7+ DC`       = c("CCR7", "CD83", "IL7R", "CXCR4", "BIRC3"), # 首次出现，全部保留
  `CD1C+ DC`       = c("CLEC10A", "HAVCR2", "MMP9"),              # 去重 ITGAX, IL1B
  `CD4+ Th`        = c("CD40LG"),                                 # 去重 IL7R, CCR4, CD3E, CD28
  `CD4+ Treg`      = c("FOXP3", "CTLA4", "IL2RA", "TIGIT", "CD27"),
  `CD8+ Tc`        = c("CD8A", "GZMA", "GZMK", "NKG7", "CCL5"),
  `CLEC9A+ DC`     = c("CLEC9A", "IRF8", "BATF3", "FLT3"),        # 去重 IDO1
  `Cyc Imm`        = c("MKI67", "CRABP1", "JAK3"),                # 去重 JCHAIN
  `LC`             = c("CD207", "CD1A", "CSF1R", "TNF"),          # 去重 CD83
  `Mac`            = c("CD163", "C1QA", "C1QB", "F13A1", "MRC1"),
  `Mast`           = c("HDC", "CPA3", "CMA1", "IL1RL1", "KIT"),
  `Mono`           = c("OSM", "CXCR1", "CXCR2"),                  # 去重 ITGAX, IL1B
  `Naive T`        = c("CD3D"),                                   # 去重 IL7R, CCR4, CD3E, CD28
  `Plasma`         = c("TCL1A", "KLRG1", "CCL22")                 # 去重 JCHAIN, IDO1
)
Immnue_major_marker_ls <- lapply(Immnue_major_marker_ls, function(genes) {
  genes[genes %in% rownames(KC_data)]
})

#KC_major_marker_ls <- c(
#  KC_basal = 'COL17A1',
#  KC_basal = 'POSTN',
#  KC_basal = 'ITGB4',
#  KC_spinous = 'KRT1',
#  KC_granular = 'KRT2',
#  KC_activated = 'KRT17',
#  KC_eccrine = 'CFTR'
#)


Idents(Immnue_data)<-'RNA_snn_res.0.6'
Immnue_data_markers <-
  FindAllMarkers(
    JoinLayers(Immnue_data),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )

Immnue_data_top5_markers <- 
  Immnue_data_markers %>% 
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

names(Immnue_data_top5_markers) <- sub('.*\\.','',names(Immnue_data_top5_markers))

# dotplot
pdf(file.path(out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 36, height = 6)
#table(unlist(KC_major_marker_ls))[table(unlist(KC_major_marker_ls)) > 1]
DotPlot(Immnue_data, group.by = "RNA_snn_res.0.6", features = Immnue_major_marker_ls)+RotatedAxis()
dev.off()


pdf(file.path(out_dir, str_glue("dotplot dims{dims} top5_marker.pdf")), width = 30, height = 6)

#table(Immnue_data$RNA_snn_res.0.2)
#levels(factor(Immnue_data$RNA_snn_res.0.2))
DotPlot(Immnue_data, group.by = "RNA_snn_res.0.6", features = Immnue_data_top5_markers)+RotatedAxis()
dev.off()

features_all <- unique(unlist(Immnue_major_marker_ls))


pdf(file.path(out_dir, str_glue("dotplot dims{dims} gene_marker.pdf")), width = 36, height = 6)
#table(unlist(KC_major_marker_ls))[table(unlist(KC_major_marker_ls)) > 1]
DotPlot(Immnue_data, group.by = "RNA_snn_res.0.6", features = features_all) +
  RotatedAxis()
dev.off()

Immnue_data$subtype <- case_when(
  Immnue_data$RNA_snn_res.0.6 == 0 ~ "Mac",
  Immnue_data$RNA_snn_res.0.6 == 1 ~ "Mac",
  Immnue_data$RNA_snn_res.0.6 == 2 ~ "CD4+ Th",
  Immnue_data$RNA_snn_res.0.6 == 3 ~ "CD8+ Tc",
  Immnue_data$RNA_snn_res.0.6 == 4 ~ "CD1C+ DC",
  Immnue_data$RNA_snn_res.0.6 == 5 ~ "LC",
  Immnue_data$RNA_snn_res.0.6 == 6 ~ "Mast",
  Immnue_data$RNA_snn_res.0.6 == 7 ~ "Mast",
  Immnue_data$RNA_snn_res.0.6 == 8 ~ "CD4+ Treg",
  Immnue_data$RNA_snn_res.0.6 == 9 ~ "Cyc Imm",
  Immnue_data$RNA_snn_res.0.6 == 10 ~ "CLEC9A+ DC",
  Immnue_data$RNA_snn_res.0.6 == 11 ~ "CCR7+ DC",
  Immnue_data$RNA_snn_res.0.6 == 12 ~ "Plasma",
  Immnue_data$RNA_snn_res.0.6 == 13 ~ "CLEC9A+ DC",
  Immnue_data$RNA_snn_res.0.6 == 14 ~ "CD1C+ DC",
  TRUE ~ "Unknown"  # 其他情况
)

immune_colors <- c(
  "Mac"              = "#9370DB",  # MediumPurple
  "CD1C+ DC"         = "#6495ED",  # CornflowerBlue
  "CLEC9A+ DC"       = "#DEB887",  # BurlyWood
  "CCR7+ DC"         = "#B0C4DE",  # LightSteelBlue
  "LC"               = "#800080",  # Purple
  "Mast"             = "#E6E682",  # Khaki
  "CD4+ Th"          = "#FF0000",  # Red
  "CD4+ Treg"        = "#006400",  # DarkGreen
  "CD8+ Tc"          = "#4B0082",  # Indigo
  "Naive T"          = "#1E90FF",  # DodgerBlue (已修复原与ORS basal冲突的#4682B4)
  "Cyc Imm"          = "#87CEFA",  # LightSkyBlue
  "Mono"             = "#FFA07A",  # LightSalmon
  "Plasma"           = "#FF69B4"  # HotPink
)
pdf(file.path(out_dir, str_glue("Dimplot dim{dims} subtype.pdf")), height = 6, width = 6)
DimPlot(Immnue_data, group.by = "subtype", label = TRUE, reduction = 'harmonyumap') +
  scale_color_manual(values = immune_colors) +
  NoLegend()
dev.off()

pdf(file.path(out_dir, str_glue("nCount_RNA.pdf")), height = 6, width = 6)
VlnPlot(Immnue_data, features = "nCount_RNA", pt.size = 0) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
####stromal_subtype
out_dir <- "stroma_subtype"  
dir.create(out_dir)
unique(all_data$annotation)
stroma_data <- subset(all_data, annotation %in% c("Fibro","EC","Peri","Schwann","Adipo","LEC"))

stroma_data <- JoinLayers(stroma_data)
## PCA --------
stroma_data <- FindVariableFeatures(stroma_data)
length(VariableFeatures(stroma_data)) #2000
stroma_data <- ScaleData(stroma_data)
stroma_data <- RunPCA(stroma_data)

### pca inspector ------
pdf(file.path(out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(stroma_data, ndims = 50)
dev.off()

pdf(file.path(out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(stroma_data, dims = 1:5)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(stroma_data, dims = 6:10)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(stroma_data, dims = 11:15)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(stroma_data, dims = 16:20)
dev.off()

## dimension reduction --------
dims <- 20 # 18

library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

stroma_data <- RunUMAP(stroma_data, dims = 1:dims,return.model = TRUE)

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 8, width = 8)
DimPlot(stroma_data, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

stroma_data <- stroma_data %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

stroma_data <- RunUMAP(
  stroma_data, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 8, width = 8)
DimPlot(stroma_data, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
stroma_data  <-
  FindNeighbors(
    stroma_data ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


stroma_data <-
  FindClusters(
    stroma_data,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 6)
DimPlot(stroma_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(stroma_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(stroma_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(stroma_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

Idents(stroma_data) <- "RNA_snn_res.0.8"

### dotplot -------
stromal_major_marker_ls <- list(
  Adipo          = c("ADIPOQ", "G0S2", "LPL", "S100B", "LTF"),
  DP             = c("COCH", "RSPO4", "TNN", "BMP4", "CRABP1"),
  DS             = c("COL11A1", "LRRC15", "TNN", "DPEP1", "POSTN"),
  HEC            = c("ACKR1", "SELE", "VWF", "FLT1", "IL3RA"),
  LEC            = c("CCL21", "TFF3", "ACKR2", "ANGPT2", "LYVE1"),
  `Papil Fib`    = c("COMP", "RSPO1", "APCDD1", "STC2", "AXIN2"),
  Peri           = c("ACTA2", "EDNRA", "RGS5", "TINAGL1", "MCAM"),
  `Perivasc Fib I` = c("CCL19", "C3", "CXCL1", "ABI3BP", "CSF1"),
  `Perivasc Fib II`= c("MYOC", "GPC3", "C3", "PDGFRA", "PLA2G2A"),
  `Retic Fib I`  = c("PRG4", "MMP2", "MMP3", "PI16", "COL1A1"),
  `Retic Fib II` = c("ANGPTL7", "COMP", "RSPO1", "DKK3", "STC2"),
  `Retic Fib III`= c("ANGPTL7", "PTGDS", "GLI1", "SFRP4", "PTCH2"),
  Schwann        = c("MPZ", "SOX2", "NRXN1", "SOX10", "S100B"),
  SM             = c("DES", "ITGA8", "MYH11", "PLN", "PCP4"),
  VEC            = c("NOTCH4", "DLL4", "EDN1", "IGFBP3", "IFI44L")
)
stroma_major_marker_ls <- list(
  Adipo          = c("ADIPOQ", "G0S2", "LPL", "S100B", "LTF"),
  DP             = c("COCH", "RSPO4", "TNN", "BMP4", "CRABP1"),
  DS             = c("COL11A1", "LRRC15", "DPEP1", "POSTN"),
  HEC            = c("ACKR1", "SELE", "VWF", "FLT1", "IL3RA"),
  LEC            = c("CCL21", "TFF3", "ACKR2", "ANGPT2", "LYVE1"),
  `Papil Fib`    = c("COMP", "RSPO1", "APCDD1", "STC2", "AXIN2"),
  Peri           = c("ACTA2", "EDNRA", "RGS5", "TINAGL1", "MCAM"),
  `Perivasc Fib I` = c("CCL19", "C3", "CXCL1", "ABI3BP", "CSF1"),
  `Perivasc Fib II`= c("MYOC", "GPC3", "PDGFRA", "PLA2G2A"),
  `Retic Fib I`  = c("PRG4", "MMP2", "MMP3", "PI16", "COL1A1"),
  `Retic Fib II` = c("ANGPTL7", "DKK3"),
  `Retic Fib III`= c("PTGDS", "GLI1", "SFRP4", "PTCH2"),
  Schwann        = c("MPZ", "SOX2", "NRXN1", "SOX10"),
  SM             = c("DES", "ITGA8", "MYH11", "PLN", "PCP4"),
  VEC            = c("NOTCH4", "DLL4", "EDN1", "IGFBP3", "IFI44L")
)

stroma_major_marker_ls <- lapply(stroma_major_marker_ls, function(genes) {
  genes[genes %in% rownames(stroma_data)]
})

#KC_major_marker_ls <- c(
#  KC_basal = 'COL17A1',
#  KC_basal = 'POSTN',
#  KC_basal = 'ITGB4',
#  KC_spinous = 'KRT1',
#  KC_granular = 'KRT2',
#  KC_activated = 'KRT17',
#  KC_eccrine = 'CFTR'
#)


Idents(stroma_data)<-'RNA_snn_res.0.8'
stroma_data_markers <-
  FindAllMarkers(
    JoinLayers(stroma_data),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )

stroma_data_top5_markers <- 
  stroma_data_markers %>% 
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

names(stroma_data_top5_markers) <- sub('.*\\.','',names(stroma_data_top5_markers))

# dotplot
pdf(file.path(out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 40, height = 6)
#table(unlist(KC_major_marker_ls))[table(unlist(KC_major_marker_ls)) > 1]
DotPlot(stroma_data, group.by = "RNA_snn_res.0.8", features = stroma_major_marker_ls)+RotatedAxis()
dev.off()


pdf(file.path(out_dir, str_glue("dotplot dims{dims} top5_marker.pdf")), width = 30, height = 6)

#table(stroma_data$RNA_snn_res.0.2)
#levels(factor(stroma_data$RNA_snn_res.0.2))
DotPlot(stroma_data, group.by = "RNA_snn_res.0.8", features = stroma_data_top5_markers)+RotatedAxis()
dev.off()

features_all <- unique(unlist(stroma_major_marker_ls))


pdf(file.path(out_dir, str_glue("dotplot dims{dims} gene_marker.pdf")), width = 36, height = 6)
#table(unlist(KC_major_marker_ls))[table(unlist(KC_major_marker_ls)) > 1]
DotPlot(stroma_data, group.by = "RNA_snn_res.0.6", features = features_all) +
  RotatedAxis()
dev.off()

stroma_data$subtype <- case_when(
  stroma_data$RNA_snn_res.0.8 == 0 ~ "Retic Fib I",
  stroma_data$RNA_snn_res.0.8 == 1 ~ "Retic Fib I",
  stroma_data$RNA_snn_res.0.8 == 2 ~ "Perivasc Fib I",
  stroma_data$RNA_snn_res.0.8 == 3 ~ "HEC",
  stroma_data$RNA_snn_res.0.8 == 4 ~ "Peri",
  stroma_data$RNA_snn_res.0.8 == 5 ~ "DS",
  stroma_data$RNA_snn_res.0.8 == 6 ~ "Papil Fib",
  stroma_data$RNA_snn_res.0.8 == 7 ~ "VEC",
  stroma_data$RNA_snn_res.0.8 == 8 ~ "VEC",
  stroma_data$RNA_snn_res.0.8 == 9 ~ "Perivasc Fib II",
  stroma_data$RNA_snn_res.0.8 == 10 ~ "DS",
  stroma_data$RNA_snn_res.0.8 == 11 ~ "Adipo",
  stroma_data$RNA_snn_res.0.8 == 12 ~ "LEC",
  stroma_data$RNA_snn_res.0.8 == 13 ~ "Schwann",
  stroma_data$RNA_snn_res.0.8 == 14 ~ "Peri",
  stroma_data$RNA_snn_res.0.8 == 15 ~ "SM",
  stroma_data$RNA_snn_res.0.8 == 16 ~ "Peri",
  TRUE ~ "Unknown"  # 其他情况
)


stroma_colors <- c(
  "Retic Fib I"      = "#DC143C",  # Crimson
  "Perivasc Fib I"   = "#20B2AA",  # LightSeaGreen
  "Perivasc Fib II"  = "#66CDAA",  # MediumAquamarine (已修复原与Diff Seb冲突的#3CB371)
  "Papil Fib"        = "#FF8C00",  # DarkOrange
  "DP"               = "#DAA520",  # GoldenRod (已修复原与Diff Inf和Peri冲突的#FFD700)
  "DS"               = "#DDA0DD",  # Plum
  "Adipo"            = "#FA8072",  # Salmon (已修复原与Mono冲突的#FFA07A)
  "HEC"              = "#4169E1",  # RoyalBlue
  "VEC"              = "#0000CD",  # MediumBlue (已修复原与HEC相近的#4169E1)
  "Peri"             = "#B8860B",  # DarkGoldenRod (已修复原与Diff Inf和DP冲突的#FFD700)
  "SM"               = "#E6E6FA",  # Lavender
  "Schwann"          = "#A0522D",  # Sienna (已修复原与Spn KC I冲突的#8B4513)
  "Perineural Fib"   = "#98FB98",  # PaleGreen
  "LEC"              = "#DB7093" 
)


pdf(file.path(out_dir, str_glue("Dimplot dim{dims} subtype.pdf")), height = 6, width = 6)
DimPlot(stroma_data, group.by = "subtype", label = TRUE, reduction = 'harmonyumap') +
  scale_color_manual(values = stroma_colors) +
  NoLegend()
dev.off()

pdf(file.path(out_dir, str_glue("nCount_RNA.pdf")), height = 6, width = 6)
VlnPlot(stroma_data, features = "nCount_RNA", pt.size = 0) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
#########
########ALL_CELL
all_data$subtype <- as.character(all_data$annotation)
unique(all_data$subtype)

a <- as.data.frame(table(all_data$subtype))
b <- a
# 2. 映射 KC 亚型
kc_cells <- colnames(KC_data)
kc_subtype <- KC_data$subtype
subtype_map_kc <- setNames(as.character(kc_subtype), kc_cells)
all_data$subtype[match(kc_cells, colnames(all_data))] <- subtype_map_kc[kc_cells]

# 3. 映射 免疫亚型（新增）
immune_cells <- colnames(Immnue_data)
immune_subtype <- Immnue_data$subtype
subtype_map_immune <- setNames(as.character(immune_subtype), immune_cells)
all_data$subtype[match(immune_cells, colnames(all_data))] <- subtype_map_immune[immune_cells]

# 4. 映射 基质亚型（新增）
stroma_cells <- colnames(stroma_data)
stroma_subtype <- stroma_data$subtype
subtype_map_stroma <- setNames(as.character(stroma_subtype), stroma_cells)
all_data$subtype[match(stroma_cells, colnames(all_data))] <- subtype_map_stroma[stroma_cells]

unique(all_data$subtype)
#############################空间细胞类型原位
all_data_1 <- subset(all_data, subset = subtype != "Others")

out_dir <- "Celltype_spatial"
dir.create(out_dir)
annotation_col <- 'subtype'
length(cell_colors)
cell_colors <- c(
  # 上皮细胞 (KC_data)
  "Spn KC I"         = "#8B4513",
  "Bas KC"           = "#CD5C5C",
  "Grn KC"           = "#D2B48C",
  "Spn KC II"        = "#2E8B57",
  "Cyc KC"           = "#F08080",
  #"Bas Inf"          = "#FF7F50",
  "Inf"         = "#FFD700",
  #"Diff Seb"         = "#3CB371",
  "Seb"          = "#90EE90",
  #"Bulge"            = "#9ACD32",
  "Melano"           = "#7CFC00",
  "Ecc duct"         = "#00FA9A",
  "Ecc gland"        = "#00FF7F",
  "Ecc Myoepi"       = "#00CED1",
  "ORS basal"        = "#4682B4",
  "ORS suprabasal"   = "#5F9EA0",
  "IRS/HS"               = "#48D1CC",
  "Others"           = "#808080",
  
  # 免疫细胞 (Immnue_data)
  "Mac"              = "#9370DB",
  "CD1C+ DC"         = "#6495ED",
  "CLEC9A+ DC"       = "#DEB887",
  "CCR7+ DC"         = "#B0C4DE",
  "LC"               = "#800080",
  "Mast"             = "#E6E682",
  "CD4+ Th"          = "#FF0000",
  "CD4+ Treg"        = "#006400",
  "CD8+ Tc"          = "#4B0082",
  "Naive T"          = "#4682B4",
  "Cyc Imm"          = "#87CEFA",
  "Mono"             = "#FFA07A",
  "Plasma"           = "#FF69B4",
  
  # 基质细胞 (stroma_data)
  "Retic Fib I"      = "#DC143C",
  "Perivasc Fib I"   = "#20B2AA",
  "Perivasc Fib II"  = "#3CB371",
  "Papil Fib"        = "#FF8C00",
  "DP"               = "#FFD700",
  "DS"               = "#DDA0DD",
  "Adipo"            = "#FFA07A",
  "HEC"              = "#1E90FF",
  "VEC"              = "#4169E1",
  "Peri"             = "#FFD700",
  "SM"               = "#E6E6FA",
  "Schwann"          = "#8B4513",
  "Perineural Fib"   = "#98FB98",
  "LEC"              = "#F08080"
)


cell_colors <- c(
  # 上皮细胞 (KC_data) - 暖色系为主
  "Spn KC I"         = "#8B4513",  # SaddleBrown
  "Bas KC"           = "#CD5C5C",  # IndianRed
  "Grn KC"           = "#D2B48C",  # Tan
  "Spn KC II"        = "#2E8B57",  # SeaGreen
  "Cyc KC"           = "#F08080",  # LightCoral
  #"Bas Inf"          = "#FF7F50",  # Coral
  "Inf"         = "#FFD700",  # Gold
  #"Diff Seb"         = "#3CB371",  # MediumSeaGreen
  "Seb"          = "#90EE90",  # LightGreen
  #"Bulge"            = "#9ACD32", # YellowGreen
  "Melano"           = "#7CFC00",  # LawnGreen
  "Ecc duct"         = "#00FA9A",  # MediumSpringGreen
  "Ecc gland"        = "#00FF7F",  # SpringGreen
  "Ecc Myoepi"       = "#00CED1",  # DarkTurquoise
  "ORS basal"        = "#4682B4",  # SteelBlue
  "ORS suprabasal"   = "#5F9EA0",  # CadetBlue
  "IRS/HS"               = "#48D1CC",  # MediumTurquoise
  "Others"           = "#808080",  # Gray
  
  # 免疫细胞 (Immune_data) - 冷色/中性系
  "Mac"              = "#9370DB",  # MediumPurple
  "CD1C+ DC"         = "#6495ED",  # CornflowerBlue
  "CLEC9A+ DC"       = "#DEB887",  # BurlyWood
  "CCR7+ DC"         = "#B0C4DE",  # LightSteelBlue
  "LC"               = "#800080",  # Purple
  "Mast"             = "#E6E682",  # Khaki
  "CD4+ Th"          = "#FF0000",  # Red
  "CD4+ Treg"        = "#006400",  # DarkGreen
  "CD8+ Tc"          = "#4B0082",  # Indigo
  "Naive T"          = "#1E90FF",  # DodgerBlue (已修复原与ORS basal冲突的#4682B4)
  "Cyc Imm"          = "#87CEFA",  # LightSkyBlue
  "Mono"             = "#FFA07A",  # LightSalmon
  "Plasma"           = "#FF69B4",  # HotPink
  
  # 基质细胞 (stroma_data) - 青绿/蓝色系
  "Retic Fib I"      = "#DC143C",  # Crimson
  "Perivasc Fib I"   = "#20B2AA",  # LightSeaGreen
  "Perivasc Fib II"  = "#66CDAA",  # MediumAquamarine (已修复原与Diff Seb冲突的#3CB371)
  "Papil Fib"        = "#FF8C00",  # DarkOrange
  "DP"               = "#DAA520",  # GoldenRod (已修复原与Diff Inf和Peri冲突的#FFD700)
  "DS"               = "#DDA0DD",  # Plum
  "Adipo"            = "#FA8072",  # Salmon (已修复原与Mono冲突的#FFA07A)
  "HEC"              = "#4169E1",  # RoyalBlue
  "VEC"              = "#0000CD",  # MediumBlue (已修复原与HEC相近的#4169E1)
  "Peri"             = "#B8860B",  # DarkGoldenRod (已修复原与Diff Inf和DP冲突的#FFD700)
  "SM"               = "#E6E6FA",  # Lavender
  "Schwann"          = "#A0522D",  # Sienna (已修复原与Spn KC I冲突的#8B4513)
  "Perineural Fib"   = "#98FB98",  # PaleGreen
  "LEC"              = "#DB7093"   # PaleVioletRed (已修复原与Cyc KC冲突的#F08080)
)
cells_position_ls <-
  cbind(
    all_data_1[['sample']],
    Embeddings(all_data_1[['spatial']])
  ) %>%
  rownames_to_column('cell_id') %>%
  group_by(sample) %>%
  nest() %>%
  deframe()


plot_dat <-
  FetchData(
    all_data_1, 
    vars = c(annotation_col)
  ) %>%
  rownames_to_column('cell_id')

p_sp_ct_ls <-
  cells_position_ls %>% # .[1] %>%
  imap(function(coords, sam){
    dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
    
    height <- 8
    width <- 8
    p <-
      dat2 %>%
      ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
      geom_tile(
        aes(fill = .data[[annotation_col]]),
        height = height, 
        width = width,
        show.legend = c(fill = TRUE, alpha = FALSE)
      ) +
      theme_void() +
      coord_equal() +
      scale_y_reverse() +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values = cell_colors) +
      theme(panel.background = element_rect(fill = "black", color = NA)) +
      ggtitle(sam)
    p
    
  })

#p_sp_ct_ls[4]
pdf(
  file.path(out_dir, "1. spatial plot of subtype in each sample.pdf"),
  height = 8,
  width = 8
)
p_sp_ct_ls %>% walk(print)
dev.off()
####HRAS
out_dir <- "HRAS_mut"
some_genes <- c("HRAS", "HRAS-ALT-37:C", "HRAS-WT-37:C")
annotation_col <- 'subtype'

# 定义每个基因的颜色
gene_colors <- list(
  "HRAS" = c("white", "blue"),           # 白→蓝
  "HRAS-ALT-37:C" = c("white", "red"),   # 白→红
  "HRAS-WT-37:C" = c("white", "yellow")  # 白→黄
)

cells_position_ls <-
  cbind(
    all_data_1[['sample']],
    Embeddings(all_data_1[['spatial']])
  ) %>%
  rownames_to_column('cell_id') %>%
  group_by(sample) %>%
  nest() %>%
  deframe()

plot_dat <-
  FetchData(
    all_data_1, 
    vars = c(annotation_col, some_genes)
  ) %>%
  rownames_to_column('cell_id')

p_sp_gene_ls <-
  cells_position_ls %>%
  imap(function(coords, sam){
    some_genes %>%
      setNames(., .) %>%
      map(function(x){
        dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
        
        stopifnot(x %in% colnames(dat2))
        
        # quantile cutoff
        q_cutoff <- 0.9
        n <- 0
        while( near(quantile(dat2[[x]], q_cutoff), 0) ){
          q_cutoff <- q_cutoff + 9/10^(n + 2)
          print(q_cutoff)
          n <- n + 1
          if(n > 3) break
        }
        
        if(n <= 3){
          dat2[[x]] <- pmin(dat2[[x]], quantile(dat2[[x]], q_cutoff))
        }
        
        val_rng <- dat2[[x]] %>% range()
        height <- 8
        width <- 8
        
        # 获取当前基因的颜色
        colors <- gene_colors[[x]]
        
        p <-
          dat2 %>%
          ggplot(aes(x = SPATIAL_1, y = SPATIAL_2)) +
          # 背景灰色点
          geom_tile(fill = "gray30", height = height, width = width) + 
          # 基因表达点
          geom_tile(
            aes(fill = .data[[x]], alpha = .data[[x]]),
            height = height, 
            width = width,
            show.legend = c(fill = TRUE, alpha = FALSE)
          ) +
          theme_void() +
          coord_equal() +
          scale_y_reverse() +
          # 使用基因特异性颜色
          scale_fill_gradient(low = colors[1], high = colors[2], limits = val_rng) +
          scale_alpha_continuous(range = c(0, 1), guide = "none") +
          theme(panel.background = element_rect(fill = "black", color = NA)) +
          ggtitle(str_glue("{sam} - {x}"))
        
        p
      })
  })

# 导出
p_sp_gene_ls %>%
  iwalk(function(p_ls, sam){
    pdf(
      file.path(out_dir, str_glue("spatial plot of HRAS mut in sample {sam}.pdf")),
      height = 8,
      width = 8
    )
    p_ls %>% walk(print)
    dev.off()
  })
##Xenium Explorer 4 CSV文件
all_meta <- all_data@meta.data
all_meta <- all_meta[,c("cell_id2","subtype")]

all_meta$color <- cell_colors[all_meta$subtype]
all_meta$color[is.na(all_meta$color)] <- "#A9A9A9" 
out_dir <- "Celltype_spatial"
colnames(all_meta) <- c("cell_id","group","color")
write.csv(all_meta,file.path(out_dir, "cell_id_color_subtype.csv"),row.names = FALSE)
unique(all_data$subtype)
length(unique(all_meta$group))
length(unique(all_meta$color))
#####
target_genes <- c("HRAS", "HRAS-ALT-37:C", "HRAS-WT-37:C")

# 初始化一个空的数据框用来存放结果
summary_result <- data.frame()

for (gene in target_genes) {
  # 1. 提取当前基因的表达数据
  gene_expr <- FetchData(all_data_1, vars = gene)[, 1]
  
  # 2. 判定表达并计算比例
  is_expressed <- gene_expr > 0
  total_cells <- ncol(all_data_1)
  expressed_cells <- sum(is_expressed)
  proportion <- round(expressed_cells / total_cells * 100, 2)
  
  # 3. 存入结果表
  temp_df <- data.frame(
    Gene = gene,
    Expressed_Cells = expressed_cells,
    Total_Cells = total_cells,
    Proportion_Percent = proportion
  )
  summary_result <- rbind(summary_result, temp_df)
}

# 查看汇总结果
print(summary_result)

target_genes <- c("HRAS", "HRAS-ALT-37:C", "HRAS-WT-37:C")

# 使用 lapply 循环处理，返回一个包含各基因结果的列表
all_gene_results <- lapply(target_genes, function(gene) {
  
  # 提取数据与判定表达
  gene_expr <- FetchData(all_data_1, vars = gene)[, 1]
  is_expressed <- gene_expr > 0
  
  # 提取表达该基因的细胞的亚群信息
  expressed_subtype <- all_data_1$subtype[is_expressed]
  
  # 如果该基因没有任何细胞表达，则跳过计算
  if (length(expressed_subtype) == 0) {
    return(data.frame(Gene = gene, Subtype = "None", Count = 0, Percentage = 0))
  }
  
  # 计算各亚群的数量和比例
  subtype_counts <- table(expressed_subtype)
  subtype_proportions <- prop.table(subtype_counts)
  
  # 整理为数据框，并加上基因名列
  result_df <- data.frame(
    Gene = gene,
    Subtype = names(subtype_counts),
    Count = as.numeric(subtype_counts),
    Percentage = round(as.numeric(subtype_proportions) * 100, 2)
  )
  
  return(result_df)
})

# 将所有基因的结果合并成一张大表格
final_table <- bind_rows(all_gene_results)
out_dir <- "HRAS_mut"
cairo_pdf(file.path(out_dir, "HRAS_allprop_barplot.pdf"), width = 10, height = 6)

# 2. 将之前的绘图代码完整放入 print() 中
print(
  ggplot(final_table, aes(x = Gene, y = Percentage, fill = Subtype)) +
    geom_bar(stat = "identity", width = 0.6, color = "white", size = 0.1) +
    
    scale_fill_manual(values = cell_colors, drop = FALSE) +
    
    theme_minimal(base_size = 14) +
    theme(
      # 🚫 去除所有背景网格线
      panel.grid = element_blank(), 
      
      # ➕ 添加黑色的 X 轴和 Y 轴线（粗细设为 0.5）
      axis.line = element_line(color = "black", size = 0.5),
      
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, vjust = 0.5),
      axis.text.y = element_text(color = "black"),
      legend.position = "right",
      legend.text = element_text(size = 10),
      legend.title = element_text(face = "bold")
    ) +
    labs(
      x = "Gene",
      y = "Proportion (%)",
      fill = "Celltype",
    ) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) 
)

# 3. 关闭设备并写入文件
dev.off()

#######等高线
target_genes <- c("HRAS", "HRAS-ALT-37:C", "HRAS-WT-37:C")
annotation_col <- 'subtype'
HRAS_position_ls <- map(cells_position_ls, function(sample_df) {
  
  # 3.1 从 Seurat 对象中提取当前样本细胞的表达量
  # 注意：这里假设 sample_df$cell_id 中的 ID 与 all_data 的 colnames 完全匹配
  # FetchData 会自动处理 cell_id 的匹配
  
  # 提取表达矩阵子集
  gene_expr_matrix <- FetchData(all_data_1, 
                                vars = c(annotation_col,target_genes), 
                                cells = sample_df$cell_id)
  
  # 3.2 将提取的矩阵转换为 data.frame 并加上 cell_id 列以便合并
  # rownames(gene_expr_matrix) 就是 cell_id
  expr_df <- gene_expr_matrix %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "cell_id")
  
  # 3.3 将表达量数据 (expr_df) 左连接合并到原始坐标数据 (sample_df) 中
  # 使用 left_join 保证原始坐标数据的顺序和行数不变
  final_df <- left_join(sample_df, expr_df, by = "cell_id")
  
  return(final_df)
})



dat =HRAS_position_ls[[1]]

kde_ls <- 
  dat %>% 
  remove_missing() %>%
  group_by(subtype) %>%
  nest() %>%
  deframe() %>%
  
  map(function(x){
    kde <- MASS::kde2d(x$SPATIAL_1, x$SPATIAL_2)
    data.frame(
      z = kde$z %>% c(),
      x = kde$x %>% rep(times = length(kde$y)),
      y = kde$y %>% rep(each  = length(kde$x))
    )
  })


kde_df <- 
  kde_ls %>%
  enframe("subtype") %>%
  unnest(value)

kde_df %>% 
  ggplot(aes(x = x, y = y, z = z, color = subtype)) +
  geom_contour() +
  theme_void()

kde_df %>% 
  ggplot(aes(x = x, y = y, z = z)) +
  geom_contour(aes(color = after_stat(nlevel))) +
  scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(Inf, "RdBu"))) +
  facet_wrap(~ subtype) +
  theme_void()
ggsave("contour plot - levels.pdf")

dat2 %>%
  ggplot(aes(x =SPATIAL_1 , SPATIAL_2)) + 
  # geom_point(aes(color = annotation3))
  geom_density_2d( ) +
  facet_wrap(~ subtype)


dat2 %>% 
  # head(100000) %>%
  ggplot(aes(x = SPATIAL_1, SPATIAL_2)) + 
  # geom_point(aes(color = annotation3)) +
  geom_tile(aes(fill = subtype, height = 10, width = 10), alpha = 0.5) +
  
  geom_density2d(aes(color = subtype), contour_var = "density", show.legend = FALSE ) +
  # ggsci::scale_color_aaas() +
  # ggsci::scale_fill_aaas() +
  coord_equal() +
  theme_void()

ggsave("contour plot - combine.pdf")
file.remove("contour plot - levels.pdf")
#####HRAS
all_cells <- dat %>% remove_missing()
hras_positive_cells <- all_cells %>% filter(HRAS > 0)

# 2. 只对 HRAS 阳性细胞进行二维核密度估计 (KDE)
kde_result <- MASS::kde2d(
  x = hras_positive_cells$SPATIAL_1, 
  y = -hras_positive_cells$SPATIAL_2, 
  n = 100            # 网格分辨率
)

# 3. 将 KDE 结果转换为 ggplot 可用的长格式数据框
kde_df <- data.frame(
  x = rep(kde_result$x, times = length(kde_result$y)),
  y = rep(kde_result$y, each = length(kde_result$x)),
  z = as.vector(kde_result$z)
)

# 4. 分层绘图
out_dir <- "contour_plot"
pdf(file.path(out_dir, "contour_plot_HRAS.pdf"), width = 10, height = 6)
print(
  ggplot() +
  # 【第一层：背景】画出所有细胞的淡淡散点，保留完整空间位点
  geom_tile(
    data = all_cells,
    aes(x = SPATIAL_1, y = -SPATIAL_2, fill = subtype), # 注意这里 y 也要取负
    height = 8,
    width = 8
  ) +
  scale_fill_manual(values = cell_colors) +
  
  # 【第二层：前景】叠加 HRAS 阳性细胞的密度等高线
  geom_contour(data = kde_df, 
               aes(x = x, y = y, z = z, color = after_stat(level)), 
               size = 0.8) +
  
  # 设置配色与主题
  scale_color_viridis_c(option = "C") +
  theme_void() +
  
  coord_equal() +
  labs(fill= "Celltype",
       color = "Density")+
  theme(
    legend.key.size = unit(0.5, 'cm'),       # 缩小图例中色块/线条图标的大小
    legend.text = element_text(size = 10),   # 调整图例标签文字大小
    legend.title = element_text(size = 12)   # 调整图例标题文字大小
  )
)
dev.off()
#ggsave("contour_plot_HRAS.pdf", plot = last_plot(), width = 10, height = 12, units = "cm", dpi = 300)

p_hras_contour_ls <- 
  HRAS_position_ls %>% 
  imap(function(coords, sam){
    
    # --- 1. 准备当前样本的绘图数据 ---
    # coords 就是当前样本的 HRAS_position_ls[[i]] (即之前的 dat)
    all_cells <- coords %>% remove_missing()
    hras_positive_cells <- all_cells %>% filter(HRAS > 0)
    
    # 如果该样本没有阳性细胞，直接返回空，防止后续报错
    if (nrow(hras_positive_cells) == 0) return(NULL)
    
    # --- 2. 计算当前样本的二维核密度估计 (KDE) ---
    kde_result <- MASS::kde2d(
      x = hras_positive_cells$SPATIAL_1, 
      y = -hras_positive_cells$SPATIAL_2, 
      n = 100            
    )
    
    # 转换为长格式数据框
    kde_df <- data.frame(
      x = rep(kde_result$x, times = length(kde_result$y)),
      y = rep(kde_result$y, each = length(kde_result$x)),
      z = as.vector(kde_result$z)
    )
    
    # --- 3. 绘制当前样本的叠加图 ---
    p <- ggplot() +
      # 【第一层：背景】细胞类型注释底图
      geom_tile(
        data = all_cells,
        aes(x = SPATIAL_1, y = -SPATIAL_2, fill = subtype), # Y轴同步翻转
        height = 8,
        width = 8
      ) +
      scale_fill_manual(values = cell_colors) +
      
      # 【第二层：前景】HRAS 阳性细胞的密度等高线
      geom_contour(data = kde_df, 
                   aes(x = x, y = y, z = z, color = after_stat(level)), 
                   size = 0.8) +
      
      # 设置配色、主题与标签
      scale_color_viridis_c(option = "C") +
      theme_void() +
      coord_equal() +
      labs(fill = "Celltype", color = "Density") +
      theme(
        panel.background = element_rect(fill = "white", color = NA), # 保持黑色背景
        legend.key.size = unit(0.5, 'cm'),       
        legend.text = element_text(size = 10),   
        legend.title = element_text(size = 12)   
      ) +
      ggtitle(sam)+ # 加上样本名称作为标题
      guides(
        fill = guide_legend(order = 1),       # 离散图例用 guide_legend
        color = guide_colorbar(order = 2)     # 连续图例用 guide_colorbar
      )
    
    return(p)
  })
pdf(
  file.path(out_dir, "contour_plot_HRAS in each sample.pdf"),
  height = 6,
  width = 10
)
p_hras_contour_ls %>% walk(print)
dev.off()

####alt
p_hras_contour_ls <- 
  HRAS_position_ls %>% 
  imap(function(coords, sam){
    
    # --- 1. 准备当前样本的绘图数据 ---
    # coords 就是当前样本的 HRAS_position_ls[[i]] (即之前的 dat)
    all_cells <- coords %>% remove_missing()
    hras_positive_cells <- all_cells %>% filter(`HRAS-WT-37:C` > 0)
    
    # 如果该样本没有阳性细胞，直接返回空，防止后续报错
    if (nrow(hras_positive_cells) == 0) return(NULL)
    
    # --- 2. 计算当前样本的二维核密度估计 (KDE) ---
    kde_result <- MASS::kde2d(
      x = hras_positive_cells$SPATIAL_1, 
      y = -hras_positive_cells$SPATIAL_2, 
      n = 100            
    )
    
    # 转换为长格式数据框
    kde_df <- data.frame(
      x = rep(kde_result$x, times = length(kde_result$y)),
      y = rep(kde_result$y, each = length(kde_result$x)),
      z = as.vector(kde_result$z)
    )
    
    # --- 3. 绘制当前样本的叠加图 ---
    p <- ggplot() +
      # 【第一层：背景】细胞类型注释底图
      geom_tile(
        data = all_cells,
        aes(x = SPATIAL_1, y = -SPATIAL_2, fill = subtype), # Y轴同步翻转
        height = 8,
        width = 8
      ) +
      scale_fill_manual(values = cell_colors) +
      
      # 【第二层：前景】HRAS 阳性细胞的密度等高线
      geom_contour(data = kde_df, 
                   aes(x = x, y = y, z = z, color = after_stat(level)), 
                   size = 0.8) +
      
      # 设置配色、主题与标签
      scale_color_viridis_c(option = "C") +
      theme_void() +
      coord_equal() +
      labs(fill = "Celltype", color = "Density") +
      theme(
        panel.background = element_rect(fill = "white", color = NA), # 保持黑色背景
        legend.key.size = unit(0.5, 'cm'),       
        legend.text = element_text(size = 10),   
        legend.title = element_text(size = 12)   
      ) +
      ggtitle(sam)+ # 加上样本名称作为标题
      guides(
        fill = guide_legend(order = 1),       # 离散图例用 guide_legend
        color = guide_colorbar(order = 2)     # 连续图例用 guide_colorbar
      )
    
    return(p)
  })
pdf(
  file.path(out_dir, "contour_plot_HRAS-WT-37_C in each sample.pdf"),
  height = 6,
  width = 10
)
p_hras_contour_ls %>% walk(print)
dev.off()
#########

KC_major_marker_ls <- c(
  # 角质形成细胞
  `Cyc KC` = "MKI67",
  `Cyc KC` = "TOP2A",
  
  `Bas KC` = "COL17A1",
  `Bas KC` = "COL7A1",
  
  `Spn KC` = "RHOV",
  `Spn KC` = "KLF5",
  
  `Grn KC` = "DSC1",
  
  HFE = "GJB6",
  HFE = "SOX9",
  
  `Activated KC`="KRT17",
  
  `Ecc Duct` = "SEMA3C",
  
  Others= "SOX2-OT",
  
  Mast = "CPA",
  Mast = "TPSAB1",
  Mast = "KIT"
)

KC_major_marker_ls <- KC_major_marker_ls[KC_major_marker_ls%in%rownames(KC_data)]
pdf(file.path(out_dir, str_glue("dotplot dims{dims} major_marker_refer2.pdf")), width = 20, height = 6)
DotPlot(KC_data, group.by = "RNA_snn_res.0.4", features = KC_major_marker_ls)+RotatedAxis()
dev.off()

file.remove("/mnt/NexData_file/home/lzl/20260508_XE2511106427-01/contour_plot/ contour_plot_HRAS in each sample.pdf")

######
dir.create("RDS_need")
saveRDS(KC_data,"/mnt/NexData_file/home/lzl/20260508_XE2511106427-01/RDS_need/KC_subtype.rds")
saveRDS(Immnue_data,"/mnt/NexData_file/home/lzl/20260508_XE2511106427-01/RDS_need/Immnue_subtype.rds")
saveRDS(stroma_data,"/mnt/NexData_file/home/lzl/20260508_XE2511106427-01/RDS_need/stroma_subtype.rds")
saveRDS(all_data,"/mnt/NexData_file/home/lzl/20260508_XE2511106427-01/RDS_need/allcell_subtype.rds")
saveRDS(cells_position_ls,"/mnt/NexData_file/home/lzl/20260508_XE2511106427-01/RDS_need/cells_position_ls.rds")

###
all_data <- readRDS("/mnt/NexData_file/home/lzl/20260508_XE2511106427-01/RDS_need/allcell_subtype.rds")
unique(all_data$subtype)

KC_data <- readRDS("/mnt/NexData_file/home/lzl/20260508_XE2511106427-01/RDS_need/KC_subtype.rds")
unique(KC_data$subtype)
Immnue_data <- readRDS("/mnt/NexData_file/home/lzl/20260508_XE2511106427-01/RDS_need/Immnue_subtype.rds")
unique(Immnue_data$subtype)
stroma_data <- readRDS("/mnt/NexData_file/home/lzl/20260508_XE2511106427-01/RDS_need/stroma_subtype.rds")
unique(stroma_data$subtype)
all_data <- readRDS('/mnt/NexData_file/home/lzl/20260508_XE2511106427-01/RDS_need/allcell_subtype.rds')

####
KC_data_cn_meta <- subset(all_data, annotation %in% c("Bas KC","Grn KC","Spn KC","Ecc Duct","HFE","Seb","Melano","Ecc Gland"))@meta.data
KC_data@meta.data <- KC_data_cn_meta
saveRDS(KC_data,"/mnt/NexData_file/home/lzl/20260508_XE2511106427-01/RDS_need/KC_subtype_CN_NS.rds")

Immnue_data_cn_meta <- subset(all_data, annotation %in% c("Lym","Mac","DC","Mast"))@meta.data
Immnue_data@meta.data <-Immnue_data_cn_meta 
saveRDS(Immnue_data,"/mnt/NexData_file/home/lzl/20260508_XE2511106427-01/RDS_need/Immnue_subtype_CN_NS.rds")

stroma_data_cn_meta <- subset(all_data, annotation %in% c("Fibro","EC","Peri","Schwann","Adipo","LEC"))@meta.data
stroma_data@meta.data <- stroma_data_cn_meta
saveRDS(stroma_data,"/mnt/NexData_file/home/lzl/20260508_XE2511106427-01/RDS_need/stroma_subtype_CN_NS.rds")
