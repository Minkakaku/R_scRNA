all_data <- readRDS('Plasma_subtype/all_data_Plasma_T_annotation.rds')


mTEC_out_dir <- "mTEC_subtype" 
dir.create(mTEC_out_dir)


library(Seurat)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Hs.eg.db)

#all_data <- readRDS(file.path("R_OUTPUT", "all_data.rds"))

# markers: https://cloud.tencent.com/developer/article/2115987

# mTECs -------------------------------------------------------------------
#这次只提取一种细胞，不用替换
all_data$Plasma_T_annotation <- str_replace(all_data$Plasma_T_annotation, 'cTEC', 'TEC')
all_data$Plasma_T_annotation <- str_replace(all_data$Plasma_T_annotation, 'mTEC', 'TEC')
#all_data$Plasma_T_annotation <- str_replace(all_data$Plasma_T_annotation, 'Basal', 'TEC')
table(all_data$Plasma_T_annotation)

mTEC <- subset(all_data, Plasma_T_annotation == 'TEC')

## PCA --------
mTEC <- FindVariableFeatures(mTEC)

options(future.globals.maxSize = 50 * 1024^3) # 把并行全局变量上限改成 50GB（足够用）
mTEC <- ScaleData(mTEC, vars.to.regress = 'nCount_RNA')
mTEC <- RunPCA(mTEC)

### pca inspector ------
pdf(file.path(mTEC_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(mTEC, ndims = 50)
dev.off()


## dimension reduction --------
dims <- 13 # 18

library(future)

#options(future.globals.maxSize = 50 * 1024^3)
#options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
#plan(multisession, workers = 10)

mTEC <- RunUMAP(mTEC, dims = 1:dims,return.model = TRUE)

pdf(file.path(mTEC_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(mTEC, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

mTEC <- mTEC %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

mTEC <- RunUMAP(
  mTEC, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(mTEC_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(mTEC, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'

#
### cluster -------
mTEC  <-
  FindNeighbors(
    mTEC ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


mTEC <-
  FindClusters(
    mTEC,
    resolution = c(0.1,0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(mTEC_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(mTEC, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(mTEC, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(mTEC, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(mTEC, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(mTEC, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(mTEC, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(mTEC, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(mTEC, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

Idents(mTEC) <- "RNA_snn_res.0.6"

### dotplot -------
mTEC_major_marker_ls <- c(
  # AIRE+ 成熟高功能 mTEC (核心耐受亚群)
  "AIRE+ mTEChigh" = "AIRE",
  "AIRE+ mTEChigh" = "CD80",
  
  # Post-AIRE mTEC (AIRE下调后的终末分化mTEC)
  "Post-AIRE mTEC" = "KRT10",
  "Post-AIRE mTEC" = "KRT1",
  "Post-AIRE mTEC" = "IVL",
  # Tuft cell (簇状细胞，mTEC终末分化特化亚群)
  "Tuft cell" = "POU2F3",
  "Tuft cell" = "CFTR",
  
  # CCL21+ mTEC low (未成熟/趋化引导型mTEC)
  "CCL21+ mTEClow" = "CCL21",
  "CCL19+ mTEClow" = "CCL19"
)



Idents(mTEC)<-'RNA_snn_res.0.6'
mTEC_markers <-
  FindAllMarkers(
    JoinLayers(mTEC),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(mTEC_markers, file = "mTEC_subtype/mTEC_all_markers.csv", row.names = TRUE)

mTEC_top5_markers <- 
  mTEC_markers %>% 
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

names(mTEC_top5_markers) <- sub('.*\\.','',names(mTEC_top5_markers))


# dotplot
# 如果已做其他细胞类型的亚群分析，这里就可以顺便画一下其他类型的细胞亚群基因，比如Fibro_major_marker_ls2，否则只用画T细胞
pdf(file.path(mTEC_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(mTEC, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(mTEC, group.by = "RNA_snn_res.0.6", features = mTEC_major_marker_ls)+RotatedAxis()+theme_90angle
#DotPlot(mTEC, group.by = "RNA_snn_res.0.6", features = mTEC_markers)+RotatedAxis()
DotPlot(mTEC, group.by = "RNA_snn_res.0.6", features = mTEC_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(mTEC, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()



pdf(file.path(mTEC_out_dir, str_glue("dotplot dims{dims} T_major_marker.pdf")), width = 25, height = 6)
DotPlot(mTEC, group.by = "RNA_snn_res.0.6", features = T_major_marker_ls)+RotatedAxis()
DotPlot(mTEC, group.by = "RNA_snn_res.0.6", features = T_major_marker_ls_1)+RotatedAxis()
dev.off()

saveRDS(mTEC ,file = "mTEC_subtype/mTEC_cluster.rds")
mTEC <- readRDS("mTEC_subtype/mTEC_cluster.rds")

some_genes <- c("CD19", "CD38", "CD138")

FeaturePlot(
  mTEC,
  features = some_genes,
  ncol = 3,          # 3个图并排
  pt.size = 0.2,     # 点大小
  label = TRUE,
  repel = TRUE
)

ggsave(
  "mTEC_UMAP_CD19_CD38_CD138.pdf",
  width = 15, height = 4
)


VlnPlot(
  mTEC,
  features = some_genes,
  ncol = 3,
  pt.size = 0      
)


ggsave(
  "mTEC_VlnPlot_CD19_CD38_CD138.pdf",
  width = 15, height = 6
)

table(mTEC_markers$cluster)

## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "CCL19+ mTEClow", # 
  '1'  = "cTEC", #
  '2'  = "CCL19+ mTEClow", #
  '3'  = "cTEC", # 
  '4'  = "cTEC", # 
  '5'  = "cTEC", # 
  '6'  = "cTEC", # 
  '7'  = "CCL19+ mTEClow", # 
  '8'  = "cTEC",
  '9'  = "Cycling TEC",
  '10'  = "AIRE+ mTEChigh",
  '11'  = "CCL19+ mTEClow",
  '12'  = "CCL19+ mTEClow"
)

Idents(mTEC) <- "RNA_snn_res.0.6"

mTEC <- RenameIdents(mTEC, cell_cluster)
mTEC$subtype <- Idents(mTEC)


pdf(file.path(mTEC_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(mTEC, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(mTEC_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 15)
DimPlot(mTEC, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


# 更改原注释
# 1. 先把 all_data$T_annotation 复制到新列 mTEC_T_annotation
all_data$mTEC_Plasma_T_annotation <- as.character(all_data$Plasma_T_annotation)

# 2. 创建映射：细胞名 → 精细T亚型
cell_subtype <- mTEC$subtype
names(cell_subtype) <- colnames(mTEC)

# 3. 只把 "T cell" 替换成精细亚型，其他细胞保持原样
all_data$mTEC_Plasma_T_annotation[colnames(mTEC)] <- as.character(cell_subtype)

# 查看结果（正确）
table(all_data$mTEC_Plasma_T_annotation)
table(all_data$Plasma_T_annotation)
table(mTEC$subtype)
saveRDS(all_data, "mTEC_subtype/all_data_mTEC_Plasma_T_annotation.rds")


#只保留细胞亚群
keep <- c(
  "AIRE+ mTEChigh",
  "CCL19+ mTEClow"
)

mTEC_clean <- subset(mTEC, subtype %in% keep)

mTEC_clean$subtype <- droplevels(mTEC_clean$subtype)
Idents(mTEC_clean) <- mTEC_clean$subtype

table(mTEC_clean$subtype)


mTEC <- mTEC_clean
mTEC <- RunUMAP(mTEC, dims = 1:dims,return.model = TRUE)
pdf(file.path(mTEC_out_dir, str_glue("Dimplot dim{dims} mTEC annotation.pdf")), height = 6, width = 8)
DimPlot(mTEC, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(mTEC_out_dir, str_glue("Dimplot dim{dims} mTEC annotation split.pdf")), height = 6, width = 15)
DimPlot(mTEC, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

saveRDS(mTEC_clean, "mTEC_subtype/mTEC_clean.rds")

mTEC_major_marker_ls_1 <- c(
  # AIRE+ 成熟高功能 mTEC (核心耐受亚群)
  "AIRE+ mTEChigh" = "AIRE",
  "AIRE+ mTEChigh" = "CD80",
  
  "CCL19+ mTEClow" = "CCL19"
)

pdf(file.path(mTEC_out_dir, str_glue("dotplot mTEC subtype major_marker.pdf")), width = 10, height = 6)
DotPlot(mTEC, group.by = "subtype", features = mTEC_major_marker_ls_1 )+RotatedAxis()+theme_90angle
dev.off()
## cell props barplot  ----------------------------------------------------------
mTEC_meta <- mTEC[[]]

mTEC_meta$sample <- factor(mTEC_meta$sample,levels=c(mTEC_meta$sample %>% unique()))

p_col_proj <-
  mTEC_meta %>%
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
  filename = file.path("mTEC_subtype", "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  mTEC_meta %>%
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
  filename = file.path(mTEC_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)

#IFN基因集打分
IFN_geneset_ls <- list(
  IFN_a = c("ADAR","CASP1","CASP8","CCRL2","CD47","CNP","CSF1","CXCL10","CXCL11","EIF2AK2","HERC6","IFI35","IFI44L","IFIH1","IFIT2","IFIT3","IFITM1","IL15","IL4R","IL7","IRF1","IRF7","LAMP3","LPAR6","MOV10","MX1","OAS1","OGFR","PLSCR1","PNPT1","PROCR","RIPK2","RNF31","RSAD2","SAMD9","SELL","SP110","STAT2","TAP1","TDRD7","TRAFD1","TRIM21","TRIM25","TRIM5"),
  IFN_g = c("ADAR","APOL6","BANK1","CASP1","CASP3","CASP7","CASP8","CCL7","CD274","CD38","CD40","CD86","CDKN1A","CFB","CFH","CIITA","CMKLR1","CSF2RB","CXCL10","CXCL11","CXCL9","EIF2AK2","FAS","FCGR1A","FPR1","GCH1","GZMA","HERC6","HIF1A","ICAM1","IDO1","IFI35","IFI44L","IFIH1","IFIT1","IFIT2","IFIT3","IFNAR2","IL10RA","IL15","IL15RA","IL18BP","IL2RB","IL4R","IL6","IL7","IRF1","IRF4","IRF5","IRF7","IRF8","ITGB7","JAK2","KLRK1","LATS2","MVP","MX1","MYD88","NCOA3","NFKB1","NLRC5","NOD1","NUP93","OAS3","OGFR","P2RY14","PDE4B","PFKP","PIM1","PLA2G4A","PLSCR1","PML","PNP","PNPT1","PTGS2","PTPN1","PTPN2","PTPN6","RAPGEF6","RBCK1","RIPK1","RIPK2","RNF31","RSAD2","SAMHD1","SELP","SLAMF7","SOCS1","SOCS3","SP110","SPPL2A","SRI","STAT1","STAT2","STAT3","STAT4","TAP1","TDRD7","TNFSF10","TRAFD1","TRIM21","TRIM25","VAMP5","VCAM1","ZBP1")
)


mTEC_score <- JoinLayers(mTEC) %>%
  AddModuleScore(
    features = IFN_geneset_ls,
    name = names(IFN_geneset_ls)
  )

mTEC_score %>% head
# 5. 绘制小提琴图：IFN-α ------------------------------------------
Idents(mTEC)<-'RNA_snn_res.0.6'
Idents(mTEC_score) <- 'RNA_snn_res.0.6'

p1 <- VlnPlot(
  mTEC_score,
  features = "IFN_a1",  # AddModuleScore自动在末尾加1
  pt.size = 0
) +
  labs(x = "Cluster", y = "IFN-alpha Score")+
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_blank()
  )

# 6. 绘制小提琴图：IFN-γ ------------------------------------------
p2 <- VlnPlot(
  mTEC_score,
  features = "IFN_g2",  # AddModuleScore自动在末尾加1
  pt.size = 0
) +
  labs(x = "Cluster", y = "IFN-gamma Score")+
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_blank()
  )

# 7. 保存图片 -----------------------------------------------------
ggsave(
  plot = p1,
  filename = paste0(mTEC_out_dir, "/IFN_alpha_score_Violin.pdf"),
  width = 8, height = 5
)
ggsave(
  plot = p2,
  filename = paste0(mTEC_out_dir, "/IFN_gamma_score_Violin.pdf"),
  width = 8, height = 5
)

# 组合图（可选）
p_combined <- p1 + p2
ggsave(
  plot = p_combined,
  filename = paste0(mTEC_out_dir, "/IFN_score_Violin.pdf"),
  width = 8, height = 8
)



Idents(mTEC_score) <- 'subtype'

p3 <- VlnPlot(
  mTEC_score,
  features = "IFN_a1",  # AddModuleScore自动在末尾加1
  pt.size = 0
) +
  labs(x = "subtype", y = "IFN-alpha Score")+
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_blank()
  )

# 6. 绘制小提琴图：IFN-γ ------------------------------------------
p4 <- VlnPlot(
  mTEC_score,
  features = "IFN_g2",  # AddModuleScore自动在末尾加1
  pt.size = 0
) +
  labs(x = "subtype", y = "IFN-gamma Score")+
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_blank()
  )



# 组合图（可选）
p_combined <- p3 + p4
ggsave(
  plot = p_combined,
  filename = paste0(mTEC_out_dir, "/IFN_score_subtype_Violin.pdf"),
  width = 8, height = 8
)

# 按照cluster
# 5. 箱线图：IFN-α ------------------------------------------
p1 <- mTEC_score@meta.data %>%
  ggplot(aes(x = RNA_snn_res.0.6, y = IFN_a1, fill = RNA_snn_res.0.6)) +
  geom_boxplot(outlier.shape = NA) +          # 箱线图
  geom_jitter(
    width = 0.2,
    height = 0,
    alpha = 0.5,
    size = 1,
    color = "steelblue"
  ) +
  ggpubr::geom_pwc(tip.length = 0, hide.ns = FALSE, label = 'p.signif') +
  labs(x = "Cluster", y = "IFN-alpha Score") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_blank(),
    legend.position = "none"
  )

# 6. 箱线图：IFN-γ ------------------------------------------
p2 <- mTEC_score@meta.data %>%
  ggplot(aes(x = RNA_snn_res.0.6, y = IFN_g2, fill = RNA_snn_res.0.6)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    width = 0.2,
    height = 0,
    alpha = 0.5,
    size = 1,
    color = "steelblue"
  ) +
  ggpubr::geom_pwc(tip.length = 0, hide.ns = FALSE, label = 'p.signif') +
  labs(x = "Cluster", y = "IFN-gamma Score") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_blank(),
    legend.position = "none"
  )

p_combined <- p1 + p2
ggsave(
  plot = p_combined,
  filename = paste0(mTEC_out_dir, "/IFN_score_cluster_Boxplot.pdf"),
  width = 8, height = 8
)


# 按照subtype
# 5. 箱线图：IFN-α ------------------------------------------
p3 <- mTEC_score@meta.data %>%
  ggplot(aes(x = subtype, y = IFN_a1, fill = subtype)) +
  geom_boxplot(outlier.shape = NA) +          # 箱线图
  geom_jitter(
    width = 0.2,
    height = 0,
    alpha = 0.5,
    size = 1,
    color = "steelblue"
  ) +
  ggpubr::geom_pwc(tip.length = 0, hide.ns = FALSE, label = 'p.signif') +
  labs(x = "Subtype", y = "IFN-alpha Score") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_blank(),
    legend.position = "none"
  )

# 6. 箱线图：IFN-γ ------------------------------------------
p4 <- mTEC_score@meta.data %>%
  ggplot(aes(x = subtype, y = IFN_g2, fill = subtype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    width = 0.2,
    height = 0,
    alpha = 0.5,
    size = 1,
    color = "steelblue"
  ) +
  ggpubr::geom_pwc(tip.length = 0, hide.ns = FALSE, label = 'p.signif') +
  labs(x = "Subtype", y = "IFN-gamma Score") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_blank(),
    legend.position = "none"
  )

p_combined <- p3 + p4
ggsave(
  plot = p_combined,
  filename = paste0(mTEC_out_dir, "/IFN_score_subtype_Boxplot.pdf"),
  width = 8, height = 8
)


# 映射到空间切片
subtype_col <- 'subtype'

cells_position_ls <-
  cbind(
    mTEC[['sample']],
    Embeddings(mTEC[['spatial']])
  ) %>%
  rownames_to_column('cell_id') %>%
  group_by(sample) %>%
  nest() %>%
  deframe()

some_genes <- c("CD19", "CD38", "CD138")

plot_dat <-
  FetchData(
    mTEC, 
    vars = c(subtype_col, some_genes)
  ) %>%
  rownames_to_column('cell_id')


# plot celltype ------------------
all_cells <-  mTEC$subtype %>% unique %>% as.character()
cell_colors <- structure(
  names = all_cells,
  ggsci:: pal_igv()(length(all_cells))
)



p_sp_ct_ls <-
  cells_position_ls %>% # .[1] %>%
  imap(function(coords, sam){
    dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
    
    height <- 5
    width  <- 5
    p <-
      dat2 %>%
      ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
      geom_tile(
        aes(fill = .data[[subtype_col]]),
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


pdf(
  file.path("mTEC_subtype/spatial plot of celltype in each sample.pdf"),
  height = 10,
  width  = 10
)
p_sp_ct_ls %>% walk(print)
dev.off()


# plot each cell in all sample -----------
p_sp_ct_ls_ls <-
  cells_position_ls %>% # .[1] %>%
  imap(function(coords, sam){
    dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
    
    height <- 5
    width  <- 5
    
    p_ls <-
      all_cells %>%
      setNames(., .) %>%
      map(function(x){
        cell_colors[names(cell_colors) != x] <- NA
        
        p <-
          dat2 %>%
          ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
          geom_tile(
            aes(fill = .data[[subtype_col]]),
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
    p_ls
  })

# export by sample
p_sp_ct_ls_ls %>%
  iwalk(function(p_ls, sample){
    pdf(
      file.path(str_glue("mTEC_subtype/spatial plot of celltype in each sample - {sample}.pdf")),
      height = 10,
      width  = 10
    )
    p_ls %>% walk(print)
    dev.off()
  })


