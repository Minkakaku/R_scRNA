library(harmony)
library(Seurat)
library(ggplot2)
library(tidyr)
library(patchwork)
library(dplyr)

rm(list = ls())
setwd("D:/Project/AE2507033447-08/")

savePlot <- function(p, output, file_name, unit_width = 2, unit_height = 2, base_w = 6, base_h = 6) {
  library(ggplot2)
  
  # 1. 健壮地提取映射变量名
  get_var_name <- function(mapping) {
    if (is.null(mapping)) return(NULL)
    # 尝试多种方式提取，兼容不同版本的 ggplot2
    name <- tryCatch(rlang::as_name(mapping), error = function(e) NULL)
    if (is.null(name)) name <- as.character(mapping)[2]
    return(name)
  }
  
  x_var <- get_var_name(p$mapping$x)
  y_var <- get_var_name(p$mapping$y)
  
  # 2. 计算分面的行列数
  f_cols <- 1
  f_rows <- 1
  # 检查是否使用了分面
  if (!is.null(p$facet) && !inherits(p$facet, "FacetNull")) {
    facet_build <- ggplot2::ggplot_build(p)$layout$layout
    f_cols <- max(facet_build$COL)
    f_rows <- max(facet_build$ROW)
  }
  
  # 3. 计算单个分面的宽度
  if (!is.null(x_var) && x_var %in% colnames(p$data)) {
    col_data <- p$data[[x_var]]
    if (is.factor(col_data) || is.character(col_data)) {
      n_levels <- length(unique(na.omit(col_data)))
      # 如果水平太多（如 >20），自动压缩单位宽度以防总宽过大
      actual_unit_w <- if(n_levels > 20) unit_width * 0.6 else unit_width
      single_width <- n_levels * actual_unit_w
    } else {
      single_width <- base_w
    }
  } else {
    single_width <- base_w
  }
  
  # 4. 计算单个分面的高度
  if (!is.null(y_var) && y_var %in% colnames(p$data)) {
    col_data <- p$data[[y_var]]
    if (is.factor(col_data) || is.character(col_data)) {
      n_levels <- length(unique(na.omit(col_data)))
      actual_unit_h <- if(n_levels > 20) unit_height * 0.6 else unit_height
      single_height <- n_levels * actual_unit_h
    } else {
      single_height <- base_h
    }
  } else {
    single_height <- base_h
  }
  
  # 5. 总尺寸计算与安全上限控制 (重要！)
  # 140 英寸太大会导致内存溢出，这里限制最大宽度为 45 英寸
  final_width <- min(max(single_width * f_cols, 6), 45)
  final_height <- min(max(single_height * f_rows, 6), 45)
  
  # 6. 保存文件
  if(!dir.exists(output)) dir.create(output, recursive = TRUE)
  file_path <- file.path(output, file_name)
  
  message(glue::glue("Saving Plot: {file_name} | Size: {round(final_width, 1)}x{round(final_height, 1)} in"))
  
  # 保存 PNG (300 DPI)
  tryCatch({
    ggsave(paste0(file_path, ".png"), plot = p, device = "png", 
           width = final_width, height = final_height, limitsize = FALSE, dpi = 300)
  }, error = function(e) {
    warning(glue::glue("PNG save failed for {file_name}: {e$message}"))
  })
  
  # 保存 PDF
  tryCatch({
    ggsave(paste0(file_path, ".pdf"), plot = p, device = "pdf", 
           width = final_width, height = final_height, limitsize = FALSE)
  }, error = function(e) {
    warning(glue::glue("PDF save failed for {file_name}: {e$message}"))
  })
}

colors <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A23","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928","grey",'royalblue')
sce <- readRDS("CD4_T/CD4_T_sub.rds")
features = c("CD3D","CD3E","CD4","CD8A","TRDC","TRGC1","CD68","C1QA","S100A8","EPCAM","NKG7","GNLY")
# 聚类
sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:26,k.param = 20)
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:26 , verbose = TRUE,n.neighbors = 20,min.dist = 0.5)
# sce <- RunTSNE(sce, reduction = "harmony", dims = 1:26, verbose = TRUE)
sce <- FindClusters(sce, resolution = 0.8)
p <- FeaturePlot(sce,features = features) 
p2 <- DotPlot(sce,features = features,group.by = "RNA_snn_res.0.8") # 这里的cluster 9,10,11是CD8T,需要删掉
ggsave(
  p,width = 16,height = 11,
  filename = "CD4_T/CD4_8_qc_feature.pdf"
)
ggsave(
  p2,width = 10,height = 11,
  filename = "CD4_T/CD4_8_qc_dot.pdf"
)
p1 <- DimPlot(sce,label = T,reduction = 'umap',group.by = 'RNA_snn_res.0.8')+scale_color_discrete(palette  = colors)
ggsave(
  p1,width = 8,height = 8,
  filename = "CD4_T/CD4_8_qc_umap.pdf"
)
# > table(sce$RNA_snn_res.0.8)
# 
# 0    1    2    3    4    5    6    7    8    9   10   11   12 
# 6677 6377 5635 5587 5062 4721 4102 2906 2539 1718 1418 1050  924 

# 这个里面簇7 是CD4 ,很少量
# sce_8 <- readRDS("CD8_T/CD8_T_sub.rds")
# p <- FeaturePlot(sce_8,features = features)
# p2 <- DotPlot(sce_8,features = features,group.by = "RNA_snn_res.0.5") 
# p1 <- DimPlot(sce_8,label =T)+scale_color_discrete(palette  = colors)
# ggsave(
#   p1,width = 8,height = 8,
#   filename = "CD8_T/CD8_4_qc_umap.pdf"
# )
# ggsave(
#   p2,width = 10,height = 11,
#   filename = "CD8_T/CD8_4_qc_dot.pdf"
# )
# ggsave(
#   p,width = 16,height = 11,
#   filename = "CD8_T/CD8_4_qc_feature.pdf"
# )
# table(sce_8$RNA_snn_res.0.5)
# 0     1     2     3     4     5     6     7 
# 13247 12969  7447  6271  6254  5144  2541  1650 
outdir1 <- "D:/Project/AE2507033447-08/CD4_T/sub_CD4/"
dir.create(outdir1)
sub = subset(sce ,subset = RNA_snn_res.0.8 %in% c(setdiff(unique(sce$RNA_snn_res.0.8),c("11","9","10"))))
DefaultAssay(sub) <- "RNA"
sub <- NormalizeData(sub, verbose = FALSE)
sub <- FindVariableFeatures(sub, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
sub <- ScaleData(sub, verbose = FALSE)
sub <- RunPCA(sub, npcs = 50, verbose = FALSE)
p_elbow <- ElbowPlot(sub, ndims = 50)
ggsave(
  filename = file.path(outdir1, "T_sub_pca_elbowplot.pdf"),
  plot = p_elbow,
  width = 6,
  height = 4
)
ggsave(
  filename = file.path(outdir1, "T_sub_pca_elbowplot.png"),
  plot = p_elbow,
  width = 6,
  height = 4,
  dpi = 300
)

pca_sample <- DimPlot(sub,reduction='pca',group.by = 'sample')
ggsave(
  filename = file.path(outdir1, "T_sub_pca_sample.pdf"),
  plot = pca_sample,
  width = 6,
  height = 4
)
ggsave(
  filename = file.path(outdir1, "T_sub_pca_sample.png"),
  plot = pca_sample,
  width = 6,
  height = 4,
  dpi = 300
)
# Harmony去批次
sub <- RunHarmony(sub, group.by.vars = "sample")
p_elbow <- ElbowPlot(sub, ndims = 50,reduction = 'harmony')
ggsave(
  filename = file.path(outdir1, "T_sub_harmony_elbowplot.pdf"),
  plot = p_elbow,
  width = 6,
  height = 4
)
ggsave(
  filename = file.path(outdir1, "T_sub_harmony_elbowplot.png"),
  plot = p_elbow,
  width = 6,
  height = 4,
  dpi = 300
)

# 聚类
sub <- FindNeighbors(sub, reduction = "harmony", dims = 1:19,k.param = 20)
sub <- RunUMAP(sub, reduction = "harmony", dims = 1:19 , verbose = TRUE,n.neighbors = 20,min.dist = 0.5)
# sub <- RunTSNE(sub, reduction = "harmony", dims = 1:19, verbose = TRUE)
sub <- FindClusters(sub, resolution = 1.0)
annof <- c(
  # Treg
  "FOXP3","CTLA4","TIGIT","IL2RA",
  # Th1 
  # "IFNG","ID2",'TBX21', 'TNF',
  # Th17
  "IL23R","RORC","CCR6","IL17F","RORA","ANXA1","CD40LG",
  # Tfh
  "CXCR5","PDCD1","TNFSF8","IL6ST",
  # Tn
  "CCR7","SELL","TCF7","LEF1",
  # Tem
  "CCL5","GZMK","KLRB1",
  # Tcm
  "IL7R","BACH2","ICAM2",
  # Tex
  "TNFRSF18","RGS1","TOX","HAVCR2"
  # Temra
  # "GZMH","PRF1","NKG7"
)

dot_data <- DotPlot(
  object = sub,
  features = annof,
  group.by = "RNA_snn_res.1",
  cols = c("lightgrey", "blue"), # 基础配色，后续会被覆盖
  dot.scale = 6
) + theme(axis.text.x = element_text(angle = 45, hjust = 1))


dot_data

Allmark <- FindAllMarkers(sub,logfc.threshold = 0.5,min.pct = 0.3,group.by = "RNA_snn_res.1",only.pos = T)
write.csv(Allmark,"CD4_T/sub_CD4/ALL_markers_0.5_0.3.csv",quote = F,row.names = F)
# 0,2,7,10,11,12Treg
# 3 Tn
# 1,4 Tem
# 6,9 Tcm
# 13 Th17
# 5 Tfh
# 8 Tex

# 根据结果注释细胞
result <- cbind(seq(0,13),c("Treg",
                            "CD4_Tem",
                            "Treg",
                            "CD4_Tn",
                            "CD4_Tem",
                            "CD4_Tfh",
                            "CD4_Tcm",
                            "Treg",
                            "CD4_Tex",
                            "CD4_Tcm",
                            "Treg",
                            "Treg",
                            "Treg",
                            "CD4_Th17"
)) %>% as.data.frame()

# 选择结果写入seurat对象
new.cluster.ids <- result$V2
sub$RNA_snn_res.1 <-  factor(sub$RNA_snn_res.1,levels = as.character(c(0:13)))
names(new.cluster.ids) <- levels(sub$RNA_snn_res.1)
Idents(sub) <- "RNA_snn_res.1"
sce2 <- RenameIdents(sub, new.cluster.ids)
sce2@meta.data[["cell_type"]] <- sce2@active.ident

p1 <- DimPlot(sce2, reduction = "umap",group.by = "RNA_snn_res.1",
              label = TRUE, pt.size = 0.5) + scale_color_discrete(palette = colors)
p2 <- DimPlot(sce2, reduction = "umap",group.by = "cell_type",
              label = TRUE, pt.size = 0.5) + scale_color_discrete(palette = colors)

p1+p2

saveRDS(sce2 ,paste0(outdir1, "CD4T_sub_harmony_umap_tsne_cell_type.rds"))


p_d5 <- DimPlot(sce2,
                reduction = 'umap',
                group.by = "cell_type",
                label = T
) + scale_color_discrete(palette = colors)

p_d6 <- DimPlot(sce2,
                reduction = 'umap',
                group.by = "cell_type",
)+ scale_color_discrete(palette = colors)


ggsave(
  filename = file.path(outdir1, "CD4T_sub_umap_celltype_label.pdf"),
  plot = p_d5,
  width = 10,
  height = 7,
)
ggsave(
  filename = file.path(outdir1, "CD4T_sub_umap_celltype_label.png"),
  plot = p_d5,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  filename = file.path(outdir1, "CD4T_sub_umap_celltype.pdf"),
  plot = p_d6,
  width = 10,
  height = 7,
)
ggsave(
  filename = file.path(outdir1, "CD4T_sub_umap_celltype.png"),
  plot = p_d6,
  width = 10,
  height = 7,
  dpi = 300
)


p_d5 <- DimPlot(sce2,
                reduction = 'umap',
                group.by = "cell_type",
                label = T,split.by = "group",ncol = 2,
) + scale_color_discrete(palette = colors)

p_d6 <- DimPlot(sce2,
                reduction = 'umap',
                group.by = "cell_type",split.by = "group",ncol = 2,
)+ scale_color_discrete(palette = colors)


ggsave(
  filename = file.path(outdir1, "CD4T_sub_umap_celltype_label_split.pdf"),
  plot = p_d5,
  width = 16,
  height = 12,
)
ggsave(
  filename = file.path(outdir1, "CD4T_sub_umap_celltype_label_split.png"),
  plot = p_d5,
  width = 16,
  height = 12,
  dpi = 300
)

ggsave(
  filename = file.path(outdir1, "CD4T_sub_umap_celltype_split.pdf"),
  plot = p_d6,
  width = 16,
  height = 12,
)
ggsave(
  filename = file.path(outdir1, "CD4T_sub_umap_celltype_split.png"),
  plot = p_d6,
  width = 16,
  height = 12,
  dpi = 300
)

# 比例图
data <- sce2@meta.data %>% as.data.frame() 
prop = data %>% select("cell_type","group") %>% table() %>% prop.table() %>%
  as.data.frame() %>% rename_with(~c("Celltype","Group", "Proportion"))
prop = prop %>% group_by(Group) %>% mutate(Proportion = Proportion / sum(Proportion))
prop1 = prop %>% mutate(Proportion = Proportion*100)
write.table(prop1,paste0(outdir1,"./celltype_proportion.xls"),sep = "\t",quote = F,row.names = F)

p1 = ggplot(prop,aes(x=Group, y=Proportion, fill=Celltype))+
  geom_bar(stat="identity",position="fill")+
  theme_classic()+
  scale_fill_manual(values=colors)+
  ggtitle("Proportional Barplot of celltype")+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        axis.text.x=element_text(color='black', angle=45, hjust=1))

p2 <- ggplot(prop,aes(x=Celltype, y=Proportion, fill=Group))+
  geom_bar(stat="identity",position="fill")+
  theme_classic()+
  scale_fill_manual(values=colors)+
  ggtitle("Proportional Barplot of group")+
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
        axis.text.x=element_text(color='black', angle=90, hjust=1,vjust = 0.5)) +
  guides(fill=guide_legend(title=NULL))+geom_hline(yintercept  = c(0.25, 0.5, 0.75), linetype = "dashed", color = "blue", size = 1)

savePlot(p1,outdir1,"celltype_barplot_proportion",unit_width = 2, unit_height = 2)
savePlot(p2,outdir1,"group_barplot_proportion",unit_width = 1, unit_height = 2)

df3 = data %>% select(sample,'cell_type') %>% rename_with(~c("sample",'celltype'))
df3 = df3 %>% group_by(sample,celltype) %>% table() %>% prop.table() %>% as.data.frame() %>% rename_with(~c("sample","celltype","Proportion"))
df3 = df3 %>% group_by(sample) %>% mutate(Proportion = Proportion / sum(Proportion))
p1 = ggplot(df3,aes(x=sample, y=Proportion, fill=celltype))+
  geom_bar(stat="identity",position="fill")+
  theme_classic()+
  scale_fill_manual(values=colors)+
  ggtitle("Proportional Barplot of celltype")+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        axis.text.x=element_text(color='black', angle=45, hjust=1))
savePlot(p1,outdir1,"celltype_barplot_proportion_sample",unit_width = 1, unit_height = 2)


Anno_final <- c(
  # Treg
  "FOXP3","CTLA4","TIGIT","IL2RA",
  # Tem
  "CCL5","GZMK","KLRB1",
  # Th17
  "IL23R","RORC","CCR6","IL17A","RORA",
  # Tfh
  "PDCD1","TNFSF8","IL6ST","TOX",
  # Tn
  "CCR7","SELL","TCF7","LEF1",
  # Tcm
  "IL7R","BACH2","KLF2",
  # Tex
  "TNFRSF18","RGS1","CXCL13","HAVCR2"
  # Temra
  # "GZMH","PRF1","NKG7"
)
dot_data <- DotPlot(
  object = sub,
  features =Anno_final,
  group.by = "RNA_snn_res.1",
  cols = c("lightgrey", "blue"), 
  dot.scale = 6
) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(outdir1, "markers_dot.pdf"),
  plot = dot_data,
  width = 16,
  height = 12,
)
ggsave(
  filename = file.path(outdir1, "markers_dot.png"),
  plot = dot_data,
  width = 16,
  height = 12,
  dpi = 300
)

dot_data <- DotPlot(
  object = sce2,
  features =Anno_final,
  group.by = "cell_type",
  cols = c("lightgrey", "blue"), 
  dot.scale = 6
) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(outdir1, "anno_dot.pdf"),
  plot = dot_data,
  width = 16,
  height = 8,
)
ggsave(
  filename = file.path(outdir1, "anno_dot.png"),
  plot = dot_data,
  width = 16,
  height = 8,
  dpi = 300
)

pp <- FeaturePlot(sce2,features = Anno_final,ncol = 6) 
ggsave(
  filename = file.path(outdir1, "markers_umap.pdf"),
  plot = pp,
  width = 25,
  height = 21,
)
ggsave(
  filename = file.path(outdir1, "markers_umap.png"),
  plot = pp,
  width = 25,
  height = 21,
  dpi = 300
)


plot_df <- FetchData(sce2, vars = c("PDCD1", "CTLA4", "group"))

# 转换成长格式
plot_long <- plot_df %>%
  pivot_longer(cols = c("PDCD1", "CTLA4"), 
               names_to = "Gene", 
               values_to = "Expression")

# 绘图（分面小提琴图）
pp <- ggplot(plot_long, aes(x = group, y = Expression, fill = group)) +
  geom_violin(trim = FALSE, scale = "width") + scale_fill_discrete(palette = colors)+
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 2) +
  theme_minimal() +
  labs(title = "Expression of PDCD1 and CTLA4 by group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(outdir1, "Treatment_vio.pdf"),
  plot = pp,
  width = 10,
  height = 6,
)
ggsave(
  filename = file.path(outdir1, "Treatment_vio.png"),
  plot = pp,
  width = 10,
  height = 6,
  dpi = 300
)

