library(Seurat)
library(tidyverse)


# all_data_sc ------------
load("scRNA annotation/Ps_RNA_annotation.RData")

all_data_sc <- test.seu

plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
          "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead", 
          "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
          "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
          "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
          "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
          "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
          "#005995", "#fa625f", "#f17d80")


major_marker_ls <- c(
  # "Basal_cell" = "KRT5", 
  "Basal_cell" = "KRT14",
  "Basal_cell" = "COL17A1", 
  "KC" = "KRT5", 
  "KC" = "CDH1", 
  # "Keratinocyte_progenitor_cell" = "ABCG2", 
  "Keratinocyte_progenitor_cell" = "CD34", 
  "Keratinocyte_progenitor_cell" = "TP63", 
  "Spinous_cell" = "KRT1", 
  "Spinous_cell" = "KRT10", 
  # "Spinous_cell" = "KRT5", 
  # "Spinous_cell" = "KRT14", 
  "Sebaceous_gland" = "MGST1", 
  "Sebaceous_gland" = "FASN", 
  "Granular_cell" = "SPINK5", 
  "Granular_cell" = "FLG", 
  # "Granular_cell" = "KRT1", 
  # "Granular_cell" = "KRT10", 
  # "Proliferating_keratinocyte_Mitotic_cell" = "UBE2C", 
  "Proliferating_keratinocyte_Mitotic_cell" = "TOP2A", 
  # "Proliferating_keratinocyte_Mitotic_cell" = "PCNA", 
  "Proliferating_keratinocyte_Mitotic_cell" = "MKI67", 
  "Hair_follicle" = "KRT6B", 
  "Hair_follicle" = "KRT17", 
  # "Hair_follicle" = "CD34", 
  # "Hair_follicle" = "RUNX1", 
  # "Sweat_Gland_Cel" = "SCGB1D2", 
  "Sweat_Gland_Cel" = "SCGB2A2", 
  # "Sweat_Gland_Cel" = "DCD", 
  "Sweat_Gland_Cel" = "AQP5",
  "Endothelial_cell" = "PECAM1", 
  "Endothelial_cell" = "SELE", 
  "Endothelial_cell" = "CLDN5", 
  "Endothelial_cell" = "CDH5", 
  # "Endothelial_cell" = "VWF", 
  "VEC" = "VWF", 
  "VEC" = "CD36", 
  "VEC" = "ACKR1", 
  # "VEC" = "FABP4", 
  "Fibroblast" = "COL1A1", 
  "Fibroblast" = "DCN", 
  "Fibroblast" = "COL1A2", 
  # "Fibroblast" = "SFRP2", 
  # "Fibroblast" = "COL6A5", 
  # "Fibroblast" = "COL18A1", 
  "LEC" = "CCL21", 
  "LEC" = "LYVE1", 
  "Langerhans_cells" = "CD1A", 
  "Langerhans_cells" = "CD207", 
  "Dendritic_cell" = "CD1C", 
  "Dendritic_cell" = "CD80", 
  "Dendritic_cell" = "CD83", 
  "Dendritic_cell" = "CD86", 
  # "Dendritic_cell" = "HLA-A", 
  # "Dendritic_cell" = "HLA-B", 
  # "Dendritic_cell" = "HLA-C", 
  # "Dendritic_cell" = "HLA-DRA", 
  # "Dendritic_cell" = "HLA-DPA1", 
  # "Dendritic_cell" = "HLA-DQA1", 
  "Dendritic_cell" = "LAMP3", 
  "T_cell" = "CD3E", 
  "T_cell" = "CD4", 
  "T_cell" = "CD8A", 
  "NK_cell" = "FCGR3A", 
  "NK_cell" = "CD2", 
  # "NK_cell" = "TNFRSF8", 
  # "NK_cell" = "CD44", 
  "Mast_cell" = "TPSAB1", 
  "Mast_cell" = "TPSB2", 
  "Mast_cell" = "CPA3", 
  "Melanocyte_cell" = "MLANA", 
  "Melanocyte_cell" = "TYRP1", 
  "Melanocyte_cell" = "PMEL", 
  "Pericyte" = "RGS5", 
  "Pericyte" = "FABP4", 
  "vSMC" = "MYL9", 
  "vSMC" = "TPM2", 
  "vSMC" = "RERGL", 
  # "Stem_cell" = "LGR5", 
  "Stem_cell" = "SOX9", 
  # "Stem_cell" = "CD34", 
  "Stem_cell" = "KRT15", 
  "Stem_cell" = "TCF4", 
  "Stem_cell" = "CTNNB1", 
  "Stem_cell" = "LRIG1", 
  "Stem_cell" = "CD200", 
  "Stem_cell" = "KRT19", 
  # "Osteogenic_precursor_like_cell" = "BMPR1B", 
  "Dermal_papilla_cell" = "PROM1", 
  # "Embryonic_stem_cell" = "NANOG", 
  # "Embryonic_stem_cell" = "POU5F1", 
  "Merkel_cell" = "KRT18", 
  # "Merkel_cell" = "KRT20", 
  "Merkel_cell" = "NCAM1", 
  "Pluripotent_stem_cell" = "ANPEP", 
  "Progenitor_cell" = "EPCAM", 
  "Progenitor_cell" = "NCAM2",
  # "Neuron" = "CAMK2A",
  # "Neuron" = "CBLN2",
  "Neuron" = "LDB2",
  # "Neuron" = "GAD1",
  "Neuron" = "LHFPL3",
  "Neuron" = "L1CAM",
  "Neuron" = "TUBB3"#,
  # "Neuron" = "PCDH15"
  # "Progenitor_cell" = "PROM1", 
  # "Treg" = "LAG3", 
  # "Treg" = "ITGA2", 
  # "Treg" = "FOXP3", 
  # "Treg" = "IKZF2", 
  # "Treg" = "NRP1", 
  # "Treg" = "NT5E", 
  # "Treg" = "IL10", 
  # "Treg" = "TGFB1", 
)



all_data_sc$sample <- all_data_sc$orig.ident1
all_data_sc$group <- gsub("_\\d+$", "", all_data_sc$sample)



p <- DimPlot(all_data_sc, group.by = "seurat_clusters", pt.size = 0.1, reduction = "umap", label = T, raster = FALSE, cols = plan) + 
  theme_bw() + #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"),
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"),
    #legend.position = "none", #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  guides(color = guide_legend(override.aes = list(size = 8), ncol = 2)) + 
  labs(title = "Cluster") + # 标题
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
ggsave(p, file = "scRNA annotation/Ps_RNA_annotation/umap.pdf", width = 8, height = 6)
ggsave(p, file = "scRNA annotation/Ps_RNA_annotation/umap.png", width = 8, height = 6, dpi = 300)



p <- DimPlot(all_data_sc, group.by = "celltype", pt.size = 0.1, reduction = "umap", raster = FALSE, cols = plan) + 
  theme_bw() + #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"), 
    axis.title.x = element_text(size = 15, face = "bold"), 
    #legend.position = "none", #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  guides(color = guide_legend(override.aes = list(size = 8), ncol = 1)) + 
  labs(title = "Cell Type") + # 标题
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
ggsave(p, file = "scRNA annotation/Ps_RNA_annotation/umap_annot.pdf", width = 7.5, height = 5)
ggsave(p, file = "scRNA annotation/Ps_RNA_annotation/umap_annot.png", width = 7.5, height = 5, dpi = 300)



p <- DimPlot(all_data_sc, group.by = "sample", pt.size = 0.1, reduction = "umap", raster = FALSE) + 
  theme_bw() + #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"), 
    axis.title.x = element_text(size = 15, face = "bold"), 
    #legend.position = "none", #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  labs(title = "Sample") + # 标题
  guides(color = guide_legend(override.aes = list(size = 8), ncol = 2)) + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
ggsave(p, file = "scRNA annotation/Ps_RNA_annotation/umap_sample.pdf", width = 8, height = 5)
ggsave(p, file = "scRNA annotation/Ps_RNA_annotation/umap_sample.png", width = 8, height = 5, dpi = 300)




p <- DimPlot(all_data_sc, group.by = "celltype", ncol = 4, split.by = "group", pt.size = 0.1, raster = FALSE, reduction = "umap", cols = plan) + 
  theme_bw() + #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"), 
    axis.title.x = element_text(size = 15, face = "bold"), 
    strip.background = element_rect(fill = "white", colour = NA), 
    strip.text.x = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5), 
    #legend.position = "none", #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  labs(title = "") + # 标题
  guides(color = guide_legend(override.aes = list(size = 8), ncol = 1)) + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
ggsave(p, file = "scRNA annotation/Ps_RNA_annotation/umap_group_split.pdf", width = 18, height = 5)
ggsave(p, file = "scRNA annotation/Ps_RNA_annotation/umap_group_split.png", width = 18, height = 5, dpi = 300)




Idents(all_data_sc) <- "celltype"
Idents(all_data_sc) <- factor(Idents(all_data_sc), levels = sort(unique(all_data_sc@meta.data$celltype)))
# gene_list <- read.table("marker.txt", header = T, sep = "\t", quote = "")

p <- VlnPlot(all_data_sc, 
             pt.size = 0, # 设置点的大小，但是如果堆叠，则强制不显示点
             stack = T, # 是否堆叠
             flip = F, # 翻转
             #adjust = 5, # 调整粗细
             fill.by = "ident", 
             cols = plan, 
             features = major_marker_ls) + 
  labs(x = "") + 
  theme(legend.position = 'none', 
        axis.text.x = element_text(size = 12, color = "black", angle = 90, vjust = 0.5),
        strip.text.y = element_text(hjust = 0))
ggsave(p, file = "scRNA annotation/Ps_RNA_annotation/vlnplot_annot.pdf", width = 15, height = 6)
ggsave(p, file = "scRNA annotation/Ps_RNA_annotation/vlnplot_annot.png", width = 15, height = 6, dpi = 300)



## 气泡图
p <- DotPlot(all_data_sc, features = major_marker_ls, dot.scale = 10, cols = c("lightgrey", "red")) + 
  labs(x = "", y = "") + 
  theme_bw() + #改变ggplot2的主题
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    axis.text.y = element_text(size = 12, color = "black"), 
    axis.text.x = element_text(size = 12, color = "black", angle = 45, vjust = 1, hjust = 1)) + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
ggsave(p, file = "scRNA annotation/Ps_RNA_annotation/dotplot_annot.pdf", width = 18, height = 5)
ggsave(p, file = "scRNA annotation/Ps_RNA_annotation/dotplot_annot.png", width = 18, height = 5, dpi=300)


saveRDS(all_data_sc, "scRNA annotation/Ps_RNA_annotation/all_data_sc.rds")







# T_cell ------------
load("scRNA annotation/Ps_T_annotation.RData")

all_data_sc <- test.seu

plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
          "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead", 
          "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
          "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
          "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
          "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
          "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
          "#005995", "#fa625f", "#f17d80")


major_marker_ls <- c(
  # "Basal_cell" = "KRT5", 
  "Basal_cell" = "KRT14",
  "Basal_cell" = "COL17A1", 
  "KC" = "KRT5", 
  "KC" = "CDH1", 
  # "Keratinocyte_progenitor_cell" = "ABCG2", 
  "Keratinocyte_progenitor_cell" = "CD34", 
  "Keratinocyte_progenitor_cell" = "TP63", 
  "Spinous_cell" = "KRT1", 
  "Spinous_cell" = "KRT10", 
  # "Spinous_cell" = "KRT5", 
  # "Spinous_cell" = "KRT14", 
  "Sebaceous_gland" = "MGST1", 
  "Sebaceous_gland" = "FASN", 
  "Granular_cell" = "SPINK5", 
  "Granular_cell" = "FLG", 
  # "Granular_cell" = "KRT1", 
  # "Granular_cell" = "KRT10", 
  # "Proliferating_keratinocyte_Mitotic_cell" = "UBE2C", 
  "Proliferating_keratinocyte_Mitotic_cell" = "TOP2A", 
  # "Proliferating_keratinocyte_Mitotic_cell" = "PCNA", 
  "Proliferating_keratinocyte_Mitotic_cell" = "MKI67", 
  "Hair_follicle" = "KRT6B", 
  "Hair_follicle" = "KRT17", 
  # "Hair_follicle" = "CD34", 
  # "Hair_follicle" = "RUNX1", 
  # "Sweat_Gland_Cel" = "SCGB1D2", 
  "Sweat_Gland_Cel" = "SCGB2A2", 
  # "Sweat_Gland_Cel" = "DCD", 
  "Sweat_Gland_Cel" = "AQP5",
  "Endothelial_cell" = "PECAM1", 
  "Endothelial_cell" = "SELE", 
  "Endothelial_cell" = "CLDN5", 
  "Endothelial_cell" = "CDH5", 
  # "Endothelial_cell" = "VWF", 
  "VEC" = "VWF", 
  "VEC" = "CD36", 
  "VEC" = "ACKR1", 
  # "VEC" = "FABP4", 
  "Fibroblast" = "COL1A1", 
  "Fibroblast" = "DCN", 
  "Fibroblast" = "COL1A2", 
  # "Fibroblast" = "SFRP2", 
  # "Fibroblast" = "COL6A5", 
  # "Fibroblast" = "COL18A1", 
  "LEC" = "CCL21", 
  "LEC" = "LYVE1", 
  "Langerhans_cells" = "CD1A", 
  "Langerhans_cells" = "CD207", 
  "Dendritic_cell" = "CD1C", 
  "Dendritic_cell" = "CD80", 
  "Dendritic_cell" = "CD83", 
  "Dendritic_cell" = "CD86", 
  # "Dendritic_cell" = "HLA-A", 
  # "Dendritic_cell" = "HLA-B", 
  # "Dendritic_cell" = "HLA-C", 
  # "Dendritic_cell" = "HLA-DRA", 
  # "Dendritic_cell" = "HLA-DPA1", 
  # "Dendritic_cell" = "HLA-DQA1", 
  "Dendritic_cell" = "LAMP3", 
  "T_cell" = "CD3E", 
  "T_cell" = "CD4", 
  "T_cell" = "CD8A", 
  "NK_cell" = "FCGR3A", 
  "NK_cell" = "CD2", 
  # "NK_cell" = "TNFRSF8", 
  # "NK_cell" = "CD44", 
  "Mast_cell" = "TPSAB1", 
  "Mast_cell" = "TPSB2", 
  "Mast_cell" = "CPA3", 
  "Melanocyte_cell" = "MLANA", 
  "Melanocyte_cell" = "TYRP1", 
  "Melanocyte_cell" = "PMEL", 
  "Pericyte" = "RGS5", 
  "Pericyte" = "FABP4", 
  "vSMC" = "MYL9", 
  "vSMC" = "TPM2", 
  "vSMC" = "RERGL", 
  # "Stem_cell" = "LGR5", 
  "Stem_cell" = "SOX9", 
  # "Stem_cell" = "CD34", 
  "Stem_cell" = "KRT15", 
  "Stem_cell" = "TCF4", 
  "Stem_cell" = "CTNNB1", 
  "Stem_cell" = "LRIG1", 
  "Stem_cell" = "CD200", 
  "Stem_cell" = "KRT19", 
  # "Osteogenic_precursor_like_cell" = "BMPR1B", 
  "Dermal_papilla_cell" = "PROM1", 
  # "Embryonic_stem_cell" = "NANOG", 
  # "Embryonic_stem_cell" = "POU5F1", 
  "Merkel_cell" = "KRT18", 
  # "Merkel_cell" = "KRT20", 
  "Merkel_cell" = "NCAM1", 
  "Pluripotent_stem_cell" = "ANPEP", 
  "Progenitor_cell" = "EPCAM", 
  "Progenitor_cell" = "NCAM2",
  # "Neuron" = "CAMK2A",
  # "Neuron" = "CBLN2",
  "Neuron" = "LDB2",
  # "Neuron" = "GAD1",
  "Neuron" = "LHFPL3",
  "Neuron" = "L1CAM",
  "Neuron" = "TUBB3"#,
  # "Neuron" = "PCDH15"
  # "Progenitor_cell" = "PROM1", 
  # "Treg" = "LAG3", 
  # "Treg" = "ITGA2", 
  # "Treg" = "FOXP3", 
  # "Treg" = "IKZF2", 
  # "Treg" = "NRP1", 
  # "Treg" = "NT5E", 
  # "Treg" = "IL10", 
  # "Treg" = "TGFB1", 
)



all_data_sc$sample <- all_data_sc$orig.ident1
all_data_sc$group <- gsub("_\\d+$", "", all_data_sc$sample)



p <- DimPlot(all_data_sc, group.by = "seurat_clusters", pt.size = 0.1, reduction = "umap", label = T, raster = FALSE, cols = plan) + 
  theme_bw() + #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"),
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"),
    #legend.position = "none", #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  guides(color = guide_legend(override.aes = list(size = 8), ncol = 2)) + 
  labs(title = "Cluster") + # 标题
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
ggsave(p, file = "scRNA annotation/Ps_T_annotation/umap.pdf", width = 8, height = 6)
ggsave(p, file = "scRNA annotation/Ps_T_annotation/umap.png", width = 8, height = 6, dpi = 300)



p <- DimPlot(all_data_sc, group.by = "celltype", pt.size = 0.1, label = T, reduction = "umap", raster = FALSE, cols = plan) + 
  theme_bw() + #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"), 
    axis.title.x = element_text(size = 15, face = "bold"), 
    #legend.position = "none", #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  guides(color = guide_legend(override.aes = list(size = 8), ncol = 1)) + 
  labs(title = "Cell Type") + # 标题
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
ggsave(p, file = "scRNA annotation/Ps_T_annotation/umap_annot.pdf", width = 7.5, height = 5)
ggsave(p, file = "scRNA annotation/Ps_T_annotation/umap_annot.png", width = 7.5, height = 5, dpi = 300)



p <- DimPlot(all_data_sc, group.by = "sample", pt.size = 0.1, reduction = "umap", raster = FALSE) + 
  theme_bw() + #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"), 
    axis.title.x = element_text(size = 15, face = "bold"), 
    #legend.position = "none", #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  labs(title = "Sample") + # 标题
  guides(color = guide_legend(override.aes = list(size = 8), ncol = 2)) + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
ggsave(p, file = "scRNA annotation/Ps_T_annotation/umap_sample.pdf", width = 8, height = 5)
ggsave(p, file = "scRNA annotation/Ps_T_annotation/umap_sample.png", width = 8, height = 5, dpi = 300)




p <- DimPlot(all_data_sc, group.by = "celltype", ncol = 4, split.by = "group", pt.size = 0.1, raster = FALSE, reduction = "umap", cols = plan) + 
  theme_bw() + #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"), 
    axis.title.x = element_text(size = 15, face = "bold"), 
    strip.background = element_rect(fill = "white", colour = NA), 
    strip.text.x = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5), 
    #legend.position = "none", #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  labs(title = "") + # 标题
  guides(color = guide_legend(override.aes = list(size = 8), ncol = 1)) + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
ggsave(p, file = "scRNA annotation/Ps_T_annotation/umap_group_split.pdf", width = 18, height = 5)
ggsave(p, file = "scRNA annotation/Ps_T_annotation/umap_group_split.png", width = 18, height = 5, dpi = 300)




Idents(all_data_sc) <- "celltype"
Idents(all_data_sc) <- factor(Idents(all_data_sc), levels = sort(unique(all_data_sc@meta.data$celltype)))
# gene_list <- read.table("marker.txt", header = T, sep = "\t", quote = "")

p <- VlnPlot(all_data_sc, 
             pt.size = 0, # 设置点的大小，但是如果堆叠，则强制不显示点
             stack = T, # 是否堆叠
             flip = F, # 翻转
             #adjust = 5, # 调整粗细
             fill.by = "ident", 
             cols = plan, 
             features = major_marker_ls) + 
  labs(x = "") + 
  theme(legend.position = 'none', 
        axis.text.x = element_text(size = 12, color = "black", angle = 90, vjust = 0.5),
        strip.text.y = element_text(hjust = 0))
ggsave(p, file = "scRNA annotation/Ps_T_annotation/vlnplot_annot.pdf", width = 15, height = 6)
ggsave(p, file = "scRNA annotation/Ps_T_annotation/vlnplot_annot.png", width = 15, height = 6, dpi = 300)



## 气泡图
p <- DotPlot(all_data_sc, features = major_marker_ls, dot.scale = 10, cols = c("lightgrey", "red")) + 
  labs(x = "", y = "") + 
  theme_bw() + #改变ggplot2的主题
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    axis.text.y = element_text(size = 12, color = "black"), 
    axis.text.x = element_text(size = 12, color = "black", angle = 45, vjust = 1, hjust = 1)) + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
ggsave(p, file = "scRNA annotation/Ps_T_annotation/dotplot_annot.pdf", width = 18, height = 5)
ggsave(p, file = "scRNA annotation/Ps_T_annotation/dotplot_annot.png", width = 18, height = 5, dpi=300)


saveRDS(all_data_sc, "scRNA annotation/Ps_T_annotation/all_data_sc.rds")









# Endothelial ------------
load("scRNA annotation/Ps_Endothelial_annotation.RData")

all_data_sc <- test.seu

plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
          "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead", 
          "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
          "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
          "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
          "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
          "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
          "#005995", "#fa625f", "#f17d80")


major_marker_ls <- c(
  # "Basal_cell" = "KRT5", 
  "Basal_cell" = "KRT14",
  "Basal_cell" = "COL17A1", 
  "KC" = "KRT5", 
  "KC" = "CDH1", 
  # "Keratinocyte_progenitor_cell" = "ABCG2", 
  "Keratinocyte_progenitor_cell" = "CD34", 
  "Keratinocyte_progenitor_cell" = "TP63", 
  "Spinous_cell" = "KRT1", 
  "Spinous_cell" = "KRT10", 
  # "Spinous_cell" = "KRT5", 
  # "Spinous_cell" = "KRT14", 
  "Sebaceous_gland" = "MGST1", 
  "Sebaceous_gland" = "FASN", 
  "Granular_cell" = "SPINK5", 
  "Granular_cell" = "FLG", 
  # "Granular_cell" = "KRT1", 
  # "Granular_cell" = "KRT10", 
  # "Proliferating_keratinocyte_Mitotic_cell" = "UBE2C", 
  "Proliferating_keratinocyte_Mitotic_cell" = "TOP2A", 
  # "Proliferating_keratinocyte_Mitotic_cell" = "PCNA", 
  "Proliferating_keratinocyte_Mitotic_cell" = "MKI67", 
  "Hair_follicle" = "KRT6B", 
  "Hair_follicle" = "KRT17", 
  # "Hair_follicle" = "CD34", 
  # "Hair_follicle" = "RUNX1", 
  # "Sweat_Gland_Cel" = "SCGB1D2", 
  "Sweat_Gland_Cel" = "SCGB2A2", 
  # "Sweat_Gland_Cel" = "DCD", 
  "Sweat_Gland_Cel" = "AQP5",
  "Endothelial_cell" = "PECAM1", 
  "Endothelial_cell" = "SELE", 
  "Endothelial_cell" = "CLDN5", 
  "Endothelial_cell" = "CDH5", 
  # "Endothelial_cell" = "VWF", 
  "VEC" = "VWF", 
  "VEC" = "CD36", 
  "VEC" = "ACKR1", 
  # "VEC" = "FABP4", 
  "Fibroblast" = "COL1A1", 
  "Fibroblast" = "DCN", 
  "Fibroblast" = "COL1A2", 
  # "Fibroblast" = "SFRP2", 
  # "Fibroblast" = "COL6A5", 
  # "Fibroblast" = "COL18A1", 
  "LEC" = "CCL21", 
  "LEC" = "LYVE1", 
  "Langerhans_cells" = "CD1A", 
  "Langerhans_cells" = "CD207", 
  "Dendritic_cell" = "CD1C", 
  "Dendritic_cell" = "CD80", 
  "Dendritic_cell" = "CD83", 
  "Dendritic_cell" = "CD86", 
  # "Dendritic_cell" = "HLA-A", 
  # "Dendritic_cell" = "HLA-B", 
  # "Dendritic_cell" = "HLA-C", 
  # "Dendritic_cell" = "HLA-DRA", 
  # "Dendritic_cell" = "HLA-DPA1", 
  # "Dendritic_cell" = "HLA-DQA1", 
  "Dendritic_cell" = "LAMP3", 
  "T_cell" = "CD3E", 
  "T_cell" = "CD4", 
  "T_cell" = "CD8A", 
  "NK_cell" = "FCGR3A", 
  "NK_cell" = "CD2", 
  # "NK_cell" = "TNFRSF8", 
  # "NK_cell" = "CD44", 
  "Mast_cell" = "TPSAB1", 
  "Mast_cell" = "TPSB2", 
  "Mast_cell" = "CPA3", 
  "Melanocyte_cell" = "MLANA", 
  "Melanocyte_cell" = "TYRP1", 
  "Melanocyte_cell" = "PMEL", 
  "Pericyte" = "RGS5", 
  "Pericyte" = "FABP4", 
  "vSMC" = "MYL9", 
  "vSMC" = "TPM2", 
  "vSMC" = "RERGL", 
  # "Stem_cell" = "LGR5", 
  "Stem_cell" = "SOX9", 
  # "Stem_cell" = "CD34", 
  "Stem_cell" = "KRT15", 
  "Stem_cell" = "TCF4", 
  "Stem_cell" = "CTNNB1", 
  "Stem_cell" = "LRIG1", 
  "Stem_cell" = "CD200", 
  "Stem_cell" = "KRT19", 
  # "Osteogenic_precursor_like_cell" = "BMPR1B", 
  "Dermal_papilla_cell" = "PROM1", 
  # "Embryonic_stem_cell" = "NANOG", 
  # "Embryonic_stem_cell" = "POU5F1", 
  "Merkel_cell" = "KRT18", 
  # "Merkel_cell" = "KRT20", 
  "Merkel_cell" = "NCAM1", 
  "Pluripotent_stem_cell" = "ANPEP", 
  "Progenitor_cell" = "EPCAM", 
  "Progenitor_cell" = "NCAM2",
  # "Neuron" = "CAMK2A",
  # "Neuron" = "CBLN2",
  "Neuron" = "LDB2",
  # "Neuron" = "GAD1",
  "Neuron" = "LHFPL3",
  "Neuron" = "L1CAM",
  "Neuron" = "TUBB3"#,
  # "Neuron" = "PCDH15"
  # "Progenitor_cell" = "PROM1", 
  # "Treg" = "LAG3", 
  # "Treg" = "ITGA2", 
  # "Treg" = "FOXP3", 
  # "Treg" = "IKZF2", 
  # "Treg" = "NRP1", 
  # "Treg" = "NT5E", 
  # "Treg" = "IL10", 
  # "Treg" = "TGFB1", 
)



all_data_sc$sample <- all_data_sc$orig.ident1
all_data_sc$group <- gsub("_\\d+$", "", all_data_sc$sample)



p <- DimPlot(all_data_sc, group.by = "seurat_clusters", pt.size = 0.1, reduction = "umap", label = T, raster = FALSE, cols = plan) + 
  theme_bw() + #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"),
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"),
    #legend.position = "none", #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  guides(color = guide_legend(override.aes = list(size = 8), ncol = 2)) + 
  labs(title = "Cluster") + # 标题
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
ggsave(p, file = "scRNA annotation/Ps_Endothelial_annotation/umap.pdf", width = 8, height = 6)
ggsave(p, file = "scRNA annotation/Ps_Endothelial_annotation/umap.png", width = 8, height = 6, dpi = 300)



p <- DimPlot(all_data_sc, group.by = "celltype", pt.size = 0.1, label = T, reduction = "umap", raster = FALSE, cols = plan) + 
  theme_bw() + #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"), 
    axis.title.x = element_text(size = 15, face = "bold"), 
    #legend.position = "none", #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  guides(color = guide_legend(override.aes = list(size = 8), ncol = 1)) + 
  labs(title = "Cell Type") + # 标题
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
ggsave(p, file = "scRNA annotation/Ps_Endothelial_annotation/umap_annot.pdf", width = 7.5, height = 5)
ggsave(p, file = "scRNA annotation/Ps_Endothelial_annotation/umap_annot.png", width = 7.5, height = 5, dpi = 300)



p <- DimPlot(all_data_sc, group.by = "sample", pt.size = 0.1, reduction = "umap", raster = FALSE) + 
  theme_bw() + #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"), 
    axis.title.x = element_text(size = 15, face = "bold"), 
    #legend.position = "none", #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  labs(title = "Sample") + # 标题
  guides(color = guide_legend(override.aes = list(size = 8), ncol = 2)) + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
ggsave(p, file = "scRNA annotation/Ps_Endothelial_annotation/umap_sample.pdf", width = 8, height = 5)
ggsave(p, file = "scRNA annotation/Ps_Endothelial_annotation/umap_sample.png", width = 8, height = 5, dpi = 300)




p <- DimPlot(all_data_sc, group.by = "celltype", ncol = 4, split.by = "group", pt.size = 0.1, raster = FALSE, reduction = "umap", cols = plan) + 
  theme_bw() + #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"), 
    axis.title.x = element_text(size = 15, face = "bold"), 
    strip.background = element_rect(fill = "white", colour = NA), 
    strip.text.x = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5), 
    #legend.position = "none", #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  labs(title = "") + # 标题
  guides(color = guide_legend(override.aes = list(size = 8), ncol = 1)) + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
ggsave(p, file = "scRNA annotation/Ps_Endothelial_annotation/umap_group_split.pdf", width = 18, height = 5)
ggsave(p, file = "scRNA annotation/Ps_Endothelial_annotation/umap_group_split.png", width = 18, height = 5, dpi = 300)




Idents(all_data_sc) <- "celltype"
Idents(all_data_sc) <- factor(Idents(all_data_sc), levels = sort(unique(all_data_sc@meta.data$celltype)))
# gene_list <- read.table("marker.txt", header = T, sep = "\t", quote = "")

p <- VlnPlot(all_data_sc, 
             pt.size = 0, # 设置点的大小，但是如果堆叠，则强制不显示点
             stack = T, # 是否堆叠
             flip = F, # 翻转
             #adjust = 5, # 调整粗细
             fill.by = "ident", 
             cols = plan, 
             features = major_marker_ls) + 
  labs(x = "") + 
  theme(legend.position = 'none', 
        axis.text.x = element_text(size = 12, color = "black", angle = 90, vjust = 0.5),
        strip.text.y = element_text(hjust = 0))
ggsave(p, file = "scRNA annotation/Ps_Endothelial_annotation/vlnplot_annot.pdf", width = 15, height = 6)
ggsave(p, file = "scRNA annotation/Ps_Endothelial_annotation/vlnplot_annot.png", width = 15, height = 6, dpi = 300)



## 气泡图
p <- DotPlot(all_data_sc, features = major_marker_ls, dot.scale = 10, cols = c("lightgrey", "red")) + 
  labs(x = "", y = "") + 
  theme_bw() + #改变ggplot2的主题
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    axis.text.y = element_text(size = 12, color = "black"), 
    axis.text.x = element_text(size = 12, color = "black", angle = 45, vjust = 1, hjust = 1)) + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
ggsave(p, file = "scRNA annotation/Ps_Endothelial_annotation/dotplot_annot.pdf", width = 18, height = 5)
ggsave(p, file = "scRNA annotation/Ps_Endothelial_annotation/dotplot_annot.png", width = 18, height = 5, dpi=300)


saveRDS(all_data_sc, "scRNA annotation/Ps_Endothelial_annotation/all_data_sc.rds")









# Fibroblast ------------
load("scRNA annotation/Ps_Fibroblast_annotation.RData")

all_data_sc <- test.seu

plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
          "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead", 
          "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
          "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
          "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
          "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
          "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
          "#005995", "#fa625f", "#f17d80")


major_marker_ls <- c(
  # "Basal_cell" = "KRT5", 
  "Basal_cell" = "KRT14",
  "Basal_cell" = "COL17A1", 
  "KC" = "KRT5", 
  "KC" = "CDH1", 
  # "Keratinocyte_progenitor_cell" = "ABCG2", 
  "Keratinocyte_progenitor_cell" = "CD34", 
  "Keratinocyte_progenitor_cell" = "TP63", 
  "Spinous_cell" = "KRT1", 
  "Spinous_cell" = "KRT10", 
  # "Spinous_cell" = "KRT5", 
  # "Spinous_cell" = "KRT14", 
  "Sebaceous_gland" = "MGST1", 
  "Sebaceous_gland" = "FASN", 
  "Granular_cell" = "SPINK5", 
  "Granular_cell" = "FLG", 
  # "Granular_cell" = "KRT1", 
  # "Granular_cell" = "KRT10", 
  # "Proliferating_keratinocyte_Mitotic_cell" = "UBE2C", 
  "Proliferating_keratinocyte_Mitotic_cell" = "TOP2A", 
  # "Proliferating_keratinocyte_Mitotic_cell" = "PCNA", 
  "Proliferating_keratinocyte_Mitotic_cell" = "MKI67", 
  "Hair_follicle" = "KRT6B", 
  "Hair_follicle" = "KRT17", 
  # "Hair_follicle" = "CD34", 
  # "Hair_follicle" = "RUNX1", 
  # "Sweat_Gland_Cel" = "SCGB1D2", 
  "Sweat_Gland_Cel" = "SCGB2A2", 
  # "Sweat_Gland_Cel" = "DCD", 
  "Sweat_Gland_Cel" = "AQP5",
  "Endothelial_cell" = "PECAM1", 
  "Endothelial_cell" = "SELE", 
  "Endothelial_cell" = "CLDN5", 
  "Endothelial_cell" = "CDH5", 
  # "Endothelial_cell" = "VWF", 
  "VEC" = "VWF", 
  "VEC" = "CD36", 
  "VEC" = "ACKR1", 
  # "VEC" = "FABP4", 
  "Fibroblast" = "COL1A1", 
  "Fibroblast" = "DCN", 
  "Fibroblast" = "COL1A2", 
  # "Fibroblast" = "SFRP2", 
  # "Fibroblast" = "COL6A5", 
  # "Fibroblast" = "COL18A1", 
  "LEC" = "CCL21", 
  "LEC" = "LYVE1", 
  "Langerhans_cells" = "CD1A", 
  "Langerhans_cells" = "CD207", 
  "Dendritic_cell" = "CD1C", 
  "Dendritic_cell" = "CD80", 
  "Dendritic_cell" = "CD83", 
  "Dendritic_cell" = "CD86", 
  # "Dendritic_cell" = "HLA-A", 
  # "Dendritic_cell" = "HLA-B", 
  # "Dendritic_cell" = "HLA-C", 
  # "Dendritic_cell" = "HLA-DRA", 
  # "Dendritic_cell" = "HLA-DPA1", 
  # "Dendritic_cell" = "HLA-DQA1", 
  "Dendritic_cell" = "LAMP3", 
  "T_cell" = "CD3E", 
  "T_cell" = "CD4", 
  "T_cell" = "CD8A", 
  "NK_cell" = "FCGR3A", 
  "NK_cell" = "CD2", 
  # "NK_cell" = "TNFRSF8", 
  # "NK_cell" = "CD44", 
  "Mast_cell" = "TPSAB1", 
  "Mast_cell" = "TPSB2", 
  "Mast_cell" = "CPA3", 
  "Melanocyte_cell" = "MLANA", 
  "Melanocyte_cell" = "TYRP1", 
  "Melanocyte_cell" = "PMEL", 
  "Pericyte" = "RGS5", 
  "Pericyte" = "FABP4", 
  "vSMC" = "MYL9", 
  "vSMC" = "TPM2", 
  "vSMC" = "RERGL", 
  # "Stem_cell" = "LGR5", 
  "Stem_cell" = "SOX9", 
  # "Stem_cell" = "CD34", 
  "Stem_cell" = "KRT15", 
  "Stem_cell" = "TCF4", 
  "Stem_cell" = "CTNNB1", 
  "Stem_cell" = "LRIG1", 
  "Stem_cell" = "CD200", 
  "Stem_cell" = "KRT19", 
  # "Osteogenic_precursor_like_cell" = "BMPR1B", 
  "Dermal_papilla_cell" = "PROM1", 
  # "Embryonic_stem_cell" = "NANOG", 
  # "Embryonic_stem_cell" = "POU5F1", 
  "Merkel_cell" = "KRT18", 
  # "Merkel_cell" = "KRT20", 
  "Merkel_cell" = "NCAM1", 
  "Pluripotent_stem_cell" = "ANPEP", 
  "Progenitor_cell" = "EPCAM", 
  "Progenitor_cell" = "NCAM2",
  # "Neuron" = "CAMK2A",
  # "Neuron" = "CBLN2",
  "Neuron" = "LDB2",
  # "Neuron" = "GAD1",
  "Neuron" = "LHFPL3",
  "Neuron" = "L1CAM",
  "Neuron" = "TUBB3"#,
  # "Neuron" = "PCDH15"
  # "Progenitor_cell" = "PROM1", 
  # "Treg" = "LAG3", 
  # "Treg" = "ITGA2", 
  # "Treg" = "FOXP3", 
  # "Treg" = "IKZF2", 
  # "Treg" = "NRP1", 
  # "Treg" = "NT5E", 
  # "Treg" = "IL10", 
  # "Treg" = "TGFB1", 
)



all_data_sc$sample <- all_data_sc$orig.ident1
all_data_sc$group <- gsub("_\\d+$", "", all_data_sc$sample)



p <- DimPlot(all_data_sc, group.by = "seurat_clusters", pt.size = 0.1, reduction = "umap", label = T, raster = FALSE, cols = plan) + 
  theme_bw() + #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"),
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"),
    #legend.position = "none", #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  guides(color = guide_legend(override.aes = list(size = 8), ncol = 2)) + 
  labs(title = "Cluster") + # 标题
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
ggsave(p, file = "scRNA annotation/Ps_Fibroblast_annotation/umap.pdf", width = 8, height = 6)
ggsave(p, file = "scRNA annotation/Ps_Fibroblast_annotation/umap.png", width = 8, height = 6, dpi = 300)



p <- DimPlot(all_data_sc, group.by = "celltype", pt.size = 0.1, label = T, reduction = "umap", raster = FALSE, cols = plan) + 
  theme_bw() + #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"), 
    axis.title.x = element_text(size = 15, face = "bold"), 
    #legend.position = "none", #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  guides(color = guide_legend(override.aes = list(size = 8), ncol = 1)) + 
  labs(title = "Cell Type") + # 标题
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
ggsave(p, file = "scRNA annotation/Ps_Fibroblast_annotation/umap_annot.pdf", width = 7.5, height = 5)
ggsave(p, file = "scRNA annotation/Ps_Fibroblast_annotation/umap_annot.png", width = 7.5, height = 5, dpi = 300)



p <- DimPlot(all_data_sc, group.by = "sample", pt.size = 0.1, reduction = "umap", raster = FALSE) + 
  theme_bw() + #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"), 
    axis.title.x = element_text(size = 15, face = "bold"), 
    #legend.position = "none", #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  labs(title = "Sample") + # 标题
  guides(color = guide_legend(override.aes = list(size = 8), ncol = 2)) + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
ggsave(p, file = "scRNA annotation/Ps_Fibroblast_annotation/umap_sample.pdf", width = 8, height = 5)
ggsave(p, file = "scRNA annotation/Ps_Fibroblast_annotation/umap_sample.png", width = 8, height = 5, dpi = 300)




p <- DimPlot(all_data_sc, group.by = "celltype", ncol = 4, split.by = "group", pt.size = 0.1, raster = FALSE, reduction = "umap", cols = plan) + 
  theme_bw() + #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"), 
    axis.title.x = element_text(size = 15, face = "bold"), 
    strip.background = element_rect(fill = "white", colour = NA), 
    strip.text.x = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5), 
    #legend.position = "none", #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  labs(title = "") + # 标题
  guides(color = guide_legend(override.aes = list(size = 8), ncol = 1)) + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
ggsave(p, file = "scRNA annotation/Ps_Fibroblast_annotation/umap_group_split.pdf", width = 18, height = 5)
ggsave(p, file = "scRNA annotation/Ps_Fibroblast_annotation/umap_group_split.png", width = 18, height = 5, dpi = 300)




Idents(all_data_sc) <- "celltype"
Idents(all_data_sc) <- factor(Idents(all_data_sc), levels = sort(unique(all_data_sc@meta.data$celltype)))
# gene_list <- read.table("marker.txt", header = T, sep = "\t", quote = "")

p <- VlnPlot(all_data_sc, 
             pt.size = 0, # 设置点的大小，但是如果堆叠，则强制不显示点
             stack = T, # 是否堆叠
             flip = F, # 翻转
             #adjust = 5, # 调整粗细
             fill.by = "ident", 
             cols = plan, 
             features = major_marker_ls) + 
  labs(x = "") + 
  theme(legend.position = 'none', 
        axis.text.x = element_text(size = 12, color = "black", angle = 90, vjust = 0.5),
        strip.text.y = element_text(hjust = 0))
ggsave(p, file = "scRNA annotation/Ps_Fibroblast_annotation/vlnplot_annot.pdf", width = 15, height = 6)
ggsave(p, file = "scRNA annotation/Ps_Fibroblast_annotation/vlnplot_annot.png", width = 15, height = 6, dpi = 300)



## 气泡图
p <- DotPlot(all_data_sc, features = major_marker_ls, dot.scale = 10, cols = c("lightgrey", "red")) + 
  labs(x = "", y = "") + 
  theme_bw() + #改变ggplot2的主题
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    axis.text.y = element_text(size = 12, color = "black"), 
    axis.text.x = element_text(size = 12, color = "black", angle = 45, vjust = 1, hjust = 1)) + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
ggsave(p, file = "scRNA annotation/Ps_Fibroblast_annotation/dotplot_annot.pdf", width = 18, height = 5)
ggsave(p, file = "scRNA annotation/Ps_Fibroblast_annotation/dotplot_annot.png", width = 18, height = 5, dpi=300)


saveRDS(all_data_sc, "scRNA annotation/Ps_Fibroblast_annotation/all_data_sc.rds")









# Myeloid ------------
load("scRNA annotation/Ps_Myeloid_annotation.RData")

all_data_sc <- test.seu

plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
          "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead", 
          "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
          "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
          "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
          "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
          "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
          "#005995", "#fa625f", "#f17d80")


major_marker_ls <- c(
  # "Basal_cell" = "KRT5", 
  "Basal_cell" = "KRT14",
  "Basal_cell" = "COL17A1", 
  "KC" = "KRT5", 
  "KC" = "CDH1", 
  # "Keratinocyte_progenitor_cell" = "ABCG2", 
  "Keratinocyte_progenitor_cell" = "CD34", 
  "Keratinocyte_progenitor_cell" = "TP63", 
  "Spinous_cell" = "KRT1", 
  "Spinous_cell" = "KRT10", 
  # "Spinous_cell" = "KRT5", 
  # "Spinous_cell" = "KRT14", 
  "Sebaceous_gland" = "MGST1", 
  "Sebaceous_gland" = "FASN", 
  "Granular_cell" = "SPINK5", 
  "Granular_cell" = "FLG", 
  # "Granular_cell" = "KRT1", 
  # "Granular_cell" = "KRT10", 
  # "Proliferating_keratinocyte_Mitotic_cell" = "UBE2C", 
  "Proliferating_keratinocyte_Mitotic_cell" = "TOP2A", 
  # "Proliferating_keratinocyte_Mitotic_cell" = "PCNA", 
  "Proliferating_keratinocyte_Mitotic_cell" = "MKI67", 
  "Hair_follicle" = "KRT6B", 
  "Hair_follicle" = "KRT17", 
  # "Hair_follicle" = "CD34", 
  # "Hair_follicle" = "RUNX1", 
  # "Sweat_Gland_Cel" = "SCGB1D2", 
  "Sweat_Gland_Cel" = "SCGB2A2", 
  # "Sweat_Gland_Cel" = "DCD", 
  "Sweat_Gland_Cel" = "AQP5",
  "Endothelial_cell" = "PECAM1", 
  "Endothelial_cell" = "SELE", 
  "Endothelial_cell" = "CLDN5", 
  "Endothelial_cell" = "CDH5", 
  # "Endothelial_cell" = "VWF", 
  "VEC" = "VWF", 
  "VEC" = "CD36", 
  "VEC" = "ACKR1", 
  # "VEC" = "FABP4", 
  "Fibroblast" = "COL1A1", 
  "Fibroblast" = "DCN", 
  "Fibroblast" = "COL1A2", 
  # "Fibroblast" = "SFRP2", 
  # "Fibroblast" = "COL6A5", 
  # "Fibroblast" = "COL18A1", 
  "LEC" = "CCL21", 
  "LEC" = "LYVE1", 
  "Langerhans_cells" = "CD1A", 
  "Langerhans_cells" = "CD207", 
  "Dendritic_cell" = "CD1C", 
  "Dendritic_cell" = "CD80", 
  "Dendritic_cell" = "CD83", 
  "Dendritic_cell" = "CD86", 
  # "Dendritic_cell" = "HLA-A", 
  # "Dendritic_cell" = "HLA-B", 
  # "Dendritic_cell" = "HLA-C", 
  # "Dendritic_cell" = "HLA-DRA", 
  # "Dendritic_cell" = "HLA-DPA1", 
  # "Dendritic_cell" = "HLA-DQA1", 
  "Dendritic_cell" = "LAMP3", 
  "T_cell" = "CD3E", 
  "T_cell" = "CD4", 
  "T_cell" = "CD8A", 
  "NK_cell" = "FCGR3A", 
  "NK_cell" = "CD2", 
  # "NK_cell" = "TNFRSF8", 
  # "NK_cell" = "CD44", 
  "Mast_cell" = "TPSAB1", 
  "Mast_cell" = "TPSB2", 
  "Mast_cell" = "CPA3", 
  "Melanocyte_cell" = "MLANA", 
  "Melanocyte_cell" = "TYRP1", 
  "Melanocyte_cell" = "PMEL", 
  "Pericyte" = "RGS5", 
  "Pericyte" = "FABP4", 
  "vSMC" = "MYL9", 
  "vSMC" = "TPM2", 
  "vSMC" = "RERGL", 
  # "Stem_cell" = "LGR5", 
  "Stem_cell" = "SOX9", 
  # "Stem_cell" = "CD34", 
  "Stem_cell" = "KRT15", 
  "Stem_cell" = "TCF4", 
  "Stem_cell" = "CTNNB1", 
  "Stem_cell" = "LRIG1", 
  "Stem_cell" = "CD200", 
  "Stem_cell" = "KRT19", 
  # "Osteogenic_precursor_like_cell" = "BMPR1B", 
  "Dermal_papilla_cell" = "PROM1", 
  # "Embryonic_stem_cell" = "NANOG", 
  # "Embryonic_stem_cell" = "POU5F1", 
  "Merkel_cell" = "KRT18", 
  # "Merkel_cell" = "KRT20", 
  "Merkel_cell" = "NCAM1", 
  "Pluripotent_stem_cell" = "ANPEP", 
  "Progenitor_cell" = "EPCAM", 
  "Progenitor_cell" = "NCAM2",
  # "Neuron" = "CAMK2A",
  # "Neuron" = "CBLN2",
  "Neuron" = "LDB2",
  # "Neuron" = "GAD1",
  "Neuron" = "LHFPL3",
  "Neuron" = "L1CAM",
  "Neuron" = "TUBB3"#,
  # "Neuron" = "PCDH15"
  # "Progenitor_cell" = "PROM1", 
  # "Treg" = "LAG3", 
  # "Treg" = "ITGA2", 
  # "Treg" = "FOXP3", 
  # "Treg" = "IKZF2", 
  # "Treg" = "NRP1", 
  # "Treg" = "NT5E", 
  # "Treg" = "IL10", 
  # "Treg" = "TGFB1", 
)



all_data_sc$sample <- all_data_sc$orig.ident1
all_data_sc$group <- gsub("_\\d+$", "", all_data_sc$sample)



p <- DimPlot(all_data_sc, group.by = "seurat_clusters", pt.size = 0.1, reduction = "umap", label = T, raster = FALSE, cols = plan) + 
  theme_bw() + #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"),
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"),
    #legend.position = "none", #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  guides(color = guide_legend(override.aes = list(size = 8), ncol = 2)) + 
  labs(title = "Cluster") + # 标题
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
ggsave(p, file = "scRNA annotation/Ps_Myeloid_annotation/umap.pdf", width = 8, height = 6)
ggsave(p, file = "scRNA annotation/Ps_Myeloid_annotation/umap.png", width = 8, height = 6, dpi = 300)



p <- DimPlot(all_data_sc, group.by = "celltype", pt.size = 0.1, label = T, reduction = "umap", raster = FALSE, cols = plan) + 
  theme_bw() + #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"), 
    axis.title.x = element_text(size = 15, face = "bold"), 
    #legend.position = "none", #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  guides(color = guide_legend(override.aes = list(size = 8), ncol = 1)) + 
  labs(title = "Cell Type") + # 标题
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
ggsave(p, file = "scRNA annotation/Ps_Myeloid_annotation/umap_annot.pdf", width = 7.5, height = 5)
ggsave(p, file = "scRNA annotation/Ps_Myeloid_annotation/umap_annot.png", width = 7.5, height = 5, dpi = 300)



p <- DimPlot(all_data_sc, group.by = "sample", pt.size = 0.1, reduction = "umap", raster = FALSE) + 
  theme_bw() + #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"), 
    axis.title.x = element_text(size = 15, face = "bold"), 
    #legend.position = "none", #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  labs(title = "Sample") + # 标题
  guides(color = guide_legend(override.aes = list(size = 8), ncol = 2)) + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
ggsave(p, file = "scRNA annotation/Ps_Myeloid_annotation/umap_sample.pdf", width = 8, height = 5)
ggsave(p, file = "scRNA annotation/Ps_Myeloid_annotation/umap_sample.png", width = 8, height = 5, dpi = 300)




p <- DimPlot(all_data_sc, group.by = "celltype", ncol = 4, split.by = "group", pt.size = 0.1, raster = FALSE, reduction = "umap", cols = plan) + 
  theme_bw() + #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
    legend.text = element_text(size = 15), #图例字体大小和间距
    axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"), 
    axis.title.x = element_text(size = 15, face = "bold"), 
    strip.background = element_rect(fill = "white", colour = NA), 
    strip.text.x = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5), 
    #legend.position = "none", #去掉图例
    #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
    plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
  ) + 
  labs(title = "") + # 标题
  guides(color = guide_legend(override.aes = list(size = 8), ncol = 1)) + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
ggsave(p, file = "scRNA annotation/Ps_Myeloid_annotation/umap_group_split.pdf", width = 18, height = 5)
ggsave(p, file = "scRNA annotation/Ps_Myeloid_annotation/umap_group_split.png", width = 18, height = 5, dpi = 300)




Idents(all_data_sc) <- "celltype"
Idents(all_data_sc) <- factor(Idents(all_data_sc), levels = sort(unique(all_data_sc@meta.data$celltype)))
# gene_list <- read.table("marker.txt", header = T, sep = "\t", quote = "")

p <- VlnPlot(all_data_sc, 
             pt.size = 0, # 设置点的大小，但是如果堆叠，则强制不显示点
             stack = T, # 是否堆叠
             flip = F, # 翻转
             #adjust = 5, # 调整粗细
             fill.by = "ident", 
             cols = plan, 
             features = major_marker_ls) + 
  labs(x = "") + 
  theme(legend.position = 'none', 
        axis.text.x = element_text(size = 12, color = "black", angle = 90, vjust = 0.5),
        strip.text.y = element_text(hjust = 0))
ggsave(p, file = "scRNA annotation/Ps_Myeloid_annotation/vlnplot_annot.pdf", width = 15, height = 6)
ggsave(p, file = "scRNA annotation/Ps_Myeloid_annotation/vlnplot_annot.png", width = 15, height = 6, dpi = 300)



## 气泡图
p <- DotPlot(all_data_sc, features = major_marker_ls, dot.scale = 10, cols = c("lightgrey", "red")) + 
  labs(x = "", y = "") + 
  theme_bw() + #改变ggplot2的主题
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    axis.text.y = element_text(size = 12, color = "black"), 
    axis.text.x = element_text(size = 12, color = "black", angle = 45, vjust = 1, hjust = 1)) + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
ggsave(p, file = "scRNA annotation/Ps_Myeloid_annotation/dotplot_annot.pdf", width = 18, height = 5)
ggsave(p, file = "scRNA annotation/Ps_Myeloid_annotation/dotplot_annot.png", width = 18, height = 5, dpi=300)


saveRDS(all_data_sc, "scRNA annotation/Ps_Myeloid_annotation/all_data_sc.rds")




