setwd('/mnt/NexData_file/home/zwl/Project/other/AE2511186665-lst')
library(ggplot2) 
library('viridis')
library(dendextend)
library(circlize)

scRNA= readRDS('./merge_raw.rds')

################################################################################
#########上皮
celltype=c('AT1', 'AT1/AT2', 'AT2', 'Basal', 'Ciliated', 'Deuterosomal', 'PNEC', 'Sox9 Epi')
all_markers <- c("Ager", "Pdpn", "Rtkn2", "Sema3b", "Muc1", "Etv5", "Abca3", "Sftpb", "Sftpc", "Lamp3", 
                 "Krt5", "Trp63", "Foxj1", "Rsph1", "Cdhr3", "Cdhr4", "Deup1", "Foxn4", "Cdc20b", "Calca", "Ascl1", "Sox9", "Id2")
subscRNA = subset(scRNA,cell_type %in% celltype)
# colors1 = colorRampPalette(c("#4575B4", "#FFFFBF", "#D73020") )(50)
data_use = Seurat::DotPlot(subscRNA, 
                           features = all_markers,
                           group.by  = c('cell_type'))$data

# PCA 均值计算
mat.pca = subscRNA@reductions$pca@cell.embeddings
mat.pca = mat.pca[rownames(subscRNA@meta.data), ]
mat.pca = as.data.frame(mat.pca)
mat.pca$celltype <- subscRNA@meta.data$cell_type

library(tidyverse)
result <- mat.pca %>%
  group_by(celltype) %>%
  summarise(across(PC_1:PC_40, mean, na.rm = TRUE))
result <- result %>%
  column_to_rownames(var = "celltype")
result <- t(result)

library(ggplot2)
library(ggdendro)

hc <- hclust(dist(t(as.matrix(result))), method = 'ward.D')
hc <- as.dendrogram(hc)
hc_dat <- dendro_data(hc, type = "triangle")
color = c("#d2981a", "#a53e1f", "#457277", "#8f657d", "#42819F", "#86AA7D", "#CBB396",
  "#E2706A", "#6A518C", "#A64C8D", "#5A9F6D", "#D97F2F", "#4E85A6", "#C95A6F", "#84A729")
# plot
p1 = ggplot() +
  geom_segment(data = hc_dat$segments, aes(x = x, y = -y, xend = xend, yend = -yend), size = 1, color = '#8c6bb1') +
  geom_text(data = hc_dat$labels, aes(x = x, y = -y + 5, label = label), hjust = 0, angle = 0) +
  geom_point(data = hc_dat$labels, aes(x = x, y = -y, color = label), size = 6) +
  geom_text(data = hc_dat$labels, aes(x = x, y = -y, label = x)) +
  theme_dendro() +
  NoLegend() +
  coord_flip() +
  expand_limits(y = max(hc_dat$segments$y)) +
  # ylim(-50, 10) 
  scale_color_manual(values = color)
p1
# 提取排序后的细胞类型名称
target_order <- labels(hc) 

# 强制让 DotPlot 的 y 轴顺序等于聚类树的顺序
data_use$id <- factor(data_use$id, levels = target_order)
p <- ggplot(data_use, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled),color='black',shape=21) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
  labs(x='',y=NULL)+
  scale_size(range = c(0,5), limits = c(0, 100), breaks = c(0,25,50,75,100))+
  theme(panel.grid.major = element_line(colour = "grey90",size=0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y=  element_blank(),
        axis.ticks.y =  element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 90,vjust=0.2,hjust=0.95))
p
library(patchwork)
P2 = p1+p+plot_layout(ncol = 2, widths = c(0.5,1))
ggsave(filename = './Epi_marker_list_bubble.png',w=9,h=4,P2,dpi = 150)
ggsave(filename = './Epi_marker_list_bubble.pdf',w=9,h=4,P2,dpi = 150)

################################################################################
#########髓系
celltype <- c("AM", "Basophil", "IM", "Mast", "Neutrophil", "cDC1", "cDC2", "iMON", "maDC")
all_markers <- c("Cybb", "Cd68", "Fcer1a", "Il18rap", "Cxcr2", "Ighm", "Ms4a2", "Gata2", "Kit",
                 "S100a8", "S100a9", "Csf3r", "Flt3", "Clec9a", "Zbtb46", "Ccr2", "Itgb7", "Cd209a",
                 "Tnfrsf1b", "Ctsb", "Lyz2", "Spi1", "Sirpa", "Ccr7", "Lad1")
subscRNA = subset(scRNA,cell_type %in% celltype)
# colors1 = colorRampPalette(c("#4575B4", "#FFFFBF", "#D73020") )(50)
data_use = Seurat::DotPlot(subscRNA, 
                           features = all_markers,
                           group.by  = c('cell_type'))$data

# PCA 均值计算
mat.pca = subscRNA@reductions$pca@cell.embeddings
mat.pca = mat.pca[rownames(subscRNA@meta.data), ]
mat.pca = as.data.frame(mat.pca)
mat.pca$celltype <- subscRNA@meta.data$cell_type

library(tidyverse)
result <- mat.pca %>%
  group_by(celltype) %>%
  summarise(across(PC_1:PC_40, mean, na.rm = TRUE))
result <- result %>%
  column_to_rownames(var = "celltype")
result <- t(result)

library(ggplot2)
library(ggdendro)

hc <- hclust(dist(t(as.matrix(result))), method = 'ward.D')
hc <- as.dendrogram(hc)
hc_dat <- dendro_data(hc, type = "triangle")
color = c("#d2981a", "#a53e1f", "#457277", "#8f657d", "#42819F", "#86AA7D", "#CBB396",
          "#E2706A", "#6A518C", "#A64C8D", "#5A9F6D", "#D97F2F", "#4E85A6", "#C95A6F", "#84A729")
# plot
p1 = ggplot() +
  geom_segment(data = hc_dat$segments, aes(x = x, y = -y, xend = xend, yend = -yend), size = 1, color = '#8c6bb1') +
  geom_text(data = hc_dat$labels, aes(x = x, y = -y + 5, label = label), hjust = 0, angle = 0) +
  geom_point(data = hc_dat$labels, aes(x = x, y = -y, color = label), size = 6) +
  geom_text(data = hc_dat$labels, aes(x = x, y = -y, label = x)) +
  theme_dendro() +
  NoLegend() +
  coord_flip() +
  expand_limits(y = max(hc_dat$segments$y)) +
  # ylim(-50, 10) 
  scale_color_manual(values = color)
p1
# 提取排序后的细胞类型名称
target_order <- labels(hc) 

# 强制让 DotPlot 的 y 轴顺序等于聚类树的顺序
data_use$id <- factor(data_use$id, levels = target_order)
p <- ggplot(data_use, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled),color='black',shape=21) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
  labs(x='',y=NULL)+
  scale_size(range = c(0,5), limits = c(0, 100), breaks = c(0,25,50,75,100))+
  theme(panel.grid.major = element_line(colour = "grey90",size=0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y=  element_blank(),
        axis.ticks.y =  element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 90,vjust=0.2,hjust=0.95))
p
library(patchwork)
P2 = p1+p+plot_layout(ncol = 2, widths = c(0.5,1));P2
ggsave(filename = './Myeloid_marker_list_bubble.png',w=9,h=4,P2,dpi = 150)
ggsave(filename = './Myeloid_marker_list_bubble.pdf',w=9,h=4,P2,dpi = 150)

################################################################################
#########淋巴
celltype=c("B", "CD4_Tcm", "CD4_Tfh", "CD4_Th1", "CD4_Th17", "CD8_Tcm", "CD8_Teff", "ILC", "NK", "Treg")
all_markers <- c("Cd79a", "Ms4a1", "Cd19", "Cd79b", "Il7r", "Tcf7", "Ccr7", "Sell", "Pdcd1", "Icos", "Cxcl13", 
                 "Bcl6", "Cxcr5", "Bhlhe40", "Il12rb2", "Tbx21", "Ifng", "Cxcr3", "Il23r", "Il17a", "Il17f",
                 "Ccr6", "Rorc", "Il7r", "Tcf7", "Ccr7", "Sell", "Gzma", "Prf1", "Gzmb", "Il7r", "Cxcr6", "Ncr1",
                 "Prf1", "Eomes", "Ikzf2", "Tigit", "Il2ra", "Ctla4", "Foxp3")
subscRNA = subset(scRNA,cell_type %in% celltype)
# colors1 = colorRampPalette(c("#4575B4", "#FFFFBF", "#D73020") )(50)
data_use = Seurat::DotPlot(subscRNA, 
                           features = unique(all_markers),
                           group.by  = c('cell_type'))$data

# PCA 均值计算
mat.pca = subscRNA@reductions$pca@cell.embeddings
mat.pca = mat.pca[rownames(subscRNA@meta.data), ]
mat.pca = as.data.frame(mat.pca)
mat.pca$celltype <- subscRNA@meta.data$cell_type

library(tidyverse)
result <- mat.pca %>%
  group_by(celltype) %>%
  summarise(across(PC_1:PC_40, mean, na.rm = TRUE))
result <- result %>%
  column_to_rownames(var = "celltype")
result <- t(result)

library(ggplot2)
library(ggdendro)

hc <- hclust(dist(t(as.matrix(result))), method = 'ward.D')
hc <- as.dendrogram(hc)
hc_dat <- dendro_data(hc, type = "triangle")
color = c("#d2981a", "#a53e1f", "#457277", "#8f657d", "#42819F", "#86AA7D", "#CBB396",
          "#E2706A", "#6A518C", "#A64C8D", "#5A9F6D", "#D97F2F", "#4E85A6", "#C95A6F", "#84A729")
# plot
p1 = ggplot() +
  geom_segment(data = hc_dat$segments, aes(x = x, y = -y, xend = xend, yend = -yend), size = 1, color = '#8c6bb1') +
  geom_text(data = hc_dat$labels, aes(x = x, y = -y + 5, label = label), hjust = 0, angle = 0) +
  geom_point(data = hc_dat$labels, aes(x = x, y = -y, color = label), size = 6) +
  geom_text(data = hc_dat$labels, aes(x = x, y = -y, label = x)) +
  theme_dendro() +
  NoLegend() +
  coord_flip() +
  expand_limits(y = max(hc_dat$segments$y)) +
  # ylim(-50, 10) 
  scale_color_manual(values = color)
p1
# 提取排序后的细胞类型名称
target_order <- labels(hc) 

# 强制让 DotPlot 的 y 轴顺序等于聚类树的顺序
data_use$id <- factor(data_use$id, levels = target_order)
p <- ggplot(data_use, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled),color='black',shape=21) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
  labs(x='',y=NULL)+
  scale_size(range = c(0,5), limits = c(0, 100), breaks = c(0,25,50,75,100))+
  theme(panel.grid.major = element_line(colour = "grey90",size=0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y=  element_blank(),
        axis.ticks.y =  element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 90,vjust=0.2,hjust=0.95))
p
library(patchwork)
P2 = p1+p+plot_layout(ncol = 2, widths = c(0.5,1))
P2
ggsave(filename = './Lymphocyte_marker_list_bubble.png',w=12,h=4,P2,dpi = 150)
ggsave(filename = './Lymphocyte_marker_list_bubble.pdf',w=12,h=4,P2,dpi = 150)


################################################################################
#########基质
celltype= c("AEC", "AF1", "AF2", "ASMC", "CAP1", "CAP2", "Chondrocyte", "EPC", "LEC", "Mesothelial", "Pericyte", "SCMF", "VEC", "VSMC")
all_markers <- c("Dkk2", "Gja5", "Bmx", "Cxcl12", "Tcf21", "Wnt2", "Slc38a5", "Mfap5", "Scara5", "Cdon",
                 "Des", "Lgr6", "Actg2", "Apln", "Il7r", "Gpihbp1", "Kit", "Plvap", "Tek", "Car4", "Apln",
                 "Ednrb", "Col2a1", "Hapln1", "Kit", "Foxf1", "Dll4", "Prox1", "Ccl21a", "Ptx3", "Hoxd9",
                 "Wt1", "Frem2", "Upk3b", "Lamc3", "Higd1b", "Cox4i2", "Notch3", "Tnc", "Acta2", "Tagln", 
                 "Mustn1", "Sdc2", "Loxl2", "Dkk3", "Igf1", "Tpm2", "Myl9", "Ephb4", "Nr2f2", "Cnn1", "Itga7",
                 "Ntrk3", "Hey2")
subscRNA = subset(scRNA,cell_type %in% celltype)
# colors1 = colorRampPalette(c("#4575B4", "#FFFFBF", "#D73020") )(50)
data_use = Seurat::DotPlot(subscRNA, 
                           features = unique(all_markers),
                           group.by  = c('cell_type'))$data

# PCA 均值计算
mat.pca = subscRNA@reductions$pca@cell.embeddings
mat.pca = mat.pca[rownames(subscRNA@meta.data), ]
mat.pca = as.data.frame(mat.pca)
mat.pca$celltype <- subscRNA@meta.data$cell_type

library(tidyverse)
result <- mat.pca %>%
  group_by(celltype) %>%
  summarise(across(PC_1:PC_40, mean, na.rm = TRUE))
result <- result %>%
  column_to_rownames(var = "celltype")
result <- t(result)

library(ggplot2)
library(ggdendro)

hc <- hclust(dist(t(as.matrix(result))), method = 'ward.D')
hc <- as.dendrogram(hc)
hc_dat <- dendro_data(hc, type = "triangle")
color = c("#d2981a", "#a53e1f", "#457277", "#8f657d", "#42819F", "#86AA7D", "#CBB396",
          "#E2706A", "#6A518C", "#A64C8D", "#5A9F6D", "#D97F2F", "#4E85A6", "#C95A6F", "#84A729")
# plot
p1 = ggplot() +
  geom_segment(data = hc_dat$segments, aes(x = x, y = -y, xend = xend, yend = -yend), size = 1, color = '#8c6bb1') +
  geom_text(data = hc_dat$labels, aes(x = x, y = -y + 5, label = label), hjust = 0, angle = 0) +
  geom_point(data = hc_dat$labels, aes(x = x, y = -y, color = label), size = 6) +
  geom_text(data = hc_dat$labels, aes(x = x, y = -y, label = x)) +
  theme_dendro() +
  NoLegend() +
  coord_flip() +
  expand_limits(y = max(hc_dat$segments$y)) +
  # ylim(-50, 10) 
  scale_color_manual(values = color)
p1
# 提取排序后的细胞类型名称
target_order <- labels(hc) 

# 强制让 DotPlot 的 y 轴顺序等于聚类树的顺序
data_use$id <- factor(data_use$id, levels = target_order)
p <- ggplot(data_use, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled),color='black',shape=21) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
  labs(x='',y=NULL)+
  scale_size(range = c(0,5), limits = c(0, 100), breaks = c(0,25,50,75,100))+
  theme(panel.grid.major = element_line(colour = "grey90",size=0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y=  element_blank(),
        axis.ticks.y =  element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 90,vjust=0.2,hjust=0.95))
p
library(patchwork)
P2 = p1+p+plot_layout(ncol = 2, widths = c(0.5,1))
P2
ggsave(filename = './Stromal_marker_list_bubble.png',w=14,h=4,P2,dpi = 150)
ggsave(filename = './Stromal_marker_list_bubble.pdf',w=14,h=4,P2,dpi = 150)

################################################################################
#########ALL
celltype=c('AT1', 'AT1/AT2', 'AT2', 'Basal', 'Ciliated', 'Deuterosomal', 'PNEC', 'Sox9 Epi',
           "AM", "Basophil", "IM", "Mast", "Neutrophil", "cDC1", "cDC2", "iMON", "maDC",
           "B", "CD4_Tcm", "CD4_Tfh", "CD4_Th1", "CD4_Th17", "CD8_Tcm", "CD8_Teff", "ILC", "NK", "Treg",
           "AEC", "AF1", "AF2", "ASMC", "CAP1", "CAP2", "Chondrocyte", "EPC", "LEC", "Mesothelial", "Pericyte", "SCMF", "VEC", "VSMC"
           )
all_markers <- c("Ager", "Pdpn", "Rtkn2", "Sema3b", "Muc1", "Etv5", "Abca3", "Sftpb", "Sftpc", "Lamp3", 
                 "Krt5", "Trp63", "Foxj1", "Rsph1", "Cdhr3", "Cdhr4", "Deup1", "Foxn4", "Cdc20b", "Calca", "Ascl1", "Sox9", "Id2",
                 "Cybb", "Cd68", "Fcer1a", "Il18rap", "Cxcr2", "Ighm", "Ms4a2", "Gata2", "Kit",
                 "S100a8", "S100a9", "Csf3r", "Flt3", "Clec9a", "Zbtb46", "Ccr2", "Itgb7", "Cd209a",
                 "Tnfrsf1b", "Ctsb", "Lyz2", "Spi1", "Sirpa", "Ccr7", "Lad1",
                 "Cd79a", "Ms4a1", "Cd19", "Cd79b", "Il7r", "Tcf7", "Ccr7", "Sell", "Pdcd1", "Icos", "Cxcl13", 
                 "Bcl6", "Cxcr5", "Bhlhe40", "Il12rb2", "Tbx21", "Ifng", "Cxcr3", "Il23r", "Il17a", "Il17f",
                 "Ccr6", "Rorc", "Il7r", "Tcf7", "Ccr7", "Sell", "Gzma", "Prf1", "Gzmb", "Il7r", "Cxcr6", "Ncr1",
                 "Prf1", "Eomes", "Ikzf2", "Tigit", "Il2ra", "Ctla4", "Foxp3",
                 "Dkk2", "Gja5", "Bmx", "Cxcl12", "Tcf21", "Wnt2", "Slc38a5", "Mfap5", "Scara5", "Cdon",
                 "Des", "Lgr6", "Actg2", "Apln", "Il7r", "Gpihbp1", "Kit", "Plvap", "Tek", "Car4", "Apln",
                 "Ednrb", "Col2a1", "Hapln1", "Kit", "Foxf1", "Dll4", "Prox1", "Ccl21a", "Ptx3", "Hoxd9",
                 "Wt1", "Frem2", "Upk3b", "Lamc3", "Higd1b", "Cox4i2", "Notch3", "Tnc", "Acta2", "Tagln", 
                 "Mustn1", "Sdc2", "Loxl2", "Dkk3", "Igf1", "Tpm2", "Myl9", "Ephb4", "Nr2f2", "Cnn1", "Itga7",
                 "Ntrk3", "Hey2")
# subscRNA = subset(scRNA,cell_type %in% celltype)
# colors1 = colorRampPalette(c("#4575B4", "#FFFFBF", "#D73020") )(50)
data_use = Seurat::DotPlot(scRNA, 
                           features = unique(all_markers),
                           group.by  = c('cell_type'))$data

# PCA 均值计算
mat.pca = scRNA@reductions$pca@cell.embeddings
mat.pca = mat.pca[rownames(scRNA@meta.data), ]
mat.pca = as.data.frame(mat.pca)
mat.pca$celltype <- scRNA@meta.data$cell_type

library(tidyverse)
result <- mat.pca %>%
  group_by(celltype) %>%
  summarise(across(PC_1:PC_40, mean, na.rm = TRUE))
result <- result %>%
  column_to_rownames(var = "celltype")
result <- t(result)

library(ggplot2)
library(ggdendro)

hc <- hclust(dist(t(as.matrix(result))), method = 'ward.D')
hc <- as.dendrogram(hc)
hc_dat <- dendro_data(hc, type = "triangle")
color = colorRampPalette(c("#4575B4", "#5A9F6D", "#FEE090", "#D73027"))(43)
# plot
p1 = ggplot() +
  geom_segment(data = hc_dat$segments, aes(x = x, y = -y, xend = xend, yend = -yend), size = 1, color = '#8c6bb1') +
  geom_text(data = hc_dat$labels, aes(x = x, y = -y + 5, label = label), hjust = 0, angle = 0) +
  geom_point(data = hc_dat$labels, aes(x = x, y = -y, color = label), size = 6) +
  geom_text(data = hc_dat$labels, aes(x = x, y = -y, label = x)) +
  theme_dendro() +
  NoLegend() +
  coord_flip() +
  expand_limits(y = max(hc_dat$segments$y)) +
  # ylim(-50, 10) 
  scale_color_manual(values = color)
p1
# 提取排序后的细胞类型名称
target_order <- labels(hc) 

# 强制让 DotPlot 的 y 轴顺序等于聚类树的顺序
data_use$id <- factor(data_use$id, levels = target_order)
p <- ggplot(data_use, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, fill = avg.exp.scaled),color='black',shape=21) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
  labs(x='',y=NULL)+
  scale_size(range = c(0,5), limits = c(0, 100), breaks = c(0,25,50,75,100))+
  theme(panel.grid.major = element_line(colour = "grey90",size=0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y=  element_blank(),
        axis.ticks.y =  element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 90,vjust=0.2,hjust=0.95))
p
library(patchwork)
P2 = p1+p+plot_layout(ncol = 2, widths = c(0.5,1))
P2
ggsave(filename = './ALL_marker_list_bubble_tree.png',w=6,h=12,p1,dpi = 150)
ggsave(filename = './ALL_marker_list_bubble_tree.pdf',w=6,h=12,p1,dpi = 150)

ggsave(filename = './ALL_marker_list_bubble_main.png',w=20,h=12,p,dpi = 150)
ggsave(filename = './ALL_marker_list_bubble_main.pdf',w=20,h=12,p,dpi = 150)
