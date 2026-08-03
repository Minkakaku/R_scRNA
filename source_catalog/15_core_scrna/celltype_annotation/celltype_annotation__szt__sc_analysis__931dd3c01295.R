library(tidyverse)
library(reticulate)

# setting dir -------
sc_output <- "scRNA" # dir.create(img_output)
dir.create(sc_output)

use_python("/home/szt/miniforge3/envs/rpy/bin/python")

reticulate::py_config()

sc <- import("scanpy")
adata <- sc$read("scRNA/scRNA.h5ad")

dup_gene_idx <- duplicated(adata$var$Gene)

counts_x <- adata$X %>% t %>% as.matrix() %>% .[!dup_gene_idx, ]
colnames(counts_x) <- rownames(adata$obs)
# rownames(counts_x) <- adata$var$Gene %>% str_remove("_.*$") %>% .[!dup_gene_idx]
rownames(counts_x) <- adata$var$Gene %>% .[!dup_gene_idx]

scRNA2 <- CreateSeuratObject(
  counts_x,
  meta.data = adata$obs
)

scRNA2 %>%head
scRNA2 %>%Features %>%head


selected_genes <- c("CD19", "CD79A", "CD79B", "PTPRC", "CD3D", "CD3E")
genes <- intersect(rownames(scRNA2), selected_genes)
selected_genesids <- c("ENSG00000177455", "ENSG00000105369", "ENSG00000007312", "ENSG00000262418", "ENSG00000167286", "ENSG00000198851")
selected_gene_ids <- c("CD19-ENSG00000177455", "CD79A-ENSG00000105369", "CD79B-ENSG00000007312", "PTPRC-ENSG00000262418", "CD3D-ENSG00000167286", "CD3E-ENSG00000198851")
genes <- intersect(rownames(scRNA2), selected_gene_ids)


saveRDS(sc_output, "scRNA.rds")



library(Seurat)
library(tidyverse)
library(dplyr)

all_data_sc <- readRDS("scRNA/scRNA.rds")
sc_out_dir <- "scRNA"

all_data_sc <- NormalizeData(all_data_sc, normalization.method = "LogNormalize", scale.factor = 10000)
all_data_sc <- FindVariableFeatures(object = all_data_sc, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(all_data_sc)

all_data_sc <- ScaleData(object = all_data_sc, features = all.genes, vars.to.regress = NULL) 

all_data_sc <- RunPCA(object = all_data_sc)

DimPlot(all_data_sc, reduction = "pca", pt.size = 0.5)
ggsave("scRNA/PCA.pdf")

pdf("scRNA/ElbowPlot.pdf")
ElbowPlot(all_data_sc, ndims = 50)
dev.off()


all_data_sc <- FindNeighbors(all_data_sc, reduction = "pca", dims = 1:18)
all_data_sc <- FindClusters(all_data_sc, resolution = c(0.8, 1.0, 1.2, 1.5), algorithm = 1)
all_data_sc <- RunUMAP(all_data_sc, dims = 1:18, n.neighbors = 30, min.dist = 0.3)


umap_1 <- DimPlot(all_data_sc, group.by = "RNA_snn_res.1.5", reduction = "umap", pt.size = 0.5, label = TRUE)

ggsave("scRNA/umap_res15.pdf", width = 10, height = 8)


saveRDS(all_data_sc, file = "scRNA/all_data_sc_clustered.rds")
all_data_sc <- readRDS("scRNA/all_data_sc_clustered.rds")




VlnPlot(object = all_data_sc, features = c("nFeature_RNA"), ncol = 1, pt.size = 0)



# dotplot ------

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


pdf(file.path(sc_output, str_glue("dotplot scRNA CellType custom_marker.pdf")), width = 40, height = 10)
DotPlot(all_data_sc, group.by = "CellType", features = major_marker_ls, dot.scale = 8) + theme_90angle
dev.off()



pdf(file.path(sc_output, str_glue("dotplot scRNA seurat_clusters major_marker_ls_1 res15.pdf")), width = 40, height = 15)
DotPlot(all_data_sc, group.by = "seurat_clusters", features = major_marker_ls_1, dot.scale = 8) + theme_90angle
dev.off()

pdf(file.path(sc_output, str_glue("dotplot scRNA seurat_clusters major_marker_ls_2 res15.pdf")), width = 40, height = 15)
DotPlot(all_data_sc, group.by = "seurat_clusters", features = major_marker_ls_2, dot.scale = 8) + theme_90angle
dev.off()



## annotation ----
Idents(all_data_sc) <- "RNA_snn_res.1.5"

cell_cluster <- c(
  '0'    = "Mg-TAM",
  '1'    = "GPC", 
  '2'    = "Mg-TAM",
  '3'    = "Ast-Tumor",
  '4'    = "Mo-TAM",
  '5'    = "Mg-TAM",
  '6'    = "Mo-TAM",
  '7'    = "Oligo",
  '8'    = "Ast-Tumor",
  '9'    = "Ast-Tumor",
  '10'   = "Mo-TAM",
  '11'   = "Mg-TAM",
  '12'   = "Ast-Tumor",
  '13'   = "Mg-TAM",
  '14'   = "Mo-TAM",
  '15'   = "Ast-Tumor",
  '16'   = "Oligo",
  '17'   = "Ast-Tumor",
  '18'   = "Endothelial",
  '19'   = "Ast-Tumor",
  '20'   = "Ast-Tumor",
  '21'   = "Mo-TAM",
  '22'   = "Mo-TAM",
  '23'   = "Astrocyte",
  '24'   = "Neuron",
  '25'   = "GPC",
  '26'   = "Astrocyte",
  '27'   = "Oli-Tumor",
  '28'   = "Mg-TAM",
  '29'   = "Oligo",
  '30'   = "B cell",
  '31'   = "Pericyte",
  '32'   = "NK, T-cells",
  '33'   = "GPC",
  '34'   = "Mo-TAM",
  '35'   = "Pericyte",
  '36'   = "Neutrophil",
  '37'   = "Endothelial",
  '38'   = "Endothelial",
  '39'   = "Pericyte",
  '40'   = "Fibroblast",
  '41'   = "Mo-TAM"
)

all_data_sc <- RenameIdents(all_data_sc, cell_cluster)
all_data_sc$annotation <- Idents(all_data_sc)



pdf(file.path(sc_output, str_glue("dotplot scRNA annotation major_marker_ls_1 res15.pdf")), width = 30, height = 10)
DotPlot(all_data_sc, group.by = "annotation", features = major_marker_ls_1, dot.scale = 8) + theme_90angle
dev.off()

pdf(file.path(sc_output, str_glue("dotplot scRNA annotation major_marker_ls_2 res15.pdf")), width = 40, height = 10)
DotPlot(all_data_sc, group.by = "annotation", features = major_marker_ls_2, dot.scale = 8) + theme_90angle
dev.off()







# scRNA dimplot
library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(scales)
library(RColorBrewer)
library(ggsci)

plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
          "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", "#de4307", 
          "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", "#fa6e57", "#f69e53", 
          "#248888", "#999999", "#E7475E", "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", 
          "#FF6D31", "#f1ac9d", "#f06966", "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", 
          "#FF7A5A", "#8ED2C9", "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", 
          "#ffc300", "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
          "#f17d80", "#737495")

# UMAP图
p4 <-DimPlot(all_data_sc, group.by = "annotation", pt.size = 0.1, reduction = "umap", label = T, raster = FALSE, cols = plan, repel = T) +
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
ggsave(p4, file = "scRNA/sc_umap_annotation.pdf", width = 8, height = 6)

# UMAP图
p5 <-DimPlot(all_data_sc, group.by = "CellType", pt.size = 0.1, reduction = "umap", label = T, raster = FALSE, cols = plan, repel = T) +
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
ggsave(p5, file = "scRNA/sc_umap_annot.pdf", width = 8, height = 6)



# p5 <- DimPlot(all_data_sc_subset, group.by = "annotation", pt.size = 0.1, reduction = "umap", label = T, raster = FALSE, cols = plan) +
#   theme_bw() + #改变ggplot2的主题
#   theme( #进一步修改主题
#     panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #去掉背景线
#     legend.text = element_text(size = 15), #图例字体大小和间距
#     axis.text.y = element_text(vjust = 0.5, size = 10, face = "bold"),
#     axis.text.x = element_text(vjust = 0.5, size = 10, face = "bold"),
#     axis.title.y = element_text(size = 15, face = "bold"),
#     axis.title.x = element_text(size = 15, face = "bold"),
#     #legend.position = "none", #去掉图例
#     #plot.title = element_text(hjust = 0.5,size=20) #改变标题位置和字体大小
#     plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5, lineheight = 0.5)
#   ) +
#   guides(color = guide_legend(override.aes = list(size = 8), ncol = 1)) +
#   labs(title = "Cell Type") + # 标题
#   theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) #边框粗细
# ggsave(p5, file = "sc_subset_umap_annot.pdf", width = 10, height = 6)













# 细胞比例柱状图
library(reshape2)
library(patchwork)
library(dplyr)
library(ggplot2)

plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
          "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", "#de4307", 
          "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", "#fa6e57", "#f69e53", 
          "#248888", "#999999", "#E7475E", "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", 
          "#FF6D31", "#f1ac9d", "#f06966", "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", 
          "#FF7A5A", "#8ED2C9", "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", 
          "#ffc300", "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
          "#f17d80", "#737495")

head(all_data_sc@meta.data)

out_gene_umi = as.matrix(all_data_sc@meta.data)
row_name = c("Cell", rownames(out_gene_umi))
out_gene_umi = rbind(colnames(out_gene_umi), out_gene_umi)
out_gene_umi = as.data.frame(cbind(as.matrix(row_name), out_gene_umi))
write.table(out_gene_umi, file = "scRNA/Filter_Basicinfo_nGene_nUMI_mito.txt", sep = "\t", quote = F, row.names = F, col.names = F)

data <- read.table("scRNA/Filter_Basicinfo_nGene_nUMI_mito.txt", header = T, sep = "\t", check.names = F)
data <- data[, c(1, 2, 3, 13, 18)]
colnames(data) <- c("Cells", "Sample", "Group", "Cluster", "Annotation")
head(data)
final_data <- data

mydata<-final_data %>% count(Sample, Annotation)
mymtrix<-dcast(mydata, Sample~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = "scRNA/Cell_type_stat.txt", quote = F, row.names = T, col.names = F, sep = "\t")

rownames(mymtrix)<-mymtrix[, "Sample"]
mymtrix<-mymtrix[, -1]
mymtrix<-as.matrix(mymtrix)
percent<-do.call(rbind, lapply(1:nrow(mymtrix), function(i){
  mymtrix[i, ]<-round(mymtrix[i, ]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent)<-rownames(mymtrix)
data4plot1<-melt(percent)
percent<-t(percent)
DATA<-cbind(rownames(percent), percent)
DATA<-rbind(c("Sample", colnames(percent)), DATA)
write.table(DATA, file = "scRNA/Cell_type_stat_percent.txt", quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot1)<-c("Sample", "Celltype", "Abund")

mydata<-final_data %>% count(Group, Annotation)
mymtrix<-dcast(mydata, Group~Annotation)
mymtrix[is.na(mymtrix)]<-0
mymtrix_out<-t(mymtrix)
write.table(as.data.frame(mymtrix_out), file = "scRNA/Cell_type_stat_group.txt", quote = F, row.names = T, col.names = F, sep = "\t")
rownames(mymtrix)<-mymtrix[, "Group"]
mymtrix<-mymtrix[, -1]
mymtrix<-as.matrix(mymtrix)
percent<-do.call(rbind, lapply(1:nrow(mymtrix), function(i){
  mymtrix[i, ]<-round(mymtrix[i, ]/rowSums(mymtrix)[i]*100, 2)
}))
rownames(percent)<-rownames(mymtrix)
data4plot2<-melt(percent)
percent<-t(percent)
DATA<-cbind(rownames(percent), percent)
DATA<-rbind(c("Group", colnames(percent)), DATA)
write.table(DATA, file = "scRNA/Cell_type_stat_percent_group.txt", quote = F, row.names = F, col.names = F, sep = "\t")
colnames(data4plot2)<-c("Group", "Celltype", "Abund")


# 限制顺序
# data4plot1$Sample <- factor(data4plot1$Sample, levels = c("M4-1", "M4-2", "M24-1", "M24-2"))
# data4plot2$Group <- factor(data4plot2$Group, levels = c("M4", "M24"))


# 柱状图
plot_list <- c()
plot_list[[1]] <- ggplot(data4plot1, aes(x = Sample, y = Abund, fill = Celltype)) +
  geom_bar(stat = "identity", size = 0.5, width = 0.7)+
  labs(x = "Sample", y = "Percent(%)", title = "")+
  scale_fill_manual(values = plan)+
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 15, colour = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 15, colour = "black"),
    axis.title = element_text(size = 15, colour = "black", face = "bold"),
    legend.title = element_text(size = 15, colour = "black", face = "bold"),
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5, lineheight = 0.5),
    strip.background = element_rect(fill = "white", size = 2))+NoLegend()

plot_list[[2]] <- ggplot(data4plot2, aes(x = Group, y = Abund, fill = Celltype)) +
  geom_bar(stat = "identity", size = 0.5, width = 0.7)+
  labs(x = "Group", y = "", title = "")+
  scale_fill_manual(values = plan)+
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 15, colour = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 15, colour = "black"),
    legend.text = element_text(size = 15, colour = "black"),
    legend.title = element_text(size = 15, colour = "black", face = "bold"),
    axis.title = element_text(size = 15, colour = "black", face = "bold"),
    plot.title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5, lineheight = 0.5),
    strip.background = element_rect(fill = "white", size = 2))

p <- wrap_plots(plot_list, ncol = 2)+
  plot_layout(widths = c(2, 1)) # 两张图各自的宽度
ggsave(p, file = "scRNA/ratio_barplot.pdf", width = 12, height = 6)













