
###########################################################

#使用harmony整合多样本，去除批次效应
#安装harmony
#library(devtools)
#install_github("immunogenomics/harmony")

###########################################################
library(Seurat)
library(harmony)
library(ggsci)
library(ggplot2)
library(tidydr)
library(cols4all)
script_path <- normalizePath(dirname(rstudioapi::getSourceEditorContext()$path))
setwd(script_path)

my_theme <- theme_void()+theme(plot.margin = margin(5.5, 15, 5.5, 5.5))+
  theme_dr(xlength = 0.2, ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, units = "inches"), ends = "last", type = "closed"))+
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 15, color = "black"))
#----------------------------------------------------------
data.filt <- readRDS("../4-数据整合查看质控/all.data.rds")

#进行数据标准化
data.filt <- NormalizeData(data.filt)

#筛选高变基因
data.filt <- FindVariableFeatures(data.filt, nfeatures = 2000)

#进行数据归一化
data.filt <- ScaleData(data.filt)

# 对周期基因进行主成分分析
s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes, match = rownames(data.filt))

g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes, match = rownames(data.filt))

PCA <- RunPCA(data.filt, features = c(s_genes, g2m_genes))
# 对主成分分析结果进行可视化
PCA_plot <- DimPlot(PCA, reduction = "pca", group.by = "Phase", shuffle = T,
                    cols = c("G1" = "#219EBC", "G2M" = "#FEB705", "S" = "#E73F74"),
                    raster = F, pt.size = 0.75) + my_theme
ggsave("PCA_plot_before_phase.png", plot = PCA_plot, width = 7, height = 6, dpi = 300)
ggsave("PCA_plot_before_phase.pdf", plot = PCA_plot, width = 7, height = 6, dpi = 300)
rm(PCA)
#----------------------------------------------------------
#重新进行数据归一化，消除细胞周期等因素的影响
data.filt <- ScaleData(data.filt, vars.to.regress = c("S.Score", "G2M.Score", 
                                                      "nFeature_RNA", "nCount_RNA", 
                                                      "percent_mito", "percent_ribo"))
PCA <- RunPCA(data.filt, features = c(s_genes, g2m_genes))
# 对主成分分析结果进行可视化
PCA_plot <- DimPlot(PCA, reduction = "pca", group.by = "Phase", shuffle = T,
                    cols = c("G1" = "#219EBC", "G2M" = "#FEB705", "S" = "#E73F74"),
                    raster = F, pt.size = 0.75) + my_theme
ggsave("PCA_plot_after_phase.png", plot = PCA_plot, width = 7, height = 6, dpi = 300) 
ggsave("PCA_plot_after_phase.pdf", plot = PCA_plot, width = 7, height = 6, dpi = 300) 
rm(PCA)
#----------------------------------------------------------
#PCA降维
data.filt <- RunPCA(data.filt, npcs = 80)#指定80个主成分


#量化拐点识别
pct <- data.filt[["pca"]]@stdev / sum(data.filt[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
pc.use <- min(which(cumu > 90 & pct < 5)[1], 
              sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1)

p <- ElbowPlot(data.filt, ndims = 80)$data %>% ggplot() + 
  geom_point(aes(x = dims, y = stdev)) + 
  geom_vline(xintercept = pc.use, color = "#E73F74") +
  theme_bw() + labs(title = "Elbow plot: quantitative approach") + 
  theme(axis.text.x = element_text(hjust = 0.5, size = 10, color = "black"),
        axis.text.y = element_text(hjust = 0.5, size = 10, color = "black"),
        axis.ticks = element_line(linewidth = 0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(size = 12, color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
ggsave(plot = p, filename = "ElbowPlot.png", width = 6, height = 6, dpi = 300)
ggsave(plot = p, filename = "ElbowPlot.pdf", width = 6, height = 6, dpi = 300)

#降维聚类
pc.num <- 1:pc.use
data.filt <- RunUMAP(data.filt, dims = pc.num)
#--------------------------------------------------
#c4a_gui()
col_sample <- c("191598A1" = "#49B1DD", "191598A2" = "#E98DAF", "191598A3" = "#E9CF91", "191598B1" = "#95C6C6",
                "191598B2" = "#FFF4F2", "191598B3" = "#2A9D8E", "191598C1" = "#FEB705", "191598C2" = "#F66F69",
                "191598C3" = "#E73F74","191598D1" = "#FFBFFD", "191598D2" = "#EFFDAF", "191598D3" = "#EF5FFF")

plot1 <- DimPlot(data.filt, group.by = "Samples", reduction = "umap", shuffle = T,
                 cols = col_sample, raster = F, pt.size = 0.5) + my_theme
ggsave("QC_umap_before_harmony.pdf", plot = plot1, width = 10, height = 8, dpi = 300)
ggsave("QC_umap_before_harmony.png", plot = plot1, width = 10, height = 8, dpi = 300)
#########################################################################################

##Harmony整合

#########################################################################################
system.time({scRNA_harmony <- RunHarmony(data.filt, group.by.vars = "orig.ident", block.size = 0.01,
                                         max.iter.harmony = 40, max.iter.cluster = 40)})

scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = pc.num)

plot3 <- DimPlot(scRNA_harmony, group.by='Samples', reduction = "umap", shuffle = T,
                 cols = col_sample, raster = F, pt.size = 0.5) + my_theme
ggsave("scRNA_umap_harmony.pdf", plot = plot3, width = 10, height = 8, dpi = 300)
ggsave("scRNA_umap_harmony.png", plot = plot3, width = 10, height = 8, dpi = 300)

#run garbage collect to free up memory
gc()
##################################################################
#保存去除批次效应后的数据
dim(scRNA_harmony)
#23291 83229
saveRDS(scRNA_harmony, "scRNA_harmony.rds")

