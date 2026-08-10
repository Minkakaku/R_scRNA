script_path <- normalizePath(dirname(rstudioapi::getSourceEditorContext()$path))
setwd(script_path)
library(Seurat)
library(openxlsx)
library(SeuratWrappers)
library(ggsci)
library(ggplot2)
library(tidydr)

options(future.globals.maxSize = 1e9)

my_theme <- theme_void()+theme(plot.margin = margin(5.5, 15, 5.5, 5.5))+
  theme_dr(xlength = 0.2, ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, units = "inches"), ends = "last", type = "closed"))+
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 15, color = "black"))
#--------------------------------------------------------
#合并所有样本
dir <- "../3-decontX去除游离RNA污染/scData/"
sample_name <- gsub(pattern = "_decontX.rds", replacement = "", x = dir(dir))

for (i in 1:length(sample_name)) {
  assign(x = sample_name[i], 
         value = readRDS(paste0(dir, sample_name[i], "_decontX.rds")))
}

#add.cell.ids参数的作用是在细胞ID前添加分组信息，但是要注意，信息的顺序需要和合并的顺序保持一致
all.data <- merge(x = sdata.191598A1, 
                  y = c(sdata.191598A2, sdata.191598A3,
                        sdata.191598B1, sdata.191598B2, 
                        sdata.191598B3, sdata.191598C1, 
                        sdata.191598C2, sdata.191598C3,
                        sdata.191598D1, sdata.191598D2,
                        sdata.191598D3), 
                  add.cell.ids = c("191598A1", "191598A2", "191598A3", "191598B1", "191598B2", 
                                   "191598B3", "191598C1", "191598C2", "191598C3", "191598D1", 
                                   "191598D2", "191598D3"))
dim(all.data)
# 23291 83229
#--------------------------------------------------------
#数据质量控制指标的可视化
all.data$Samples <- factor(x = all.data$orig.ident,
                           levels = c("191598A1", "191598A2", "191598A3", "191598B1", "191598B2", 
                                     "191598B3", "191598C1", "191598C2", "191598C3", "191598D1", 
                                      "191598D2", "191598D3"))

#添加其它的metadata
col_sample <- c("191598A1" = "#49B1DD", "191598A2" = "#E98DAF", "191598A3" = "#E9CF91", "191598B1" = "#95C6C6",
                "191598B2" = "#FFF4F2", "191598B3" = "#2A9D8E", "191598C1" = "#FEB705", "191598C2" = "#F66F69",
                "191598C3" = "#E73F74","191598D1" = "#FFBFFD", "191598D2" = "#EFFDAF", "191598D3" = "#EF5FFF")
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo")

p <- VlnPlot(all.data, features = feats, group.by = "orig.ident", 
             ncol = 2, pt.size = 0, raster = F, cols = col_sample) + NoLegend()
ggsave(plot = p, filename = "Samples.png", width = 12, height = 8, dpi = 300)
ggsave(plot = p, filename = "Samples.pdf", width = 12, height = 8, dpi = 300)
#--------------------------------------------------------
#计算细胞周期
#提取g2m特征向量
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes, match = rownames(all.data))

#提取s期特征向量
s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes, match = rownames(all.data))

#对细胞周期阶段进行评分
all.data <- CellCycleScoring(object = all.data, g2m.features = g2m_genes, s.features = s_genes)

saveRDS(object = all.data, file = "all.data.rds")
