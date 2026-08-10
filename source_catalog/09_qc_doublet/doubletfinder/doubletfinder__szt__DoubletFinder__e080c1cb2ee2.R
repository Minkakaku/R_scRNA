
#进一步去除双胞
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
###########################################################
library(Seurat)
packageVersion("Seurat")
library(DoubletFinder)
library(ggplot2)
library(tidydr)
script_path <- normalizePath(dirname(rstudioapi::getSourceEditorContext()$path))
setwd(script_path)

my_theme <- theme_void()+theme(plot.margin = margin(5.5, 15, 5.5, 5.5))+
  theme_dr(xlength = 0.2, ylength = 0.2, arrow = grid::arrow(length = unit(0.1, units = "inches"), ends = "last", type = "closed"))+
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 20, color = "black"),
        axis.title = element_text(size = 15, color = "black"))
###########################################################
dir <- "../1-数据读取和预处理/scData/"
sample_name <- gsub(pattern = ".rds", replacement = "", x = dir(dir))

for (i in 1:length(sample_name)) {
  #i <- 1
  data.filt <- readRDS(paste0(dir, sample_name[i], ".rds"))

  #进行数据归一化
  data.filt <- NormalizeData(data.filt)
  data.filt <- FindVariableFeatures(data.filt, selection.method = "vst", nfeatures = 2000)
  data.filt <- ScaleData(data.filt)
  
  #PCA降维
  data.filt <- RunPCA(data.filt, npcs = 80)
  
  pc.num <- 1:40
  data.filt <- RunUMAP(data.filt, dims = pc.num)
  data.filt <- FindNeighbors(data.filt, dims = pc.num)
  data.filt <- FindClusters(data.filt, resolution = 1)
  
  #寻找最优pK值
  sweep.res <- paramSweep(data.filt, PCs = pc.num, sct = F)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
  
  #define the expected number of doublet cells.
  homotypic.prop <- modelHomotypic(data.filt$seurat_clusters) #最好提供celltype
  DoubletRate <- ncol(data.filt)*8*1e-6
  nExp_poi <- round(DoubletRate*ncol(data.filt))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  #使用doubletFinder_v3函数预测双细胞
  data.filt <- doubletFinder(data.filt, PCs = pc.num, pN = 0.25, pK = pK_bcmvn,
                                nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
  DF.name <- colnames(data.filt@meta.data)[grepl("DF.classification", colnames(data.filt@meta.data))]
  print(table(data.filt@meta.data[,DF.name]))
  
  #将双细胞检测结果可视化
  p <- DimPlot(data.filt, reduction = "umap", group.by = DF.name, 
               cols = c("Singlet" = "#219EBCCC", "Doublet" = "#E73F74CC"), 
               shuffle = T, order = c("Doublet", "Singlet"),
               raster = F, pt.size = 0.75) + my_theme + NoLegend() + 
    labs(title = paste("DoubletFinder", sample_name[i]))
  ggsave(plot = p, filename = paste0(sample_name[i], "_Singlet_umap.png"), width = 6.5, height = 6.5, dpi = 300)
  
  #移除这些双细胞
  print(dim(data.filt))
  data.filt <- data.filt[,data.filt@meta.data[,DF.name] == "Singlet"]
  print(dim(data.filt))
  ##################################################################
  #保存质控过滤后的数据
  epi_mat <- data.filt@assays$RNA@counts #提取处理后的基因表达矩阵。
  epi_phe <- data.filt@meta.data[,c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo")]#提取处理后的元数据。
  data.filt <- CreateSeuratObject(counts = epi_mat, meta.data = epi_phe)
  
  saveRDS(data.filt, paste0(sample_name[i], "_Singlet.rds"))
}

