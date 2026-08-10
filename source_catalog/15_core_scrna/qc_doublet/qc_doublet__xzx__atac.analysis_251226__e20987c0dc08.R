.libPaths(c("/home/szt/miniforge3/envs/rpy/lib/R/library",
            "/home/tt/R/x86_64-conda-linux-gnu-library/4.3"
))


library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(reshape2)
library(ggplot2)
library(dplyr)
# load the RNA and ATAC data

input = "/mnt/NexData_file/nas/outputdata/cellranger_out/FWSC20251628"
outpath = "/mnt/NexData_file/home/xzx/project/FWSC20251628"
if(!dir.exists(outpath)){dir.create(outpath)}
setwd(outpath)

#allsample = list.dirs(input,recursive = F,full.names = F)
allsample = c('CHL722','LY0410','LYW25112','YT0702')
# gg = list()
# for(sample in allsample){
#   counts <- Read10X_h5(paste0(input,"/",sample,"/outs/filtered_feature_bc_matrix.h5"))
#   fragpath <- paste0(input,"/",sample,"/outs/atac_fragments.tsv.gz")
#   gg[[sample]] <- CreateSeuratObject(
#     counts = counts$`Gene Expression`,
#     assay = "RNA",project = sample
#   )
#   gg[[sample]][["percent.mt"]] <- PercentageFeatureSet(gg[[sample]], pattern = "^MT-")##线粒体
#   gg[[sample]][["percent.rb"]] <- PercentageFeatureSet(gg[[sample]], pattern = "^RP[SL]")##核糖体
  
#   # create ATAC assay and add it to the object
#   grange.counts <- StringToGRanges(rownames(counts$Peaks),sep=c(":","-"))
#   grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
#   atac_counts = counts$Peaks[as.vector(grange.use),]
  
#   annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
#   seqlevelsStyle(annotation) <- "UCSC"
#   genome(annotation)<- "hg38"
  
#   gg[[sample]][["ATAC"]] <- CreateChromatinAssay(
#     counts = atac_counts,
#     sep = c(":", "-"),
#     genome = "hg38",
#     fragments = fragpath,
#     annotation = annotation,
#     min.cells = 10
#   )
  
  
#   DefaultAssay(gg[[sample]]) <- "ATAC"
#   gg[[sample]] <- NucleosomeSignal(gg[[sample]])
#   gg[[sample]] <- TSSEnrichment(gg[[sample]],fast=FALSE)
  
  
#   ###2.BeforFilter
#   if(!dir.exists("02.Data_QC")){dir.create("02.Data_QC")}
#   dir.create(paste0("02.Data_QC/",sample))
#   p = DensityScatter(gg[[sample]], x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
#   ggsave(paste0("02.Data_QC/",sample,"/DensityScatter.pdf"),p,width = 8,height = 8)
#   ggsave(paste0("02.Data_QC/",sample,"/DensityScatter.png"),p,width = 8,height = 8)
  
#   p = VlnPlot(
#     object = gg[[sample]],
#     features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.rb", "nCount_ATAC","nFeature_ATAC", "TSS.enrichment", "nucleosome_signal"),
#     ncol = 8,
#     pt.size = 0
#   )
  
#   ggsave(paste0("02.Data_QC/",sample,"/ViolinPlot.beforQC.pdf"),p,width = 24,height = 10)
#   ggsave(paste0("02.Data_QC/",sample,"/ViolinPlot.beforQC.png"),p,width = 24,height = 10)
  
#   gg[[sample]]$high.tss <- ifelse(gg[[sample]]$TSS.enrichment > 4 , 'High', 'Low')
#   TSS<-TSSPlot(gg[[sample]], group.by = 'high.tss') + NoLegend()
#   ggsave(paste0("02.Data_QC/",sample,"/TSSplot.pdf"), plot = TSS, width = 10, height = 8)
#   ggsave(paste0("02.Data_QC/",sample,"/TSSplot.png"), plot = TSS, width = 10, height = 8)
  
#   gg[[sample]]$nucleosome_group <- ifelse(gg[[sample]]$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
#   Frag<-FragmentHistogram(object = gg[[sample]], group.by = 'nucleosome_group')
#   ggsave(paste0("02.Data_QC/",sample,"/Nucleosome banding pattern diagram.pdf"), plot = Frag, width = 10, height = 8)
#   ggsave(paste0("02.Data_QC/",sample,"/Nucleosome banding pattern diagram.png"), plot = Frag, width = 10, height = 8)
  
#   plot1 <- FeatureScatter(object = gg[[sample]], feature1 = "nCount_RNA", feature2 = "percent.mt")
#   plot2 <- FeatureScatter(object = gg[[sample]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
#   p = CombinePlots(plots = list(plot1, plot2),ncol=2)
#   ggsave(paste0("02.Data_QC/",sample,"/Counts_correlation.png"),p,width = 15,height = 7,bg="white")
#   ggsave(paste0("02.Data_QC/",sample,"/Counts_correlation.pdf"),p,width = 15,height = 7)
  
# }


# # create a Seurat object containing the RNA adata

# pbmc = merge(gg[[1]],gg[-1])


# ###2.BeforFilter
# if(!dir.exists("02.Data_QC")){dir.create("02.Data_QC")}

# p = VlnPlot(
#   object = pbmc,
#   features = c("nCount_ATAC","nFeature_ATAC", "TSS.enrichment", "nucleosome_signal"),
#   ncol = 4,
#   pt.size = 0
# )

# ggsave("02.Data_QC/ViolinPlot.ATAC.beforQC.pdf",p,width = 20,height = 10)
# ggsave("02.Data_QC/ViolinPlot.ATAC.beforQC.png",p,width = 20,height = 10)

# p = VlnPlot(
#   object = pbmc,
#   features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.rb"),
#   ncol = 4,
#   pt.size = 0
# )

# ggsave("02.Data_QC/ViolinPlot.scRNA.beforQC.pdf",p,width = 20,height = 10)
# ggsave("02.Data_QC/ViolinPlot.scRNA.beforQC.png",p,width = 20,height = 10)


# plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# p = CombinePlots(plots = list(plot1, plot2),ncol=2)
# ggsave("02.Data_QC/Counts_correlation.png",p,width = 15,height = 7,bg="white")
# ggsave("02.Data_QC/Counts_correlation.pdf",p,width = 15,height = 7)

# write.table(pbmc@meta.data,"02.Data_QC/CellInfos.Raw.xls",sep = "\t",quote = F)

# saveRDS(pbmc,"02.Data_QC/result.raw.rds")
pbmc = readRDS(paste0(outpath, "/02.Data_QC/result.raw.rds"))
pbmc <- subset(
  x = pbmc,
  subset = nCount_RNA < 15000 &
    nCount_RNA > 500 &
    nFeature_RNA < 5000 &
    percent.mt < 30 &
    percent.rb < 15 &
    nCount_ATAC < 40000 &
    nCount_ATAC > 600 &
    nucleosome_signal < 2 &
    TSS.enrichment > 4
)

###3.AfterFilter
if(!dir.exists("03.Cellfilter")){dir.create("03.Cellfilter")}
p = VlnPlot(
  object = pbmc,
  features = c("nCount_ATAC","nFeature_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

ggsave("03.Cellfilter/ViolinPlot.ATAC.afterQC.pdf",p,width = 20,height = 10)
ggsave("03.Cellfilter/ViolinPlot.ATAC.afterQC.png",p,width = 20,height = 10)

p = VlnPlot(
  object = pbmc,
  features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.rb"),
  ncol = 4,
  pt.size = 0
)

ggsave("03.Cellfilter/ViolinPlot.scRNA.afterQC.pdf",p,width = 20,height = 10)
ggsave("03.Cellfilter/ViolinPlot.scRNA.afterQC.png",p,width = 20,height = 10)


plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

p = CombinePlots(plots = list(plot1, plot2),ncol=2)
ggsave("03.Cellfilter/Counts_correlation.png",p,width = 15,height = 7,bg="white")
ggsave("03.Cellfilter/Counts_correlation.pdf",p,width = 15,height = 7)



for(sample in allsample){
  dir.create(paste0("03.Cellfilter/",sample))
  p = DensityScatter(pbmc, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
  ggsave(paste0("03.Cellfilter/",sample,"/DensityScatter.pdf"),p,width = 8,height = 8)
  ggsave(paste0("03.Cellfilter/",sample,"/DensityScatter.png"),p,width = 8,height = 8)
  
  # pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 4 , 'High', 'Low')
  # TSS<-TSSPlot(pbmc, group.by = 'high.tss')
  # ggsave(paste0("03.Cellfilter/",sample,"/TSSplot.pdf"), plot = TSS, width = 8, height = 8)
  # ggsave(paste0("03.Cellfilter/",sample,"/TSSplot.png"), plot = TSS, width = 8, height = 8)
  # 
  Frag<-FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')
  ggsave(paste0("03.Cellfilter/",sample,"/Nucleosome banding pattern diagram.pdf"), plot = Frag, width = 10, height = 8)
  ggsave(paste0("03.Cellfilter/",sample,"/Nucleosome banding pattern diagram.png"), plot = Frag, width = 10, height = 8)
}



write.table(data.frame(Cells=rownames(pbmc@meta.data),pbmc@meta.data),"03.Cellfilter/CellInfos.Filter.xls",sep="\t",quote=F,row.names=F)

saveRDS(pbmc,"03.Cellfilter/result.QC.rds")

###4.降维聚类
##单细胞降维
if(!dir.exists("04.DimensionReduction_Cluster")){dir.create("04.DimensionReduction_Cluster")}
DefaultAssay(pbmc) <- "RNA"
pbmc <- JoinLayers(object = pbmc,assay = "RNA")
pbmc <- SCTransform(pbmc)
pbmc <- RunPCA(pbmc)

DefaultAssay(pbmc) <- "SCT"
top10 <- head(VariableFeatures(pbmc), 10)

plot1 <- VariableFeaturePlot(pbmc)
p <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
ggsave("04.DimensionReduction_Cluster/VariableGenesPlot.pdf",p,width = 10,height = 8)
ggsave("04.DimensionReduction_Cluster/VariableGenesPlot.png",p,width = 10,height = 8,bg="white")

gg = VizDimLoadings(pbmc, dims = 1:20, reduction = "pca",balanced = TRUE,nfeatures = 20)
ggsave("04.DimensionReduction_Cluster/VizPCAPlot.pdf",gg,width = 12,height = 18)
ggsave("04.DimensionReduction_Cluster/VizPCAPlot.png",gg,width = 12,height = 18)

gg = DimHeatmap(pbmc, dims = 1:20, cells = min(500,length(colnames(pbmc))),nfeatures = 20, balanced = TRUE, fast = FALSE,ncol = 4)
ggsave("04.DimensionReduction_Cluster/PCHeatmap.pdf",gg,width = 12,height = 18)
ggsave("04.DimensionReduction_Cluster/PCHeatmap.png",gg,width = 12,height = 18)

gg = ElbowPlot(pbmc,ndims = 50)
ggsave("04.DimensionReduction_Cluster/PCElbowPlot.pdf",gg,width = 8,height = 8)
ggsave("04.DimensionReduction_Cluster/PCElbowPlot.png",gg,width = 8,height = 8,bg="white")

gg = DimPlot(pbmc, reduction = "pca")
ggsave("04.DimensionReduction_Cluster/PCAPlot.pdf",gg,width = 10,height = 8)
ggsave("04.DimensionReduction_Cluster/PCAPlot.png",gg,width = 10,height = 8)

pbmc <- RunUMAP(pbmc,reduction = "pca",dims=1:15, n.neighbors = 15,reduction.name = "rna.umap",reduction.key = 'rnaUMAP_')

##ATAC降维
DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc,reduction = "lsi",dims = 2:10,n.neighbors = 15,reduction.name = "atac.umap",reduction.key = "atacUMAP_")


pbmc <- FindMultiModalNeighbors(pbmc,reduction.list = list("pca","lsi") ,dims.list = list(1:15,2:10))
pbmc <- RunUMAP(pbmc,nn.name = "weighted.nn",reduction.name = "wnn.umap",reduction.key = "wnnUMAP_")
pbmc <- FindClusters(pbmc,graph.name = "wsnn", algorithm = 3, verbose = F,resolution = 0.5)

p = DimPlot(pbmc, reduction = "wnn.umap",group.by='seurat_clusters', label = TRUE, label.size = 4.5, repel = TRUE) + ggtitle("WNN")
ggsave("04.DimensionReduction_Cluster/cluster_umap.png",p,width = 10,height = 10)
saveRDS(pbmc,"04.DimensionReduction_Cluster/result.dim.rds")



###5.markergene
if(!dir.exists("05.MarkerGene")){dir.create("05.MarkerGene")}
DefaultAssay(pbmc)="SCT"

allmarkers <- FindAllMarkers(object = pbmc, only.pos = T, min.pct = 0.1,logfc.threshold = 0.25, test.use = "wilcox",base=2)
print(head(allmarkers))
write.table(data.frame(Genes = allmarkers$gene,allmarkers[,1:6]),file="05.MarkerGene/AllMarkerGenes.xls",sep="\t",row.names=F,quote=F)

top10 <- allmarkers %>% dplyr::filter(p_val_adj <0.05) %>% group_by(cluster) %>% top_n(10, avg_log2FC)

p = DoHeatmap(object = pbmc,top10$gene)
ggsave("05.MarkerGene/Markers.top10.Heatmap.png",p,width=20,height=20)
ggsave("05.MarkerGene/Markers.top10.Heatmap.pdf",p,width=20,height=20)

# top3 <- allmarkers %>% dplyr::filter(p_val_adj <0.05) %>% group_by(cluster) %>% top_n(3, avg_log2FC)
# top3_markers = top3$gene
# names(top3_markers) = top3$cluster
# p = DotPlot(pbmc, group.by = "seurat_clusters", features = top3_markers, dot.scale = 10,cols=c("lightgrey","red"))+theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))
# ggsave("05.MarkerGene/Markers.top3.DotPlot.png",p,width = 20,height = 10,bg="white")
# ggsave("05.MarkerGene/Markers.top3.DotPlot.pdf",p,width = 20,height = 10)

if(!dir.exists("05.MarkerGene/MarkerUMAP")){dir.create("05.MarkerGene/MarkerUMAP")}
allmarkers = dplyr::filter(allmarkers,p_val_adj<0.05)
for(j in levels(pbmc)){
  markers <- allmarkers[allmarkers$cluster==j,]
  plotmarkers = markers$gene[1:min(20,nrow(markers))]
  ncols = ceiling(sqrt(min(20,length(plotmarkers))))
  nrows = floor(length(plotmarkers)/ncols)
  gg = FeaturePlot(pbmc, plotmarkers, cols = c("lightgrey", "red"),reduction = "wnn.umap",pt.size =0.5,ncol = ncols)
  ggsave(paste0("05.MarkerGene/MarkerUMAP/MarkersUMAPPlot_Cluster",j,".png"),gg,width=1.25*ncols*4,height=nrows*4)
  ggsave(paste0("05.MarkerGene/MarkerUMAP/MarkersUMAPPlot_Cluster",j,".pdf"),gg,width=1.25*ncols*4,height=nrows*4)
}

###6.Celltype_annotation
if(!dir.exists("06.Celltype_annotation")){dir.create("06.Celltype_annotation")}

table = table(Idents(pbmc),pbmc@meta.data$orig.ident)
table = dcast(data.frame(table),Var1~Var2)
colnames(table)[1]="Cluster"
write.table(table,file="06.Celltype_annotation/CellNumber.Statistics.xls",sep="\t",row.names=F,quot=F)
data = data.frame(Cell = rownames(pbmc@meta.data),Cluster=pbmc@meta.data$seurat_clusters)
write.table(data,file="06.Celltype_annotation/CellIdent.Summary.xls",sep="\t",row.names=F,quot=F)


p = FeaturePlot(pbmc,reduction = "wnn.umap",c("PTPRC","CD3D","CD3E","CD79A","MS4A1","CD19","JCHAIN","LYZ","CSF1R","C1QA","CSF3R","HLA-DRA","PECAM1","CDH5","ACTA2","TAGLN","DCN","COL1A1","RGS5","KCNJ8","PAX7","MYF5","EPCAM","KRT6B","KRT6C","KRT14","KRT7","SLPI"),ncol = 5)
ggsave("06.Celltype_annotation/Feature.UMAPPlot.png",p,width = 25,height = 28)
ggsave("06.Celltype_annotation/Feature.UMAPPlot.pdf",p,width = 25,height = 28)

write.table(data.frame(Cells=rownames(pbmc@meta.data),pbmc@meta.data),"06.Celltype_annotation/CellInfos.final.xls",sep="\t",quote=F,row.names=F)

saveRDS(pbmc,"06.Celltype_annotation/result.rds")

# type = read.table("D:\\单细胞项目\\30.FWSC20251628\\result_multi\\06.Celltype_annotation\\type.txt",header = T,row.names = 1,sep = "\t")
# pbmc = AddMetaData(pbmc,type,col.name = "cell_type")
# rank = names(sort(-table(pbmc@meta.data$cell_type)))
# pbmc@meta.data$cell_type = factor(pbmc@meta.data$cell_type,levels = rank)
# p1 <- DimPlot(pbmc, reduction = "rna.umap",group.by='cell_type', label = TRUE, label.size = 4.5, repel = TRUE) + ggtitle("RNA")
# p2 <- DimPlot(pbmc, reduction = "atac.umap",group.by='cell_type', label = TRUE, label.size = 4.5, repel = TRUE) + ggtitle("ATAC")
# p3 <- DimPlot(pbmc, reduction = "wnn.umap",group.by='cell_type', label = TRUE, label.size = 4.5, repel = TRUE) + ggtitle("WNN")
# 
# p = (p1 & NoLegend()) + (p2 & NoLegend()) + p3 & theme(plot.title = element_text(hjust = 0.5))
# ggsave("06.Celltype_annotation/Celltype.UMAPPlot.png",p,width = 24,height = 8)
# ggsave("06.Celltype_annotation/Celltype.UMAPPlot.pdf",p,width = 24,height = 8)
# 
# p1 <- DimPlot(pbmc, reduction = "rna.umap",group.by='seurat_clusters', label = TRUE, label.size = 4.5, repel = TRUE) + ggtitle("RNA")
# p2 <- DimPlot(pbmc, reduction = "atac.umap",group.by='seurat_clusters', label = TRUE, label.size = 4.5, repel = TRUE) + ggtitle("ATAC")
# p3 <- DimPlot(pbmc, reduction = "wnn.umap",group.by='seurat_clusters', label = TRUE, label.size = 4.5, repel = TRUE) + ggtitle("WNN")
# 
# p = (p1 & NoLegend()) + (p2 & NoLegend()) + p3 & theme(plot.title = element_text(hjust = 0.5))
# ggsave("06.Celltype_annotation/Cluster.UMAPPlot.png",p,width = 24,height = 8)
# ggsave("06.Celltype_annotation/Cluster.UMAPPlot.pdf",p,width = 24,height = 8)
# 
# p1 <- DimPlot(pbmc, reduction = "rna.umap",group.by='orig.ident', label = TRUE, label.size = 4.5, repel = TRUE) + ggtitle("RNA")
# p2 <- DimPlot(pbmc, reduction = "atac.umap",group.by='orig.ident', label = TRUE, label.size = 4.5, repel = TRUE) + ggtitle("ATAC")
# p3 <- DimPlot(pbmc, reduction = "wnn.umap",group.by='orig.ident', label = TRUE, label.size = 4.5, repel = TRUE) + ggtitle("WNN")
# 
# p = (p1 & NoLegend()) + (p2 & NoLegend()) + p3 & theme(plot.title = element_text(hjust = 0.5))
# ggsave("06.Celltype_annotation/Sample.UMAPPlot.png",p,width = 24,height = 8)
# ggsave("06.Celltype_annotation/Sample.UMAPPlot.pdf",p,width = 24,height = 8)

# dat = as.matrix(table(pbmc@meta.data$cell_type))
# df1 = data.frame(cell_type = rownames(dat),number = dat[,1])
# df1$prop = df1$number/sum(df1$number)
# df1$result = paste0(round(df1$prop*100, 2),"% ",df1$cell_type)
# 
# df1 = df1[rank,]
# df1$result=factor(df1$result,levels = df1$result)
# library(scales)
# p = ggplot(df1,aes(x=1,y=prop,fill = result))+scale_x_continuous(limits=c(0,1.45))+
#   geom_col(width = 0.9)+theme_classic()+coord_polar("y",direction = -1)+labs(title="Cell Frequency",fill=sample)+scale_fill_manual(values = hue_pal()(nrow(df1)))+
#   theme(axis.line = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),plot.title  = element_text(hjust=0.5))
# ggsave("06.Celltype_annotation/PropPlot_frequency.png",p,width = 5.5,height = 5)
# ggsave("06.Celltype_annotation/PropPlot_frequency.pdf",p,width = 5.5,height = 5)
