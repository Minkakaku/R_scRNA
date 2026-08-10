setwd("/mnt/NexData_file/home/xzx/project/AE2509195142/shouhou4/")

library(monocle)
library(Seurat)
library(ggplot2)
#library(argparse)
library(dplyr)
library(tidyverse)
library(ggsci)
library(patchwork)
library(stringr)
library(getopt)

#读取参数和文件
command = matrix(c('rds','r',1,"character",
                   'type','t',2,"character",
                   'downcluster','D',0,"logical",
                   'maxcells','m',2,"numeric",
                   'cluster','c',2,"character",
                   'ptsize','p',2,"numeric",
                   'BEAM','B',0,"logical",
                   'output','o',1,"character")
                 ,byrow=TRUE, ncol=4)

args=getopt(command)


#### ===分析1 基因集评分=== ####
pbmc <- readRDS('./stromal_anno_v3.rds')
DimPlot(pbmc,group.by = 'leiden')

library(Seurat)
library(ggplot2)
color = read.delim("../shouhou3/stromal/leiden_color.xls")
celltype_color <- color$leiden.1
names(celltype_color) <- color$leiden
celltype_color
genesets_addmodulescore <- function(gene_sets,gene_name,scRNA,group_name,out_dir){
  genes_list <- list(gene_sets)
  
  scRNA_fin <- scRNA
  scRNA_fin <- AddModuleScore(
    object = scRNA_fin,
    features = genes_list,
    ctrl = 200,
    #nbin = 5,
    name = gene_name
  )
  # 提取需要的数据（细胞ID、分组、模块分数）
  plot_data <- data.frame(
    cell_type = sc@meta.data[[group_name]],  # 分组信息
    score = scRNA_fin@meta.data[[paste0(gene_name,"1")]]  # 模块分数
  )
  
  #plot_data$cell_type <- factor(plot_data$cell_type, levels = cell_type_order)
  p2 <- ggplot(plot_data, aes(x = cell_type, y = score, fill = cell_type)) +
    geom_violin(scale = "width") +
    scale_fill_manual(
      values = celltype_color
    )+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.background = element_rect(fill = 'white'),
          panel.grid = element_blank()) +
    labs(title = paste0(gene_name," Module Score"), x = '', y = "Score",fill=NULL)
  #+
  #  scale_fill_manual(values = cell_type_colors)
  p2
  ggsave(plot = p2,filename = paste0(out_dir,gene_name,"_",group_name,"_addmodulescore.png"),
         width = 9,height = 6,dpi = 300,bg = "white")
  ggsave(plot = p2,filename = paste0(out_dir,gene_name,"_",group_name,"_addmodulescore.pdf"),
         width = 9,height = 6,dpi = 300,bg = "white")
  
  return(scRNA_fin)
}

sc <- pbmc

Pericyte <- c('PDGFRB', 'RGS5', 'CSPG4','MCAM','ABCC9','KCNJ8','NOTCH3','DES','ANGPT1')
sc <- genesets_addmodulescore(gene_sets = Pericyte,gene_name = "Pericyte",scRNA = sc,group_name = "leiden",
                              out_dir = "./modulescore/")
# sc <- genesets_addmodulescore(gene_sets = Pericyte,gene_name = "Pericyte",scRNA = sc,group_name = "cell_type",
#                               out_dir = "./stromal/")

SMC <- c('ACTA2','TAGLN','MYH11','CNN1','TPM2','LMOD1','MYL9','CALD1','ACTG2','SMTN')
sc <- genesets_addmodulescore(gene_sets = SMC,gene_name = "SMC",scRNA = sc,group_name = "leiden",
                              out_dir = "./modulescore/")
# sc <- genesets_addmodulescore(gene_sets = SMC,gene_name = "SMC",scRNA = sc,group_name = "cell_type",
#                               out_dir = "./stromal/")


Drift <- c("VIMFN1","ITGA5","VCAN","POSTN","SPARC","CTGF","TGFBI","SERPINE1","TIMP1","MMP2","MMP14","LGALS3","SPP1")
sc <- genesets_addmodulescore(gene_sets = Drift,gene_name = "Drift",scRNA = sc,group_name = "leiden",
                              out_dir = "./modulescore/")

SASP <- c("IL6", "CXCL8", "CCL2", "CXCL1", "CXCL2", "CCL11", "CCL19", "GDF15", "B2M", "IL1B", "CSF2", "NFKBIA", "NFKB1", "RELA")
sc <- genesets_addmodulescore(gene_sets = SASP,gene_name = "SASP",scRNA = sc,group_name = "leiden",
                              out_dir = "./modulescore/")

ECM_drift <- c("FN1","VCAN","COL1A1","COL1A2","COL3A1","POSTN","SPARC","CTGF","TGFBI","THBS2","LRRC15","TNC","MMP2","TIMP1","SERPINE1")
sc <- genesets_addmodulescore(gene_sets = ECM_drift,gene_name = "ECM_drift",scRNA = sc,group_name = "leiden",
                              out_dir = "./modulescore/")

ECM_homeostatic <- c("DCN","LUM","COL14A1","MFAP4","FBLN1","COL6A1","COL6A2","COL6A3","PI16","CFD","DPT","C7")
sc <- genesets_addmodulescore(gene_sets = ECM_homeostatic,gene_name = "ECM_homeostatic",scRNA = sc,group_name = "leiden",
                              out_dir = "./modulescore/")
saveRDS(sc,"./modulescore/modulesc.rds")

################################################################################
#### === 分析2 monocle === ####
Idents(sc) <- sc$leiden
pbmc <- subset(sc,idents=c('3','6','8'))
DimPlot(pbmc, group.by = c('leiden','cell_type'))
Idents(pbmc) <- pbmc$cell_type
#table(pbmc$leiden)
#table(pbmc$cell_type_new)
table(pbmc$cell_type)

#pbmc <- subset(pbmc,ident=c('Epithelial_0', 'Epithelial_1', 'Epithelial_2', 'Epithelial_3', 'Epithelial_5'))
clustercolname <- 'cell_type'

sample_cells_per_cluster <- function(seurat_object,colnames="leiden", num_cells = 3000) {
  # 获取cluster信息
  Idents(seurat_object)=colnames
  clusters <- Idents(seurat_object)
  cluster_ids <- levels(clusters)  # 获取所有cluster的名称
  
  sampled_cell_ids <- c() # 用于存储所有抽样的细胞ID
  
  for (cluster_id in cluster_ids) {
    # 获取当前cluster的所有细胞ID
    cell_ids_in_cluster <- names(clusters[clusters == cluster_id])
    
    # 确定抽样数量
    num_cells_to_sample <- min(num_cells, length(cell_ids_in_cluster))  # 如果cluster细胞数小于1000，则抽取所有
    
    # 随机抽样
    sampled_ids <- sample(cell_ids_in_cluster, size = num_cells_to_sample, replace = FALSE)
    sampled_cell_ids <- c(sampled_cell_ids, sampled_ids) # 将抽样的细胞ID添加到列表中
  }
  subseurat = subset(seurat_object,cells=sampled_cell_ids)
  return(subseurat)
}

downcluster = TRUE
if(downcluster){
  if(!is.null(args$maxcells)){maxcell = args$maxcells}else{maxcell = 3000}
  pbmc = sample_cells_per_cluster(pbmc,num_cells = maxcell,colnames = clustercolname)
  print("Cell number of each cluster:")
  print(table(pbmc@meta.data[,clustercolname]))
}

p_data <- pbmc@meta.data
print("Meta data:")
head(p_data)

f_data <- data.frame(gene_short_name=row.names(pbmc),row.names=row.names(pbmc))
pd <- new('AnnotatedDataFrame',data = p_data)
fd <- new('AnnotatedDataFrame',data = f_data)
cds <- newCellDataSet(pbmc@assays$RNA@counts,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())#count适用于negbinomial.size,data适用于tobit且后续无需运行estimateDispersions
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds,min_expr=0.1)
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed>10))


diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr ="~cell_type" )
deg <- subset(diff,qval<0.01)
deg <- deg %>% arrange(qval,desc(num_cells_expressed))
deg <- deg[order(deg$qval,decreasing = F),]
head(deg)
deg_output <- data.frame(Genes = rownames(deg),deg[,c("pval","qval","num_cells_expressed")])
write.table(deg_output,file="./monocle/monocle.DEG.xls",col.names=T,row.names=F,sep="\t",quote=F)
ordergene <- rownames(deg)[1:min(1000,nrow(deg))]
length(ordergene)
cds <- setOrderingFilter(cds,ordergene[1:min(length(ordergene),2000)])
#cds <- setOrderingFilter(cds,ordergene)


pdf("./monocle/monocle.ordergenes.pdf")
plot_ordering_genes(cds)
dev.off()
png("./monocle/monocle.ordergenes.png")
plot_ordering_genes(cds)
dev.off()

cds <- reduceDimension(cds,max_components = 2,reduction_method='DDRTree')#,residualModelFormulaStr = "~sample"
cds <- orderCells(cds,reverse = TRUE)#重设起始点使用参数root_state=1
#cds <- orderCells(cds,root_state = 2)#重设起始点使用参数root_state=1


monocleplot <- function(cds,colnames,ptsize=1){
  
  p = plot_cell_trajectory(cds,color_by=colnames,size=ptsize,show_backbone=TRUE)
  ggsave(paste0("./monocle/monocle.trajectory.",colnames,".png"),p,width = 7,height = 7.5)
  ggsave(paste0("./monocle/monocle.trajectory.",colnames,".pdf"),p,width = 7,height = 7.5)
  
  ncols = ceiling(sqrt(length(unique(cds@phenoData@data[,colnames]))))
  nrows = ceiling(length(unique(cds@phenoData@data[,colnames]))/ncols)
  p = plot_cell_trajectory(cds,color_by=colnames,size=ptsize,show_backbone=TRUE)+facet_wrap(vars(!!sym(colnames)),nrow=nrows)
  ggsave(paste0("./monocle/monocle.trajectory.",colnames,".faceted.png"),p,width = 4*ncols,height = 4.2*nrows)
  ggsave(paste0("./monocle/monocle.trajectory.",colnames,".faceted.pdf"),p,width = 4*ncols,height = 4.2*nrows)
  
}

plotcol = c("group","State","cell_type","Pseudotime")
#if(!is.null(args$type)){plotcol = c(plotcol,typeuse)}
#if("cell_type" %in% colnames(pbmc@meta.data)){plotcol = c(plotcol,"cell_type")}
for(i in plotcol){
  if(i %in% colnames(cds@phenoData@data)){
    monocleplot(cds,colnames = i,ptsize=1)
  }else{
    print(paste0(i," is not exist in colnames"))
  }
  
}

#clustercolname <- 'group'


p = plot_complex_cell_trajectory(cds,x=1,y=2,color_by = clustercolname)
p
ggsave("./monocle/monocle.trajectory.tree.pdf",p,width = 7,height = 7.5)
ggsave("./monocle/monocle.trajectory.tree.png",p,width = 7,height = 7.5)

ncols = ceiling(sqrt(length(unique(cds@phenoData@data[,clustercolname]))))
nrows = ceiling(length(unique(cds@phenoData@data[,clustercolname]))/ncols)
p = plot_complex_cell_trajectory(cds,x=1,y=2,color_by = clustercolname)+facet_wrap(vars(!!sym(clustercolname)),nrow=nrows)
p
ggsave("./monocle/monocle.trajectory.tree.faceted.pdf",p,width = 4*ncols,height = 4.2*nrows)
ggsave("./monocle/monocle.trajectory.tree.faceted.png",p,width = 4*ncols,height = 4.2*nrows)

'Drift1','ECM_drift1','ECM_homeostatic1','SASP1','SMC1','Pericyte1'

cds@phenoData@data$ECM_drift <- cds@phenoData@data$ECM_drift1
cds@phenoData@data$Drift <- cds@phenoData@data$Drift1
cds@phenoData@data$ECM_homeostatic <- cds@phenoData@data$ECM_homeostatic1
cds@phenoData@data$SASP <- cds@phenoData@data$SASP1
cds@phenoData@data$SMC <- cds@phenoData@data$SMC1
cds@phenoData@data$Pericyte <- cds@phenoData@data$Pericyte1

#'Drift','ECM_drift','ECM_homeostatic',
col_name <- c('SASP','SMC','Pericyte','Drift','ECM_drift','ECM_homeostatic')

for (col in col_name) {
  p = plot_cell_trajectory(cds,x=1,y=2,color_by = col,cell_size = .5)
  ggsave(paste0("./monocle/monocle_pseudotime_celltype_",col,".png"),p,width = 7,height = 7.5)
  ggsave(paste0("./monocle/monocle_pseudotime_celltype_",col,".pdf"),p,width = 7,height = 7.5)
  
  p = plot_cell_trajectory(cds,x=1,y=2,color_by = col,cell_size = .5)+ 
    facet_wrap(~cell_type)+scale_color_gradientn(colours = c('grey','#ff7f0e','firebrick'), values = c(0,0.5,1,1.5))
  ggsave(paste0("./monocle/monocle_pseudotime_celltype_",col,"_faceted.png"),p,width = 9,height = 5)
  ggsave(paste0("./monocle/monocle_pseudotime_celltype_",col,"_faceted.pdf"),p,width = 9,height = 5)
}


pseudotime_df = cds@phenoData@data[,c("sample",clustercolname,"Size_Factor","num_genes_expressed","Pseudotime","State",'group')]
df = data.frame(Cells=rownames(pseudotime_df),pseudotime_df)
write.table(df,"./monocle/monocle.metadata.xls",sep = "\t",quote = F,row.names = F)


new_varMetadata <- data.frame(
  labelDescription = colnames(pData(cds)),
  row.names = colnames(pData(cds))
)

# 更新cds对象的varMetadata
varMetadata(cds) <- new_varMetadata
outpath = '/mnt/NexData_file/home/xzx/project/AE2509195142/shouhou4/monocle/'
for(branch_point in 1:((length(levels(cds@phenoData@data[,"State"]))-1)/2)){
  outbeam = paste0(outpath,"/BEAM_analysis.Node",branch_point)
  print(paste0("Analysis Beam Node ",branch_point))
  if(!dir.exists(outbeam)){
    dir.create(outbeam)
  }
  
  BEAM_res = BEAM(cds,branch_point = branch_point,cores = 4, progenitor_method = "duplicate")
  BEAM_res = BEAM_res[,c("gene_short_name","pval","qval")]
  BEAM_res = BEAM_res[order(BEAM_res$qval),]
  write.csv(BEAM_res,paste0(outbeam,"/BEAM_res.csv"),row.names = F)
  selectgenes = row.names(subset(BEAM_res,qval<0.05))
  p = plot_genes_branched_heatmap(cds[selectgenes[1:min(100,length(selectgenes))],],
                                  branch_point = branch_point,
                                  num_clusters = 4,
                                  cores=1,
                                  use_gene_short_name = T,
                                  show_rownames = T,return_heatmap = T)
  ggsave(paste0(outbeam,"/BEAM_heatmap.png"),p$ph_res,width=7,height = 8)
  ggsave(paste0(outbeam,"/BEAM_heatmap.pdf"),p$ph_res,width=7,height = 8)
  
  p=plot_genes_branched_pseudotime(cds[selectgenes[1:min(20,length(selectgenes))],],branch_point = branch_point,color_by = "State",ncol = 4)
  ggsave(paste0(outbeam,"/BEAM_pseudotime_top20.png"),p,width=15,height = 12)
  ggsave(paste0(outbeam,"/BEAM_pseudotime_top20.pdf"),p,width=15,height = 12)
  
  p=plot_genes_branched_pseudotime(cds[selectgenes[1:min(20,length(selectgenes))],],branch_point = branch_point,color_by = "cell_type",ncol = 4)
  ggsave(paste0(outbeam,"/BEAM_pseudotime_cell_type_top20.png"),p,width=15,height = 12)
  ggsave(paste0(outbeam,"/BEAM_pseudotime_cell_type_top20.pdf"),p,width=15,height = 12)
  
  p=plot_genes_branched_pseudotime(cds[selectgenes[1:min(20,length(selectgenes))],],branch_point = branch_point,color_by = "group",ncol = 4)
  ggsave(paste0(outbeam,"/BEAM_pseudotime_group_top20.png"),p,width=15,height = 12)
  ggsave(paste0(outbeam,"/BEAM_pseudotime_group_top20.pdf"),p,width=15,height = 12)
  
  
}



# 按拟时间做差异分析（检验基因表达是否随拟时间变化）
diff_test_res <- differentialGeneTest(cds[ordergene,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- diff_test_res %>% dplyr::filter(num_cells_expressed >10)%>% #防止小类干扰
  dplyr::arrange(qval) %>% head(100) %>% row.names() 
write.csv(file = './monocle/Pseudotime_ordergene.csv',diff_test_res)
p = plot_pseudotime_heatmap(cds[sig_gene_names,],
                            num_clusters = 5,
                            cores = 1,
                            show_rownames = T,return_heatmap = TRUE)
p
ggsave(filename = "./monocle/monocle.top100gene_cortime.png",p,h=8,w=6,dpi = 150)
ggsave(filename = "./monocle/monocle.top100gene_cortime.pdf",p,h=8,w=6,dpi = 150)

save.image("./monocle.Rdata")

# 将 Pseudotime 信息添加到 Seurat 对象的 meta.data 中
pbmc@meta.data$Pseudotime <- cds@phenoData@data$Pseudotime


# 提取 UMAP 和 Pseudotime 信息
umap_data <- FetchData(pbmc, vars = c("UMAP_1", "UMAP_2", "Pseudotime"))

# 绘制 Pseudotime 映射到 UMAP 图
p <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = Pseudotime)) +
  geom_point(size = 1) +
  #scale_color_gradientn(colors = c('blue', 'cyan', 'green', 'yellow', 'orange', 'red')) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))
ggsave(filename = "./monocle/monocle.pseudotime_UMAP.png",p,h=8,w=8,dpi = 300)
ggsave(filename = "./monocle/monocle.pseudotime_UMAP.pdf",p,h=8,w=8,dpi = 300)


################################################################################
####SDaiNC图表
#### Pericyte score与SMC score在cluster层面的统计:
df <- sc@meta.data

dt <- df %>%
  group_by(leiden) %>%
  summarise(
    Pericyte_mean = mean(Pericyte1),
    Pericyte_sd   = sd(Pericyte1),
    Pericyte_med  = median(Pericyte1),
    Pericyte_IQR  = IQR(Pericyte1),
    
    SMC_mean      = mean(SMC1),
    SMC_sd        = sd(SMC1),
    SMC_med       = median(SMC1),
    SMC_IQR       = IQR(SMC1),
    
    Drift_mean = mean(Drift1),
    Drift_sd   = sd(Drift1),
    Drift_med  = median(Drift1),
    Drift_IQR  = IQR(Drift1),
    
    SASP_mean      = mean(SASP1),
    SASP_sd        = sd(SASP1),
    SASP_med       = median(SASP1),
    SASP_IQR       = IQR(SASP1),
    
    ECM_drift_mean = mean(ECM_drift1),
    ECM_drift_sd   = sd(ECM_drift1),
    ECM_drift_med  = median(ECM_drift1),
    ECM_drift_IQR  = IQR(ECM_drift1),
    
    ECM_homeostatic_mean = mean(ECM_homeostatic1),
    ECM_homeostatic_sd   = sd(ECM_homeostatic1),
    ECM_homeostatic_med  = median(ECM_homeostatic1),
    ECM_homeostatic_IQR  = IQR(ECM_homeostatic1),
  )

write.csv(dt,"./summary_ScoreinLeiden.csv",row.names = FALSE)
