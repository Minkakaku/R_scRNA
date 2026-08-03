all_data2 <- readRDS("result0202/all_data.rds")

result_dir <- "result0325"
dir.create(result_dir)

all_data2@meta.data$group %>% table()
#OE、EAOC分析------
all_data_OE_EAOC <- subset(all_data2, group %in% c("EAOC","OE"))
all_data_OE_EAOC@meta.data$group %>% table()

#图1d、e、总的ump  分组umap  注释marker热图------
all_data_OE_EAOC@meta.data$subtype %>% table()
all_data_OE_EAOC@meta.data$celltype %>% table()

all_data_OE_EAOC@meta.data$celltype %>% unique()
all_data_OE_EAOC@meta.data$celltype <- as.character(all_data_OE_EAOC@meta.data$celltype)

celltype_level <- c("Endothelial", "Pericyte", "Stromal", "SMC", "Myeloid", "Mast", "T_NK" ,"Epithelial")
all_data_OE_EAOC@meta.data$celltype <- factor(all_data_OE_EAOC@meta.data$celltype,levels = celltype_level)

saveRDS(all_data_OE_EAOC,"result0325/all_data_OE_EAOC.rds")

all_data_OE_EAOC <- readRDS("result0325/all_data_OE_EAOC.rds")

plan <- c("#EA7070","#FDC4B6","#E59572","#2694AB","#96CEB4","#ED5485",
          "#D9534F","#FFAD60","#57D1C9","#FFEEAD","#7F95D1","#FFE869",
          "#DE4307","#F29C2B","#F6D04D","#8BC24C","#4695D6","#FED95C",
          "#FA6E57","#F69E53","#248888","#999999","#E7475E","#F0D879",
          "#29A2C6","#FFCB18","#73B66B","#FF6D31","#F1AC9D","#F06966",
          "#DEE2D1","#6ABE83","#FCF4D9","#FFB85F","#FF7A5A","#8ED2C9",
          "#A5DFF9","#EF5285","#60C5BA","#FEEE7D","#E8A0B8","#FFC300",
          "#BCCF3D","#02C9C9","#FBF579","#005995","#FA625F","#600473",
          "#F17D80","#737495")

plan <- c("#005995",
          "#96CEB4",
          "#248888",
          "#73B66B",
          "#FF6D31",
          #"#FDC4B6",
          "#EF5285",
          "#FFC300",
          "#900475"
)

pdf(file.path(result_dir, str_glue("Dimplot celltype umap.pdf")), height = 6, width = 8)
DimPlot(all_data_OE_EAOC, 
        group.by = "celltype", 
        label = TRUE, 
        cols = plan,
        reduction = "harmonyumap")

dev.off()

pdf(file.path(result_dir, str_glue("Dimplot celltype umap group.pdf")), height = 6, width = 12)
DimPlot(all_data_OE_EAOC, 
        group.by = "celltype", 
        label = TRUE, 
        cols = plan,
        reduction = "harmonyumap",
        split.by = "group"
)
dev.off()


# major_marker_ls <- c(
#   "Endothelial" = "CDH5",
#   "Endothelial" = "PECAM1",
#   "Endothelial" = "CD34",
# 
#   "Pericyte" = "RGS5",
#   "Pericyte" = "PDGFRB",
#   "Pericyte" = "MCAM",
#   
#   "Stromal" = "OSR2",
#   "Stromal" = "PDGFRA",
#   
#   "SMC" = "MYLK",
#   "SMC" = "CSRP1",
# 
#   "Myeloid" = "CD68",
#   "Myeloid" = "CSF1R",
#   "Myeloid" = "SPI1",
#   "Myeloid" = "CD14",
#   
#   "Mast" = "KIT",
#   "Mast" = "MS4A2",
# 
#   "T_NK" = "CD3E",
#   "T_NK" = "CD8A",
# 
#   "Normal Epi" = "EPCAM",
#   "Normal Epi" = "PAX8",
# 
#   "Tumor Epi" = "EPCAM",
#   "Tumor Epi" = "PAX8"
# )

major_marker_ls <- c(
  "CDH5", "PECAM1", "CD34",     # Endothelial
  "RGS5", "PDGFRB", "MCAM",      # Pericyte
  "OSR2", "PDGFRA",              # Stromal
  "MYLK", "CSRP1",               # SMC
  "CD68", "CSF1R", "SPI1", "CD14",  # Myeloid
  "KIT", "MS4A2",                # Mast
  "CD3E", "CD8A",                # T_NK
  "EPCAM", "PAX8"#,               # Normal Epi
  #"EPCAM", "PAX8"                # Tumor Epi
)


available_genes <- major_marker_ls[major_marker_ls %in% rownames(all_data_OE_EAOC)]
missing_genes <- setdiff(major_marker_ls, rownames(all_data_OE_EAOC))

all_data_OE_EAOC <- JoinLayers(all_data_OE_EAOC,layers = "data")

avg_expression <- AverageExpression(
  all_data_OE_EAOC,
  assays =DefaultAssay(all_data_OE_EAOC),
  features = available_genes,
  group.by = "celltype",
  slot = "data"
)[[1]]

head(avg_expression)
colnames(avg_expression)[colnames(avg_expression) == "T-NK"] <- "T_NK"

celltype_order <- c(
  "Endothelial",
  "Pericyte", 
  "Stromal",
  "SMC",        
  "Myeloid",    
  "Mast",       
  "T_NK",  
  "Epithelial"
  # "Normal Epi", 
  # "Tumor Epi"   
)

avg_expression_ordered <- avg_expression[, intersect(celltype_order, colnames(avg_expression))]
avg_expression_ordered <- avg_expression_ordered[intersect(major_marker_ls,rownames(avg_expression)), ]

library(pheatmap)
library(viridis)

pdf(file.path(result_dir, "marker_celltype_heatmap.pdf"), width = 6, height = 8)
p <- pheatmap(
  avg_expression_ordered,
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = viridis(100),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = "black",
  #annotation_row = gene_groups,
  #annotation_colors = annotation_colors,
  main = "",
  #gaps_col = c(4, 6)  # 在特定列之间添加间隙
)

print(p)

dev.off()


#图1f、注释结果空间映射------
all_data2@meta.data$project %>% unique()

# result_dir <- 'result' 
# dir.create(result_dir)

##add FOV画ImageDimPlot----
source("/mnt/NexData_file/home/hyj/h5ad2rds(1).R")
library(Seurat)
library(SeuratDisk)
library(sceasy)
library(reticulate)
library(dplyr)
library(ggplot2)
library(argparser)
library(patchwork)
library(RColorBrewer)
library(arrow)

all_data3 <-  JoinLayers(all_data2,layers = 'count')
colnames(all_data3@meta.data)

all_data3@meta.data$cell_id %>% head()

###改行名
all_data3@meta.data$cell_id1 <- gsub("^.+?_(?=[a-z]{8}-\\d+$)", "", 
                                     rownames(all_data3@meta.data), 
                                     perl = TRUE)
T1<- data_frame(rownames(all_data3@meta.data),all_data3@meta.data[c("orig.ident", "cell_id","cell_id1")])

all_data3@meta.data$cell_id2 <- rownames(all_data3@meta.data)

rownames(all_data3@meta.data) <- paste0(all_data3@meta.data$sample,"_",all_data3@meta.data$cell_id1)

T1<- data_frame(rownames(all_data3@meta.data),all_data3@meta.data[c("orig.ident", "cell_id","cell_id1","cell_id2")])

rownames(all_data3@assays$RNA@layers$count) <- rownames(all_data3)
colnames(all_data3@assays$RNA@layers$count) <- colnames(all_data3)
xenium <- CreateSeuratObject(
  all_data3@assays$RNA@layers$count,
  assay = "Xenium")

xenium@meta.data <- all_data3@meta.data

umap <- all_data3@reductions$umap@cell.embeddings
rownames(umap) <- rownames(xenium@meta.data)
metadata <- xenium@meta.data
xenium <- Add_meta_umapData_Xenium(object = xenium,umapfile = umap,metedata = metadata)
colnames(xenium@assays$Xenium) <-  colnames(xenium)
row.names(xenium@reductions$umap@cell.embeddings) <- colnames(xenium)

xenium <- RenameCells(xenium,new.names = colnames(xenium))
label <- unique(xenium@meta.data$sample)

#options(future.globals.maxSize = 11 * 1024^3)
#xenium <- NormalizeData(xenium)

xenium <- JoinLayers(xenium,layers = "data")


###多个文件夹，一个文件夹内包含多个fov###----

##加入1个测试----
library(arrow)
# coord1 <- Read_fov_spatial(data.dir = "/data/xenium/20241107__051016__20241107_FWSP20241358/output-XETG00305__0041445__23_03967__20241107__051512",assay = 'Xenium',
#                            mols.qv.threshold = 20,label = "CYS_23_03967",object = xenium,cov = NULL)
# xenium@images[["CYS-23-03967"]] = coord1

#options(future.globals.maxSize = 10 * 1024^3)
#xenium <- NormalizeData(xenium)
DefaultBoundary(xenium@images[["EAOC-3"]]) <- "segmentation"
p <- ImageFeaturePlot(xenium, fov = "EAOC-3",features = c("CD3E"),
                      max.cutoff = "q90",axes = FALSE,
                      border.color = "white", 
                      border.size = 0.1,
                      cols = c('#003366', '#FFD700'))
ggsave(p,filename = paste0("result0227/","_feature_spatial.png"),width = 12 ,height = 12,dpi = 300)


##加入16个样本FOV----
coord1 <- Read_fov_spatial(data.dir = "/data/xenium/XE2511116456-01/output-XETG00305__0076863__E2324414-1_E2437436-13_S2415628-14__20251206__030955",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "OE-3",object = xenium,cov = NULL)
xenium@images[["OE-3"]] = coord1

coord2 <- Read_fov_spatial(data.dir = "/data/xenium/XE2511116456-01/output-XETG00305__0076863__E2324414-1_E2437436-13_S2415628-14__20251206__030955",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "OE-4",object = xenium,cov = NULL)
xenium@images[["OE-4"]] = coord2

coord3 <- Read_fov_spatial(data.dir = "/data/xenium/XE2511116456-01/output-XETG00305__0076863__E2324414-1_E2437436-13_S2415628-14__20251206__030955",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "EAOC_Eu-1",object = xenium,cov = NULL)
xenium@images[["EAOC_Eu-1"]] = coord3

coord4 <- Read_fov_spatial(data.dir = "/data/xenium/XE2511116456-01/output-XETG00305__0076863__E2409744-6__20251206__030955",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "EAOC-3",object = xenium,cov = NULL)
xenium@images[["EAOC-3"]] = coord4

coord5 <- Read_fov_spatial(data.dir = "/data/xenium/XE2511116456-01/output-XETG00305__0076863__s2415628-3_5__20251206__030955",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "EAOC-1",object = xenium,cov = NULL)
xenium@images[["EAOC-1"]] = coord5

coord6 <- Read_fov_spatial(data.dir = "/data/xenium/XE2511116456-01/output-XETG00305__0076863__s2415628-3_5__20251206__030955",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "EAOC-2",object = xenium,cov = NULL)
xenium@images[["EAOC-2"]] = coord6

coord7 <- Read_fov_spatial(data.dir = "/data/xenium/XE2511116456-01/output-XETG00305__0076863__S2424653-6_ES2504309-6__20251206__030955",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "EAOC-4",object = xenium,cov = NULL)
xenium@images[["EAOC-4"]] = coord7

coord8 <- Read_fov_spatial(data.dir = "/data/xenium/XE2511116456-01/output-XETG00305__0076863__S2424653-6_ES2504309-6__20251206__030955",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "EAOC-6",object = xenium,cov = NULL)
xenium@images[["EAOC-6"]] = coord8

coord9 <- Read_fov_spatial(data.dir = "/data/xenium/XE2511116456-01/output-XETG00305__0076863__S2424653-7_E2303678-2_E2421767-6__20251206__030955",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "OE-1",object = xenium,cov = NULL)
xenium@images[["OE-1"]] = coord9

coord10 <- Read_fov_spatial(data.dir = "/data/xenium/XE2511116456-01/output-XETG00305__0076863__S2424653-7_E2303678-2_E2421767-6__20251206__030955",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "OE-2",object = xenium,cov = NULL)
xenium@images[["OE-2"]] = coord10

coord11 <- Read_fov_spatial(data.dir = "/data/xenium/XE2511116456-01/output-XETG00305__0076863__S2424653-7_E2303678-2_E2421767-6__20251206__030955",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "EAOC-5",object = xenium,cov = NULL)
xenium@images[["EAOC-5"]] = coord11

coord12 <- Read_fov_spatial(data.dir = "/data/xenium/XE2511116456-01/output-XETG00305__0076863__TMA__20251206__030955",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "EAOC_Eu-3",object = xenium,cov = NULL)
xenium@images[["EAOC_Eu-3"]] = coord12

coord13 <- Read_fov_spatial(data.dir = "/data/xenium/XE2511116456-01/output-XETG00305__0076863__TMA__20251206__030955",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "OE_Eu-2",object = xenium,cov = NULL)
xenium@images[["OE_Eu-2"]] = coord13

coord14 <- Read_fov_spatial(data.dir = "/data/xenium/XE2511116456-01/output-XETG00305__0076863__TMA__20251206__030955",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "EAOC_Eu-2",object = xenium,cov = NULL)
xenium@images[["EAOC_Eu-2"]] = coord14

coord15 <- Read_fov_spatial(data.dir = "/data/xenium/XE2511116456-01/output-XETG00305__0076863__TMA__20251206__030955",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "OE_Eu-3",object = xenium,cov = NULL)
xenium@images[["OE_Eu-3"]] = coord15

coord16 <- Read_fov_spatial(data.dir = "/data/xenium/XE2511116456-01/output-XETG00305__0076863__TMA__20251206__030955",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "OE_Eu-1",object = xenium,cov = NULL)
xenium@images[["OE_Eu-1"]] = coord16


saveRDS(xenium,"xenium_FOV.rds")
xenium <- readRDS("xenium_FOV.rds")

###Imagedimplot----
xenium@meta.data$celltype %>% unique()
celltype_level <- c("Endothelial", "Pericyte", "Stromal", "SMC", "Myeloid", "Mast", "T_NK" ,"Epithelial")
xenium@meta.data$celltype <- factor(xenium@meta.data$celltype,levels = celltype_level)

plan <- c("Endothelial" = "#005995",
          "Pericyte" = "#96CEB4",
          "Stromal" = "#248888",
          "SMC" = "#73B66B",
          "Myeloid" = "#FF6D31",
          "Mast" = "#EF5285",
          "T_NK" ="#FFC300",
          "Epithelial" = "#900475" 
)

table(xenium$sample,xenium$celltype)

# sample_label <- unique(xenium@meta.data$sample)
sample_label <- c("OE-1","OE-2","OE-3","OE-4",
                  "EAOC-1","EAOC-2","EAOC-3","EAOC-4","EAOC-5","EAOC-6")

for( i in sample_label){
  #j = gsub("_","-",i)
  DefaultBoundary(xenium@images[[i]]) <- "segmentation"
  p <- ImageDimPlot(xenium,
                    fov = i,
                    nmols = 20000,
                    dark.background = F,
                    border.color = "black", 
                    border.size = 0.1,
                    cols = plan,
                    #shuffle.cols = F, #随机打乱颜色顺序
                    #size = 1.2,
                    group.by = "celltype"
                    #molecules = c("INHBA"),
                    #mols.size = 1,
                    #mols.cols = "red",
                    #mols.alpha = 0.5, #分子的透明度值
                    
  )
  ggsave(p,filename = paste0("result0325//Imagedimplot/",i,"_anno_spatial.png"),width = 12 ,height = 12,dpi = 300)
  ggsave(p,filename = paste0("result0325//Imagedimplot/",i,"_anno_spatial.pdf"),width = 12 ,height = 12)
}

#图1g、EAOC上皮marker空间映射-------
EAOCmarker <- c("WT1","MSLN","HNF1B","PIK3CA","ARID1A") 

Imagefeature_dir <- "result0227/ImageFeaturePlot/"
dir.create(Imagefeature_dir)


for( i in sample_label){ 
  for(gene in EAOCmarker){
    DefaultBoundary(xenium@images[[i]]) <- "segmentation"
    p <- ImageFeaturePlot(xenium, fov = i,
                          features = gene,
                          max.cutoff = "q90",
                          dark.background = F,
                          axes = FALSE,
                          border.color = "black", 
                          border.size = 0.1,
                          cols = c('#003366', '#FFD700')
    )
    ggsave(p,filename = file.path(Imagefeature_dir,paste0(i, "_", gene, "_feature_spatial.png")), width = 12 ,height = 12,dpi = 300)
    ggsave(p,filename = file.path(Imagefeature_dir,paste0(i, "_", gene, "_feature_spatial.pdf")), width = 12 ,height = 12)
    
  }
}


plan <- c("Endothelial" = "#005995",
          "Pericyte" = "#96CEB4",
          "Stromal" = "#248888",
          "SMC" = "#73B66B",
          "Myeloid" = "#FF6D31",
          "Mast" = "#FDC4B6",
          "T_NK" = "#EF5285",
          "Normal Epi" = "#FFC300",
          "Tumor Epi" =  "#900475"             
)


pdf(file.path(result_dir, str_glue("Dimplot Epi subtype umap.pdf")), height = 6, width = 8)

DimPlot(all_data_OE_EAOC[,all_data_OE_EAOC$celltype %in% c("Epithelial")], 
        group.by = "subtype", 
        label = F, 
        cols = plan,
        reduction = "harmonyumap")

dev.off()

pdf(file.path(result_dir, str_glue("Dimplot Epi subtype umap group.pdf")), height = 6, width = 12)

DimPlot(all_data_OE_EAOC[,all_data_OE_EAOC$celltype %in% c("Epithelial")], 
        group.by = "subtype", 
        label = F, 
        cols = plan,
        reduction = "harmonyumap",
        split.by = "group"
)
dev.off()



#图1h、分组小提琴-------
all_data_OE_EAOC@meta.data$group %>% table()
all_data_OE_EAOC@meta.data$sample %>% table()

all_data_OE_EAOC@meta.data$group <- as.character(all_data_OE_EAOC@meta.data$group)
table(all_data_OE_EAOC@meta.data$group,all_data_OE_EAOC@meta.data$subtype)

all_data_EAOC_Epi <- subset(all_data_OE_EAOC, group %in% c("EAOC")&subtype %in% c("Normal Epi","Tumor Epi"))
all_data_EAOC_Epi@meta.data$group %>% table()
table(all_data_EAOC_Epi@meta.data$sample,all_data_EAOC_Epi@meta.data$subtype)

ViolinPlot_dir <- "result0325/ViolinPlot"
dir.create(ViolinPlot_dir)

pdf(file.path(ViolinPlot_dir, "EAOC_markers_by_sample.pdf"), width = 16, height = 5)
p2 <- VlnPlot(all_data_OE_EAOC,
              features = EAOCmarker,
              group.by = "sample", 
              pt.size = 0,
              split.by = "subtype",  
              #split.plot = TRUE,
              ncol = length(EAOCmarker)#,
              #cols = c("#E69F00", "#56B4E9"),
              #same.y.lims = T
)
print(p2)
dev.off()

p3_ls <-
  EAOCmarker %>%
  setNames(., .) %>% 
  map(function(x){
    idx <- FetchData(all_data_OE_EAOC, vars = x) %>% {rownames(.)[.[[1]] != 0] }
    seu <- all_data_OE_EAOC[, idx]
    
    p <- VlnPlot(seu,
                 features = x,
                 group.by = "sample",  # 按样本分组
                 pt.size = 0
    )
    p
  })

p3 <- patchwork::wrap_plots(p3_ls, ncol = 5)

p3

##数据处理-----
expr_data <- FetchData(
  object = all_data_EAOC_Epi,
  vars = c(EAOCmarker, "sample", "subtype") 
)

#长格式，适配ggplot2
expr_long <- expr_data %>%
  pivot_longer(
    cols = all_of(EAOCmarker),
    names_to = "gene",         
    values_to = "expression"   
  ) %>%
  filter(expression != 0) %>%
  mutate(
    sample = factor(sample, levels = c("EAOC-1","EAOC-2","EAOC-3","EAOC-4","EAOC-5","EAOC-6","OE-1","OE-2","OE-3","OE-4")),
    subtype = factor(subtype)
  )

group_colors <- c("Normal Epi" = "#A23B72", "Tumor Epi" = "#2E86AB")

##所有marker----
p <- ggplot(expr_long, aes(x = sample, y = expression, fill = subtype)) +
  geom_violin(scale = "width", alpha = 0.8, position = position_dodge(width = 0.8)) +
  #geom_boxplot(width = 0.1, alpha = 1, position = position_dodge(width = 0.8)) +
  facet_wrap(~gene, scales = "free_y", ncol = 1) + 
  scale_fill_manual(values = group_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 10, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Sample",
    y = "Gene Expression",
    fill = "subtype",
    title = ""
  )


ggsave(
  filename = "result0325/ViolinPlot/EAOC_marker_violin_plot.pdf",
  plot = p,
  width = 8,  
  height = 8, 
  dpi = 300,   
  device = "pdf"
)

ggsave(
  filename = "result0325/ViolinPlot/EAOC_marker_violin_plot.png",
  plot = p,
  width = 8,
  height = 8,
  dpi = 300
)

##单独marker-sample-----
for (i in EAOCmarker) {
  expr <- expr_long %>% 
    filter(gene == i)
  
  p <- ggplot(expr, aes(x = sample, y = expression, fill = subtype)) +
    geom_violin(scale = "width", alpha = 0.8, position = position_dodge(width = 0.8)) +
    #geom_boxplot(width = 0.1, alpha = 1, position = position_dodge(width = 0.8)) +
    # facet_wrap(~gene, scales = "free_y", ncol = 1) + 
    scale_fill_manual(values = group_colors) +
    theme_bw() +
    theme(
      
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 10, face = "bold"),
      
      strip.text = element_text(size = 10, face = "bold"),
      
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 8),
      
      panel.grid = element_blank()
    ) +
    labs(
      x = "Sample",
      y = "Gene Expression",
      fill = "subtupe",
      title = ""
    )
  
  
  ggsave(
    filename = paste0("result0325/ViolinPlot/",i,"_violin_plot.pdf"),
    plot = p,
    width = 8,  
    height = 3, 
    dpi = 300,   
    device = "pdf"
  )
  
  ggsave(
    filename = paste0("result0325/ViolinPlot/",i,"_violin_plot.png"),
    plot = p,
    width = 8,
    height = 3,
    dpi = 300
  )
}





# library(ggpubr)
##单独marker-group-----
# for (i in EAOCmarker) {
#   expr <- expr_long %>% 
#     filter(gene == i)
#   
#   p <- ggplot(expr, aes(x = group, y = expression, fill = subtype)) +
#     geom_violin(scale = "width", alpha = 0.8, position = position_dodge(width = 0.8)) +
#     geom_boxplot(width = 0.1, alpha = 1, position = position_dodge(width = 0.8)) +
#     # facet_wrap(~gene, scales = "free_y", ncol = 1) + 
#     scale_fill_manual(values = group_colors) +
#     stat_compare_means(method = "t.test",size = 8, label = "p.signif", label.y = 2.2, label.x = 1.4
#     )+
#     theme_bw() +
#     theme(
#       
#       axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
#       axis.text.y = element_text(size = 8),
#       axis.title = element_text(size = 10, face = "bold"),
#       
#       strip.text = element_text(size = 10, face = "bold"),
#       
#       legend.title = element_text(size = 10, face = "bold"),
#       legend.text = element_text(size = 8),
#       
#       panel.grid = element_blank()
#     ) +
#     labs(
#       x = "Sample",
#       y = "Gene Expression",
#       fill = "subtype",
#       title = ""
#     )
#   
#   
#   ggsave(
#     filename = paste0("result0227/ViolinPlot/",i,"_group_violin_plot.pdf"),
#     plot = p,
#     width = 5,  
#     height = 5, 
#     dpi = 300,   
#     device = "pdf"
#   )
#   
#   ggsave(
#     filename = paste0("result0227/ViolinPlot/",i,"_group_violin_plot.png"),
#     plot = p,
#     width = 5,
#     height = 5,
#     dpi = 300
#   )
# }

#图2a、OE-EAOC 细胞比例变化图------
library(tidyverse)
library(Seurat)
library(ggalluvial)
library(ggthemes)
library(ggpubr)

dat_props_immune <- all_data_OE_EAOC[[]] %>% 
  group_by(group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

dat_props_immune %>% head
# # A tibble: 6 × 4
# # Groups:   group [1]
# group annotation  count   freq
# <fct> <fct>       <int>  <dbl>
#   1 Ctrl  B cell       1437 0.153 
# 2 Ctrl  Plasma cell  1515 0.162 
# 3 Ctrl  T cell       4222 0.451 
# 4 Ctrl  NK cell       152 0.0162
# 5 Ctrl  Macrophage   1192 0.127 
# 6 Ctrl  Neutrophil    253 0.0270

# p_props_immune  <-
#   dat_props_immune %>%
#   ggplot(aes(x = group, y = freq, fill = subtype)) +
#   geom_col(width = 0.95) +
#   theme_classic() +
#   theme(
#     axis.line = element_blank()
#   ) +
#   scale_x_discrete(expand = c(0, 0)) +
#   scale_fill_manual(values = plan) +
#   scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
#   labs(
#     x = NULL,
#     y = NULL
#   )

new_row <- tibble(
  group   = "OE",
  subtype = factor("Tumor Epi", levels = levels(dat_props_immune$subtype)),
  count   = 1L,
  freq    = 0.00001
)
dat_props_immune <- bind_rows(dat_props_immune, new_row)

p_props_immune_flow <- 
  dat_props_immune %>% 
  # mutate(id = 1:n()) %>%
  
  ggplot(aes(x = group, y = freq, stratum = subtype, alluvium = subtype, fill = subtype, label = subtype)) +
  geom_stratum(alpha = 1, color = NA) +
  geom_flow(alpha = 0.5) +
  #geom_text(stat = "stratum") +
  theme_classic2() +
  theme(
    axis.line = element_blank(),
    axis.text = element_text(size = 12)
  ) +
  scale_fill_manual(values = plan) +
  #scale_x_discrete(expand = c(0.01, 0.01)) +
  #scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(
    x = "group",
    y = "cell proportion"
  )

ggsave(file = "result0325/props_flow.pdf", p_props_immune_flow,width = 6,height = 6)
ggsave(file = "result0325/props_flow.png", p_props_immune_flow,width = 6,height = 6)




#图2d、OE_上皮 → EAOC_ Normal → EAOC_ Tumor 三组周围25微米范围内细胞比例统计（邻近度）------
all_data_OE_EAOC <- readRDS("result0325/all_data_OE_EAOC.rds")
all_data_OE_EAOC@meta.data %>% head()
rownames(all_data_OE_EAOC@meta.data) %>% head()
all_data_OE_EAOC@meta.data$subtype %>% table()

group_subtype_level <- c("OE_Normal Epi","EAOC_Normal Epi","EAOC_Tumor Epi","Endothelial", "Pericyte", "Stromal", "SMC", "Myeloid", "Mast", "T_NK")

all_data_OE_EAOC@meta.data <- all_data_OE_EAOC@meta.data %>% 
  mutate(
    group_subtype = case_when(
      subtype %in% c("Normal Epi","Tumor Epi") ~ paste0(group, "_", subtype),
      TRUE ~ subtype  
    ),
    group_subtype = factor(group_subtype,levels=group_subtype_level)
  )

unique(all_data_OE_EAOC$group_subtype)

table(all_data_OE_EAOC$group, all_data_OE_EAOC$subtype)
table(all_data_OE_EAOC$group, all_data_OE_EAOC$group_subtype)


cells_boundary_ls <- readRDS(file.path("~/R_OUTPUT/data_for_analysis", "cells_boundary_ls.rds"))
cells_position_ls <- readRDS(file.path("~/R_OUTPUT/data_for_analysis", "cells_position_ls.rds"))

library(Seurat)
library(tidyverse)

source("spatial analysis.R")

# proximity_out_dir <- "proximity_out_dir"
# dir.create(proximity_out_dir)

###sample名称和slide名称需一致，slide名称是cells_position_ls各FOV的名字
## Proximity analysis -------
## 1. input data ------------
cells_position_ls %>% head() 

# all_data_cn_meta
all_data_meta_pos <-
  cells_position_ls %>%
  imap(function(cells_position, slide){ # browser()
    
    all_data_meta <- all_data_OE_EAOC[[]]  
    
    proj_meta_sub <-
      all_data_meta %>%  
      dplyr::select(group_subtype, sample, "group") %>% 
      rownames_to_column("cell_id") %>%
      dplyr::filter(sample == slide)
    
    cells_position[, 1:3] %>%
      inner_join(proj_meta_sub, by = "cell_id") %>%  
      mutate(group_subtype = factor(group_subtype,levels=group_subtyep_level)) %>%
      dplyr::rename(Classification1 = group_subtype,
                    X = x_centroid,
                    Y = y_centroid,
                    Slide_ID = sample,
                    Group = group
      )
  }) %>%
  enframe() %>%
  unnest(cols = value) %>%
  dplyr::select(-name) %>%
  column_to_rownames("cell_id")

all_data_meta_pos %>% head()
all_data_meta_pos$Slide_ID %>% table()
all_data_meta_pos$Group %>% table()
all_data_meta_pos$Classification1 %>% table()

all_data_meta_pos$cell_id <- rownames(all_data_meta_pos)
all_data_meta_pos <- all_data_meta_pos %>%
  as.data.frame() %>%       
  remove_rownames()    

write.csv(all_data_meta_pos,"result0325/proximity analysis/all_data_meta_pos.csv")

all_data_meta_pos <- read.csv("result0325/proximity analysis/all_data_meta_pos.csv")

# all_data_meta_pos_subset <- readRDS('data/all_data_meta_pos_subset.rds')
# all_data_meta_pos_subset %>% head

#        cell_id Classification1   Group           Slide_ID        X        Y
# 1 3_aaaacdgo-1     Fibroblasts   MUT    20231128__R01-001 2228.621 860.3357
# 2 3_aaaanalm-1   Myeloid cells   MUT    20231128__R01-001 2248.922 850.1516
# 3 3_aaaaonpl-1   Myeloid cells   MUT    20231128__R01-001 2270.447 856.8892
# 4 3_aaabciio-1         B cells   MUT    20231128__R01-001 2279.562 838.6568
# 5 3_aaabdkea-1     Fibroblasts   MUT    20231128__R01-001 2292.319 858.3909
# 6 3_aaacdgip-1         B cells   MUT    20231128__R01-001 2288.924 868.3226
# ......


## 2. proximity analysis -------
# 参数介绍
# x,y 指定细胞的坐标
# group.by 一般是切片或者样本编号，代表分别计算最近邻关系
# classification 一般是细胞的注释
# k 计算的最近邻的细胞数
# searchtype 邻近度使用radius

# radius是只查看当前细胞的多少um内的细胞
# mpp 坐标与um的转换比例，xenium是1，celldive是0.325
# returntype 可以每个细胞多少radius内的细胞计数count，也可以返回细胞频率freq
# cell_enough_ratio_cutoff 为了降低计算时间，k值不建议设置太高，
# 如果k值未能完全覆盖全部的邻近细胞，则cell_enough_ratio_cutoff会自动加k值，
# 直到默认90%的细胞均已经计算了全部邻近细胞

### 2.1 proximity analysis --------
res_pa <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
    k = 750, # 'auto',
    searchtype = "radius",
    radius = 25, #  50um半径
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )


# return res_pa
# rows: all celltypes
# cols: query cell
res_pa_colsum <- colSums(res_pa)
res_pa <- sweep(res_pa, 2, res_pa_colsum, FUN = "/")
if(any(res_pa_colsum == 0)){
  res_pa[, res_pa_colsum == 0] <- 0
}

# res_pa <- readRDS('data/res_pa.rds')


# summarized by Group and celltype
# pa_groups <-
#   with(
#     all_data_meta_pos,
#     str_c(Classification1, Group, Slide_ID, sep = '..')
#   )


pa_groups <-
  all_data_meta_pos %>%
  group_by(Slide_ID) %>% # compute pa_groups per Slide_ID
  nest() %>%
  deframe() %>%
  imap(function(dat, Slide_ID){
    with(
      dat,
      str_c(Classification1, Group, Slide_ID, sep = '..')
    )
  }) %>%
  unlist(use.names = FALSE)

# 分组内取细胞比例的中位数
# apply对每行操作时会转置结果
res_pa_df <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups, FUN = "median")})

res_pa_df  <-
  res_pa_df %>%
  as.data.frame() %>%
  rownames_to_column('group') %>%
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
  separate(group, into = c("Target", "Group", "project"), sep = "\\.\\.") 

# plot boxplot
# Target:中心细胞
# Others: 其他细胞占比

p_res_pa <-
  res_pa_df %>%
  # dplyr::filter(Target == "Epithelial cells") %>%
  ggplot(aes(x = Others, y = cell_props, fill = Group)) +
  geom_boxplot(width = 0.35, outlier.colour = NA) +
  geom_jitter(width = 0.25) +
  facet_wrap(~ Target, scales = "free_y",ncol = 5)+
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# export plot
pdf(
  file.path("result0325/proximity analysis", "proximity analysis.pdf"),
  height = 10,
  width  = 20
)
p_res_pa %>% print()
dev.off()

table(res_pa_df$Target,res_pa_df$Group)


###循环OE_上皮 → EAOC_ Normal → EAOC_ Tumor-------
cell_list <- c("OE_Normal Epi", "EAOC_Normal Epi", "EAOC_Tumor Epi")
lapply(cell_list,function(cell){
  p_res_pa <-
    res_pa_df %>%
    dplyr::filter(Target == cell) %>%
    ggplot(aes(x = Others, y = cell_props, fill = Group)) +
    geom_boxplot(width = 0.35, outlier.colour = NA) +
    geom_jitter(width = 0.25) +
    facet_wrap(~ Target, scales = "free_y",ncol = 5)+
    ggpubr::theme_pubclean() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  # export plot
  pdf(
    file.path("result0325/proximity analysis", paste0(cell,"_proximity analysis.pdf")),
    height = 5,
    width  = 5
  )
  p_res_pa %>% print()
  dev.off()
})


##统计25μm范围内的细胞数量------
res_pa <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
    k = 750, # 'auto',
    searchtype = "radius",
    radius = 25, #  50um半径
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )

pa_groups <-
  all_data_meta_pos %>%
  group_by(Slide_ID) %>% # compute pa_groups per Slide_ID
  nest() %>%
  deframe() %>%
  imap(function(dat, Slide_ID){
    with(
      dat,
      str_c(Classification1, Group, Slide_ID, sep = '..')
    )
  }) %>%
  unlist(use.names = FALSE)

res_pa_df <- apply(res_pa, 1, function(x) {
  tapply(x, INDEX = pa_groups, FUN = median, na.rm = TRUE)  # 或用 mean
})

res_pa_df  <-
  res_pa_df %>%
  as.data.frame() %>%
  rownames_to_column('group') %>%
  pivot_longer(-group, names_to = "Others", values_to = "cell_counts") %>%  # 列名改为 cell_counts
  separate(group, into = c("Target", "Group", "project"), sep = "\\.\\.") 

res_pa_df$Target = factor(res_pa_df$Target,levels=group_subtype_level)
res_pa_df$Others = factor(res_pa_df$Others,levels=group_subtype_level)

p_res_pa <-
  res_pa_df %>%
  ggplot(aes(x = Others, y = cell_counts, fill = Group)) +  # 注意：y = cell_counts
  geom_boxplot(width = 0.35, outlier.colour = NA) +
  geom_jitter(width = 0.25, alpha = 0.6) +
  facet_wrap(~ Target, scales = "free_y", ncol = 5) +
  ggpubr::theme_pubclean() +
  labs(y = "Number") +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(
  file.path("result0325/proximity analysis", "proximity_analysis_cell_counts.pdf"),
  plot = p_res_pa,
  height = 10,
  width = 20
)


### 循环绘制特定上皮细胞-----
cell_list <- c("OE_Normal Epi", "EAOC_Normal Epi", "EAOC_Tumor Epi")
lapply(cell_list, function(cell) {
  p_sub <- res_pa_df %>%
    dplyr::filter(Target == cell) %>%
    ggplot(aes(x = Others, y = cell_counts, fill = Group)) +
    geom_boxplot(width = 0.35, outlier.colour = NA) +
    geom_jitter(width = 0.25, alpha = 0.6) +
    facet_wrap(~ Target, scales = "free_y", ncol = 5) +
    ggpubr::theme_pubclean() +
    labs(y = "Number of cells within 25 μm") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  ggsave(
    file.path("result0325/proximity analysis", paste0(cell, "_proximity_analysis_counts.pdf")),
    plot = p_sub,
    height = 5,
    width = 5
  )
})

###三种细胞合并展示------
res_pa_df$Target %>% unique()
res_pa_df$Others %>% unique()
plan <- c("Endothelial" = "#005995",
          "Pericyte" = "#96CEB4",
          "Stromal" = "#248888",
          "SMC" = "#73B66B",
          "Myeloid" = "#FF6D31",
          "Mast" = "#FDC4B6",
          "T_NK" = "#EF5285",
          "OE_Normal Epi" =  "#FBF579",
          "EAOC_Normal Epi" = "#FFC300",
          "EAOC_Tumor Epi" =  "#900475"             
)

p_res_pa <-
  res_pa_df %>%
  dplyr::filter(Target %in% cell_list) %>%
  ggplot(aes(x = Target, y = cell_counts, fill = Others)) + 
  geom_col(position = "stack") +  # 堆叠柱状图
  #geom_boxplot(width = 0.35, outlier.colour = NA) +
  #geom_jitter(width = 0.25, alpha = 0.6) +
  #facet_wrap(~ Target, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = plan) +
  ggpubr::theme_pubclean() +
  labs(y = "Number", fill = "celltype") +  
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1#, #vjust = 0.5
                                   ),
        legend.position = "right" )
 
ggsave(
  file.path("result0325/proximity analysis", "proximity_analysis_3type_cell_counts.pdf"),
  plot = p_res_pa,
  height = 5,
  width = 5
)



###三种细胞合并展示(桑基图)------
library(tidyverse)
library(Seurat)
library(ggalluvial)
library(ggthemes)
library(ggpubr)

dat_props <- res_pa_df %>%
  dplyr::filter(Target %in% cell_list) %>%
  mutate(freq = cell_counts/sum(cell_counts),
         celltype = Others,
) %>% 
  group_by(Target, celltype) %>%
  summarise(cell_counts = sum(cell_counts), .groups = "drop")

dat_props_immune %>% head
# # A tibble: 6 × 4
# # Groups:   group [1]
# group annotation  count   freq
# <fct> <fct>       <int>  <dbl>
#   1 Ctrl  B cell       1437 0.153 
# 2 Ctrl  Plasma cell  1515 0.162 
# 3 Ctrl  T cell       4222 0.451 
# 4 Ctrl  NK cell       152 0.0162
# 5 Ctrl  Macrophage   1192 0.127 
# 6 Ctrl  Neutrophil    253 0.0270


p_props_flow <- 
  dat_props %>% 
  # mutate(id = 1:n()) %>%
  
  ggplot(aes(x = Target, y = cell_counts, stratum = celltype, alluvium = celltype, fill = celltype, label = celltype)) +
  geom_stratum(alpha = 1, color = NA) +
  geom_flow(alpha = 0.5) +
  #geom_text(stat = "stratum") +
  theme_classic2() +
  theme(
    axis.line = element_blank(),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1#, #vjust = 0.5
    ),
    legend.position = "right")+
  scale_fill_manual(values = plan) +
  #scale_x_discrete(expand = c(0.01, 0.01)) +
  #scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(
    x = "group",
    y = "cell proportion"
  )

ggsave(file = "result0325/proximity analysis/3type_cell_counts_flow.pdf", p_props_flow,width = 6,height = 6)
ggsave(file = "result0325/proximity analysis/3type_cell_counts_flow.png", p_props_flow,width = 6,height = 6)






#图2e、f、g、h、邻域分析（加umap）、空间映射、CN种细胞比例变化统计-------
##邻域分析-----
all_data_OE_EAOC@meta.data %>% head()
CN_data <- all_data_OE_EAOC

library(Seurat)
library(tidyverse)
source("spatial analysis.R")
#####data ------
CN_out_dir <- "result0227/CN analysis"
dir.create(CN_out_dir)

#all_data <- readRDS("data_for_analysis/all_data_annotation_rctd.rData")
Idents(CN_data) <- "subtype"

cells_position_ls<-readRDS(file.path("~/R_OUTPUT/data_for_analysis", "cells_position_ls.rds"))
cells_position_df <- cells_position_ls %>% bind_rows()

#meta.data$cell_id需要和cells_position_ls里的cell_id一致

# CN_data@meta.data$cell_id2 <- rownames(CN_data@meta.data)
# 
# CN_data@meta.data$cell_id1 <- CN_data@meta.data$cell_id
# 
# CN_data@meta.data$cell_id<- paste(
#   CN_data$project,
#   CN_data$cell_id1,
#   sep = "_"
# )
# 
# rownames(CN_data@meta.data) <- CN_data@meta.data$cell_id



##### CN analysis ---------
# options(warn = 0)
# Classification1 must be a factor
# add a column: obj_id

all_data_submeta_for_cn_analysis <- 
  CN_data[[c("sample", "group", "subtype")]] %>% 
  rownames_to_column("cell_id") %>%
  left_join(cells_position_df[c('cell_id', 'x_centroid', 'y_centroid')], by = "cell_id")

colnames(all_data_submeta_for_cn_analysis)[colnames(all_data_submeta_for_cn_analysis) %in% c("x_centroid", "y_centroid")] <- c("X", "Y")

all_data_submeta_for_cn_analysis <-all_data_submeta_for_cn_analysis %>% 
  mutate(
    subtype = factor(subtype), # import!
    obj_id = cell_id # import!
  )

saveRDS(all_data_submeta_for_cn_analysis,"result0227/CN analysis/all_data_submeta_for_cn_analysis.rds")


#####CN k10 ---------
res_cn_10 <- CN_compute_CN(
  all_data_submeta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "sample",
  classification = "subtype",
  k.nn           = 10,
  k.cluster      = 10
)

saveRDS(res_cn_10, file.path(CN_out_dir, "res_cn_10.rds"))

res_cn_10 <- readRDS(file = file.path("result0227/CN analysis/","res_cn_10.rds"))

#####heatmap ----------------------
# FoldChange cluster/allCell
cn10_cluster_mean <- res_cn_10$cluster_centers
cn10_baseMean     <- res_cn_10$base_mean

# log2(FC)
cn10_cluster_FC <- log2(1 + sweep(cn10_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn10_baseMean))

# heatmap
p_cn10_heat <-
  pheatmap::pheatmap(
    cn10_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn10_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_10_heatmap.pdf"),
  height = 6,
  width  = 10
)
cowplot::plot_grid(p_cn10_heat$gtable)
dev.off()


####CN annotaton --------------------
# # CN k10
# cn_annotation <- c(
#   "CN1"  = "CAF",
#   "CN2"  = "myelo_TB",
#   "CN3"  = "CAF_Stroma_myelo",
#   "CN4"  = "Tumor_Lum",
#   "CN5"  = "TB",
#   "CN6"  = "Plasma_T",
#   "CN7"  = "stroma",
#   "CN8"  = "Tumor_immune",
#   "CN9"  = "Tumor",
#   "CN10" = "Epi"
# )
# 
# p_cn_heat2 <-
#   cn10_cluster_FC %>%
#   pheatmap::pheatmap(
#     color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
#     border_color = NA,
#     treeheight_row = 0,
#     labels_row = cn_annotation,
#     # clustering_method = 'ward.D2',
#     cluster_rows = F,
#     scale = "column",
#     silent = TRUE
#   )
# 
# pdf(
#   file.path(CN_out_dir, "CN_heatmap_annotation.pdf"),
#   height = 6,
#   width  = 11
# )
# cowplot::plot_grid(p_cn_heat2$gtable)
# dev.off()


#####CN10 spatial plot -----------------
cn10_sp_out_dir <- file.path(CN_out_dir, 'CN10 spatial plot') 
dir.create(cn10_sp_out_dir)

cells_boundary_ls<-readRDS(file.path("~/R_OUTPUT/data_for_analysis", "cells_boundary_ls.rds"))
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

CN_data@meta.data$cell_id <- NULL

All_type_meta_cn10 <- 
  CN_data[[]] %>%
  rownames_to_column('cell_id') %>%
  # left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id") %>%
  # mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
  left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id")

saveRDS(All_type_meta_cn10, 'result0227/CN analysis/all_data_meta_cn10.rds')



xenium@meta.data %>% head()

xenium@meta.data <- xenium@meta.data %>%
  rownames_to_column("cell_id_rowname") %>% 
  left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id2"), by = "cell_id2")

rownames(xenium@meta.data) <- xenium@meta.data$cell_id_rowname

xenium@meta.data$CN %>% table()

###Imagedimplot----
plan <- c("CN1" = "#005995",
          "CN2" = "#96CEB4",
          "CN3" = "#248888",
          "CN4" = "#73B66B",
          "CN5" = "#FF6D31",
          "CN6" = "#FDC4B6",
          "CN7" = "#EF5285",
          "CN8" = "#FFC300",
          "CN9" = "#900475",
          "CN10" = "gray90"
)



# sample_label <- unique(xenium@meta.data$sample)
sample_label <- c("OE-1","OE-2","OE-3","OE-4",
                  "EAOC-1","EAOC-2","EAOC-3","EAOC-4","EAOC-5","EAOC-6")

for( i in sample_label){
  #j = gsub("_","-",i)
  DefaultBoundary(xenium@images[[i]]) <- "segmentation"
  p <- ImageDimPlot(xenium,
                    fov = i,
                    nmols = 20000,
                    dark.background = F,
                    border.color = "black", 
                    border.size = 0.1,
                    cols = plan,
                    #shuffle.cols = F, #随机打乱颜色顺序
                    #size = 1.2,
                    group.by = "CN"
                    #molecules = c("INHBA"),
                    #mols.size = 1,
                    #mols.cols = "red",
                    #mols.alpha = 0.5, #分子的透明度值
                    
  )
  ggsave(p,filename = paste0("result0227/CN analysis/CN10 spatial plot/",i,"_anno_spatial.png"),width = 12 ,height = 12,dpi = 300)
  ggsave(p,filename = paste0("result0227/CN analysis/CN10 spatial plot/",i,"_anno_spatial.pdf"),width = 12 ,height = 12)
}


#####cell props ------
All_type_cn10_count_df <-
  All_type_meta_cn10 %>% 
  group_by(sample, group, CN) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

# All_type_cn10_count_df <- factor

p_cn10_props <-
  All_type_cn10_count_df %>%
  
  ggplot(aes(x = sample, y = freq, fill = CN)) +
  geom_col() +
  # ggsci::scale_fill_aaas() +
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN10_barplot.pdf"),
  height = 6,
  width  = 10
)
p_cn10_props %>% print
dev.off()

Table1 <- data_frame(All_type_meta_cn10$CN,All_type_meta_cn10$sample,All_type_meta_cn10$group)


All_type_cn10_count_df <-
  All_type_meta_cn10 %>% 
  group_by(sample, group, CN) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))


p_cn10_props_group <-
  All_type_cn10_count_df %>%
  
  ggplot(aes(x = group, y = freq, fill = group)) +
  ###统计
  stat_compare_means(label = "p.format",size=3,aes(group = group))+
  geom_boxplot() +
  ggsci::scale_fill_aaas(
    name = "Group"#,  # 图例标题
  )+
  # theme(axis.text.x = element_blank()
  # )+
  # labs(
  #   x = NULL
  # )+
  facet_wrap(~CN, nrow = 1)


pdf(
  file.path(CN_out_dir, "CN10_boxplot_group.pdf"),
  height = 6,
  width  = 15
)

p_cn10_props_group %>% print
dev.off()

##CN_umap-------
all_data_OE_EAOC@meta.data %>% head()
all_data_OE_EAOC@meta.data$cell_id %>% head()

all_data_OE_EAOC1 <- all_data_OE_EAOC
all_data_OE_EAOC1@meta.data <- all_data_OE_EAOC1@meta.data %>%
  rownames_to_column("cell_id2") %>% 
  left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id2"), by = "cell_id2") %>% 
  column_to_rownames("cell_id2")

all_data_OE_EAOC1@meta.data$CN %>% table()
is.na(all_data_OE_EAOC1@meta.data$CN) %>% table()
all_data_OE_EAOC1@meta.data %>% head()


plan <- c("CN1" = "#005995",
          "CN2" = "#96CEB4",
          "CN3" = "#248888",
          "CN4" = "#73B66B",
          "CN5" = "#FF6D31",
          "CN6" = "#FDC4B6",
          "CN7" = "#EF5285",
          "CN8" = "#FFC300",
          "CN9" = "#900475",
          "CN10" = "gray90"
)


pdf(file.path(CN_out_dir, str_glue("Dimplot CN umap.pdf")), height = 6, width = 8)
DimPlot(all_data_OE_EAOC1, 
        group.by = "CN", 
        label = TRUE, 
        cols = plan,
        reduction = "harmonyumap")

dev.off()


pdf(file.path(CN_out_dir, str_glue("Dimplot CN umap group.pdf")), height = 6, width = 12)
DimPlot(all_data_OE_EAOC1, 
        group.by = "CN", 
        label = TRUE, 
        cols = plan,
        reduction = "harmonyumap",
        split.by = "group"
)
dev.off()


#图2i、去掉EU、样本细胞比例统--------
p_col_proj <-
  all_data_OE_EAOC@meta.data %>%
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
  filename = file.path("result0227/fig2i_cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)




#20260415analysis-----
#图2d、OE_上皮 → EAOC_ Normal → EAOC_ Tumor 三组周围25微米范围内细胞比例统计（邻近度）------
all_data_OE_EAOC <- readRDS("result0325/all_data_OE_EAOC.rds")
all_data_OE_EAOC@meta.data %>% head()
rownames(all_data_OE_EAOC@meta.data) %>% head()
all_data_OE_EAOC@meta.data$subtype %>% table()

group_subtype_level <- c("OE_Normal Epi","EAOC_Normal Epi","EAOC_Tumor Epi","Endothelial", "Pericyte", "Stromal", "SMC", "Myeloid", "Mast", "T_NK")

all_data_OE_EAOC@meta.data <- all_data_OE_EAOC@meta.data %>% 
  mutate(
    group_subtype = case_when(
      subtype %in% c("Normal Epi","Tumor Epi") ~ paste0(group, "_", subtype),
      TRUE ~ subtype  
    ),
    group_subtype = factor(group_subtype,levels=group_subtype_level)
  )

unique(all_data_OE_EAOC$group_subtype)

table(all_data_OE_EAOC$group, all_data_OE_EAOC$subtype)
table(all_data_OE_EAOC$group, all_data_OE_EAOC$group_subtype)


cells_boundary_ls <- readRDS(file.path("~/R_OUTPUT/data_for_analysis", "cells_boundary_ls.rds"))
cells_position_ls <- readRDS(file.path("~/R_OUTPUT/data_for_analysis", "cells_position_ls.rds"))

library(Seurat)
library(tidyverse)

source("spatial analysis.R")

proximity_out_dir <- "result0415/proximity_out_dir"
dir.create(proximity_out_dir)

###sample名称和slide名称需一致，slide名称是cells_position_ls各FOV的名字
## Proximity analysis -------
## 1. input data ------------
cells_position_ls %>% head() 

# all_data_cn_meta
all_data_meta_pos <-
  cells_position_ls %>%
  imap(function(cells_position, slide){ # browser()
    
    all_data_meta <- all_data_OE_EAOC[[]]  
    
    proj_meta_sub <-
      all_data_meta %>%  
      dplyr::select(group_subtype, sample, "group") %>% 
      rownames_to_column("cell_id") %>%
      dplyr::filter(sample == slide)
    
    cells_position[, 1:3] %>%
      inner_join(proj_meta_sub, by = "cell_id") %>%  
      mutate(group_subtype = factor(group_subtype,levels=group_subtype_level)) %>%
      dplyr::rename(Classification1 = group_subtype,
                    X = x_centroid,
                    Y = y_centroid,
                    Slide_ID = sample,
                    Group = group
      )
  }) %>%
  enframe() %>%
  unnest(cols = value) %>%
  dplyr::select(-name) %>%
  column_to_rownames("cell_id")

all_data_meta_pos %>% head()
all_data_meta_pos$Slide_ID %>% table()
all_data_meta_pos$Group %>% table()
all_data_meta_pos$Classification1 %>% table()

all_data_meta_pos$cell_id <- rownames(all_data_meta_pos)
all_data_meta_pos <- all_data_meta_pos %>%
  as.data.frame() %>%       
  remove_rownames()    

write.csv(all_data_meta_pos,file.path(proximity_out_dir,"all_data_meta_pos.csv"))

all_data_meta_pos <- read.csv("result0415/proximity analysis/all_data_meta_pos.csv")

# all_data_meta_pos_subset <- readRDS('data/all_data_meta_pos_subset.rds')
# all_data_meta_pos_subset %>% head

#        cell_id Classification1   Group           Slide_ID        X        Y
# 1 3_aaaacdgo-1     Fibroblasts   MUT    20231128__R01-001 2228.621 860.3357
# 2 3_aaaanalm-1   Myeloid cells   MUT    20231128__R01-001 2248.922 850.1516
# 3 3_aaaaonpl-1   Myeloid cells   MUT    20231128__R01-001 2270.447 856.8892
# 4 3_aaabciio-1         B cells   MUT    20231128__R01-001 2279.562 838.6568
# 5 3_aaabdkea-1     Fibroblasts   MUT    20231128__R01-001 2292.319 858.3909
# 6 3_aaacdgip-1         B cells   MUT    20231128__R01-001 2288.924 868.3226
# ......


## 2. proximity analysis -------
# 参数介绍
# x,y 指定细胞的坐标
# group.by 一般是切片或者样本编号，代表分别计算最近邻关系
# classification 一般是细胞的注释
# k 计算的最近邻的细胞数
# searchtype 邻近度使用radius

# radius是只查看当前细胞的多少um内的细胞
# mpp 坐标与um的转换比例，xenium是1，celldive是0.325
# returntype 可以每个细胞多少radius内的细胞计数count，也可以返回细胞频率freq
# cell_enough_ratio_cutoff 为了降低计算时间，k值不建议设置太高，
# 如果k值未能完全覆盖全部的邻近细胞，则cell_enough_ratio_cutoff会自动加k值，
# 直到默认90%的细胞均已经计算了全部邻近细胞

### 2.1 proximity analysis --------
res_pa <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
    k = 750, # 'auto',
    searchtype = "radius",
    radius = 25, #  50um半径
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )


# return res_pa
# rows: all celltypes
# cols: query cell
res_pa_colsum <- colSums(res_pa)
res_pa <- sweep(res_pa, 2, res_pa_colsum, FUN = "/")
if(any(res_pa_colsum == 0)){
  res_pa[, res_pa_colsum == 0] <- 0
}

# res_pa <- readRDS('data/res_pa.rds')


# summarized by Group and celltype
# pa_groups <-
#   with(
#     all_data_meta_pos,
#     str_c(Classification1, Group, Slide_ID, sep = '..')
#   )


pa_groups <-
  all_data_meta_pos %>%
  group_by(Slide_ID) %>% # compute pa_groups per Slide_ID
  nest() %>%
  deframe() %>%
  imap(function(dat, Slide_ID){
    with(
      dat,
      str_c(Classification1, Group, Slide_ID, sep = '..')
    )
  }) %>%
  unlist(use.names = FALSE)

# 分组内取细胞比例的中位数
# apply对每行操作时会转置结果
res_pa_df <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups, FUN = "median")})



res_pa_df  <-
  res_pa_df %>%
  as.data.frame() %>%
  rownames_to_column('group') %>%
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
  separate(group, into = c("Target", "Group", "project"), sep = "\\.\\.") 


res_pa_df$Target = factor(res_pa_df$Target,levels=group_subtype_level)
res_pa_df$Others = factor(res_pa_df$Others,levels=group_subtype_level)
# plot boxplot
# Target:中心细胞
# Others: 其他细胞占比

p_res_pa <-
  res_pa_df %>%
  # dplyr::filter(Target == "Epithelial cells") %>%
  ggplot(aes(x = Others, y = cell_props, fill = Group)) +
  geom_boxplot(width = 0.35, outlier.colour = NA) +
  geom_jitter(width = 0.25) +
  facet_wrap(~ Target, scales = "free_y",ncol = 5)+
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# export plot
pdf(
  file.path(proximity_out_dir, "proximity analysis.pdf"),
  height = 10,
  width  = 20
)
p_res_pa %>% print()
dev.off()

table(res_pa_df$Target,res_pa_df$Group)


####循环OE_上皮 → EAOC_ Normal → EAOC_ Tumor-------
# cell_list <- c("OE_Normal Epi", "EAOC_Normal Epi", "EAOC_Tumor Epi")
# lapply(cell_list,function(cell){
#   p_res_pa <-
#     res_pa_df %>%
#     dplyr::filter(Target == cell) %>%
#     ggplot(aes(x = Others, y = cell_props, fill = Group)) +
#     geom_boxplot(width = 0.35, outlier.colour = NA) +
#     geom_jitter(width = 0.25) +
#     facet_wrap(~ Target, scales = "free_y",ncol = 5)+
#     ggpubr::theme_pubclean() +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
#   
#   # export plot
#   pdf(
#     file.path(proximity_out_dir, paste0(cell,"_proximity analysis.pdf")),
#     height = 5,
#     width  = 5
#   )
#   p_res_pa %>% print()
#   dev.off()
# })


####三种细胞合并展示------
res_pa_df$Target %>% unique()
res_pa_df$Others %>% unique()
plan <- c("Endothelial" = "#005995",
          "Pericyte" = "#96CEB4",
          "Stromal" = "#248888",
          "SMC" = "#73B66B",
          "Myeloid" = "#FF6D31",
          "Mast" = "#FDC4B6",
          "T_NK" = "#EF5285",
          "OE_Normal Epi" =  "#FBF579",
          "EAOC_Normal Epi" = "#FFC300",
          "EAOC_Tumor Epi" =  "#900475"             
)

p_res_pa <-
  res_pa_df %>%
  dplyr::filter(Target %in% cell_list) %>%
  ggplot(aes(x = Target, y = cell_props, fill = Others)) + 
  geom_col(position = "stack") +  # 堆叠柱状图
  #geom_boxplot(width = 0.35, outlier.colour = NA) +
  #geom_jitter(width = 0.25, alpha = 0.6) +
  #facet_wrap(~ Target, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = plan) +
  ggpubr::theme_pubclean() +
  labs(y = "Prop", fill = "celltype") +  
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1#, #vjust = 0.5
  ),
  legend.position = "right" )

ggsave(
  file.path(proximity_out_dir, "proximity_analysis_3type_cell_prop.pdf"),
  plot = p_res_pa,
  height = 5,
  width = 5
)



res_pa_df %>% head()


####三种细胞合并展示(桑基图)------
library(tidyverse)
library(Seurat)
library(ggalluvial)
library(ggthemes)
library(ggpubr)

cell_list <- c("OE_Normal Epi", "EAOC_Normal Epi", "EAOC_Tumor Epi")
celltype_order <- c(
  "Endothelial",
  "Pericyte", 
  "Stromal",
  "SMC",        
  "Myeloid",    
  "Mast",       
  "T_NK"#,  
  #"Epithelial"
  # "Normal Epi", 
  # "Tumor Epi"   
)

dat_props <- res_pa_df %>%
  dplyr::filter(Target %in% cell_list) %>%
  dplyr::filter(!Others %in% cell_list) %>% 
  mutate(freq = cell_props,
         celltype = Others,
  ) %>% 
  group_by(Target, celltype)

dat_props %>% head

dat_props_fixed <- dat_props %>%
  group_by(Target, celltype) %>%
  summarise(freq = sum(freq), .groups = "drop") %>% 
  mutate(Target = factor(Target,level = cell_list),
         celltype = factor(celltype,level = celltype_order))

# dat_props_immune <- all_data_OE_EAOC[[]] %>% 
#   group_by(group, subtype) %>%
#   summarise(count = n()) %>%
#   mutate(freq = count/sum(count))
# 
# dat_props_immune %>% head
# # A tibble: 6 × 4
# # Groups:   group [1]
# group annotation  count   freq
# <fct> <fct>       <int>  <dbl>
#   1 Ctrl  B cell       1437 0.153 
# 2 Ctrl  Plasma cell  1515 0.162 
# 3 Ctrl  T cell       4222 0.451 
# 4 Ctrl  NK cell       152 0.0162
# 5 Ctrl  Macrophage   1192 0.127 
# 6 Ctrl  Neutrophil    253 0.0270


# # A tibble: 6 × 4
# # Groups:   group [1]
# group annotation  count   freq
# <fct> <fct>       <int>  <dbl>
#   1 Ctrl  B cell       1437 0.153 
# 2 Ctrl  Plasma cell  1515 0.162 
# 3 Ctrl  T cell       4222 0.451 
# 4 Ctrl  NK cell       152 0.0162
# 5 Ctrl  Macrophage   1192 0.127 
# 6 Ctrl  Neutrophil    253 0.0270


p_props_flow <- 
  dat_props_fixed %>% 
  # mutate(id = 1:n()) %>%
  
  ggplot(aes(x = Target, y = freq, stratum = celltype, alluvium = celltype, fill = celltype, label = celltype)) +
  geom_stratum(alpha = 1, color = NA) +
  geom_flow(alpha = 0.5) +
  #geom_text(stat = "stratum") +
  theme_classic2() +
  theme(
    axis.line = element_blank(),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1#, #vjust = 0.5
    ),
    legend.position = "right")+
  scale_fill_manual(values = plan) +
  #scale_x_discrete(expand = c(0.01, 0.01)) +
  #scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(
    x = "group",
    y = "cell proportion"
  )

ggsave(file = file.path(proximity_out_dir,"3type_cell_prop_flow.pdf"), p_props_flow,width = 6,height = 6)
ggsave(file = file.path(proximity_out_dir,"3type_cell_prop_flow.png"), p_props_flow,width = 6,height = 6)




###2.2统计25μm范围内的细胞数量------
res_pa <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
    k = 750, # 'auto',
    searchtype = "radius",
    radius = 25, #  50um半径
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )

pa_groups <-
  all_data_meta_pos %>%
  group_by(Slide_ID) %>% # compute pa_groups per Slide_ID
  nest() %>%
  deframe() %>%
  imap(function(dat, Slide_ID){
    with(
      dat,
      str_c(Classification1, Group, Slide_ID, sep = '..')
    )
  }) %>%
  unlist(use.names = FALSE)

res_pa_df <- apply(res_pa, 1, function(x) {
  tapply(x, INDEX = pa_groups, FUN = median, na.rm = TRUE)  # 或用 mean
})

res_pa_df  <-
  res_pa_df %>%
  as.data.frame() %>%
  rownames_to_column('group') %>%
  pivot_longer(-group, names_to = "Others", values_to = "cell_counts") %>%  # 列名改为 cell_counts
  separate(group, into = c("Target", "Group", "project"), sep = "\\.\\.") 

res_pa_df$Target = factor(res_pa_df$Target,levels=group_subtype_level)
res_pa_df$Others = factor(res_pa_df$Others,levels=group_subtype_level)

p_res_pa <-
  res_pa_df %>%
  ggplot(aes(x = Others, y = cell_counts, fill = Group)) +  # 注意：y = cell_counts
  geom_boxplot(width = 0.35, outlier.colour = NA) +
  geom_jitter(width = 0.25, alpha = 0.6) +
  facet_wrap(~ Target, scales = "free_y", ncol = 5) +
  ggpubr::theme_pubclean() +
  labs(y = "Number") +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(
  file.path(proximity_out_dir, "proximity_analysis_cell_counts.pdf"),
  plot = p_res_pa,
  height = 10,
  width = 20
)


#### 循环绘制特定上皮细胞-----
cell_list <- c("OE_Normal Epi", "EAOC_Normal Epi", "EAOC_Tumor Epi")
lapply(cell_list, function(cell) {
  p_sub <- res_pa_df %>%
    dplyr::filter(Target == cell) %>%
    ggplot(aes(x = Others, y = cell_counts, fill = Group)) +
    geom_boxplot(width = 0.35, outlier.colour = NA) +
    geom_jitter(width = 0.25, alpha = 0.6) +
    facet_wrap(~ Target, scales = "free_y", ncol = 5) +
    ggpubr::theme_pubclean() +
    labs(y = "Number of cells within 25 μm") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  ggsave(
    file.path(proximity_out_dir, paste0(cell, "_proximity_analysis_counts.pdf")),
    plot = p_sub,
    height = 5,
    width = 5
  )
})

####三种细胞合并展示------
res_pa_df$Target %>% unique()
res_pa_df$Others %>% unique()
plan <- c("Endothelial" = "#005995",
          "Pericyte" = "#96CEB4",
          "Stromal" = "#248888",
          "SMC" = "#73B66B",
          "Myeloid" = "#FF6D31",
          "Mast" = "#FDC4B6",
          "T_NK" = "#EF5285",
          "OE_Normal Epi" =  "#FBF579",
          "EAOC_Normal Epi" = "#FFC300",
          "EAOC_Tumor Epi" =  "#900475"             
)

p_res_pa <-
  res_pa_df %>%
  dplyr::filter(Target %in% cell_list) %>%
  ggplot(aes(x = Target, y = cell_counts, fill = Others)) + 
  geom_col(position = "stack") +  # 堆叠柱状图
  #geom_boxplot(width = 0.35, outlier.colour = NA) +
  #geom_jitter(width = 0.25, alpha = 0.6) +
  #facet_wrap(~ Target, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = plan) +
  ggpubr::theme_pubclean() +
  labs(y = "Number", fill = "celltype") +  
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1#, #vjust = 0.5
  ),
  legend.position = "right" )

ggsave(
  file.path(proximity_out_dir, "proximity_analysis_3type_cell_counts.pdf"),
  plot = p_res_pa,
  height = 5,
  width = 5
)



####三种细胞合并展示(桑基图)------
library(tidyverse)
library(Seurat)
library(ggalluvial)
library(ggthemes)
library(ggpubr)

cell_list <- c("OE_Normal Epi", "EAOC_Normal Epi", "EAOC_Tumor Epi")
celltype_order <- c(
  "Endothelial",
  "Pericyte", 
  "Stromal",
  "SMC",        
  "Myeloid",    
  "Mast",       
  "T_NK"#,  
  #"Epithelial"
  # "Normal Epi", 
  # "Tumor Epi"   
)

dat_props <- res_pa_df %>%
  dplyr::filter(Target %in% cell_list) %>%
  dplyr::filter(!Others %in% cell_list) %>% 
  mutate(freq = cell_counts/sum(cell_counts),
         celltype = Others,
  ) %>% 
  group_by(Target, celltype) %>%
  summarise(cell_counts = sum(cell_counts), .groups = "drop") %>% 
  mutate(Target = factor(Target,level = cell_list),
         celltype = factor(celltype,level = celltype_order))

dat_props %>% head
# # A tibble: 6 × 4
# # Groups:   group [1]
# group annotation  count   freq
# <fct> <fct>       <int>  <dbl>
#   1 Ctrl  B cell       1437 0.153 
# 2 Ctrl  Plasma cell  1515 0.162 
# 3 Ctrl  T cell       4222 0.451 
# 4 Ctrl  NK cell       152 0.0162
# 5 Ctrl  Macrophage   1192 0.127 
# 6 Ctrl  Neutrophil    253 0.0270


p_props_flow <- 
  dat_props %>% 
  # mutate(id = 1:n()) %>%
  
  ggplot(aes(x = Target, y = cell_counts, stratum = celltype, alluvium = celltype, fill = celltype, label = celltype)) +
  geom_stratum(alpha = 1, color = NA) +
  geom_flow(alpha = 0.5) +
  #geom_text(stat = "stratum") +
  theme_classic2() +
  theme(
    axis.line = element_blank(),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1#, #vjust = 0.5
    ),
    legend.position = "right")+
  scale_fill_manual(values = plan) +
  #scale_x_discrete(expand = c(0.01, 0.01)) +
  #scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(
    x = "group",
    y = "cell count"
  )

ggsave(file = file.path(proximity_out_dir,"3type_cell_counts_flow.pdf"), p_props_flow,width = 6,height = 6)
ggsave(file = file.path(proximity_out_dir,"3type_cell_counts_flow.png"), p_props_flow,width = 6,height = 6)






#20260428analysis-----
#图2d、OE_上皮 → EAOC_ Normal → EAOC_ Tumor 三组周围25微米范围内细胞比例统计（邻近度）------
all_data_OE_EAOC <- readRDS("result0325/all_data_OE_EAOC.rds")
all_data_OE_EAOC@meta.data %>% head()
rownames(all_data_OE_EAOC@meta.data) %>% head()
all_data_OE_EAOC@meta.data$subtype %>% table()

group_subtype_level <- c("OE_Normal Epi","EAOC_Normal Epi","EAOC_Tumor Epi","Endothelial", "Pericyte", "Stromal", "SMC", "Myeloid", "Mast", "T_NK")

all_data_OE_EAOC@meta.data <- all_data_OE_EAOC@meta.data %>% 
  mutate(
    group_subtype = case_when(
      subtype %in% c("Normal Epi","Tumor Epi") ~ paste0(group, "_", subtype),
      TRUE ~ subtype  
    ),
    group_subtype = factor(group_subtype,levels=group_subtype_level)
  )

unique(all_data_OE_EAOC$group_subtype)

table(all_data_OE_EAOC$group, all_data_OE_EAOC$subtype)
table(all_data_OE_EAOC$group, all_data_OE_EAOC$group_subtype)


cells_boundary_ls <- readRDS(file.path("~/R_OUTPUT/data_for_analysis", "cells_boundary_ls.rds"))
cells_position_ls <- readRDS(file.path("~/R_OUTPUT/data_for_analysis", "cells_position_ls.rds"))

library(Seurat)
library(tidyverse)

source("spatial analysis.R")

proximity_out_dir <- "result0428/proximity_out_dir"
dir.create(proximity_out_dir)

###sample名称和slide名称需一致，slide名称是cells_position_ls各FOV的名字
## Proximity analysis -------
## 1. input data ------------
cells_position_ls %>% head() 

# all_data_cn_meta
all_data_meta_pos <-
  cells_position_ls %>%
  imap(function(cells_position, slide){ # browser()
    
    all_data_meta <- all_data_OE_EAOC[[]]  
    
    proj_meta_sub <-
      all_data_meta %>%  
      dplyr::select(group_subtype, sample, "group") %>% 
      rownames_to_column("cell_id") %>%
      dplyr::filter(sample == slide)
    
    cells_position[, 1:3] %>%
      inner_join(proj_meta_sub, by = "cell_id") %>%  
      mutate(group_subtype = factor(group_subtype,levels=group_subtype_level)) %>%
      dplyr::rename(Classification1 = group_subtype,
                    X = x_centroid,
                    Y = y_centroid,
                    Slide_ID = sample,
                    Group = group
      )
  }) %>%
  enframe() %>%
  unnest(cols = value) %>%
  dplyr::select(-name) %>%
  column_to_rownames("cell_id")

all_data_meta_pos %>% head()
all_data_meta_pos$Slide_ID %>% table()
all_data_meta_pos$Group %>% table()
all_data_meta_pos$Classification1 %>% table()

all_data_meta_pos$cell_id <- rownames(all_data_meta_pos)
all_data_meta_pos <- all_data_meta_pos %>%
  as.data.frame() %>%       
  remove_rownames()    

write.csv(all_data_meta_pos,file.path(proximity_out_dir,"all_data_meta_pos.csv"))

all_data_meta_pos <- read.csv("result0428/proximity analysis/all_data_meta_pos.csv")

# all_data_meta_pos_subset <- readRDS('data/all_data_meta_pos_subset.rds')
# all_data_meta_pos_subset %>% head

#        cell_id Classification1   Group           Slide_ID        X        Y
# 1 3_aaaacdgo-1     Fibroblasts   MUT    20231128__R01-001 2228.621 860.3357
# 2 3_aaaanalm-1   Myeloid cells   MUT    20231128__R01-001 2248.922 850.1516
# 3 3_aaaaonpl-1   Myeloid cells   MUT    20231128__R01-001 2270.447 856.8892
# 4 3_aaabciio-1         B cells   MUT    20231128__R01-001 2279.562 838.6568
# 5 3_aaabdkea-1     Fibroblasts   MUT    20231128__R01-001 2292.319 858.3909
# 6 3_aaacdgip-1         B cells   MUT    20231128__R01-001 2288.924 868.3226
# ......


## 2. proximity analysis -------
# 参数介绍
# x,y 指定细胞的坐标
# group.by 一般是切片或者样本编号，代表分别计算最近邻关系
# classification 一般是细胞的注释
# k 计算的最近邻的细胞数
# searchtype 邻近度使用radius

# radius是只查看当前细胞的多少um内的细胞
# mpp 坐标与um的转换比例，xenium是1，celldive是0.325
# returntype 可以每个细胞多少radius内的细胞计数count，也可以返回细胞频率freq
# cell_enough_ratio_cutoff 为了降低计算时间，k值不建议设置太高，
# 如果k值未能完全覆盖全部的邻近细胞，则cell_enough_ratio_cutoff会自动加k值，
# 直到默认90%的细胞均已经计算了全部邻近细胞



### 2.1 proximity analysis --------
radius_list <- c(10,25,30,60)

for (r in radius_list) {
  res_pa <-
    ProximutyAnalysis_NN_multi(
      all_data_meta_pos,
      query = all_data_meta_pos,
      x = "X",
      y = "Y",
      group.by = "Slide_ID",
      classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
      k = 750, # 'auto',
      searchtype = "radius",
      radius = r, #  50um半径
      mpp = 1, # !important, 1 for xenium
      returntype = "count",
      cell_enough_ratio_cutoff = 90
    )
  
  
  # return res_pa
  # rows: all celltypes
  # cols: query cell
  res_pa_colsum <- colSums(res_pa)
  res_pa <- sweep(res_pa, 2, res_pa_colsum, FUN = "/")
  if(any(res_pa_colsum == 0)){
    res_pa[, res_pa_colsum == 0] <- 0
  }
  
  # res_pa <- readRDS('data/res_pa.rds')
  
  
  # summarized by Group and celltype
  # pa_groups <-
  #   with(
  #     all_data_meta_pos,
  #     str_c(Classification1, Group, Slide_ID, sep = '..')
  #   )
  
  
  pa_groups <-
    all_data_meta_pos %>%
    group_by(Slide_ID) %>% # compute pa_groups per Slide_ID
    nest() %>%
    deframe() %>%
    imap(function(dat, Slide_ID){
      with(
        dat,
        str_c(Classification1, Group, Slide_ID, sep = '..')
      )
    }) %>%
    unlist(use.names = FALSE)
  
  # 分组内取细胞比例的中位数
  # apply对每行操作时会转置结果
  res_pa_df <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups, FUN = "median")})
  
  
  
  res_pa_df  <-
    res_pa_df %>%
    as.data.frame() %>%
    rownames_to_column('group') %>%
    pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
    separate(group, into = c("Target", "Group", "project"), sep = "\\.\\.") 
  
  
  res_pa_df$Target = factor(res_pa_df$Target,levels=group_subtype_level)
  res_pa_df$Others = factor(res_pa_df$Others,levels=group_subtype_level)
  # plot boxplot
  # Target:中心细胞
  # Others: 其他细胞占比
  
  saveRDS(res_pa_df,file = file.path(proximity_out_dir,paste0(r,"um_res_pa_df.rds")))
  
  p_res_pa <-
    res_pa_df %>%
    # dplyr::filter(Target == "Epithelial cells") %>%
    ggplot(aes(x = Others, y = cell_props, fill = Group)) +
    geom_boxplot(width = 0.35, outlier.colour = NA) +
    geom_jitter(width = 0.25) +
    facet_wrap(~ Target, scales = "free_y",ncol = 5)+
    ggpubr::theme_pubclean() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  # export plot
  pdf(
    file.path(proximity_out_dir,paste0(r, "um_proximity analysis.pdf")),
    height = 10,
    width  = 20
  )
  p_res_pa %>% print()
  dev.off()
  
  table(res_pa_df$Target,res_pa_df$Group)
  
  
  ####三种细胞合并展示------
  # res_pa_df$Target %>% unique()
  # res_pa_df$Others %>% unique()
  # plan <- c("Endothelial" = "#005995",
  #           "Pericyte" = "#96CEB4",
  #           "Stromal" = "#248888",
  #           "SMC" = "#73B66B",
  #           "Myeloid" = "#FF6D31",
  #           "Mast" = "#FDC4B6",
  #           "T_NK" = "#EF5285",
  #           "OE_Normal Epi" =  "#FBF579",
  #           "EAOC_Normal Epi" = "#FFC300",
  #           "EAOC_Tumor Epi" =  "#900475"             
  # )
  # 
  # p_res_pa <-
  #   res_pa_df %>%
  #   dplyr::filter(Target %in% cell_list) %>%
  #   ggplot(aes(x = Target, y = cell_props, fill = Others)) + 
  #   geom_col(position = "stack") +  # 堆叠柱状图
  #   #geom_boxplot(width = 0.35, outlier.colour = NA) +
  #   #geom_jitter(width = 0.25, alpha = 0.6) +
  #   #facet_wrap(~ Target, scales = "free_y", ncol = 5) +
  #   scale_fill_manual(values = plan) +
  #   ggpubr::theme_pubclean() +
  #   labs(y = "Prop", fill = "celltype") +  
  #   theme_classic2()+
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1#, #vjust = 0.5
  #   ),
  #   legend.position = "right" )
  # 
  # ggsave(
  #   file.path(proximity_out_dir,paste0(r, "um_proximity_analysis_3type_cell_prop.pdf")),
  #   plot = p_res_pa,
  #   height = 5,
  #   width = 5
  # )
  
  res_pa_df %>% head()
  
  ###三种细胞合并展示(桑基图)------
  library(tidyverse)
  library(Seurat)
  library(ggalluvial)
  library(ggthemes)
  library(ggpubr)
  
  cell_list <- c("OE_Normal Epi", "EAOC_Normal Epi", "EAOC_Tumor Epi")
  celltype_order <- c(
    "Endothelial",
    "Pericyte", 
    "Stromal",
    "SMC",        
    "Myeloid",    
    "Mast",       
    "T_NK"#,  
    #"Epithelial"
    # "Normal Epi", 
    # "Tumor Epi"   
  )
  
  dat_props <- res_pa_df %>%
    dplyr::filter(Target %in% cell_list) %>%
    dplyr::filter(!Others %in% cell_list) %>% 
    mutate(freq = cell_props,
           celltype = Others,
    ) %>% 
    group_by(Target, celltype)
  
  dat_props %>% head
  
  dat_props_fixed <- dat_props %>%
    group_by(Target, celltype) %>%
    summarise(freq = sum(freq), .groups = "drop") %>% 
    mutate(Target = factor(Target,level = cell_list),
           celltype = factor(celltype,level = celltype_order)) %>% 
    group_by(Target) %>%
    mutate(freq = freq / sum(freq)) %>%
    ungroup() 
  
  # dat_props_immune <- all_data_OE_EAOC[[]] %>% 
  #   group_by(group, subtype) %>%
  #   summarise(count = n()) %>%
  #   mutate(freq = count/sum(count))
  # 
  # dat_props_immune %>% head
  # # A tibble: 6 × 4
  # # Groups:   group [1]
  # group annotation  count   freq
  # <fct> <fct>       <int>  <dbl>
  #   1 Ctrl  B cell       1437 0.153 
  # 2 Ctrl  Plasma cell  1515 0.162 
  # 3 Ctrl  T cell       4222 0.451 
  # 4 Ctrl  NK cell       152 0.0162
  # 5 Ctrl  Macrophage   1192 0.127 
  # 6 Ctrl  Neutrophil    253 0.0270
  
  
  p_props_flow <- 
    dat_props_fixed %>% 
    # mutate(id = 1:n()) %>%
    
    ggplot(aes(x = Target, y = freq, stratum = celltype, alluvium = celltype, fill = celltype, label = celltype)) +
    geom_stratum(alpha = 1, color = NA) +
    geom_flow(alpha = 0.5) +
    #geom_text(stat = "stratum") +
    theme_classic2() +
    theme(
      axis.line = element_blank(),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1#, #vjust = 0.5
      ),
      legend.position = "right")+
    scale_fill_manual(values = plan) +
    #scale_x_discrete(expand = c(0.01, 0.01)) +
    #scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    labs(
      x = "group",
      y = "cell proportion"
    )
  
  ggsave(file = file.path(proximity_out_dir,paste0(r, "um_3type_cell_prop_flow.pdf")), 
         p_props_flow,
         width = 6,height = 6)
  ggsave(file = file.path(proximity_out_dir,paste0(r, "um_3type_cell_prop_flow.png")), 
         p_props_flow,
         width = 6,height = 6)

}

# 
###2.2统计25μm范围内的细胞数量------
# res_pa <-
#   ProximutyAnalysis_NN_multi(
#     all_data_meta_pos,
#     query = all_data_meta_pos,
#     x = "X",
#     y = "Y",
#     group.by = "Slide_ID",
#     classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
#     k = 750, # 'auto',
#     searchtype = "radius",
#     radius = 25, #  50um半径
#     mpp = 1, # !important, 1 for xenium
#     returntype = "count",
#     cell_enough_ratio_cutoff = 90
#   )
# 
# pa_groups <-
#   all_data_meta_pos %>%
#   group_by(Slide_ID) %>% # compute pa_groups per Slide_ID
#   nest() %>%
#   deframe() %>%
#   imap(function(dat, Slide_ID){
#     with(
#       dat,
#       str_c(Classification1, Group, Slide_ID, sep = '..')
#     )
#   }) %>%
#   unlist(use.names = FALSE)
# 
# res_pa_df <- apply(res_pa, 1, function(x) {
#   tapply(x, INDEX = pa_groups, FUN = median, na.rm = TRUE)  # 或用 mean
# })
# 
# res_pa_df  <-
#   res_pa_df %>%
#   as.data.frame() %>%
#   rownames_to_column('group') %>%
#   pivot_longer(-group, names_to = "Others", values_to = "cell_counts") %>%  # 列名改为 cell_counts
#   separate(group, into = c("Target", "Group", "project"), sep = "\\.\\.") 
# 
# res_pa_df$Target = factor(res_pa_df$Target,levels=group_subtype_level)
# res_pa_df$Others = factor(res_pa_df$Others,levels=group_subtype_level)
# 
# p_res_pa <-
#   res_pa_df %>%
#   ggplot(aes(x = Others, y = cell_counts, fill = Group)) +  # 注意：y = cell_counts
#   geom_boxplot(width = 0.35, outlier.colour = NA) +
#   geom_jitter(width = 0.25, alpha = 0.6) +
#   facet_wrap(~ Target, scales = "free_y", ncol = 5) +
#   ggpubr::theme_pubclean() +
#   labs(y = "Number") +  
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# 
# ggsave(
#   file.path(proximity_out_dir, "proximity_analysis_cell_counts.pdf"),
#   plot = p_res_pa,
#   height = 10,
#   width = 20
# )
# 
# 
#### 循环绘制特定上皮细胞-----
# cell_list <- c("OE_Normal Epi", "EAOC_Normal Epi", "EAOC_Tumor Epi")
# lapply(cell_list, function(cell) {
#   p_sub <- res_pa_df %>%
#     dplyr::filter(Target == cell) %>%
#     ggplot(aes(x = Others, y = cell_counts, fill = Group)) +
#     geom_boxplot(width = 0.35, outlier.colour = NA) +
#     geom_jitter(width = 0.25, alpha = 0.6) +
#     facet_wrap(~ Target, scales = "free_y", ncol = 5) +
#     ggpubr::theme_pubclean() +
#     labs(y = "Number of cells within 25 μm") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
#   
#   ggsave(
#     file.path(proximity_out_dir, paste0(cell, "_proximity_analysis_counts.pdf")),
#     plot = p_sub,
#     height = 5,
#     width = 5
#   )
# })
# 
####三种细胞合并展示------
# res_pa_df$Target %>% unique()
# res_pa_df$Others %>% unique()
# plan <- c("Endothelial" = "#005995",
#           "Pericyte" = "#96CEB4",
#           "Stromal" = "#248888",
#           "SMC" = "#73B66B",
#           "Myeloid" = "#FF6D31",
#           "Mast" = "#FDC4B6",
#           "T_NK" = "#EF5285",
#           "OE_Normal Epi" =  "#FBF579",
#           "EAOC_Normal Epi" = "#FFC300",
#           "EAOC_Tumor Epi" =  "#900475"             
# )
# 
# p_res_pa <-
#   res_pa_df %>%
#   dplyr::filter(Target %in% cell_list) %>%
#   ggplot(aes(x = Target, y = cell_counts, fill = Others)) + 
#   geom_col(position = "stack") +  # 堆叠柱状图
#   #geom_boxplot(width = 0.35, outlier.colour = NA) +
#   #geom_jitter(width = 0.25, alpha = 0.6) +
#   #facet_wrap(~ Target, scales = "free_y", ncol = 5) +
#   scale_fill_manual(values = plan) +
#   ggpubr::theme_pubclean() +
#   labs(y = "Number", fill = "celltype") +  
#   theme_classic2()+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1#, #vjust = 0.5
#   ),
#   legend.position = "right" )
# 
# ggsave(
#   file.path(proximity_out_dir, "proximity_analysis_3type_cell_counts.pdf"),
#   plot = p_res_pa,
#   height = 5,
#   width = 5
# )
# 
# 
# 
####三种细胞合并展示(桑基图)------
# library(tidyverse)
# library(Seurat)
# library(ggalluvial)
# library(ggthemes)
# library(ggpubr)
# 
# cell_list <- c("OE_Normal Epi", "EAOC_Normal Epi", "EAOC_Tumor Epi")
# celltype_order <- c(
#   "Endothelial",
#   "Pericyte", 
#   "Stromal",
#   "SMC",        
#   "Myeloid",    
#   "Mast",       
#   "T_NK"#,  
#   #"Epithelial"
#   # "Normal Epi", 
#   # "Tumor Epi"   
# )
# 
# dat_props <- res_pa_df %>%
#   dplyr::filter(Target %in% cell_list) %>%
#   dplyr::filter(!Others %in% cell_list) %>% 
#   mutate(freq = cell_counts/sum(cell_counts),
#          celltype = Others,
#   ) %>% 
#   group_by(Target, celltype) %>%
#   summarise(cell_counts = sum(cell_counts), .groups = "drop") %>% 
#   mutate(Target = factor(Target,level = cell_list),
#          celltype = factor(celltype,level = celltype_order))
# 
# dat_props %>% head
# # # A tibble: 6 × 4
# # # Groups:   group [1]
# # group annotation  count   freq
# # <fct> <fct>       <int>  <dbl>
# #   1 Ctrl  B cell       1437 0.153 
# # 2 Ctrl  Plasma cell  1515 0.162 
# # 3 Ctrl  T cell       4222 0.451 
# # 4 Ctrl  NK cell       152 0.0162
# # 5 Ctrl  Macrophage   1192 0.127 
# # 6 Ctrl  Neutrophil    253 0.0270
# 
# 
# p_props_flow <- 
#   dat_props %>% 
#   # mutate(id = 1:n()) %>%
#   
#   ggplot(aes(x = Target, y = cell_counts, stratum = celltype, alluvium = celltype, fill = celltype, label = celltype)) +
#   geom_stratum(alpha = 1, color = NA) +
#   geom_flow(alpha = 0.5) +
#   #geom_text(stat = "stratum") +
#   theme_classic2() +
#   theme(
#     axis.line = element_blank(),
#     axis.text = element_text(size = 12),
#     axis.text.x = element_text(angle = 45, hjust = 1#, #vjust = 0.5
#     ),
#     legend.position = "right")+
#   scale_fill_manual(values = plan) +
#   #scale_x_discrete(expand = c(0.01, 0.01)) +
#   #scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
#   labs(
#     x = "group",
#     y = "cell count"
#   )
# 
# ggsave(file = file.path(proximity_out_dir,"3type_cell_counts_flow.pdf"), p_props_flow,width = 6,height = 6)
# ggsave(file = file.path(proximity_out_dir,"3type_cell_counts_flow.png"), p_props_flow,width = 6,height = 6)
# 
# 


