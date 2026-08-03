library(tidyverse)
library(Seurat)


all_data_sc2 <- readRDS('all_data_sc2.rds')

out_dir <- "20260605 dat analysis"
dir.create(out_dir)


# 1. dotplot -------------
marker_ls <- c(
  "AT1" = "Ager",
  "AT2" = "Sftpa1",
  "AT2" = "Sftpd",
  
  "Ciliated cell" = "Cdhr3",
  "Ciliated cell" = "Cetn2",
  "Ciliated cell" = "Dnah11",
  # "Club cell" = "Scgb1a1",
  "Club cell" = "Scgb3a1",
  "Club cell" = "Scgb3a2",
  
  "Endothelial" = "Pecam1",
  "Endothelial" = "Vwf",
  
  "Fibroblasts" = "Pdgfra",
  "Fibroblasts" = "Col1a1",
  "Fibroblasts" = "Dcn",

  
  "Smooth muscle cells" = "Tagln",
  "Smooth muscle cells" = "Myh11",
  
  "B/Plasma cell" = "Cd19",
  "B/Plasma cell" = "Cd79a", 
  "B/Plasma cell" = "Ms4a1",
  
  "B/Plasma cell" = "Mzb1",
  "B/Plasma cell" = "Jchain",
  
  "NK/T cell" = "Cd3e", 
  "NK/T cell" = "Cd4",
  "NK/T cell" = "Cd8a",
  
  "NK/T cell" = "Ncam1",
  "NK/T cell" = "Gnly",
  "NK/T cell" = "Nkg7",
  
  "DC" = "Lyz1",
  "DC" = "Il3ra",
  "DC" = "Clec9a",
  "DC" = "Irf4",
  # "DC" = "Xcr1",
  # "DC" = "Itgax",
  # "DC" = "Bst2",
  # "DC" = "Siglech",
  
  "Monocyte" = "Cd14",
  "Monocyte" = "Fcgr3",
  "Monocyte" = "Vcan",
  # "Monocyte" = "Fcna",
  "Monocyte" = "Csf1r",
  # "Monocyte" = "Ly6c1", 
  
  
  "Alveolar Macrophage" = "Apoc1",
  "Alveolar Macrophage" = "Gpnmb",
  # "Alveolar Macrophage" = "Inhba",
  # "Alveolar Macrophage" = "Fcgr1", ### CD64
  "Alveolar Macrophage" = "Siglecf",

  
  "Interstitial Macrophages" = "Adgre1", ### F4/80
  "Interstitial Macrophages" = "Mrc1", 
  
  "Neutrophil" = "Csf3r",
  "Neutrophil" = "S100a8",
  "Neutrophil" = "S100a9"
  # "Neutrophil"  = "Ly6g", ### 
  # "Neutrophil"  = "Itgam" # CD11b

)

all_data_sc2$annotation <- factor(all_data_sc2$annotation, levels = anno_levels)
anno_levels <- c(
  'AT1',
  'AT2',
  'Ciliated/Club cell',
  'Endothelial',
  'Fibroblasts',
  'Smooth muscle cells',
  'B/Plasma',
  'NK/T',
  'DC',
  'Monocyte',
  'Alveolar mac.',
  'Interstitial mac.',
  'Neutrophil'
)

# export
pdf(file.path(out_dir, str_glue("1. dotplot of cell marker.pdf")), width = 25, height = 6)
  DotPlot(all_data_sc2, group.by = "annotation", features = marker_ls) + 
    RotatedAxis() + 
    theme(strip.clip = "off") +
    scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(Inf, "RdBu")))
dev.off()


# 2. Ro/e ---------------
# library(tidyverse)

# all cells
roe_out_dir <- out_dir

anno_column <- "annotation"
grp_column  <- "project"

all_data_meta_roe <-
  all_data_sc2[[]] %>%
  dplyr::filter(.data[[anno_column]]  != "doublet") %>%
  mutate({{anno_column}} := fct_drop(  .data[[anno_column]] ))


# Cell Preference
library(ComplexHeatmap)
# 生成细胞数目分布频数统计表observe.data
observe.data <- unclass(table(all_data_meta_roe[[anno_column]], all_data_meta_roe[[grp_column]]))
stopifnot(is.matrix(observe.data))

expected <-
  function(x){
    if(!is.matrix(x)){stop("Must be a matrix")}
    rtot <- margin.table(x, 1)
    ctot <- margin.table(x, 2)
    tot  <- margin.table(x)
    outer(rtot, ctot, "*")/tot
  }
expected.data <- expected(observe.data) # 生成期望频数表
plot.data <- observe.data/expected.data # 计算实际细胞数目与期望细胞数目的比值
# 绘图数据的列顺序
# col.order <- ggplot2:::`%||%`(
#   levels(all_data_meta_roe[[grp_column]]),
#   unique(all_data_meta_roe[[grp_column]])
# ) 
col.order <- c(
  "P3_Normoxia",
  "P7_Normoxia",
  "P14_Normoxia",
  
  "P3_Hyperoxia",  
  "P7_Hyperoxia",  
  "P14_Hyperoxia"
)
plot.data <- plot.data[, col.order] # 对绘图数据列进行排序
plot.data <- plot.data[order(apply(plot.data, 1, which.max), decreasing = T), ] # 对绘图数据行进行排序
# 如果想深色色块在主对角线使用which.min，在副对角线则使用which.max

# 热图颜色设定
roe_cols <- structure(
  c("#FEE6CE", "#FDC08C", "#F5904B", "#E6550D"),
  names = c("1", "1.5", "2", ">3")
)

cell_fun = function(j, i, x, y, width, height, fill) { # browser()
  
  # print(c(j, i, x, y, width, height, fill))
  
  if(plot.data[i,j] < as.numeric(names(roe_cols)[1])) {
    grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = roe_cols[1], fill = roe_cols[1]))
    grid::grid.text("\u00B1", x = x, y = y)
  } else if(plot.data[i,j] < as.numeric(names(roe_cols)[2])) {
    grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = roe_cols[2], fill = roe_cols[2]))
    grid::grid.text("+", x = x, y = y)
  } else if(plot.data[i,j] < as.numeric(names(roe_cols)[3])) {
    grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = roe_cols[3], fill = roe_cols[3]))
    grid::grid.text("++", x = x, y = y)
  } else {
    grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = roe_cols[4], fill = roe_cols[4]))
    grid::grid.text("+++", x = x, y = y)
  }
  
}

# 生成热图
roi_hm <-
  ComplexHeatmap::Heatmap(
    plot.data, 
    cell_fun = cell_fun,
    cluster_rows = F, 
    cluster_columns = F,
    width = unit(5, "cm"), 
    height = unit(15, "cm"), 
    show_heatmap_legend = F
  )

# stat comparasion
chi2.T <- sum((observe.data - expected.data)^2/expected.data)
p_chisq <- pchisq(q = chi2.T, df = (nrow(observe.data)-1)*(ncol(observe.data)-1), lower.tail = FALSE)
print(p_chisq)

# 生成图例
lgd <-
  ComplexHeatmap::Legend(
    labels = names(roe_cols), 
    title = expression(R[o/e]), 
    legend_gp = grid::gpar(fill = roe_cols)
  )

pdf(
  file.path(roe_out_dir, "2. cell props boxplot - major annotation - Roe plot.pdf"), 
  width = 8, 
  height = 12
)
ComplexHeatmap::draw(roi_hm, annotation_legend_list = lgd)
dev.off()


# 2. Ro/e -----------------
ROIE <- function(crosstab){
  ## Calculate the Ro/e value from the given crosstab
  ##
  ## Args:
  #' @crosstab: the contingency table of given distribution
  ##
  ## Return:
  ## The Ro/e matrix 
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  return(roie)
}

divMatrix <- function(m1, m2){
  ## Divide each element in turn in two same dimension matrixes
  ##
  ## Args:
  #' @m1: the first matrix
  #' @m2: the second matrix
  ##
  ## Returns:
  ## a matrix with the same dimension, row names and column names as m1. 
  ## result[i,j] = m1[i,j] / m2[i,j]
  dim_m1 <- dim(m1)
  dim_m2 <- dim(m2)
  if( sum(dim_m1 == dim_m2) == 2 ){
    div.result <- matrix( rep(0,dim_m1[1] * dim_m1[2]) , nrow = dim_m1[1] )
    row.names(div.result) <- row.names(m1)
    colnames(div.result) <- colnames(m1)
    for(i in 1:dim_m1[1]){
      for(j in 1:dim_m1[2]){
        div.result[i,j] <- m1[i,j] / m2[i,j]
      }
    }   
    return(div.result)
  }
  else{
    warning("The dimensions of m1 and m2 are different")
  }
}

roe_out_dir <- out_dir

df = all_data_sc2[[]]
summary <- table(df[,c('annotation','project')])
roe <- as.data.frame(ROIE(summary))

long_df <- roe %>%
  rownames_to_column('cell_type') %>%
  pivot_longer(
    -cell_type,
    names_to = "group",
    values_to = "Ro/e"
  )
write.csv(roe, file = file.path(roe_out_dir,"Roe.xls"),sep = "\t",quote=F)
# roe = read.csv('/lvdata/wzb/scRNA/FW2024-213/fig/fig1/fig1c/Roe.xls',row.names = 1)
# roe$celltype = rownames(roe)




# 依据Roe值进行定义
long_df <-
  long_df %>% 
  mutate(class = case_when(
    `Ro/e` >= 1.05 ~ "Enrichment",
    `Ro/e` <= 0.95 ~ "Depletion",
    TRUE ~ "No Change",
  )) %>%
  mutate(class = factor(class, levels = c("Depletion","No Change","Enrichment")))


row.order <- c(
  "P3_Normoxia",
  "P7_Normoxia",
  "P14_Normoxia",
  
  "P3_Hyperoxia",  
  "P7_Hyperoxia",  
  "P14_Hyperoxia"
)
if(!is.null(row.order)){
  long_df$group <- factor(long_df$group, levels = row.order)
}

# 创建预设的样式
p_roe <-
  long_df %>%
  ggplot(aes(x = group, y = cell_type))+
  geom_point(aes(color = class, size = `Ro/e`))+
  scale_color_manual(values = c(
    "Depletion" = '#4a6fe3',
    "No Change" = 'grey',
    "Enrichment" = '#8e063b'
  ))+
  scale_size(
    range  = c(2, 8),
    limits = long_df$`Ro/e` %>% range() %>% {c(floor(.[1]), ceiling(.[2]))},
    breaks = long_df$`Ro/e` %>% range() %>% {seq(floor(.[1]), ceiling(.[2]), length.out = 5)}
  ) +
  coord_flip() +
  theme_bw() +
  theme(
    plot.margin = margin(.5, .5, .5, .5,'cm'),
    panel.grid.major=element_line(colour=NA),
    axis.text.x=element_text(hjust=1,vjust=1,angle = 60,size=10)
  ) +
  labs(x = "Group",y = "Celltype",color="",size="Ro/e")

ggsave(p_roe, filename = file.path(roe_out_dir, "2. Roe bubbleplot.png"), width = 24, height = 14, units = "cm")
ggsave(p_roe, filename = file.path(roe_out_dir, "2. Roe bubbleplot.pdf"), width = 24, height = 14, units = "cm")

p_roi_heatmap <-
  pheatmap::pheatmap(
    roe[, row.order], 
    color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "Reds"))(10),
    cluster_rows = F, 
    cluster_cols = F, 
    display_numbers = TRUE,
    number_color = "black",
    # cellwidth = 30,
    # cellheight = 16, 
    fontsize = 10,                    # fontsize字号大小
    border_color = '#ffffff',
    angle_col='45'
  )


pdf(file.path(roe_out_dir, "2. Roe heatmap.pdf"), width = ncol(roe)*1, height = nrow(roe)*0.5)
p_roi_heatmap$gtable %>% cowplot::plot_grid(scale = 0.9)
dev.off()

# png(file.path(roe_out_dir, "2. Roe heatmap.png"), width = ncol(roe)*1, height = nrow(roe)*0.5)
# p_roi_heatmap$gtable %>% cowplot::plot_grid(scale = 0.9)
# dev.off()


# 3. hallmark ----------------
h_geneset <- msigdbr::msigdbr(category = "H", species = "Mus musculus")

h_geneset$gs_name %>%
  unique() %>% 
  str_to_lower() %>%
  str_remove("^hallmark_")

h_geneset_ls <-
  h_geneset[c('gs_name', 'gene_symbol')] %>%
  with(split(gene_symbol, gs_name)) %>%
  {set_names(
    .,
    names(.) %>% 
      str_to_lower() %>%
      str_remove("^hallmark_")
  )}

all_data_sc2_score <-
  all_data_sc2 %>%
  JoinLayers() %>%
  AddModuleScore(
    features = h_geneset_ls,
  )

all_data_sc2_score %>% head

# msigdb anno
h_all_data_sc2_score <- 
  all_data_sc2_score[[]] %>%
  dplyr::select(project, Cluster1:Cluster50) %>%
  group_by(project) %>%
  summarise(across(everything(), mean))

p_hallmark_1 <-
  pheatmap::pheatmap(
    h_all_data_sc2_score %>% column_to_rownames("project") %>% t,
    scale = 'row',
    labels_row = names(h_geneset_ls),
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(Inf, "RdBu")))(50),
    border_color = NA,
    silent = TRUE
  )

# export 
pdf(
  file.path(out_dir, "3. Hallmark pathway score of all samples - all cells.pdf"), 
  height = 12, 
  width  = 5
)
p_hallmark_1$gtable %>% cowplot::plot_grid() %>% print()
dev.off()







