
setwd('~/')

library(Seurat)
library(tidyverse)
library(patchwork)
library(babelgene)
library(openxlsx)
library(VennDiagram)
library(grid)

common_markers <- read.table('DEGs_by_group/common_markers.txt', header = F) %>% unlist
PCDpath <- openxlsx::read.xlsx('DEGs_by_group/18PCDgenes.xlsx')
PCDpath_mouse <- PCDpath


for(col in colnames(PCDpath)){
  human_genes <- as.character(PCDpath[[col]])
  human_genes <- human_genes[!is.na(human_genes)]
  human_genes <- human_genes[human_genes != ""]
  human_genes <- trimws(human_genes) # 去首尾空格
  if(length(human_genes)==0) next
  
  # 转换
  orth <- orthologs(genes = human_genes, species = "mouse")
  mouse_genes <- orth$symbol
  
  # 补齐长度
  final <- rep(NA, nrow(PCDpath))
  if(length(mouse_genes)>0) {
    final[seq_along(mouse_genes)] <- mouse_genes
  }
  
  PCDpath_mouse[[col]] <- final
}


PCDgenes <- unlist(PCDpath_mouse)[!(unlist(PCDpath_mouse) %>% is.na())]
#common_markers[common_markers %in% str_to_title(PCDgenes)]
heart_markers <- read.csv('DEGs_by_group/heart_DEGs.csv') %>% column_to_rownames('X')
liver_markers <- read.csv('DEGs_by_group/liver_DEGs.csv') %>% column_to_rownames('X')
gut_markers <- read.csv('DEGs_by_group/gut_DEGs.csv') %>% column_to_rownames('X')
lung_markers <- read.csv('DEGs_by_group/lung_DEGs.csv') %>% column_to_rownames('X')
kidney_markers <- read.csv('DEGs_by_group/kidney_DEGs.csv') %>% column_to_rownames('X')
rownames(heart_markers)[rownames(heart_markers) %in% PCDgenes]
rownames(liver_markers)[rownames(liver_markers) %in% PCDgenes]
rownames(gut_markers)[rownames(gut_markers) %in% PCDgenes]
rownames(lung_markers)[rownames(lung_markers) %in% PCDgenes]
rownames(kidney_markers)[rownames(kidney_markers) %in% PCDgenes]


heart_obj <- readRDS("heart_annotation/heart.rds")
liver_obj <- readRDS("liver_annotation/liver.rds")
gut_obj   <- readRDS("gut_annotation/gut.rds")
lung_obj  <- readRDS("lung_annotation/lung.rds")
kidney_obj<- readRDS("kidney_annotation/kidney.rds")

seurat_list <- list(
  heart = heart_obj,
  liver = liver_obj,
  gut   = gut_obj,
  lung  = lung_obj,
  kidney= kidney_obj
)
table(heart_obj$annotation)
table(heart_obj$subtype)
table(heart_obj$group)
table(liver_obj$annotation)
table(liver_obj$subtype_re)
table(liver_obj$group)
table(gut_obj$annotation)
table(gut_obj$subtype)
table(gut_obj$group)
colnames(gut_obj@meta.data)
table(lung_obj$annotation)
table(lung_obj$subtype)
table(lung_obj$group)
table(kidney_obj$annotation)
table(kidney_obj$subtype)
table(kidney_obj$group)




tissues <- c("heart", "liver", "gut", "lung", "kidney")

#tissues <- c("liver")

# 2. 创建输出文件夹
dir.create("PCD_analysis_results", showWarnings = F)
dir.create("PCD_analysis_results/1_PCDgenes_degs_venn_plots", showWarnings = F)
dir.create("PCD_analysis_results/2_score_violin_plots", showWarnings = F)
dir.create("PCD_analysis_results/3_median_score_DEGs", showWarnings = F)

# 3. 核心循环：5个组织独立分析
for(tissue in tissues){  browser()
  #取出对应组织的Seurat对象
  seurat_obj <- seurat_list[[tissue]]
  
  #取出对应组织的差异基因 markers
  markers_var <- paste0(tissue, "_markers")
  tissue_markers <- get(markers_var)
  
  #需求1：【每个组织独立取交集】markers ∩ PCDgenes + 韦恩图
  tissue_markers_genes <- rownames(tissue_markers)
  tissue_PCD_genes <- intersect(tissue_markers_genes, PCDgenes)
  markers_genes <- setdiff(tissue_markers_genes, PCDgenes)
  PCD_genes <- setdiff(PCDgenes, tissue_markers_genes)
  
  cat(tissue, "组织筛选到PCD交集基因：", length(tissue_PCD_genes), "个\n")
  print(tissue_PCD_genes)
  
  #保存该组织的PCD交集基因
  write.table(tissue_PCD_genes,
              paste0("PCD_analysis_results/1_PCDgenes_degs_venn_plots/", tissue, "_PCD_intersect_genes.txt"),
              row.names = F, col.names = F, quote = F)
  
  #保存韦恩图基因列表
  venn_gene_list <- data.frame(
    gene = c(markers_genes, PCD_genes, tissue_PCD_genes),
    type = c(
      rep("markers_genes", length(markers_genes)),
      rep("PCD_genes", length(PCD_genes)),
      rep("Intersection", length(tissue_PCD_genes))
    )
  )
  write.csv(venn_gene_list,
            paste0("PCD_analysis_results/1_PCDgenes_degs_venn_plots/", tissue, "_PCD_markers_venn_list.csv"),
            row.names = F)
  
  #画韦恩图
  pdf(paste0("PCD_analysis_results/1_PCDgenes_degs_venn_plots/", tissue, "_PCD_markers_venn.pdf"),
      width = 6, height = 6)
  plot.new()
  venn.plot <- venn.diagram(
    x = list(
      Tissue_Markers = tissue_markers_genes,
      PCD_Genes = PCDgenes
    ),
    category = c(paste0(tissue, "_markers"), "PCD_genes"),
    fill = c("#E41A1C", "#377EB8"),
    alpha = 0.5,
    cat.cex = 1.5,
    cex = 2,
    filename = NULL,
    cat.pos = c(350, 10)
  )
  grid.draw(venn.plot)
  dev.off()
  
  #如果没有基因，跳过
  if(length(tissue_PCD_genes) == 0){
    cat(tissue, "无PCD交集基因，跳过\n\n")
    next
  }
  
  
  # 需求2：只取 HS_组织 样本 → AddModuleScore → 小提琴图
  if(tissue == "gut"){
    #肠道：4组合并
    seurat_hs <- subset(seurat_obj, group %in% c("HS0_gut", "HS6_gut", "HS12_gut", "HS_gut"))
  } else {
    #其他组织：正常分组
    hs_group <- paste0("HS_", tissue)
    seurat_hs <- subset(seurat_obj, group == hs_group)
  }
  
  #打分
  seurat_hs <- AddModuleScore(
    object = seurat_hs,
    features = list(tissue_PCD_genes),
    name = "PCD_score"
  )
  
  if(tissue == "liver"){
    celltype_col <- "subtype_re"
  } else {
    celltype_col <- "subtype"
  }
  
  #小提琴图
  p <- VlnPlot(
    seurat_hs,
    features = "PCD_score1",
    group.by = celltype_col
  ) + ggtitle(paste0(tissue, " - HS_group_PCD_genes_score"))
  
  ggsave(paste0("PCD_analysis_results/2_score_violin_plots/", tissue, "_PCD_score_violin.png"),
         p, width = 16, height = 6, dpi = 300)
  ggsave(paste0("PCD_analysis_results/2_score_violin_plots/", tissue, "_PCD_score_violin.pdf"),
         p, width = 16, height = 6)
  
  
  # 需求3：每个细胞类型独立分析 → 组织内所有细胞类型DEGs取交集 + 韦恩图
  all_celltype_degs <- list()
  
  celltypes <- as.character(seurat_hs[[celltype_col]][,1])
  
  for (ct in unique(celltypes)) {  
    
    cat("  ↳ 分析细胞类型：", ct, "\n")
    
    ct_cells <- colnames(seurat_hs)[celltypes == ct]
    seurat_ct <- subset(seurat_hs, cells = ct_cells)
    
    ct_score <- seurat_ct$PCD_score1
    ct_median <- median(ct_score)
    
    seurat_ct$cell_death_status <- ifelse(ct_score >= ct_median, "death_related", "non_death")
    
    Idents(seurat_ct) <- "cell_death_status"
    degs_ct <- FindMarkers(
      object = seurat_ct,
      ident.1 = "death_related",
      ident.2 = "non_death",
      test.use = "wilcox",
      min.pct = 0.1,
      logfc.threshold = 0.25,
      only.pos = F
    )
    
    write.table(
      degs_ct,
      paste0("PCD_analysis_results/3_median_score_DEGs/", tissue, "_", gsub("/","_",ct), "_DEGs.txt"),
      sep = "\t", quote = F, row.names = T
    )
    
    #保存显著基因（单独也存一份）
    degs_ct_sig <- rownames(degs_ct)[degs_ct$p_val_adj < 0.05]
    write.table(
      degs_ct_sig,
      paste0("PCD_analysis_results/3_median_score_DEGs/", tissue, "_", gsub("/","_",ct), "_sig_DEGs.txt"),
      row.names = F, col.names = F, quote = F
    )
    
    if (length(degs_ct_sig) == 0) {
      cat("    - 细胞类型", ct, "无显著DEGs\n")
      next
    }
    
    all_celltype_degs[[ct]] <- degs_ct_sig
  }
  
  #取交集（修复空列表报错）
  if(length(all_celltype_degs) > 0){
    final_intersect <- Reduce(intersect, all_celltype_degs)
  } else {
    final_intersect <- c()
  }
  
  #保存组织最终交集基因
  write.table(
    final_intersect,
    paste0("PCD_analysis_results/3_median_score_DEGs/", tissue, "_all_celltypes_DEG_intersection.txt"),
    row.names = F, col.names = F, quote = F
  )
  
  cat(tissue, "组织最终交集基因数量：", length(final_intersect), "✅\n\n")
}


# 单独重跑一下小提琴图
tissues <- c("heart", "liver", "gut", "lung", "kidney")

#tissues <- c("liver")

for(tissue in tissues){
  
  # 读取之前保存的交集基因
  tissue_PCD_genes <- read.table(
    paste0("PCD_analysis_results/1_PCDgenes_degs_venn_plots/", tissue, "_PCD_intersect_genes.txt"),
    header=F, sep="\t", stringsAsFactors=F
  )[,1]
  
  if(length(tissue_PCD_genes) == 0){
    cat(tissue, "无PCD交集基因，跳过\n")
    next
  }
  
  # 取出Seurat对象
  seurat_obj <- seurat_list[[tissue]]
  
  # 提取HS组
  if(tissue == "gut"){
    seurat_hs <- subset(seurat_obj, group %in% c("HS0_gut", "HS6_gut", "HS12_gut", "HS_gut"))
  } else {
    hs_group <- paste0("HS_", tissue)
    seurat_hs <- subset(seurat_obj, group == hs_group)
  }
  
  # 打分
  seurat_hs <- AddModuleScore(
    object = seurat_hs,
    features = list(tissue_PCD_genes),
    name = "PCD_score"
  )
  
  # 细胞类型列
  celltype_col <- ifelse(tissue == "liver", "subtype_re", "subtype")
  
  # 小提琴图（无点）
  p <- VlnPlot(
    seurat_hs,
    features = "PCD_score1",
    group.by = celltype_col,
    pt.size = 0
  ) + ggtitle(paste0(tissue, " - HS_group_PCD_genes_score"))
  
  ggsave(paste0("PCD_analysis_results/2_score_violin_plots/", tissue, "_PCD_score_violin.png"),
         p, width = 16, height = 6, dpi = 300)
  ggsave(paste0("PCD_analysis_results/2_score_violin_plots/", tissue, "_PCD_score_violin.pdf"),
         p, width = 16, height = 6)
}




# DSS ------------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3) # 100Gb


seurat_ls <- 
  seurat_list %>% 
  imap(function(seurat_obj, tissue){ # browser()
    # 取出对应组织的差异基因 markers
    markers_var <- paste0(tissue, "_markers")
    tissue_markers <- get(markers_var)
    
    # markers ∩ PCDgenes 取交集
    tissue_markers_genes <- rownames(tissue_markers)
    tissue_PCD_genes <- intersect(tissue_markers_genes, PCDgenes)
    # markers_genes <- setdiff(tissue_markers_genes, PCDgenes)
    # PCD_genes <- setdiff(PCDgenes, tissue_markers_genes)
    
    # 只取 HS_组织 样本 → AddModuleScore
    if(tissue == "gut"){
      #肠道：4组合并
      seurat_hs <- subset(seurat_obj, group %in% c("HS0_gut", "HS6_gut", "HS12_gut", "HS_gut"))
    } else {
      #其他组织：正常分组
      hs_group <- paste0("HS_", tissue)
      seurat_hs <- subset(seurat_obj, group == hs_group)
    }
    
    #打分
    seurat_hs <- AddModuleScore(
      object = seurat_hs,
      features = list(tissue_PCD_genes),
      name = "PCD_score"
    )
    
    meta_DSS <- 
      seurat_hs@meta.data %>% 
      rownames_to_column('cell_id') %>% 
      group_by(group) %>% 
      nest() %>% 
      deframe() %>% 
      map(function(meta){
        meta$DSS <- ifelse(meta$PCD_score1 > quantile(meta$PCD_score1, probs = 0.95), TRUE, FALSE)
        # print(meta$DSS %>% table())
        meta
      }) %>% 
      bind_rows() %>% 
      column_to_rownames('cell_id')
      
    seurat_hs <- AddMetaData(seurat_hs, meta_DSS['DSS'])    
    
    seurat_hs
  })

seurat_ls

seurat_ls$heart@meta.data %>% head()
seurat_list$heart@meta.data %>% head()

seurat_ls$heart$DSS %>% table()
seurat_ls$liver$DSS %>% table()
seurat_ls$gut$DSS %>% table()
seurat_ls$lung$DSS %>% table()
seurat_ls$kidney$DSS %>% table()

saveRDS(seurat_ls, 'R_OUTPUT/seurat_DSS_HSgroup_ls.rds')


## spatial plot -------------
spatial_theme <- 
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_rect(color = "black", size = 1, fill = NA), #边框
        #axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.title = element_text(size = 15),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.title = element_text(hjust = 0.5, size = 20))


allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#8B4513","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#DEB887")
cells_boundary_ls <- readRDS('R_OUTPUT/data_for_analysis/cells_boundary_ls.rds')


all_data_meta <- 
  seurat_ls %>% 
  map(function(seurat_obj){
    seurat_obj@meta.data %>% 
      rownames_to_column('cell_id')
  }) %>% 
  bind_rows()

all_data_meta$DSS %>% table()

p_dimPlot_ls <- 
  cells_boundary_ls %>%
  imap(function(cells_boundary, sample) { # browser()
    if (!(sample %in% unique(all_data_meta$sample))) return(NULL)
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) + 
      geom_polygon(aes(fill = DSS, group = cell_id)) + 
      ggtitle(sample) + 
      # ggrastr::geom_point_rast(size = 0.05 , alpha = 1, raster.dpi = 300)  +  
      scale_fill_manual(values = c('lightgrey', 'red')) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.7))) + 
      coord_equal() +
      scale_y_reverse() + 
      spatial_theme
    
    p
  }) %>% 
  purrr::discard(function(x) length(x) < 1)

p_dimPlot_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path('3. spatial plot',  
                  str_glue("spatial dim plot {grp}_DSS.pdf")), 
        height = 10, 
        width = 10)
    print(p)
    dev.off()
  })  




## density plot ------------
out_dir <- 'DSS_density_plot'
dir.create(out_dir)


# test
dat <- # coord and annotation
  cells_boundary_ls$HS01gut %>% 
  left_join(all_data_meta[c('cell_id', 'DSS')], by = 'cell_id')

dat2 <- dat %>% remove_missing()
kde <- MASS::kde2d(dat2$vertex_x, dat2$vertex_y)

?geom_contour

data.frame(
  z = kde$z %>% c(),
  x = kde$x %>% rep(times = length(kde$y)),
  y = kde$y %>% rep(each  = length(kde$x))
) %>%
  ggplot(aes(x = x, y = y, z = z)) +
  geom_contour() +
  theme_void()


kde_ls <- 
  dat %>% 
  remove_missing() %>%
  group_by(DSS) %>%
  nest() %>%
  deframe() %>%
  
  map(function(x){
    kde <- MASS::kde2d(x$vertex_x, x$vertex_y)
    data.frame(
      z = kde$z %>% c(),
      x = kde$x %>% rep(times = length(kde$y)),
      y = kde$y %>% rep(each  = length(kde$x))
    )
  })


kde_df <- 
  kde_ls %>%
  enframe("DSS") %>%
  unnest(value)

kde_df %>% 
  ggplot(aes(x = x, y = y, z = z, color = DSS)) +
  geom_contour() +
  theme_void()

kde_df %>% 
  ggplot(aes(x = x, y = -y, z = z)) +
  geom_contour(aes(color = after_stat(nlevel))) +
  scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(Inf, "RdBu"))) +
  facet_wrap(~ DSS) +
  theme_void()
ggsave(file.path(out_dir, "contour plot - levels.pdf"), width = 20, height = 8)

dat2 %>%
  ggplot(aes(x = vertex_x, -vertex_y)) + 
  # geom_point(aes(color = DSS))
  geom_density_2d( ) +
  facet_wrap(~ DSS)


dat2 %>% 
  arrange(desc(DSS)) %>%
  # head(100000) %>%
  ggplot(aes(x = vertex_x, -vertex_y, fill = DSS)) + 
  # geom_point(aes(color = DSS)) +
  # geom_tile(aes(fill = DSS, height = 10, width = 10), alpha = 0.5) + 
  geom_polygon(aes(fill = DSS, group = cell_id), alpha = 0.5) +
  scale_fill_manual(values = c(
    'TRUE' = 'red',
    'FALSE' = 'lightgrey')
  ) +
  
  geom_density2d(
    data = ~ dplyr::filter(., DSS == "TRUE"), 
    aes(color = DSS),
    contour_var = "density",
    show.legend = FALSE,
    # n = 250
    # h = 500
  ) +
  coord_equal() +
  theme_void() 
  # facet_wrap(~ DSS)

ggsave(file.path(out_dir, "contour plot - combine.pdf"), width = 10, height = 8)
# test




contour_plot_ls <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){
    if (!(sample %in% unique(all_data_meta$sample))) return(NULL)
    
    dat <- # coord and annotation
      cells_boundary %>% 
      left_join(all_data_meta[c('cell_id', 'DSS')], by = 'cell_id')
    
    dat2 <- dat %>% remove_missing()
    
    dat2 %>% 
      arrange(desc(DSS)) %>%
      # head(100000) %>%
      ggplot(aes(x = vertex_x, -vertex_y, fill = DSS)) + 
      # geom_point(aes(color = DSS)) +
      # geom_tile(aes(fill = DSS, height = 10, width = 10), alpha = 0.5) + 
      geom_polygon(aes(fill = DSS, group = cell_id), alpha = 0.5) +
      scale_fill_manual(values = c(
        'TRUE' = 'red',
        'FALSE' = 'lightgrey')
      ) +
      
      geom_density2d(
        data = ~ dplyr::filter(., DSS == "TRUE"), 
        aes(color = DSS),
        contour_var = "density",
        show.legend = FALSE,
        # n = 250
        # h = 500
      ) +
      coord_equal() +
      theme_void() 
  }) %>% 
  purrr::discard(function(x) length(x) < 1)

contour_plot_ls %>% 
  iwalk(function(contour_plot, sample){
    pdf(file.path(out_dir, str_glue("contour plot - combine - {sample}.pdf")), width = 10, height = 8)
    print(contour_plot)
    dev.off()
  })


contour_plot_ls2 <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){
    if (!(sample %in% unique(all_data_meta$sample))) return(NULL)
    
    dat <- # coord and annotation
      cells_boundary %>% 
      left_join(all_data_meta[c('cell_id', 'DSS', 'PCD_score1')], by = 'cell_id') %>% 
      dplyr::rename(PCD_score = PCD_score1)
    
    dat2 <- dat %>% remove_missing()
    
    dat2 %>% 
      arrange(desc(DSS)) %>%
      # head(100000) %>%
      ggplot(aes(x = vertex_x, -vertex_y, fill = PCD_score, color = DSS)) + 
      # geom_point(aes(color = DSS)) +
      # geom_tile(aes(fill = DSS, height = 10, width = 10), alpha = 0.5) + 
      geom_polygon(aes(fill = PCD_score, color = DSS, group = cell_id), size = 0.1, alpha = 0.5) + 
      scale_fill_gradient2(
        low = "white",             
        high = "red",              # 高表达
        midpoint = 0,  # 中间值
        na.value = "lightgrey"     # NA值的颜色
      ) +
      scale_color_manual(values = c(
        'TRUE' = 'red',
        'FALSE' = 'lightgrey'),
        guide = guide_legend(override.aes = list(linewidth = 1, fill = NA))
      ) +
      
      geom_density2d(
        data = ~ dplyr::filter(., DSS == "TRUE"), 
        aes(color = DSS),
        contour_var = "density",
        show.legend = FALSE,
        # n = 250
        # h = 500
      ) +
      coord_equal() +
      theme_void() 
  }) %>% 
  purrr::discard(function(x) length(x) < 1)

contour_plot_ls2 %>% 
  iwalk(function(contour_plot, sample){
    pdf(file.path(out_dir, str_glue("contour plot - combine - {sample} - score.pdf")), width = 10, height = 8)
    print(contour_plot)
    dev.off()
  })


contour_plot_ls3 <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){
    if (!(sample %in% unique(all_data_meta$sample))) return(NULL)
    
    dat <- # coord and annotation
      cells_boundary %>% 
      left_join(all_data_meta[c('cell_id', 'DSS', 'PCD_score1')], by = 'cell_id') %>% 
      dplyr::rename(PCD_score = PCD_score1)
    
    dat2 <- dat %>% remove_missing()
    
    dat2 %>% 
      arrange(desc(DSS)) %>%
      # head(100000) %>%
      ggplot(aes(x = vertex_x, -vertex_y, fill = PCD_score, color = DSS)) + 
      # geom_point(aes(color = DSS)) +
      # geom_tile(aes(fill = DSS, height = 10, width = 10), alpha = 0.5) + 
      geom_polygon(aes(fill = PCD_score, color = DSS, group = cell_id), size = 0.1, alpha = 0.5) + 
      scale_fill_gradient2(
        low = "white",             
        high = "red",              # 高表达
        midpoint = 0,  # 中间值
        na.value = "lightgrey"     # NA值的颜色
      ) +
      scale_color_manual(values = c(
        'TRUE' = 'red',
        'FALSE' = 'lightgrey'),
        guide = guide_legend(override.aes = list(linewidth = 1, fill = NA))
      ) +
      
      # geom_density2d(
      #   data = ~ dplyr::filter(., DSS == "TRUE"), 
      #   aes(color = DSS),
      #   contour_var = "density",
      #   show.legend = FALSE,
      #   # n = 250
      #   # h = 500
      # ) +
      coord_equal() +
      theme_void() 
  }) %>% 
  purrr::discard(function(x) length(x) < 1)

contour_plot_ls3 %>% 
  iwalk(function(contour_plot, sample){
    pdf(file.path(out_dir, str_glue("DSS - {sample} - score.pdf")), width = 10, height = 8)
    print(contour_plot)
    dev.off()
  })




## boxplot / annotation --------------
boxplot_DSS_prop_annotation <- 
  all_data_meta %>% 
  select(DSS, annotation, sample, group) %>% 
  filter(DSS == 'TRUE') %>% 
  group_by(group) %>% 
  nest() %>% 
  deframe() %>% 
  imap(function(dat, group){
    dat %>% 
      group_by(sample) %>% 
      nest() %>% 
      deframe() %>% 
      map(function(dat2){
        dat2 %>% 
          group_by(annotation) %>% 
          summarise(count = n()) %>% 
          mutate(freq = count/sum(count))
      }) %>% 
      bind_rows() %>% 
      
      ggplot(aes(x = annotation, y = freq, fill = annotation)) + 
      geom_boxplot(outliers = F) + 
      geom_jitter(
        width = 0.2, 
        height = 0, 
        alpha = 0.5, 
        size = 1, 
        color = "steelblue"
      ) + 
      # ggpubr::geom_pwc(tip.length = 0, hide.ns = F, label = 'p.signif') + 
      scale_fill_manual(values = allcolour) + 
      theme_bw() + 
      theme(
        strip.text = element_text(size = 12), 
        axis.text.x = element_text(angle = 90)
      ) + 
      ggtitle(str_glue('Cell type of DSS - {group}'))
  }) %>% 
  wrap_plots(n = 3)
  
ggsave(boxplot_DSS_prop_annotation, filename = file.path(out_dir, 'boxplot_DSS_prop_by_annotation.pdf'), width = 15, height = 15)  
  

## boxplot / cell_type --------------
boxplot_DSS_prop_cell_type <- 
  all_data_meta %>% 
  select(DSS, subtype, sample, group) %>% 
  dplyr::rename(cell_type = subtype) %>% 
  filter(DSS == 'TRUE') %>% 
  group_by(group) %>% 
  nest() %>% 
  deframe() %>% 
  imap(function(dat, group){
    dat %>% 
      group_by(sample) %>% 
      nest() %>% 
      deframe() %>% 
      map(function(dat2){
        dat2 %>% 
          group_by(cell_type) %>% 
          summarise(count = n()) %>% 
          mutate(freq = count/sum(count))
      }) %>% 
      bind_rows() %>% 
      
      ggplot(aes(x = cell_type, y = freq, fill = cell_type)) + 
      geom_boxplot(outliers = F) + 
      geom_jitter(
        width = 0.2, 
        height = 0, 
        alpha = 0.5, 
        size = 1, 
        color = "steelblue"
      ) + 
      # ggpubr::geom_pwc(tip.length = 0, hide.ns = F, label = 'p.signif') + 
      scale_fill_manual(values = allcolour) + 
      theme_bw() + 
      theme(
        strip.text = element_text(size = 12), 
        axis.text.x = element_text(angle = 90)
      ) + 
      ggtitle(str_glue('Cell type of DSS - {group}'))
  }) %>% 
  wrap_plots(n = 3)

ggsave(boxplot_DSS_prop_cell_type, filename = file.path(out_dir, 'boxplot_DSS_prop_by_cell_type.pdf'), width = 15, height = 15)  


## boxplot / subtype --------------
boxplot_DSS_prop_subtype <- 
  all_data_meta %>% 
  select(DSS, subtype_re, sample, group) %>% 
  dplyr::rename(subtype = subtype_re) %>% 
  filter(DSS == 'TRUE') %>% 
  filter(!is.na(subtype)) %>% 
  group_by(group) %>% 
  nest() %>% 
  deframe() %>% 
  imap(function(dat, group){ browser()
    dat %>% 
      group_by(sample) %>% 
      nest() %>% 
      deframe() %>% 
      map(function(dat2){
        dat2 %>% 
          group_by(subtype) %>% 
          summarise(count = n()) %>% 
          mutate(freq = count/sum(count))
      }) %>% 
      bind_rows() %>% 
      
      ggplot(aes(x = subtype, y = freq, fill = subtype)) + 
      geom_boxplot(outliers = F) + 
      geom_jitter(
        width = 0.2, 
        height = 0, 
        alpha = 0.5, 
        size = 1, 
        color = "steelblue"
      ) + 
      # ggpubr::geom_pwc(tip.length = 0, hide.ns = F, label = 'p.signif') + 
      scale_fill_manual(values = allcolour) + 
      theme_bw() + 
      theme(
        strip.text = element_text(size = 12), 
        axis.text.x = element_text(angle = 90)
      ) + 
      ggtitle(str_glue('Cell type of DSS - {group}'))
  })

ggsave(boxplot_DSS_prop_subtype$HS_liver, filename = file.path(out_dir, 'boxplot_DSS_prop_by_subtype.pdf'), width = 5, height = 5)  




## spatial plot / cell_type -------------
DSS_spatial_plot_ls <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){ # browser()
    if (!(sample %in% unique(all_data_meta$sample))) return(NULL)
    
    dat <- # coord and annotation
      cells_boundary %>% 
      left_join(all_data_meta[c('cell_id', 'DSS', 'subtype')], by = 'cell_id') %>% 
      mutate(
        cell_type = ifelse(DSS == 'TRUE', as.character(subtype), 'non-DSS'), 
        cell_type = factor(cell_type, levels = c(unique(as.character(subtype)), 'non-DSS'))
      )
    
    dat2 <- dat %>% remove_missing()
    
    dat2 %>% 
      arrange(desc(DSS)) %>%
      # head(100000) %>%
      ggplot(aes(x = vertex_x, -vertex_y, fill = cell_type, color = DSS)) + 
      # geom_point(aes(color = DSS)) +
      # geom_tile(aes(fill = DSS, height = 10, width = 10), alpha = 0.5) + 
      geom_polygon(aes(fill = cell_type, color = DSS, group = cell_id), size = 0.1, alpha = 0.5) + 
      scale_fill_manual(values = c(allcolour[1:length(levels(dat2$cell_type)) - 1], 'lightgrey'), drop = FALSE) +
      scale_color_manual(values = c(
        'TRUE' = 'red',
        'FALSE' = 'grey'),
        guide = guide_legend(override.aes = list(linewidth = 1, fill = NA))
      ) +
      
      # geom_density2d(
      #   data = ~ dplyr::filter(., DSS == "TRUE"), 
      #   aes(color = DSS),
      #   contour_var = "density",
      #   show.legend = FALSE,
      #   # n = 250
      #   # h = 500
      # ) +
      coord_equal() +
      theme_void() 
  }) %>% 
  purrr::discard(function(x) length(x) < 1)

DSS_spatial_plot_ls %>% 
  iwalk(function(contour_plot, sample){
    pdf(file.path(out_dir, str_glue("DSS - {sample} - cell_type.pdf")), width = 10, height = 8)
    print(contour_plot)
    dev.off()
  })
 

## spatial plot / subtype -------------
DSS_spatial_plot_subtype_ls <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){
    if (!(sample %in% unique(all_data_meta$sample))) return(NULL)
    
    dat <- # coord and annotation
      cells_boundary %>% 
      left_join(all_data_meta[c('cell_id', 'DSS', 'subtype_re')], by = 'cell_id') %>% 
      filter(!is.na(subtype_re)) %>% 
      mutate(
        subtype = 
          ifelse(DSS == 'TRUE', as.character(subtype_re), 'non-DSS') %>% 
          str_replace_all(c('Netrophil' = 'Neutrophil')), 
        subtype = factor(subtype, levels = c('B cell', 'Cholangiocyte', 'Endothelial', 'Fibroblast', 'Macrophage', 
                                             'Neutrophil', 'Pericentral-hepatocyte', 'Periportal-hepatocyte', 
                                             'Pericyte', 'non-DSS'))
        )
    
    if (dim(dat)[1] == 0) return(NULL)
    
    # browser()
    dat2 <- dat %>% remove_missing()
    
    dat2 %>% 
      arrange(desc(DSS)) %>%
      # head(100000) %>%
      ggplot(aes(x = vertex_x, -vertex_y, fill = subtype, color = DSS)) + 
      # geom_point(aes(color = DSS)) +
      # geom_tile(aes(fill = DSS, height = 10, width = 10), alpha = 0.5) + 
      geom_polygon(aes(fill = subtype, color = DSS, group = cell_id), size = 0.1, alpha = 0.5) + 
      scale_fill_manual(values = c(allcolour[1:9], 'lightgrey'), drop = FALSE) +
      scale_color_manual(values = c(
        'TRUE' = 'red',
        'FALSE' = 'grey'),
        guide = guide_legend(override.aes = list(linewidth = 1, fill = NA))
      ) +
      
      # geom_density2d(
      #   data = ~ dplyr::filter(., DSS == "TRUE"), 
      #   aes(color = DSS),
      #   contour_var = "density",
      #   show.legend = FALSE,
      #   # n = 250
      #   # h = 500
      # ) +
      coord_equal() +
      theme_void() 
  }) %>% 
  purrr::discard(function(x) length(x) < 1)

DSS_spatial_plot_subtype_ls %>% 
  iwalk(function(contour_plot, sample){
    pdf(file.path(out_dir, str_glue("DSS - {sample} - subtype.pdf")), width = 10, height = 8)
    print(contour_plot)
    dev.off()
  })




# DSS nearest neighbor --------------
source("spatial analysis.R")

## input data ------------
cells_position_ls <- readRDS("R_OUTPUT/data_for_analysis/cells_position_ls.rds")

all_data <- merge(seurat_ls[[1]], seurat_ls[-1])
# all_data_cn_meta
all_data_meta_pos <-
  cells_position_ls %>%
  imap(function(cells_position, slide){#browser()
    
    all_data_meta <- all_data[[]]  
    
    proj_meta_sub <-
      all_data_meta %>%  
      dplyr::select(DSS, sample) %>% 
      rownames_to_column("cell_id") %>%
      dplyr::filter(sample == slide)
    
    cells_position[, 1:3] %>%
      inner_join(proj_meta_sub, by = "cell_id") %>%  
      mutate(DSS = factor(DSS)) %>%
      dplyr::rename(Classification1 = DSS,
                    X = x_centroid,
                    Y = y_centroid,
                    Slide_ID = sample)
  }) %>%
  enframe() %>%
  unnest(cols = value) %>%
  dplyr::select(-name)  
  # column_to_rownames("cell_id")


all_data_meta_pos_subset <- all_data_meta_pos
all_data_meta_pos_subset %>% head

# cell_id             X      Y Classification1 Slide_ID
# <chr>           <dbl>  <dbl> <fct>           <chr>   
#   1 TMA2_aajjehak-1  983. 15512. FALSE           HS1heart
# 2 TMA2_aajjenoh-1 1038. 15420. FALSE           HS1heart
# 3 TMA2_aajjhlmf-1  859. 15497. FALSE           HS1heart
# 4 TMA2_aajjibme-1  890. 15497. FALSE           HS1heart
# 5 TMA2_aajjilgk-1  893. 15518. FALSE           HS1heart
# 6 TMA2_aajjkjab-1  879. 15496. FALSE           HS1heart
# ......


## nn ----------
res_nn_ls <-
  all_data_meta_pos_subset %>%
  group_by(Slide_ID) %>%
  nest() %>%
  deframe() %>% #.[1] %>%
  imap(function(dat, name){
    num_nn_area <- 6
    
      
    # 1. find all nn 5 cells
    dat_query <- dat %>% dplyr::filter(Classification1 == "TRUE")
    print(dim(dat_query))
    res_nn <- RANN::nn2(
      dat[c('X', 'Y')],
      query = dat_query[c('X', 'Y')],
      k = 10
    )
    
    # 2. remove duplicated cells between 5 distances
    nn_idx  <- res_nn$nn.idx[, (1:num_nn_area)]  
    nn_dist <- res_nn$nn.dists[, (1:num_nn_area)] 
    
    nn_df <-
      data.frame(
        query_cell_id = rep(dat_query$cell_id, times = num_nn_area),
        row = rep(1:nrow(nn_idx), times = num_nn_area),
        col = rep(1:ncol(nn_idx), each = nrow(nn_idx)),
        idxs = c(nn_idx),
        dists = c(nn_dist),
        Slide_ID = name
      )

    nn_df_filter <-
      nn_df %>%
      group_by(idxs) %>%
      dplyr::filter(dists == min(dists))
    nn_df_filter$nn_cell_id <- dat$cell_id[nn_df_filter$idxs]
    
    # res_nn$query_cell <- dat_query$cell_id
    # res_nn$all_cell   <- dat$cell_id
    # res_nn
    nn_df_filter
  })

DSS_neighbor_group <- 
  res_nn_ls %>% 
  enframe(name = 'sample') %>% 
  unnest() %>% 
  mutate(DSS_neighbor_group = str_c('G', col)) %>% 
  dplyr::select(nn_cell_id, DSS_neighbor_group) %>% 
  dplyr::rename(cell_id = nn_cell_id)
write.csv(DSS_neighbor_group, file.path(out_dir, 'DSS_neighbor_group.csv'), row.names = F)


# dup_idx <- nn_df$idxs %>% .[duplicated(.)] %>% unique()
# if(length(dup_idx) > 0){
#   nn_df_dup <- 
#     nn_df %>%
#     dplyr::filter(dists %in% dup_idx) %>%
#     group_by(idxs) %>%
#     dplyr::filter(dists == min(dists))
#   
#   nn_df_filter <-
#     nn_df %>%
#     dplyr::filter(!dists %in% dup_idx) %>%
#     rbind(nn_df_dup)
# }else{
#   nn_df_filter <- nn_df
# }




# # 初始化空向量
# G1 <- c()
# G2 <- c()
# G3 <- c()
# G4 <- c()
# G5 <- c()
# 
# result_lists <- 
#   res_nn_ls %>% 
#   walk(function(res_nn){ # browser()
#     dat <- 
#       res_nn %>% 
#       ungroup() %>% 
#       select(col, nn_cell_id) %>% 
#       group_by(col) %>% 
#       nest() %>% 
#       deframe() %>% 
#       map(function(x){
#         x$nn_cell_id %>% unlist()
#       })
#     
#     # 使用 assign 赋值到全局环境
#     assign("G1", c(G1, dat[[1]]), envir = .GlobalEnv)
#     assign("G2", c(G2, dat[[2]]), envir = .GlobalEnv)
#     assign("G3", c(G3, dat[[3]]), envir = .GlobalEnv)
#     assign("G4", c(G4, dat[[4]]), envir = .GlobalEnv)
#     assign("G5", c(G5, dat[[5]]), envir = .GlobalEnv)
#   })
# write.csv(G1, file.path(out_dir, 'G1_cell_id.csv'))  
# write.csv(G2, file.path(out_dir, 'G2_cell_id.csv'))
# write.csv(G3, file.path(out_dir, 'G3_cell_id.csv'))
# write.csv(G4, file.path(out_dir, 'G4_cell_id.csv'))
# write.csv(G5, file.path(out_dir, 'G5_cell_id.csv'))


