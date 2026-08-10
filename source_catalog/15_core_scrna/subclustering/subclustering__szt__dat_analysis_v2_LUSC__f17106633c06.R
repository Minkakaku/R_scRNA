library(tidyverse)
library(Seurat)

out_dir_analysis2 <- 'LUSC_out_dir_20250922'
dir.create(out_dir_analysis2)

all_data_subcluster2 <- readRDS('all_data_subcluster2.rds')
# 
# cell_cols <- c(
#   'Epithelials'     = "#5050FF",
#   
#   "Tumor C1 SOX2"   = "#003399",
#   "Tumor C2 SDC4"   = "#1A0099",
#   "Tumor C3 MKI67"  = "#0000FF",
#   
#   'AT2 cell'        = "#4775FF",
#   'Endothelial'     = "#5DB1DD",
#   'Pericyte'        = "#0099CC",
#   
#   
#   'B'               = "#F0E685",
#   'plasma'          = "#FFD147",
#   
#   'CD8-Tem'         = "#802268",
#   'CD8-Temra'       = "#990080",
#   'CD4-Treg'        = "#AE1F63",
#   'CD4-Tfh'         = "#D60047",
#   'Cycling T cell'  = "#FF1463",
#   
#   
#   'NF'              = "#99CC00",
#   'dCAF'            = "#6BD76B",
#   'infCAF'          = "#00CC33",
#   'mCAF'            = "#33CC00",
#   'CDH1 tCAF'       = "#00CC99",
#   'SOX2 tCAF'       = "#0091AF",
#   
#   'cDC1'            = "#AE1F63",
#   'cDC2'            = "#FF1463",
#   'Mast cell'                = "#D60047",
#   'Proliferating Macrophage' = "#7A65A5",
#   'CXCL10 Macrophage'        = "#802268",
#   'CYP27A1 Macrophage'       = "#660099",
#   'Macro Macrophage'         = "#990080"
# )

cell_cols <-
  all_data_subcluster2$annotation2 %>%
  unique() %>%
  structure(
    names = .,
    ggsci::pal_igv()(length(.))
  )

all_data_cn <- readRDS("CN_analysis - LUSC/all_data_meta_cn.rds")

# 0. major group ------------
gs_h <- msigdbr::msigdbr(
  db_species = "HS",
  species = "Homo sapiens",
  collection = "H"
) %>%
  with(split(gene_symbol, gs_name))

all_data_score <- all_data_subcluster2 %>% AddModuleScore(gs_h)

p_hallmark_heat_major <-
  all_data_score[[]] %>%
  dplyr::select(
    group,
    starts_with("Cluster")
  ) %>%
  group_by(group) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames("group") %>%
  {pheatmap::pheatmap(
    t(.),
    # annotation_col = .['group'],
    # scale = "row",
    border_color = NA,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    labels_row = names(gs_h),
    silent = TRUE
  )}

pdf(
  file.path(out_dir_analysis2, "2 hallmark heatmap - epithelial.pdf"), 
  height = 10, 
  width  = 7
)
p_hallmark_heat_major$gtable %>%  cowplot::plot_grid()
dev.off()


p_dot_hallmark <-
  all_data_score[[]] %>%
  mutate(project = str_c(group, project, sep = ": ")) %>% 
  dplyr::select(
    project,
    starts_with("Cluster")
  ) %>%
  group_by(project) %>%
  summarise(across(everything(), mean)) %>%
  mutate(across(-project, scale)) %>%
  pivot_longer(
    -project,
    names_to  = "geneset",
    values_to = "score"
  ) %>%
  mutate(geneset = fct_recode(
    geneset,
    !!!structure(
      names = names(gs_h),
      str_c("Cluster", 1:50)
    )
  )) %>%
  ggplot(aes(x = project, y = geneset, color = score, alpha = score, size = score)) +
  geom_point() +
  theme_bw() +
  scale_color_gradientn(colours = colorRampPalette((RColorBrewer::brewer.pal(11, "Reds")))(100))

pdf(
  file.path(out_dir_analysis2, "2 hallmark dotplot - epithelial.pdf"), 
  height = 12, 
  width  = 7
)
p_dot_hallmark %>%  cowplot::plot_grid()
dev.off()




## Ro/e -------

roe_out_dir <- out_dir_analysis2

anno_column <- "annotation2"
grp_column  <- "project"

all_data_meta_roe <- all_data_subcluster2[[]]

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
col.order <- ggplot2:::`%||%`(
  levels(all_data_meta_roe[[grp_column]]),
  unique(all_data_meta_roe[[grp_column]])
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
    height = unit(11, "cm"), 
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
  file.path(roe_out_dir, "0 cell props boxplot - Roe plot - all cells.pdf"), 
  width = 8, 
  height = 6
)
ComplexHeatmap::draw(roi_hm, annotation_legend_list = lgd)
dev.off()

## Ro/e per CN -----------
all_data_cn %>%
  group_by(CN2) %>%
  nest() %>%
  deframe() %>%
  imap(function(df, name){
    roe_out_dir <- out_dir_analysis2
    
    anno_column <- "annotation2"
    grp_column  <- "project"
    
    all_data_meta_roe <- df
    
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
    col.order <- ggplot2:::`%||%`(
      levels(all_data_meta_roe[[grp_column]]),
      unique(all_data_meta_roe[[grp_column]])
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
        height = unit(11, "cm"), 
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
      file.path(roe_out_dir, str_glue("0 cell props boxplot - Roe plot - {name}.pdf")), 
      width = 8, 
      height = 6
    )
    ComplexHeatmap::draw(roi_hm, annotation_legend_list = lgd)
    dev.off()
    
  })



# 1. umap -------------------
p_umap_major <- 
  DimPlot(
    all_data_subcluster2, 
    group.by = "annotation_full",
    reduction = "harmonyumap"
  ) +
  ggsci::scale_color_aaas()


p_umap_sub <- 
  DimPlot(
    all_data_subcluster2, 
    group.by = "annotation2",
    reduction = "harmonyumap"
  ) +
  scale_color_manual(values = cell_cols)

p_umap_grp <- 
  DimPlot(
    all_data_subcluster2, 
    group.by = "group",
    reduction = "harmonyumap"
  ) 

p_umap_project <- 
  DimPlot(
    all_data_subcluster2, 
    group.by = "project",
    reduction = "harmonyumap"
  ) 

# export
ggsave(
  p_umap_major,
  filename = file.path(out_dir_analysis2, "1. umap plot - major annotation.pdf"),
  height = 6,
  width = 8
)
ggsave(
  p_umap_sub,
  filename = file.path(out_dir_analysis2, "1. umap plot - final annotation.pdf"),
  height = 6,
  width = 12
)
ggsave(
  p_umap_grp,
  filename = file.path(out_dir_analysis2, "1. umap plot - group.pdf"),
  height = 6,
  width = 8
)
ggsave(
  p_umap_project,
  filename = file.path(out_dir_analysis2, "1. umap plot - project.pdf"),
  height = 6,
  width = 8
)

## bar plot 
p_col_major <- 
  all_data_subcluster2[[]] %>%
  group_by(disease, group, annotation_full) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = annotation_full)) +
  geom_col() +
  # facet_wrap(~ disease) +
  ggsci::scale_fill_aaas()

p_col_sub <- 
  all_data_subcluster2[[]] %>%
  group_by(disease, group, annotation2) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = annotation2)) +
  geom_col() +
  # facet_wrap(~ disease) +
  scale_fill_manual(values = cell_cols)

# export
ggsave(
  p_col_major,
  filename = file.path(out_dir_analysis2, "1. bar plot - major - group.pdf"),
  height = 6,
  width  = 8
)

ggsave(
  p_col_sub,
  filename = file.path(out_dir_analysis2, "1. bar plot - final - group.pdf"),
  height = 6,
  width  = 8
)


## mcrophage ---------------

### bar plot -------
p_col_mac <- 
  sp_mac [[]] %>% 
  dplyr::filter(
    str_detect(annotation_mac, "like|low", negate = TRUE)
  ) %>% 
  group_by(project, annotation_mac) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = project, y = freq, fill = annotation_mac)) +
  geom_col() +
  # facet_wrap(~ disease) +
  ggsci::scale_fill_aaas()

# export
ggsave(
  p_col_mac,
  filename = file.path(out_dir_analysis2, "1. bar plot - macrophage - project.pdf"),
  height = 6,
  width  = 8
)



# sp_mac <- readRDS('sp_subcluster/sp_mac.rds')
sp_mac$disease <-
  ifelse(
    sp_mac$sample %in% c("A_2", "A_14", "A_39", "B_2"),
    'LUSC',
    'LUAD'
  )

p_umap_mac <- 
  DimPlot(
    sp_mac, 
    group.by = "annotation_mac",
    reduction = "harmonyumap"
  ) +
  ggsci::scale_color_igv()

# export
ggsave(
  p_umap_mac,
  filename = file.path(out_dir_analysis2, "1. umap plot - macrophage - anno.pdf"),
  height = 6,
  width  = 8
)

# bar plot
p_col_mac <- 
  sp_mac[[]] %>%
  # dplyr::filter(!str_detect(annotation_mac, "like|low")) %>%
  group_by(disease, group, annotation_mac) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = annotation_mac)) +
  geom_col() +
  # facet_wrap(~ disease) +
  ggsci::scale_fill_igv()


p_col_mac2 <- 
  sp_mac[[]] %>%
  group_by(disease, group, annotation_mac) %>%
  summarise(count = n()) %>%
  dplyr::filter(
    annotation_mac != "low quality cells",
    !str_detect(annotation_mac, "-like")
  ) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = annotation_mac)) +
  geom_col() +
  facet_wrap(~ disease) +
  ggsci::scale_fill_aaas()

# export
ggsave(
  p_col_mac,
  filename = file.path(out_dir_analysis2, "1. bar plot - macrophage - group.pdf"),
  height = 6,
  width  = 8
)

ggsave(
  p_col_mac2,
  filename = file.path(out_dir_analysis2, "1. bar plot - macrophage - group - remove others cells.pdf"),
  height = 6,
  width  = 8
)





### de analysis ----------
sp_mac2 <- sp_mac[, !str_detect(sp_mac$annotation_mac, "low|-like")]



## dotplot
p_dot_mac <-
  DotPlot(
    sp_mac2,
    features = c(
      "CD14", "CD68",
      "CCL18", "MMP9",  "MMP12", "CXCL10", "TOP2A", "MKI67",
      "CLEC9A", "IDO1", "LAMP3", "CLEC10A", "FCER1A",  "CD1C",
      "MS4A2", "KIT"
    ),
    group.by = "annotation_mac"
  ) +
  RotatedAxis()

ggsave(
  p_dot_mac,
  filename = file.path(out_dir_analysis2, "1. dotplot - macrohage.pdf"),
  height = 4, width = 12
)

#### 1.1 major group de analysis ---------

all_marker_mac <- 
  sp_mac2 %>%
  JoinLayers() %>%
  FindAllMarkers(
    group.by = "annotation_mac",
    return.thresh = 1.1
  )

# export
all_marker_mac %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir_analysis2, "1 de analysis - macrophage - table.xlsx")
  )

# scatter
p_volcanoplot_mac_ls <- 
  list(Myeloid = all_marker_mac) %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 1
    dodge_width <- 0.1
    
    df_scatter <-
      dat %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      slice_sample(prop = 1) %>%
      group_by(cluster) %>%
      
      mutate(
        x_delta = 1:n(),
        x_delta = scales::rescale(x_delta, to = c(dodge_width, 1 - dodge_width)),
        xx = x_base - 1 + x_delta,
        # y = ifelse(avg_log2FC > 0, avg_log2FC + rect_height, avg_log2FC - rect_height)
        y = avg_log2FC
      ) %>%
      mutate(
        group = ifelse(p_val_adj > 0.01, "p > 0.01", "p < 0.01"),
        group = factor(group, levels = c("p > 0.01", "p < 0.01"))
      ) %>%
      dplyr::arrange(group)
    
    # df_label 
    df_label_up <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_max(avg_log2FC, n = 10)
    
    df_label_down <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_min(avg_log2FC, n = 10)
    
    
    
    # grp rect  
    df_grp_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise() %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1, xmax = x_base) %>%
      mutate(ymin = -rect_height, ymax = rect_height)
    
    
    # background rect
    df_bg_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise(ymax = max(avg_log2FC) + rect_height, ymin = min(avg_log2FC) - rect_height) %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1 + dodge_width, xmax = x_base - dodge_width)
    
    
    
    df_scatter %>%
      ggplot(aes(x = xx, y = y)) +
      geom_rect(
        data = df_bg_rect, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        inherit.aes = FALSE
      ) +
      
      geom_point(aes(color = group), size = 1) +
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_text(
        data = df_grp_rect %>% mutate(xx = 1/2 * (xmin + xmax)),
        aes(x = xx, y = 0, label = cluster),
        inherit.aes = FALSE
      ) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
      ggsci::scale_fill_igv() +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text.x = element_blank()
      ) +
      labs(
        x = NULL,
        y = "Log2FC",
        title = name
      )
    
  })


pdf(
  file.path(out_dir_analysis2, "1 de analysis - macrophage.pdf"), 
  height = 5, 
  width  = 15
)
p_volcanoplot_mac_ls %>%  walk(print)
dev.off()

#### 1.2 enrichment analysis -----------
all_gene_id <- 
  clusterProfiler::bitr(
    rownames(sp_mac2),
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = "org.Hs.eg.db"
  ) %>%
  .$ENTREZID

all_geneid_mac_ls <-
  list(Myeloid = all_marker_mac ) %>%
  map(function(dat){
    gene <-
      dat %>%
      dplyr::filter(avg_log2FC > 0.25, p_val_adj < 0.05, pct.1 > 0.1) %>%
      group_by(cluster) %>%
      nest %>%
      deframe() %>%
      map(~{
        clusterProfiler::bitr(
          .x$gene,
          fromType = "SYMBOL",
          toType   = "ENTREZID",
          OrgDb    = "org.Hs.eg.db"
        ) %>%
          .$ENTREZID
      })
    
    gene  
  })

all_geneid_mac_ls <- all_geneid_mac_ls$Myeloid


all_enrich_mac_ls <- 
  all_geneid_mac_ls %>%
  map(function(x){
    clusterProfiler::enrichGO(
      gene = x,
      ont = "BP",
      OrgDb = "org.Hs.eg.db",
      universe = all_gene_id,
      pvalueCutoff  = 1, # 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff  = 1 # 0.2
    )
  })

# enrichment 
all_enrich_mac_ls <-
  all_enrich_mac_ls %>%
  map(function(x){ # browser()
    xx <- x@result$geneID %>% str_split("/")
    
    all_genes <- unlist(xx)
    gene_db <- 
      clusterProfiler::bitr(
        all_genes,
        fromType = "ENTREZID",
        toType   = "SYMBOL",
        OrgDb    = "org.Hs.eg.db"
      ) 
    
    xx2 <-
      xx %>%
      map_chr(function(g){
        str_c(gene_db$SYMBOL[match(g, gene_db$ENTREZID)], collapse = ",")
      })
    x@result$geneSymbol <- xx2
    x
  })


# enrichment plot
p_all_enrich_mac_ls <-
  all_enrich_mac_ls %>%
  imap(function(enrich, name){ # browser()
    enrichplot::dotplot(
      enrich,
      x = "GeneRatio",
      color = "p.adjust",
      showCategory = 10,
      label_format = 40
    ) + 
      # scale_fill_gradient2(low = "gray", high = "red") +
      ggtitle(name)
  })


# export
all_enrich_mac_ls %>%
  map('result') %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir_analysis2, "1 enrich analysis - macrophage - table - v2 - universe.xlsx")
  )


pdf(
  file.path(out_dir_analysis2, "1 enrich analysis - macropage - v2 - universe.pdf"), 
  height = 6.5, 
  width  = 8
)
p_all_enrich_mac_ls %>%  walk(print)
dev.off()


#### 1.3 macrophage score ------------------
mac_pathway_ls <- list(
  M1 = c("IL12A", "IL23A", "TNF", "IL6", "CD86", "HLA-DMA", "HLA-DMB", 
         "HLA-DOA", "HLA-DOB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", 
         "HLA-DQB1", "HLA-DQB2", "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "IL1B", 
         "MARCO", "NOS2", "FCGR1A", "CD80", "CXCL9", "CXCL10", "CXCL11", 
         "CD86", "IL1A", "IL1B", "IL6", "TNF", "CCL5", "IRF5", "IRF1", 
         "CD40", "IDO1", "KYNU", "CCR7", "PTPRC", "CD68", "CSF1R", "LY75", 
         "CD14"),
  M2 = c("ARG1", "ARG2", "IL10", "FCGR2A", "CD163", "FCER2", "CD200R1", 
         "PDCD1LG2", "CD274", "MARCO", "CSF1R", "MRC1", "IL1RN", "IL1R2", 
         "IL4R", "CCL4", "CCL13", "CCL20", "CCL17", "CCL18", "CCL22", 
         "CCL24", "LYVE1", "VEGFA", "VEGFB", "VEGFC", "VEGFD", "EGF", 
         "CTSA", "CTSB", "CST3", "CTSD", "TGFB1", "TGFB2", "TGFB3", "MMP14", 
         "MMP19", "MMP9", "CLEC7A", "WNT7B", "FASLG", "TNFSF12", "TNFSF8", 
         "CD276", "VTCN1", "MSR1", "FN1", "IRF4", "PTPRC", "CD68", "CSF1R", 
         "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "LY75", "CD14"),
  Phagocytosis = c("MRC1", "CD163", "MERTK", "C1QB"),
  Angiogenesis = c("CCND2", "CCNE1", "CD44", "CXCR4", "E2F3", "EDN1", "EZH2", 
                   "FGF18", "FGFR1", "FYN", "HEY1", "ITGAV", "JAG1", "JAG2", "MMP9", 
                   "NOTCH1", "PDGFA", "PTK2", "SPP1", "STC1", "TNFAIP6", "TYMP", 
                   "VAV2", "VCAN", "VEGFA")
  
)

sp_mye_score <- 
  sp_mac2 %>%
  AddModuleScore(features = mac_pathway_ls) 
colnames(sp_mye_score@meta.data)[match(paste0("Cluster", 1:4), colnames(sp_mye_score@meta.data))] <- names(mac_pathway_ls)


pdf(
  file.path(out_dir_analysis2, "1 pathway score - macrophage pathway score feature plot.pdf"), 
  height = 4, 
  width  = 24
)
dim_path_ls <- 
  FeaturePlot(
    sp_mye_score, 
    features = c("M1", "M2", "Phagocytosis", "Angiogenesis"), 
    order = T, 
    reduction = "harmonyumap",
    combine = FALSE) %>%
  map(function(p) 
    p + 
      scale_color_gradientn(
        colours = Seurat:::SpatialColors(50), 
        # transform = "log1p"
      )
  )

dim_path_ls[['CellType']] <- 
  DimPlot(
    sp_mye_score,
    reduction = "harmonyumap",
    group.by = "annotation_mac"
  ) & 
  ggsci::scale_color_aaas()
patchwork::wrap_plots(dim_path_ls, nrow = 1)
dev.off()

## Vlnplot
VlnPlot(
  sp_mye_score,
  features = c("M1", "M2"),
  group.by = "annotation_mac",
  pt.size = 0
)

VlnPlot(
  sp_mye_score,
  features = c('Phagocytosis', 'Angiogenesis'),
  group.by = "annotation_mac",
  pt.size = 0
)
#### 1.4 Ro/e ---------------

roe_out_dir <- out_dir_analysis2

anno_column <- "annotation_mac"
grp_column  <- "project"

all_data_meta_roe <- 
  sp_mac[[]] %>% 
  dplyr::filter(str_detect(annotation_mac, "like|low", negate = TRUE)) %>%
  mutate(annotation_mac = fct_drop(annotation_mac))


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
col.order <- ggplot2:::`%||%`(
  levels(all_data_meta_roe[[grp_column]]),
  unique(all_data_meta_roe[[grp_column]])
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
    height = unit(11, "cm"), 
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
  file.path(roe_out_dir, "1 cell props boxplot - Roe plot - macrophage.pdf"), 
  width = 8, 
  height = 6
)
ComplexHeatmap::draw(roi_hm, annotation_legend_list = lgd)
dev.off()

## epi ------------------
sp_epi <- readRDS('sp_subcluster/sp_epi.rds')


p_umap_epi <- 
  DimPlot(
    sp_epi, 
    group.by = "annotation_epi",
    reduction = "harmonyumap"
  ) +
  ggsci::scale_color_igv()

# export
ggsave(
  p_umap_epi,
  filename = file.path(out_dir_analysis2, "1. umap plot - epi - anno.pdf"),
  height = 6,
  width  = 8
)

# bar plot
p_col_epi <- 
  sp_epi[[]] %>%
  group_by(disease, group, annotation_epi) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = annotation_epi)) +
  geom_col() +
  facet_wrap(~ disease) +
  ggsci::scale_fill_aaas()


# export
ggsave(
  p_col_epi,
  filename = file.path(out_dir_analysis2, "1. bar plot - epi - group.pdf"),
  height = 6,
  width  = 8
)


## dotplot
p_dot_epi <-
  DotPlot(
    sp_epi,
    features = c(
      'EPCAM',
      'APOL1', 'SLC2A1', 'TYMP', 'APOL2','ALG11', 
      'BIRC5', 'RRM1', 'SMC4', 'TUBB', 'NTS',
      'CEBPA', 'CES1', 'CLDN1', 'NDRG4', 'GJA1', 'FGFBP2',
      'CDKN1A', 'PLAU', 'SDC4', 'A2ML1', 'FUT3', 'IL1RN'
    ),
    group.by = "annotation_epi"
  ) +
  RotatedAxis()

ggsave(
  p_dot_epi,
  filename = file.path(out_dir_analysis2, "1. dotplot - epithelial.pdf"),
  height = 4, 
  width = 12
) 


### 2.1 epi hallmark ---------
gs_h <- msigdbr::msigdbr(
    db_species = "HS",
    species = "Homo sapiens",
    collection = "H"
  ) %>%
  with(split(gene_symbol, gs_name))

sp_epi_score <- sp_epi %>% AddModuleScore(gs_h)

p_hallmark_heat <-
  sp_epi_score[[]] %>%
  dplyr::select(
    annotation_epi,
    starts_with("Cluster")
  ) %>%
  group_by(annotation_epi) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames("annotation_epi") %>%
  {pheatmap::pheatmap(
    t(.),
    # annotation_col = .['group'],
    scale = "row",
    border_color = NA,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    labels_row = names(gs_h),
    silent = TRUE
  )}

pdf(
  file.path(out_dir_analysis2, "2 hallmark heatmap - epithelial.pdf"), 
  height = 10, 
  width  = 7
)
p_hallmark_heat$gtable %>%  cowplot::plot_grid()
dev.off()

# remove CEBPA Epithelial
p_hallmark_heat2 <-
  sp_epi_score[[]] %>%
  dplyr::select(
    annotation_epi,
    starts_with("Cluster")
  ) %>%
  dplyr::filter(annotation_epi != "CEBPA Epithelial") %>% 
  group_by(annotation_epi) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames("annotation_epi") %>%
  {pheatmap::pheatmap(
    t(.),
    # annotation_col = .['group'],
    scale = "row",
    border_color = NA,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    labels_row = names(gs_h),
    silent = TRUE
  )}

pdf(
  file.path(out_dir_analysis2, "2 hallmark heatmap - epithelial - remove CEBPA Epithelial.pdf"), 
  height = 10, 
  width  = 7
)
p_hallmark_heat2$gtable %>%  cowplot::plot_grid()
dev.off()


# group  hallmark 
p_hallmark_heat_grp <-
  sp_epi_score[[]] %>%
  dplyr::select(
    project,
    starts_with("Cluster")
  ) %>%
  # dplyr::filter(annotation_epi != "CEBPA Epithelial") %>% 
  group_by(project) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames("project") %>%
  {pheatmap::pheatmap(
    t(.),
    # annotation_col = .['group'],
    scale = "row",
    border_color = NA,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    labels_row = names(gs_h),
    silent = TRUE
  )}

pdf(
  file.path(out_dir_analysis2, "2 hallmark heatmap - group.pdf"), 
  height = 10, 
  width  = 7
)
p_hallmark_heat_grp$gtable %>%  cowplot::plot_grid()
dev.off()

# group  hallmark 
p_hallmark_heat_grp2 <-
  sp_epi_score[[]] %>% 
  dplyr::filter(annotation_epi != "CEBPA Epithelial") %>%
  dplyr::select(
    project,
    starts_with("Cluster")
  ) %>%
  group_by(project) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames("project") %>%
  {pheatmap::pheatmap(
    t(.),
    # annotation_col = .['group'],
    scale = "row",
    border_color = NA,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    labels_row = names(gs_h),
    silent = TRUE
  )}

pdf(
  file.path(out_dir_analysis2, "2 hallmark heatmap - group - remove CEBPA Epithelial.pdf"), 
  height = 10, 
  width  = 7
)
p_hallmark_heat_grp2$gtable %>%  cowplot::plot_grid()
dev.off()

### 2.2  major group de analysis ---------

all_marker_epi <- 
  sp_epi %>%
  JoinLayers() %>%
  FindAllMarkers(
    group.by = "annotation_epi",
    return.thresh = 1.1
  )

# export
all_marker_epi %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir_analysis2, "2 de analysis - epithelial - table.xlsx")
  )

# scatter
p_volcanoplot_epi_ls <- 
  list(Epithelial = all_marker_epi) %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 1
    dodge_width <- 0.1
    
    df_scatter <-
      dat %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      slice_sample(prop = 1) %>%
      group_by(cluster) %>%
      
      mutate(
        x_delta = 1:n(),
        x_delta = scales::rescale(x_delta, to = c(dodge_width, 1 - dodge_width)),
        xx = x_base - 1 + x_delta,
        # y = ifelse(avg_log2FC > 0, avg_log2FC + rect_height, avg_log2FC - rect_height)
        y = avg_log2FC
      ) %>%
      mutate(
        group = ifelse(p_val_adj > 0.01, "p > 0.01", "p < 0.01"),
        group = factor(group, levels = c("p > 0.01", "p < 0.01"))
      ) %>%
      dplyr::arrange(group)
    
    # df_label 
    df_label_up <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_max(avg_log2FC, n = 10)
    
    df_label_down <- df_scatter %>%
      dplyr::filter(as.character(group) == "p < 0.01") %>%
      group_by(cluster) %>% 
      slice_min(avg_log2FC, n = 10)
    
    
    
    # grp rect  
    df_grp_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise() %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1, xmax = x_base) %>%
      mutate(ymin = -rect_height, ymax = rect_height)
    
    
    # background rect
    df_bg_rect <-
      dat %>%
      group_by(cluster) %>%
      summarise(ymax = max(avg_log2FC) + rect_height, ymin = min(avg_log2FC) - rect_height) %>%
      mutate(x_base = as.numeric(factor(cluster))) %>%
      mutate(xmin = x_base - 1 + dodge_width, xmax = x_base - dodge_width)
    
    
    
    df_scatter %>%
      ggplot(aes(x = xx, y = y)) +
      geom_rect(
        data = df_bg_rect, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        inherit.aes = FALSE
      ) +
      
      geom_point(aes(color = group), size = 1) +
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_text(
        data = df_grp_rect %>% mutate(xx = 1/2 * (xmin + xmax)),
        aes(x = xx, y = 0, label = cluster),
        inherit.aes = FALSE
      ) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
      ggsci::scale_fill_igv() +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.text.x = element_blank()
      ) +
      labs(
        x = NULL,
        y = "Log2FC",
        title = name
      )
    
  })


pdf(
  file.path(out_dir_analysis2, "2 de analysis - epithelial.pdf"), 
  height = 5, 
  width  = 15
)
p_volcanoplot_epi_ls %>%  walk(print)
dev.off()

### 2.3 enrichment analysis -----------
# all_gene_id <- 
#   clusterProfiler::bitr(
#     rownames(sp_mac2),
#     fromType = "SYMBOL",
#     toType   = "ENTREZID",
#     OrgDb    = "org.Hs.eg.db"
#   ) %>%
#   .$ENTREZID

all_geneid_epi_ls <-
  list(Epithelial = all_marker_epi ) %>%
  map(function(dat){
    gene <-
      dat %>%
      dplyr::filter(avg_log2FC > 0.25, p_val_adj < 0.05, pct.1 > 0.1) %>%
      group_by(cluster) %>%
      nest %>%
      deframe() %>%
      map(~{
        clusterProfiler::bitr(
          .x$gene,
          fromType = "SYMBOL",
          toType   = "ENTREZID",
          OrgDb    = "org.Hs.eg.db"
        ) %>%
          .$ENTREZID
      })
    
    gene  
  })

all_geneid_epi_ls <- all_geneid_epi_ls$Epithelial


all_enrich_epi_ls <- 
  all_geneid_epi_ls %>%
  map(function(x){
    clusterProfiler::enrichGO(
      gene = x,
      ont = "BP",
      OrgDb = "org.Hs.eg.db",
      universe = all_gene_id,
      pvalueCutoff  = 1, # 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff  = 1 # 0.2
    )
  })

# enrichment 
all_enrich_epi_ls <-
  all_enrich_epi_ls %>%
  map(function(x){ # browser()
    xx <- x@result$geneID %>% str_split("/")
    
    all_genes <- unlist(xx)
    gene_db <- 
      clusterProfiler::bitr(
        all_genes,
        fromType = "ENTREZID",
        toType   = "SYMBOL",
        OrgDb    = "org.Hs.eg.db"
      ) 
    
    xx2 <-
      xx %>%
      map_chr(function(g){
        str_c(gene_db$SYMBOL[match(g, gene_db$ENTREZID)], collapse = ",")
      })
    x@result$geneSymbol <- xx2
    x
  })


# enrichment plot
p_all_enrich_epi_ls <-
  all_enrich_epi_ls %>%
  imap(function(enrich, name){ # browser()
    enrichplot::dotplot(
      enrich,
      x = "GeneRatio",
      color = "p.adjust",
      showCategory = 10,
      label_format = 40
    ) + 
      # scale_fill_gradient2(low = "gray", high = "red") +
      ggtitle(name)
  })


# export
all_enrich_epi_ls %>% 
  map('result') %>% 
  openxlsx::write.xlsx(
    file = file.path(out_dir_analysis2, "2 enrich analysis - epithelial - table - v2 - universe.xlsx")
  )

pdf(
  file.path(out_dir_analysis2, "2 enrich analysis - epithelial - v2 - universe.pdf"), 
  height = 6.5, 
  width  = 8
)
p_all_enrich_epi_ls %>% walk(print)
dev.off()


### 2.4 Ro/e ---------------

roe_out_dir <- out_dir_analysis2

anno_column <- "annotation_epi"
grp_column  <- "project"

all_data_meta_roe <- sp_epi[[]]

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
col.order <- ggplot2:::`%||%`(
  levels(all_data_meta_roe[[grp_column]]),
  unique(all_data_meta_roe[[grp_column]])
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
    height = unit(11, "cm"), 
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
  file.path(roe_out_dir, "2 cell props boxplot - Roe plot - epi cells.pdf"), 
  width = 8, 
  height = 6
)
ComplexHeatmap::draw(roi_hm, annotation_legend_list = lgd)
dev.off()


### 2.5 dotplot of wnt5a and mcam --------
sp_epi[, sp_epi$annotation_epi != "CEBPA Epithelial"] %>%
  DotPlot(
    group.by = "annotation_epi",
    split.by = "group",
    feature = c("WNT5A")
  )





## fibro ---------------------
# sp_fibro <- readRDS("sp_subcluster/sp_fibro.rds")


p_umap_fibro <- 
  DimPlot(
    sp_fibro, 
    group.by = "annotation_fibro",
    reduction = "harmonyumap"
  ) +
  ggsci::scale_color_igv()

# export
ggsave(
  p_umap_fibro,
  filename = file.path(out_dir_analysis2, "1. umap plot - fibro - anno.pdf"),
  height = 6,
  width  = 8
)

# bar plot
p_col_fibro <- 
  sp_fibro[[]] %>%
  group_by(disease, group, annotation_fibro) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = annotation_fibro)) +
  geom_col() +
  facet_wrap(~ disease) +
  ggsci::scale_fill_aaas()


p_col_fibro2 <- 
  sp_fibro[[]] %>%
  group_by(disease, group, annotation_fibro) %>%
  summarise(count = n()) %>%
  dplyr::filter(
    annotation_fibro != "low quality cells",
    !str_detect(annotation_fibro, "-like|Others")
  ) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = annotation_fibro)) +
  geom_col() +
  facet_wrap(~ disease) +
  ggsci::scale_fill_aaas()

# export
ggsave(
  p_col_fibro,
  filename = file.path(out_dir_analysis2, "1. bar plot - fibro - group.pdf"),
  height = 6,
  width  = 8
)

ggsave(
  p_col_fibro2,
  filename = file.path(out_dir_analysis2, "1. bar plot - fibro - group - remove others cells.pdf"),
  height = 6,
  width  = 8
)



## t ------------------
sp_cd4t <- readRDS("sp_subcluster_LUSC/sp_cd4t.rds")
sp_cd8t <- readRDS("sp_subcluster_LUSC/sp_cd8t.rds")

sp_t_merge <- merge(sp_cd4t, sp_cd8t)

sp_t_merge$annotation_t <- 
  ifelse(
    is.na(sp_t_merge$annotation_cd4t),
    sp_t_merge$annotation_cd8t,
    sp_t_merge$annotation_cd4t
  )



DotPlot(
  all_data_subcluster2[, all_data_subcluster2$disease == 'LUSC'],
  features = c('THBS1',	'CD47'),
  group.by = "annotation2"
)
 # p_umap_t_merge <- 
#   DimPlot(
#     sp_t_merge, 
#     group.by = "annotation_t",
#     reduction = "harmonyumap"
#   ) +
#   ggsci::scale_color_igv()
# 
# # export
# ggsave(
#   p_umap_t_merge,
#   filename = file.path(out_dir_analysis2, "1. umap plot - fibro - anno.pdf"),
#   height = 6,
#   width  = 8
# )

### 4.1 bar plot ------------
p_col_t_merge <- 
  sp_t_merge[[]] %>%
  group_by(disease, group, annotation_t) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = annotation_t)) +
  geom_col() +
  facet_wrap(~ disease) +
  ggsci::scale_fill_aaas()


p_col_t_merge2 <- 
  sp_t_merge[[]] %>%
  group_by(disease, group, annotation_t) %>%
  summarise(count = n()) %>%
  dplyr::filter(
    annotation_t != "low quality cells",
    !str_detect(annotation_t, "-like|Other")
  ) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = annotation_t)) +
  geom_col() +
  facet_wrap(~ disease) +
  ggsci::scale_fill_aaas()

# export
ggsave(
  p_col_t_merge,
  filename = file.path(out_dir_analysis2, "1. bar plot - T - group.pdf"),
  height = 6,
  width  = 8
)

ggsave(
  p_col_t_merge2,
  filename = file.path(out_dir_analysis2, "1. bar plot - T - group - remove others cells.pdf"),
  height = 6,
  width  = 8
)

# barplot - project
p_col_t_merge_proj <- 
  sp_t_merge[[]] %>%
  group_by(disease, project, annotation_t) %>%
  summarise(count = n()) %>%
  dplyr::filter(
    annotation_t != "low quality cells",
    !str_detect(annotation_t, "-like|Other")
  ) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = project, y = freq, fill = annotation_t)) +
  geom_col() +
  # facet_wrap(~ disease) +
  ggsci::scale_fill_aaas()


ggsave(
  p_col_t_merge_proj,
  filename = file.path(out_dir_analysis2, "4. bar plot - T - project - remove others cells.pdf"),
  height = 6,
  width  = 8
)
### 4.2 Ro/e ---------------

roe_out_dir <- out_dir_analysis2

anno_column <- "annotation_t"
grp_column  <- "project"

all_data_meta_roe <- 
  sp_t_merge[[]] %>%
  dplyr::filter(str_detect(annotation_t, "like|low|Others", negate = TRUE)) %>%
  mutate(annotation_t = fct_drop(annotation_t))

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
col.order <- ggplot2:::`%||%`(
  levels(all_data_meta_roe[[grp_column]]),
  unique(all_data_meta_roe[[grp_column]])
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
    height = unit(11, "cm"), 
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
  file.path(roe_out_dir, "4 cell props boxplot - Roe plot - t cells.pdf"), 
  width = 8, 
  height = 6
)
ComplexHeatmap::draw(roi_hm, annotation_legend_list = lgd)
dev.off()



# all_data ---------

## hallmark ---------
gs_h <- msigdbr::msigdbr(
  db_species = "HS",
  species = "Homo sapiens",
  collection = "H"
) %>%
  with(split(gene_symbol, gs_name))

all_data_subcluster2 <- all_data_subcluster2 %>% AddModuleScore(gs_h)

p_hallmark_heat_all_grp <-
  all_data_subcluster2[[]] %>%
  dplyr::select(
    project,
    starts_with("Cluster")
  ) %>%
  group_by(project) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames("project") %>%
  {pheatmap::pheatmap(
    t(.),
    # annotation_col = .['group'],
    scale = "row",
    border_color = NA,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    labels_row = names(gs_h),
    silent = TRUE
  )}

pdf(
  file.path(out_dir_analysis2, "hallmark heatmap - all cells.pdf"), 
  height = 10, 
  width  = 7
)
p_hallmark_heat_all_grp$gtable %>%  cowplot::plot_grid()
dev.off()





DotPlot(
  all_data_subcluster2,
  features = c("CXCL16", "CXCR6"),
  group.by = "annotation2"
)


# CCL18_CCR8
DotPlot(
  all_data_subcluster2,
  features = c("CCL18", "CCR8"),
  group.by = "annotation2"
)


# CCL18_CCR8
DotPlot(
  all_data_subcluster2,
  features = c("CCL18", "CCR8"),
  group.by = "project"
)



DotPlot(
  all_data_subcluster2,
  features = c('LGALS9',	'P4HB'),
  group.by = "annotation2"
)











