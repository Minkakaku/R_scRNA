library(tidyverse)
library(Seurat)

# define default-------------
out_dir <- "DEGs_by_group" 
dir.create(out_dir)

allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#8B4513","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#DEB887")

all_data <- readRDS("R_OUTPUT/data_for_analysis/all_data_filtered2.rds")
all_data_cn_meta <- readRDS("CN_result_all_10CN/all_data_meta_cn.rds")
all_data_cn_meta %>% dim()

all_data_meta <- 
  all_data[[]] %>% 
  rownames_to_column('cell_id') %>% 
  left_join(all_data_cn_meta[c('cell_id', 'CN')], by = 'cell_id') %>% 
  column_to_rownames('cell_id')
all_data_meta %>% head()
all_data <- AddMetaData(all_data, all_data_meta['CN'])
all_data@meta.data %>% head()


# setdiff DEGs ------------
P_markers <- 
  read.csv(file.path(out_dir, 'P_tumor_area_DEGs.csv')) %>% 
  column_to_rownames('X')
P_markers %>% head()
(P_markers$p_val > 0.05) %>% table()

A_markers <- 
  read.csv(file.path(out_dir, 'A_tumor_area_DEGs.csv')) %>% 
  column_to_rownames('X')
A_markers %>% head()
(A_markers$p_val > 0.05) %>% table()


A_minus_P <- setdiff(rownames(A_markers), rownames(P_markers))
'CTTN' %in% rownames(A_markers)
'CTTN' %in% rownames(P_markers)

P_minus_A <- setdiff(rownames(P_markers), rownames(A_markers))




# scoring ------- 
pathway<-list(A_minus_P, P_minus_A)
names(pathway)<-c('A_minus_P', 'P_minus_A')

n_na_idx <- pathway %>% map_lgl(function(x) all(is.na(x))) %>% {! .}%>%  which
pathway <- pathway[n_na_idx] %>% map(function(x) x[!is.na(x)])

for ( i in seq_along(names(pathway))){
  all_data <- AddModuleScore(JoinLayers(all_data), pathway[i], name =  names(pathway)[i],nbin=5,ctrl=10)
}

p<-VlnPlot(all_data,
           features = c('A_minus_P1','P_minus_A1'),
           group.by='CN',
           pt.size = 0)


ggsave(
  p,
  filename = file.path(out_dir, '/vlnplot_disease_associated_gene_score.pdf'),
  height = 10,
  width  = 15
) 


all_data_meta <- all_data@meta.data %>% rownames_to_column('cell_id')
score_col <- c('A_minus_P1','P_minus_A1')
boxplot <- 
  all_data_meta %>% 
  select(CN, all_of(score_col)) %>% 
  pivot_longer(all_of(score_col), names_to = 'gene', values_to = 'expr') %>% 
  # filter(expr > 0) %>% 
  ggplot(aes(x = CN, y = expr, fill = CN)) + 
  geom_boxplot(outliers = F) + 
  # geom_jitter(
  #   width = 0.2, 
  #   height = 0, 
  #   alpha = 0.5, 
  #   size = 2, 
  #   color = "steelblue"
  # ) + 
  ggpubr::geom_pwc(tip.length = 0, hide.ns = F, label = 'p.signif') + 
  facet_wrap(~ gene, scale = "free_y", ncol = 2) + 
  scale_fill_manual(values = allcolour) + 
  labs(y = "Expression Level") +
  theme_bw() + 
  theme(
    strip.text = element_text(size = 12), 
    axis.text.x = element_text(angle = 90)
  )

ggsave(file.path(out_dir, 'Boxplot_disease_associated_gene_score.pdf'), width = 10, height = 15)


# A_minus_P DEGs boxplot ----------
all_data_meta <- 
  all_data[[]] %>% 
  rownames_to_column('cell_id')
data_expr <-  
  all_data %>% 
  FetchData(vars = unlist(A_minus_P), layer = 'data') %>% 
  rownames_to_column('cell_id')
stopifnot(all(all_data_meta$cell_id == data_expr$cell_id))

data_expr_full <- left_join(all_data_meta, data_expr, by = 'cell_id')

boxplot <- 
  data_expr_full %>% 
  select(CN, all_of(A_minus_P)) %>% 
  pivot_longer(all_of(A_minus_P), names_to = 'gene', values_to = 'expr') %>% 
  # filter(expr > 0) %>% 
  ggplot(aes(x = CN, y = expr, fill = CN)) + 
  geom_boxplot(outliers = F) + 
  # geom_jitter(
  #   width = 0.2, 
  #   height = 0, 
  #   alpha = 0.5, 
  #   size = 2, 
  #   color = "steelblue"
  # ) + 
  # ggpubr::geom_pwc(tip.length = 0, hide.ns = F, label = 'p.signif') + 
  facet_wrap(~ gene, scale = "free_y", ncol = 5) + 
  scale_fill_manual(values = allcolour) + 
  labs(y = "Expression Level") +
  theme_bw() + 
  theme(
    strip.text = element_text(size = 12), 
    axis.text.x = element_text(angle = 90)
  )

ggsave(file.path(out_dir, 'Boxplot_A_minus_P_genes.pdf'), width = 15, height = 190, limitsize = FALSE)




# P_minus_A DEGs boxplot ----------
all_data_meta <- 
  all_data[[]] %>% 
  rownames_to_column('cell_id')
data_expr <-  
  all_data %>% 
  FetchData(vars = unlist(P_minus_A), layer = 'data') %>% 
  rownames_to_column('cell_id')
stopifnot(all(all_data_meta$cell_id == data_expr$cell_id))

data_expr_full <- left_join(all_data_meta, data_expr, by = 'cell_id')

boxplot <- 
  data_expr_full %>% 
  select(CN, all_of(P_minus_A)) %>% 
  pivot_longer(all_of(P_minus_A), names_to = 'gene', values_to = 'expr') %>% 
  # filter(expr > 0) %>% 
  ggplot(aes(x = CN, y = expr, fill = CN)) + 
  geom_boxplot(outliers = F) + 
  # geom_jitter(
  #   width = 0.2, 
  #   height = 0, 
  #   alpha = 0.5, 
  #   size = 2, 
  #   color = "steelblue"
  # ) + 
  # ggpubr::geom_pwc(tip.length = 0, hide.ns = F, label = 'p.signif') + 
  facet_wrap(~ gene, scale = "free_y", ncol = 5) + 
  scale_fill_manual(values = allcolour) + 
  labs(y = "Expression Level") +
  theme_bw() + 
  theme(
    strip.text = element_text(size = 12), 
    axis.text.x = element_text(angle = 90)
  )

ggsave(file.path(out_dir, 'Boxplot_P_minus_A_genes.pdf'), width = 15, height = 190, limitsize = FALSE)
