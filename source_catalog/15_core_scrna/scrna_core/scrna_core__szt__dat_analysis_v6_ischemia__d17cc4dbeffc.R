library(tidyverse)
library(Seurat)


all_data <- readRDS("all_data_projected.rds")
all_data <- all_data %>% JoinLayers()
DefaultAssay(all_data) <- "RNA"


out_dir_v6 <- "20260303 dat_anlysis"
dir.create(out_dir_v6)



path_ls <-
  readxl::excel_sheets("Supplementary Table 2_41586_2024_7806_MOESM4_ESM.xlsx") %>%
  setNames(., str_trim(str_extract(., "\\s[A-Za-z]*\\s"))) %>%
  setNames(., str_replace(names(.), "Ischemic", "IZ")) %>%
  map(function(x){
    dat <- readxl::read_excel("Supplementary Table 2_41586_2024_7806_MOESM4_ESM.xlsx", sheet = x)
    dat$gene
  }) %>%
  map(str_to_title)

path_ls <- path_ls[c('BZ', 'IZ')]


# signature of papers
path_ls2 <-
  readxl::read_excel("44161_2022_160_MOESM2_ESM.xlsx", sheet = 1) %>%
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map('gene') %>%
  setNames(., fct_recode(names(.), RZ = 'Myh6', BZ1 = 'Ankrd1', BZ2 = 'Xirp2'))

path_ls <- c(path_ls, path_ls2)



# signature from chatgpt
path_ls$Hypoxia  <- c(
  'Vegfa', 'Slc2a1', 'Pdk1', 'Adm', 'Bnip3', 'Egln3', 'Ca9',
  'Ldha', 'Pgk1', 'Pfkfb3', 'Hk2', 'Eno1', 'Tpi1'
)

path_ls$metabolic_ischemia <- c(
  'Hk1', 'Hk2', 'Pfkl', 'Pfkm', 'Pkm', 'Ldha'
)

path_ls$FAO <- c(
  'Cpt1b', 'Acadl', 'Acox1', 'Ppara'
)


# add module score
all_data_score_zone <- 
  all_data %>%
  AddModuleScore(path_ls)
colnames(all_data_score_zone@meta.data)[
  str_detect(colnames(all_data_score_zone@meta.data), 'Cluster')
] <- 
  names(path_ls)




# 1 nn analysis ----------
cell_coords_ls <- readRDS("cells_position_ls.rds")
cell_coords_ls2 <- cell_coords_ls %>% map(~  .[.$cell_id %in% colnames(all_data), ])

## find nn point -----------
k_nn <- 25
res_nn_ls <- 
  cell_coords_ls2 %>%
  map(function(dat){ # browser()
    cn_nn <- 
      RANN::nn2(
        dat[c('x', 'y')], 
        k = k_nn
      )
    
    list(
      nn = cn_nn,
      cell_id = dat$cell_id
    ) 
  })




## combined nn expression -----------
dat_zone_sum_ls <- 
  res_nn_ls %>%
  map(function(res_nn){ # browser()
    # dat
    print('subset data')
    dat_used <- FetchData(
      all_data_score_zone,
      vars = names(path_ls),
      cells = res_nn$cell_id
    )
    dat_used <- dat_used[res_nn$cell_id, ]
    
    stopifnot(all(rownames(dat_used) == res_nn$cell_id))
    
    
    # idx
    nn_idx <- res_nn$nn$nn.idx %>% t %>% c()
    group_idx <- rep(res_nn$cell_id, each = k_nn)
    
    # sum
    print('sum of nn cells')
    dat_merged <- dat_used %>% apply(2, function(x) tapply(x[nn_idx], group_idx, sum) )
    
    print('done.')
    
    # return
    dat_merged
  })



## kmeans cluster -------------
res_kmeans_used_ls <- 
  dat_zone_sum_ls %>%
  map(function(dat){
    set.seed(123)
    kmeans(dat, centers = 4) 
  })


p_heat_ls <-
  res_kmeans_used_ls %>%
  map('centers') %>%
  map(function(x){
    pheatmap::pheatmap(
      t(x),
      scale = "row",
      cluster_cols = F,
      # cluster_rows = F,
      border_color = NA,
      silent = T
    )
  })


pdf(
  file = file.path(out_dir_v6, "1 kmeans cluster IZ and BZ - heatmap.pdf"),
  height = 5,
  width  = 5
)
p_heat_ls %>% map('gtable') %>% map(cowplot::plot_grid)
dev.off()


df_kmeans_cluster <-
  res_kmeans_used_ls %>%
  map('cluster') %>%
  map(enframe, 'cell_id', 'cluster') %>%
  imap(function(dat, name){
    # cluster annotation by heatmap
    if(name == 'A02994A6'){
      dat <- dat %>%
        mutate(
          cluster = str_c('cluster', cluster),
          # cluster = fct_recode(
          #   cluster,
          #   "RZ_c1" = 'cluster1',
          #   "RZ_c2" = 'cluster2',
          #   "BZ/IZ_c3" = 'cluster3',
          #   "BZ/IZ_c4" = 'cluster4'
          # )
        )
    }else{
      dat <- dat %>%
        mutate(
          cluster = str_c('cluster', cluster),
          # cluster = fct_recode(
          #   cluster,
          #   "BZ1_c1" = 'cluster1',
          #   "IZ_c2" = 'cluster2',
          #   "RZ_c3" = 'cluster3',
          #   "BZ2_c4" = 'cluster4'
          # )
        )
    }
  }) %>%
  bind_rows() %>%
  column_to_rownames("cell_id")



p_sp_iz_used_ls <-
  cell_coords_ls2 %>%
  imap(function(coord, name){
    coord %>%
      left_join(
        rownames_to_column(df_kmeans_cluster, "cell_id"),
        by = "cell_id"
      ) %>%

      slice_sample(n = 50000) %>%

      ggplot(aes(x = x, y = -y, color = as.character(cluster))) +
      geom_point(size = 0.3) +
      coord_equal() +
      theme_minimal() +
      # scale_color_manual(values = c(
      #   '1' = "red",
      #   '2' = 'blue'
      # )) +
      guides(color = guide_legend(override.aes = list(size = 2))) +
      labs(
        color = 'cluster',
        title = name
      )
  })


pdf(
  file = file.path(out_dir_v6, "1 kmeans cluster IZ and BZ - spatial plot.pdf"),
  height = 20,
  width  = 20
)
p_sp_iz_used_ls %>% walk(print)
dev.off()



## vlnplot of IZ and BZ ------------
all_data_score_zone <-
  all_data_score_zone %>%
  AddMetaData(df_kmeans_cluster)

p_vln_zone_ls <-
  path_ls %>% 
  names() %>%
  map(function(x){
    all_data_score_zone %>%
      VlnPlot(
        features = x,
        group.by = "project",
        split.by = "cluster",
        pt.size = 0
      )  +
      geom_boxplot(width = 0.15, outliers = F, position = position_dodge(width = 0.9)) +
      ggtitle(x)
    
  })



pdf(
  file = file.path(out_dir_v5, "1 kmeans cluster IZ and BZ - vln plot.pdf"),
  height = 6,
  width  = 12
)
p_vln_zone_ls %>% walk(print)
dev.off()

## 20260507 kmeans cluster by ROI -------------
roi_ls <- 
  all_data[[]] %>% 
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  nest() %>%
  deframe() %>%
  map(~ with(.x, split(cell_id, roi)))
  


res_kmeans_used_ls_ls <- 
  roi_ls %>%
  imap(function(roi, proj){
    roi %>%
      map(function(r){
        dat <- dat_zone_sum_ls[[proj]]
        dat <- dat[rownames(dat) %in% r, ]
        set.seed(123)
        kmeans(dat, centers = 4) 
      })
  }) 



p_heat_ls_ls <-
  res_kmeans_used_ls_ls %>%
  imap(function(res_kmeans_used_ls, proj){
    res_kmeans_used_ls %>%
      map('centers') %>%
      imap(function(x, roi_name){
        pheatmap::pheatmap(
          t(x),
          scale = "row",
          cluster_cols = F,
          # cluster_rows = F,
          border_color = NA,
          silent = T,
          main = str_glue("{proj} - {roi_name}")
        )
      })
  })
   


pdf(
  file = file.path(out_dir_v6, "1 kmeans cluster IZ and BZ - heatmap - v2.pdf"),
  height = 5,
  width  = 5
)
p_heat_ls_ls %>% flatten() %>% map('gtable') %>% map(cowplot::plot_grid)
dev.off()


# spatial plot 
df_kmeans_cluster2 <-
  res_kmeans_used_ls_ls %>%
  purrr::flatten() %>%
  map('cluster') %>%
  map(enframe, 'cell_id', 'cluster') %>%
  imap(function(dat, name){
    dat %>%
      mutate(cluster = str_c('cluster', cluster))
  }) %>%
  bind_rows() %>%
  column_to_rownames("cell_id")
          
          
p_sp_iz_used_ls2 <-
  cell_coords_ls2 %>%
  imap(function(coord, name){
    coord %>%
      left_join(
        rownames_to_column(df_kmeans_cluster2, "cell_id"), 
        by = "cell_id"
      ) %>%
      mutate(cluster = as.character(cluster)) %>% 
      # slice_sample(n = 50000) %>%
      
      ggplot(aes(x = x, y = -y, color = cluster)) +
      # geom_point(size = 0.3) +
      geom_tile(aes(fill = cluster), height = 10, width = 10) +
      coord_equal() +
      theme_minimal() +
      # scale_color_manual(values = c(
      #   '1' = "red",
      #   '2' = 'blue'
      # )) +
      # guides(color = guide_legend(override.aes = list(size = 2))) +
      labs(
        color = 'cluster',
        title = name
      )
  })


pdf(
  file = file.path(out_dir_v6, "1 kmeans cluster IZ and BZ - spatial plot - v2.pdf"),
  height = 20,
  width  = 20
)
p_sp_iz_used_ls2 %>% walk(print)
dev.off()

