library(tidyverse)
library(Seurat)
library(clusterProfiler)


#### flatDoubleList
flatDoubleList <- function(x, sep = "---"){
  stopifnot(is.vector(x))
  
  res_ls <- list()
  
  # layer 1
  name_1 <- names(x)
  id_1   <- seq_along(x)
  if(is.null(name_1)) name_1 <- id_1
  
  for(i in id_1){
    y <- x[[i]]
    stopifnot(is.vector(y))
    
    # layer 2
    name_2 <- names(y)
    id_2   <- seq_along(y)
    if(is.null(name_2)) name_2 <- id_2
    
    for(j in id_2){
      res  <- y[[j]]
      name <- paste0(name_1[i], sep, name_2[j])
      res_ls[[name]] <- res
    }
  }
  res_ls
}

##### nestSingleList
nestSingleList <- function(x, sep = "---"){
  stopifnot(is.vector(x))
  
  name <- names(x)
  id_ <- seq_along(name)
  stopifnot(all(str_detect(name, sep)))
  name_1 <- str_split(name, pattern = sep) %>% map_chr(~ .x[1])
  name_2 <- str_split(name, pattern = sep) %>% map(~ .x[2:length(.x)]) %>% map_chr(str_c, collapse = sep)
  
  res_ls <- list()
  for(i in id_){
    n1 <- name_1[i]
    n2 <- name_2[i]
    res <- x[[i]]
    res_ls[[n1]][[n2]] <- res
  }
  res_ls
}



# out dir
out_dir_analysis <- "LUSC_out_dir_20250920"
dir.create(out_dir_analysis)



all_data_subcluster2 <- readRDS('all_data_subcluster2.rds')
all_data_subcluster2 <- all_data_subcluster2 %>% JoinLayers()

# 最关键的就是R组和S组中不同细胞类型中基因表达的差异。
# 特别关注上皮细胞，T细胞和巨噬细胞；
# 差异基因还需要进行聚类分析，注明哪些通路差异特别显著。

# 1. cell props ------------

# major annotation
p_col_major_proj <-
  all_data_subcluster2[[]] %>%
  group_by(project, group, annotation_full) %>%
  summarise(count = n()) %>%
  group_by(project, group) %>%
  mutate(freq = count / sum(count)) %>%
  
  ggplot(aes(x = project, y = freq, fill = annotation_full)) +
  geom_col() +
  ggsci::scale_fill_igv()

p_col_major_grp <-
  all_data_subcluster2[[]] %>%
  group_by(project, group, annotation_full) %>%
  summarise(count = n()) %>%
  group_by(group) %>%
  mutate(freq = count / sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = annotation_full)) +
  geom_col() +
  ggsci::scale_fill_igv()

# export 
ggsave(
  p_col_major_proj,
  filename = file.path(out_dir_analysis, '1. cell props - major - project.pdf'),
  height = 6,
  width  = 6
)
ggsave(
  p_col_major_grp,
  filename = file.path(out_dir_analysis, '1. cell props - major - group.pdf'),
  height = 6,
  width  = 6
)

# sub annotation
p_col_sub_proj <-
  all_data_subcluster2[[]] %>%
  group_by(project, group, annotation2) %>%
  summarise(count = n()) %>%
  group_by(project, group) %>%
  mutate(freq = count / sum(count)) %>%
  
  ggplot(aes(x = project, y = freq, fill = annotation2)) +
  geom_col() +
  ggsci::scale_fill_igv()

p_col_sub_grp <-
  all_data_subcluster2[[]] %>%
  group_by(project, group, annotation2) %>%
  summarise(count = n()) %>%
  group_by(group) %>%
  mutate(freq = count / sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = annotation2)) +
  geom_col() +
  ggsci::scale_fill_igv()

# export 
ggsave(
  p_col_sub_proj,
  filename = file.path(out_dir_analysis, '1. cell props - sub - project.pdf'),
  height = 6,
  width  = 6
)
ggsave(
  p_col_sub_grp,
  filename = file.path(out_dir_analysis, '1. cell props - sub - group.pdf'),
  height = 6,
  width  = 6
)


# 2. major de analysis -------------
gene_description <- read_csv("all_data_analysis_V1_combined /data/XeniumPrimeHuman5Kpan_tissue_pathways_metadata.csv")

all_markers_major <- 
  FindAllMarkers(
    all_data_subcluster2,
    group.by = "group",
    logfc.threshold = 0,
    return.thresh = 1.1
  )
  
# add gene description
all_markers_major <- 
  all_markers_major %>%
  # rownames_to_column('gene') %>%
  left_join(
    gene_description,
    by = c(gene = "gene_name")
  )
  

# expot de table
write.csv(
  all_markers_major,
  file = file.path(out_dir_analysis, str_glue("2. de analysis table - all - group.csv")),
  row.names = FALSE
)


# 3. de analysis --------------
all_data_subcluster2 %>% head
all_data_subcluster2@assays$RNA@data[1:5, 1:50]


table(all_data_subcluster2$group, all_data_subcluster2$annotation2)

all_marker_ls <-
  all_data_subcluster2$annotation2 %>%
  unique() %>%
  setNames(., .) %>%
  
  map(function(x){
    all_data_subcluster2 %>% 
      {Idents(.) <- "annotation2"; .} %>% 
      FindMarkers(
        ident.1 = 'R',
        ident.2 = 'S',
        group.by = "group",
        subset.ident = x,
        logfc.threshold = 0
      )
  })


# add gene description
all_marker_ls <- 
  all_marker_ls %>%
  map(function(de){
    de %>%
      rownames_to_column('gene') %>%
      left_join(
        gene_description,
        by = c(gene = "gene_name")
      )
  })

# expot de table
all_marker_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(de, name){
    
    write.csv(
      de,
      file = file.path(out_dir_analysis, str_glue("3. de analysis table - {name} - group.csv")),
      row.names = FALSE
    )
  })


all_marker <- 
  all_marker_ls %>%
  # map(rownames_to_column, "gene") %>% 
  enframe("cluster") %>%
  unnest("value")



all_de_genes_for_enrichmenet <-
  all_marker %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map(function(x){
    up_gene <- 
      x %>%
      dplyr::filter(
        p_val_adj < 0.05,
        # pct.1 - pct.2 > 0.1,
        avg_log2FC > 1
      ) %>%
      .$gene
    
    down_gene <- 
      x %>% 
      dplyr::filter(
        p_val_adj < 0.05,
        # pct.1 - pct.2 > 0.1,
        avg_log2FC < -1
      ) %>%
      .$gene
    
    # return
    list(
      R = up_gene,
      S = down_gene
    )
  }) %>%
  flatDoubleList()


## multi volcano plot ------------ 
# all_markers_fib_ls$`S3 fibroblasts` %>% head
# p_val avg_log2FC pct.1 pct.2     p_val_adj cluster       gene
# MTRNR2L12  1.259519e-114  2.0699800 0.421 0.139 3.364427e-110    CD-F  MTRNR2L12
# FP236383.3 1.647682e-110  0.5967926 1.000 1.000 4.401288e-106    CD-F FP236383.3
# FP236383.5  4.720774e-79  0.9386205 0.750 0.558  1.261013e-74    CD-F FP236383.5
# FP236383.4  8.018641e-69  0.9385583 0.673 0.449  2.141939e-64    CD-F FP236383.4
# FP671120.7  3.710999e-63  0.6419110 0.977 0.949  9.912819e-59    CD-F FP671120.7
# SLCO5A1     1.146982e-47  0.9227888 0.460 0.258  3.063819e-43    CD-F    SLCO5A1

set.seed(1234)
rect_height <- 1
dodge_width <- 0.1
p_val_cutoff <- 0.05


# scatter
p_multi_volcano_ls <-
  list(
    all_marker_ls[1:15],
    all_marker_ls[16:length(all_marker_ls)]
  ) %>%
    map(function(dat_ls){
      df_scatter <-
        dat_ls %>%
        # map(rownames_to_column, "gene") %>%
        enframe('cluster') %>%
        unnest(value) %>% 
        dplyr::filter(!is.na(p_val_adj)) %>% 
        mutate(cluster = factor(cluster), x_base = as.numeric(cluster)) %>%
        slice_sample(prop = 1) %>%
        group_by(cluster) %>%
        
        mutate(
          x_delta = 1:n(),
          x_delta = scales::rescale(x_delta, to = c(dodge_width, 1 - dodge_width)),
          xx = x_base - 1 + x_delta,
          y = ifelse(avg_log2FC > 0, avg_log2FC + rect_height, avg_log2FC - rect_height)
        ) %>%
        mutate(
          group = ifelse(p_val_adj > p_val_cutoff, str_glue("p > {p_val_cutoff}"), str_glue("p < {p_val_cutoff}") ),
          group = factor(group, levels = c(str_glue("p > {p_val_cutoff}"), str_glue("p < {p_val_cutoff}") ))
        ) %>%
        dplyr::arrange(group)
      
      # df_label 
      df_label_up <- df_scatter %>%
        dplyr::filter(as.character(group) == str_glue("p < {p_val_cutoff}") ) %>%
        group_by(cluster) %>% 
        slice_max(avg_log2FC, n = 10)
      
      df_label_down <- df_scatter %>%
        dplyr::filter(as.character(group) == str_glue("p < {p_val_cutoff}") ) %>%
        group_by(cluster) %>% 
        slice_min(avg_log2FC, n = 10)
      
      
      
      # grp rect  
      df_grp_rect <-
        df_scatter %>%
        group_by(cluster) %>%
        summarise() %>%
        mutate(x_base = as.numeric(factor(cluster))) %>%
        mutate(xmin = x_base - 1, xmax = x_base) %>%
        mutate(ymin = -rect_height, ymax = rect_height)
      
      
      # background rect
      df_bg_rect <-
        df_scatter %>%
        group_by(cluster) %>%
        summarise(ymax = max(avg_log2FC) + rect_height, ymin = min(avg_log2FC) - rect_height) %>%
        mutate(x_base = as.numeric(factor(cluster))) %>%
        mutate(xmin = x_base - 1 + dodge_width, xmax = x_base - dodge_width)
      
      
      # text label
      df_text_grp <-
        df_scatter %>%
        group_by(cluster) %>%
        summarise() %>%
        mutate(x_base = as.numeric(factor(cluster))) %>%
        mutate(x = x_base - 1) %>%
        mutate(x = (x_base + x)/2)
      
      
      p_multi_volcano <-
        df_scatter %>%
        ggplot(aes(x = xx, y = y, fill = cluster)) +
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
        # geom_tile(
        #   data = ~ distinct(., cluster, .keep_all = TRUE),
        #   aes(x = x_base, y = 0),
        #   width = 1
        # ) +
        geom_text(
          data = df_text_grp,
          aes(x = x, y = 0, label = cluster)
        ) +
        
        ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
        ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
        scale_color_manual(values = structure(
          c(
            "gray70", "red"
          ),
          names = c(
            str_glue("p > {p_val_cutoff}"), str_glue("p < {p_val_cutoff}")
          )
        )) +
        ggsci::scale_fill_igv() +
        theme_minimal() +
        theme(
          panel.grid = element_blank(),
          axis.line.y = element_line(color = "black"),
          axis.text.x = element_blank()
        ) +
        labs(
          title = "R vs S",
          x = NULL,
          y = "Log2FC"
        )
      
      p_multi_volcano
    })

pdf(
  file  = file.path(out_dir_analysis, "3 de analysis - multi volcano plot.pdf"),
  width = 40,
  height = 8,
)
p_multi_volcano_ls %>% walk(print)
dev.off()

# ggsave(
#   p_multi_volcano,
#   filename = file.path(out_dir_analysis, "3 de analysis - multi volcano plot.pdf"),
#   width = 50,
#   height = 8,
#   limitsize = FALSE
# )

## heatmap -----
all_de_genes_ls <- all_de_genes_for_enrichmenet %>% nestSingleList()  %>% map(enframe, "group", "gene") %>% map(~ unnest(.x, gene) )

avg_expr <- AverageExpression(all_data_subcluster2, group.by = "sample") %>% .$RNA %>% as.matrix()

anno_df_sample_group <-
  all_data_subcluster2[[c("project", "group")]] %>% 
  distinct() %>% 
  mutate(project = str_replace(project, "_", "-")) %>%
  remove_rownames() %>%
  column_to_rownames("project") %>%
  .[colnames(avg_expr), ,drop=FALSE]

all_heatmap_ls <-
  all_de_genes_ls %>%
  imap(function(de, name){
    if(nrow(de) < 1) return()
    p <- 
      avg_expr[de$gene, ] %>%
      t %>% 
      scale %>% 
      t %>% 
      ComplexHeatmap::Heatmap(
        col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
        name = "expression",
        column_split = anno_df_sample_group$group,
        top_annotation = ComplexHeatmap::HeatmapAnnotation(df = anno_df_sample_group, col = list(group = c(R = "red", S = "green"))),
        row_split = de$group,
        row_gap =  unit(2, "mm"),
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2",
        show_row_names = FALSE,
        column_title = name
      )
  })

# epxort
pdf(
  file = file.path(out_dir_analysis, "3 de heatmap of sample.pdf"),
  height = 8,
  width  = 6
)
all_heatmap_ls %>% walk(print)
dev.off()



# 4. enrichemt analysis -----
all_xenium_gene <- 
  clusterProfiler::bitr(
    rownames(all_data_subcluster2),
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = "org.Hs.eg.db"
  )  %>%
  .$ENTREZID



dir_enrich <- file.path(out_dir_analysis, "enrichment analysis")
dir.create(dir_enrich)

all_de_genes_id_for_enrichmenet <-
  all_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::bitr(
      x,
      OrgDb = "org.Hs.eg.db",
      fromType = "SYMBOL",
      toType = "ENTREZID"
    ) %>%
      .$ENTREZID
  }) %>%
  discard(function(x) length(x) < 1)

gene_db <- 
  clusterProfiler::bitr(
    unique(unlist(all_de_genes_id_for_enrichmenet)),
    fromType = "ENTREZID",
    toType   = "SYMBOL",
    OrgDb    = "org.Hs.eg.db"
  ) 

## KEGG ------------
enrich_ls <-
  all_de_genes_id_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichKEGG(
      x,
      organism = "hsa", 
      keyType = "ncbi-geneid", 
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH",
      universe = all_xenium_gene,
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
    )
  })

enrich_ls <-
  enrich_ls %>%
  map(function(x){ # browser()
    xx <- x@result$geneID %>% str_split("/")
    
    # gene symbol
    xx2 <-
      xx %>%
      map(function(g){ # browser()
        gene_name_vec <- gene_db$SYMBOL[match(g, gene_db$ENTREZID)] 
        gene_name <- str_c(str_replace_na(gene_name_vec), collapse = ",")
        
        gene_loc_vec <- gene_description$location[match(gene_name_vec, gene_description$gene_name)]
        gene_loc  <- str_c(str_replace_na(gene_loc_vec), collapse = ",")
        
        list(
          gene_name,
          gene_loc
        )
      })
    
    x@result$geneSymbol <- xx2 %>% map(1) %>% unlist(use.names = FALSE)
    x@result$geneLocation <- xx2 %>% map(2) %>% unlist(use.names = FALSE)
    x
  })

# saveRDS(enrich_ls, file.path(dir_enrich, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    
    write.csv(
      x@result,
      file = file.path(dir_enrich, str_glue("4.1 KEGG - enriment table - {name}.csv")),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  imap(function(x, name){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    ) +
      ggtitle(name)
  })

# export dotplot
enrich_dot_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(dir_enrich, str_glue("4.1 dotplot - {name}.pdf")),
      height = 6,
      width = 8
    )
  })

## GOBP -------
gobp_enrich_ls <-
  all_de_genes_id_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichGO(
      x,
      OrgDb = 'org.Hs.eg.db',
      ont = "BP",
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      universe = all_xenium_gene,
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2
    )
  })

gobp_enrich_ls <-
  gobp_enrich_ls %>%
  map(function(x){ # browser()
    xx <- x@result$geneID %>% str_split("/")
    
    # gene symbol
    xx2 <-
      xx %>%
      map(function(g){ # browser()
        gene_name_vec <- gene_db$SYMBOL[match(g, gene_db$ENTREZID)] 
        gene_name <- str_c(str_replace_na(gene_name_vec), collapse = ",")
        
        gene_loc_vec <- gene_description$location[match(gene_name_vec, gene_description$gene_name)]
        gene_loc  <- str_c(str_replace_na(gene_loc_vec), collapse = ",")
        
        list(
          gene_name,
          gene_loc
        )
      })
    
    x@result$geneSymbol <- xx2 %>% map(1) %>% unlist(use.names = FALSE)
    x@result$geneLocation <- xx2 %>% map(2) %>% unlist(use.names = FALSE)
    x
  })
# saveRDS(gobp_enrich_ls, file.path(dir_enrich, 'gobp_enrich_ls.rds'))

# expot enrichment table
gobp_enrich_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    
    write.csv(
      x@result,
      file = file.path(dir_enrich, str_glue("4.2 GOBP - enriment table - {name}.csv")),
      row.names = FALSE
    )
  })

# dotplot
gobp_enrich_dot_ls <-
  gobp_enrich_ls %>%
  imap(function(x, name){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    ) +
      ggtitle(name)
  })

# export dotplot
gobp_enrich_dot_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(dir_enrich, str_glue("4.2 BOBP - dotplot - {name}.pdf")),
      height = 6,
      width = 8
    )
  })








