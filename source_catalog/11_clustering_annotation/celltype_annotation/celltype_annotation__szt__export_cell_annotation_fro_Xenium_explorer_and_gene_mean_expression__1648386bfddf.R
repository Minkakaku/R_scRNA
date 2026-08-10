library(tidyverse)
library(Seurat)

all_data <- readRDS("report_all_samples/output/data/adata_celltype.rds")

all_data_epi <- readRDS("all_data_epi.rds")
all_data_cn  <- readRDS("20251022 epi nn analysis/all_data_cn.rds")


out_dir_temp <- "20251031"
dir.create(out_dir_temp)

# 1. gene expr of ------------
h_geneset <- msigdbr::msigdbr(category = "H")

h_geneset_ls <- 
  h_geneset %>%
  group_by(gs_name) %>%
  nest() %>%
  deframe() %>%
  map('gene_symbol')

avg_exprt_h_ls <-
  h_geneset_ls %>%
  map(function(x){
    dat <- all_data %>%
      FetchData(
        vars = c('group', x)
      )
    apply(
      dat[-1],
      2, 
      function(x)tapply(x, dat$group, mean)
    )
  })

names(avg_exprt_h_ls) <-
  names(avg_exprt_h_ls) %>% 
  str_sub(start = 10) %>% 
  str_to_lower() %>% 
  str_sub(end = 31)

avg_exprt_h_ls %>%
  map(as.data.frame) %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir_temp, "mean expression of each gene in hallmark geneset.xlsx"),
    asTable = T,
    withFilter = F
  )
  
# 2. export cell annotation -------
out_dir_temp_nn <- file.path(out_dir_temp, "nn cluster2")
dir.create(out_dir_temp_nn)



all_data[[c('fov', 'cell_id', 'annotation2')]] %>%
  dplyr::rename(group = annotation2) %>%
  group_by(fov) %>%
  nest() %>%
  deframe() %>%
  iwalk(function(dat, name){
    write.csv(
      dat, 
      file = file.path(out_dir_temp_nn, str_glue("{name}.csv")),
      row.names = FALSE
    )
  })


out_dir_temp_gpc3 <- file.path(out_dir_temp, "gpc3 cell2")
dir.create(out_dir_temp_gpc3)

unique( all_data$annotation[!is.na(all_data$AUC_GPC3_Group)] )

all_data[[c('fov', 'cell_id', 'annotation', 'AUC_GPC3_Group')]] %>%
  dplyr::mutate(
    annotation2 = ifelse(is.na(AUC_GPC3_Group), annotation, AUC_GPC3_Group)
  ) %>%
  dplyr::select(fov, cell_id, group = annotation2) %>%
  group_by(fov) %>%
  nest() %>%
  deframe() %>%
  iwalk(function(dat, name){
    write.csv(
      dat, 
      file = file.path(out_dir_temp_gpc3, str_glue("{name}.csv")),
      row.names = FALSE
    )
  })


# 20251105 WNT gene expression ---------------
wnt_path <- KEGGREST::keggGet('hsa04310')

wnt_gene <- 
  wnt_path[[1]]$GENE %>% matrix(ncol = 2, byrow = T) %>%
  as.data.frame() %>%
  setNames(c("gene_id", "gene_name"))

wnt_gene %>%
  openxlsx::write.xlsx(
    file = "gene of WNT pathway.xlsx"
  )

gene_db <- 
  clusterProfiler::bitr(
    wnt_gene$gene_id,
    fromType = "ENTREZID",
    toType = "SYMBOL",
    OrgDb = "org.Hs.eg.db"
  )

wnt_gene_full <- 
  wnt_gene %>%
  left_join(gene_db, by = c(gene_id = 'ENTREZID'))

wnt_gene_full2 <-
  wnt_gene_full %>%
  dplyr::filter(SYMBOL %in% rownames(all_data))

avg_exprt_wnt_ls <-
  list(WNT = wnt_gene_full2$SYMBOL) %>%
  map(function(x){
    dat <- all_data %>%
      FetchData(
        vars = c('group', x)
      )
    apply(
      dat[-1],
      2, 
      function(x)tapply(x, dat$group, mean)
    )
  })

wnt_gene_full2_add_expr <-
  wnt_gene_full2 %>%
  left_join(
    rownames_to_column(as.data.frame(t(avg_exprt_wnt_ls$WNT)), "SYMBOL"),
    by = "SYMBOL"
  )

openxlsx::write.xlsx(
  wnt_gene_full2_add_expr,
  file = file.path(out_dir_temp, "mean expression of each gene in KEGG WNT pathway.xlsx"),
  asTable = T,
  withFilter = F
)



