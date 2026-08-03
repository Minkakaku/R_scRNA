################################
######### using sceasy #########
################################

# terminal
cd /mnt/NexData_file/home/szt/others/20250908_heart/data_raw
conda activate tools

R

# R env
library(reticulate)
use_python("/home/szt/miniforge3/envs/tools/bin/python")

`%>%` <- magrittr::`%>%`

all_files <- list.files(pattern = "h5ad$", full.names = TRUE)
names(all_files) <- all_files %>% basename() %>% stringr::str_replace("h5ad", "rds")

# others bin
all_files[! stringr::str_detect(all_files, "cellbin")]  %>%
  purrr::iwalk(function(adata, rds){
    print(rds)
    if(file.exists(rds)) return()
    adata_raw <-
      sceasy::convertFormat(
        adata,
        from = "anndata",
        to = "seurat",
        main_layer = "data"
      )
    saveRDS(adata_raw, file = rds)
  })



# cellbin
all_files[stringr::str_detect(all_files, "cellbin")] %>%
  purrr::iwalk(function(h5ad, rds){
    print(rds)
    if(file.exists(rds)) return()
    
    sc <- import('scanpy')
    adata <- sc$read(h5ad)
    
    saveRDS(adata$obsm['cell_border'], stringr::str_replace(rds, "rds", "cellborder.rds"))
    adata$obsm['cell_border'] <- NULL
    
    adata <- adata$raw$to_adata()
    
    # var; obs
    var_df <- adata$var
    obs_df <- adata$obs
    
    # count 
    count <- adata$X %>% as.matrix() %>% t
    
    # count unique 
    dup_gene <- var_df$real_gene_name %>% .[duplicated(.)]

    count_uniq <- NULL
    for(gene in dup_gene){
      gene_idx <- var_df$real_gene_name %in% gene
      count_uniq <- rbind(count_uniq, colSums(count[gene_idx, ,drop=FALSE]))
    }
    rownames(count_uniq) <- dup_gene
    
    # combine count
    dup_idx <- var_df$real_gene_name %in% dup_gene
    non_dup_count <- count[!dup_idx, ]
    rownames(non_dup_count) <- var_df$real_gene_name[!dup_idx] 
    
    count_final <- rbind(non_dup_count, count_uniq)
    colnames(count_final) <- rownames(obs_df)
    
    # seurat object
    adata_raw <- SeuratObject::CreateSeuratObject(
      counts = count_final,
      meta.data = obs_df
    )
    # sc$write(adata = adata, filename = stringr::str_replace(h5ad, ".h5ad", "2.h5ad"))
    # adata_raw <-
    #   sceasy::convertFormat(
    #     stringr::str_replace(h5ad, ".h5ad", "2.h5ad"),
    #     from = "anndata",
    #     to = "seurat",
    #     main_layer = "data"
    #   )
    saveRDS(adata_raw, file = rds)
  })


# bin20 combined duplicated gene 
all_files[stringr::str_detect(all_files, "bin20")]  %>%
  purrr::iwalk(function(h5ad, rds){
    
    rds <- stringr::str_replace(rds, "rds", "dropdup.rds")
    print(rds)
    
    if(file.exists(rds)) return()
    
    sc <- import('scanpy')
    adata <- sc$read(h5ad)
    adata <- adata$raw$to_adata()
    
    # var; obs
    var_df <- adata$var
    obs_df <- adata$obs
    
    # count 
    count <- adata$X %>% as.matrix() %>% t
    
    # count unique 
    dup_gene <- var_df$real_gene_name %>% .[duplicated(.)]
    
    count_uniq <- NULL
    for(gene in dup_gene){
      gene_idx <- var_df$real_gene_name %in% gene
      count_uniq <- rbind(count_uniq, colSums(count[gene_idx, ,drop=FALSE]))
    }
    rownames(count_uniq) <- dup_gene
    
    # combine count
    dup_idx <- var_df$real_gene_name %in% dup_gene
    non_dup_count <- count[!dup_idx, ]
    rownames(non_dup_count) <- var_df$real_gene_name[!dup_idx] 
    
    count_final <- rbind(non_dup_count, count_uniq)
    colnames(count_final) <- rownames(obs_df)
    
    # seurat object
    adata_raw <- SeuratObject::CreateSeuratObject(
      counts = count_final,
      meta.data = obs_df
    )
    
    saveRDS(adata_raw, file = rds)
  })



# adata_celltype <-
#   sceasy::convertFormat(
#     "adata_celltype.h5ad", 
#     from = "anndata", 
#     to = "seurat",
#     main_layer = "scale.data" # <----- 
#   )
# 
# adata_celltype_anno <-
#   sceasy::convertFormat(
#     "adata_celltype_anno.h5ad", 
#     from = "anndata", 
#     to = "seurat",
#     main_layer = "scale.data" # <----- 
#   )
# 
# adata_raw <-
#   sceasy::convertFormat(
#     "adata_raw.h5ad", 
#     from = "anndata", 
#     to = "seurat" 
#   )
# 
# saveRDS(adata_celltype, 'adata_celltype.rds')
# saveRDS(adata_celltype_anno, 'adata_celltype_anno.rds')
# saveRDS(adata_raw, 'adata_raw.rds')


######### export gene_info #########

library(reticulate)
use_python("/home/szt/miniforge3/envs/tools/bin/python")

scanpy <- import("scanpy")
`%>%` <- magrittr::`%>%`

all_files <- list.files('data_raw', pattern = "h5ad$", full.names = TRUE)
names(all_files) <- all_files %>% basename() %>% stringr::str_replace("h5ad", "rds")

adata <- scanpy$read("data_raw/A02994A6.cellbin_1.0.adjusted.h5ad")
adata$var %>%
  rownames_to_column('gene_id') %>%
  write.csv(
    file = "data_raw/gene annotation.csv",
    row.names = FALSE
  )

