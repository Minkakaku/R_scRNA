library(Seurat)
library(tidyverse)
library(CellChat)

source("spatial analysis.R")

out_dir_260413 <- "2026413 dat_analysis"
dir.create(out_dir_260413)

# read data  -------
# sp
all_data_full <- readRDS("all_data_full_v6.rds")

# sc
all_data_sc <- readRDS("data/gastric_sc_st_datas_20251117/scRNA/adata_raw.rds")
all_data_sc[['RNA']] <- as(all_data_sc[['RNA']], "Assay5")
all_data_sc <- NormalizeData(all_data_sc)
all_data_sc$celltype <- all_data_sc$major_celltype
all_data_sc$project <- all_data_sc$sample_id



# all_data_full[[]] %>%
#   group_by(CN_cellcharter, major_celltype) %>%
#   summarise(count = n()) %>%
#   mutate(freq = count / sum(count)) %>%
#   
#   ggplot(aes(x = CN_cellcharter, y = freq, fill = major_celltype)) +
#   geom_col() +
#   ggsci::scale_fill_igv()
# 
# all_data_full[[]] %>%
#   group_by(CN_cellcharter, major_celltype) %>%
#   summarise(count = n()) %>%
#   mutate(freq = count / sum(count)) %>%
#   
#   dplyr::filter(freq > 0.01) %>% 
#   arrange(CN_cellcharter) %>%
#   summarise(count = n())
# 
# 
# c(
#   'sub_celltype',
#   'sub_celltype2',
#   'sub_celltype3'
# ) %>%
#   map(function(x){
#     all_data_full[[]] %>%
#       group_by(CN_cellcharter, .data[[x]]) %>%
#       summarise(count = n()) %>%
#       mutate(freq = count / sum(count)) %>%
#       
#       ggplot(aes(x = CN_cellcharter, y = freq, fill = .data[[x]])) +
#       geom_col() +
#       ggsci::scale_fill_igv() +
#       ggtitle(x)
#   })


#  cn celltypes ---------
cn_ct_ls <-
  all_data_full[[]] %>%
  group_by(CN_cellcharter, sub_celltype2) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  dplyr::filter(CN_cellcharter %in% c("xCN1", "xCN6", "xCN7")) %>%
  dplyr::filter(freq > 0.01) %>%
  group_by(CN_cellcharter) %>%
  nest() %>%
  deframe() %>%
  map("sub_celltype2") %>%
  map(as.character) %>%
  map(sort)


all_data_sc$sub_celltype2 %>%  unique %>% as.character() %>% sort


cn_ct_to_sc_ls <-
  cn_ct_ls %>%
  map(function(x){
    # map to single cell celltype
    x %>%
      fct_recode(
        caf = "CAF", 
        CD4T = "CD4T",
        CD8T = "CD8T",
        enterocytes = "enterocytes", 
        goblet_cells = "goblet_cells", 
        Macrophage = "Macrophage", 
        Monocyte = "Monocytes", 
        NK = "NK",
        Pericyte = "Pericyte",
        tumor_cells = "tumor_cells", 
        normal_epi = "pmc",
        normal_epi = "gmc", 
        normal_epi = "enteroendocrine_cells"
        # Neutrophils 
      ) %>%
      as.character() %>%
      setdiff('Neutrophils')
  })


# cellcaht ---------
all_data_meta_all     <- all_data_sc[[]] 
all_data_norm_all     <- LayerData(JoinLayers(all_data_sc), 'data') 
stopifnot(all(rownames(all_data_meta_all) == colnames(all_data_norm_all)))

# only used 50,000 cells
set.seed(12345)
subset_cell_idx <- 
  all_data_meta_all %>%
  rownames_to_column("cell_id") %>%
  slice_sample(prop = 1) #  signif(50000/nrow(all_data_meta_all)))
all_data_meta <- all_data_meta_all[subset_cell_idx$cell_id, ]
all_data_norm <- all_data_norm_all[, subset_cell_idx$cell_id]


all_grp_ls <-
  cn_ct_to_sc_ls %>%
  map(function(x){
    all_data_meta %>%
      rownames_to_column("cell_id") %>%
      dplyr::filter(celltype %in% x) %>%
      dplyr::pull("cell_id")
  })



# using CellChat db --------
all_cellchat_ls <-
  all_grp_ls %>%
  map(function(cells){ # browser()
    stopifnot(all(c('celltype', 'project') %in% colnames(all_data_sc[[]])))
    
    all_data_meta_all$celltype <- as.character(all_data_meta_all$celltype)
    all_data_meta$celltype     <- as.character(all_data_meta$celltype)
    
    meta <- all_data_meta[cells, ] 
    dat  <- all_data_norm[, cells]
    
    # table
    cells <- all_data_meta_all$celltype  %>% unique()
    if(!all(cells %in% meta$celltype )){
      missing_cells <- cells %>% setdiff(unique(meta$celltype ))
      missing_meta <-
        all_data_meta_all %>%
        rownames_to_column('cell_id') %>%
        dplyr::filter(celltype  %in% missing_cells) %>%
        group_by(celltype ) %>%
        slice_head(n = 1) %>%
        column_to_rownames("cell_id")
      meta <- rbind(meta, missing_meta)
      print(missing_meta)
      dat <- cbind(dat, all_data_norm_all[, rownames(missing_meta)])
    }
    
    meta$samples <- factor(meta$project)
    
    cellchat <- createCellChat(object = dat, meta = meta, group.by = "celltype")
    
    # cellchat db
    CellChatDB.use <- 
      CellChatDB.human %>% 
      subsetDB(
        search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"), 
        key    = "annotation"
      )
    
    cellchat@DB <- CellChatDB.use
    
    
    # preprocess
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    
    # optinal - project PPI data
    # cellchat = projectData(cellchat, PPI.human)
    cellchat <- smoothData(cellchat, adj = PPI.human)
    
    # Compute the communication probability
    cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = TRUE, population.size = TRUE)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    
    # Infer the cell-cell communication at a signaling pathway level
    cellchat <- computeCommunProbPathway(cellchat)
    
    # aggregated cell-cell communication network
    cellchat <- aggregateNet(cellchat)
    # cellchat <- netAnalysis_computeCentrality(cellchat)
    
    cellchat
  })


# saveRDS(all_cellchat_ls, file.path(out_dir_260413, "all_cellchat_ls2.rds"))



cellchat <- mergeCellChat(all_cellchat_ls, add.names = names(all_cellchat_ls), cell.prefix = TRUE)


cellchat %>% subsetCommunication()

all_weigth_ls <-
  all_cellchat_ls %>%
  imap(function(x, name){
    x@netP$prob %>%  
      apply(3, sum) %>%
      data.frame(weight = .) %>%
      setNames(name) %>%
      rownames_to_column("Path")
  }) %>%
  purrr::reduce(full_join, by = "Path")
  
all_count_ls <-
  all_cellchat_ls %>%
  imap(function(x, name){
    x %>%
      subsetCommunication() %>%
      group_by(pathway_name) %>%
      summarise(count = n()) %>%
      dplyr::rename(Path = pathway_name, {{name}} := count)
  }) %>%
  purrr::reduce(full_join, by = "Path")


# plot
p_cn_ct_ls <-
  all_cellchat_ls %>%
  names() %>%
  setNames(., .) %>%
  map(function(x){
    
    weight <- all_weigth_ls[c("Path", x)] %>% dplyr::rename(weight = 2)
    count  <- all_count_ls[c("Path", x)] %>% dplyr::rename(count = 2)
    
    plot_dat <- weight %>% left_join(count, by = "Path")
    plot_dat$weight[is.na(plot_dat$weight)] <- 0
    plot_dat$count[is.na(plot_dat$count)] <- 0
    
    plot_dat %>%
      ggplot(aes(x = reorder(Path, -count), y = count, fill = weight)) +
      geom_col() +
      theme_classic() +
      RotatedAxis() +
      theme(
        
      ) +
      labs(
        x = NULL,
        y = "Number of interations",
        fill = "Interaction strength",
        title = str_glue("scRNA-seq {x} CellChat")
      )
  })


# export
pdf(
  file = file.path(out_dir_260413, "CN cellchat in single cell.pdf"),
  height = 6,
  width = 14
)
p_cn_ct_ls %>% walk(print)
dev.off()


