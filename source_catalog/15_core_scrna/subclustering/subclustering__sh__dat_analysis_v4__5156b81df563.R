# conda activate rpy
library(Seurat)
library(tidyverse)

# all_data_sc <- readRDS('all_data_sc_doublet.rds')
all_data_only_prior <- readRDS("/data/NAS02/sh/Projects/FWSC2024-192/00.Mission/all_data_only_prior.rds")


out_dir <- "/data/NAS02/sh/Projects/FWSC2024-192/00.Mission/result_24.12.18" 
dir.create(out_dir)


sample_grp_info <- read_tsv("/data/NAS02/sh/Projects/FWSC2024-192/00.Mission/sample_info.tsv")

# read gene list
gene_ls <-
  readxl::read_excel("/data/NAS02/sh/Projects/FWSC2024-192/00.Mission/41591_2020_944_MOESM2_ESM.xlsx", sheet = "Supplementary Table 25", skip = 4) %>%
  map(function(x){
    x[!is.na(x)]
  })
names(gene_ls) <- 
  names(gene_ls) %>%
  str_replace("\\r\\n", "")


# 1. add more cytokine analysis ----------------------
p_dotplot_only_prior_more_cytekine_ls <-  
  gene_ls %>%
  imap(function(genes, name){ # browser()
    DotPlot(
      all_data_only_prior,
      features = genes,
      # cols = c("white", "lightgrey", "blue"),
      cols = c("RdYlBu"),
      group.by = "annotation",
      split.by = "group"
    ) +
    geom_hline(
      yintercept = all_data_only_prior$annotation %>% unique() %>% as.character() %>% seq_along() %>% {3 * .[-length(.)] + 0.5},
      linetype = 2
    ) +
    scale_color_gradient2(low = "white", mid = "lightgrey", high = "blue") +
    RotatedAxis() +
    ggtitle(name)
  })

  
# export
p_dotplot_only_prior_more_cytekine_ls %>% 
  iwalk(function(p, name){ 
    name <- name %>% str_remove("\\+")
    len_gene <- p_dotplot_only_prior_more_cytekine_ls$ISG$data$features.plot %>% n_distinct()
    pdf(
      file.path(out_dir, str_glue("1. dotplot of cytokine - group - {name}.pdf")),
      height = 10,
      width  = 4 + len_gene / 6
    )
    print(p)
    dev.off()
  })



  
# dot_plot_dat <-
#   all_data_only_prior %>%
#   FetchData(
#     vars = gene_ls[[1]],
#     layer = "data"
#   ) %>%
#   cbind(
#     all_data_only_prior[[c("annotation", "group")]]
#   ) %>%
#   group_by(annotation, group) %>%
#   summarise(across(everything(), list(
#     mean.expr = function(x) mean(expm1(x), na.rm = T),
#     pct.exp   = function(x) sum(x > 0)/length(x)
#   ))) %>%
#   ungroup() %>% 
#   pivot_longer(
#     -c(annotation, group),
#     names_to  = c('name'),
#     values_to = "value"
#   ) %>%
#   separate(col = name, into = c("gene", "var"), sep = "_") %>%
#   pivot_wider(
#     names_from = 'var',
#     values_from = "value"
#   ) %>%
#   group_by(gene) %>%
#   mutate(mean.expr = as.numeric(scale(log1p(mean.expr))))
# 
# 
# dot_plot_dat %>%
#   ggplot() +
#   geom_point(aes(x = mean.expr, y = pct.exp, color = group)) +
#   facet_wrap(~ group)

# 1.2 add more cytokine - split by project ---------

p_dotplot_only_prior_more_cytekine_byProject_ls <-  
  gene_ls %>%
  imap(function(genes, name){
    # porject 
    all_data_only_prior$project %>%
      unique() %>%
      set_names(., .) %>%
    map(function(proj){
      Idents(all_data_only_prior) <- 'project'
      all_data_only_prior %>%
        subset(idents = proj) %>%
        DotPlot(
          features = genes,
          # cols = c("white", "lightgrey", "blue"),
          cols = c("RdGy"),
          group.by = "annotation",
          split.by = "group"
        ) +
        geom_hline(
          yintercept = all_data_only_prior$annotation %>% unique() %>% as.character() %>% seq_along() %>% {3 * .[-length(.)] + 0.5},
          linetype = 2
        ) +
        scale_color_gradient2(low = "white", mid = "lightgrey", high = "blue") +
        RotatedAxis() +
        ggtitle(str_glue("{proj} - {name}"))
    })
    

  })


# export
p_dotplot_only_prior_more_cytekine_byProject_ls %>% 
  iwalk(function(p_ls, name){
    name <- name %>% str_remove("\\+")
    len_gene <- p_dotplot_only_prior_more_cytekine_ls$ISG$data$features.plot %>% n_distinct()
    pdf(
      file.path(out_dir, str_glue("1. dotplot of cytokine - project - {name}.pdf")),
      height = 10,
      width  = 4 + len_gene / 6
    )
    print(p_ls)
    dev.off()
  })


# 2. IFNG CXCL9 IL10 IFNGR1  IFNGR2  个案呈现 --------
p_dotplot_only_prior_5genes_byProject_ls <-
  all_data_only_prior$project %>%
  unique() %>%
  set_names(., .) %>%
  map(function(proj){
    Idents(all_data_only_prior) <- 'project'
    all_data_only_prior %>%
      subset(idents = proj) %>%
      DotPlot(
        features = c('IFNG', 'CXCL9', 'IL10', 'IFNGR1', 'IFNGR2'),
        # cols = c("white", "lightgrey", "blue"),
        cols = c("RdGy"),
        group.by = "annotation",
        split.by = "group"
      ) +
      geom_hline(
        yintercept = all_data_only_prior$annotation %>% unique() %>% as.character() %>% seq_along() %>% {3 * .[-length(.)] + 0.5},
        linetype = 2
      ) +
      scale_color_gradient2(low = "white", mid = "lightgrey", high = "blue") +
      RotatedAxis() +
      ggtitle(proj)
  })

# export
pdf(
  file.path(out_dir, str_glue("2. dotplot of gene CXCL9 IFNG IL10 IFNGR - split by project.pdf")),
  height = 6,
  width  = 8
)
print(p_dotplot_only_prior_5genes_byProject_ls)
dev.off()



# 3. IFN打分  ---------

ifn_geneset <- c(
  'CXCL10', 'GBP1', 'CXCL9', 'GBP5', 'GBP4', 'FAM26F', 
  'STAT1', 'IRF1', 'SLAMF7', 'WARS', 'NFKBIA', 'C3'
)

undebug(AddModuleScore)

all_data_only_prior2 <- 
  all_data_only_prior %>%
  AddModuleScore(
    list(IFN = ifn_geneset), 
    nbin = 24,
    ctrl = 100,
    name = 'IFN'
  ) 


all_data_only_prior2$IFN <- all_data_only_prior2$IFN1

p_vln_ifn_pathway <-
  all_data_only_prior2 %>% {
    .$group <- factor(.$group, levels = c("NC", "NLAHS", "LAHS"))
    .
  } %>%
  VlnPlot(
    features = "IFN", 
    group.by = "annotation",
    split.by = "group",
    pt.size = 0,
    log = TRUE
  ) + 
  geom_vline(
    xintercept = all_data_only_prior2$annotation %>% unique() %>% as.character() %>% seq_along() %>% {1 * .[-length(.)] + 0.5},
    linetype = 2
  ) +
  ggsci::scale_fill_aaas()

p_vln_ifn_pathway2 <-
  all_data_only_prior2 %>% {
    .$group <- factor(.$group, levels = c("NC", "NLAHS", "LAHS"))
    .$project <- str_c(.$group, .$project, sep = "_")
    .$project <- reorder(.$project, as.numeric(.$group))
    .
  } %>%
  VlnPlot(
    features = "IFN", 
    group.by = "annotation",
    split.by = "project",
    pt.size = 0,
    log = TRUE
  )  +
  geom_vline(
    xintercept = all_data_only_prior2$annotation %>% unique() %>% as.character() %>% seq_along() %>% {1 * .[-length(.)] + 0.5},
    linetype = 2
  )  +
  ggsci::scale_fill_igv()

# export
pdf(
  file.path(out_dir, str_glue("3. Vlnplot of IFN pathway score - split by group.pdf")),
  height = 6,
  width  = 10
)
print(p_vln_ifn_pathway)
dev.off()

pdf(
  file.path(out_dir, str_glue("3. Vlnplot of IFN pathway score - split by project.pdf")),
  height = 6,
  width  = 15
)
print(p_vln_ifn_pathway2)
dev.off()



