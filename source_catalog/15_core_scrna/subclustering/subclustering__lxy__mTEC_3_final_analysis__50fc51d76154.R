
mTEC_out_dir_3 <- "mTEC_subtype_3" 
dir.create(mTEC_out_dir_3)

#clean后如果run_umap时图形形状不理想，重新findcluster
mTEC <-
  FindNeighbors(
    mTEC ,
    reduction = 'pca',
    dims = 1:13,
    k.param = 20
  )


mTEC <-
  FindClusters(
    mTEC,
    resolution = c(0.1,0.2, 0.4, 0.6, 0.8)
  )

mTEC <-
  FindClusters(
    mTEC,
    resolution = 0.25
  )

# dimplot
pdf(file.path(mTEC_out_dir_3, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(mTEC, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(mTEC, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(mTEC, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(mTEC, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(mTEC, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'umap') + NoLegend()
# DimPlot(mTEC, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(mTEC, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(mTEC, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

saveRDS(mTEC, "mTEC_subtype_3/mTEC_clean_1.rds")


Idents(mTEC)<-'RNA_snn_res.0.4'
mTEC_markers <-
  FindAllMarkers(
    JoinLayers(mTEC),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(mTEC_markers, file = "mTEC_subtype_3/mTEC_all_markers_0.4.csv", row.names = TRUE)

mTEC_top5_markers <- 
  mTEC_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 5) %>%
  nest() %>%
  deframe() %>% 
  map(pull, 'gene') %>%
  imap(function(genes, name){
    re <- genes %>%
      map(function(gene){
        # browser()
        names(gene) = str_glue("cluster_{name}")
        gene
      })
    re
  }) %>%
  unlist() %>% 
  .[!duplicated(.)]

names(mTEC_top5_markers) <- sub('.*\\.','',names(mTEC_top5_markers))


# dotplot
# 如果已做其他细胞类型的亚群分析，这里就可以顺便画一下其他类型的细胞亚群基因，比如Fibro_major_marker_ls2，否则只用画T细胞
pdf(file.path(mTEC_out_dir_3, str_glue("dotplot dims{dims} res.0.4 major_marker.pdf")), width = 25, height = 6)
DotPlot(mTEC, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(mTEC, group.by = "RNA_snn_res.0.4", features = mTEC_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(mTEC, group.by = "RNA_snn_res.0.4", features = mTEC_top5_markers)+RotatedAxis()+theme_90angle
dev.off()

pdf(file.path(mTEC_out_dir_3, str_glue("dotplot dims{dims} T_major_marker.pdf")), width = 25, height = 6)
DotPlot(mTEC, group.by = "RNA_snn_res.0.4", features = T_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(mTEC, group.by = "RNA_snn_res.0.4", features = T_major_marker_ls_1)+RotatedAxis()+theme_90angle
DotPlot(mTEC, group.by = "RNA_snn_res.0.4", features = Plasma_major_marker_ls)+RotatedAxis()+theme_90angle
dev.off()

saveRDS(mTEC ,file = "mTEC_subtype_3/mTEC_cluster.rds")
mTEC <- readRDS("mTEC_subtype_3/mTEC_cluster.rds")

table(mTEC_markers$cluster)




cell_cluster <- c(
  '0'  = "TEC precursor", # 
  '1'  = "mTEC_1", #
  '2'  = "mTEC_2", #
  '3'  = "mTEC_3", # 
  '4'  = "cTEC_Naive_T", # 
  '5'  = "mTEC_4", # 
  '6'  = "mTEC_5", # 
  '7'  = "cTEC", # 
  '8'  = "cTEC_Naive_T", # 
  '9'  = "mTEC_6", # 
  '10'  = "mTEC_7"
)

Idents(mTEC) <- "RNA_snn_res.0.4"

mTEC <- RenameIdents(mTEC, cell_cluster)
mTEC$subtype <- Idents(mTEC)


pdf(file.path(mTEC_out_dir_3, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(mTEC, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(mTEC_out_dir_3, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 15)
DimPlot(mTEC, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()


# 更改原注释
# 1. 先把 all_data$T_annotation 复制到新列 mTEC_T_annotation
all_data$mTEC_3_Plasma_T_annotation <- as.character(all_data$mTEC_2_Plasma_T_annotation)

# 2. 创建映射：细胞名 → 精细T亚型
cell_subtype <- mTEC$subtype
names(cell_subtype) <- colnames(mTEC)

# 3. 只把 "T cell" 替换成精细亚型，其他细胞保持原样
all_data$mTEC_3_Plasma_T_annotation[colnames(mTEC)] <- as.character(cell_subtype)

# 查看结果（正确）
table(all_data$mTEC_3_Plasma_T_annotation)
table(all_data$mTEC_2_Plasma_T_annotation)
table(mTEC$subtype)
saveRDS(all_data, "mTEC_subtype_3/all_data_mTEC_3_Plasma_T_annotation.rds")


#只保留细胞亚群
keep <- c(
  "mTEC_1",
  "mTEC_2",
  "mTEC_3",
  "mTEC_4",
  "mTEC_5",
  "mTEC_6",
  "mTEC_7"
)

mTEC_clean <- subset(mTEC, subtype %in% keep)

mTEC_clean$subtype <- droplevels(mTEC_clean$subtype)
Idents(mTEC_clean) <- mTEC_clean$subtype

table(mTEC_clean$subtype)


mTEC <- mTEC_clean
mTEC <- RunUMAP(mTEC, dims = 1:13,n.neighbors = 50,return.model = TRUE)



pdf(file.path(mTEC_out_dir_3, str_glue("Dimplot dim{dims} mTEC annotation.pdf")), height = 6, width = 8)
DimPlot(mTEC, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()





pdf(file.path(mTEC_out_dir_3, str_glue("Dimplot dim{dims} mTEC annotation split.pdf")), height = 6, width = 15)
DimPlot(mTEC, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

saveRDS(mTEC, "mTEC_subtype_3/mTEC_clean.rds")



mTEC_markers <-
  FindAllMarkers(
    JoinLayers(mTEC),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(mTEC_markers, file = "mTEC_subtype_3/mTEC_all_markers_subtype.csv", row.names = TRUE)

mTEC_top5_markers <- 
  mTEC_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 5) %>%
  nest() %>%
  deframe() %>% 
  map(pull, 'gene') %>%
  imap(function(genes, name){
    re <- genes %>%
      map(function(gene){
        # browser()
        names(gene) = str_glue("cluster_{name}")
        gene
      })
    re
  }) %>%
  unlist() %>% 
  .[!duplicated(.)]

names(mTEC_top5_markers) <- sub('.*\\.','',names(mTEC_top5_markers))


# dotplot
# 如果已做其他细胞类型的亚群分析，这里就可以顺便画一下其他类型的细胞亚群基因，比如Fibro_major_marker_ls2，否则只用画T细胞
pdf(file.path(mTEC_out_dir_3, str_glue("dotplot dims{dims} major_subtype_marker.pdf")), width = 25, height = 6)
DotPlot(mTEC, group.by = "subtype", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(mTEC, group.by = "subtype", features = mTEC_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(mTEC, group.by = "subtype", features = mTEC_top5_markers)+RotatedAxis()+theme_90angle
dev.off()


## cell props barplot  ----------------------------------------------------------
mTEC_meta <- mTEC[[]]

mTEC_meta$sample <- factor(mTEC_meta$sample,levels=c(mTEC_meta$sample %>% unique()))

p_col_proj <-
  mTEC_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path("mTEC_subtype_3", "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  mTEC_meta %>%
  group_by(group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path(mTEC_out_dir_3, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)

#IFN基因集打分
IFN_geneset_ls <- list(
  IFN_a = c("ADAR","CASP1","CASP8","CCRL2","CD47","CNP","CSF1","CXCL10","CXCL11","EIF2AK2","HERC6","IFI35","IFI44L","IFIH1","IFIT2","IFIT3","IFITM1","IL15","IL4R","IL7","IRF1","IRF7","LAMP3","LPAR6","MOV10","MX1","OAS1","OGFR","PLSCR1","PNPT1","PROCR","RIPK2","RNF31","RSAD2","SAMD9","SELL","SP110","STAT2","TAP1","TDRD7","TRAFD1","TRIM21","TRIM25","TRIM5"),
  IFN_g = c("ADAR","APOL6","BANK1","CASP1","CASP3","CASP7","CASP8","CCL7","CD274","CD38","CD40","CD86","CDKN1A","CFB","CFH","CIITA","CMKLR1","CSF2RB","CXCL10","CXCL11","CXCL9","EIF2AK2","FAS","FCGR1A","FPR1","GCH1","GZMA","HERC6","HIF1A","ICAM1","IDO1","IFI35","IFI44L","IFIH1","IFIT1","IFIT2","IFIT3","IFNAR2","IL10RA","IL15","IL15RA","IL18BP","IL2RB","IL4R","IL6","IL7","IRF1","IRF4","IRF5","IRF7","IRF8","ITGB7","JAK2","KLRK1","LATS2","MVP","MX1","MYD88","NCOA3","NFKB1","NLRC5","NOD1","NUP93","OAS3","OGFR","P2RY14","PDE4B","PFKP","PIM1","PLA2G4A","PLSCR1","PML","PNP","PNPT1","PTGS2","PTPN1","PTPN2","PTPN6","RAPGEF6","RBCK1","RIPK1","RIPK2","RNF31","RSAD2","SAMHD1","SELP","SLAMF7","SOCS1","SOCS3","SP110","SPPL2A","SRI","STAT1","STAT2","STAT3","STAT4","TAP1","TDRD7","TNFSF10","TRAFD1","TRIM21","TRIM25","VAMP5","VCAM1","ZBP1")
)


mTEC_score <- JoinLayers(mTEC) %>%
  AddModuleScore(
    features = IFN_geneset_ls,
    name = names(IFN_geneset_ls)
  )

mTEC_score %>% head

# 按照subtype
# 5. 箱线图：IFN-α ------------------------------------------
p3 <- mTEC_score@meta.data %>%
  ggplot(aes(x = subtype, y = IFN_a1, fill = subtype)) +
  geom_boxplot(outlier.shape = NA) +          # 箱线图
  geom_jitter(
    width = 0.2,
    height = 0,
    alpha = 0.5,
    size = 1,
    color = "steelblue"
  ) +
  ggpubr::geom_pwc(tip.length = 0, hide.ns = FALSE, label = 'p.signif') +
  labs(x = "Subtype", y = "IFN-alpha Score") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_blank(),
    legend.position = "none"
  )

# 6. 箱线图：IFN-γ ------------------------------------------
p4 <- mTEC_score@meta.data %>%
  ggplot(aes(x = subtype, y = IFN_g2, fill = subtype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    width = 0.2,
    height = 0,
    alpha = 0.5,
    size = 1,
    color = "steelblue"
  ) +
  ggpubr::geom_pwc(tip.length = 0, hide.ns = FALSE, label = 'p.signif') +
  labs(x = "Subtype", y = "IFN-gamma Score") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_blank(),
    legend.position = "none"
  )

p_combined <- p3 + p4
ggsave(
  plot = p_combined,
  filename = paste0(mTEC_out_dir_3, "/IFN_score_subtype_Boxplot.pdf"),
  width = 8, height = 8
)



#补体基因集打分
Complement_geneset_ls <- list(
  Complement = c("ACTN2","ADAM9","ADRA2B","ANXA5","APOBEC3F","APOBEC3G","C9","CALM3","CASP1","CASP10","CASP3","CASP7","CASP9","CBLB","CD36","CD40LG","CD46","CD55","CD59","CDA","CDH13","CFB","CFH","COL4A2","CP","CR1","CR2","CSRP1","CTSC","CTSH","CTSL","DOCK9","DPP4","DUSP5","DUSP6","F2","F5","F7","F8","FCN1","FDX1","FYN","GATA3","GMFB","GNAI2","GP1BA","GRB2","GZMA","GZMB","GZMK","HNF4A","IL6","IRF1","IRF7","ITGAM","ITIH1","JAK2","KCNIP2","KIF2A","LAMP2","LCK","LGMN","LIPA","LRP1","LTA4H","LYN","MAFF","ME1","MMP12","MMP13","MMP14","MMP8","NOTCH4","OLR1","PCSK9","PDGFB","PHEX","PIK3CA","PIK3CG","PIM1","PLA2G4A","PLA2G7","PLAT","PLAUR","PLG","PLSCR1","PRCP","PREP","PRKCD","PSEN1","RAF1","RASGRP1","SERPINB2","SERPINC1","SERPINE1","SIRT6","SRC","TFPI2","TIMP2","TMPRSS6","WAS","XPNPEP1","ZEB1")
)


mTEC_score <- JoinLayers(mTEC) %>%
  AddModuleScore(
    features = Complement_geneset_ls,
    name = names(Complement_geneset_ls)
  )

mTEC_score %>% head

# 按照subtype
# 5. 箱线图：Complement ------------------------------------------
p1 <- mTEC_score@meta.data %>%
  ggplot(aes(x = subtype, y = Complement1, fill = subtype)) +
  geom_boxplot(outlier.shape = NA) +          # 箱线图
  geom_jitter(
    width = 0.2,
    height = 0,
    alpha = 0.5,
    size = 1,
    color = "steelblue"
  ) +
  ggpubr::geom_pwc(tip.length = 0, hide.ns = FALSE, label = 'p.signif') +
  labs(x = "Subtype", y = "IFN-alpha Score") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_blank(),
    legend.position = "none"
  )

ggsave(
  plot = p1,
  filename = paste0(mTEC_out_dir_3, "/Complement_score_subtype_Boxplot.pdf"),
  width = 5, height = 8
)






# 映射到空间切片
subtype_col <- 'subtype'

cells_position_ls <-
  cbind(
    mTEC[['sample']],
    Embeddings(mTEC[['spatial']])
  ) %>%
  rownames_to_column('cell_id') %>%
  group_by(sample) %>%
  nest() %>%
  deframe()

some_genes <- c("CD19", "CD38", "CD138")

plot_dat <-
  FetchData(
    mTEC, 
    vars = c(subtype_col, some_genes)
  ) %>%
  rownames_to_column('cell_id')


# plot celltype ------------------
all_cells <-  mTEC$subtype %>% unique %>% as.character()
cell_colors <- structure(
  names = all_cells,
  ggsci:: pal_igv()(length(all_cells))
)



p_sp_ct_ls <-
  cells_position_ls %>% # .[1] %>%
  imap(function(coords, sam){
    dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
    
    height <- 5
    width  <- 5
    p <-
      dat2 %>%
      ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
      geom_tile(
        aes(fill = .data[[subtype_col]]),
        height = height, 
        width = width,
        show.legend = c(fill = TRUE, alpha = FALSE)
      ) +
      theme_void() +
      coord_equal() +
      scale_y_reverse() +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values = cell_colors) +
      theme(panel.background = element_rect(fill = "black", color = NA)) +
      ggtitle(sam)
    p
    
  })


pdf(
  file.path("mTEC_subtype_3/spatial plot of celltype in each sample.pdf"),
  height = 10,
  width  = 10
)
p_sp_ct_ls %>% walk(print)
dev.off()


# plot each cell in all sample -----------
p_sp_ct_ls_ls <-
  cells_position_ls %>% # .[1] %>%
  imap(function(coords, sam){
    dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
    
    height <- 5
    width  <- 5
    
    p_ls <-
      all_cells %>%
      setNames(., .) %>%
      map(function(x){
        cell_colors[names(cell_colors) != x] <- NA
        
        p <-
          dat2 %>%
          ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
          geom_tile(
            aes(fill = .data[[subtype_col]]),
            height = height, 
            width = width,
            show.legend = c(fill = TRUE, alpha = FALSE)
          ) +
          theme_void() +
          coord_equal() +
          scale_y_reverse() +
          # ggsci::scale_fill_igv() +
          scale_fill_manual(values = cell_colors) +
          theme(panel.background = element_rect(fill = "black", color = NA)) +
          ggtitle(sam)
        p
      })
    p_ls
  })

# export by sample
p_sp_ct_ls_ls %>%
  iwalk(function(p_ls, sample){
    pdf(
      file.path(str_glue("mTEC_subtype_3/spatial plot of celltype in each sample - {sample}.pdf")),
      height = 10,
      width  = 10
    )
    p_ls %>% walk(print)
    dev.off()
  })



## enrichment ---------------------------------------
mTEC_de_genes_for_enrichmenet <- 
  mTEC_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.1) %>%
  dplyr::filter(abs(avg_log2FC) > 0.25) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")

mTEC_er_output<-'mTEC_subtype_3/enrichment' 
dir.create(mTEC_er_output)
enrich_ls <-
  mTEC_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichGO(
      x,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
      ont='ALL',
      pAdjustMethod = 'BH',
      pvalueCutoff = 0.05, 
      qvalueCutoff = 0.05,
      keyType = 'SYMBOL'
    )
  })


# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(mTEC_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(mTEC_er_output, name, "Go enriment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    # if(nrow(as.data.frame(x)) < 1){
    #   return(NULL)
    # }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

# export dotplot
enrich_dot_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(mTEC_er_output, name, "go dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(mTEC_er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })


mTEC_de_genes_id_for_enrichmenet <-
  mTEC_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::bitr(
      x,
      OrgDb = "org.Hs.eg.db",
      fromType = "SYMBOL",
      toType = "ENTREZID"
    ) %>%
      .$ENTREZID
  }) %>%
  purrr::discard(function(x) length(x) < 1)


enrich_ls <-
  mTEC_de_genes_id_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichKEGG(
      x,
      organism = "hsa", 
      keyType = "ncbi-geneid", 
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
    )
  })

# saveRDS(mTEC_er_output, 'enrich_ls.rds')

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(mTEC_er_output, name))
    
    write.csv(
      x@result,
      file = file.path(mTEC_er_output, name, "KEGG enriment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

enrich_bar_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    barplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

# export dotplot
enrich_dot_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(mTEC_er_output, name, "kegg dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

enrich_bar_ls %>%
  purrr::discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    # cluster_name <- paste0("cluster_", name)
    # dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(mTEC_er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })

# saveRDS(mTEC, 'mTEC_subtype_3/all_mTEC.rds')