#Fibroblast亚类加入all_data2--------
stromal_obs <- read.csv("stromal_obs.csv")
stromal_obs$cell_id_sample %>% head()

stromal_subtype <- stromal_obs %>% 
  dplyr::select("cell_id_sample","subtype","subtype_leiden") %>%
  mutate(subtype_stromal = subtype,
         leiden_stromal = subtype_leiden
         ) %>% 
  dplyr::select(-"subtype",-"subtype_leiden")
  
ttt <- all_data2@meta.data$cell_id_sample %>% duplicated()
all_data2@meta.data$cell_id_sample %>% head()
all_data2@meta.data$subtype %>% table()

# 检查重复
if (any(duplicated(stromal_subtype$cell_id_sample))) {
  print("stromal_subtype中有重复的cell_id_sample")
  stromal_subtype <- stromal_subtype[!duplicated(stromal_subtype$cell_id_sample), ]
}

all_data2@meta.data <- all_data2@meta.data %>%
  left_join(stromal_subtype, by = "cell_id_sample")

all_data2@meta.data$subtype_stromal %>% table()

# 对于非基质细胞，保持原来的subtype
all_data2@meta.data <- all_data2@meta.data %>%
  mutate(
    subtype_stromal = ifelse(is.na(subtype_stromal), as.character(subtype), subtype_stromal)
  )

print(table(all_data2@meta.data$subtype_stromal, useNA = "always"))

rownames(all_data2@meta.data) <- paste0(all_data2@meta.data$project,"_",all_data2@meta.data$cell_id)
ttt <- data_frame(rownames(all_data2@meta.data),all_data2@meta.data$project,all_data2@meta.data$cell_id)

##加入stromal亚类的all_data------
saveRDS(all_data2,"result0127/all_data.rds")

#5、EAOC、OE两组样本基质细胞比例统计（整体、亚类）-------
metadata <- all_data2[[]]

# metadata$sample <- factor(metadata$sample, 
#                             levels=c(metadata$sample %>% unique())
# )
p_col_proj <-
  metadata %>%
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
  filename = file.path("result0127/cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)


p_col_group <-
  metadata %>%
  group_by(group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_group,
  filename = file.path("result0127/cell_freq_group.pdf"),
  height = 6,
  width  = 8
)

##分面细胞比例统计------
p_col_group <-
  metadata %>%
  filter(group %in% c("EAOC","OE")) %>% 
  # filter(subtype %in% c("Stromal")) %>% 
  group_by(group,subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = group)) + 
  facet_grid(~subtype)+
  geom_col() +
  #ggsci::scale_fill_igv() +
  scale_fill_manual(values = c("#4DBBD5", "#E64B35"))+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_group,
  filename = file.path("result0127/cell_freq_group_OE_EAOC_Stromal.pdf"),
  height = 6,
  width  = 12
)

library(Seurat)
library(tidyverse)
library(BPCells)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tibble)
library(data.table)

table(all_data2$group, all_data2$subtype)
table(all_data2$celltype)

#1、EAOC_Tumor再分群--------
EAOC_Tumor <- subset(all_data2, subtype =='Tumor Epi')
#EAOC_Tumor@images <- list()

## PCA --------
DefaultAssay(EAOC_Tumor) <- "RNA"
EAOC_Tumor <- FindVariableFeatures(EAOC_Tumor)
EAOC_Tumor <- ScaleData(EAOC_Tumor, vars.to.regress = 'nCount_RNA')
EAOC_Tumor <- RunPCA(EAOC_Tumor)

EAOC_Tumor_out_dir<-'result0127/EAOC_Tumor_subtype/' 
dir.create(EAOC_Tumor_out_dir)


## pca inspector ------
pdf(file.path(EAOC_Tumor_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(EAOC_Tumor, ndims = 50)
dev.off()

pdf(file.path(EAOC_Tumor_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(EAOC_Tumor, dims = 1:5)
dev.off()
pdf(file.path(EAOC_Tumor_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(EAOC_Tumor, dims = 6:10)
dev.off()
pdf(file.path(EAOC_Tumor_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(EAOC_Tumor, dims = 11:15)
dev.off()
pdf(file.path(EAOC_Tumor_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(EAOC_Tumor, dims = 16:20)
dev.off()


## dimension reduction --------
dims <- 20

library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


EAOC_Tumor <- RunUMAP(EAOC_Tumor, dims = 1:dims,return.model = TRUE)

pdf(file.path(EAOC_Tumor_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(EAOC_Tumor, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

EAOC_Tumor <- EAOC_Tumor %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

EAOC_Tumor <- RunUMAP(
  EAOC_Tumor, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(EAOC_Tumor_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(EAOC_Tumor, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'

### cluster -------
EAOC_Tumor  <-
  FindNeighbors(
    EAOC_Tumor ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


EAOC_Tumor <-
  FindClusters(
    EAOC_Tumor,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(EAOC_Tumor_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(EAOC_Tumor, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(EAOC_Tumor, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(EAOC_Tumor, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(EAOC_Tumor, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


### dotplot -------
EAOC_Tumor_major_marker_ls <- c(
  # AT2：肺泡 2 型细胞
  "AT2" = "Sftpc",
  "AT2" = "Sftpd",
  "AT2" = "Pgc",
  "AT2" = "Napsa",
  "AT2" ="Wif1", 
  "AT2" ="Etv5", 
  "AT2" ="Abca3",
  
  # Basal：基底细胞
  "Basal" = "Krt5",
  "Basal" = "Krt15",
  "Basal" = "Krt17",
  "Basal" = "Trp63",
  
  # AT1：肺泡 1 型细胞
  "AT1" = "Rtkn2",
  "AT1" = "Ager",
  "AT1" = "Col4a3",
  
  # Transitional AT2：过渡态肺泡 2 型细胞
  "Transitional AT2" = "Lamp3",
  
  # Ciliated：纤毛细胞
  "Ciliated" = "C20orf85",  # 保留原名（无标准小鼠同源基因）
  "Ciliated" = "Foxj1",
  "Ciliated" = "Scgb1a1",
  "Ciliated" = "Trp73",
  
  # Secretory：分泌细胞
  "Secretory" = "Muc5b",
  "Secretory" = "Wfdc2",
  "Secretory" = "Mmp7",
  "Secretory" = "Scgb3a2",
  
  # PNEC：肺神经内分泌细胞
  "PNEC" = "Scg2",
  "PNEC" = "Chgb",
  "PNEC" = "Calca",
  "PNEC" = "Ascl1",
  
  # Mast：肥大细胞
  "Mast" = "Cpa3",
  "Mast" = "Tpsab1",
  "Mast" = "Kit",
  
  # Proliferating：增殖细胞
  "Proliferating" = "Top2a",
  "Proliferating" = "Mki67",
  "Proliferating" = "Cdk1",
  
  # 未明确分类细胞
  "Unclassified" = "Nkx2-1",
  "Unclassified" = "Bmp4",
  "Unclassified" = "Spink1"
)
valid_EAOC_Tumorgenes <- EAOC_Tumor_major_marker_ls[EAOC_Tumor_major_marker_ls %in% rownames(EAOC_Tumor)]
cat("有效基因:", length(valid_EAOC_Tumorgenes), "/", length(EAOC_Tumor_major_marker_ls))

# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_gene <- setdiff(EAOC_Tumor_major_marker_ls, valid_EAOC_Tumorgenes)
missing_gene


# find markers
Idents(EAOC_Tumor) <- "RNA_snn_res.0.4"
EAOC_Tumor_markers <-
  FindAllMarkers(
    JoinLayers(EAOC_Tumor),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

EAOC_Tumor_top10_markers <- 
  EAOC_Tumor_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10) %>%
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

names(EAOC_Tumor_top10_markers) <- sub('.*\\.','',names(EAOC_Tumor_top10_markers))


#提取deepseek
EAOC_Tumor_top10_markers1 <-
  EAOC_Tumor_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10)


# dotplot
pdf(file.path(EAOC_Tumor_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 35, height = 6)
DotPlot(EAOC_Tumor, group.by = "RNA_snn_res.0.4", features = major_marker)+RotatedAxis()
DotPlot(EAOC_Tumor, group.by = "RNA_snn_res.0.4", features = EAOC_Tumor_major_marker_ls)+RotatedAxis()
DotPlot(EAOC_Tumor, group.by = "RNA_snn_res.0.4", features = EAOC_Tumor_top10_markers)+RotatedAxis()
dev.off()


for (gene in valid_EAOC_Tumorgenes) {
  p <- FeaturePlot(EAOC_Tumor, features = gene, reduction = 'harmonyumap')+RotatedAxis()
  ggsave(
    filename = file.path(EAOC_Tumor_out_dir, paste0("featureplot/",gene, "_feature.png")),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}


## annotation -----------
cell_cluster <- c(
  '0' = "AT2", # SFTPC SFTPB SFTPD
  '1' = "Ciliated", # FOXJ1 TPPP3 DNAH5
  '2' = "AT1", # AGER PDPN CLDN4
  '3' = "AT2", # SFTPC SFTPB SFTPD
  '4' = "Ciliated", # FOXJ1 TPPP3 DNAH5
  '5' = "AT1", # AGER PDPN CLDN4
  '6' = "AT2", # SFTPC SFTPB SFTPD
  '7' = "AT1", # AGER PDPN CLDN4
  '8' = "AT1", # AGER PDPN CLDN4
  '9' = "AT2", # SFTPC SFTPB SFTPD
  '10' = "Inf_EAOC_Tumor", # IFITM3 ISG15 IFI44L
  '11' = "delete" # no marker
)

Idents(EAOC_Tumor) <- "RNA_snn_res.0.4"

EAOC_Tumor <- RenameIdents(EAOC_Tumor, cell_cluster)
EAOC_Tumor$subtype <- Idents(EAOC_Tumor)

pdf(file.path(EAOC_Tumor_out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(EAOC_Tumor, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(EAOC_Tumor_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 7, width = 18)
DimPlot(EAOC_Tumor, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

## cell props barplot  --------
EAOC_Tumor_meta <- EAOC_Tumor[[]]

EAOC_Tumor_meta$sample <- factor(EAOC_Tumor_meta$sample, 
                          levels=c(EAOC_Tumor_meta$sample %>% unique())
)
p_col_proj <-
  EAOC_Tumor_meta %>%
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
  filename = file.path("result0127/EAOC_Tumor_subtype/cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  EAOC_Tumor_meta %>%
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
  filename = file.path("result0127/EAOC_Tumor_subtype/cell_freq_group.pdf"),
  height = 6,
  width  = 4
)


## volcano-------
EAOC_Tumor_markers <-
  FindAllMarkers(
    JoinLayers(EAOC_Tumor),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(EAOC_Tumor_markers,'result0127/EAOC_Tumor_subtype/EAOC_Tumor_diff_gene.csv')

library(EnhancedVolcano)
type<-unique(EAOC_Tumor_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- EAOC_Tumor_markers[EAOC_Tumor_markers$cluster == type, ]
    p<-EnhancedVolcano(markers,
                       lab = rownames(markers),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = str_glue('{name} volcano plot'),
                       pCutoff = 10e-2,
                       FCcutoff = 0,
                       pointSize = 3.0,
                       labSize = 6.0)
    p
  })
volcano_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path('result0127/EAOC_Tumor_subtype', str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })

EAOC_Tumor_markers_ls <- list(EAOC_Tumor_markers)
names(EAOC_Tumor_markers_ls) <- 'EAOC_Tumor'

# scatter
p_volcanoplot_major_ls <- 
  EAOC_Tumor_markers_ls %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 0.5
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
        y = ifelse(avg_log2FC > 0, avg_log2FC + rect_height, avg_log2FC - rect_height)
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
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_point(aes(color = group), size = 1) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
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

p_volcanoplot_major_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path('result0127/EAOC_Tumor_subtype', str_glue("{name} diffgene dot.pdf")),
      height = 10,
      width = 20
    )
  })


## enrichment ----------------------------
EAOC_Tumor_de_genes_for_enrichmenet <- 
  EAOC_Tumor_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.15) %>%
  dplyr::filter(abs(avg_log2FC) > 0.5) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")


enrich_ls <- 
  EAOC_Tumor_de_genes_for_enrichmenet %>%
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


er_output<-'result0127/EAOC_Tumor_subtype/enrichment' 
dir.create(er_output)

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(er_output, name))
    
    write.csv(
      x@result,
      file = file.path(er_output, name, "Go enrichment table.csv"),
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
      filename = file.path(er_output, name, "go dotplot.pdf"),
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
      filename = file.path(er_output, name, "go barplot.pdf"),
      height = 6,
      width = 8
    )
  })

EAOC_Tumor_de_genes_id_for_enrichmenet <-
  EAOC_Tumor_de_genes_for_enrichmenet %>%
  map(function(x){ # browser()
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
  EAOC_Tumor_de_genes_id_for_enrichmenet %>%
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

# saveRDS(enrich_ls, file.path(er_output, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>%
  # discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    dir.create(file.path(er_output, name))
    
    write.csv(
      x@result,
      file = file.path(er_output, name, "KEGG enriment table.csv"),
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
      filename = file.path(er_output, name, "kegg dotplot.pdf"),
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
      filename = file.path(er_output, name, "kegg barplot.pdf"),
      height = 6,
      width = 8
    )
  })

saveRDS(EAOC_Tumor, 'result0127/EAOC_Tumor_subtype/EAOC_Tumor.rds')


#2. 上皮细胞：计算 EAOC_All_Normal v.s. OE_ALL 差异表达基因-------
##转到203分析
table(all_data2$group, all_data2$subtype)

EPI_oe_eaoc <- subset(all_data2, subtype %in% c("Normal Epi")&group %in% c("EAOC","OE"))

table(EPI_oe_eaoc$group, EPI_oe_eaoc$subtype)

saveRDS(EPI_oe_eaoc,"result0127/EPI_oe_eaoc/EPI_oe_eaoc.rds")

EPI_oe_eaoc <- readRDS("result0127/EPI_oe_eaoc/EPI_oe_eaoc.rds")
library(tidyverse)
library(Seurat)
library(patchwork)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(RColorBrewer)
library(EnhancedVolcano)
library(dplyr)




output_dir <- "result0127/EPI_oe_eaoc/FindMarkers"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

celltype <- "Normal_Epi"

subset_obj <- JoinLayers(EPI_oe_eaoc)

Idents(subset_obj) <- 'group'

diffmarkers <- FindMarkers(
  subset_obj,
  ident.1 = 'EAOC',
  ident.2 = 'OE',
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

nrow(diffmarkers)

diffmarkers_up <- diffmarkers %>% 
  filter(avg_log2FC > 0) %>% 
  arrange(desc(avg_log2FC))
diffmarkers_down <- diffmarkers %>% 
  filter(avg_log2FC < 0) %>% 
  arrange(desc(avg_log2FC))

nrow(diffmarkers_up)
nrow(diffmarkers_down)

write.csv(diffmarkers, file.path(output_dir, "Normal_Epi_group_diffmarker.csv"))
write.csv(diffmarkers_up, file.path(output_dir, "Normal_Epi_group_diffmarker_up.csv"))
write.csv(diffmarkers_down, file.path(output_dir, "Normal_Epi_group_diffmarker_down.csv"))

# 富集分析 - 上调基因
# GO 上调
enrichGO_up_ls <- enrichGO(
  rownames(diffmarkers_up),
  OrgDb = org.Hs.eg.db, 
  ont = 'ALL',
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05,
  keyType = 'SYMBOL'
)

write.csv(enrichGO_up_ls@result, 
          file.path(output_dir, "Normal_Epi_GO_up_enrichment_table.csv"), 
          row.names = FALSE)

# GO 上调点图
enrichGO_up_dot <- enrichplot::dotplot(
    enrichGO_up_ls,
    showCategory = 10,
    label_format = 50
  )
ggsave(
    enrichGO_up_dot,
    filename = file.path(output_dir, "Normal_Epi_GO_up_dotplot.pdf"),
    height = 6,
    width = 8
  )


# GO 上调条形图

enrichGO_up_bar <- barplot(
    enrichGO_up_ls,
    showCategory = 10,
    label_format = 50
  )
ggsave(
    enrichGO_up_bar,
    filename = file.path(output_dir, "Normal_Epi_GO_up_barplot.pdf"),
    height = 6,
    width = 8
  )


# KEGG 上调
de_up_genes_id_for_enrichment <- 
  clusterProfiler::bitr(
    rownames(diffmarkers_up),
    OrgDb = "org.Hs.eg.db",
    fromType = "SYMBOL",
    toType = "ENTREZID"
  ) %>%
  .$ENTREZID


enrichKEGG_up_ls <- enrichKEGG(
    de_up_genes_id_for_enrichment,
    organism = "hsa", 
    keyType = "ncbi-geneid",
    pvalueCutoff = 0.05, 
    pAdjustMethod = "BH", 
    minGSSize = 10,
    maxGSSize = 500, 
    qvalueCutoff = 0.2,
    use_internal_data = F
)
  

enrichKEGG_up_ls <- setReadable(enrichKEGG_up_ls, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
write.csv(enrichKEGG_up_ls@result, 
          file.path(output_dir, "Normal_Epi_KEGG_up_enrichment_table.csv"), 
          row.names = FALSE)

# KEGG 上调点图
enrichKEGG_up_dot <- enrichplot::dotplot(
  enrichKEGG_up_ls,
  showCategory = 10,
  label_format = 50
)
ggsave(
  enrichKEGG_up_dot,
  filename = file.path(output_dir, "Normal_Epi_KEGG_up_dotplot.pdf"),
  height = 6,
  width = 8
)

# KEGG 上调条形图

enrichKEGG_up_bar <- barplot(
  enrichKEGG_up_ls,
  showCategory = 10,
  label_format = 50
)
ggsave(
  enrichKEGG_up_bar,
  filename = file.path(output_dir, "Normal_Epi_KEGG_up_barplot.pdf"),
  height = 6,
  width = 8
)


# 富集分析 - 下调基因

cat("进行下调基因GO富集分析...\n")

# GO 下调
enrichGO_down_ls <- enrichGO(
  rownames(diffmarkers_down),
  OrgDb = org.Hs.eg.db, 
  ont = 'ALL',
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05,
  keyType = 'SYMBOL'
)

write.csv(enrichGO_down_ls@result, 
          file.path(output_dir, "Normal_Epi_GO_down_enrichment_table.csv"), 
          row.names = FALSE)

# GO 下调点图

enrichGO_down_dot <- enrichplot::dotplot(
    enrichGO_down_ls,
    showCategory = 10,
    label_format = 50
  )
ggsave(
    enrichGO_down_dot,
    filename = file.path(output_dir, "Normal_Epi_GO_down_dotplot.pdf"),
    height = 6,
    width = 8
  )


# GO 下调条形图

enrichGO_down_bar <- barplot(
    enrichGO_down_ls,
    showCategory = 10,
    label_format = 50
  )
ggsave(
    enrichGO_down_bar,
    filename = file.path(output_dir, "Normal_Epi_GO_down_barplot.pdf"),
    height = 6,
    width = 8
  )


# KEGG 下调
de_down_genes_id_for_enrichment <- 
  clusterProfiler::bitr(
    rownames(diffmarkers_down),
    OrgDb = "org.Hs.eg.db",
    fromType = "SYMBOL",
    toType = "ENTREZID"
  ) %>%
  .$ENTREZID


enrichKEGG_down_ls <- enrichKEGG(
  de_down_genes_id_for_enrichment,
  organism = "hsa", 
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05, 
  pAdjustMethod = "BH", 
  minGSSize = 10,
  maxGSSize = 500, 
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
)

  
enrichKEGG_down_ls <- setReadable(enrichKEGG_down_ls, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
write.csv(enrichKEGG_down_ls@result, 
          file.path(output_dir, "Normal_Epi_KEGG_down_enrichment_table.csv"), 
          row.names = FALSE)

# KEGG 下调点图
enrichKEGG_down_dot <- enrichplot::dotplot(
  enrichKEGG_down_ls,
  showCategory = 10,
  label_format = 50
)
ggsave(
  enrichKEGG_down_dot,
  filename = file.path(output_dir, "Normal_Epi_KEGG_down_dotplot.pdf"),
  height = 6,
  width = 8
)

# KEGG 下调条形图

enrichKEGG_down_bar <- barplot(
  enrichKEGG_down_ls,
  showCategory = 10,
  label_format = 50
)
ggsave(
  enrichKEGG_down_bar,
  filename = file.path(output_dir, "Normal_Epi_KEGG_down_barplot.pdf"),
  height = 6,
  width = 8
)


# 火山图
p_volcano <- EnhancedVolcano(
  diffmarkers,
  lab = rownames(diffmarkers),
  x = 'avg_log2FC',
  y = 'p_val_adj',
  title = '',
  pCutoff = 0.05,
  FCcutoff = 0.25,
  pointSize = 2.0,
  labSize = 3.0
)
ggsave(
  p_volcano,
  filename = file.path(output_dir, "Normal_Epi_volcano.pdf"),
  height = 8,
  width = 10
)



#3、上述差异基因在Eu中是否有体现?----------
result3_dir <- "result0127/prevDiff_EU_show"
dir.create(result3_dir)
table(all_data2$group, all_data2$subtype)

EPI_Eu <- subset(all_data2, subtype %in% c("Normal Epi")&group %in% c("EAOC_Eu","OE_Eu"))

pdf(file.path(result3_dir, str_glue("Dimplot subtype group - harmony.pdf")), height = 6, width = 8)
DimPlot(EPI_Eu, group.by = "subtype", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(EPI_Eu, group.by = "subtype", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()




## 筛选上调下调基因------
# EAOC_Tumor_top10_markers <- 
#   EAOC_Tumor_markers %>% 
#   group_by(cluster) %>%
#   dplyr::filter(p_val_adj < 0.05) %>%
#   arrange(desc(avg_log2FC)) %>%
#   slice_head(n = 30) %>%
#   arrange(desc(pct.1)) %>%
#   slice_head(n = 10) %>%
#   nest() %>%
#   deframe() %>% 
#   map(pull, 'gene') %>%
#   imap(function(genes, name){
#     re <- genes %>%
#       map(function(gene){
#         # browser()
#         names(gene) = str_glue("cluster_{name}")
#         gene
#       })
#     re
#   }) %>%
#   unlist() %>% 
#   .[!duplicated(.)]
# 
# names(EAOC_Tumor_top10_markers) <- sub('.*\\.','',names(EAOC_Tumor_top10_markers))



# diffmarkers_up <- diffmarkers %>% 
#   filter(avg_log2FC > 0) %>% 
#   arrange(desc(avg_log2FC))
# diffmarkers_down <- diffmarkers %>% 
#   filter(avg_log2FC < 0) %>% 
#   arrange(desc(avg_log2FC))
# 
# cat("上调基因数:", nrow(diffmarkers_up), "\n")
# cat("下调基因数:", nrow(diffmarkers_down), "\n")


diff_up20 <- diffmarkers_up %>% 
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 20) #%>%
  #rownames()
write.csv(diff_up20, file.path(result3_dir, "Normal_Epi_diff_up20.csv"))


for (gene in diff_up20) {
  p <- FeaturePlot(EPI_Eu, features = gene, reduction = 'harmonyumap')+RotatedAxis()
  ggsave(
    filename = file.path(result3_dir, paste0("featureplot/diff_up20_", gene, "_feature.png")),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}


diff_down20 <- diffmarkers_down %>% 
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 20) #%>%
  #rownames()

write.csv(diff_down20, file.path(result3_dir, "Normal_Epi_diff_down20.csv"))

for (gene in diff_down20) {
  p <- FeaturePlot(EPI_Eu, features = gene, reduction = 'harmonyumap')+RotatedAxis()
  ggsave(
    filename = file.path(result3_dir, paste0("featureplot/diff_down20_", gene, "_feature.png")),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}



#4、上皮细胞：三组梯度分析--------
#梯度分析：OE_上皮 → EAOC_ Normal → EAOC_ Tumor；三组细胞各自findallmarker，挑出线性升高的基因
table(all_data2$group, all_data2$subtype)

TumorEpi_OE_EAOC <- subset(all_data2, subtype %in% c("Normal Epi","Tumor Epi")&group %in% c("EAOC","OE"))

table(TumorEpi_OE_EAOC$group, TumorEpi_OE_EAOC$subtype)
table(TumorEpi_OE_EAOC$subtype)

TumorEpi_OE_EAOC@meta.data %>% head()

TumorEpi_OE_EAOC@meta.data$group_subtype <- paste0(TumorEpi_OE_EAOC@meta.data$group,"_",TumorEpi_OE_EAOC@meta.data$subtype)

table(TumorEpi_OE_EAOC$group, TumorEpi_OE_EAOC$group_subtype)
table(TumorEpi_OE_EAOC$group_subtype)

Idents(TumorEpi_OE_EAOC) <- "group_subtype"


result4_dir <- "result0127/TumorEpi_OE_EAOC/"
dir.create(result4_dir)

##FindAllMarker-------
### volcano-------
TumorEpi_OE_EAOC_markers <-
  FindAllMarkers(
    JoinLayers(TumorEpi_OE_EAOC),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )
write.csv(TumorEpi_OE_EAOC_markers,file.path(result4_dir, 'TumorEpi_OE_EAOC_diff_gene.csv'))

library(EnhancedVolcano)
type<-unique(TumorEpi_OE_EAOC_markers$cluster)
volcano_ls <- type %>% 
  setNames(.,.) %>%
  imap(function(type,name){
    markers<- TumorEpi_OE_EAOC_markers[TumorEpi_OE_EAOC_markers$cluster == type, ]
    p<-EnhancedVolcano(markers,
                       lab = rownames(markers),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = str_glue('{name} volcano plot'),
                       pCutoff = 10e-2,
                       FCcutoff = 0,
                       pointSize = 3.0,
                       labSize = 6.0)
    p
  })
volcano_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(result4_dir, str_glue("{name} diffgene valcano.pdf")),
      height = 6,
      width = 8
    )
  })

TumorEpi_OE_EAOC_markers_ls <- list(TumorEpi_OE_EAOC_markers)
names(TumorEpi_OE_EAOC_markers_ls) <- 'TumorEpi_OE_EAOC'

# scatter
p_volcanoplot_major_ls <- 
  TumorEpi_OE_EAOC_markers_ls %>%
  imap(function(dat, name){
    
    set.seed(1234)
    rect_height <- 0.5
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
        y = ifelse(avg_log2FC > 0, avg_log2FC + rect_height, avg_log2FC - rect_height)
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
      geom_rect(
        data = df_grp_rect,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
        inherit.aes = FALSE
      ) +
      geom_point(aes(color = group), size = 1) +
      
      ggrepel::geom_text_repel(data = df_label_up, aes(label = gene)) +
      ggrepel::geom_text_repel(data = df_label_down, aes(label = gene)) +
      scale_color_manual(values = c(
        'p < 0.01' = "red",
        'p > 0.01' = "gray70"
      )) +
      # scale_fill_manual(values = grp_cols_sc) +
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

p_volcanoplot_major_ls %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(result4_dir, str_glue("{name} diffgene dot.pdf")),
      height = 10,
      width = 20
    )
  })

###top30特征基因dotplot----
TumorEpi_OE_EAOC_top10_markers <-
  TumorEpi_OE_EAOC_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  # arrange(desc(pct.1)) %>%
  # slice_head(n = 10) #%>%
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

names(TumorEpi_OE_EAOC_top10_markers) <- sub('.*\\.','',names(TumorEpi_OE_EAOC_top10_markers))
names(TumorEpi_OE_EAOC_top10_markers) <- sub('cluster_','',names(TumorEpi_OE_EAOC_top10_markers))

# dotplot
pdf(file.path(result4_dir, str_glue("TumorEpi_OE_EAOC_top10_markers dotplot .pdf")), width = 35, height = 6)
DotPlot(TumorEpi_OE_EAOC, 
        #group.by = "RNA_snn_res.0.4", 
        features = TumorEpi_OE_EAOC_top10_markers)+RotatedAxis()
dev.off()


###top100的三种细胞共有基因dotplot-----
TumorEpi_OE_EAOC_top10_markers111 <-
  TumorEpi_OE_EAOC_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 100) %>%
  # arrange(desc(pct.1)) %>%
  # slice_head(n = 10) #%>%
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
  .[duplicated(.)]

names(TumorEpi_OE_EAOC_top10_markers111) <- sub('.*\\.','',names(TumorEpi_OE_EAOC_top10_markers111))
names(TumorEpi_OE_EAOC_top10_markers111) <- sub('cluster_','',names(TumorEpi_OE_EAOC_top10_markers111))

pdf(file.path(result4_dir, str_glue("TumorEpi_OE_EAOC_common_markers dotplot .pdf")), width = 35, height = 6)
DotPlot(TumorEpi_OE_EAOC, 
        #group.by = "RNA_snn_res.0.4", 
        features = TumorEpi_OE_EAOC_top10_markers111)+RotatedAxis()
dev.off()



# 
# allgene <- TumorEpi_OE_EAOC_markers %>% 
#   arrange(desc(avg_log2FC)) %>%
#   slice_head(n = 200) %>%
#   rownames()
# 
# 
# pdf(file.path(result4_dir, str_glue("TumorEpi_OE_EAOC_all_diff_gene200 dotplot .pdf")), width = 35, height = 6)
# DotPlot(TumorEpi_OE_EAOC, 
#         #group.by = "RNA_snn_res.0.4", 
#         features = allgene)+RotatedAxis()
# dev.off()

###top20gene热图--------
anno <- data.frame(celltype=TumorEpi_OE_EAOC$group_subtype,
                   group=TumorEpi_OE_EAOC$group,
                   #cluster=sce$seurat_clusters, 
                   row.names = colnames(TumorEpi_OE_EAOC))
head(anno)
anno <- anno[order(anno$celltype), ]

data <- GetAssayData(JoinLayers(TumorEpi_OE_EAOC), slot = "scale.data")
data <- data[, rownames(anno)]


available_genes <- intersect(TumorEpi_OE_EAOC_top10_markers, rownames(data))
missing_genes <- setdiff(TumorEpi_OE_EAOC_top10_markers, rownames(data))

data <- data[available_genes, ]


# Specify colors
ann_colors = list(
  celltype=c("OE_Normal Epi" = "#117733",
             "EAOC_Normal Epi" = "#cc61b0",
             "EAOC_Tumor Epi" = "#86cde4"
    # 'Naive CD4 T'="#117733",'Memory CD4 T'="#cc61b0",'CD14+ Mono'="#99c945",
    #          'B'="#2f8ac4",'CD8 T'="#8d59ea",'FCGR3A+ Mono'="#dddcd8",'NK'="#d12125",
    #          'DC'="#1daa7e",'Platelet'="#86cde4")
))

ann_colors


library(pheatmap)

p <- pheatmap(data,
              scale = "row",
              treeheight_row = 0,
              treeheight_col = 0,
              show_colnames = FALSE,
              annotation_col = anno[, 1, drop = FALSE],
              annotation_names_col = FALSE,
              cluster_cols = FALSE,
              cluster_rows = FALSE,
              annotation_colors = ann_colors,
              col = colorRampPalette(c("#ceb0ff", "white", "#ff4628"))(100),
              fontsize_row = 8
)

pdf(file = file.path(result4_dir, "group_subtype_top30.pdf"), width = 9,height = 10)
print(p)
dev.off()










