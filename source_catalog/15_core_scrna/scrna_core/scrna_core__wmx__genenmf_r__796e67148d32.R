# BiocManager::install("UCell")
# install.packages('GeneNMF',type = 'binary')
# remotes::install_github("carmonalab/GeneNMF")  # or from GitHub
.libPaths("D:/Softwares/Anaconda/envs/r_env/Lib/R/library") # 选取lib路径
install.packages("S7")
library(GeneNMF)
# data(sampleObj) # test 一下
# sampleObj <- runNMF(sampleObj, k=5)
# options(device = NULL)
library(patchwork)
library(Seurat)
library(ggplot2)
library(UCell)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(tidyverse)

setwd("D:/Project/geneNMF/")

# 读入数据 
seu <- readRDS("Tumor_result.rds")

# DimPlot(seu, group.by="group") + theme(aspect.ratio = 1)
# DimPlot(seu,reduction='umap',group.by='cell_type',raster=FALSE)
# DefaultAssay(seu) <- "RNA"

# seu.list <- SplitObject(seu, split.by = "sample")
# mt_genes <- rownames(seu)[grepl("^MT-", rownames(seu))]
# ribo <- rownames(seu)[grepl("^(RPS[0-9]+|RPL[0-9]+)", rownames(seu))]
# delete_genes = c(mt_genes,ribo)
# 加载已经完成的NMF结果
load("geneNMF.Rdata")

### 核心运行代码
# geneNMF.programs <- multiNMF(seu.list, assay="RNA", k=4:12,hvg.blocklist = mt_genes)
# 循环保存结果，查看图片，确定最后的选择
# geneNMF.programs <- multiNMF(seu.list, 
#                              assay = "RNA", 
#                              k = 4:6,
#                              min.cells.per.sample = 10,
#                              hvg=FALSE,
#                              hvg.blocklist = delete_genes, # 删除在高变基因内的线粒体基因的影响
#                              min.exp = 0.1, # 删除低表达基因的影响，Minimum average log-expression value for retaining genes
#                              seed = 2026
#                             )

# mplist = c(6,8,10)
# splist = c(2,5)
# weiglist = c(0.6,0.8)
# maxgenelist = c(150,200)
# melist = c("cosine", "jaccard")

# for(mp in mplist){
#   for(sp in splist){
#     for(weig in weiglist){
#       for(maxg in maxgenelist){
#         for(me in melist){
          
#           # 尝试执行，捕获错误
#           tryCatch({
#             geneNMF.metaprograms <- getMetaPrograms(
#               geneNMF.programs,
#               metric = me,
#               hclust.method = "ward.D2",
#               remove.empty = TRUE,
#               max.genes = maxg,
#               specificity.weight = sp,
#               weight.explained = weig,
#               nMP = mp
#             )
            
#             pdf_file <- paste0("test/","MP_",mp,"_Specificity_",sp,"_Weight_",weig,"_maxgene_",maxg,"_metric_",me,"_heatmap.pdf")
#             pdf(pdf_file, width = 15, height = 15)
#             plotMetaPrograms(geneNMF.metaprograms,
#                              show_rownames = FALSE,
#                              similarity.cutoff = c(0.1,1))
#             dev.off()
            
#             cat("Success:", pdf_file, "\n")
            
#           }, error = function(e) {
#             # 打印出错参数和错误信息
#             cat("ERROR for parameters: ",
#                 "mp=", mp, 
#                 ", sp=", sp, 
#                 ", weig=", weig, 
#                 ", maxg=", maxg, 
#                 ", me=", me,
#                 " - ", e$message, "\n")
#             # 确保图形设备关闭（如果已经打开）
#             if(dev.cur() != 1) dev.off()
#           })
          
#         }
#       }
#     }
#   }
# }
### 确定最终参数,MP= 6-8 合适 ,另外cosine的算法下聚类更清晰，specificity.weight = 5 效果会比2更好一些 。max.genes，weight.explained的影响不大
# 最终挑选的是MP = 6
# geneNMF.metaprograms <- getMetaPrograms(
#                                         geneNMF.programs,
#                                         metric = "cosine",
#                                         hclust.method = "ward.D2",
#                                         remove.empty = TRUE, # Whether to remove meta-programs with no genes above confidence threshold
#                                         max.genes = 200,
#                                         specificity.weight = 5,# A parameter controlling how specific gene should be for each program. 'specificity.weight=0' no constraint on specificity, and positive values impose increasing specificity.
#                                         weight.explained = 0.8, # Fraction of NMF weights explained by selected genes. For example if weight.explained=0.5, all genes that together account for 50% of NMF weights are used to return program signatures.
#                                         nMP= 6
#                                       )

# geneNMF.metaprograms$metaprograms.metrics
# plotMetaPrograms(geneNMF.metaprograms,show_rownames = F,
#                  similarity.cutoff = c(0.1,1))

### 聚类热图绘制
# dir.create('./geneNMF/')
# png(filename = "./geneNMF/heatmap_MP.png",width = 10,height = 10,res = 300,units = 'in')
# plotMetaPrograms(geneNMF.metaprograms,show_rownames = F,
#                  similarity.cutoff = c(0.1,1))
# dev.off()

# pdf(file = "./geneNMF/heatmap_MP.pdf",width = 10,height = 10)
# plotMetaPrograms(geneNMF.metaprograms,show_rownames = F,
#                  similarity.cutoff = c(0.1,1))
# dev.off()
# save(geneNMF.programs,geneNMF.metaprograms,file = "geneNMF.Rdata")
### 对应MPs的基因热图绘制
# all_w <- lapply(geneNMF.programs, function(x) x$w)
# all_w_renamed <- lapply(names(all_w), function(df_name) {
#   # 提取当前数据框
#   df <- all_w[[df_name]]
#   # 将列名中的"pattern"替换为数据框的名称（df_name）
#   colnames(df) <- gsub("pattern", df_name, colnames(df))
#   # 返回修改后的数据框
#   return(df)
# })
# # 保持列表的名称与原all_w一致

# names(all_w_renamed) <- names(all_w)
# combined_df <- do.call(cbind, all_w_renamed)
# sorted_labels <- geneNMF.metaprograms$programs.tree$labels[geneNMF.metaprograms$programs.tree$order]

# combined_df_sorted <- combined_df[, sorted_labels]
# all_genes <- unlist(geneNMF.metaprograms$metaprograms.genes)

# top5_genes_per_mp <- lapply(geneNMF.metaprograms$metaprograms.genes, function(genes) {
#   head(genes, 5)  
# })


# all_top5_genes <- unlist(top5_genes_per_mp)

# unique_top5_genes <- unique(all_top5_genes)

# combined_df_fin <- combined_df_sorted[rownames(combined_df_sorted) %in% unique_top5_genes,]
# combined_df_fin <- as.data.frame(combined_df_fin)

# combined_df_fin <- t(scale(t(combined_df_fin)))
# pheatmap::pheatmap(combined_df_fin,show_colnames = F,cluster_cols = F,cluster_rows = F)


# geneNMF.programs %>% map('w') %>% do.call(cbind, .) %>% .[geneNMF.metaprograms$metaprograms.genes %>% map(head, 5) %>% unlist, ] %>%  pheatmap::pheatmap(scale = 'row', cluster_rows = F, cluster_cols = T)
# 
#  

### 相关表格输出
# library(tidyverse)
# # 1. 提取基因和权重列表
# genes_list <- geneNMF.metaprograms[["metaprograms.genes"]]
# weights_list <- geneNMF.metaprograms[["metaprograms.genes.weights"]]

# # 2. 检查列表长度是否一致（避免配对错误）
# if (length(genes_list) != length(weights_list)) {
#   stop("基因列表和权重列表长度不一致，请检查数据！")
# }

# # 3. 按元程序整合基因-权重对
# # 初始化空列表存储每个元程序的结果
# meta_program_df_list <- list()

# # 遍历每个元程序（i 为元程序编号）
# for (i in seq_along(genes_list)) {
#   # 提取当前元程序的基因和权重
#   genes <- genes_list[[i]]
#   weights <- weights_list[[i]]
  
#   # 构建当前元程序的临时数据框
#   temp_df <- tibble(
#     meta_program = paste0("MP", i),  # 元程序名称（可自定义）
#     gene = genes,                              # 基因名称
#     weight = weights                           # 对应权重
#   )
  
#   # 存入列表
#   meta_program_df_list[[i]] <- temp_df
# }

# # 4. 合并所有元程序的结果为一个总表格
# final_table <- bind_rows(meta_program_df_list)

# # 5. 输出表格（可选：保存为CSV/Excel，或打印查看）
# # 打印前10行预览
# head(final_table, 10)

# # 保存为CSV文件
# write.csv(final_table, 
#           file = "./geneNMF/geneNMF_metaprograms_genes_weights.csv", 
#           row.names = FALSE,  # 不保存行名
#           quote = FALSE)      # 避免基因名被引号包裹

# head(geneNMF.metaprograms$metaprograms.metrics,10)

# nmf_metrics <- geneNMF.metaprograms$metaprograms.metrics
# write.csv(nmf_metrics, 
#           file = "./geneNMF/geneNMF_metaprograms_metrics.csv", 
#           row.names = TRUE,  # 不保存行名
#           quote = FALSE)

### MP基因富集分析
library(msigdbr)
library(fgsea)
library(data.table)

# 修改GeneNMF的函数，输出fouldenrinchment
runGSEA_mod <- function(genes, universe=NULL,
                    category="H", subcategory=NULL,
                    species="Homo sapiens",
                    custom.db=NULL,
                    pval.thr=0.05) {
  
  if (!requireNamespace("fgsea", quietly = TRUE) |
      !requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Function 'runGSEA' requires the 'fgsea' and 'msigdbr' packages.
         Please install them.", call. = FALSE)
  }  
  
  if (any(duplicated(genes))) {
    genes <- genes[!duplicated(genes)]
  }
  
  if (!is.null(custom.db)) {
    if (!is.list(custom.db) || is.null(names(custom.db))) {
      stop("custom.db should be a named list", call. = FALSE)
    }
    DB_list <- custom.db
  } else {
    msig_df <- msigdbr::msigdbr(
      species = species,
      collection = category,
      subcollection = subcategory
    )
    DB_list <- split(x = msig_df$gene_symbol, f = msig_df$gs_name)
  }
  
  # 如果 universe 没给，默认用所有通路基因作为背景
  if (is.null(universe)) {
    universe <- unique(unlist(DB_list))
  } else {
    universe <- unique(universe)
  }
  
  # 去掉不在背景中的输入基因
  genes <- unique(genes[genes %in% universe])
  
  fgRes <- fgsea::fora(
    pathways = DB_list,
    genes = genes,
    universe = universe
  )
  
  # 增加期望重叠数和富集倍数
  fgRes$expected_overlap <- length(genes) * fgRes$size / length(universe)
  fgRes$foldEnrichment <- fgRes$overlap / fgRes$expected_overlap
  
  
  fgRes <- fgRes[fgRes$pval <= pval.thr, ]
  fgRes <- fgRes[order(fgRes$pval, -fgRes$foldEnrichment), ]
  
  return(fgRes)
}

gmt_file <- "D:/gmt/h.all.v2026.1.Hs.symbols.gmt"
hallmark_list <- fgsea::gmtPathways(gmt_file)
top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
  runGSEA_mod(genes = program, 
          universe = rownames(seu), 
          custom.db = hallmark_list,   # 关键：传入本地列表
          pval.thr = 0.05)
})

gsea_df <- rbindlist(
  top_p,
  idcol = "MetaProgram"
)
gsea_df$overlapGenes <- vapply(
  gsea_df$overlapGenes,
  function(x) paste(x, collapse = ","),
  character(1)
)

write.csv(
  gsea_df,
  file = "./geneNMF/geneNMF_metaprograms_pathway_H.csv",
  row.names = FALSE,
  quote = FALSE
)

get_top5_pathway <- function(gsea_res, mp_name) {
  # 过滤显著通路（可选，p.adjust < 0.05），无显著则取top5
  sig_res <- gsea_res %>% 
    filter(padj < 0.05) %>% # 可选：仅保留显著通路
    arrange(desc(abs(foldEnrichment))) # foldEnrichment降序（优先核心通路）
  # 若显著通路不足5个，补充所有显著通路，再用非显著补到5个
  if (nrow(sig_res) < 5) {
    non_sig <- gsea_res %>% 
      filter(!pathway %in% sig_res$pathway) %>% 
      arrange(desc(abs(foldEnrichment))) %>% 
      head(5 - nrow(sig_res))
    top5 <- bind_rows(sig_res, non_sig)
  } else {
    top5 <- head(sig_res, 5)
  }
  # 添加MP名称列
  top5$MetaProgram <- mp_name
  return(top5)
}

# 2. 遍历所有MP，提取top5通路并整合
# 给top_p命名（MP1/MP2/...，或用原始名称）
#names(top_p) <- paste0("MP", 1:length(top_p)) 
# 批量提取top5
top5_all <- map2_dfr(top_p, names(top_p), get_top5_pathway)

# 3. 查看整合后的数据结构
head(top5_all)

top5_all$MetaProgram <- factor(
  top5_all$MetaProgram,
  levels = paste0("MP", 1:6)
)

top5_all <- top5_all %>%
  arrange(MetaProgram, desc(foldEnrichment)) %>%   # 组内排序标准可换
  mutate(pathway = factor(pathway, levels = unique(pathway)))

p <- ggplot(
  top5_all,
  aes(
    x = MetaProgram,
    y = pathway,
    size = overlap,      # 或 size
    color = -log10(padj)
  )
) +
  geom_point(alpha = 0.9) +
  scale_size(range = c(2, 8)) +
  scale_color_gradient(
    low = "#08306B",
    high = "red",
    name = "padj"
  ) +
  labs(
    x = "MetaProgram",
    y = "Pathway",
    title = "Top Pathways of Each MetaProgram (GSEA)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

p
ggsave(plot = p,"./geneNMF/dotplot_metaprograms_pathway_H.png",width = 10,height = 10,dpi = 300)
ggsave(plot = p,"./geneNMF/dotplot_metaprograms_pathway_H.pdf",width = 10,height = 10,dpi = 300)


## KEGG
gmt_file <- "D:/gmt/c2.cp.kegg_medicus.v2026.1.Hs.symbols.gmt"
kegg_list <- fgsea::gmtPathways(gmt_file)
top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
  runGSEA_mod(program, universe=rownames(seu), custom.db = kegg_list)
})
top_p_kegg <- lapply(top_p, function(gsea_res) {
  # 转换为数据框（确保格式统一）
  res_df <- as.data.frame(gsea_res)
  # 筛选通路名包含"KEGG_"的行（MSigDB中KEGG通路均以KEGG_开头）
  kegg_res <- res_df[grepl("^KEGG_", res_df$pathway) | grepl("^KEGG_", res_df$pathway), ]
  return(kegg_res)
})

gsea_df <- rbindlist(
  top_p_kegg,
  idcol = "MetaProgram"
)
gsea_df$overlapGenes <- vapply(
  gsea_df$overlapGenes,
  function(x) paste(x, collapse = ","),
  character(1)
)

write.csv(
  gsea_df,
  file = "./geneNMF/geneNMF_metaprograms_pathway_KEGG.csv",
  row.names = FALSE,
  quote = FALSE
)

get_top5_pathway <- function(gsea_res, mp_name) {
  # 过滤显著通路（可选，p.adjust < 0.05），无显著则取top5
  sig_res <- gsea_res %>% 
    filter(padj < 0.05) %>% # 可选：仅保留显著通路
    arrange(desc(abs(foldEnrichment))) # foldEnrichment降序（优先核心通路）
  # 若显著通路不足5个，补充所有显著通路，再用非显著补到5个
  if (nrow(sig_res) < 5) {
    non_sig <- gsea_res %>% 
      filter(!pathway %in% sig_res$pathway) %>% 
      arrange(desc(abs(foldEnrichment))) %>% 
      head(5 - nrow(sig_res))
    top5 <- bind_rows(sig_res, non_sig)
  } else {
    top5 <- head(sig_res, 5)
  }
  # 添加MP名称列
  top5$MetaProgram <- mp_name
  return(top5)
}

# 2. 遍历所有MP，提取top5通路并整合
# 给top_p命名（MP1/MP2/...，或用原始名称）
#names(top_p) <- paste0("MP", 1:length(top_p)) 
# 批量提取top5
top5_all <- map2_dfr(top_p_kegg, names(top_p_kegg), get_top5_pathway)

# 3. 查看整合后的数据结构
head(top5_all)

top5_all$MetaProgram <- factor(
  top5_all$MetaProgram,
  levels = paste0("MP", 1:10)
)

top5_all <- top5_all %>%
  arrange(MetaProgram, desc(foldEnrichment)) %>%   # 组内排序标准可换
  mutate(pathway = factor(pathway, levels = unique(pathway)))


p <- ggplot(
  top5_all,
  aes(
    x = MetaProgram,
    y = pathway,
    size = overlap,      # 或 size
    color = -log10(padj)
  )
) +
  geom_point(alpha = 0.9) +
  scale_size(range = c(2, 8)) +
  scale_color_gradient(
    low = "#08306B",
    high = "red",
    name = "padj"
  ) +
  labs(
    x = "MetaProgram",
    y = "Pathway",
    title = "Top Pathways of Each MetaProgram (GSEA)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

p
ggsave(plot = p,"./geneNMF/dotplot_metaprograms_pathway_kegg.png",width = 15,height = 10,dpi = 300)
ggsave(plot = p,"./geneNMF/dotplot_metaprograms_pathway_kegg.pdf",width = 15,height = 10,dpi = 300)


## GO
gmt_file <- "D:/gmt/c5.go.v2026.1.Hs.symbols.gmt"
go_list <- fgsea::gmtPathways(gmt_file)
top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
  runGSEA_mod(program, universe=rownames(seu), custom.db = go_list)
})
top_p_GO <- lapply(top_p, function(gsea_res) {
  # 转换为数据框（确保格式统一）
  res_df <- as.data.frame(gsea_res)
  # 筛选通路名包含"KEGG_"的行（MSigDB中KEGG通路均以KEGG_开头）
  kegg_res <- res_df[grepl("^GO", res_df$pathway) | grepl("^GO", res_df$pathway), ]
  return(kegg_res)
})


library(data.table)

gsea_df <- rbindlist(
  top_p_GO,
  idcol = "MetaProgram"
)
gsea_df$overlapGenes <- vapply(
  gsea_df$overlapGenes,
  function(x) paste(x, collapse = ","),
  character(1)
)

# immune_keywords <- c("IMMUNE", "INFLAMMATORY", "CYTOKINE", "CHEMOKINE", 
#                      "INTERFERON", "ANTIGEN", "MHC")
# immune_results <- gsea_df %>%
#   filter(grepl(paste(immune_keywords, collapse = "|"), 
#                pathway, ignore.case = TRUE))
write.csv(
  gsea_df,
  file = "./geneNMF/geneNMF_metaprograms_pathway_GO.csv",
  row.names = FALSE,
  quote = FALSE
)

# write.csv(
#   immune_results,
#   file = "D:/项目/AE2511046228-02/06.Mission/geneNMF/geneNMF_metaprograms_pathway_GO.csv",
#   row.names = FALSE,
#   quote = FALSE
# )


get_top5_pathway <- function(gsea_res, mp_name) {
  # 过滤显著通路（可选，p.adjust < 0.05），无显著则取top5
  sig_res <- gsea_res %>% 
    filter(padj < 0.05) %>% # 可选：仅保留显著通路
    arrange(desc(abs(foldEnrichment))) # foldEnrichment降序（优先核心通路）
  # 若显著通路不足5个，补充所有显著通路，再用非显著补到5个
  if (nrow(sig_res) < 5) {
    non_sig <- gsea_res %>% 
      filter(!pathway %in% sig_res$pathway) %>% 
      arrange(desc(abs(foldEnrichment))) %>% 
      head(5 - nrow(sig_res))
    top5 <- bind_rows(sig_res, non_sig)
  } else {
    top5 <- head(sig_res, 5)
  }
  # 添加MP名称列
  top5$MetaProgram <- mp_name
  return(top5)
}

# 2. 遍历所有MP，提取top5通路并整合
# 给top_p命名（MP1/MP2/...，或用原始名称）
#names(top_p) <- paste0("MP", 1:length(top_p)) 
# 批量提取top5
top5_all <- map2_dfr(top_p_GO, names(top_p_GO), get_top5_pathway)

# 3. 查看整合后的数据结构
head(top5_all)

top5_all$MetaProgram <- factor(
  top5_all$MetaProgram,
  levels = paste0("MP", 1:10)
)

top5_all <- top5_all %>%
  arrange(MetaProgram, desc(foldEnrichment)) %>%   # 组内排序标准可换
  mutate(pathway = factor(pathway, levels = unique(pathway)))

# immune_results$abs_foldEnrich <- abs(immune_results$foldEnrichment)
# # 提取MP1-5中每个组abs最大的5行
# top5_by_mp <- immune_results %>%
#   mutate(abs_foldEnrich = abs(foldEnrichment)) %>%
#   filter(MetaProgram %in% paste0("MP", 1:6)) %>%
#   group_by(MetaProgram) %>%
#   slice_max(abs_foldEnrich, n = 5, with_ties = FALSE) %>%
#   ungroup() %>%
#   dplyr::select(-abs_foldEnrich)

p <- ggplot(
  top5_all,
  aes(
    x = MetaProgram,
    y = pathway,
    size = overlap,      # 或 size
    color = -log10(padj)
  )
) +
  geom_point(alpha = 0.9) +
  scale_size(range = c(2, 8)) +
  scale_color_gradient(
    low = "#08306B",
    high = "red",
    name = "padj"
  ) +
  labs(
    x = "MetaProgram",
    y = "Pathway",
    title = "Top Pathways of Each MetaProgram (GSEA)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

p
ggsave(plot = p,"./geneNMF/dotplot_metaprograms_pathway_GO.png",width = 15,height = 10,dpi = 300)
ggsave(plot = p,"./geneNMF/dotplot_metaprograms_pathway_GO.pdf",width = 15,height = 10,dpi = 300)



#### Signature scores for gene programs

# mp.genes <- geneNMF.metaprograms$metaprograms.genes
# # geneNMF.metaprograms$metaprograms.genes.weights$MP7
# seu <- AddModuleScore_UCell(seu, features = mp.genes, ncores=4, name = "")
# mp_cols <- paste0("MP", 1:6)
# seu$MP_group <- colnames(seu@meta.data[, mp_cols])[max.col(seu@meta.data[, mp_cols], ties.method = "first")]
# dir.create("./geneNMF/SC/")
# metadata <- seu@meta.data
# write.csv(metadata,"./geneNMF/SC/metaprograms.csv")

# p1 <- FeaturePlot(seu, features = names(mp.genes),ncol=3)&
#   scale_color_viridis(option="B") &
#   theme(aspect.ratio = 1, axis.text=element_blank(), axis.ticks=element_blank())
# ggsave(plot = p1,"./geneNMF/SC/UMAP_metaprograms.png",width = 16,height = 9,dpi = 300)
# ggsave(plot = p1,"./geneNMF/SC/UMAP_metaprograms.pdf",width = 16,height = 9,dpi = 300)


# genes_to_plot <- names(mp.genes)  
# plot_df <- FetchData(seu, vars = c(genes_to_plot, "cell_type"))

# plot_long <- pivot_longer(plot_df, 
#                           cols = all_of(genes_to_plot),
#                           names_to = "Gene", 
#                           values_to = "Expression")

# p2 <- ggplot(plot_long, aes(x = cell_type, y = Expression, fill = cell_type)) +
#   geom_violin(trim = FALSE, scale = "width") +
#   geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
#   facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
#   theme_minimal() +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(title = NULL)
# p2
# ggsave(plot = p2,"./geneNMF/SC/vlnplot_metaprograms_celltype.png",width = 16,height = 9,dpi = 300)
# ggsave(plot = p2,"./geneNMF/SC/vlnplot_metaprograms_celltype.pdf",width = 16,height = 9,dpi = 300)

# plot_df <- FetchData(seu, vars = c(genes_to_plot, "group"))
# # 转换为长格式
# plot_long <- plot_df %>%
#   pivot_longer(cols = all_of(genes_to_plot),
#                names_to = "Gene",
#                values_to = "Expression")

# # 绘制小提琴图（按基因分面）
# p3 <- ggplot(plot_long, aes(x = group, y = Expression, fill = group)) +
#   geom_violin(trim = FALSE, scale = "width") +
#   geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
#   facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
#   theme_minimal() +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1))
# p3
# ggsave(plot = p3,"./geneNMF/SC/vlnplot_metaprograms_group.png",width = 16,height = 9,dpi = 300)
# ggsave(plot = p3,"./geneNMF/SC/vlnplot_metaprograms_group.pdf",width = 16,height = 9,dpi = 300)

# plot_df <- FetchData(seu, vars = c(genes_to_plot, "sample"))
# # 转换为长格式
# plot_long <- plot_df %>%
#   pivot_longer(cols = all_of(genes_to_plot),
#                names_to = "Gene",
#                values_to = "Expression")

# # 绘制小提琴图（按基因分面）
# p4 <- ggplot(plot_long, aes(x = sample, y = Expression, fill = sample)) +
#   geom_violin(trim = FALSE, scale = "width") +
#   geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
#   facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
#   theme_minimal() +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1))
# p4
# ggsave(plot = p4,"./geneNMF/SC/vlnplot_metaprograms_sample.png",width = 21,height = 9,dpi = 300)
# ggsave(plot = p4,"./geneNMF/SC//vlnplot_metaprograms_sample.pdf",width = 21,height = 9,dpi = 300)

# # p5 <- VlnPlot(seu, features=names(mp.genes), group.by = "leiden",
# #               pt.size = 0, ncol=5) & theme(axis.title = element_blank())
# # p5
# # ggsave(plot = p5,"./geneNMF/vlnplot_metaprograms_sample.png",width = 16,height = 9,dpi = 300)
# # ggsave(plot = p5,"./geneNMF/vlnplot_metaprograms_sample.pdf",width = 16,height = 9,dpi = 300)

# save.image("./geneNMF_R.Rdata")

# p <- DimPlot(
#   seu,
#   group.by = "MP_group",
#   reduction = "umap",
#   pt.size = 0.2
# ) +
#   ggtitle("Cell Type") +
#   theme_classic() +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 16),
#     axis.title = element_text(size = 12, color = "black"),
#     axis.text = element_blank(),         # 和例图一样，不显示坐标刻度数字
#     axis.ticks = element_blank(),
#     axis.line = element_line(color = "black"),
#     legend.title = element_blank(),
#     legend.text = element_text(size = 11),
#     legend.position = "right",
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
#   ) +
#   labs(x = "UMAP1", y = "UMAP2")

# p

# ggsave(plot = p,"./geneNMF/SC/dimplot_metaprograms.png",width = 7,height = 6.5,dpi = 300)
# ggsave(plot = p,"./geneNMF/SC//dimplot_metaprograms.pdf",width = 7,height = 6.5,dpi = 300)





# # matrix <- seu@meta.data[,names(mp.genes)]
# # 
# # #dimred <- scale(matrix)
# # dimred <- as.matrix(matrix)
# # 
# # colnames(dimred) <- paste0("MP_",seq(1, ncol(dimred)))
# # #New dim reduction
# # seu@reductions[["MPsignatures"]] <- new("DimReduc",
# #                                         cell.embeddings = dimred,
# #                                         assay.used = "RNA",
# #                                         key = "MP_",
# #                                         global = FALSE)
# # 
# # set.seed(123)
# # seu <- RunUMAP(seu, reduction="MPsignatures", dims=1:length(seu@reductions[["MPsignatures"]]),
# #                metric = "euclidean", reduction.name = "umap_MP")
# # 
# # DimPlot(seu, group.by = c('sample', 'group'), reduction = 'umap_MP')
# # 
# # 
# # 
# # set.seed(123)
# # seu <- 
# #   seu %>%
# #   FindNeighbors(reduction = 'MPsignatures') %>%
# #   FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5))
# # 
# # DimPlot(seu, group.by = c('RNA_snn_res.0.1'), reduction = 'umap_MP')
# # 
# # 
# # 
# # umap<-adata@reductions$umap@cell.embeddings
#   # write.csv(umap,"D:/项目/AE2511046228-02/03.Mission/CD4T/umapreductions.csv")