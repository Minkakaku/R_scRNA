# === 头部参数接收修改开始 ===
args <- commandArgs(trailingOnly = TRUE)
# 动态接收参数
case <- args[1]     # 实验组 (例如 YHY)
control <- args[2]  # 对照组 (例如 PDQ)
specie <- args[3]   # 物种 (human 或 mouse)
qsdata <- args[4]   # qs 数据文件路径
cpu <- as.numeric(args[5]) # CPU 核心数

source("/mnt/NexData_file/home/lyh/Script/R/source.R")
if (!exists("output_dir")) {output_dir <- getwd()}
setup_results_dir(output_dir,dir_name = "DE_Analysis")
pbmc_filt <- qread(qsdata)

## Data
library(scRNAtoolVis)
library(Seurat)

# parameter
colors <- color3
# case <- "PmMET"
# control <- "RT"
# specie <- "human"  # 选项: "human" 或 "mouse" (根据你的数据修改这里!)

# === 2. 物种配置 ===
if (tolower(specie) == "human") {
  library(org.Hs.eg.db);org_db <- org.Hs.eg.db;kegg_org <- "hsa";pathway_specie="Homo_sapiens";org_pkg <- "org.Hs.eg.db"
} else if (tolower(specie) == "mouse") {
  library(org.Mm.eg.db);org_db <- org.Mm.eg.db;kegg_org <- "mmu";pathway_specie="Mus_musculus";org_pkg <- "org.Mm.eg.db"
} else {stop("specie 参数错误，请设置为 'human' 或 'mouse'")}

# 取子集
print(paste("Subsetting data for groups:", case, "and", control))
pbmc_filt$Group %>% unique() %>% print_vec()
pbmc_filt <- subset(pbmc_filt,subset = Group %in% c(case, control))


# PARTⅠ: DEG
## Analysis
save_dir <- file.path(output_dir, "DE_Analysis", paste0(case, "_vs_", control, "_DEG"))
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
setwd(save_dir)

DefaultAssay(pbmc_filt)  <- "RNA"
library(future)
plan("multicore", workers = cpu)
options(future.globals.maxSize = 100 * 1024^3)

pbmc_filt@meta.data %>% with(table(Celltype,Group))

Idents(pbmc_filt) <- "Group"
DE_DF <- FindMarkers(pbmc_filt, ident.1 = case, ident.2 = control, logfc.threshold = 0, min.pct = 0, only.pos = FALSE) %>%
  rownames_to_column("GeneName") %>% 
  mutate(Celltype = "all_cell")
for(ct in unique(pbmc_filt$Celltype)){
  sub_obj <- subset(pbmc_filt, subset = Celltype == ct)
  if(sum(sub_obj$Group == case) < 3 | sum(sub_obj$Group == control) < 3) next
  marker_DF <- FindMarkers(sub_obj, ident.1 = case, ident.2 = control, logfc.threshold = 0, min.pct = 0, only.pos = FALSE) %>%
    rownames_to_column("GeneName") %>% 
    mutate(Celltype = ct)
  DE_DF <- bind_rows(DE_DF, marker_DF)
}
DE_DF <- DE_DF %>% dplyr::rename(Pvalue=p_val,Qvalue=p_val_adj,log2FC=avg_log2FC)
openxlsx::write.xlsx(DE_DF,file = paste0(case, "vs", control, "_DEG_inEachCelltype.xlsx"),rowNames = FALSE)

if (!exists("DE_DF")) {DE_DF <- read.xlsx(paste0(case, "vs", control, "_DEG_inEachCelltype.xlsx"))}
DE_DF %>% str

## Visualization(DEG)
### Volcano
markers <- DE_DF %>% filter(Celltype!="all_cell") %>% dplyr::rename(gene=GeneName,avg_log2FC=log2FC,p_val=Pvalue,cluster=Celltype,p_val_adj=Qvalue)

n_clust <- length(unique(markers$cluster))
width <- max(6, min(n_clust, 8) * 1.5 + 1);height <- 9
fig.size(width,height)
jjVolcano(markers, tile.col=colors,size = 3,celltypeSize = 3.5,legend.position = c(0.9, 1))
ggsave1("all_cluster_volcano",width=width,height=height)

n_clust <- length(unique(markers$cluster));calc_nrow <- ceiling(n_clust / 10)
width <- max(6, min(n_clust, 8) * 1.5 + 1);height <- calc_nrow * 4 + 2
fig.size(width,height)
markerVolcano(markers,labelCol = colors,nrow=calc_nrow,pforce=3)
ggsave1("all_cluster_volcano_pct", width = width, height = height)

library(dplyr); library(ggplot2); library(ggrepel); library(ggnewscale)

# === 预定义主题 ===
my_theme <- theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(), legend.position = 'right',
        legend.justification = c(0, 1), legend.key.size = unit(0.4, 'cm'),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(face = "bold", size = 14),
        panel.border = element_rect(linewidth = 1))

# === 循环绘图 ===
for(ct in unique(DE_DF$Celltype)){
  message(paste("Processing:", ct))
  
  # 1. 数据准备与分层
  sub_data <- DE_DF %>% filter(Celltype == ct) %>%
    mutate(group = case_when(Pvalue < 0.01 & log2FC > 0.25 ~ "Up",
                             Pvalue < 0.01 & log2FC < -0.25 ~ "Down",
                             TRUE ~ "No"))
  
  # 2. 筛选 Top 10 标注基因 (Top 5 Up + Top 5 Down by Pvalue)
  candidates <- sub_data %>% 
    filter(group != "No", !grepl("^RPS|^RPL|^MT-", GeneName)) %>%
    arrange(Pvalue)
  
  top_genes <- c(head(filter(candidates, group == "Up")$GeneName, 10),
                 head(filter(candidates, group == "Down")$GeneName, 10))
  
  plot_data <- sub_data %>% mutate(is_label = ifelse(GeneName %in% top_genes, "Yes", "No"))
  
  # 3. 拆分图层数据
  bg_no   <- filter(plot_data, group == "No")
  bg_down <- filter(plot_data, group == "Down", is_label == "No")
  bg_up   <- filter(plot_data, group == "Up",   is_label == "No")
  lab_down<- filter(plot_data, group == "Down", is_label == "Yes")
  lab_up  <- filter(plot_data, group == "Up",   is_label == "Yes")
  lab_all <- filter(plot_data, is_label == "Yes")
  
  # 4. 动态坐标轴范围
  lim_max <- max(abs(plot_data$log2FC), na.rm = TRUE) * 1.1
  
  # 5. 绘图
  p <- ggplot() +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "gray60", linewidth = 0.4) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "gray60", linewidth = 0.4) +
    
    # 背景灰点
    geom_point(data = bg_no, aes(x = log2FC, y = -log10(Pvalue)), color = "grey85", size = 1, alpha = 0.5) +
    
    # 蓝色系 (Down)
    geom_point(data = bg_down, aes(x = log2FC, y = -log10(Pvalue), color = log2FC), size = 2, alpha = 0.5) +
    geom_text_repel(data = lab_down, aes(x = log2FC, y = -log10(Pvalue), label = GeneName, color = log2FC),
                    size = 3.5, fontface = "bold.italic", bg.color = "white", bg.r = 0.15,
                    max.overlaps = Inf, seed = 123, segment.size = 0.3, segment.color = "black") +
    scale_color_gradient(name = "Down", low = "#08306B", high = "#4292C6", guide = guide_colorbar(order = 1)) +
    
    new_scale_color() +
    
    # 红色系 (Up)
    geom_point(data = bg_up, aes(x = log2FC, y = -log10(Pvalue), color = log2FC), size = 2, alpha = 0.5) +
    geom_text_repel(data = lab_up, aes(x = log2FC, y = -log10(Pvalue), label = GeneName, color = log2FC),
                    size = 3.5, fontface = "bold.italic", bg.color = "white", bg.r = 0.15,
                    max.overlaps = Inf, seed = 123, segment.size = 0.3, segment.color = "black") +
    scale_color_gradient(name = "Up", low = "#FC9272", high = "#CB181D", guide = guide_colorbar(order = 2)) +
    
    new_scale_fill() +
    
    # 高亮金点
    geom_point(data = lab_all, aes(x = log2FC, y = -log10(Pvalue), fill = abs(log2FC)),
               shape = 21, color = "black", stroke = 0.6, size = 3.5) +
    scale_fill_gradient(name = "Highlighted", low = "#FFEB3B", high = "#FF9800", guide = guide_colorbar(order = 3)) +
    
    my_theme +
    labs(title = ct, x = expression(log[2] ~ "Fold Change"), y = expression(-log[10] ~ italic(P) ~ "-value")) +
    scale_x_continuous(limits = c(-lim_max, lim_max))
  
  # 6. 保存 (文件名处理了特殊字符)
  safe_name <- gsub("[/|:?]", "_", ct)
  ggsave1(paste0("Volcano_", safe_name), plot = p, width = 6, height = 4.5)
}

### DEG_Number
DE_DF %>% str

deg <- DE_DF %>% filter(Pvalue < 0.05, abs(log2FC) > 0.25) %>% filter(Celltype != "all_cell") %>%
  mutate(Direction = ifelse(log2FC > 0, "Up", "Down")) %>%group_by(Celltype, Direction) %>% summarise(Num = n(), .groups = 'drop')
# 指定因子水平，确保 Up 在上或在特定顺序（可选）
deg$Direction <- factor(deg$Direction, levels = c("Up", "Down"))
deg %>% write.csv("DEG_Number.csv",row.names=F)
# 数据转换：将 Down 的数值改为负数
deg_mirror <- deg %>% mutate(Num_plot = ifelse(Direction == "Down", -Num, Num))
fig.size(5,6)
# 镜像作图
ggplot(deg_mirror, aes(x = Celltype, y = Num_plot, fill = Direction)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_y_continuous(labels = abs) + 
  scale_fill_manual(values = c("Up" = "#BC3C29", "Down" = "#0072B5")) + 
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) + # 添加 0 轴线
  labs(x = NULL, y = "Number of DEGs", title = "DEG Distribution") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(), # 去除背景网格
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5)
  ) +
  geom_text(aes(label = Num), position = position_stack(vjust = 0.5), size = 3, color = "white") # 添加数值标签
ggsave1("Number_of_DEGs",width=5,height=6)

### Heatmap
library(grid)
library(scop)

if (!exists("DE_DF")) {DE_DF <- read.xlsx(paste0(case, "vs", control, "_DEG_inEachCelltype.xlsx"))}
deg <- DE_DF %>% filter(Pvalue < 0.05, abs(log2FC) > 0.25) %>% filter(Celltype != "all_cell") %>% 
    group_by(Celltype) %>% slice_max(order_by=log2FC,n=5) %>% ungroup()

# scop 分组+细胞类切面热图
target_genes <- deg$GeneName
# 核心：利用 deg 中的 Celltype 自动生成切分向量
split_vec <- factor(deg$Celltype, levels = unique(deg$Celltype))
# 动态计算尺寸 (根据基因数和细胞类型数自适应)
base_w <- 4; base_h <- 2 
n_rows <- length(target_genes)
# 注意：若无 Group 分组需求，n_cols 计算可简化
n_cols <- length(unique(pbmc_filt$Celltype)) * length(unique(pbmc_filt$Group)) 
calc_w <- base_w + (n_cols * 0.2)
calc_h <- base_h + (n_rows * 0.15) # 稍微增加行高权重防止基因名重叠
p <- scop::GroupHeatmap(
  srt = pbmc_filt, 
  features = target_genes, 
  group.by = "Celltype",
  split.by = "Group", # 如果不需要按样本分列，可注释此行
  nlabel = 0,
  # --- 分区设置 (自动化) ---
  feature_split = split_vec,      # 传入自动构建的向量
  cluster_row_slices = FALSE,     # 关键：禁止重排板块，保持 deg 顺序
  row_title_rot = 0,              # 标题水平便于阅读
  feature_split_palcolor = colors,#左侧分面条颜色
  # --- 基础与美化 ---
  show_row_names = TRUE,ow_names_side = 'right',
  cluster_rows = FALSE,cluster_columns = FALSE,
  cell_split_palette = "Dark2",
  group_palcolor = list(colors), 
  cell_split_palcolor = group_color,
  # --- 细节调整 ---
  ht_params = list(
    row_gap = unit(2, "mm"),
    row_names_gp = gpar(fontsize = 10, fontface = "italic"), 
    row_title_gp = gpar(fontsize = 11, fontface = "bold"),
    rect_gp = gpar(col = "white", lwd = 0.5) #以此增加格线区分度
  )
)
ggsave1('DEG_top5_Heatmap', calc_w, calc_h, plot = p$plot)


# PARTⅡ: Enrichment
## Analysis
library(clusterProfiler)
# 定义富集分析结果的主目录及子目录
setup_results_dir(output_dir,dir_name = "Enrichment")
save_dir <- file.path(output_dir, "Enrichment", paste0(case, "_vs_", control, "_Enrichment"))
dir.create(file.path(save_dir, "RDS"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(save_dir, "Tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(save_dir, "Enrichment_Plots"), recursive = TRUE, showWarnings = FALSE)

# === 准备差异基因数据 ===
if (!exists("DE_DF")) {DE_DF <- read.xlsx(paste0(case, "vs", control, "_DEG_inEachCelltype.xlsx"))}
sig_deg <- DE_DF %>% filter(Pvalue < 0.05,abs(log2FC) > 0.25) %>% dplyr::rename(p_val=Pvalue,p_val_adj=Qvalue,avg_log2FC=log2FC)
run_species <- if(exists("specie") && tolower(specie) == "mouse") "Mus_musculus" else "Homo_sapiens" # 确定物种参数
celltypes <- unique(sig_deg$Celltype)

library(parallel)
library(doParallel)
cl <- makeCluster(cpu)
registerDoParallel(cl)

message(paste0(">>> 启动并行计算 (Cores: ", cpu, ")..."))

export_vars <- c("ggsave1", "kegg_org", "save_dir", "sig_deg", "org_pkg", "specie")

all_enrich_results_list <- foreach(ct = celltypes, 
                                   .packages = c("dplyr", "clusterProfiler", "ggplot2", "stringr", "openxlsx"),
                                   .export = export_vars) %dopar% {
  
  library(org_pkg, character.only = TRUE)
  org_db <- get(org_pkg) 
  
  safe_ct <- gsub("[^[:alnum:]]", "_", ct)
  genes <- sig_deg$GeneName[sig_deg$Celltype == ct]
  ct_res_list <- list()
  
  # --- 步骤 1: ID 转换 ---
  eg_df <- clusterProfiler::bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org_db)
  eg <- eg_df$ENTREZID
  
  if(length(eg) == 0) {
    stop(paste0("细胞类型 ", ct, " ID转换后没有得到有效的 Entrez ID"))
  }

  # --- 步骤 2: 定义分析任务 ---
  tasks <- list(
    GO_BP = list(fun=enrichGO, args=list(ont="BP", OrgDb=org_db, keyType="ENTREZID")),
    GO_CC = list(fun=enrichGO, args=list(ont="CC", OrgDb=org_db, keyType="ENTREZID")),
    GO_MF = list(fun=enrichGO, args=list(ont="MF", OrgDb=org_db, keyType="ENTREZID")),
    KEGG  = list(fun=enrichKEGG, args=list(organism=kegg_org))
  )
  
  # --- 步骤 3: 循环运行 ---
  for (db in names(tasks)) {
    res <- do.call(tasks[[db]]$fun, c(list(gene=eg, pvalueCutoff=0.05, qvalueCutoff=0.2), tasks[[db]]$args))
    
    if (!is.null(res) && nrow(res@result) > 0) {
      res <- setReadable(res, OrgDb=org_db, keyType="ENTREZID")
      ct_res_list[[db]] <- res
      
      saveRDS(res, file.path(save_dir, "RDS", paste0(safe_ct, "_", db, ".rds")))
      write.csv(filter(as.data.frame(res), p.adjust < 0.05), 
                file.path(save_dir, "Tables", paste0(safe_ct, "_", db, ".csv")), row.names=FALSE)

      p1 <- dotplot(res, showCategory=20) + ggtitle(paste(ct, "-", db)) + theme(axis.text.y=element_text(size=10))
      ggsave1(file.path(save_dir, "Enrichment_Plots", paste0(safe_ct, "_", db, "_Dot")), plot=p1, width=8, height=7)
    }
  }
  
  return(ct_res_list)
}

stopCluster(cl)

# === 5. 结果整理 ===
names(all_enrich_results_list) <- celltypes
all_enrich_results <- all_enrich_results_list
saveRDS(all_enrich_results, file.path(save_dir, "RDS", "All_Results.rds"))

message(paste0("分析完成，结果已归档至: ", save_dir))

## Visualization(Enrichment)
### scop plot
# === 将 clusterProfiler 结果转换为 scop 格式 ===
extract_enrich_df <- function(res_obj, celltype, database) {
  if (is.null(res_obj)) return(NULL)
  df <- tryCatch(as.data.frame(res_obj), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df$ID <- as.character(df$ID)
  df$Description <- as.character(df$Description)
  df$geneID <- as.character(df$geneID)
  df$Groups <- celltype
  df$Database <- database 
  return(df)
}
df_list <- list()
# 遍历 names(all_enrich_results) 确保取到名字
for (ct in names(all_enrich_results)) {
  ct_data <- all_enrich_results[[ct]]
  for (db in names(ct_data)) {
    temp_df <- extract_enrich_df(ct_data[[db]], celltype = ct, database = db)
    if (!is.null(temp_df)) {
      df_list[[length(df_list) + 1]] <- temp_df
    }
  }
}
combined_df <- dplyr::bind_rows(df_list)
# 构造 scop 格式对象
scop_formatted_res <- list(enrichment = combined_df)
message(paste0(">>> 转换完成，合并了 ", nrow(combined_df), " 行富集结果。"))

sel_celltype <- celltypes[1]
sel_id <- scop_formatted_res$enrichment %>% filter(if("Groups" %in% names(.)) Groups == sel_celltype else TRUE) %>% filter(if("Count" %in% names(.)) Count < 50 else TRUE) %>% arrange(p.adjust)  %>% slice_head(n = 4) %>% pull(ID)
# === A. 绘制 Network 图 (展示基因与通路的关系) ===
p_network <- EnrichmentPlot(
  res = scop_formatted_res,plot_type = "network",padjustCutoff = 0.01,group_use = sel_celltype,
  db = c("GO_BP","GO_CC","GO_MF"), topTerm = 3, network_layout = "fr", id_use=sel_id
) + ggtitle(paste0(sel_celltype, " - GO Network"))
ggsave1(file.path(save_dir, "Enrichment_Plots", "Network_Plot"), plot=p_network, width=12, height=12)

# # === B. 绘制 Enrichmap (展示通路间的相似性聚类) ===
# p_enrichmap <- EnrichmentPlot(
#   res = scop_formatted_res,
#   plot_type = "enrichmap",
#   group_use = sel_celltype,     # 指定细胞类型
#   db = c("GO_BP","GO_CC","GO_MF"),
#   topTerm = 30,                 # 聚类通常需要较多通路
#   enrichmap_cluster = "fast_greedy",
#   enrichmap_label = "term"      # 标签显示通路名
# ) + ggtitle(paste0(sel_celltype, " - GO Enrichmap"))

# print(p_enrichmap)
# # ggsave("Enrichmap_Plot.pdf", p_enrichmap, width = 10, height = 8)

# === C. 绘制 Comparison 图 (适合展示所有细胞类型的富集) ===
p_compare <- EnrichmentPlot(res = scop_formatted_res,plot_type = "comparison",compare_only_sig=T,db = c("GO_BP"),topTerm = 5)
p_compare <- EnrichmentPlot(
  res = scop_formatted_res,plot_type = "comparison",compare_only_sig=T,
  db = c("GO_BP"),#,"GO_CC","GO_MF","KEGG"
  topTerm = 5,aspect.ratio=nrow(p_compare$data)/length(celltypes)/length(celltypes)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),)
ggsave1(file.path(save_dir, "Enrichment_Plots", "Comparison"), plot=p_compare, width = length(celltypes)/1.2, height = nrow(p_compare$data)/length(celltypes)/3)

# === D. 绘制 Lollipop 图 (棒棒糖图) ===
p_lollipop <- EnrichmentPlot(
  res = scop_formatted_res,plot_type = "lollipop",group_use = sel_celltype,
  db = c("GO_BP"), topTerm = 10, aspect.ratio=2.5
)
ggsave1(file.path(save_dir, "Enrichment_Plots", "Lollipop"), plot=p_lollipop, width = 7, height = 4)

### combined bar
# 1. 数据准备
prepare_composite_data_v4 <- function(enrich_list, celltype, top_n = 5, gene_width = 77) { # 基因字符宽度
  ct_res <- enrich_list[[celltype]]
  if(is.null(ct_res)) stop(paste("未找到细胞类型:", celltype))
  # 定义顺序
  db_order <- c("GO_BP", "GO_CC", "GO_MF", "KEGG", "DO")
  db_label_map <- c("GO_BP" = "BP", "GO_CC" = "CC", "GO_MF" = "MF", "KEGG" = "KEGG", "DO" = "DO")
  plot_df <- data.frame()
  for (db_name in db_order) {
    if (!is.null(ct_res[[db_name]]) && nrow(ct_res[[db_name]]@result) > 0) {
      res_df <- as.data.frame(ct_res[[db_name]]) %>% filter(p.adjust < 0.05)
      if (nrow(res_df) > 0) {
        top_df <- res_df %>% 
          arrange(p.adjust) %>% 
          head(top_n) %>%
          mutate(
            ONTOLOGY = db_label_map[[db_name]], 
            LogP = -log10(p.adjust),
            # 截断基因字符串
            Genes_Formatted = str_trunc(geneID, width = gene_width, side = "right")
          ) %>%
          select(Description, p.adjust, LogP, Count, Genes_Formatted, ONTOLOGY)
        
        plot_df <- rbind(plot_df, top_df)
      }
    }
  }
  if(nrow(plot_df) == 0) return(NULL)
  # 排序逻辑：先按 ONTOLOGY 分组，再按 LogP 排序
  valid_levels <- intersect(db_label_map[db_order], unique(plot_df$ONTOLOGY))
  plot_df$ONTOLOGY <- factor(plot_df$ONTOLOGY, levels = valid_levels)
  
  plot_df <- plot_df %>% 
    arrange(desc(ONTOLOGY), LogP) %>% 
    mutate(Description = factor(Description, levels = Description))
  return(plot_df)
}
# 辅助函数：计算分类标签的中心位置
get_category_labels <- function(df) {
  df %>% group_by(ONTOLOGY) %>%mutate(rn = row_number(), n = n()) %>%filter(rn == ceiling(n / 2)) %>%select(ONTOLOGY, Description) %>% ungroup()
}

# 绘图主程序
target_ct <- names(all_enrich_results_list)[1]
plot_data <- prepare_composite_data_v4(all_enrich_results_list, target_ct, top_n = 5)
cat_labels <- get_category_labels(plot_data)
my_colors <- c("BP"="#F98400", "CC"="#CD534C", "MF"="#009E73", "KEGG"="#5BBCD6", "DO"="#E64B35")
# --- 距离控制参数 (保持您的紧凑设置) ---
max_p <- max(plot_data$LogP);strip_width <- max_p * 0.04      
strip_x <- -max_p * 0.08;dot_x <- -max_p * 0.03
main_title <- paste0("Enrichment Analysis: ", target_ct)
p <- ggplot(plot_data, aes(y = Description)) +
  # 1. 左侧整根分类 Bar 
  # [修复]: show.legend = FALSE (防止它生成多余的图注)
  geom_tile(aes(x = strip_x, fill = ONTOLOGY), 
            width = strip_width, height = 1, color = NA, show.legend = FALSE) +
  # 2. 左侧分类文字
  # [修复]: show.legend = FALSE (去除图注里的 'a')
  geom_text(data = cat_labels, 
            aes(x = strip_x, label = ONTOLOGY), 
            angle = 90, color = "white", fontface = "bold", size = 4, show.legend = FALSE) +
  # 3. 中间圆点 (负责生成 Count 图注)
  geom_point(aes(x = dot_x, color = ONTOLOGY, size = Count), 
             shape = 21, fill = "white", stroke = 1.2) +
  # 4. 圆点内数字
  geom_text(aes(x = dot_x, label = Count), 
            size = 3, color = "black", show.legend = FALSE) +
  # 5. 右侧数据 Bar (负责生成 Ontology 图注)
  geom_col(aes(x = LogP, fill = ONTOLOGY), width = 0.45, alpha = 0.8) +
  # 6. 通路名称
  geom_text(aes(x = 0.2, label = Description), 
            hjust = 0, vjust = 0.2, size = 4.5, fontface = "bold", color = "black", show.legend = FALSE) +  
  # 7. 基因列表
  geom_text(aes(x = 0.2, label = Genes_Formatted, color = ONTOLOGY), 
            hjust = 0, vjust = 2.5, size = 3.5, fontface = "italic", show.legend = FALSE) +
  # --- 样式与标尺 ---
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  scale_size_continuous(range = c(4, 8)) + # 控制圆点大小范围
  # --- 图注控制 ---
  guides(
    fill = guide_legend(title = "Ontology", order = 1, override.aes = list(alpha = 1)),
    size = guide_legend(title = "Count", order = 2, override.aes = list(color = "black", fill = "white")), 
    color = "none"     # 3. 隐藏多余的线条颜色图注
  ) +
  # X轴设置
  scale_x_continuous(
    expand = expansion(mult = c(0.01, 0.05)), limits = c(strip_x - strip_width/2 - (max_p * 0.02), NA), 
    breaks = seq(0, max_p + 10, by = 5),labels = function(x) ifelse(x < 0, "", x)
  ) +
  facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y", switch = "y") +
  theme_classic() +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), axis.title.y = element_blank(),
    strip.text = element_blank(), strip.background = element_blank(),axis.title.x = element_text(size = 16),
    axis.line.x = element_line(color = "black"), axis.text.x = element_text(color = "black", size = 10),
    plot.title = element_text(hjust = 0.5,face = "bold",size = 18,margin = margin(b = 10)),
    # 图注位置
    legend.position = "right",legend.box = "vertical",legend.title = element_text(face = "bold", size = 10),
    panel.spacing = unit(0.1, "lines")
  ) +
  labs(x = "-log10(p.adjust)",title = main_title)
ggsave1(file.path(save_dir, "Enrichment_Plots", "Complex_combined_bar"), plot=p, width = 10, height = 10)
fig.size(10,10)
print(p)

# PARTⅢ: GSEA
## Analysis
# output_dir <- "/public23/4_After_Sales_Personal/FuChenXuan/0_sc_pipeline/1_methods/DEG_Enrichment"
setup_results_dir(output_dir,dir_name = "Enrichment")
# 设置具体的保存子目录
save_dir <- file.path(output_dir, "Enrichment", paste0(case, "_vs_", control, "_GSEA"))
dir.create(file.path(save_dir, "RDS"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(save_dir, "Tables"), recursive = TRUE, showWarnings = FALSE)

library(clusterProfiler)
library(stringr)
library(GseaVis)
if (!exists("DE_DF")) {DE_DF <- read.xlsx(paste0(case, "vs", control, "_DEG_inEachCelltype.xlsx"))}

library(parallel)
library(doParallel)
library(foreach)

# === 1. 数据准备 ===
if (!exists("DE_DF")) DE_DF <- read.xlsx(paste0(case, "vs", control, "_DEG_inEachCelltype.xlsx"))
if (!exists("org_pkg")) org_pkg <- org_db@packageName 
df_list <- split(DE_DF, DE_DF$Celltype)
celltypes <- names(df_list)

# === 2. 启动集群 (PSOCK模式) ===
cl <- makeCluster(cpu)
registerDoParallel(cl)
message(paste0(">>> 启动并行 GSEA (Cores: ", cpu, ")... | ", case, " vs ", control))
export_vars <- c("df_list", "kegg_org", "save_dir", "case", "control", "org_pkg")

# === 3. 并行循环 ===
results <- foreach(ct = celltypes, 
                   .packages = c("dplyr", "stringr", "openxlsx", "clusterProfiler", "DOSE"),
                   .export = export_vars) %dopar% {
  library(org_pkg, character.only = TRUE)
  org_db_local <- get(org_pkg)
  
  safe_ct <- str_replace_all(ct, "[/ :]", "_")
  prefix  <- paste0(case, "_vs_", control, "_", safe_ct)
  sub_df  <- df_list[[ct]]
  gl_sym <- sort(setNames(sub_df$log2FC, sub_df$GeneName), decreasing = TRUE)
  ids <- tryCatch({
    bitr(names(gl_sym), fromType="SYMBOL", toType="ENTREZID", OrgDb=org_db_local)
  }, error = function(e) NULL)
  gl_ent <- NULL
  if(!is.null(ids)){
    dedup <- ids[!duplicated(ids$SYMBOL), ]
    gl_ent <- sort(setNames(gl_sym[dedup$SYMBOL], dedup$ENTREZID), decreasing = TRUE)
  }
  # --- Step C: 运行 GSEA GO ---
  try({
    go <- gseGO(
      geneList     = gl_sym,
      OrgDb        = org_db_local,
      keyType      = "SYMBOL",
      ont          = "ALL",
      minGSSize    = 10, maxGSSize = 1000, pvalueCutoff = 1,
      verbose      = FALSE, seed = 123
    )
    if(!is.null(go) && nrow(go) > 0){
      saveRDS(go, file.path(save_dir, "RDS", paste0("GSEA_GO_", prefix, ".rds")))
      write.xlsx(go@result, file.path(save_dir, "Tables", paste0("GSEA_GO_", prefix, ".xlsx")))
    }
  })
  # --- Step D: 运行 GSEA KEGG ---
  if(!is.null(gl_ent) && length(gl_ent) > 0){
    try({
      kegg <- gseKEGG(
        geneList     = gl_ent,
        organism     = kegg_org,
        minGSSize    = 10, maxGSSize = 1000, pvalueCutoff = 1,
        verbose      = FALSE, seed = 123
      )
      if(!is.null(kegg) && nrow(kegg) > 0){
        saveRDS(kegg, file.path(save_dir, "RDS", paste0("GSEA_KEGG_", prefix, ".rds")))
        kegg_read <- setReadable(kegg, OrgDb=org_db_local, keyType="ENTREZID")
        write.xlsx(kegg_read@result, file.path(save_dir, "Tables", paste0("GSEA_KEGG_", prefix, ".xlsx")))
      }
    })
  }
  return(ct)
}
# === 4. 结束清理 ===
stopCluster(cl)
message(paste0("所有分析完成，结果已归档至: ", save_dir))

## Visualization(GSEA)
#1.柱状图 top10 ; 2.峰峦图 top10 （enrich_ridge_plot）; 3.GSAEvis top5 （gseaNb）; 4.折线图 top5 （scplotter::GSEAPlot）
# === 1. 定义核心绘图函数 ===
plot_gsea_pipeline <- function(category, rds_dir, root_save_dir) {
  
  message(paste0(">>> 正在处理分类: ", category))
  
  # 路径准备
  files <- list.files(rds_dir, pattern = paste0("GSEA_", category, "_.*\\.rds"))
  cell_types <- str_remove_all(files, paste0("GSEA_", category, "_|\\.rds"))
  plot_root <- file.path(root_save_dir, category)
  if(!dir.exists(plot_root)) dir.create(plot_root, recursive = TRUE)
  
  for(ct in cell_types){
    rds_file <- file.path(rds_dir, paste0("GSEA_", category, "_", ct, ".rds"))
    gsea_obj <- readRDS(rds_file)
    gsea_obj@readable=FALSE
    if(nrow(gsea_obj) == 0) next
    
    # 统一转换 Symbol
    gsea_sym <- tryCatch(setReadable(gsea_obj, OrgDb = org_db, keyType = "ENTREZID"), 
                         error = function(e) gsea_obj)
    
    # 准备输出目录
    out_dir <- file.path(plot_root, ct)
    if(!dir.exists(out_dir)) dir.create(out_dir)
    
    # --- 1. 双向柱状图 (Top 10 Up/Down) ---
    res <- gsea_sym@result
    res$Group <- ifelse(res$NES > 0, "Up", "Down")
    top_up <- res %>% filter(Group == "Up") %>% arrange(p.adjust) %>% head(10)
    top_down <- res %>% filter(Group == "Down") %>% arrange(p.adjust) %>% head(10)
    plot_data <- rbind(top_up, top_down)
    
    if(nrow(plot_data) > 0){
      p1 <- ggplot(plot_data, aes(reorder(Description, NES), NES, fill = Group)) +
        geom_col() + theme_bw() + coord_flip() + geom_hline(yintercept = 0) +
        geom_text(data = plot_data[plot_data$NES > 0,], aes(y = -0.1, label = Description), hjust = 1, size = 3, fontface = "bold") +
        geom_text(data = plot_data[plot_data$NES < 0,], aes(y = 0.1, label = Description), hjust = 0, size = 3, fontface = "bold") +
        labs(x = NULL, y = "NES", title = paste(ct, category, "Top Pathways")) +
        scale_fill_manual(values = c("Up" = "#D73027", "Down" = "#313695")) +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
              panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(color = 'grey60'))
      ggsave1(file.path(out_dir, "1_Barplot"), width = 7, height = 5, plot = p1)
    }
    
    # --- 2. 峰峦图 (Ridge Plot) Top 10 ---
    sig_terms <- c(head(top_up$ID, 5), head(top_down$ID, 5)) # Top 5 Up + Top 5 Down
    
    if(length(sig_terms) > 0){
      ridge_data <- gsea_sym@result %>% filter(ID %in% sig_terms) %>%
        mutate(Description = factor(Description, levels = Description[order(NES)]))
      
      fc_vals <- gsea_obj@geneList
      min_fc <- min(fc_vals, na.rm=TRUE); max_fc <- max(fc_vals, na.rm=TRUE)
      point_x <- min_fc 
      x_lims <- c(point_x - 0.1, max_fc) # 右侧严格去白
      
      p2 <- enrich_ridge_plot(gsea_obj, terms_ID = sig_terms) + # 用 Entrez 对象绘图
        geom_point(data = ridge_data, aes(x = point_x + 0.15, y = Description, color = NES, size = setSize), alpha = 0.85) +
        scale_color_gradient2(name = "NES", low = "blue", mid = "#f7f7f7", high = "red", midpoint = 0) +
        scale_size_continuous(name = "Count", range = c(3, 8)) +
        scale_fill_gradientn(name = "-Log10(p)", colours = c("#FFD700", "#FF8C00", "#DA70D6", "#8A2BE2")) +
        scale_x_continuous(limits = x_lims, expand = expansion(mult = c(0.01, 0))) +
        theme_minimal(base_size = 14) +
        theme(panel.grid = element_blank(), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
              axis.text.y = element_text(size = 12, color = "black", margin = margin(r = 5)),
              plot.title = element_text(hjust = 0.5, face = "bold", size = 16)) +
        labs(x = "Log2 Fold Change", y = NULL, title = paste(ct, "- Ridge"))
      ggsave1(file.path(out_dir, "2_RidgePlot"), width = 10, height = 6, plot = p2)
    }
    
    # --- 3. GseaVis (Top 5 联合) ---
    top_terms_id <- res %>% arrange(desc(abs(NES))) %>% head(5) %>% pull(ID)
    if(length(top_terms_id) > 0){
      cols <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD")[1:length(top_terms_id)]
      p3 <- gseaNb(object = gsea_sym, geneSetID = top_terms_id, newGsea = FALSE, addPval = TRUE, pvalY = 0.85, curveCol = cols)
      ggsave1(file.path(out_dir, "3_GseaVis"), width = 10, height = 6, plot = p3)
    }
    
    # --- 4. scplotter 折线图 (Combined) ---
    if(length(top_terms_id) > 0){
      tryCatch({
        p4 <- scplotter::GSEAPlot(
          data = gsea_sym, 
          gene_ranks = gsea_sym@geneList, 
          gene_sets = gsea_sym@geneSets, 
          gs = top_terms_id[1:4], 
          in_form = "dose", combine = TRUE, ncol = 2, 
        )
        ggsave1(file.path(out_dir, "4_GSEAPlot_Combined"), 12, 16, plot = p4)
      }, error = function(e) message(paste("  [scplotter Error]", ct, ":", e$message)))
    }
  }
}

# === 2. 执行绘图 ===
rds_dir_path <- file.path(save_dir, "RDS")
plot_base_dir <- file.path(save_dir, "GSEA_Plots")

# 分别运行 GO 和 KEGG
plot_gsea_pipeline("GO", rds_dir_path, plot_base_dir)
plot_gsea_pipeline("KEGG", rds_dir_path, plot_base_dir)


# 拼图
library(magick)
# 拼接差异基因
de_path <- file.path(output_dir, "DE_Analysis", paste0(case, "_vs_", control, "_DEG"))
img1 <- image_read(file.path(de_path,"all_cluster_volcano_pct.png")) %>% image_scale("5000")
img2 <- image_read(file.path(de_path,paste0("Volcano_",celltypes[1],".png"))) %>% image_scale("x2000")
img3 <- image_read(file.path(de_path,"DEG_top5_Heatmap.png")) %>% image_scale("x2000")
img4 <- image_read(file.path(de_path,"Number_of_DEGs.png")) %>% image_scale("x2000")
# 进行拼接布局
row2 <- image_append(c(img2, img3, img4)) %>% image_scale("5000")
combined_grid1 <- image_append(c(img1, row2), stack = TRUE)
image_write(combined_grid1, path = file.path(output_dir,"DE_Analysis", paste0(case, "_vs_", control, "_combined_plot.png")), format = "png")

# 拼接Enrichment
enrich_path <- file.path(output_dir, "Enrichment", paste0(case, "_vs_", control, "_Enrichment"),"Enrichment_Plots")
img2 <- image_read(file.path(enrich_path,"Network_Plot.png")) %>% image_scale("x2000")
img1 <- image_read(file.path(enrich_path,paste0(celltypes[1],"_KEGG_Dot.png"))) %>% image_scale("x2000")
img3 <- image_read(file.path(enrich_path,"Comparison.png")) %>% image_scale("x2000")
img4 <- image_read(file.path(enrich_path,"Complex_combined_bar.png")) %>% image_scale("x2000")
# 进行拼接布局
row1 <- image_append(c(img1, img2)) %>% image_scale("5000")
row2 <- image_append(c(img3, img4)) %>% image_scale("5000")
combined_grid2 <- image_append(c(row1, row2), stack = TRUE)
image_write(combined_grid2, path = file.path(output_dir,"Enrichment", paste0(case, "_vs_", control, "_combined_plot_Enrichment.png")), format = "png")

# 拼接GSEA
gsea_path <- file.path(output_dir,"Enrichment",paste0(case, "_vs_", control, "_GSEA"),"GSEA_Plots","KEGG" ,paste0(case, "_vs_", control, "_",celltypes[1]))
img1 <- image_read(file.path(gsea_path,"1_Barplot.png")) %>% image_scale("x2000")
img2 <- image_read(file.path(gsea_path,"2_RidgePlot.png")) %>% image_scale("x2000")
img3 <- image_read(file.path(gsea_path,"3_GseaVis.png")) %>% image_scale("x2000")
img4 <- image_read(file.path(gsea_path,"4_GSEAPlot_Combined.png")) %>% image_scale("x2000")
# 进行拼接布局
row1 <- image_append(c(img1, img2)) %>% image_scale("5000")
row2 <- image_append(c(img3, img4)) %>% image_scale("5000")
combined_grid2 <- image_append(c(row1, row2), stack = TRUE)
image_write(combined_grid2, path = file.path(output_dir,"Enrichment", paste0(case, "_vs_", control, "_combined_plot_GSEA.png")), format = "png")