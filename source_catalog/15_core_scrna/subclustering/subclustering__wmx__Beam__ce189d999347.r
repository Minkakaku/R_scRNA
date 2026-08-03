rm(list = ls())
library('igraph', lib.loc = '/mnt/NexData_file/home/wmx/envs/scRNA/lib/R/library/') # 加载2.0.3的igraph不然ordercells会报错
.libPaths("/mnt/NexData_file/share/scRNA-tools/envs/pip_monocle2/lib/R/library/") 
library(monocle)
library(ggthemes,lib.loc = '/mnt/NexData_file/home/wmx/envs/scRNA/lib/R/library/')
library(Seurat)
library(BiocGenerics)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(reshape2)
library(cowplot)
library(pheatmap)
library(scales)
library(plyr)
library(VGAM)


getgene = function(seuCDS,cds_subset,genelist,qvalueCutOff,TopGeneNum){
  if (is.null(genelist)) {
    # 使用Beam直接做速度巨慢
    genelist = seuCDS@featureData@data$gene_short_name[seuCDS@featureData@data$use_for_ordering]
    BEAM_res = BEAM(seuCDS[genelist,],progenitor_method = "duplicate",
                branch_point = branchPoint,cores = 1)
    # BEAM_res = BEAM(seuCDS,branch_point = branchPoint,cores = 32)
    row.names(BEAM_res) = BEAM_res$gene_short_name
    ingene = intersect(row.names(BEAM_res),row.names(seuCDS))
    BEAM_res = BEAM_res[ingene,] %>% arrange(qval)
    # qvalueCutOff  =2
    subBEAM_res = subset(BEAM_res,qval<qvalueCutOff)
    subBEAM_res = subBEAM_res[subBEAM_res$use_for_ordering =='TRUE',]
    write.table(as.data.frame(BEAM_res),
                "AllBEAMGene.txt",
                sep="\t",row.names=F,quote=F)
    # BEAM_res = read.table("AllBEAMGene.txt",sep = "\t",header = T)
    # TopGeneNum=500
    genes = subBEAM_res[1:TopGeneNum,"gene_short_name"]
    genes = na.omit(genes)
    # genes = subBEAM_res$gene_short_name
  }else{
    genes = unique(as.character(read.table(genelist,sep = "\t",header = T)[,1]))
    genes = intersect(rownames(cds_subset),genes)
  }
  print("getgene done")
  return(genes)
}

beamheatmap = function(seuCDS,branchPoint,genes,numClusters){
  new_cds <- buildBranchCellDataSet(seuCDS[genes,], 
                                    branch_states=NULL, 
                                    branch_point=branchPoint, 
                                    progenitor_method = 'duplicate',
                                    )
  new_cds@dispFitInfo <- seuCDS@dispFitInfo
  branch_states = NULL
  if(is.null(branch_states)) {
    progenitor_state <- subset(pData(seuCDS), Pseudotime == 0)[, 'State']
    branch_states <- setdiff(pData(seuCDS)$State, progenitor_state)
  }
  col_gap_ind <- 101
  numClusters = 3
  listmap = plot_genes_branched_heatmap(seuCDS[genes[1:20],],
                                        branch_point = branchPoint,
                                        num_clusters = numClusters,
                                        cores = 1,
                                        use_gene_short_name = T,
                                        show_rownames = T,
                                        return_heatmap=TRUE)
  
  ggsave(paste0("Branch",branchPoint,".genes_branched_heatmap.pdf"), listmap$ph_res)
  ggsave(paste0("Branch",branchPoint,".genes_branched_heatmap.png"), listmap$ph_res)
  gene2Cluster = data.frame(Gene=rownames(listmap$annotation_row),listmap$annotation_row)
  order = listmap$ph$tree_row$order
  gene2Cluster = gene2Cluster[order,]
  write.table(gene2Cluster,file=paste0("Branch.",branchPoint,".Gene2Cluster.txt"),sep="\t",quote=F,row.names=F)
  getwd()
  print("beamheatmap done")
}

plot_genes_branched_pseudotime_new = function (cds_subset, 
                                                branch_states = NULL, 
                                                branch_point=1,
                                                branch_labels = NULL,
                                                method = "fitting", 
                                                min_expr = NULL, 
                                                cell_size = 0.75,
                                                nrow = NULL, 
                                                ncol = 1, 
                                                panel_order = NULL, 
                                                color_by = "State",
                                                expression_curve_linetype_by = "Branch", 
                                                trend_formula = "~ sm.ns(Pseudotime, df=3) * Branch", 
                                                reducedModelFormulaStr = NULL, 
                                                label_by_short_name = TRUE,
                                                relative_expr = TRUE,
                                                logy = TRUE,
                                                #gene_pairs = NULL,
                                                ...){
  Branch <- NA  
  if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
    integer_expression <- TRUE
  }
  else {
    integer_expression <- FALSE
  }
  if (integer_expression) {
    CM <- exprs(cds_subset)
    if (relative_expr){
      if (is.null(sizeFactors(cds_subset))) {
        stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
      }
      CM <- Matrix::t(Matrix::t(CM)/sizeFactors(cds_subset))
    }
    cds_exprs <- reshape2::melt(round(as.matrix(CM)))
  }
  else {
    cds_exprs <- reshape2::melt(exprs(cds_subset))
  }
  if (is.null(min_expr)) {
    min_expr <- cds_subset@lowerDetectionLimit
  }
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_pData <- pData(cds_subset)
  
  cds_fData <- fData(cds_subset)
  cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
  cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
  if (integer_expression) {
    cds_exprs$adjusted_expression <- round(cds_exprs$expression)
  }
  else {
    cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
  }
  if (label_by_short_name == TRUE) {
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
      cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    }
    else {
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  }
  else {
    cds_exprs$feature_label <- cds_exprs$f_id
  }
  cds_exprs$feature_label <- as.factor(cds_exprs$feature_label)
  # trend_formula <- paste("adjusted_expression", trend_formula,
  #     sep = "")
  cds_exprs$Branch <- as.factor(cds_exprs$Branch) 
  
  new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime, Branch = pData(cds_subset)$Branch)
  
  full_model_expectation <- genSmoothCurves(cds_subset, cores=1, trend_formula = trend_formula, 
                                            relative_expr = T, new_data = new_data)
  colnames(full_model_expectation) <- colnames(cds_subset)
  
  cds_exprs$full_model_expectation <- apply(cds_exprs,1, function(x) full_model_expectation[x[2], x[1]])
  if(!is.null(reducedModelFormulaStr)){
    reduced_model_expectation <- genSmoothCurves(cds_subset, cores=1, trend_formula = reducedModelFormulaStr,
                                                 relative_expr = T, new_data = new_data)
    colnames(reduced_model_expectation) <- colnames(cds_subset)
    cds_exprs$reduced_model_expectation <- apply(cds_exprs,1, function(x) reduced_model_expectation[x[2], x[1]])
  }
  
  # FIXME: If you want to show the bifurcation time for each gene, this function
  # should just compute it. Passing it in as a dataframe is just too complicated
  # and will be hard on the user. 
  # if(!is.null(bifurcation_time)){
  #     cds_exprs$bifurcation_time <- bifurcation_time[as.vector(cds_exprs$gene_short_name)]
  # }
  if (method == "loess")
    cds_exprs$expression <- cds_exprs$expression + cds@lowerDetectionLimit
  if (label_by_short_name == TRUE) {
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
      cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    }
    else {
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  }
  else {
    cds_exprs$feature_label <- cds_exprs$f_id
  }
  cds_exprs$feature_label <- factor(cds_exprs$feature_label)
  if (is.null(panel_order) == FALSE) {
    cds_exprs$feature_label <- factor(cds_exprs$feature_label,
                                      levels = panel_order)
  }
  cds_exprs$expression[is.na(cds_exprs$expression)] <- min_expr
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_exprs$full_model_expectation[is.na(cds_exprs$full_model_expectation)] <- min_expr
  cds_exprs$full_model_expectation[cds_exprs$full_model_expectation < min_expr] <- min_expr
  
  if(!is.null(reducedModelFormulaStr)){
    cds_exprs$reduced_model_expectation[is.na(cds_exprs$reduced_model_expectation)] <- min_expr
    cds_exprs$reduced_model_expectation[cds_exprs$reduced_model_expectation < min_expr] <- min_expr
  }
  
  cds_exprs$State <- as.factor(cds_exprs$State)
  cds_exprs$Branch <- as.factor(cds_exprs$Branch)
  
  q <- ggplot(aes(Pseudotime, expression), data = cds_exprs)

  if (is.null(color_by) == FALSE) {
    q <- q + geom_point(aes_string(color = color_by), size = I(cell_size))
  }
  if(logy){
    q <- q + scale_y_log10()
  }
  
  if (is.null(reducedModelFormulaStr) == FALSE)
    q <- q +  facet_wrap(~feature_label +
                           pval, nrow = nrow, ncol = ncol, scales = "free_y")
  else q <- q +  facet_wrap(~feature_label,
                            nrow = nrow, ncol = ncol, scales = "free_y")
  if (method == "loess")
    q <- q + stat_smooth(aes(fill = Branch, color = Branch),
                         method = "loess")
  else if (method == "fitting") {
    q <- q + geom_line(aes_string(x = "Pseudotime", y = "full_model_expectation",
                                  linetype = "Branch"), data = cds_exprs) #+ scale_color_manual(name = "Type", values = c(colour_cell, colour), labels = c("Pre-branch", "AT1", "AT2", "AT1", "AT2")
  }
  
  if(!is.null(reducedModelFormulaStr)) {
    q <- q + geom_line(aes_string(x = "Pseudotime", y = "reduced_model_expectation"),
                       color = 'black', linetype = 2, data =  cds_exprs)   
  }
  
  q <- q + ylab("Expression") + xlab("Pseudotime (stretched)")
  
  q <- q + monocle_theme_opts()
  q + expand_limits(y = min_expr)
}
monocle_theme_opts <- function(){
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}
beamgene = function(cds_subset,genes,branchPoint){
  dir.create("GeneBranched")
  for (i in genes) {
    gg = plot_genes_branched_pseudotime_new(cds_subset[i,],
                                            branch_point = branch_point,
                                            color_by = "Cluster",
                                            ncol = 1,
                                            logy=TRUE,
                                            cell_size=1)
    ggsave(paste0("GeneBranched/",i,".branched.png"),gg,width=6,height = 3)
    ggsave(paste0("GeneBranched/",i,".branched.pdf"),gg,width=6,height = 3)
  }
  zip("GeneBranched.zip","./GeneBranched")
  
  file.remove(paste0("./GeneBranched/",dir("./GeneBranched")))
  file.remove("./GeneBranched/")
}
mainfun = function(rds,branchPoint,genelist,numClusters,qvalueCutOff,TopGeneNum){
  # rds
  seuCDS = readRDS(rds)
  print("Rds load finish")
  # seuCDS <- estimateDispersions(seuCDS)
  # print("despersions finish")
  # 构建branch
  cds_subset = buildBranchCellDataSet(seuCDS,branch_point = branchPoint,progenitor_method ="duplicate")
  # dir.create("BEAM")
  # setwd("./BEAM")
  print("BranchCellData finish")
  genes = getgene(seuCDS,cds_subset,genelist,qvalueCutOff,TopGeneNum)
  # genes  = read.csv('/lvdata/wzb/scRNA/FW2024-213/2024.05.20/NK/BEAM1/AllBEAMGene.txt',sep = "\t")
  # genes = genes$gene_short_name[1:100]
  print(genes)
  print("genes Cal finish")
  cds_subset = cds_subset[genes,]
  beamheatmap(seuCDS,branchPoint,genes,numClusters)  
  print("Beam heatmap finish")
  # 这一部分先不画图
  # if (length(genes)<=1000) {
  #   beamgene(cds_subset,genes,branchPoint)
  # }
  
}


rds = "/mnt/NexData_file/home/wmx/Project/AE2508064183/shouhou5/monocle_real_count/Pseudotime/Pseudotime.monocle.rds"
branchPoint = 1
qvalueCutOff = 0.05
TopGeneNum = 100
genelist = NULL
numClusters = 3
outpath = '/mnt/NexData_file/home/wmx/Project/AE2508064183/shouhou5/monocle_real_count/BEAM_1'
dir.create(outpath)
setwd(outpath)
mainfun(rds,branchPoint,genelist,numClusters,qvalueCutOff,TopGeneNum)

branchPoint = 2 
outpath = '/mnt/NexData_file/home/wmx/Project/AE2508064183/shouhou5/monocle_real_count/BEAM_2'
dir.create(outpath)
setwd(outpath)
mainfun(rds,branchPoint,genelist,numClusters,qvalueCutOff,TopGeneNum)

# 补充 热图基因的表达分面图
dir_paths <- c(
  "/mnt/NexData_file/home/wmx/Project/AE2508064183/shouhou5/monocle_real_count/BEAM_1/",
  "/mnt/NexData_file/home/wmx/Project/AE2508064183/shouhou5/monocle_real_count/BEAM_2/"
)

ttf <- unlist(lapply(dir_paths, function(d) {
  list.files(d, pattern = "\\.Gene2Cluster\\.txt$", full.names = TRUE)
}))
print(ttf)
cds <- readRDS(rds)

for (f in ttf) {
    # 读取文件
    gene_data <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    selectgenes <- as.character(gene_data[, 1])
    
    # 根据文件路径自动判断 branch_point
    # 如果路径包含 BEAM_1 则为 1，包含 BEAM_2 则为 2
    bp <- ifelse(grepl("BEAM_1", f), 1, 2)
    
    message("Drawing: ", basename(f), " | Branch Point: ", bp)
    
    # 限制基因数量并过滤不存在于 cds 中的基因
    genes_to_plot <- selectgenes[1:min(100, length(selectgenes))]
    genes_to_plot <- intersect(genes_to_plot, rownames(cds))
    
    if(length(genes_to_plot) > 0) {
        # 绘图
        p = plot_genes_branched_pseudotime(cds[genes_to_plot,],
                                           branch_point = bp,
                                           color_by = "group",
                                           ncol = 4)
        ggsave(p,
               filename = paste0("/mnt/NexData_file/home/wmx/Project/AE2508064183/shouhou5/monocle_real_count/BEAM_",bp,"/Branch_gene_plot.pdf"),
                width = 20,height = 25)
        ggsave(p,
               filename = paste0("/mnt/NexData_file/home/wmx/Project/AE2508064183/shouhou5/monocle_real_count/BEAM_",bp,"/Branch_gene_plot.png"),
                width = 20,height = 25)         
    } else {
        message("No matching genes found in cds for file: ", basename(f))
    }
}

# 补充 拟时序相关基因的热图
fData(cds)$num_cells_expressed <- rowSums(exprs(cds) > 0)
ordergene = cds@featureData@data$gene_short_name[cds@featureData@data$use_for_ordering]
diff_test_res <- differentialGeneTest(cds[ordergene,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- diff_test_res %>% dplyr::filter(num_cells_expressed >10)%>% #防止小类干扰
  dplyr::arrange(qval) %>% head(100) %>% row.names() 

write.csv(file = paste0('/mnt/NexData_file/home/wmx/Project/AE2508064183/shouhou5/monocle_real_count/','Pseudotime_ordergene.csv'),diff_test_res)
p1 = plot_pseudotime_heatmap(cds[sig_gene_names,],
                            num_clusters = 5,
                            cores = 1,
                            show_rownames = T,return_heatmap = TRUE)
ggsave(filename = paste0('/mnt/NexData_file/home/wmx/Project/AE2508064183/shouhou5/monocle_real_count/',"monocle.top100gene_cortime.png"),p1,h=8,w=6,dpi = 300)
ggsave(filename = paste0('/mnt/NexData_file/home/wmx/Project/AE2508064183/shouhou5/monocle_real_count/',"monocle.top100gene_cortime.pdf"),p1,h=8,w=6,dpi = 300)


# test
# 1) 只抽少量基因测试，别动细胞
# set.seed(123)

# genes_pool <- rownames(subset(fData(seuCDS), use_for_ordering == TRUE))
# if (length(genes_pool) == 0) genes_pool <- rownames(seuCDS)
# genes_use <- sample(genes_pool, min(50, length(genes_pool)))

# cds_test <- seuCDS[genes_use, ]
# # cds_test <- estimateDispersions(cds_test)

# BEAM_res_test <- BEAM(
#   cds_test,
#   branch_point = 1,
#   progenitor_method = "duplicate",
#   cores = 1
# )

# BEAM_res_test <- BEAM_res_test[order(BEAM_res_test$qval), ]
# head(BEAM_res_test[, c("gene_short_name", "pval", "qval", "status")], 20)
# table(BEAM_res_test$status, useNA = "ifany")
