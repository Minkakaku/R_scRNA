#######################################################################################
#######################################################################################
### 【注释1.5】最大得分大于次大得分两倍以上，按照预设的注释将细胞分群，输入obj_patient，输出obj_patient_annotated，并保存
#######################################################################################
#######################################################################################

# 定义特征基因列表  
# 定义特征基因列表
feature_genes <- list(
  Tumor_cell = c("EPCAM", "MUC1", "CLDN1", "CDH1", "TP63",  #general
                 "KRAS", "SDC4", "NKX2-1","NTRK2","MUC16" #sup
                 ),
  CAF = c("COL5A1", "COL5A2", "COL11A1", "POSTN", "MMP11", "MMP13", "FAP",
          "TNXB","PRG4"
          ),
  T_cell = c("CD3E", "CD2", "CD247", "CD4", "CD8A", "CD5", "CD6", "CD7", "THEMIS", "ITK", #pan
             "NCAM1", "PRF1", "GZMA", "GZMB", "GZMK", "GZMH", "PRF1", #TCL
             "HAVCR2", "PDCD1", "TIGIT", "LAG3", "CD96", #Tex
             "FOXP3", "IL2RA", "CCR8", "CTLA4", #Treg
             "CXCR3", "IFNG", "CCR4", "IL10", #Th
             "CXCR5", "CXCL13", #Tfh
             "TCF7","CD28" #Tn
             ),
  Endothelial_cell = c("PECAM1", "EPAS1","CD34","NOTCH4", "ENG", "PLVAP","S1PR1","ESM1", #Pan
                       "GJA5", #art
                       "SELE","SELP", #vein
                       "CA4", #cap
                       "PDPN","PROX1" #lymphatic
                       ),
  Mural_cell = c("COL4A1", "COL4A2", "PDGFRB", "NOTCH3","RGS5","MCAM","KCNJ8"),
  B_cell = c("CD19", "CD79A", "MS4A1", "CD79B", "BLK","PAX5","BNAK1","CD22"),
  Plasma_cell = c("CD38", "MZB1", "POU2AF1","XBP1","PIM2","TNFRSF17","FCRL5","IRF4","SLAMF7"),
  Mast_cell = c("HDC", "KIT", "IL1RL1","SLC18A2","MS4A2"),
  Macrophage = c("FCGR2A", "CD14", "CYBB", "PTAFR", "ALOX5", "C5AR1", "CLEC7A", "TLR2", "CD86","TREM2", #M1  
                 "MRC1", "CD163", "STAB1", "SIGLEC1", "LAIR1", "CCL18", "VSIG4", "TREM2","CXCL9","FOLR2", #M2
                 "MSR1", "LGMN", "LIPA", "CSF1R", "CPVL", "MERTK", "MARCO", "C3AR1", #TRM
                 "CTSL", "CD68", "ITGAX", "CLEC5A", "ITGAM", "CD86", "CXCL16" #general
                 ),
  Neutrophil = c("CXCR1", "CSF3R","CXCR2","CEACAM3","PADI4","DEFA1","FPR2"),
  DC = c("CLEC9A","CXR1","CLNK","ENPP1", #DC1
         "CD1C","FCER1A","CLEC10A", #DC2
         "LAMP3","CCR7","FSCN1","CCL22" #DC3
         )
)

# 定义注释函数  
annotate_cells <- function(seurat_obj) {  
  # 获取标准化数据  
  data <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")  
  
  # 检查基因存在性  
  missing_genes <- sapply(feature_genes, function(genes) {  
    setdiff(genes, rownames(data))  
  })  
  if (any(sapply(missing_genes, length) > 0)) {  
    warning("Missing genes: ", paste(unique(unlist(missing_genes)), collapse = ", "))  
  }  
  
  # 计算特征得分矩阵  
  score_matrix <- sapply(names(feature_genes), function(cell_type) {  
    genes <- intersect(feature_genes[[cell_type]], rownames(data))  
    if (length(genes) == 0) return(rep(0, ncol(data)))  
    Matrix::colMeans(data[genes, , drop = FALSE])  
  })  
  
  # 确定最大得分和次大得分  
  max_score <- apply(score_matrix, 1, max)  
  second_max_score <- apply(score_matrix, 1, function(x) {  
    sort(x, decreasing = TRUE)[2]  
  })  
  
  # 判断是否满足分类条件  
  classification_condition <- max_score > 1.5 * second_max_score  
  
  # 保留满足条件的细胞  
  seurat_obj <- subset(seurat_obj, cells = colnames(data)[classification_condition])  
  
  # 获取满足条件细胞的最终注释  
  final_annot <- apply(score_matrix[classification_condition, ], 1, function(x) {  
    names(which.max(x))  
  })  
  
  # 设置细胞标识  
  seurat_obj$celltype <- factor(final_annot, levels = names(feature_genes))  
  Idents(seurat_obj) <- "celltype"  
  
  # 强制保留所有类别  
  levels(Idents(seurat_obj)) <- names(feature_genes)  
  
  return(seurat_obj)  
}  

# 批量处理所有样本  
obj_patient_annotated <- lapply(obj_patient, function(obj) {  
  tryCatch({  
    annotate_cells(obj)  
  }, error = function(e) {  
    message("Error processing ", obj@project.name, ": ", e$message)  
    return(NULL)  
  })  
})  

# 过滤处理失败的样本  
obj_patient_annotated <- obj_patient_annotated[!sapply(obj_patient_annotated, is.null)]  

# 验证结果  
sapply(obj_patient_annotated, function(obj) {  
  c(  
    Cells = ncol(obj),  
    Celltypes = length(unique(Idents(obj)))  
  )  
})  

# 保存最终结果  
saveRDS(obj_patient_annotated, file = "obj_patient_annotated.rds")