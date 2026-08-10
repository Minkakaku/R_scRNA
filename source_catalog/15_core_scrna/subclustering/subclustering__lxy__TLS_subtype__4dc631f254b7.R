all_data <- readRDS('R_OUTPUT/all_data.rds')

all_data$annotation %>% table()
all_data$sample %>% table()

TLS_out_dir <- "TLS_subtype" 
dir.create(TLS_out_dir)

# 1. 读取你的CSV文件（细胞ID列表）
#sed -i '1,2d' *.csv


# 定义CSV文件目录
csv_dir <- "/mnt/NexData_file/home/lxy/case/XE2603162231-02/20260616_XE2603162231-02/TLS_subtype/ROI"

# 获取所有CSV文件
csv_files <- list.files(csv_dir, pattern = "*_cells_stats\\.csv$", full.names = TRUE)

# 初始化sample_TLS列
all_data$sample_TLS <- NA

# 批量添加标签
for(csv_file in csv_files) { #browser()
  base_name <- str_remove(basename(csv_file), "_cells_stats\\.csv$")
  cell_data <- read.csv(csv_file, stringsAsFactors = FALSE)
  csv_cells <- cell_data$Cell.ID
  
  csv_cells_prefixed <- paste0("TMA-1_", csv_cells)
  
  # 匹配并添加标签
  cells_to_label <- intersect(rownames(all_data@meta.data), csv_cells_prefixed)
  if(length(cells_to_label) > 0) {
    all_data@meta.data[cells_to_label, "sample_TLS"] <- base_name
    message(sprintf("标记 %d 个细胞为 %s", length(cells_to_label), base_name))
  }
}

table(all_data$sample_TLS)

# 提取有标签的细胞
cells_with_label <- rownames(all_data@meta.data)[!is.na(all_data$sample_TLS)]
all_data_TLS <- subset(all_data, cells = cells_with_label)

# 保存

saveRDS(all_data_TLS, file.path(TLS_out_dir, "all_data_sampleTLS.rds"))


#TLS <- readRDS('TLS_subtype/all_data_sampleTLS.rds')
TLS <- all_data_TLS


## PCA --------
TLS <- FindVariableFeatures(TLS)

TLS <- ScaleData(TLS, vars.to.regress = 'nCount_RNA')
TLS <- RunPCA(TLS)

### pca inspector ------
pdf(file.path(TLS_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(TLS, ndims = 50)
dev.off()


## dimension reduction --------
dims <- 22 # 18

#library(future)
#options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
#plan(multisession, workers = 10)

TLS <- RunUMAP(TLS, dims = 1:dims,return.model = TRUE)

pdf(file.path(TLS_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(TLS, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

TLS <- TLS %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

TLS <- RunUMAP(
  TLS, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(TLS_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(TLS, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
TLS  <-
  FindNeighbors(
    TLS ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


TLS <-
  FindClusters(
    TLS,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(TLS_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(TLS, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(TLS, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(TLS, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(TLS, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(TLS, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(TLS, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(TLS, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

Idents(TLS) <- "RNA_snn_res.0.6"

### dotplot -------
TLS_major_marker_ls <- c(
  # Tfh 滤泡辅助T细胞
  "Tfh" = "Cd4",
  "Tfh" = "Cxcr5",
  "Tfh" = "Icos",
  "Tfh" = "Pdcd1",
  #"Tfh" = "Bcl6",
  "Tfh" = "Il21",
  "Tfh" = "Maf",
  "Tfh" = "Sh2d1a",
  "Tfh" = "Batf",
  "Tfh" = "Tox2",
  "Tfh" = "Ascl2",
  "Tfh" = "Slamf6",
  "Tfh" = "Cd40lg",
  "Tfh" = "Cxcr4",
  
  # Pan-B 泛B细胞
  "Pan-B" = "Cd79a",
  "Pan-B" = "Cd79b",
  "Pan-B" = "Cd19",
  "Pan-B" = "Ms4a1",
  "Pan-B" = "Cd74",
  "Pan-B" = "H2-Aa",
  "Pan-B" = "H2-Ab1",
  "Pan-B" = "Cd37",
  "Pan-B" = "Cd22",
  "Pan-B" = "Bank1",
  "Pan-B" = "Pax5",
  
  # GC B 生发中心B细胞
  "GC B" = "Aicda",
  "GC B" = "Bcl6",
  "GC B" = "Fas",
  "GC B" = "S1pr2",
  "GC B" = "Rgs13",
  "GC B" = "Gcsam",
  #"GC B" = "Cxcr4", # 可选标记
  
  # Cycling B 增殖B细胞
  "Cycling B" = "Mki67",
  "Cycling B" = "Top2a",
  "Cycling B" = "Ccnb1",
  "Cycling B" = "Ccna2",
  "Cycling B" = "Mcm5",
  "Cycling B" = "Stmn1",
  
  # Plasmablast/Plasma 浆母/浆细胞
  "Plasmablast/Plasma" = "Jchain",
  "Plasmablast/Plasma" = "Sdc1",
  "Plasmablast/Plasma" = "Xbp1",
  "Plasmablast/Plasma" = "Prdm1",
  "Plasmablast/Plasma" = "Irf4",
  "Plasmablast/Plasma" = "Mzb1",
  #"Plasmablast/Plasma" = "Tnfrsf17", # 可选标记
  
  # FRC-like 成纤维网状细胞
  "FRC-like" = "Pdpn",
  "FRC-like" = "Pdgfra",
  "FRC-like" = "Ccl19",
  "FRC-like" = "Ccl21a",
  "FRC-like" = "Il7",
  "FRC-like" = "Ltbr",
  "FRC-like" = "Icam1",
  "FRC-like" = "Vcam1",
  "FRC-like" = "Thy1",
  "FRC-like" = "Col1a1",
  "FRC-like" = "Col1a2",
  "FRC-like" = "Col3a1",
  "FRC-like" = "Col6a1",
  "FRC-like" = "Dcn",
  "FRC-like" = "Lum",
  "FRC-like" = "Fn1",
  
  # HEV-like 高内皮微静脉内皮细胞
  "HEV-like" = "Pecam1",
  "HEV-like" = "Cdh5",
  "HEV-like" = "Kdr",
  "HEV-like" = "Vwf",
  "HEV-like" = "Emcn",
  "HEV-like" = "Plvap",
  "HEV-like" = "Glycam1",
  "HEV-like" = "Cd34",
  "HEV-like" = "Podxl",
  "HEV-like" = "Madcam1",
  "HEV-like" = "Chst4",
  "HEV-like" = "Chst2",
  "HEV-like" = "Fut7",
  "HEV-like" = "Fut4",
  "HEV-like" = "Gcnt1",
  "HEV-like" = "B3gnt3",
  "HEV-like" = "St3gal4",
  "HEV-like" = "St3gal6"
)


B_major_marker_ls <- c(
  'B cell' = 'Cd19', 
  'B cell' = 'Cd79a', 
  'B cell' = 'Ms4a1', 
  
  'naive B' = 'Cd20', 
  'naive B' = 'Cd24', 
  'naive B' = 'Ighm', 
  'naive B' = 'Bank1', 
  'naive B' = 'Fcrl1', 
  'naive B' = 'Hla-Dqb1', 
  'naive B' = 'Hla-Dqa1', 
  'naive B' = 'Hvcn1', 
  'naive B' = 'Ighd', 
  'naive B' = 'Fcer2', 
  'naive B' = 'Il4r', 
  
  'Follicular B' = 'Cd40', 
  'Follicular B' = 'Bcl6', 
  
  'memory B' = 'Cd27', 
  'memory B' = 'Fcrl4', 
  'memory B' = 'Taci', 
  'memory B' = 'Aim2', 
  'memory B' = 'Txnip', 
  'memory B' = 'Gpr183', 
  'memory B' = 'Tnfrsf13b', 
  'memory B' = 'Cd22', 
  'memory B' = 'Fcrl2', 
  'memory B' = 'Adam28', 
  'memory B' = 'Igll5', 
  'memory B' = 'P2rx5', 
  'memory B' = 'Blk', 
  
  #'Marginal zone B' = 'Cr2', 
  'Marginal zone B' = 'Cd1c', 
  
  'cycling B' = 'Mki67', 
  'cycling B' = 'Stmn1', 
  'cycling B' = 'Hmgb2', 
  'cycling B' = 'Birc5', 
  
  'Plasma' = 'Cd138', 
  'Plasma' = 'Prdm1', 
  'Plasma' = 'Blimp1', 
  'Plasma' = 'Igha1', 
  'Plasma' = 'Igha2', 
  'Plasma' = 'Ighg1', 
  'Plasma' = 'Ighg2', 
  'Plasma' = 'Ighg3', 
  'Plasma' = 'Ighg4', 
  'Plasma' = 'Igkc', 
  'Plasma' = 'Tnfrsf17', 
  'Plasma' = 'Jchain', 
  'Plasma' = 'Mzb1', 
  'Plasma' = 'Xbp1', 
  'Plasma' = 'Sdc1', 
  
  'Plasmablast' = 'Cd38', 
  'Plasmablast' = 'Irf4', 
  
  'Breg' = 'Il10', 
  'Breg' = 'Cd1d', 
  'Breg' = 'Tgfb1', 
  
  'GC B' = 'S1pr2', 
  'GC B' = 'Rgs13', 
  'GC B' = 'Tcl1a', 
  'GC B' = 'Lrmp', 
  'GC B' = 'Sugct', 
  'GC B' = 'Mme', 
  'GC B' = 'Aicda',
  
  "FDC" = "Cr2",
  "FDC" = "Cr1"
  #"FDC" = "FCER2",
)

T_major_marker_ls <- c(
  'CD8' = 'Cd8a', 
  'CD8' = 'Cd8b', 
  
  'CD4' = 'Cd4', 
  'CD4' = 'Trdc', 
  'CD4' = 'Il2', 
  
  'Naive' = 'Lef1', 
  'Naive' = 'Sell', 
  'Naive' = 'Ybx3', 
  'Naive' = 'Tgfb1', 
  'Naive' = 'Cxcr3', 
  
  'Memory' = 'Gpr183', 
  'Memory' = 'Il7r', 
  'Memory' = 'Fos', 
  'Memory' = 'Jun', 
  'Memory' = 'Klrb1', 
  
  'TRM' = 'Cd89', 
  'TRM' = 'Cxcr8', 
  
  'TEM' = 'Ocl4', 
  'TEM' = 'Ocl5', 
  'TEM' = 'Jak3', 
  'TEM' = 'Gzmk', 
  'TEM' = 'Gzma', 
  'TEM' = 'Nkg7', 
  
  'TEMRA' = 'Klrg1', 
  
  'TCM' = 'Tcf7', 
  'TCM' = 'Ido1', 
  'TCM' = 'Skiv2l', 
  'TCM' = 'Traf2', 
  'TCM' = 'Slc2a1', 
  'TCM' = 'Cd40lg', 
  'TCM' = 'Dusp2', 
  'TCM' = 'Lyar', 
  'TCM' = 'Il32', 
  'TCM' = 'Ccr7', 
  
  'TRM' = 'Znf683', 
  'TRM' = 'Cxcr6', 
  
  'Effector' = 'Fgfbp2', 
  'Effector' = 'Tyrobp', 
  'Effector' = 'Klrc2', 
  'Effector' = 'Cx3cr1', 
  'Effector' = 'Fcgr3a', 
  'Effector' = 'Gnly', 
  
  'cytotoxic' = 'Xcl1', 
  'cytotoxic' = 'Ifng', 
  'cytotoxic' = 'Ccl4', 
  'cytotoxic' = 'Ccl5', 
  
  'Exhausted' = 'Layn', 
  'Exhausted' = 'Lag3', 
  'Exhausted' = 'Havcr2', 
  'Exhausted' = 'Tnfrsf9', 
  'Exhausted' = 'Ctla4', 
  'Exhausted' = 'Tigit', 
  'Exhausted' = 'Tox', 
  
  'Tfh' = 'Cd200', 
  'Tfh' = 'Cd83', 
  'Tfh' = 'Cxcl3', 
  'Tfh' = 'Icos', 
  'Tfh' = 'Batf', 
  'Tfh' = 'Pdcd1', 
  'Tfh' = 'Cxcr5', 
  
  'proliferation' = 'Mki67', 
  'proliferation' = 'H4c3', 
  'proliferation' = 'Tubb', 
  'proliferation' = 'Stmn1', 
  'proliferation' = 'Top1', 
  'proliferation' = 'Anapc11', 
  'proliferation' = 'Mcm7', 
  'proliferation' = 'Dnmt1', 
  'proliferation' = 'Top2a', 
  
  'CD4_Treg' = 'Tnfrsf4', 
  'CD4_Treg' = 'Foxp3', 
  'CD4_Treg' = 'Il2ra', 
  'CD4_Treg' = 'Cd25', 
  'CD4_Treg' = 'Stat1', 
  
  'CD4_Tem' = 'Capg', 
  
  'CD4_Th1' = 'Anxa1', 
  'CD4_Th1' = 'Cxcl13', 
  'CD4_Th1' = 'Bhlhe40', 
  
  'CD4_Tfh' = 'Tox2', 
  
  'CD4_Th17' = 'Il17a', 
  'CD4_Th17' = 'Il17f', 
  'CD4_Th17' = 'Il23r', 
  'CD4_Th17' = 'Rorc', 
  'CD4_Th17' = 'Ccr6', 
  
  'CD4_Tcm' = 'S1pr1', 
  'CD4_Tcm' = 'Txk', 
  'CD4_Tcm' = 'Icam2', 
  'CD4_Tcm' = 'Cdc25b', 
  'CD4_Tcm' = 'Ptger2', 
  
  'CD8_Tex' = 'Enpd1', 
  
  'CD8_Tem' = 'Gzmb', 
  'CD8_Tem' = 'Prf1', 
  
  'CD8_MAIT' = 'Ncr3', 
  'CD8_MAIT' = 'Slc4a10', 
  
  'CD8_Trm' = 'Myadm', 
  'CD8_Trm' = 'Cd69', 
  'CD8_Trm' = 'Itgae', 
  'CD8_Trm' = 'Itga1', 
  'CD8_Tcm' = 'Cd28'
)

Macrophage_major_marker_ls <- c(
  # M1 经典促炎巨噬
  'M1' = 'Nos2', 
  'M1' = 'Il12b', 
  'M1' = 'Fcgr1', 
  'M1' = 'Cd80', 
  # Cxcr10 小鼠不存在，删除
  'M1' = 'Il23a', 
  'M1' = 'Cxcl9', 
  'M1' = 'Cxcl10', 
  'M1' = 'Cxcl11', 
  'M1' = 'Cd86', 
  'M1' = 'Il1a', 
  'M1' = 'Il1b', 
  'M1' = 'Il6', 
  'M1' = 'Tnf', 
  'M1' = 'H2-Ab1', # 替代人MHCII Mhcii
  'M1' = 'Ccl5', 
  'M1' = 'Irf5', 
  'M1' = 'Irf1', 
  'M1' = 'Cd40', 
  'M1' = 'Kynu', 
  'M1' = 'Ccr7', 
  'M1' = 'Vcan', 
  'M1' = 'Cd300e', 
  'M1' = 'Inhba', 
  'M1' = 'Tnip3', 
  
  # M2 抗炎/修复巨噬
  'M2' = 'Arg1', 
  'M2' = 'Il10', 
  'M2' = 'Fcgr2b', # Cd32 对应小鼠Fcgr2b
  'M2' = 'Cd163', 
  'M2' = 'Fcer2a', 
  'M2' = 'Cd200r1', 
  'M2' = 'Pdcd1lg2', 
  'M2' = 'Cd274', 
  'M2' = 'Marco', 
  'M2' = 'Csf1r', 
  'M2' = 'Mrc1', 
  'M2' = 'Il1rn', 
  'M2' = 'Il1r2', 
  'M2' = 'Il4r', 
  'M2' = 'Ccl4', 
  'M2' = 'Ccl13', 
  'M2' = 'Ccl20', 
  'M2' = 'Ccl17', 
  'M2' = 'Ccl22', 
  'M2' = 'Ccl24', 
  'M2' = 'Lyve1', 
  'M2' = 'Vegfa', 
  'M2' = 'Vegfb', 
  'M2' = 'Vegfc', 
  'M2' = 'Vegfd', 
  'M2' = 'Egf', 
  'M2' = 'Ctsa', 
  'M2' = 'Ctsb', 
  'M2' = 'Cstc', 
  'M2' = 'Ctsd', 
  'M2' = 'Tgfb1', 
  'M2' = 'Cd209a', 
  'M2' = 'Cd163l1', 
  'M2' = 'Folr2', 
  
  # M0 静息巨噬通用标记
  'M0' = 'Cd68', 
  'M0' = 'C1qa', 
  'M0' = 'C1qb', 
  'M0' = 'Selenop', 
  'M0' = 'Lyz2', # 小鼠主要溶菌酶Lyz2，Lyz几乎不表达
  'M0' = 'H2-Dpb1', # 替换人Hla-Dpb1
  'M0' = 'H2-Dpa1', # 替换人Hla-Dpa1
  'M0' = 'Aif1', 
  
  # TAM 肿瘤相关巨噬
  'TAM' = 'Spp1', 
  'TAM' = 'C1qc', # 修正原Cqc笔误
  
  # 树突细胞 DC
  "DC" = "Cd83",
  
  "pDC" = "Siglech", # 小鼠pDC核心标记，Lilra4人特有小鼠无
  "pDC" = "Irf4", 
  "pDC" = "Irf7", 
  "pDC" = "Irf8", 
  
  "cDC1" = "Xcr1",
  "cDC1" = "Clec9a",
  
  "cDC2" = "Cd1c", 
  "cDC2" = "Fcer1a", 
  "cDC2" = "Clec10a",
  
  "mregDC" = "Lamp3"
  #"mregDC" = "Ccr7"
)

Fibro_major_marker_ls <- c(
  'NF' = 'Igfbp6', 
  #'NF' = 'Pi16', 
  'NF' = 'Gpx3', 
  'NF' = 'Tnxb', 
  
  'CAF' = 'Fap', 
  'CAF' = 'Pdgfra', 
  'CAF' = 'Pdgfrb', 
  'CAF' = 'Hsf1', 
  'CAF' = 'Stat3', 
  'CAF' = 'Myc', 
  'CAF' = 'Cd10', 
  'CAF' = 'Gpr77', 
  
  'mCAF' = 'Mmp11', 
  'mCAF' = 'Myl9', 
  'mCAF' = 'Postn', 
  'mCAF' = 'Tpm2', 
  'mCAF' = 'Tagln', 
  'mCAF' = 'Hopx', 
  'mCAF' = 'Col1a2', 
  
  'iCAF' = 'Pla2g2a', 
  'iCAF' = 'Cfd', 
  'iCAF' = 'Cd34', 
  'iCAF' = 'Cxcl14', 
  #'iCAF' = 'Pi16', 
  'iCAF' = 'Prg4', 
  'iCAF' = 'Apod', 
  'iCAF' = 'Dpt', 
  'iCAF' = 'Cxcl12', 
  'iCAF' = 'Cxcl1', 
  'iCAF' = 'Il6', 
  'iCAF' = 'Ccl2', 
  
  'vCAF' = 'Notch3', 
  'vCAF' = 'Col18a1', 
  'vCAF' = 'Acta2', 
  'vCAF' = 'Mcam', 
  'vCAF' = 'Rgs5', 
  'vCAF' = 'Myh11', 
  'vCAF' = 'Adirf', 
  'vCAF' = 'Rergl', 
  
  'tCAF' = 'Pdpn', 
  'tCAF' = 'Mme', 
  'tCAF' = 'Ndrg1', 
  'tCAF' = 'Tmem158', 
  'tCAF' = 'Vegfa', 
  'tCAF' = 'Eno1', 
  'tCAF' = 'Gapdh', 
  'tCAF' = 'Caix', 
  'tCAF' = 'Hsph1', 
  'tCAF' = 'Hsp90aa1', 
  
  'ifnCAF' = 'Il32', 
  'ifnCAF' = 'Ido1', 
  'ifnCAF' = 'Cxcl9', 
  'ifnCAF' = 'Cxcl10', 
  'ifnCAF' = 'Cxcl11', 
  
  'apCAF' = 'Hla-Dra', 
  'apCAF' = 'Hla-Drb1', 
  'apCAF' = 'Cd74', 
  'apCAF' = 'Iglc2', 
  'apCAF' = 'Igkc', 
  'apCAF' = 'Cxcr4', 
  
  'rCAF' = 'Ccl19', 
  'rCAF' = 'Ccl21', 
  
  'dCAF' = 'Mki67', 
  'dCAF' = 'Tuba1b', 
  
  'meCAF' = 'Slc2a1', 
  'meCAF' = 'Bnip3', 
  'meCAF' = 'Adm', 
  'meCAF' = 'Angptl4'
)


Idents(TLS)<-'RNA_snn_res.0.6'
TLS_markers <-
  FindAllMarkers(
    JoinLayers(TLS),
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = TRUE
  )
write.csv(TLS_markers, file = "TLS_subtype/TLS_all_markers.csv", row.names = TRUE)

TLS_top5_markers <- 
  TLS_markers %>% 
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

names(TLS_top5_markers) <- sub('.*\\.','',names(TLS_top5_markers))


# dotplot
# 如果已做其他细胞类型的亚群分析，这里就可以顺便画一下其他类型的细胞亚群基因，比如Fibro_major_marker_ls2，否则只用画T细胞
pdf(file.path(TLS_out_dir, str_glue("dotplot dims{dims} major_marker 0.8.pdf")), width = 25, height = 6)
DotPlot(TLS, group.by = "RNA_snn_res.0.8", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(TLS, group.by = "RNA_snn_res.0.8", features = TLS_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(TLS, group.by = "RNA_snn_res.0.8", features = B_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(TLS, group.by = "RNA_snn_res.0.8", features = T_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(TLS, group.by = "RNA_snn_res.0.8", features = Macrophage_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(TLS, group.by = "RNA_snn_res.0.8", features = Fibro_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(TLS, group.by = "RNA_snn_res.0.8", features = TLS_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(TLS, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()

pdf(file.path(TLS_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(TLS, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(TLS, group.by = "RNA_snn_res.0.6", features = TLS_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(TLS, group.by = "RNA_snn_res.0.6", features = B_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(TLS, group.by = "RNA_snn_res.0.6", features = T_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(TLS, group.by = "RNA_snn_res.0.6", features = Macrophage_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(TLS, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls)+RotatedAxis()+theme_90angle
DotPlot(TLS, group.by = "RNA_snn_res.0.6", features = TLS_top5_markers)+RotatedAxis()+theme_90angle
#DotPlot(TLS, group.by = "RNA_snn_res.0.6", features = Fibro_major_marker_ls2)+RotatedAxis()
dev.off()


saveRDS(TLS ,file = "TLS_subtype/TLS_cluster.rds")
TLS <- readRDS("TLS_subtype/TLS_cluster.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "Pan_B", # 
  '1'  = "Tfh", #
  '2'  = "FRC_like", # 
  '3'  = "Pan_B", # 
  '4'  = "GC_B", # 
  '5'  = "Pan_B",#  
  '6'  = "Epithelial", #  
  '7'  = "M2", #  
  '8'  = "TCM", #  
  '9'  = "HEV_like", #  
  '10'  = "Plasma", #  
  '11'  = "HEV_like", #  
  '12'  = "cDC1", #  
  '13'  = "Smooth muscle cell", #  
  '14'  = "Th17", #  
  '15'  = "Tfh"
)

Idents(TLS) <- "RNA_snn_res.0.6"




TLS <- RenameIdents(TLS, cell_cluster)
TLS$subtype <- Idents(TLS)


pdf(file.path(TLS_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(TLS, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

pdf(file.path(TLS_out_dir, str_glue("Dimplot dim{dims} annotation split.pdf")), height = 6, width = 15)
DimPlot(TLS, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

saveRDS(TLS ,file = "TLS_subtype/TLS_allannotation.rds")



#画一个注释完的点图
TLS_major_marker_ls_2 <- c(
  # Tfh 滤泡辅助T细胞
  "Tfh" = "Cd4",
  "Tfh" = "Icos",
  "Tfh" = "Pdcd1",
  "Tfh" = "Il21",
  "Tfh" = "Maf",
  "Tfh" = "Sh2d1a",
  "Tfh" = "Batf",
  "Tfh" = "Slamf6",
  "Tfh" = "Cd40lg",
  
  # Pan-B 泛B细胞
  "Pan-B" = "Cd79a",
  "Pan-B" = "Cd79b",
  "Pan-B" = "Cd19",
  "Pan-B" = "Ms4a1",
  "Pan-B" = "Cd37",
  "Pan-B" = "Cd22",
  "Pan-B" = "Pax5",
  
  # GC B 生发中心B细胞
  "GC B" = "Aicda",
  "GC B" = "Bcl6",
  "GC B" = "Fas",
  "GC B" = "S1pr2",
  #"GC B" = "Cxcr4", # 可选标记
  
  
  # Plasmablast/Plasma 浆母/浆细胞
  "Plasma" = "Jchain",
  "Plasma" = "Sdc1",
  "Plasma" = "Prdm1",
  "Plasma" = "Irf4",
  "Plasma" = "Mzb1",
  #"Plasmablast/Plasma" = "Tnfrsf17", # 可选标记
  
  # FRC-like 成纤维网状细胞
  "FRC-like" = "Pdpn",
  "FRC-like" = "Pdgfra",
  "FRC-like" = "Ccl19",
  "FRC-like" = "Col1a1",
  "FRC-like" = "Col1a2",
  "FRC-like" = "Lum",
  
  # HEV-like 高内皮微静脉内皮细胞
  "HEV-like" = "Pecam1",
  "HEV-like" = "Cdh5",
  "HEV-like" = "Kdr",
  "HEV-like" = "Vwf",
  "HEV-like" = "Plvap",
  "HEV-like" = "Cd34",
  "HEV-like" = "Podxl",
  "HEV-like" = "Madcam1",
  
  
  'TCM' = 'Tcf7', 
  'TCM' = 'Ccr7', 
  
  'Th17' = 'Il23r', 
  'Th17' = 'Rorc', 
  
  'Epithelial' = 'Epcam',
  'Epithelial' = 'Krt20',
  
  'Smooth muscle cell' = 'Myh11',
  'Smooth muscle cell' = 'Tagln',
  
  'M2' = 'Cd68', 
  'M2' = 'Csf1r', 
  'M2' = 'Mrc1', 

  
  "cDC1" = "Xcr1",
  "cDC1" = "Clec9a"
)


pdf(file.path(TLS_out_dir, str_glue("dotplot TLS subtype major_marker.pdf")), width = 20, height = 6)
DotPlot(TLS, group.by = "subtype", features = TLS_major_marker_ls_2)+RotatedAxis()+theme_90angle
dev.off()

## cell props barplot  ----------------------------------------------------------
TLS_meta <- TLS[[]]

# 新分组显示顺序（左→右：Normal, Mild, Severe）
group_levels <- c("CON", "CRCFMT", "HDFMT")
# 对应配色（保持你指定色值不变）
group_col <- c(
  "CON"  = "#2E86AB",
  "CRCFMT"    = "#A23B72",
  "HDFMT"  = "#F18F01"
)


# 确保meta里的group是因子，顺序和颜色对应
#TLS_meta$Group <- recode(TLS_meta$group, !!!group_map)
TLS_meta$group <- factor(TLS_meta$group, levels = group_levels)

TLS_meta$sample <- factor(TLS_meta$sample,levels=c(TLS_meta$sample %>% unique()))

p_col_proj <-
  TLS_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  #scale_fill_brewer(palette = "Set2", name = "TLS Subtype") +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(TLS_out_dir, "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  TLS_meta %>%
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
  filename = file.path(TLS_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)


p_col_proj <-
  TLS_meta %>%
  group_by(sample_TLS, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample_TLS, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path(TLS_out_dir, "cell_freq_TLS.pdf"),
  height = 6,
  width  = 10
)

#改图
# 提取meta
TLS_meta <- TLS[[]]

# 新分组显示顺序（左→右：Normal, Mild, Severe）
group_levels <- c("CON", "CRCFMT", "HDFMT")
# 对应配色（保持你指定色值不变）
group_col <- c(
  "CON"  = "#2E86AB",
  "CRCFMT"    = "#A23B72",
  "HDFMT"  = "#F18F01"
)


# 确保meta里的group是因子，顺序和颜色对应
#TLS_meta$Group <- recode(TLS_meta$group, !!!group_map)
TLS_meta$group <- factor(TLS_meta$group, levels = group_levels)


# ========== 2. 计算样本水平百分比（用于箱线图 & 均值堆叠图） ==========
# 每个样本内各细胞亚群百分比
data_long <- TLS_meta %>%
  group_by(sample, group, subtype) %>%
  summarise(Count = n(), .groups = "drop_last") %>%
  group_by(sample) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

# 保存原始百分比数据
write.csv(data_long, 
          file.path("TLS_subtype", "TLS_sample_percentage.csv"), 
          row.names = FALSE)

# ==================== 按细胞类型分组的箱线图 ====================
comp_list <- list(c("CON","HDFMT"), c("CON","CRCFMT"), c("HDFMT","CRCFMT"))
p <- ggplot(data_long, aes(x = group, y = Percentage, fill = group)) +
  geom_boxplot(width = 0.7, alpha = 0.8, outlier.shape = 19, outlier.size = 1.5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7, aes(color = group)) +
  facet_wrap(~subtype, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_col) +
  scale_color_manual(values = group_col) +
  labs(title = "TLS - Cell Type Percentages",
       x = "Group", y = "Percentage (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))+
  stat_compare_means(
    comparisons = comp_list,
    method = "wilcox.test",
    label = "p.format",
    size = 3
  )

# 保存+展示
ggsave(file.path("TLS_subtype","TLS_Boxplot.pdf"), p, width = 12, height = 14, dpi = 300)
print(p)


# 映射到空间切片
# 批量赋值 annotation1
all_data$annotation1 <- "Others"

all_data$annotation1 %>% table()
#all_data <- readRDS("output/data/all_data.rds")

# 创建映射：细胞名 → 精细TLS亚型
cell_subtype <- TLS$subtype
names(cell_subtype) <- colnames(TLS)

# 只把TLS替换成精细亚型，其他细胞保持原样
all_data$annotation1[colnames(TLS)] <- as.character(cell_subtype)

# 查看结果（正确）
table(all_data$annotation1)
table(all_data$annotation)
table(TLS$subtype)


TLS_spatial_out_dir <- "/mnt/NexData_file/home/lxy/case/XE2603162231-02/20260616_XE2603162231-02/TLS_subtype/spatial" 
dir.create(TLS_spatial_out_dir)

cells_boundary_ls <- readRDS("R_OUTPUT/data_for_analysis/cells_boundary_ls.rds")

head(cells_boundary_ls)
tail(cells_boundary_ls)


all_data_meta <- 
  all_data[[]] %>% 
  rownames_to_column('cell_id') 
head(all_data_meta)
tail(all_data_meta)

# 1. 定义需要展示的细胞类型（对应table里非Others）
keep <- c("cDC1","M2", "Epithelial", "FRC_like","Smooth muscle cell", "GC_B", "HEV_like", "Pan_B", "Plasma", "TCM", "Tfh", "Th17")

# 2. 提取IGV调色板前10个颜色，匹配10种细胞
igv_col <- ggsci::pal_igv()(length(keep))
# 命名映射：细胞名 = IGV颜色
cell_color_map <- setNames(igv_col, keep)
# 追加Others灰色
cell_color_map[["Others"]] <- "#dcdcdc"

cn_sp_p_ls <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){ # browser()
    anno_data   <- all_data_meta[c("sample", "group")] %>% distinct()
    group_info  <- anno_data$Group[anno_data$sample == sample]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    dat <- dat %>%
      mutate(plot_group = ifelse(annotation1 %in% keep, annotation1, "Others"))
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation1, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() + 
      # 自定义颜色：为keep列表中的细胞分配颜色，"Others"设为灰色
      scale_fill_manual(
        values = cell_color_map,
        limits = c(keep, "Others")
      ) +
      theme_void() +
      labs(fill = NULL, title = sample)
    p
  })

cn_sp_p_ls %>% 
  iwalk(function(p, name){
    if(all(is.na(p$data$vertex_x))) return(NULL)
    if(any(is.na(p$data$vertex_x))) warning("vertex_x has NAs ...", .immediate = T)
    pdf(
      file.path(TLS_spatial_out_dir, str_glue("TLS spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


library(patchwork)

# 方法1：将所有图拼接成一个大的网格
combined_plot <- wrap_plots(cn_sp_p_ls, ncol = 4)  # 每行4张图

# 保存
ggsave(
  file.path(TLS_out_dir, "TLS spatial plot.pdf"),
  plot = combined_plot,
  width = 40,   # 根据子图数量调整
  height = 40,
  limitsize = FALSE
)



TLS_spatial_out_dir_2 <- "/mnt/NexData_file/home/lxy/case/XE2603162231-02/20260616_XE2603162231-02/TLS_subtype/TLS_spatial" 
dir.create(TLS_spatial_out_dir_2)

TLS_meta <- 
  TLS[[]] %>% 
  rownames_to_column('cell_id') 
head(TLS_meta)
tail(TLS_meta)


cn_sp_p_ls_2 <- 
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){ # browser()
    anno_data   <- TLS_meta[c("sample", "group")] %>% distinct()
    group_info  <- anno_data$Group[anno_data$sample == sample]
    
    dat <-
      cells_boundary %>%
      right_join(TLS_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() + 
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL, title = sample)
    p
  })

cn_sp_p_ls_2 %>% 
  iwalk(function(p, name){
    if(all(is.na(p$data$vertex_x))) return(NULL)
    if(any(is.na(p$data$vertex_x))) warning("vertex_x has NAs ...", .immediate = T)
    pdf(
      file.path(TLS_spatial_out_dir_2, str_glue("TLS filter spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })

combined_plot_2 <- wrap_plots(cn_sp_p_ls_2, ncol = 4)  # 每行4张图

# 保存
ggsave(
  file.path(TLS_out_dir, "TLS filter spatial plot.pdf"),
  plot = combined_plot_2,
  width = 40,   # 根据子图数量调整
  height = 40,
  limitsize = FALSE
)

# 组间差异分析和通路富集
#nohup Rscript /mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/deg_gokegg/deg_enrich_volcano.R -i /mnt/NexData_file/home/lxy/case/XE2603162231-02/20260616_XE2603162231-02/TLS_subtype/TLS_allannotation.rds --species mouse -g group -c subtype --ident1 CRCFMT --ident2 CON -o CRCFMT_vs_CON_DEG &
#nohup Rscript /mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/deg_gokegg/deg_enrich_volcano.R -i /mnt/NexData_file/home/lxy/case/XE2603162231-02/20260616_XE2603162231-02/TLS_subtype/TLS_allannotation.rds --species mouse -g group -c subtype --ident1 HDFMT --ident2 CON -o HDFMT_vs_CON_DEG &
#nohup Rscript /mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/deg_gokegg/deg_enrich_volcano.R -i /mnt/NexData_file/home/lxy/case/XE2603162231-02/20260616_XE2603162231-02/TLS_subtype/TLS_allannotation.rds --species mouse -g group -c subtype --ident1 HDFMT --ident2 CRCFMT -o HDFMT_vs_CRCFMT_DEG &
  
#聚类，划分TLS成熟度
library(dplyr)
library(tibble)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

# 读取TLS元数据
TLS_meta <- TLS@meta.data %>% filter(!is.na(subtype))
subtype_order <- c(
  "Pan_B", "Tfh", "FRC_like", "GC_B", "Epithelial", 
  "M2", "TCM", "HEV_like", "Plasma", "cDC1", 
  "Smooth muscle cell", "Th17"
)
TLS_meta$subtype <- factor(TLS_meta$subtype, levels = subtype_order)

# 统计每个sample_TLS内各subtype细胞数量、占比
tls_count_df <- 
  TLS_meta %>% 
  group_by(sample_TLS, subtype) %>%
  summarise(count = n()) %>% # 不写.groups="drop"，保留sample_TLS分组
  mutate(freq = count/sum(count))

tls_count_df %>% head()


# 长表转宽矩阵，行=单个TLS微环境，列=各细胞亚型占比
tls_count_matrix <- 
  tls_count_df %>% 
  as.data.frame() %>% 
  dplyr::select(sample_TLS, freq, subtype) %>% 
  pivot_wider(names_from = subtype, values_from = freq) %>% 
  column_to_rownames('sample_TLS')

# 无该细胞亚型的TLS填充0
na_idx <- is.na(tls_count_matrix)
tls_count_matrix[na_idx] <- 0
tls_count_matrix

tls_count_matrix %>%
  rownames_to_column("sample_TLS") %>%
  write.csv("/mnt/NexData_file/home/lxy/case/XE2603162231-02/20260616_XE2603162231-02/TLS_subtype/tls_subtype_freq_matrix.csv", row.names = FALSE)

# 读取csv还原矩阵
#tls_count_matrix <- read.csv("tls_subtype_freq_matrix.csv", row.names = 1)

# Kmeans聚类，分2类
set.seed(1234)
tls_cluster <- kmeans(tls_count_matrix, centers = 2, iter.max = 100)
tls_cluster$cluster

anno_row <- data.frame(
  K_cluster = factor(tls_cluster$cluster)
)
rownames(anno_row) <- names(tls_cluster$cluster)

# 分组配色
anno_color <- list(K_cluster = c("1" = "#E64B35", "2" = "#4DBBD5"))

# 行注释
row_anno <- rowAnnotation(
  K_cluster = anno_row$K_cluster,
  col = anno_color,
  show_legend = TRUE
)
# 颜色映射
col_fun <- colorRamp2(c(0, 1), c("white", "red"))

pdf("/mnt/NexData_file/home/lxy/case/XE2603162231-02/20260616_XE2603162231-02/TLS_subtype/TLS_subtype_ComplexHeatmap_k2.pdf", width = 15, height = 11)
Heatmap(
  matrix = tls_count_matrix,
  name = "Proportion",
  col = col_fun,
  # 坐标轴
  show_row_names = TRUE,    # Y轴 sample_TLS
  show_column_names = TRUE, # X轴 subtype
  left_annotation = row_anno,
  split = anno_row$K_cluster, # 按Kmeans结果把样本分成上下两大块
  # 聚类
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_title = "TLS subtype proportion heatmap",
  use_raster = FALSE
)
dev.off()

# Kmeans聚类，分3类
set.seed(1234)
tls_cluster <- kmeans(tls_count_matrix, centers = 3, iter.max = 100)
tls_cluster$cluster

anno_row <- data.frame(
  K_cluster = factor(tls_cluster$cluster)
)
rownames(anno_row) <- names(tls_cluster$cluster)

# 分组配色
anno_color <- list(K_cluster = c("1" = "#E64B35", "2" = "#4DBBD5", "3" = "#f5f829"))

# 行注释
row_anno <- rowAnnotation(
  K_cluster = anno_row$K_cluster,
  col = anno_color,
  show_legend = TRUE
)
# 颜色映射
col_fun <- colorRamp2(c(0, 1), c("white", "red"))

pdf("/mnt/NexData_file/home/lxy/case/XE2603162231-02/20260616_XE2603162231-02/TLS_subtype/TLS_subtype_ComplexHeatmap_k3.pdf", width = 15, height = 11)
Heatmap(
  matrix = tls_count_matrix,
  name = "Proportion",
  col = col_fun,
  # 坐标轴
  show_row_names = TRUE,    # Y轴 sample_TLS
  show_column_names = TRUE, # X轴 subtype
  left_annotation = row_anno,
  split = anno_row$K_cluster, # 按Kmeans结果把样本分成上下两大块
  # 聚类
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_title = "TLS subtype proportion heatmap",
  use_raster = FALSE
)
dev.off()


#给Seurat对象TLS新增两列meta信息 ==========
# 1. 提取样本名-数字cluster对应表
cluster_df <- enframe(tls_cluster$cluster, name = "sample_TLS", value = "K_cluster")
# 2. 数字转TLS成熟分型
type_map <- c("1" = "S-TLS", "2" = "E-TLS", "3" = "P-TLS")
cluster_df$TLS_type <- type_map[as.character(cluster_df$K_cluster)]

# 3. 匹配sample_TLS，新增两列到TLS@meta.data
TLS@meta.data <- TLS@meta.data %>%
  left_join(cluster_df, by = "sample_TLS")

# 验证新增列
head(TLS@meta.data[, c("sample_TLS", "K_cluster", "TLS_type")])
table(TLS$K_cluster)
table(TLS$TLS_type)


## cell props barplot  ----------------------------------------------------------
TLS_meta <- TLS[[]]

# 新分组显示顺序（左→右：Normal, Mild, Severe）
group_levels <- c("CON", "CRCFMT", "HDFMT")
# 对应配色（保持你指定色值不变）
group_col <- c(
  "CON"  = "#2E86AB",
  "CRCFMT"    = "#A23B72",
  "HDFMT"  = "#F18F01"
)


# 确保meta里的group是因子，顺序和颜色对应
#TLS_meta$Group <- recode(TLS_meta$group, !!!group_map)
TLS_meta$group <- factor(TLS_meta$group, levels = group_levels)

TLS_meta$sample <- factor(TLS_meta$sample,levels=c(TLS_meta$sample %>% unique()))

p_col_proj <-
  TLS_meta %>%
  group_by(sample, TLS_type) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = TLS_type)) + 
  geom_col() +
  #scale_fill_brewer(palette = "Set2", name = "TLS type") +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(TLS_out_dir, "TLS_type_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  TLS_meta %>%
  group_by(group, TLS_type) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = TLS_type)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_col_proj,
  filename = file.path(TLS_out_dir, "TLS_type_freq_group.pdf"),
  height = 6,
  width  = 4
)

saveRDS(TLS ,file = "TLS_subtype/TLS_type_allannotation.rds")
# 组间差异分析和通路富集(by TLS_type)
#nohup Rscript /mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/deg_gokegg/deg_enrich_volcano.R -i /mnt/NexData_file/home/lxy/case/XE2603162231-02/20260616_XE2603162231-02/TLS_subtype/TLS_type_allannotation.rds --species mouse -g group -c TLS_type --ident1 CRCFMT --ident2 CON -o CRCFMT_vs_CON_TLS_type_DEG &
#nohup Rscript /mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/deg_gokegg/deg_enrich_volcano.R -i /mnt/NexData_file/home/lxy/case/XE2603162231-02/20260616_XE2603162231-02/TLS_subtype/TLS_type_allannotation.rds --species mouse -g group -c TLS_type --ident1 HDFMT --ident2 CON -o HDFMT_vs_CON_TLS_type_DEG &
#nohup Rscript /mnt/NexData_file/home/lxy/case/XE2603262474-01/stRNA_2/deg_gokegg/deg_enrich_volcano.R -i /mnt/NexData_file/home/lxy/case/XE2603162231-02/20260616_XE2603162231-02/TLS_subtype/TLS_type_allannotation.rds --species mouse -g group -c TLS_type --ident1 HDFMT --ident2 CRCFMT -o HDFMT_vs_CRCFMT_TLS_type_DEG &

