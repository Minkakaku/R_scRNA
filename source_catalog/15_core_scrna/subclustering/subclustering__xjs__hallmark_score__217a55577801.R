library(Seurat)        
library(msigdbr)      
library(tidyverse)    
library(pheatmap)     
library(RColorBrewer)  
library(cowplot)  

All_T_anno <- readRDS('T_subtype/All_T_anno_remove_Tumor3.rds')
All_T_anno$subtype %>% table()
DimPlot(All_T_anno, group.by = "subtype",split.by = 'group', label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()

# Hallmark打分 -----------
All_T_anno$group %>% table()
# hallmark打分
# 输出目录
out_dir <- 'T_Hallmark_test'
dir.create(out_dir, recursive = T, showWarnings = F)

## All_T_anno Tumor hallmark打分 --------------
h_geneset_ls <- list(
  Glycolysis = c("ABCB6","ADORA2B","ALDH7A1","ALDOA","ANGPTL4","ANKZF1","AURKA","B3GAT1","BIK","CACNA1H","CASP6","CD44","CDK1","CENPA","COL5A1","CTH","CXCR4","DLD","ECD","EGFR","EGLN3","ENO2","FUT8","G6PD","GALE","GALK1","GFPT1","GOT1","GUSB","GYS1","HAX1","HK2","HMMR","HS6ST2","IDH1","IDUA","IL13RA1","IRS2","KIF20A","KIF2A","LDHA","ME1","MERTK","MET","MPI","NT5E","P4HA1","PAM","PC","PFKP","PGK1","PKM","PKP2","PMM2","PPFIA4","PRPS1","PYGB","PYGL","RBCK1","SDC2","SDC3","SLC25A10","SLC25A13","SLC37A4","SOD1","SOX9","SRD5A3","STC2","STMN1","TGFA","TKTL1","TPBG","VEGFA","VLDLR"),
  Oxidative_phosphorylation = c("ACAA2","ACAT1","AIFM1","ATP5F1A","ATP5F1C","ATP6V1H","BAX","CASP7","CPT1A","CYCS","DLD","ECHS1","FDX1","FH","GLUD1","HADHA","HADHB","HSD17B10","HSPA9","HTRA2","IDH1","IDH2","IMMT","LDHA","LRPPRC","MAOB","MFN2","MTRF1","NDUFA3","NDUFAB1","NQO2","OGDH","PDHA1","PDHB","PHB2","POR","RHOT1","SDHA","SDHB","SLC25A11","SLC25A4","SUPV3L1","TIMM50","VDAC1","VDAC2"),
  Fatty_acid_metabolism = c("ACAA2","ACSL4","ADIPOR2","ALAD","ALDH3A2","ALDOA","AOC3","APEX1","AQP7","BLVRA","BMPR1B","CA4","CBR1","CD1D","CD36","CPT1A","CRYZ","CYP1A1","DHCR24","DLD","ECHS1","ENO2","EPHX1","FABP2","FASN","FH","GAD2","GSTZ1","HADHB","HPGD","HSD17B10","IDH1","INMT","LDHA","MAOA","ME1","NBN","PCBD1","PDHA1","PDHB","PPARA","PTPRG","REEP6","SDHA","SLC22A5","TDO2","UGDH")
  #Ferroptosis = c("ACSL4","ACACA","FADS2","GPX4","GCLC","GCLM","NFE2L2","NQO1","SLC7A11","SLC3A2","TP53","ACSF2","AIFM2","ATF3","ATF4","CARS","CDO1","CHAC1","DHODH","DPP4","EMC2","EMC3","EMC6","FTH1","FTHL17","HMGCR","HMOX1","HSPB1","IREB2","LPCAT3","MDM2","MT1G","MT1X","NCOA4","PCBP1","PCBP2","PRNP","PTGS2","RPL8","SAT1","SCL7A5","SIRT1","SLC1A5","SLC40A1","STEAP3","TFR1")
)


sc_All_T_anno_score <-
  JoinLayers(All_T_anno) %>%
  AddModuleScore(
    features = h_geneset_ls, #传入3个Hallmark基因集
    name = names(h_geneset_ls)
  )
#合并Seurat多层数据 → 给每个细胞计算3个Hallmark通路的得分


sc_All_T_anno_score %>% head
#查看打分结果

# msigdb anno
metadata <- sc_All_T_anno_score[[]]

# 获取所有的 group 类型
groups <- unique(metadata$group)

# 循环每个 group，分别计算均值 + 画热图 + 输出PDF
for (g in groups) {#browser()
  
  # 1. 筛选当前 group 的细胞
  current_data <- metadata %>%
    dplyr::filter(group == g) %>%
    dplyr::select(subtype, starts_with(names(h_geneset_ls)))
  
  # 2. 按 subtype 计算通路平均分
  h_score_anno <- current_data %>%
    group_by(subtype) %>%
    summarise(across(everything(), mean), .groups = "drop") %>%
    column_to_rownames("subtype") %>%
    t()
  
  # 3. 生成热图（关键：用grid.draw输出，避免重叠）
  p <- pheatmap(
    h_score_anno,
    scale = 'row',
    labels_row = names(h_geneset_ls),
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(Inf, "RdBu")))(50),
    border_color = NA,
    cluster_cols = F,
    cluster_rows = F,
    silent = TRUE
  )
  
  # 4. 正确写入PDF（每个group一个独立文件，避免设备混乱）
  pdf_file <- file.path(out_dir, paste0("AddModuleScore_subtype_group_", g, ".pdf"))
  pdf(pdf_file, height = 6, width = 8)
  grid::grid.draw(p$gtable)  # 用grid.draw直接输出gtable，避免plot_grid累积
  dev.off()  # 立即关闭设备，保证文件完整性
}





# ACSL4 vlnplot -----------
all_data <- readRDS("~/all_data_remove_Tumor3_annotation_full_pure.rds")
DimPlot(all_data, group.by = "annotation", label = TRUE, repel = TRUE) + NoLegend()

Myeloid <- readRDS("~/Myeloid_subtype/Myeloid_pure.rds")
DimPlot(Myeloid, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()


VlnPlot(All_T_anno, features = 'ACSL4', group.by = 'subtype') 
ggsave(file.path(out_dir, 'Vlnplot_ACSL4_T.pdf'), width = 6, height = 6)
VlnPlot(subset(JoinLayers(All_T_anno), ACSL4 > 0), features = 'ACSL4', group.by = 'subtype') #去掉表达量 = 0 的细胞，再画一次小提琴图
ggsave(file.path(out_dir, 'Vlnplot_ACSL4_T_move0.pdf'), width = 6, height = 6)

VlnPlot(all_data, features = 'ACSL4', group.by = 'annotation') 
ggsave(file.path(out_dir, 'Vlnplot_ACSL4.pdf'), width = 6, height = 6)
VlnPlot(subset(JoinLayers(all_data), ACSL4 > 0), features = 'ACSL4', group.by = 'annotation') #去掉表达量 = 0 的细胞，再画一次小提琴图
ggsave(file.path(out_dir, 'Vlnplot_ACSL4_move0.pdf'), width = 6, height = 6)

VlnPlot(Myeloid, features = 'ACSL4', group.by = 'subtype') 
ggsave(file.path(out_dir, 'Vlnplot_ACSL4_Myeloid.pdf'), width = 8, height = 6)
VlnPlot(subset(JoinLayers(Myeloid), ACSL4 > 0), features = 'ACSL4', group.by = 'subtype') #去掉表达量 = 0 的细胞，再画一次小提琴图
ggsave(file.path(out_dir, 'Vlnplot_ACSL4_Myeloid_move0.pdf'), width = 8, height = 6)

