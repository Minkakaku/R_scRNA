dir_path <- 'Clonoscope_visualization'
dir.create(dir_path)

# ADPC_subclone_detection ---------
Cov_obj <- readRDS('/mnt/NexData_file/home/xjs/FWSP20242434/clonalscope_result/ADPC_subclone_detection/ADPC_Cov_obj.rds')

## Step3. visiualization --------------
clustering= Cov_obj$result_final$clustering
clustering2= Cov_obj$result_final$clustering2
result=Cov_obj$result_final$result
Zest=result$Zest
table(result$Zest)

# 调整resolution，mincell越大cluster数越少
# result=AssignCluster(clustering2, mincell = 1000)
# Zest=result$Zest
# table(result$Zest)
Zest <- readRDS(file.path('/mnt/NexData_file/home/xjs/FWSP20242434/clonalscope_result/ADPC_subclone_detection/Zest_ADPC_subclone_detection.rds'))
Zest %>% table()
Zest %>% head()
     
     
## Step4. Integrate with Spatial information (Optional) ------------------
library(Seurat)
vs <- readRDS('/mnt/NexData_file/home/szt/server203/20250318_FWSP20242434/result_8um/ADPC7_annotation.rds')

vs$orig.ident %>% table()

Zest2 <- 
  Zest %>% 
  as.data.frame() %>% 
  rownames_to_column('cell_id')
colnames(Zest2) <- c('cell_id', 'Zest')

vs_meta <- 
  vs@meta.data %>% 
  rownames_to_column('cell_id') %>% 
  left_join(Zest2, by = 'cell_id') %>% 
  mutate(
    subclone = ifelse(
      annotation_full != 'Luminal cell', 
      annotation_full, 
      Zest)
  ) %>% 
  column_to_rownames('cell_id')

vs <- AddMetaData(vs, metadata = vs_meta[c('Zest', 'subclone')])

Idents(vs) <- 'subclone'
library(gplots)
cols <- c(colors()[c(536, 652, 611, 463, 498, 71, 258, 84, 56, 26, 154, 59, 134, 78, 116, 85, 20, 259)])[1:length(Idents(vs) %>% levels())]
names(cols) <- Idents(vs) %>% levels()
pdf(file.path(dir_path, 'ADPC_subclone_detection_ct.pdf'), width = 8, height = 6)
SpatialDimPlot(vs, label.size = 3, pt.size.factor = 10, image.scale = "hires", cols=cols)
dev.off()    

vs$subclone2 <- case_when(
  vs$subclone == '1' ~ 'subclone_1', 
  vs$subclone == '16' ~ 'subclone_16', 
  vs$subclone == '28' ~ 'subclone_28', 
  vs$subclone == '40' ~ 'subclone_40', 
  TRUE ~ 'Others'
  )
vs$subclone2 %>% table()
vs$subclone %>% table()

Idents(vs) <- 'subclone2'
library(gplots)
cols <- c(colors()[c(536, 652, 611, 463, 498, 71, 258, 84, 56, 26, 154, 59, 134, 78, 116, 85, 20, 259)])[1:length(Idents(vs) %>% levels())]
names(cols) <- Idents(vs) %>% levels()
cols <- str_replace_all(cols, c('yellow' = 'lightgrey'))
names(cols) <- Idents(vs) %>% levels()
SpatialDimPlot(vs, label.size = 3, pt.size.factor = 10, image.scale = "hires", cols=cols, group.by = 'subclone2', image.alpha = 0.1)
ggsave(file.path(dir_path, 'ADPC_subclone_detection.pdf'), width = 8, height = 6)



# Annotation info output -------------------------------------------------------
vs$subclone %>% table()
cellType_anno_vs <- 
  vs@meta.data[c('project', 'subclone')] %>% 
  rownames_to_column('cell_id') %>% 
  mutate(
    Barcode = str_replace(cell_id, '^.*?_', ''), 
    group = subclone
  ) %>% 
  dplyr::select(Barcode, group)

cellType_anno_vs %>% 
  write.csv(file.path(dir_path, 'ADPC_subclone_ct.csv'), row.names = F)




# ADPC_subclone_detection with WES---------
Cov_obj <- readRDS('/mnt/NexData_file/home/xjs/FWSP20242434/clonalscope_result/ADPC_subclone_detection_WES/ADPC_Cov_obj_WES.rds')

## Step3. visiualization --------------
clustering= Cov_obj$result_final$clustering
clustering2= Cov_obj$result_final$clustering2
result=Cov_obj$result_final$result
Zest=result$Zest
table(result$Zest)

# 调整resolution，mincell越大cluster数越少
# result=AssignCluster(clustering2, mincell = 1000)
# Zest=result$Zest
# table(result$Zest)
Zest <- readRDS(file.path('/mnt/NexData_file/home/xjs/FWSP20242434/clonalscope_result/ADPC_subclone_detection_WES/Zest_ADPC_subclone_detection_WES.rds'))
Zest %>% table()
Zest %>% head()


## Step4. Integrate with Spatial information (Optional) ------------------
library(Seurat)
vs <- readRDS('/mnt/NexData_file/home/szt/server203/20250318_FWSP20242434/result_8um/ADPC7_annotation.rds')

vs$orig.ident %>% table()

Zest2 <- 
  Zest %>% 
  as.data.frame() %>% 
  rownames_to_column('cell_id')
colnames(Zest2) <- c('cell_id', 'Zest')

vs_meta <- 
  vs@meta.data %>% 
  rownames_to_column('cell_id') %>% 
  left_join(Zest2, by = 'cell_id') %>% 
  mutate(
    subclone = ifelse(
      annotation_full != 'Luminal cell', 
      annotation_full, 
      Zest)
  ) %>% 
  column_to_rownames('cell_id')

vs <- AddMetaData(vs, metadata = vs_meta[c('Zest', 'subclone')])

Idents(vs) <- 'subclone'
library(gplots)
cols <- c(colors()[c(536, 652, 611, 463, 498, 71, 258, 84, 56, 26, 154, 59, 134, 78, 116, 85, 20, 259)])[1:length(Idents(vs) %>% levels())]
names(cols) <- Idents(vs) %>% levels()
# pdf(file.path(dir_path, 'ADPC_subclone_detection_ct.pdf'), width = 8, height = 6)
# SpatialDimPlot(vs, label.size = 3, pt.size.factor = 10, image.scale = "hires", cols=cols)
# dev.off()    

vs$subclone2 <- case_when(
  vs$subclone == '1' ~ 'subclone_1', 
  vs$subclone == '19' ~ 'subclone_19', 
  vs$subclone == '23' ~ 'subclone_23', 
  vs$subclone == '53' ~ 'subclone_53', 
  vs$subclone == '55' ~ 'subclone_55', 
  TRUE ~ 'Others'
)
vs$subclone2 %>% table()
vs$subclone %>% table()

Idents(vs) <- 'subclone2'
library(gplots)
cols <- c(colors()[c(536, 652, 611, 463, 498, 71, 258, 84, 56, 26, 154, 59, 134, 78, 116, 85, 20, 259)])[1:length(Idents(vs) %>% levels())]
names(cols) <- Idents(vs) %>% levels()
cols <- str_replace_all(cols, c('yellow' = 'lightgrey'))
names(cols) <- Idents(vs) %>% levels()
SpatialDimPlot(vs, label.size = 3, pt.size.factor = 10, image.scale = "hires", cols=cols, group.by = 'subclone2', image.alpha = 0.1)
ggsave(file.path(dir_path, 'ADPC_subclone_detection_WES.pdf'), width = 8, height = 6)



# Annotation info output -------------------------------------------------------
vs$subclone %>% table()
cellType_anno_vs <- 
  vs@meta.data[c('project', 'subclone')] %>% 
  rownames_to_column('cell_id') %>% 
  mutate(
    Barcode = str_replace(cell_id, '^.*?_', ''), 
    group = subclone
  ) %>% 
  dplyr::select(Barcode, group)

cellType_anno_vs %>% 
  write.csv(file.path(dir_path, 'ADPC_subclone_ct_WES.csv'), row.names = F)




# scoring -------------
## ADPC 4cluster ---------- 
Zest <- readRDS(file.path('/mnt/NexData_file/home/xjs/FWSP20242434/clonalscope_result/ADPC_subclone_detection/Zest_ADPC_subclone_detection.rds'))
Zest %>% table()
Zest %>% head()


vs <- readRDS('/mnt/NexData_file/home/szt/server203/20250318_FWSP20242434/result_8um/ADPC7_annotation.rds')

vs$orig.ident %>% table()

Zest2 <- 
  Zest %>% 
  as.data.frame() %>% 
  rownames_to_column('cell_id')
colnames(Zest2) <- c('cell_id', 'Zest')

vs_meta <- 
  vs@meta.data %>% 
  rownames_to_column('cell_id') %>% 
  left_join(Zest2, by = 'cell_id') %>% 
  mutate(
    subclone = ifelse(
      annotation_full != 'Luminal cell', 
      'Other cell type', 
      Zest)
  ) %>% 
  column_to_rownames('cell_id')

vs <- AddMetaData(vs, metadata = vs_meta[c('Zest', 'subclone')])

vs$subclone2 <- str_replace_all(vs$subclone, '^', 'subclone_')
vs$subclone2 <- str_replace_all(vs$subclone2, 'subclone_Other cell type', 'Other cell type')

vs$subclone2 %>% table()
vs$subclone %>% table()

basal <- c('KRT5', 'KRT14' ,'TP63' ,'ITGA6','KRT17', 'KRT15' ,'DKK1' ,'CAV1', 'CAV2' ,'CYR61')
luminal <- c('NKX3-1','KRT8','KRT18','CD24','PSCA','CLDN3','AMACR','ERG','ALDH1A3')
ne <- c('ASCL1','CHGA','CHGB','SYP','FOXA2','NKX2-1','ENO2','MYCN','POU3F2','NCAM1','INSM1','EZH2','SOX2','SIAH2')
club <- c('MMP7', 'PIGR', 'LTF', 'CP', 'KRT7', 'SLPI')

pathway<-list(basal, luminal, ne, club)
names(pathway)<-c('basal', 'luminal', 'ne', 'club')

n_na_idx <- pathway %>% map_lgl(function(x) all(is.na(x))) %>% {! .}%>%  which
pathway <- pathway[n_na_idx] %>% map(function(x) x[!is.na(x)])


for ( i in seq_along(names(pathway))){
  vs <- AddModuleScore(JoinLayers(vs), pathway[i], name =  names(pathway)[i],nbin=5,ctrl=10)
}

p <- VlnPlot(vs,
             features = c('basal1', 'luminal1', 'ne1', 'club1'),
             group.by='subclone2', 
             ncol = 4, 
             pt.size = 0)
p

ggsave(
  p,
  filename = file.path(dir_path, 'ADPC subclone score - 4cluster.pdf'),
  height = 5,
  width  = 15
) 




## ADPC 19cluster ---------- 
Zest <- readRDS(file.path('/mnt/NexData_file/home/xjs/FWSP20242434/clonalscope_result/ADPC_subclone_detection/Zest_ADPC_subclone_detection_mincell_500_cluster19.rds'))
Zest %>% table()
Zest %>% head()


vs <- readRDS('/mnt/NexData_file/home/szt/server203/20250318_FWSP20242434/result_8um/ADPC7_annotation.rds')

vs$orig.ident %>% table()

Zest2 <- 
  Zest %>% 
  as.data.frame() %>% 
  rownames_to_column('cell_id')
colnames(Zest2) <- c('cell_id', 'Zest')

vs_meta <- 
  vs@meta.data %>% 
  rownames_to_column('cell_id') %>% 
  left_join(Zest2, by = 'cell_id') %>% 
  mutate(
    subclone = ifelse(
      annotation_full != 'Luminal cell', 
      'Other cell type', 
      Zest)
  ) %>% 
  column_to_rownames('cell_id')

vs <- AddMetaData(vs, metadata = vs_meta[c('Zest', 'subclone')])

vs$subclone2 <- str_replace_all(vs$subclone, '^', 'subclone_')
vs$subclone2 <- str_replace_all(vs$subclone2, 'subclone_Other cell type', 'Other cell type')

vs$subclone2 %>% table()
vs$subclone %>% table()

for ( i in seq_along(names(pathway))){
  vs <- AddModuleScore(JoinLayers(vs), pathway[i], name =  names(pathway)[i],nbin=5,ctrl=10)
}

p <- VlnPlot(vs,
             features = c('basal1', 'luminal1', 'ne1', 'club1'),
             group.by='subclone2', 
             ncol = 4, 
             pt.size = 0)
p

ggsave(
  p,
  filename = file.path(dir_path, 'ADPC subclone score - 19cluster.pdf'),
  height = 5,
  width  = 20
) 




## ADPC WES 6cluster ---------- 
Zest <- readRDS(file.path('/mnt/NexData_file/home/xjs/FWSP20242434/clonalscope_result/ADPC_subclone_detection_WES/Zest_ADPC_subclone_detection_mincell_800_cluster6.rds'))
Zest %>% table()
Zest %>% head()


vs <- readRDS('/mnt/NexData_file/home/szt/server203/20250318_FWSP20242434/result_8um/ADPC7_annotation.rds')

vs$orig.ident %>% table()

Zest2 <- 
  Zest %>% 
  as.data.frame() %>% 
  rownames_to_column('cell_id')
colnames(Zest2) <- c('cell_id', 'Zest')

vs_meta <- 
  vs@meta.data %>% 
  rownames_to_column('cell_id') %>% 
  left_join(Zest2, by = 'cell_id') %>% 
  mutate(
    subclone = ifelse(
      annotation_full != 'Luminal cell', 
      'Other cell type', 
      Zest)
  ) %>% 
  column_to_rownames('cell_id')

vs <- AddMetaData(vs, metadata = vs_meta[c('Zest', 'subclone')])

vs$subclone2 <- str_replace_all(vs$subclone, '^', 'subclone_')
vs$subclone2 <- str_replace_all(vs$subclone2, 'subclone_Other cell type', 'Other cell type')

vs$subclone2 %>% table()
vs$subclone %>% table()

for ( i in seq_along(names(pathway))){
  vs <- AddModuleScore(JoinLayers(vs), pathway[i], name =  names(pathway)[i],nbin=5,ctrl=10)
}

p <- VlnPlot(vs,
             features = c('basal1', 'luminal1', 'ne1', 'club1'),
             group.by='subclone2', 
             ncol = 4, 
             pt.size = 0)
p

ggsave(
  p,
  filename = file.path(dir_path, 'ADPC subclone score - 6cluster - WES.pdf'),
  height = 5,
  width  = 16
) 




## ADPC WES 6cluster ---------- 
Zest <- readRDS(file.path('/mnt/NexData_file/home/xjs/FWSP20242434/clonalscope_result/ADPC_subclone_detection_WES/Zest_ADPC_subclone_detection_mincell_500_cluster15.rds'))
Zest %>% table()
Zest %>% head()


vs <- readRDS('/mnt/NexData_file/home/szt/server203/20250318_FWSP20242434/result_8um/ADPC7_annotation.rds')

vs$orig.ident %>% table()

Zest2 <- 
  Zest %>% 
  as.data.frame() %>% 
  rownames_to_column('cell_id')
colnames(Zest2) <- c('cell_id', 'Zest')

vs_meta <- 
  vs@meta.data %>% 
  rownames_to_column('cell_id') %>% 
  left_join(Zest2, by = 'cell_id') %>% 
  mutate(
    subclone = ifelse(
      annotation_full != 'Luminal cell', 
      'Other cell type', 
      Zest)
  ) %>% 
  column_to_rownames('cell_id')

vs <- AddMetaData(vs, metadata = vs_meta[c('Zest', 'subclone')])

vs$subclone2 <- str_replace_all(vs$subclone, '^', 'subclone_')
vs$subclone2 <- str_replace_all(vs$subclone2, 'subclone_Other cell type', 'Other cell type')

vs$subclone2 %>% table()
vs$subclone %>% table()

for ( i in seq_along(names(pathway))){
  vs <- AddModuleScore(JoinLayers(vs), pathway[i], name =  names(pathway)[i],nbin=5,ctrl=10)
}

p <- VlnPlot(vs,
             features = c('basal1', 'luminal1', 'ne1', 'club1'),
             group.by='subclone2', 
             ncol = 4, 
             pt.size = 0)
p

ggsave(
  p,
  filename = file.path(dir_path, 'ADPC subclone score - 15cluster - WES.pdf'),
  height = 5,
  width  = 20
) 
