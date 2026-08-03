setwd('/home/zwl/Project/FWAE20242404/require2/')
library(ggplot2)
library(dplyr)
library(Seurat)
rm(list = ls())
gc()

Seuset <- readRDS(file = './Dibro_Macro_celltype.rds')
Seurat::DimPlot(Seuset,reduction = 'umap',group.by = 'sub_celltype')


gene_exp = FetchData(object = Seuset, vars = c('Ecm1','Lrp1'))
data  = cbind(gene_exp,Seuset@meta.data[rownames(gene_exp),c('sample','sub_celltype')])

# 删除包含0的行
data_B = data %>%
  filter(!if_any(everything(), ~ . == 0))
p = ggplot(data = data_B,aes(x=Ecm1,y=Lrp1,fill = sub_celltype,colour = sub_celltype,shape = sample))+geom_point()+
  facet_wrap(sample~sub_celltype,nrow = 2)+
  geom_smooth(method = "lm", col = "blue")
ggsave(filename = './scarplot_Ecm1_Lrp1.png',p,w=10,h=5)

library(rstatix)
# data_B %>% group_by(sample) %>% rstatix::t_test(MARCKSL1~CD69)
res = data.frame(P_value = NULL ,r = NULL)
for (i in unique(data_B$sample)) {
  for (j in unique(data_B$sub_celltype)) {
    data_B1 = data_B %>% dplyr::filter(sample == i,sub_celltype == j)
    if (dim(data_B1)[1]>2) {
      print(dim(data_B1))
      correlation_test  = cor.test(data_B1$Ecm1,data_B1$Lrp1,method = 'pearson')
      res = rbind(res,c(correlation_test$p.value,correlation_test$estimate,i,j))
      print(paste0(i,' ',j,': p =',correlation_test$p.value,'  r:',correlation_test$estimate))
    }
  }
}
colnames(res)<- c('p_value','r_cor','sample','celltype')
write.csv(file = 'Ecm1_Lrp1_cor_pvalue.csv',res)
#######################################将SMC Ecm1和Macro Lrp1进行相关性####################################
# 将Ecm1和Lrp1按照 sample和 celltype 分类，然后随机合并取平均得到30个bin
set.seed(123)
# 设置要得到的组数
num_groups <- 30
all_avg = data.frame()
for (ni in unique(data$sample)) {
  for (nj in unique(data$sub_celltype)) {
    original_data = data %>% dplyr::filter(sample == ni,sub_celltype == nj)%>% select(Ecm1,Lrp1)
    # 随机打乱原始数据的顺序
    shuffled_data <- sample(original_data)
    
    # 计算每组应该包含的数字个数（可能不是整数，所以取整数部分）
    group_size <- dim(shuffled_data)[1] %/% num_groups  # 整除得到每组的基本大小
    extra_numbers <- dim(shuffled_data)[1] %% num_groups  # 计算额外的数字个数
    
    # 初始化一个列表来存储分组的数据
    groups <- split(shuffled_data, rep(1:num_groups, each = group_size, length.out = dim(shuffled_data)[1] ))
    
    # 如果存在额外的数字，则随机分配到各组中
    if (extra_numbers > 0) {
      extra_indices <- sample(1:num_groups, extra_numbers)
      for (i in 1:extra_numbers) {
        groups[[extra_indices[i]]] <- c(groups[[extra_indices[i]]], shuffled_data[(group_size * num_groups + i),])
      }
    }
    
    calculate_averages_df <- function(df) {
      avg1 <- mean(df$Ecm1, na.rm = TRUE)  # 假设可能存在 NA 值，使用 na.rm = TRUE 来忽略它们
      avg2 <- mean(df$Lrp1, na.rm = TRUE)
      return(c(avg1, avg2))
    }
    
    # 计算每组的平均值
    group_averages <- sapply(groups, calculate_averages_df) %>% t() %>% data.frame()
    colnames(group_averages) <- c('Ecm1','Lrp1')
    group_averages$sample = ni
    group_averages$celltype = nj
    
    all_avg = rbind(group_averages,all_avg)
}}

for (si in c('Air','O2')) {
  for (Ei in c("Lipofibroblast","Ecm1 Myofib","Pericyte")) {
    Edata = all_avg %>% dplyr::filter(sample == si,celltype == Ei)
    for (Lj in c("Mono","Alveolar Macro","interstitial Macro","M1 Macro")) {
      Fdata = all_avg %>% dplyr::filter(sample == si,celltype == Lj)
      #不填充零了，合并试试
      xE = Edata$Ecm1
      yL = Fdata$Lrp1
      P_data = data.frame(x = xE, y = yL)
      ggplot(data = P_data,aes(x=x,y=y))+geom_point()+
        geom_smooth(method = "lm", col = "blue")
      correlation_test  = cor.test(xE,yL,method = 'pearson')
      if (correlation_test$p.value<0.05){
        print(paste0(si,' ',Ei,' vs ',Lj,': p =',correlation_test$p.value,'  r:',correlation_test$estimate))
      }
      
    }
  }
}


# for (si in c('Air','O2')) {
#   for (Ei in c("Lipofibroblast","Ecm1 Myofib","Pericyte")) {
#     Edata = data %>% dplyr::filter(sample == si,sub_celltype == Ei)
#     for (Lj in c("Mono","Alveolar Macro","interstitial Macro","M1 Macro")) {
#       Fdata = data %>% dplyr::filter(sample == si,sub_celltype == Lj)
#       #不填充零了，合并试试
#       if (dim(Edata)[1]>dim(Fdata)[1]){
#         xE = Edata$Ecm1
#         yL = c(Fdata$Lrp1,rep(0,dim(Edata)[1]-dim(Fdata)[1]))
#       }else if (dim(Edata)[1]<dim(Fdata)[1]){
#         xE = c(Edata$Ecm1,rep(0,dim(Fdata)[1]-dim(Edata)[1]))
#         yL = Fdata$Lrp1
#       }else{
#         xE = Edata$Ecm1
#         yL = Fdata$Lrp1
#       }
#       P_data = data.frame(x = xE, y = yL)
#       ggplot(data = P_data,aes(x=x,y=y))+geom_point()+
#         geom_smooth(method = "lm", col = "blue")
#       # correlation_test  = cor.test(xE,yL,method = 'pearson')
#       print(paste0(si,' ',Ei,' vs ',Lj,': p =',correlation_test$p.value,'  r:',correlation_test$estimate))
#     }
#   }
# }



#####################################################总结##################################################################
#上面方法1：同一细胞类型计算的没有意义
#    方法2：不同细胞类型的合成bin的方法不好解释
#下面尝试cellphonedb：先将小鼠转成人类基因


obj <- readRDS(file = './Dibro_Macro_celltype.rds')
Seurat::DimPlot(obj,reduction = 'umap',group.by = 'sub_celltype')

library(homologene)
gene = mouse2human(rownames(obj))[,1:2]%>% distinct(mouseGene, .keep_all = TRUE) %>% tibble::column_to_rownames('mouseGene')
gene_martix = as.matrix(FetchData(object = obj, vars = gene$mouseGene))
#inner = intersect(colnames(gene_martix),gene$mouseGene)
colnames(gene_martix) <- gene[colnames(gene_martix),]

write.table(obj$sub_celltype, 'cellphonedb_meta.txt', sep='\t', quote=F,col.names = F,row.names=T)
write.table(gene_martix, 'cellphonedb_count.txt', sep='\t', quote=F) 



