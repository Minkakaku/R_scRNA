
library(Seurat)
library(tidyverse)

all_data@meta.data %>% head()

p1 <- DimPlot(all_data, 
              group.by = "annotation_full", 
              split.by = "project", 
              cols = plan,
              raster = FALSE, 
              ncol = 3, 
              label = T
              )

ggsave(p1, file = "2. dimension reduction and cluster/Ratio/umap_split.png", width = 12, height = 6, dpi = 300)




# 写出csv导入loupe browser
meta_data <- read.table("2. dimension reduction and cluster/Ratio/meta_data.txt", header = T, sep = "\t", check.names = F)
colnames(meta_data)
# meta_data$Cell <- gsub("^[1]\\d+_", "", meta_data$Cell)
meta_data$Cell <- sub("^[^_]*_", "", meta_data$Cell)

pre1_meta_data <- meta_data[meta_data$project == "Pre_1", ]
pre2_meta_data <- meta_data[meta_data$project == "Pre_2", ]
pre3_meta_data <- meta_data[meta_data$project == "Pre_3", ]
post1_meta_data <- meta_data[meta_data$project == "Post_1", ]
post2_meta_data <- meta_data[meta_data$project == "Post_2", ]
post3_meta_data <- meta_data[meta_data$project == "Post_3", ]

pre1_meta_data <- pre1_meta_data[, c(1, 18)]
pre2_meta_data <- pre2_meta_data[, c(1, 18)]
pre3_meta_data <- pre3_meta_data[, c(1, 18)]
post1_meta_data <- post1_meta_data[, c(1, 18)]
post2_meta_data <- post2_meta_data[, c(1, 18)]
post3_meta_data <- post3_meta_data[, c(1, 18)]

colnames(pre1_meta_data) <- c("cell_id","group")
colnames(pre2_meta_data) <- c("cell_id","group")
colnames(pre3_meta_data) <- c("cell_id","group")
colnames(post1_meta_data) <- c("cell_id","group")
colnames(post2_meta_data) <- c("cell_id","group")
colnames(post3_meta_data) <- c("cell_id","group")



write.table(pre1_meta_data, file = "2. dimension reduction and cluster/Ratio/pre1_annotation.csv", sep = ",", row.names = F, col.names = T, quote = F)
write.table(pre2_meta_data, file = "2. dimension reduction and cluster/Ratio/pre2_annotation.csv", sep = ",", row.names = F, col.names = T, quote = F)
write.table(pre3_meta_data, file = "2. dimension reduction and cluster/Ratio/pre3_annotation.csv", sep = ",", row.names = F, col.names = T, quote = F)
write.table(post1_meta_data, file = "2. dimension reduction and cluster/Ratio/post1_annotation.csv", sep = ",", row.names = F, col.names = T, quote = F)
write.table(post2_meta_data, file = "2. dimension reduction and cluster/Ratio/post2_annotation.csv", sep = ",", row.names = F, col.names = T, quote = F)
write.table(post3_meta_data, file = "2. dimension reduction and cluster/Ratio/post3_annotation.csv", sep = ",", row.names = F, col.names = T, quote = F)




