cds <- readRDS("/mnt/NexData_file/home/wmx/Project/SC2605263939-01/shouhou3/Monocle_result/mac/Pseudotime/Pseudotime.monocle.rds")

library(dplyr)
library(ggplot2)

coord <- t(reducedDimS(cds))
coord <- as.data.frame(coord)
colnames(coord) <- c("Component_1", "Component_2")
coord$cell <- rownames(coord)

meta <- pData(cds) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell")

plot_df <- left_join(coord, meta, by = "cell")

head(plot_df)

outlier_df <- plot_df %>%
  filter(Component_2 < -10)

table(outlier_df$cell_type)
table(outlier_df$group)
# > table(outlier_df$cell_type)

# Kuffer cells     Mac_Ccr2     Mac_Fth1   Mac_Hspa1a    Mac_Mmp13     Mac_Saa3 
#            0            2           12            1            4            0 
# > table(outlier_df$group)

#     PD1 PD1_Pen 
#      16       3 

outlier_cells <- plot_df %>%
  filter(Component_2 < -10) %>%
  pull(cell)

sce = readRDS("/mnt/NexData_file/home/wmx/Project/SC2605263939-01/shouhou3/Data/mac_fixed.rds")
sce_clean <- sce[, !colnames(sce) %in% outlier_cells]

saveRDS(sce_clean,"/mnt/NexData_file/home/wmx/Project/SC2605263939-01/shouhou3/Data/mac_fixed_clean.rds")
