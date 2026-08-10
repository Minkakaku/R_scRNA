
dir.create("result_0703")
B_cell_out_dir_2 <- "result_0703/B_cell_subtype" 
dir.create(B_cell_out_dir_2)

B_cell <- readRDS("result_0622/B_cell_subtype/B_cell_cluster_1.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "Classical_memory_B", # 
  '1'  = "Naive_B", #
  '2'  = "Plasmablast", # 
  '3'  = "Immature_B", # 
  '4'  = "GC_centroblast", # 
  '5'  = "Immature_B",#  
  '6'  = "Immature_B", #  
  '7'  = "Immature_B", #  
  '8'  = "Atypical_memory_B", #  
  '9'  = "FDC", #  
  '10'  = "GC_centrocyte", #  
  '11'  = "Immature_B"
)
Idents(B_cell) <- "RNA_snn_res.0.6"
B_cell <- RenameIdents(B_cell, cell_cluster)
B_cell$subtype <- Idents(B_cell)



pdf(file.path(B_cell_out_dir_2, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(B_cell, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


keep <- c(
  "Plasmablast",
  "Naive_B",
  "Classical_memory_B",
  "Atypical_memory_B",
  "GC_centroblast",
  "GC_centrocyte",
  "FDC"
)

B_cell_clean <- subset(B_cell, subtype %in% keep)

B_cell_clean$subtype <- droplevels(B_cell_clean$subtype)
Idents(B_cell_clean) <- B_cell_clean$subtype

table(B_cell_clean$subtype)
saveRDS(B_cell_clean, "result_0703/B_cell_subtype/B_cell_clean.rds")

B_cell <- B_cell_clean

pdf(file.path(B_cell_out_dir_2, str_glue("Dimplot dim{dims} B_cell annotation_harmony.pdf")), height = 6, width = 8)
DimPlot(B_cell, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()


B_cell <- RunUMAP(B_cell, dims = 1:30,n.neighbors = 50,return.model = TRUE)

pdf(file.path(B_cell_out_dir_2, str_glue("Dimplot dim30 B_cell annotation.pdf")), height = 6, width = 8)
DimPlot(B_cell, group.by = "subtype", label = TRUE, reduction = 'umap',repel = TRUE) + NoLegend()
dev.off()

B_cell_major_marker_ls_2 <- c(
  # 未成熟B Immature B
  "B_cell" = "CD19",
  #"B_cell" = "CD79A",
  "B_cell" = "MS4A1",
  # 初始B naive B
  "Naive_B" = "BANK1",
  "Naive_B" = "FCRL1",
  "Naive_B" = "FCER2",
  "Naive_B" = "IL4R", 
  # 记忆B memory B
  "Classical_Memory_B"= "AIM2",
  "Classical_Memory_B"= "TNFRSF13B",
  #"Classical_Memory_B"= "CD22",
  "Classical_Memory_B"= "FCRL2",
  "Classical_Memory_B"= "CD27",
  #"memory B"= "ADAM28",
  #"memory B"= "BLK", 
  "Atypical_memory_B"= "ITGAX",
  "Atypical_memory_B"= "FCRL5",
  "Atypical_memory_B"= "BATF",
  
  "GC_centrocyte"="TCL1A",
  "GC_centrocyte"= "MME",
  "GC_centrocyte"="AICDA",
  "GC_centroblast"="MKI67",
  "GC_centroblast"="STMN1",
  
  "Plasmablast" = "IRF4",
  "Plasmablast" = "PRDM1",
  "Plasmablast" = "CD38",
  
  "FDC" = "CR1",
  "FDC" = "CR2"
)

pdf(file.path(B_cell_out_dir_2, str_glue("dotplot B_cell subtype major_marker.pdf")), width = 12, height = 6)
DotPlot(B_cell, group.by = "subtype", features = B_cell_major_marker_ls_2)+RotatedAxis()+theme_90angle
dev.off()




## cell props barplot  ----------------------------------------------------------
B_cell_meta <- B_cell[[]]

B_cell_meta$sample <- factor(B_cell_meta$sample,levels=c(B_cell_meta$sample %>% unique()))

p_col_proj <-
  B_cell_meta %>%
  group_by(sample, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = subtype)) + 
  geom_col() +
  scale_fill_brewer(palette = "Set2", name = "B_cell Subtype") +
  #ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(B_cell_out_dir_2, "cell_freq_sample.pdf"),
  height = 6,
  width  = 8
)

p_col_proj <-
  B_cell_meta %>%
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
  filename = file.path(B_cell_out_dir_2, "cell_freq_group.pdf"),
  height = 6,
  width  = 4
)

B_cell$subtype %>% table()
# 映射到空间切片
subtype_col <- 'subtype'

cells_position_ls <-
  cbind(
    B_cell[['sample']],
    Embeddings(B_cell[['spatial']])
  ) %>%
  rownames_to_column('cell_id') %>%
  group_by(sample) %>%
  nest() %>%
  deframe()


plot_dat <-
  FetchData(
    B_cell, 
    vars = c(subtype_col)
  ) %>%
  rownames_to_column('cell_id')



# plot celltype ------------------
all_cells <-  B_cell$subtype %>% unique %>% as.character()
cell_colors <- structure(
  names = all_cells,
  ggsci:: pal_igv()(length(all_cells))
)



p_sp_ct_ls <-
  cells_position_ls %>% # .[1] %>%
  imap(function(coords, sam){
    dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
    
    height <- 5
    width  <- 5
    p <-
      dat2 %>%
      ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
      geom_tile(
        aes(fill = .data[[subtype_col]]),
        height = height, 
        width = width,
        show.legend = c(fill = TRUE, alpha = FALSE)
      ) +
      theme_void() +
      coord_equal() +
      scale_y_reverse() +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values = cell_colors) +
      theme(panel.background = element_rect(fill = "black", color = NA)) +
      ggtitle(sam)
    p
    
  })


pdf(
  file.path("result_0703/B_cell_subtype/spatial plot of celltype in each sample.pdf"),
  height = 10,
  width  = 10
)
p_sp_ct_ls %>% walk(print)
dev.off()


# plot each cell in all sample -----------
p_sp_ct_ls_ls <-
  cells_position_ls %>% # .[1] %>%
  imap(function(coords, sam){
    dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
    
    height <- 5
    width  <- 5
    
    p_ls <-
      all_cells %>%
      setNames(., .) %>%
      map(function(x){
        cell_colors[names(cell_colors) != x] <- NA
        
        p <-
          dat2 %>%
          ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
          geom_tile(
            aes(fill = .data[[subtype_col]]),
            height = height, 
            width = width,
            show.legend = c(fill = TRUE, alpha = FALSE)
          ) +
          theme_void() +
          coord_equal() +
          scale_y_reverse() +
          # ggsci::scale_fill_igv() +
          scale_fill_manual(values = cell_colors) +
          theme(panel.background = element_rect(fill = "black", color = NA)) +
          ggtitle(sam)
        p
      })
    p_ls
  })

# export by sample
p_sp_ct_ls_ls %>%
  iwalk(function(p_ls, sample){
    pdf(
      file.path(str_glue("result_0703/B_cell_subtype/spatial plot of celltype in each sample - {sample}.pdf")),
      height = 10,
      width  = 10
    )
    p_ls %>% walk(print)
    dev.off()
  })

# rds文件保存
B_cell_1 <- readRDS("result_0622/B_cell_subtype/B_cell_cluster.rds")
cell_cluster <- c(
  '0'  = "Naive_B", # 
  '1'  = "T_cell-like_B_cell", #
  '2'  = "Memory_B", # 
  '3'  = "Immature_B", # 
  '4'  = "Plasmablast", # 
  '5'  = "Immature_B",#  
  '6'  = "GC_B", #  
  '7'  = "T_cell-like_B_cell", #  
  '8'  = "Immature_B", #  
  '9'  = "Immature_B", #  
  '10'  = "T_cell-like_B_cell", #  
  '11'  = "Immature_B", # 
  '12'  = "Naive_B", # 
  '13'  = "GC_B", #  
  '14'  = "Plasmablast", # 
  '15'  = "Plasmablast"
)

Idents(B_cell_1) <- "RNA_snn_res.0.8"
B_cell_1 <- RenameIdents(B_cell_1, cell_cluster)
B_cell_1$subtype <- Idents(B_cell_1)


B_cell <- readRDS("result_0622/B_cell_subtype/B_cell_cluster_1.rds")
## annotation ----------------------------------------------------------------
cell_cluster <- c(
  '0'  = "Classical_memory_B", # 
  '1'  = "Naive_B", #
  '2'  = "Plasmablast", # 
  '3'  = "Immature_B", # 
  '4'  = "GC_centroblast", # 
  '5'  = "Immature_B",#  
  '6'  = "Immature_B", #  
  '7'  = "Immature_B", #  
  '8'  = "Atypical_memory_B", #  
  '9'  = "FDC", #  
  '10'  = "GC_centrocyte", #  
  '11'  = "Immature_B"
)
Idents(B_cell) <- "RNA_snn_res.0.6"
B_cell <- RenameIdents(B_cell, cell_cluster)
B_cell$subtype <- Idents(B_cell)

B_cell_1$subtype2 <- as.character(B_cell_1$subtype)
cell_subtype <- B_cell$subtype
names(cell_subtype) <- colnames(B_cell)
B_cell_1$subtype2[colnames(B_cell)] <- as.character(cell_subtype)
table(B_cell_1$subtype2)
table(B_cell_1$subtype)
table(B_cell$subtype)
saveRDS(B_cell_1, "result_0703/B_cell_subtype/B_cell_allannotation.rds")

#合并到all_data_T_DC_annotation.rds中
all_data$annotation3 <- as.character(all_data$annotation2)
cell_subtype <- B_cell_1$subtype2
names(cell_subtype) <- colnames(B_cell_1)
all_data$annotation3[colnames(B_cell_1)] <- as.character(cell_subtype)
table(all_data$annotation3)
table(all_data$annotation2)
table(B_cell_1$subtype2)
saveRDS(all_data, "result_0703/B_cell_subtype/all_data_T_DC_B_annotation.rds")
