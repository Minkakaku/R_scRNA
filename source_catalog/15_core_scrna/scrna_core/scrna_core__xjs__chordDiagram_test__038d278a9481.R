library(CellChat)
library(circlize)
library(ComplexHeatmap)


pathways <- cellchat@netP$pathways
netVisual_chord_cell(cellchat, signaling = pathways)

all_cellchat_new_ls %>% 
  iwalk(function(cellchat, name){ # browser()
    pathways <- cellchat@netP$pathways
    p <- netVisual_chord_cell(
      cellchat, 
      signaling = pathways, 
      lab.cex = 0.8, 
      title.name = 'THBS1-mCAF to Epi', 
      # title.name = paste0(pathways, " signaling network"), 
      link.visible = T, 
      sources.use = 'THBS1-mCAF', 
      targets.use = c('Epi_0', 'Epi_2', 'Epi_11')
    )
    print(p)
  })

# 指定要显示的细胞类型
sources.use <- 'THBS1-mCAF'
targets.use <- c('Epi_0', 'Epi_2', 'Epi_11')
pathway.show <- c('CXCL', 'COLLAGEN', 'LAMININ', 'THBS', 'MPZ')

all_cellchat_new_ls %>% 
  iwalk(function(cellchat, name){
    pathway.show %>% 
      map(function(pathway){ # browser()
        p <- netVisual_chord_cell(
          cellchat,
          signaling = pathway,
          sources.use = sources.use,
          targets.use = targets.use,
          lab.cex = 0.8,
          title.name = paste0(pathway, " Signaling Network")
        )
        print(p)
      })
  })

cellchat@net$prob %>% dim()
cellchat@net$prob[1, , ]
cellchat@net$prob[, 1, ]
