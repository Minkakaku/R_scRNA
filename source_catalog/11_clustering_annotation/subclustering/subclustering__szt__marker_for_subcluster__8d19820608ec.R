major_marker_for_subcluster <- c(
  
  Endothelial = "PECAM1",
  Endothelial = "VWF",
  Endothelial = "AQP1",
  
  Epithelial  = "EPCAM", 
  Epithelial  = "KRT5",
  
  "Fibroblasts" = "SFRP4", 
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "COL1A1",
  "Fibroblasts" = "FAP",
  "Fibroblasts" = "RGS5",
  
  "Smooth muscle cells" = "TAGLN", 
  "Smooth muscle cells" = "MYH11", 
  "Smooth muscle cells" = "ACTA2",
  
  
  Immune = 'PTPRC',
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  "Plasma" = "TNFRSF17",
  "Plasma" = "JCHAIN",
  "Plasma" = "MZB1",
  
  "T cells" = "CD3E",
  "T cells" = "CD3G",
  
  "NK" = "GNLY",
  "NK" = "NKG7",
  "NK" = "NCAM1",
  
  "Myeloid" = "ITGAM",
  "Macrophage" = "CD68",
  "Macrophage" = "CD163",
  
  "Mast cell" = "CPA3",
  "Mast cell" = "CTSG",
  "Mast cell" = "TPSAB1",
  
  "Monocytes" = "CD14", 
  "Monocytes" = "FCGR3A",
  
  "Dendritic cell" = "CD83",
  "Dendritic cell" = "FCER1A",
  'Dendritic cell' = "LILRA4",  
  'Dendritic cell' = "CLEC9A",
  
  Neutrophil = "G0S2",
  Neutrophil = "CSF3R",
  Neutrophil = "S100A8",
  Neutrophil = "S100A9"
)


b_markers <- c(
  'EPCAM',
  "CD3D",
  "CD3E",
  "CCDC80",
  "PECAM1",
  "CD68",
  "CD14",
  "LILRA4",
  "CD1C",
  "IDO1",
  "XCR1",
  "CLEC9A",
  "CLEC10A",
  
  "CD19",
  "MS4A1",
  "CD79A",
  
  'CD27',
  
  
  "TNFRSF17",
  "JCHAIN",
  "MZB1",
  "XBP1"
)


mye_markers <- c(
  'EPCAM',
  "CD3D",
  "CD3E",
  "CCDC80",
  "PECAM1",
  
  "MS4A1",
  "MZB1",
  "XBP1",
  'CD27',
  
  "CD68",
  "CD14",
  "FCGR3A",
  "C1QA",
  "C1QC",
  'CSF1R', # macro
  
  'FCN1', # Mono
  
  'XCR1', # cDC1
  'IDO1', # cDC2
  'LAMP3', # mature DC
  'LILRA4', #pDC
  'CLEC9A', # cDC1
  
  'MS4A2', # mast cell
  'KIT',
  
  'CSF3R', # netutrophil
  'S100A8', # mono?
  'S100A9' # mono?
)



t_markers <- c(
  'EPCAM',
  
  "CD3D",
  "CD3E",
  "TRAC",
  "IL7R",
  "CD4",
  "CD8A",
  
  "IL2RA",
  "FOXP3",
  
  'TRDC',
  
  'EOMES',
  'PRF1',
  'NKG7',
  "NCAM1",
  "KLRF1",
  "KLRD1",
  
  'TCF7', 'SELL', 'CCR7', #Tn
  'GZMK', # Tem
  'CXCR6', # Trm
  'KLRG1', # TEMRA
  'PDCD1', 'TOX', 'TIGIT', 'CTLA4', # Tex
  'TOX2', 'CXCR5', # Tfh
  
  "FAP",
  "PDPN",
  "CCDC80",
  'PDGFRA',
  "PECAM1",
  
  "MS4A1",
  "MZB1",
  "XBP1",
  'CD27',
  
  "CD68",
  "CD14",
  "FCGR3A",
  "C1QA",
  "C1QC",
  'CSF1R', # macro
  
  'FCN1', # Mono
  
  'XCR1', # cDC1
  'IDO1', # cDC2
  'LAMP3', # mature DC
  'LILRA4', #pDC
  'CLEC9A', # cDC1
  
  'MS4A2', # mast cell
  'KIT',
  
  'CSF3R', # netutrophil
  'S100A8', # mono?
  'S100A9', # mono?
  
  "MKI67",
  "TOP2A"
)