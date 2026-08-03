# Find all markers using FindConservedMarkers
FindAllConservedMarkers <-
  function(obj, 
           ident.1 = NULL, 
           grouping.var,
           ...) {
    
  stopifnot(is(obj, "Seurat"))
  
  all_ident <- unique(Idents(obj))
  stopifnot(length(all_ident) > 1)
  
  
  if(is.null(ident.1)){
    message("ident.1 is null, will find all markers of each ident in all idents.")
  }else{
    all_ident <-  ident.1
  }
  
  all_ident %>%
    set_names(., .) %>%
    map(function(ident){
      FindConservedMarkers(
        object = obj,
        ident.1 = ident,
        grouping.var = grouping.var,
        ...
      )
    })
}




# function to normalize

sc_normalize <- function(sc, vars.to.regress = 'nCount_RNA', method = c("logNormalize", "SCT"), scale.factor = c('default', 'median')){
  
  method <- match.arg(method)
  scale.factor <- match.arg(scale.factor)
  
  if(scale.factor == 'default'){
    scale.factor <- 10000
  }else{
    scale.factor.vars <- vars.to.regress %>% str_subset("nCount")
    if(length(scale.factor.vars) < 1){
      scale.factor.vars <- "nCount_RNA"
    }else{
      scale.factor.vars <- scale.factor.vars[1]
    }
    scale.factor <- median(unlist(sc[[scale.factor.vars]]))
  }
  
  # normalize 
  if(method == "logNormalize"){
    sc <- NormalizeData(sc, scale.factor = scale.factor)
    sc <- FindVariableFeatures(sc)
    sc <- ScaleData(sc, vars.to.regress = vars.to.regress)
  }else{
    sc <- SCTransform(sc, vars.to.regress = vars.to.regress)
  }
  
  
  sc
}


# function to  pca
sc_pca <- function(
    sc, 
    ndims = 30,
    sc_out_dir = 'sc.1 pca', 
    dir.create.force = FALSE
  ){
  
  # dimension reduction
  sc <- RunPCA(sc)
  
  if(is.null(sc_out_dir)){
    ElbowPlot(sc, ndims = ndims)
    return(sc)
  }else if(!file.exists(sc_out_dir)){
    if(dir.create.force) {
      dir.create(sc_out_dir, recursive = TRUE)
    }else{
      stop(paste0(sc_out_dir, ' not exist.'))
    }
  }
  
  # dim heatmap
  pdf(file.path(sc_out_dir, 'pca plot.pdf'), height = 6, width = 8)
  print(ElbowPlot(sc, ndims = ndims) )
  dev.off()
  
  pdf(file.path(sc_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
  DimHeatmap(sc, dims = 1:5)
  dev.off()
  pdf(file.path(sc_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
  DimHeatmap(sc, dims = 6:10)
  dev.off()
  pdf(file.path(sc_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
  DimHeatmap(sc, dims = 11:15)
  dev.off()
  pdf(file.path(sc_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
  DimHeatmap(sc, dims = 16:20)
  dev.off()
  
  sc
}

# dr 
sc_dr <- function(
    sc,
    ndims = 30,
    run_harmony = TRUE,
    group.by = "project",
    return.model = TRUE,
    harmony.reduction.save = "harmony",
    harmony.reduction.name = "harmonyumap", 
    harmony.reduction.key  = "harmonyUMAP_",
    sc_out_dir = 'sc.2 dr', 
    dir.create.force = FALSE
  ){
  
  # create directory
  if(is.null(sc_out_dir)){
    message(paste0(sc_out_dir, 'is null.'))
    return(sc)
  }else if(!file.exists(sc_out_dir)){
    if(dir.create.force) {
      dir.create(sc_out_dir, recursive = TRUE)
    }else{
      stop(paste0(sc_out_dir, ' not exist.'))
    }
  }
  
  
  # umap
  sc <- RunUMAP(sc, dims = 1:ndims)
  
  pdf(file.path(sc_out_dir, str_glue("Dimplot dim{ndims} {group.by}.pdf")), height = 6, width = 8)
  p <- DimPlot(sc, group.by = group.by, label = TRUE) + NoLegend()
  print(p)
  dev.off()
  
  # harmony
  if(run_harmony){
    sc <- sc %>% 
      harmony::RunHarmony(
        dims.use = 1:ndims,
        group.by.vars = group.by,
        reduction.save = harmony.reduction.save
      )
    
    sc <- RunUMAP(
      sc, 
      dims = 1:ndims, 
      reduction = harmony.reduction.save, 
      reduction.name = harmony.reduction.name, 
      reduction.key  = harmony.reduction.key,
      return.model   = return.model
    )
    
    
    pdf(file.path(sc_out_dir, str_glue("Dimplot dim{ndims} {group.by} - {harmony.reduction.save}.pdf")), height = 6, width = 8)
    p <- DimPlot(sc, group.by = group.by, label = TRUE, reduction = harmony.reduction.name) + NoLegend()
    print(p)
    dev.off()
  }

  
  sc
}


# cluster
sc_cluster <- function(
    sc,
    ndims = 30,
    k.param = 15,
    reduction =  "harmony",
    harmony.reduction.name = "harmonyumap",
    sc_out_dir = 'sc.3 cluster', 
    dir.create.force = FALSE,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  ){
  
  
  # create directory
  if(is.null(sc_out_dir)){
    message(paste0(sc_out_dir, 'is null.'))
    return(sc)
  }else if(!file.exists(sc_out_dir)){
    if(dir.create.force) {
      dir.create(sc_out_dir, recursive = TRUE)
    }else{
      stop(paste0(sc_out_dir, ' not exist.'))
    }
  }
  
  # cluster
  sc <-
    FindNeighbors(
      sc,
      reduction = reduction,
      dims = 1:ndims,
      k.param = k.param
    )
  sc <-
    FindClusters(
      sc,
      resolution = resolution
    )
  
  assay_name <- DefaultAssay(sc)
  group_var <- str_glue("{assay_name}_snn_res.{resolution}")
  
  
  # dimplot
  pdf(file.path(sc_out_dir, str_glue("Dimplot dim{ndims} cluster - harmony.pdf")), height = 6, width = 8)
    walk(group_var, function(group.by){
        p <- DimPlot(sc, group.by = group.by, label = TRUE, reduction = "umap") + NoLegend()
        print(p)
    })
    # harmony
    if(reduction == "harmony"){
      walk(group_var, function(group.by){
        p <- DimPlot(sc, group.by = group.by, label = TRUE, reduction = harmony.reduction.name) + NoLegend()
        print(p)
      })
    }
  dev.off()
  
  sc
}


# sc_dot
sc_dotplot <- function(
    sc,
    ndims = 30,
    marker = NULL,
    sc_out_dir = 'sc.4 dotplot', 
    resolution = c(0.2, 0.4, 0.6, 0.8),
    dir.create.force = FALSE,
    width = 25, 
    height = 6
  ){
  if(is.null(marker)){
    return(sc)
  }
  
  # create directory
  if(is.null(sc_out_dir)){
    message(paste0(sc_out_dir, 'is null.'))
    return(sc)
  }else if(!file.exists(sc_out_dir)){
    if(dir.create.force) {
      dir.create(sc_out_dir, recursive = TRUE)
    }else{
      stop(paste0(sc_out_dir, ' not exist.'))
    }
  }
  
  # variable 
  assay_name <- DefaultAssay(sc)
  group_var <- str_glue("{assay_name}_snn_res.{resolution}")
  
  
  # export
  pdf(file.path(sc_out_dir, str_glue("dotplot dims{ndims} marker.pdf")), width = width, height = height)
    walk(group_var, function(group.by){
      p <- DotPlot(sc, group.by = group.by, features = marker) + RotatedAxis()
      print(p)
    })
  dev.off()
  
  sc
}

# sc_findMarkes
sc_findmarker <- function(
    sc,
    idents = "RNA_snn_res.0.2",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE,
    
    all_markers = NULL,
    p_val_adj = 0.05,
    pct.1 = 0.25,
    sort_by = c("avg_log2FC", "pct.1"),
    n_head = 10,
    
    method = c("marker", "conserved"),
    grouping.var = 'project'
){
  
  sort_by <- match.arg(sort_by)
  method  <- match.arg(method)
  
  if(is.null(all_markers)){
      Idents(sc) <- idents
      sc <- tryCatch(JoinLayers(sc), error = function(e) sc)
      
      if(DefaultAssay(sc) == "SCT") sc <- PrepSCTFindMarkers(sc)
      
      if(method == "marker"){
        all_markers <-
          FindAllMarkers(
            sc,
            test.use = test.use,
            logfc.threshold = logfc.threshold,
            min.pct = min.pct,
            only.pos = only.pos
          )
      }else{
        
        all_markers <-
          FindAllConservedMarkers(
            sc,
            grouping.var = grouping.var,
            test.use = test.use,
            logfc.threshold = logfc.threshold,
            assay = DefaultAssay(sc),
            min.pct = min.pct,
            only.pos = only.pos
          )
        all_markers <- 
          all_markers %>%
          map(rownames_to_column, "gene") %>%
          enframe("cluster") %>%
          unnest(value) %>%
          rowwise() %>%
          mutate(
            avg_log2FC = mean(c_across(ends_with("avg_log2FC")), na.rm = TRUE),
            pct.1      = mean(c_across(ends_with("pct.1")), na.rm = TRUE),
            p_val_adj  = minimump_p_val           
          )
          
      }
      
  }


    
    # print
    all_markers %>%
      group_by(cluster) %>%
      dplyr::filter(p_val_adj < {{p_val_adj}}) %>%
      arrange(desc(.data[[sort_by]]), .by_group = TRUE) %>%
      slice_head(n = n_head) %>%
      nest() %>%
      deframe() %>%
      map(pull, 'gene') %>%
      imap(~{
        xx <- str_c(.x, collapse = ",")
        str_glue("cluster_{.y}:{xx}")
      }) %>%
      walk(cat, "\n")
    
    invisible(all_markers)
}

  
sc_mark_heatmap <- function(
    sc,
    markers,
    idents = "RNA_snn_res.0.2",
    p_val_adj = 0.05,
    pct.1 = 0.25,
    n_head = 10,
    sort_by = c("avg_log2FC", "pct.1"),
    scale = 'column',
    log_trans = FALSE,
    row_cluster = FALSE
  ){ 
  
  # browser()
  
    Idents(sc) <- idents
    sc <- tryCatch(JoinLayers(sc), error = function(e) sc)
    avg_expr <- AverageExpression(sc)[['RNA']]
    
    sort_by <- match.arg(sort_by)
    
    
    de_genes <-
      markers %>% 
      group_by(cluster) %>%
      dplyr::filter(p_val_adj < {{p_val_adj}}) %>%
      dplyr::filter(pct.1 > {{pct.1}}) %>%
      arrange(desc(.data[[sort_by]]), .by_group = TRUE) %>%
      slice_head(n = n_head) %>%
      pull("gene")
    
    if(length(de_genes) < 1) {
      stop("No de genes found.")
    }
    
    if(row_cluster){
      row_cluster <- avg_expr[VariableFeatures(sc), ] %>% t %>% dist %>% hclust
    }
    if(log_trans){
      avg_expr <- log1p(avg_expr)
    }
    p <-
      avg_expr[de_genes, ] %>% 
      as.matrix() %>% 
      as.data.frame() %>% 
      # setNames(., paste('cluster', (1:ncol(.) - 1) )) %>% 
      t %>%
      pheatmap::pheatmap(
        scale = scale,
        cluster_row = row_cluster,
        cluster_col = FALSE,
        color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
        border_col = NA,
        silent = TRUE
      )
    cowplot::plot_grid(p$gtable)
}

      

  

