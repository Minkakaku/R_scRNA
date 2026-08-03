# CN analysis ----------------------------------
# sc_filter <- sc_filter %>% AddMetaData(cell_annotation_df_tidy)


CN_compute_cellfreq <- function(
    dat,
    x = "X",
    y = "Y",
    classification = "Classification1",
    return.vars = c("freq", "count"),
    k.nn = 10
){
  # browser()
  
  return.vars <- match.arg(return.vars)
  
  # x, y
  if(any(!c(x, y) %in% colnames(dat))){
    warning(str_glue("{x} or {y} not found."))
    others_x_y <- c('XMin', 'XMax',  'YMin',  'YMax') %in% colnames(dat)

    if(all(others_x_y)){
      warning("used 'XMin', 'XMax',  'YMin',  'YMax' instead.")
      dat[[x]] <- 1/2 * (as.numeric(dat[['XMin']]) + as.numeric(dat[['XMax']]))
      dat[[y]] <- 1/2 * (as.numeric(dat[['YMin']]) + as.numeric(dat[['YMax']]))
    }else{
      stop(str_glue("{x} or {y} should exist."))
    }
  }

  stopifnot(is.numeric(dat[[x]]))
  stopifnot(is.numeric(dat[[y]]))
  stopifnot(classification %in% colnames(dat))

  # NN
  cn_nn <- RANN::nn2(dat[c(x, y)], k = k.nn)
  cn_nn_index2cell <- as.character(dat[[classification]])[cn_nn$nn.idx]

  dim(cn_nn_index2cell)      <- dim(cn_nn$nn.idx)
  rownames(cn_nn_index2cell) <- dat$obj_id

  # get cell frequency table
  # add all celltypes, make sure all cell exist in every neighborhood
  #all_celltypes <- unique(dat[[classification]])
  all_celltypes <- levels(unlist(dat[[classification]]))
  cn_all_celltypes <- matrix(
    all_celltypes,
    nrow = nrow(cn_nn_index2cell),
    ncol = length(all_celltypes),
    byrow = TRUE
  )
  cn_nn_cell2 <- cbind(cn_nn_index2cell, cn_all_celltypes)

  # all cell have been added by 1, so minus 1
  # table will sort the celltype, so all celltype have identical order
  cn_nn_cellcount <- apply(cn_nn_cell2, 1, table) - 1
  cn_nn_cellfreq  <- sweep(cn_nn_cellcount, 2, FUN = "/", colSums(cn_nn_cellcount))

  
  if(return.vars == 'count'){
    return(cn_nn_cellcount)
  }
  
  cn_nn_cellfreq
}


CN_compute_CN <- function(
    dat,
    x = "X",
    y = "Y",
    cache = TRUE, # TODO：缓存机制
    group.by = "Slide_ID",
    classification = "Classification1",
    k.nn = 10,
    k.cluster = 10,
    seed=1234,
    iter.max = 1000
){
  set.seed(seed)
  stopifnot(group.by %in% colnames(dat))

  # browser()
  cn_nn_cellfreq_ls <-
    dat %>%
    group_by(.data[[group.by]]) %>%
    nest() %>%
    deframe() %>%
    map(function(sub_dat){

      CN_compute_cellfreq(
        sub_dat,
        x = x,
        y = y,
        classification = classification,
        k.nn = k.nn
      )
    })

  # cn_nn_cellfreq
  # row: CellType; column: cell_id
  cn_nn_cellfreq <- do.call(cbind, cn_nn_cellfreq_ls)

  # construct CN by all slide

  # cn_cluster_ls <- map(2:20, ~ kmeans(t(cn_nn_cellfreq), centers = .x, iter.max = 1000))
  # cn_cluster_vrc <- map(cn_cluster_ls, ~ fpc::calinhara(t(cn_nn_cellfreq), .x$cluster)) # compute 方差比

  # plot(
  #   x = 2:20,
  #   y = unlist(cn_cluster_vrc),
  #   ylab = "Calinski-Harabasz index",
  #   xlab = "Number of Clusters" k.cluster
  # )
  cn_cluster     <- kmeans(t(cn_nn_cellfreq), centers = k.cluster,iter.max = iter.max)

  cn_cluster_tidy <-
    cn_cluster$cluster %>%
    enframe("obj_id", "CN") %>%
    column_to_rownames("obj_id") %>%
    mutate(CN = str_c("CN", CN)) %>%
    mutate(CN = factor(CN, levels = str_c("CN", sort(unique(cn_cluster$cluster)))))

  # baseline # expressiong by CN
  cluster_centers <- cn_cluster$centers
  base_mean       <- rowMeans(cn_nn_cellfreq)

  res <-
    list(
      cluster_centers = cluster_centers,
      base_mean = base_mean,
      cn_nn_cellfreq = cn_nn_cellfreq,
      cn_cluster = cn_cluster_tidy
    )
  res
}

# # implement CN
# options(warn = 0)
# classfication <- "Classification1"
#
# res_cn <- CN_compute_CN(
#   sc_filter[[]],
#   x = "X",
#   y = "Y",
#   group.by = "Slide_ID",
#   classification = classfication,
#   k.nn = 10,
#   k.cluster = 10
# )

CN_compute_CN2 <- function(
    dat,
    query = NULL,
    x = "X",
    y = "Y",
    cache = TRUE, # TODO：缓存机制
    group.by = "Slide_ID",
    classification = "Classification1",
    k.nn = 10,
    k.cluster = 10,
    seed=1234,
    iter.max = 1000
){
  set.seed(seed)
  stopifnot(group.by %in% colnames(dat))
  
  # browser()
  cn_nn_cellfreq_ls <-
    dat %>%
    group_by(.data[[group.by]]) %>%
    nest() %>%
    deframe() %>%
    map(function(sub_dat){
      
      CN_compute_cellfreq(
        sub_dat,
        x = x,
        y = y,
        classification = classification,
        k.nn = k.nn
      )
    })
  
  # cn_nn_cellfreq
  # row: CellType; column: cell_id
  cn_nn_cellfreq <- do.call(cbind, cn_nn_cellfreq_ls)
  
  if(!is.null(query)){
    stopifnot( all(query %in% dat[[classification]]) )
    
    query_idx <- dat[[classification]] %in% query
    stopifnot(sum(query_idx) > 1)
    cn_nn_cellfreq <- cn_nn_cellfreq[, query_idx]
  }
  
  
  # construct CN by all slide
  
  # cn_cluster_ls <- map(2:20, ~ kmeans(t(cn_nn_cellfreq), centers = .x, iter.max = 1000))
  # cn_cluster_vrc <- map(cn_cluster_ls, ~ fpc::calinhara(t(cn_nn_cellfreq), .x$cluster)) # compute 方差比
  
  # plot(
  #   x = 2:20,
  #   y = unlist(cn_cluster_vrc),
  #   ylab = "Calinski-Harabasz index",
  #   xlab = "Number of Clusters" k.cluster
  # )
  cn_cluster     <- kmeans(t(cn_nn_cellfreq), centers = k.cluster,iter.max = iter.max)
  
  cn_cluster_tidy <-
    cn_cluster$cluster %>%
    enframe("obj_id", "CN") %>%
    column_to_rownames("obj_id") %>%
    mutate(CN = str_c("CN", CN)) %>%
    mutate(CN = factor(CN, levels = str_c("CN", sort(unique(cn_cluster$cluster)))))
  
  # baseline # expressiong by CN
  cluster_centers <- cn_cluster$centers
  base_mean       <- rowMeans(cn_nn_cellfreq)
  
  res <-
    list(
      cluster_centers = cluster_centers,
      base_mean = base_mean,
      cn_nn_cellfreq = cn_nn_cellfreq,
      cn_cluster = cn_cluster_tidy
    )
  res
}

# Proximuty analysis -----------------
ProximutyAnalysis_NN <- function(
    data,
    query = data,
    x = "X",
    y = "Y",
    classification = 'Classification1',
    celltypes = NULL,
    k = 'auto',
    searchtype = "radius",
    radius = 100, # micron
    mpp = 0.325, # micronPerPixel: 每像素的大小
    returntype = c("count", "freq"),
    cell_enough_ratio_cutoff = 90 # 保证90%以上的细胞都已经计算了全部的临近细胞
){
  # x, y
  if(any(!c(x, y) %in% colnames(data))){
    warning(str_glue("{x} or {y} not found."))
    others_x_y <- c('XMin', 'XMax',  'YMin',  'YMax') %in% colnames(data)

    if(all(others_x_y)){
      warning("used 'XMin', 'XMax',  'YMin',  'YMax' instead.")
      data[[x]] <- 1/2 * (as.numeric(data[['XMin']]) + as.numeric(data[['XMax']]))
      data[[y]] <- 1/2 * (as.numeric(data[['YMin']]) + as.numeric(data[['YMax']]))
    }else{
      stop(str_glue("{x} or {y} should exist."))
    }
  }

  # x, y
  if(any(!c(x, y) %in% colnames(query))){
    warning(str_glue("{x} or {y} not found."))
    others_x_y <- c('XMin', 'XMax',  'YMin',  'YMax') %in% colnames(query)

    if(all(others_x_y)){
      warning("used 'XMin', 'XMax',  'YMin',  'YMax' instead.")
      query[[x]] <- 1/2 * (as.numeric(query[['XMin']]) + as.numeric(query[['XMax']]))
      query[[y]] <- 1/2 * (as.numeric(query[['YMin']]) + as.numeric(query[['YMax']]))
    }else{
      stop(str_glue("{x} or {y} should exist."))
    }
  }
  stopifnot(is.numeric(data[[x]]))
  stopifnot(is.numeric(data[[y]]))
  stopifnot(is.numeric(query[[x]]))
  stopifnot(is.numeric(query[[y]]))
  stopifnot(is.numeric(radius) && radius >= 0)
  
  returntype <- match.arg(returntype)
  
  # 10 microns/cell
  if(k == "auto"){
    k <- (radius %/% 10 + 1)^2 * 2
  }else if(k == Inf){
    k <- nrow(data) %/% 2
  }
  stopifnot(is.numeric(k))
  if(k < 10) k <- 10

  # find NN
  res <-
    RANN::nn2(
      data[c(x, y)],
      query = query[c(x, y)],
      k = k,
      searchtype = searchtype,
      radius = radius/mpp
    )

  # make sure k is enough
  cell_enough <- sum(rowSums(res$nn.idx == 0) > 0)
  cell_enough_ratio <- cell_enough/nrow(res$nn.idx) * 100
  
  loop_count <- 1
  while(cell_enough_ratio < cell_enough_ratio_cutoff){
    
    k <- k + 100 * loop_count
    
    res <-
      RANN::nn2(
        data[c(x, y)],
        query = query[c(x, y)],
        k = k,
        searchtype = searchtype,
        radius = radius/mpp
      )
    
    # make sure k is enough
    cell_enough <- sum(rowSums(res$nn.idx == 0) > 0)
    cell_enough_ratio <- cell_enough/nrow(res$nn.idx) * 100
    
    # 最多加 100 * sum(1:10) = 5500个细胞就会跳出循环
    loop_count <- loop_count + 1
    if(loop_count > 10) break
  }
  
  # omit cell which is not a neighbor
  all_zeros_idx <- colSums(res$nn.idx != 0) == 0
  res$nn.idx    <- res$nn.idx[, !all_zeros_idx]
  res$nn.dists  <- res$nn.dists[, !all_zeros_idx]

  message(str_glue("k: {k} is enough for {signif(cell_enough_ratio, 1)}% cells"))

  # no classification, return
  if(is.null(classification)){
    res$query_cell_id <- if(is.null(query$cell_id)) rownames(query) else query$cell_id
    res$total_cell_id <- if(is.null(data$cell_id))  rownames(data)  else data$cell_id
    return(res)
  }

  # dat tidy
  # If there are no neighbours, nn.idx is zero
  stopifnot(is.character(classification))

  if(is.null(celltypes)){
    celltypes <- unique(data[[classification]])
  }
  
  if(is.factor(celltypes)){
    celltypes <- as.character(celltypes)
  }


  index2cell <- as.character(data[[classification]])[t(res$nn.idx)]
  index2index <- (res$nn.idx != 0) * (1:nrow(res$nn.idx))
  index2index <- c(t(index2index)) %>% .[ . != 0]

  nn_cellcount  <-
    tapply(
      index2cell,
      INDEX = index2index,
      FUN = function(dat) unclass(table(c(dat, celltypes)) - 1),
      simplify = TRUE
    )
  nn_cellcount <- t(do.call(rbind, unname(nn_cellcount)))
  if(!is.null(query$obj_id)){
    colnames(nn_cellcount) <- query$obj_id
  }
  
  # return count
  if(returntype == "count") return(nn_cellcount)
  
  # return freqency
  nn_cellfreq  <- sweep(nn_cellcount, 2, FUN = "/", colSums(nn_cellcount))
  nn_cellfreq
}


ProximutyAnalysis_NN_multi <- function(
    data,
    query = data,
    x = "X",
    y = "Y",
    group.by = c("Slide_ID", "ID"),
    classification = 'Classification1',
    k = 'auto',
    searchtype = "radius",
    radius = 100, # micron
    mpp = 0.325, # micronPerPixel: 每像素的大小
    returntype = c("count", "freq"),
    cell_enough_ratio_cutoff = 90 # 保证90%以上的细胞都已经计算了全部的临近细胞
    
){

  group.by <- match.arg(group.by)
  stopifnot(group.by %in% colnames(data))
  
  returntype <- match.arg(returntype)
  
  # browser()

  if(!is.null(classification)){
    celltypes <- unique(data[[classification]])
  }
  all_grps  <- unique(data[[group.by]])

  cellfreq_ls <-
    all_grps %>%
    set_names(.,.) %>%
    map(function(grp){
      dat <- data  %>%  dplyr::filter(.data[[group.by]] == grp)
      quy <- query %>%  dplyr::filter(.data[[group.by]] == grp)

      ProximutyAnalysis_NN(
        data = dat,
        query = quy,
        x = x,
        y = y,
        classification = classification,
        celltypes = celltypes,
        k = k,
        searchtype = searchtype,
        radius = radius,
        mpp = mpp,
        returntype = returntype,
        cell_enough_ratio_cutoff = cell_enough_ratio_cutoff
      )
    })

  if(is.null(classification)) return(cellfreq_ls)

  # cn_nn_cellfreq
  # row: CellType; column: cell_id
  cellfreq <- do.call(cbind, cellfreq_ls)

  cellfreq
}


# undebug(ProximutyAnalysis_NN)
# undebug(ProximutyAnalysis_NN_multi)
#
# res_pa <-
#   ProximutyAnalysis_NN_multi(
#     sc[[]],
#     sc[[]],
#     x = "X",
#     y = "Y",
#     classification = "Classification1",
#     k = Inf, # 'auto',
#     searchtype = "radius",
#     radius = 100,
#     mpp = 0.325
#   )
