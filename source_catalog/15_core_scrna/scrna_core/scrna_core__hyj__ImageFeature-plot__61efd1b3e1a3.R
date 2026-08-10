


SingleImagePlot_wyt <- function (data,feature, col.by = NA, col.factor = TRUE, cols = NULL, 
                                 shuffle.cols = FALSE, size = 0.1, molecules = NULL, mols.size = 0.1, 
                                 mols.cols = NULL, mols.alpha = 1, alpha = molecules %iff% 
                                   0.3 %||% 0.6, border.color = "white", border.size = NULL, cells.use =NULL,
                                 na.value = "grey50", border.by = border.by, dark.background = TRUE, ...) 
{
  if(!is.null(data) && is.data.frame(data) && nrow(data) > 0) {
    if (!all(c("x", "y", "cell", "boundary") %in% colnames(x = data))) {
      stop("Invalid data coordinates")
    }
    
    if (!is.na(x = col.by)) {
      if (!col.by %in% colnames(x = data)) {
        warning("Cannot find 'col.by' ('", col.by, "') in data coordinates", 
                immediate. = TRUE)
        col.by <- NA
      }
      else if (isTRUE(x = col.factor) && !is.factor(x = data[[col.by]])) {
        data[[col.by]] <- factor(x = data[[col.by]], 
                                 levels = unique(x = data[[col.by]]))
      }
      else if (isFALSE(x = col.factor) && is.factor(x = data[[col.by]])) {
        data[[col.by]] <- as.vector(x = data[[col.by]])
      }
    }

    if (!is.null(cells.use)) {
      data[[border.by]] <- as.character(data[[border.by]])
      data[[border.by]][!(data[[border.by]] %in% cells.use)] <- "other_cells"
      data[[col.by]][!(data[[border.by]] %in% cells.use)] <- NA
      data[[border.by]] <- factor(data[[border.by]], 
                                  levels = c(unique(data[[border.by]])))
      border.color <- c(border.color, "other_cells" = na.value)
    }
    
    if (is.na(x = col.by) && !is.null(x = cols)) {
      col.by <- RandomName(length = 7L)
      data[[col.by]] <- TRUE
    }
    if (!is.factor(x = data$boundary)) {
      data$boundary <- factor(x = data$boundary, levels = unique(x = data$boundary))
    }
    if (is.na(x = alpha)) {
      alpha <- 1L
    }
    alpha.min <- ifelse(test = alpha < 1L, yes = 1 * (10^Seurat:::.FindE(x = alpha)), 
                        no = 0.1)
    if (alpha.min == alpha) {
      alpha.min <- 1 * (10^(Seurat:::.FindE(x = alpha) - 1))
    }
    alphas <- Seurat:::.Cut(min = alpha.min, max = alpha, n = length(x = levels(x = data$boundary)))
  }
  plot <- ggplot(data = data %NA% NULL, mapping = aes_string(x = "y", 
                                                             y = "x", alpha = "boundary", 
                                                             fill = col.by %NA% NULL)) + 
    scale_alpha_manual(values = alphas, guide = "none")
  
  if(!is.null(data) && is.data.frame(data) && nrow(data) > 0) {
    if (anyDuplicated(x = data$cell)) {
      if (is.null(border.size)) {
        border.size <- 0.1
      }
      if (!is.null(border.by) && border.by %in% colnames(data)) {
        plot <- plot + geom_polygon(mapping = aes_string(group = "cell", color = border.by), 
                                    size = border.size) +
          scale_color_manual(values = border.color, guide = guide_legend(override.aes = list(fill = NA, linewidth = 3)))
      } else {
        plot <- plot + geom_polygon(mapping = aes_string(group = "cell"), 
                                    color = border.color[1], size = border.size)
      }
      
      
      if (!is.null(border.by) && border.by %in% colnames(data)) {
        if (!is.null(cells.use)) {
          fill_colors <- ifelse(levels(data[[border.by]]) == "other_cells", na.value, NA)
          plot <- plot + geom_polygon(mapping = aes_string(group = "cell", color = border.by),
                                      size = border.size) +
            scale_color_manual(
              values = border.color,
              guide = guide_legend(
                override.aes = list(
                  fill = fill_colors,
                  linewidth = 3
                )
              )
            )
        }else{
          plot <- plot + geom_polygon(mapping = aes_string(group = "cell", color = border.by), 
                                      size = border.size) +
            scale_color_manual(values = border.color, guide = guide_legend(override.aes = list(fill = NA, linewidth = 3)))
        }
        
      } else {
        plot <- plot + geom_polygon(mapping = aes_string(group = "cell"),
                                    color = border.color[1], size = border.size)
      }
      
    }
  }
  else {
    # 点的情况 - 类似处理
    if (is.null(border.size)) {
      border.size <- 0
    }
    if (is.character(border.color) && length(border.color) == 1) {
      geom_point(shape = 21, color = border.color, 
                 stroke = border.size, size = size)
    } else if (is.character(border.color) && length(border.color) > 1) {
      # 点图的边界颜色映射
      if (!is.null(border.by) && border.by %in% colnames(data)) {
        geom_point(shape = 21, mapping = aes_string(color = border.by),
                   stroke = border.size, size = size) +
          scale_color_manual(values = border.color)
      } else {
        geom_point(shape = 21, color = border.color[1], 
                   stroke = border.size, size = size)
      }
    }
    
    if (!is.null(x = cols)) {
      plot <- plot + if (is.numeric(x = cols) || cols[1L] %in% 
                         rownames(x = brewer.pal.info)) {
        palette <- brewer.pal(n = length(x = levels(x = data[[col.by]])), 
                              cols)
        if (length(palette) < length(x = levels(x = data[[col.by]]))) {
          num.blank <- length(x = levels(x = data[[col.by]])) - 
            length(palette)
          palette <- c(palette, rep(na.value, num.blank))
        }
        if (isTRUE(shuffle.cols)) {
          palette <- sample(palette)
        }
        scale_fill_manual(values = palette, na.value = na.value)
      }
      else if (cols[1] %in% DiscretePalette(n = NULL)) {
        scale_fill_manual(values = DiscretePalette(n = length(x = levels(x = data[[col.by]])), 
                                                   palette = cols, shuffle = shuffle.cols), na.value = na.value)
      }
      else {
        if (isTRUE(shuffle.cols)) {
          cols <- sample(cols)
        }
        scale_fill_manual(values = cols, na.value = na.value)
      }
    }
    if (length(x = levels(x = data$boundary)) == 1L) {
      plot <- plot + guides(alpha = "none")
    }
    if (isTRUE(x = col.factor) && length(x = levels(x = data[[col.by]])) <= 
        1L) {
      plot <- plot + guides(fill = "none")
    }
  }
  if (is.data.frame(x = molecules)) {
    if (all(c("x", "y", "molecule") %in% colnames(x = molecules))) {
      if (!is.factor(x = molecules$molecule)) {
        molecules$molecule <- factor(x = molecules$molecule, 
                                     levels = unique(x = molecules$molecule))
      }
      plot <- plot + geom_point(mapping = aes_string(fill = NULL, 
                                                     alpha = NULL, color = "molecule"), data = molecules, 
                                size = mols.size, alpha = mols.alpha, show.legend = c(color = TRUE, 
                                                                                      fill = FALSE, alpha = FALSE)) + guides(color = guide_legend(override.aes = list(size = 3L))) + 
        scale_color_manual(name = "Molecules", values = mols.cols %||% 
                             brewer.pal(n = length(x = levels(x = molecules$molecule)), 
                                        name = "Set1"), guide = guide_legend())
    }
    else {
      warning("Invalid molecule coordinates", immediate. = TRUE)
    }
  }
  if (isTRUE(dark.background)) {
    plot <- plot + DarkTheme()
  }
  return(plot)
}


ImageFeaturePlot_wyt <- function (object, features, fov = NULL, boundaries = NULL, cols = if (isTRUE(x = blend)) {
  c("lightgrey", "#ff0000", "#00ff00")} else {c("lightgrey", "firebrick1")}, 
  size = 0.5, min.cutoff = NA, max.cutoff = NA, split.by = NULL, molecules = NULL, mols.size = 0.1, mols.cols = NULL, nmols = 20000, 
  alpha = 1, border.color = "white", border.size = NULL, dark.background = TRUE, cells.use =NULL,noncells.color = "grey",
  blend = FALSE, blend.threshold = 0.5, crop = FALSE, cells = NULL, scale = c("feature", "all", "none"), overlap = FALSE, axes = FALSE, combine = TRUE, coord.fixed = TRUE, border.by = NULL) {
  cells <- cells %||% Cells(x = object)
  scale <- scale[[1L]]
  scale <- match.arg(arg = scale)
  no.right <- theme(axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
                    axis.text.y.right = element_blank(), axis.title.y.right = element_text(face = "bold", 
                                                                                           size = 14, margin = margin(r = 7)))
  fov <- fov %||% DefaultFOV(object = object)
  fov <- Filter(f = function(x) {
    return(x %in% Images(object = object) && inherits(x = object[[x]], 
                                                      what = "FOV"))
  }, x = fov)
  if (!length(x = fov)) {
    stop("No compatible spatial coordinates present")
  }
  boundaries <- boundaries %||% sapply(X = fov, FUN = function(x) {
    return(DefaultBoundary(object = object[[x]]))
  }, simplify = FALSE, USE.NAMES = TRUE)
  boundaries <- Seurat:::.BoundariesByImage(object = object, fov = fov, 
                                            boundaries = boundaries)
  fov <- names(x = boundaries)
  border.by <- boundaries %!NA% border.by %||% "ident"
  if (isTRUE(x = blend) || !is.null(x = split.by)) {
    type <- ifelse(test = isTRUE(x = "blend"), yes = "Blended", 
                   no = "Split")
    if (length(x = fov) != 1L) {
      fov <- fov[1L]
      warning(type, " image feature plots can only be done on a single image, using \"", 
              fov, "\"", call. = FALSE, immediate. = TRUE)
    }
    if (any(!overlap) && length(x = boundaries[[fov]]) > 
        1L) {
      warning(type, " image feature plots require overlapped segmentations", 
              call. = FALSE, immediate. = TRUE)
    }
    overlap <- TRUE
  }
  overlap <- rep_len(x = overlap, length.out = length(x = fov))
  crop <- rep_len(x = crop, length.out = length(x = fov))
  names(x = crop) <- names(x = overlap) <- fov
  if (isTRUE(x = blend)) {
    if (length(x = features) != 2L) {
      stop("Blended feature plots only works with two features")
    }
    default.colors <- eval(expr = formals(fun = ImageFeaturePlot)$cols)
    cols <- switch(EXPR = as.character(x = length(x = cols)), 
                   `0` = {
                     warning("No colors provided, using default colors", 
                             immediate. = TRUE)
                     default.colors
                   }, `1` = {
                     warning("Only one color provided, assuming specified is double-negative and augmenting with default colors", 
                             immediate. = TRUE)
                     c(cols, default.colors[2:3])
                   }, `2` = {
                     warning("Only two colors provided, assuming specified are for features and augmenting with '", 
                             default.colors[1], "' for double-negatives", 
                             immediate. = TRUE)
                     c(default.colors[1], cols)
                   }, `3` = cols, {
                     warning("More than three colors provided, using only first three", 
                             immediate. = TRUE)
                     cols[1:3]
                   })
  }
  if (!is.null(border.by)){
    var = c(features,border.by ,split.by[1L])
  } else{
    var = c(features,split.by[1L])
  }
  md <- FetchData(object = object, vars = var, 
                  cells = cells)
  # if (!is.null(cells.use)){
  #   md[[features]][!(md[[border.by]] %in% cells.use)] <- 0.000
  # }
  split.by <- intersect(x = split.by, y = colnames(x = md))
  if (!length(x = split.by)) {
    split.by <- NULL
  }
  imax <- ifelse(test = is.null(x = split.by), yes = ncol(x = md), 
                 no = ncol(x = md) - length(x = split.by))
  features <- colnames(x = md)[1:imax-1]
  min.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = min(md[[feature]]), 
                  no = cutoff))
  }, cutoff = min.cutoff, feature = features)
  max.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = max(md[[feature]]), 
                  no = cutoff))
  }, cutoff = max.cutoff, feature = features)
  check.lengths <- unique(x = vapply(X = list(features, min.cutoff, 
                                              max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
  if (length(x = check.lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  brewer.gran <- ifelse(test = length(x = cols) == 1, yes = brewer.pal.info[cols, 
  ]$maxcolors, no = length(x = cols))
  for (i in seq_along(along.with = features)) {
    f <- features[[i]]
    data.feature <- md[[f]]
    min.use <- SetQuantile(cutoff = min.cutoff[i], data = data.feature)
    max.use <- SetQuantile(cutoff = max.cutoff[i], data = data.feature)
    data.feature[data.feature < min.use] <- min.use
    data.feature[data.feature > max.use] <- max.use
    if (brewer.gran != 2) {
      data.feature <- if (all(data.feature == 0)) {
        rep_len(x = 0, length.out = length(x = data.feature))
      }
      else {
        as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.feature), 
                                         breaks = brewer.gran)))
      }
    }
    md[[f]] <- data.feature
  }
  if (is.null(x = split.by)) {
    split.by <- RandomName()
    md[[split.by]] <- factor(x = split.by)
  }
  if (!is.factor(x = md[[split.by]])) {
    md[[split.by]] <- factor(x = md[[split.by]])
  }
  if (isTRUE(x = blend)) {
    md <- lapply(X = levels(x = md[[split.by]]), FUN = function(x) {
      df <- md[as.character(x = md[[split.by]]) == x, , 
               drop = FALSE]
      no.expression <- features[colMeans(x = df[, features]) == 
                                  0]
      if (length(x = no.expression)) {
        stop("The following features have no value: ", 
             paste(no.expression, collapse = ", "))
      }
      return(cbind(df[, split.by, drop = FALSE], BlendExpression(data = df[, 
                                                                           features])))
    })
    md <- do.call(what = "rbind", args = md)
    features <- setdiff(x = colnames(x = md), y = split.by)
  }
  pnames <- unlist(x = lapply(X = seq_along(along.with = fov), 
                              FUN = function(i) {
                                return(if (isTRUE(x = overlap[i])) {
                                  fov[i]
                                } else {
                                  paste(fov[i], boundaries[[i]], sep = "_")
                                })
                              }))
  pdata <- vector(mode = "list", length = length(x = pnames))
  names(x = pdata) <- pnames
  for (i in names(x = pdata)) {
    ul <- unlist(x = strsplit(x = i, split = "_"))
    if (length(ul) > 1) {
      img <- paste(ul[1:length(ul) - 1], collapse = "_")
      lyr <- ul[length(ul)]
    }
    else if (length(ul) == 1) {
      img <- ul[1]
      lyr <- "centroids"
    }
    else {
      stop("the length of ul is 0. please check.")
    }
    if (is.na(x = lyr)) {
      lyr <- boundaries[[img]]
    }
    pdata[[i]] <- lapply(X = lyr, FUN = function(l) {
      df <- fortify(model = object[[img]][[l]])
      df <- df[df$cell %in% cells, , drop = FALSE]
      if (!is.null(x = md)) {
        df <- merge(x = df, y = md, by.x = "cell", by.y = 0, 
                    all.x = TRUE)
      }
      df$cell <- paste(l, df$cell, sep = "_")
      df$boundary <- l
      return(df)
    })
    pdata[[i]] <- if (!is.na(x = pdata[[i]])) {
      do.call(what = "rbind", args = pdata[[i]])
    }
    else {
      unlist(x = pdata[[i]])
    }
  }
  if (!is.null(x = molecules)) {
    molecules <- .MolsByFOV(object = object, fov = fov, molecules = molecules)
    mdata <- vector(mode = "list", length = length(x = fov))
    names(x = mdata) <- fov
    for (img in names(x = mdata)) {
      idata <- object[[img]]
      if (!img %in% names(x = molecules)) {
        mdata[[img]] <- NULL
        next
      }
      if (isTRUE(x = crop[img])) {
        idata <- Overlay(x = idata, y = idata)
      }
      imols <- gsub(pattern = paste0("^", Key(object = idata)), 
                    replacement = "", x = molecules[[img]])
      mdata[[img]] <- FetchData(object = idata, vars = imols, 
                                nmols = nmols)
    }
  }
  else {
    mdata <- NULL
  }
  if (isTRUE(x = blend)) {
    ncol <- 4
    color.matrix <- BlendMatrix(two.colors = cols[2:3], col.threshold = blend.threshold, 
                                negative.color = cols[1])
    cols <- cols[2:3]
    colors <- list(color.matrix[, 1], color.matrix[1, ], 
                   as.vector(x = color.matrix))
    blend.legend <- BlendMap(color.matrix = color.matrix)
  }
  limits <- switch(EXPR = scale, all = range(unlist(x = md[, 
                                                           features])), NULL)
  plots <- vector(mode = "list", length = length(x = levels(x = md[[split.by]])))
  names(x = plots) <- levels(x = md[[split.by]])
  for (i in seq_along(along.with = levels(x = md[[split.by]]))) {
    ident <- levels(x = md[[split.by]])[i]
    plots[[ident]] <- vector(mode = "list", length = length(x = pdata))
    names(x = plots[[ident]]) <- names(x = pdata)
    if (isTRUE(x = blend)) {
      blend.key <- suppressMessages(expr = blend.legend + 
                                      scale_y_continuous(sec.axis = dup_axis(name = ifelse(test = length(x = levels(x = md[[split.by]])) > 
                                                                                             1, yes = ident, no = "")), expand = c(0, 0)) + 
                                      labs(x = features[1L], y = features[2L], title = if (i == 
                                                                                           1L) {
                                        paste("Color threshold:", blend.threshold)
                                      }
                                      else {
                                        NULL
                                      }) + no.right)
    }
    for (j in seq_along(along.with = pdata)) {
      key <- names(x = pdata)[j]
      img <- unlist(x = strsplit(x = key, split = "_"))[1L]  
      plots[[ident]][[key]] <- vector(mode = "list", length = length(x = features) + 
                                        ifelse(test = isTRUE(x = blend), yes = 1L, no = 0L))
      data.plot <- pdata[[j]][as.character(x = pdata[[j]][[split.by]]) == 
                                ident, , drop = FALSE]
      for (y in seq_along(along.with = features)) {
        feature <- features[y]
        cols.use <- if (isTRUE(x = blend)) {
          cc <- as.numeric(x = as.character(x = data.plot[, 
                                                          feature])) + 1
          colors[[y]][sort(unique(x = cc))]
        }
        else {
          NULL
        }
        colnames(data.plot) <- gsub("-", "_", colnames(data.plot))
        data.plot$boundary <- factor(data.plot$boundary)
        
        p <- SingleImagePlot_wyt(data = data.plot, cells.use = cells.use ,col.by = gsub("-", "_", feature), size = size, col.factor = blend, 
                                 cols = cols.use, molecules = mdata[[img]], 
                                 mols.size = mols.size, mols.cols = mols.cols, 
                                 alpha = alpha, border.color = border.color,  border.by =  border.by,
                                 border.size = border.size, dark.background = dark.background) + 
          CenterTitle() + labs(fill = feature)
        if (isTRUE(x = blend)) {
          p <- p + guides(fill = "none")
        }
        if (isTRUE(coord.fixed)) {
          p <- p + coord_fixed()
        }
        if (!isTRUE(x = axes)) {
          p <- p + NoAxes(panel.background = element_blank())
        }
        else if (isTRUE(x = blend) || length(x = levels(x = md[[split.by]])) > 
                 1L) {
          if (y != 1L) {
            p <- p + theme(axis.line.y = element_blank(), 
                           axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
                           axis.title.y.left = element_blank())
          }
          if (i != length(x = levels(x = md[[split.by]]))) {
            p <- p + theme(axis.line.x = element_blank(), 
                           axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
                           axis.title.x = element_blank())
          }
        }
        if (!isTRUE(x = blend)) {
          if (length(x = cols) == 1L) {
            p <- p + scale_fill_brewer(palette = cols)
          }
          else {
            cols.grad <- cols
            fexp <- data.plot[data.plot[[split.by]] == 
                                ident, feature, drop = TRUE]
            fexp <- unique(x = fexp)
            if (length(x = fexp) == 1L) {
              warning("All cells have the same value (", 
                      fexp, ") of ", feature, call. = FALSE, 
                      immediate. = TRUE)
              if (fexp == 0) {
                cols.grad <- cols.grad[1L]
              }
            }
            if (scale == "feature") {
              limits <- range(pdata[[j]][[feature]])
            }
            p <- p + ggplot2::scale_fill_gradientn(colors = cols.grad, 
                                                   guide = "colorbar", limits = limits)
          }
        }
        p <- p + if (i == 1L) {
          ggplot2::labs(title = feature)
        }
        else {
          ggplot2::labs(title = NULL)
        }
        plots[[ident]][[key]][[y]] <- p
      }
      if (isTRUE(x = blend)) {
        plots[[ident]][[key]][[length(x = plots[[ident]][[key]])]] <- blend.key
      }
      else if (length(x = levels(x = md[[split.by]])) > 
               1L) {
        plots[[ident]][[key]][[y]] <- suppressMessages(expr = plots[[ident]][[key]][[y]] + 
                                                         scale_y_continuous(sec.axis = dup_axis(name = ident)) + 
                                                         no.right)
      }
    }
    plots[[ident]] <- unlist(x = plots[[ident]], recursive = FALSE, 
                             use.names = FALSE)
  }
  plots <- unlist(x = plots, recursive = FALSE, use.names = FALSE)
  if (isTRUE(x = combine)) {
    if (isTRUE(x = blend) || length(x = levels(x = md[[split.by]])) > 
        1L) {
      plots <- patchwork::wrap_plots(plots, ncol = ifelse(test = isTRUE(x = blend), 
                                                          yes = 4L, no = length(x = features)), nrow = length(x = levels(x = md[[split.by]])), 
                                     guides = "collect")
    }
    else {
      plots <- patchwork::wrap_plots(plots)
    }
  }
  return(plots)
}

colors_48 <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628",
  "#FFFF33", "#F781BF", "#999999", "#E5D8BD", "#B3CDE3", "#CCEBC5",
  "#FED9A6", "#FBB4AE", "#8DD3C7", "#BEBADA", "#80B1D3", "#B3DE69",
  "#FCCDE5", "#BC80BD", "#FFED6F", "#8DA0CB", "#E78AC3", "#E5C494",
  "#CCCCCC", "#FB9A99", "#800080", "#CAB2D6", "#6A3D9A", "#B15928",
  "#66C2A5","#FC8D62","#A6D854", "#FFD92F", "#E7298A","#008B8B", "#FF00FF",
  "#8B4513", "#2E8B57", "#DA70D6", "#00CED1",
  "#FF4500", "#9400D3", "#32CD32", "#8A2BE2"
)

# spatial <- readRDS("/mnt/NexData_file/home/wyt/S203_project/project/xenium_FWSP20255467_human_skin/anno_A_NC/total_anno/xenium_1_human_skin_anno_drop_low.rds")
# spatial <- subset(spatial,sample == "P3-12")
# DefaultBoundary(spatial@images[['P3-12']]) <- 'segmentation'
# Idents(spatial) <- "celltype_fac"
# names(colors_48) <- levels(spatial)
# 
# p <- ImageFeaturePlot_wyt(spatial, features = "DMKN", size = 0.75, cols = c('#003366', '#FFD700'),
#                  fov ="P3-12",nmols = 20000,border.color = colors_48, border.size = 1,border.by = "celltype_fac",
#                  cells.use = c("FSCN1+basal_KC","CDSN+Keratinocytes","FSCN1+Keratinocytes"))
# ggsave(p,filename = "/mnt/NexData_file/home/wyt/S203_project/project/anno_1.png",
#        width = 20,
#        height = 20,
#        units = "in",
#        dpi = 300)
# 
# p <- ImageFeaturePlot_wyt(spatial, features = "DMKN", size = 0.75, cols = c('#003366', '#FFD700'),
#                           fov ="P3-12",nmols = 20000,border.color = colors_48, border.size = 1,border.by = "celltype_fac")
# 
# ggsave(p,filename = "/mnt/NexData_file/home/wyt/S203_project/project/anno_2.png",
#        width = 20,
#        height = 20,
#        units = "in",
#        dpi = 300)
# 
# 
# data.plot <- readRDS("/mnt/NexData_file/home/wyt/S203_project/project/data_plot_2.rds")
# debug(SingleImagePlot_wyt)
# undebug(SingleImagePlot_wyt)
# p <- SingleImagePlot_wyt(data = data.plot, col.by = "DMKN", size = 0.75, col.factor = FALSE,
#                      cols = NULL, molecules = NULL,
#                      mols.size = 0.1, mols.cols = NULL,cells.use = c("FSCN1+basal_KC","CDSN+Keratinocytes","FSCN1+Keratinocytes"),
#                      alpha = 1, border.color = colors_48,  border.by =  "celltype_fac",
#                      border.size = 1, dark.background = TRUE) +
#   CenterTitle() + labs(fill = "DMKN")


