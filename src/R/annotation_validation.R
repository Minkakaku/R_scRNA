# 比较已有标签/聚类的描述性指标；不替代注释复核，也不重做聚类。

# 计算各计数项包含的无序对数。
pair_count <- function(x) {
  x <- as.numeric(x)
  sum(x * (x - 1) / 2)
}

# 计算两套分区的调整 Rand 指数。
adjusted_rand_index <- function(reference, clusters) {
  stopifnot(length(reference) == length(clusters), !anyNA(reference), !anyNA(clusters))
  if (!length(reference)) return(NA_real_)
  tab <- table(reference, clusters)
  n <- sum(tab)
  if (n < 2) return(NA_real_)
  index <- pair_count(tab)
  row_pairs <- pair_count(rowSums(tab))
  col_pairs <- pair_count(colSums(tab))
  total_pairs <- n * (n - 1) / 2
  expected <- row_pairs * col_pairs / total_pairs
  maximum <- (row_pairs + col_pairs) / 2
  if (maximum == expected) return(1)
  (index - expected) / (maximum - expected)
}

# 计算按几何平均熵归一化的互信息。
normalized_mutual_information <- function(reference, clusters) {
  stopifnot(length(reference) == length(clusters), !anyNA(reference), !anyNA(clusters))
  if (!length(reference)) return(NA_real_)
  tab <- table(reference, clusters)
  probability <- tab / sum(tab)
  row_probability <- rowSums(probability)
  col_probability <- colSums(probability)
  expected <- outer(row_probability, col_probability)
  keep <- probability > 0
  mutual_information <- sum(probability[keep] * log(probability[keep] / expected[keep]))
  reference_entropy <- -sum(row_probability[row_probability > 0] * log(row_probability[row_probability > 0]))
  cluster_entropy <- -sum(col_probability[col_probability > 0] * log(col_probability[col_probability > 0]))
  denominator <- sqrt(reference_entropy * cluster_entropy)
  if (denominator == 0) return(if (reference_entropy == 0 && cluster_entropy == 0) 1 else 0)
  mutual_information / denominator
}

# 以参考标签计算聚类纯度。
cluster_purity <- function(reference, clusters) {
  stopifnot(length(reference) == length(clusters), !anyNA(reference), !anyNA(clusters))
  if (!length(reference)) return(NA_real_)
  tab <- table(reference, clusters)
  sum(apply(tab, 2, max)) / sum(tab)
}

# 在给定坐标上计算各标签对自身与其他质心的相对距离。
centroid_separation_score <- function(coordinates, labels) {
  stopifnot(is.matrix(coordinates), nrow(coordinates) == length(labels))
  labels <- as.character(labels)
  keep <- complete.cases(coordinates) & !is.na(labels) & nzchar(labels)
  coordinates <- coordinates[keep, , drop = FALSE]
  labels <- labels[keep]
  label_levels <- sort(unique(labels))
  if (length(label_levels) < 2) return(NA_real_)
  centroids <- do.call(rbind, lapply(label_levels, function(label) {
    colMeans(coordinates[labels == label, , drop = FALSE])
  }))
  rownames(centroids) <- label_levels
  distances <- vapply(seq_len(nrow(centroids)), function(index) {
    sqrt(rowSums((coordinates - matrix(
      centroids[index, ], nrow(coordinates), ncol(coordinates), byrow = TRUE
    ))^2))
  }, numeric(nrow(coordinates)))
  if (is.null(dim(distances))) distances <- matrix(distances, ncol = 1)
  own_index <- match(labels, label_levels)
  own_distance <- distances[cbind(seq_len(nrow(distances)), own_index)]
  distances[cbind(seq_len(nrow(distances)), own_index)] <- Inf
  other_distance <- apply(distances, 1, min)
  mean((other_distance - own_distance) / pmax(other_distance, own_distance, .Machine$double.eps))
}
