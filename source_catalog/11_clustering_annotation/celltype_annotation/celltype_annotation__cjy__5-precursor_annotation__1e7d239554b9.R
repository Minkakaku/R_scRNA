# ============================================================
# Step 5: 前体蛋白功能注释 (基于 DB 文件)
# - GO 分类 (BP/MF/CC) — from GO_annotation.tsv
# - KEGG pathway 分类 — from KEGG_annotation.tsv
# - 亚细胞定位 — from Subcell_location.xls
# - 合并注释表输出
# ============================================================

source("scripts/utils.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

cat("========== Step 5: Precursor Protein Annotation (DB-based) ==========\n")

cfg <- load_config(snakemake@params[["config"]])

# ---------- 加载 DB 注释 ----------
db <- load_db_annotation(cfg$annotation$db_dir)

# ---------- 读取多肽信息，提取唯一蛋白 ----------
pep_info <- fread(snakemake@input[["peptide_info"]])

# 从 Protein.Ids 列提取 UniProt accession
all_proteins <- unique(na.omit(pep_info$Protein.Ids))
all_proteins <- all_proteins[all_proteins != ""]
# 处理可能的多 ID 分号分隔
all_proteins <- unique(unlist(strsplit(all_proteins, ";")))
all_proteins <- trimws(all_proteins)
all_proteins <- all_proteins[all_proteins != ""]

all_protein_groups <- unique(na.omit(pep_info$Protein.Group))
all_protein_groups <- all_protein_groups[all_protein_groups != ""]

cat("  Unique proteins (UniProt): ", length(all_proteins), "\n")
cat("  Unique Protein.Group: ", length(all_protein_groups), "\n")

# ============================================================
# Fig10: GO 功能分类 (Level 2)
# ============================================================
cat("[1/4] GO classification from DB...\n")

go_df <- NULL
if (!is.null(db$go)) {
  # 筛选本项目鉴定到的蛋白
  go_matched <- db$go[db$go$`Protein accession` %in% all_proteins, ]
  cat("  GO matched records: ", nrow(go_matched), "\n")

  if (nrow(go_matched) > 0) {
    # 按 Category + Terms 统计蛋白数
    go_stats <- go_matched[, .(Count = uniqueN(`Protein accession`)),
      by = .(Category, Terms)
    ]
    setnames(go_stats, c("Ontology", "Description", "Count"))

    # 映射 Category 名称
    go_stats$Ontology <- gsub("Biological process", "BP", go_stats$Ontology)
    go_stats$Ontology <- gsub("Molecular function", "MF", go_stats$Ontology)
    go_stats$Ontology <- gsub("Cellular component", "CC", go_stats$Ontology)
    go_stats <- go_stats[go_stats$Ontology %in% c("BP", "MF", "CC"), ]

    # 每个 Ontology 取 top 10
    go_split <- split(go_stats, go_stats$Ontology)
    go_df <- do.call(rbind, lapply(go_split, function(x) {
      head(x[order(-x$Count), ], 10)
    }))
  }
}

if (!is.null(go_df) && nrow(go_df) > 0) {
  go_df$Ontology <- factor(go_df$Ontology, levels = c("BP", "CC", "MF"))

  n_terms <- nrow(go_df)
  fig_h <- adaptive_height(n_terms, per_item = 0.35, min_h = 6, max_h = 14)

  p10 <- ggplot(go_df, aes(x = reorder(Description, Count), y = Count, fill = Ontology)) +
    geom_bar(stat = "identity", alpha = 0.85, width = 0.7) +
    geom_text(aes(label = Count), hjust = -0.15, size = 4, fontface = "bold") +
    scale_fill_manual(values = ONT_COLORS) +
    coord_flip() +
    facet_grid(Ontology ~ ., scales = "free_y", space = "free_y") +
    labs(
      title = "GO Classification of Precursor Proteins (Level 2)",
      x = NULL, y = "Number of Proteins"
    ) +
    theme_sci() +
    theme(
      legend.position = "none",
      strip.text.y = element_text(angle = 0, size = 15, face = "bold"),
      axis.text.y = element_text(size = 14, face = "bold")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12)))
  save_plot(p10, snakemake@output[["fig_go"]], width = 10, height = fig_h)
} else {
  save_empty_plot(snakemake@output[["fig_go"]], "No GO classification results")
}

# ============================================================
# Fig11: KEGG 通路分类
# ============================================================
cat("[2/4] KEGG classification from DB...\n")

kegg_df <- NULL
if (!is.null(db$kegg)) {
  kegg_matched <- db$kegg[db$kegg$`Protein accession` %in% all_proteins, ]
  cat("  KEGG matched records: ", nrow(kegg_matched), "\n")

  if (nrow(kegg_matched) > 0) {
    kegg_stats <- kegg_matched[, .(Count = uniqueN(`Protein accession`)),
      by = .(Terms, `Pathway level 1`, `Pathway level 2`)
    ]
    setnames(kegg_stats, c("Pathway", "Level1", "Level2", "Count"))

    # 按 Level2 汇总
    kegg_l2 <- kegg_stats[, .(Count = sum(Count)), by = .(Level1, Level2)]
    kegg_l2 <- kegg_l2[order(-Count), ]
    kegg_df <- head(kegg_l2, 25)
  }
}

if (!is.null(kegg_df) && nrow(kegg_df) > 0) {
  # 按 Level1 分面展示，每个 facet 内 Level2 按 Count 排序
  kegg_df$Level1 <- factor(kegg_df$Level1)
  fig_h <- adaptive_height(nrow(kegg_df), per_item = 0.42, min_h = 6, max_h = 16)

  # 构建 Level1 颜色映射
  l1_levels <- levels(kegg_df$Level1)
  l1_colors <- KEGG_L1_COLORS[l1_levels]
  missing <- is.na(l1_colors)
  if (any(missing)) {
    extra <- SCI_COLORS_NPG[seq_len(sum(missing))]
    l1_colors[missing] <- extra
  }
  names(l1_colors) <- l1_levels

  p11 <- ggplot(kegg_df, aes(x = reorder(Level2, Count), y = Count, fill = Level1)) +
    geom_bar(stat = "identity", alpha = 0.85, width = 0.7) +
    geom_text(aes(label = Count), hjust = -0.15, size = 4, fontface = "bold") +
    scale_fill_manual(values = l1_colors) +
    coord_flip() +
    facet_grid(Level1 ~ ., scales = "free_y", space = "free_y") +
    labs(
      title = "KEGG Pathway Classification of Precursor Proteins",
      x = NULL, y = "Number of Proteins"
    ) +
    theme_sci() +
    theme(
      legend.position = "none",
      strip.text.y = element_text(angle = 0, size = 13, face = "bold"),
      axis.text.y = element_text(size = 13, face = "bold")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12)))
  save_plot(p11, snakemake@output[["fig_kegg"]], width = 12, height = fig_h)
} else {
  save_empty_plot(snakemake@output[["fig_kegg"]], "No KEGG classification results")
}

# ============================================================
# Fig12: 亚细胞定位
# ============================================================
cat("[3/4] Subcellular localization from DB...\n")

subcell_df <- NULL
if (!is.null(db$subcell)) {
  sc_matched <- db$subcell[db$subcell$`Protein accession` %in% all_proteins, ]
  cat("  Subcell matched records: ", nrow(sc_matched), "\n")

  if (nrow(sc_matched) > 0) {
    sc_stats <- sc_matched[, .(Count = .N), by = .(`Subcellular location`)]
    setnames(sc_stats, c("Location", "Count"))
    sc_stats <- sc_stats[order(-Count), ]
    subcell_df <- head(sc_stats, 12)
    subcell_df$Percent <- round(subcell_df$Count / sum(subcell_df$Count) * 100, 1)
  }
}

if (!is.null(subcell_df) && nrow(subcell_df) > 0) {
  n_loc <- nrow(subcell_df)
  loc_colors <- colorRampPalette(SCI_COLORS_NPG)(n_loc)
  subcell_df$LocationLabel <- vapply(
    subcell_df$Location,
    function(x) paste(strwrap(x, width = 28), collapse = "\n"),
    character(1)
  )
  subcell_xmax <- max(subcell_df$Count, na.rm = TRUE) * 1.38

  p12 <- ggplot(subcell_df, aes(x = reorder(LocationLabel, Count), y = Count, fill = Location)) +
    geom_bar(stat = "identity", alpha = 0.85, width = 0.7) +
    geom_text(aes(label = paste0(Count, "  (", Percent, "%", ")")),
      hjust = -0.05, size = 4.2, fontface = "bold"
    ) +
    scale_fill_manual(values = setNames(loc_colors, subcell_df$Location[order(-subcell_df$Count)])) +
    coord_flip() +
    labs(
      title = "Subcellular Localization of Precursor Proteins",
      x = NULL, y = "Number of Proteins"
    ) +
    theme_sci() +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 13, face = "bold", lineheight = 0.9),
      plot.margin = margin(10, 18, 10, 18)
    ) +
    scale_y_continuous(limits = c(0, subcell_xmax), expand = expansion(mult = c(0, 0.02)))
  save_plot(p12, snakemake@output[["fig_subcell"]],
    width = 10,
    height = adaptive_height(n_loc, 0.4, 5, 10)
  )
} else {
  save_empty_plot(snakemake@output[["fig_subcell"]], "No subcellular localization results")
}

# ============================================================
# Table01: 合并注释表（含 Peptide_Count / GO / KEGG / Subcell）
# ============================================================
cat("[4/4] Building annotation table...\n")

annot_table <- data.frame(
  Protein = all_proteins,
  Gene = NA_character_,
  stringsAsFactors = FALSE
)

# ── 基因名（从 master annotation）──────────────────────────
if (!is.null(db$master)) {
  master <- as.data.frame(db$master)
  gene_map <- master[, intersect(c("Entry", "Gene names"), colnames(master))]
  gene_map <- gene_map[!duplicated(gene_map$Entry), ]
  annot_table <- merge(annot_table, gene_map, by.x = "Protein", by.y = "Entry", all.x = TRUE)
  if ("Gene names" %in% colnames(annot_table)) {
    annot_table$Gene <- ifelse(!is.na(annot_table$`Gene names`) & annot_table$`Gene names` != "",
      annot_table$`Gene names`, annot_table$Gene
    )
    annot_table$`Gene names` <- NULL
  }
}

# Protein description from product.tsv
if (!is.null(db$product)) {
  product_df <- as.data.frame(db$product)
  if (all(c("Protein accession", "Protein description") %in% colnames(product_df))) {
    product_df <- product_df[!duplicated(product_df$`Protein accession`), ]
    annot_table <- merge(
      annot_table,
      product_df[, c("Protein accession", "Protein description")],
      by.x = "Protein",
      by.y = "Protein accession",
      all.x = TRUE
    )
    colnames(annot_table)[colnames(annot_table) == "Protein description"] <- "Description"
  }
}

# ── 多肽数目 ──────────────────────────────────────────────
pep_count_map <- table(trimws(unlist(strsplit(as.character(pep_info$Protein.Ids), ";"))))
annot_table$Peptide_Count <- as.integer(pep_count_map[annot_table$Protein])
annot_table$Peptide_Count[is.na(annot_table$Peptide_Count)] <- 0L

# ── GO 注释（BP / MF / CC，条目用 "; " 拼接）─────────────
if (!is.null(db$go)) {
  go_df <- as.data.frame(db$go)
  go_df <- go_df[go_df$`Protein accession` %in% all_proteins, ]
  collapse_terms <- function(prot, ont) {
    sub <- go_df[go_df$`Protein accession` == prot &
      grepl(ont, go_df$Category, ignore.case = TRUE), "Terms"]
    if (length(sub) == 0) {
      return(NA_character_)
    }
    paste(unique(sub), collapse = "; ")
  }
  annot_table$GO_BP <- sapply(annot_table$Protein, collapse_terms, ont = "Biological")
  annot_table$GO_MF <- sapply(annot_table$Protein, collapse_terms, ont = "Molecular")
  annot_table$GO_CC <- sapply(annot_table$Protein, collapse_terms, ont = "Cellular")
}

# ── KEGG 通路（通路名用 "; " 拼接）────────────────────────
if (!is.null(db$kegg)) {
  kegg_df <- as.data.frame(db$kegg)
  kegg_matched <- kegg_df[kegg_df$`Protein accession` %in% all_proteins, ]
  kegg_map <- tapply(
    kegg_matched$Terms, kegg_matched$`Protein accession`,
    function(v) paste(unique(v), collapse = "; ")
  )
  annot_table$KEGG <- kegg_map[annot_table$Protein]
}

# ── 亚细胞定位 ────────────────────────────────────────────
if (!is.null(db$subcell)) {
  sc <- as.data.frame(db$subcell)
  sc_dedup <- sc[!duplicated(sc$`Protein accession`), ]
  annot_table <- merge(annot_table, sc_dedup[, c("Protein accession", "Subcellular location")],
    by.x = "Protein", by.y = "Protein accession", all.x = TRUE
  )
  colnames(annot_table)[colnames(annot_table) == "Subcellular location"] <- "Subcellular"
}

# ── NOG 功能类别 ──────────────────────────────────────────
if (!is.null(db$nog)) {
  nog_df <- as.data.frame(db$nog)
  if ("Protein accession" %in% colnames(nog_df) && "Terms" %in% colnames(nog_df)) {
    nog_matched <- nog_df[nog_df$`Protein accession` %in% all_proteins, ]
    nog_map <- tapply(
      nog_matched$Terms, nog_matched$`Protein accession`,
      function(v) paste(unique(v), collapse = "; ")
    )
    annot_table$NOG <- nog_map[annot_table$Protein]
  }
}

# 列重排
col_order <- c(
  "Protein", "Gene", "Peptide_Count",
  if ("Description" %in% colnames(annot_table)) "Description",
  if ("GO_BP" %in% colnames(annot_table)) c("GO_BP", "GO_MF", "GO_CC"),
  if ("KEGG" %in% colnames(annot_table)) "KEGG",
  if ("Subcellular" %in% colnames(annot_table)) "Subcellular",
  if ("NOG" %in% colnames(annot_table)) "NOG"
)
col_order <- col_order[col_order %in% colnames(annot_table)]
annot_table <- annot_table[, col_order]

write.table(annot_table, snakemake@output[["table_annot"]],
  sep = "\t", row.names = FALSE, quote = FALSE, na = ""
)

# ============================================================
# 分类注释表：GO / KEGG / 亚细胞定位（保留原始列，补充 Gene / Description）
# ============================================================

# 构建 Protein accession -> Gene / Description 映射
meta_map <- data.frame(
  `Protein accession` = annot_table$Protein,
  Gene                = annot_table$Gene,
  Description         = if ("Description" %in% colnames(annot_table)) annot_table$Description else NA_character_,
  stringsAsFactors    = FALSE,
  check.names         = FALSE
)

# 辅助函数：筛选鉴定蛋白 + join Gene/Description，保留全部原始列
join_meta <- function(raw_df, meta) {
  raw_df <- as.data.frame(raw_df)
  raw_df <- raw_df[raw_df[["Protein accession"]] %in% all_proteins, , drop = FALSE]
  merged <- merge(meta, raw_df, by = "Protein accession", all.y = TRUE)
  # 列顺序：Protein accession | Gene | Description | 其余原始列
  orig_extra <- setdiff(colnames(raw_df), "Protein accession")
  col_order <- c("Protein accession", "Gene", "Description", orig_extra)
  col_order <- col_order[col_order %in% colnames(merged)]
  merged[, col_order, drop = FALSE]
}

# -- GO 原始注释表 --
if (!is.null(snakemake@output[["table_go"]]) && !is.null(db$go)) {
  go_tbl <- join_meta(db$go, meta_map)
  write.table(go_tbl, snakemake@output[["table_go"]],
    sep = "\t", row.names = FALSE, quote = FALSE, na = ""
  )
  cat("  GO annotation table rows:", nrow(go_tbl), "\n")
} else if (!is.null(snakemake@output[["table_go"]])) {
  write.table(data.frame(), snakemake@output[["table_go"]],
    sep = "\t", row.names = FALSE, quote = FALSE, na = ""
  )
}

# -- KEGG 原始注释表 --
if (!is.null(snakemake@output[["table_kegg"]]) && !is.null(db$kegg)) {
  kegg_tbl <- join_meta(db$kegg, meta_map)
  write.table(kegg_tbl, snakemake@output[["table_kegg"]],
    sep = "\t", row.names = FALSE, quote = FALSE, na = ""
  )
  cat("  KEGG annotation table rows:", nrow(kegg_tbl), "\n")
} else if (!is.null(snakemake@output[["table_kegg"]])) {
  write.table(data.frame(), snakemake@output[["table_kegg"]],
    sep = "\t", row.names = FALSE, quote = FALSE, na = ""
  )
}

# -- 亚细胞定位原始注释表 --
if (!is.null(snakemake@output[["table_subcell"]]) && !is.null(db$subcell)) {
  sc_tbl <- join_meta(db$subcell, meta_map)
  write.table(sc_tbl, snakemake@output[["table_subcell"]],
    sep = "\t", row.names = FALSE, quote = FALSE, na = ""
  )
  cat("  Subcell annotation table rows:", nrow(sc_tbl), "\n")
} else if (!is.null(snakemake@output[["table_subcell"]])) {
  write.table(data.frame(), snakemake@output[["table_subcell"]],
    sep = "\t", row.names = FALSE, quote = FALSE, na = ""
  )
}

cat("  Annotation table rows:", nrow(annot_table), "\n")
cat("  Columns:", paste(colnames(annot_table), collapse = ", "), "\n")
cat("========== Step 5 Complete ==========\n")
