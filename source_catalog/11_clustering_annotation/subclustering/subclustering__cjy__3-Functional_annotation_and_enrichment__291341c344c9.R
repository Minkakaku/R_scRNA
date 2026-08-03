library(data.table)
library(dplyr)
library(tidyverse)
library(ggsci)
library(ggh4x)
library(RColorBrewer)
library(corto)
library(ggrepel)
library(ggplot2)
library(readxl)
library(igraph)
library(ggraph)
library(ggpubr)
library(reshape2)
library(jsonlite)
############################## Parameters ######################################
# workflow = "Metabolomics" # Must be one of "Metabolomics", "Lipidomics"
# organism = "hsa" # Must be organism IDs from KEGG
# sample_info = "input/samples.csv"
# comparison = "input/comparison.txt"
# db_dir = "DB/HMDB_20250114"
# results_dir = "Results"

diff_file <- snakemake@input[["dep_result"]]
workflow <- snakemake@params[["workflow"]]
organism <- snakemake@params[["organism"]]
sample_info <- snakemake@input[["sample_info"]]
db_dir1 <- snakemake@params[["db_dir"]]
results_dir <- snakemake@params[["results_dir"]]

############################ Functions defined #################################
diff.classy.plot1 <- function(data, font = 14, color.pal = "Set3", angle = 75, margin = 15, margin2 = 10, type = "Metabolomics", label.font = 4) {
  if (type == "Metabolomics") {
    n <- brewer.pal.info[color.pal, "maxcolors"]
    if (length(unique(data$HMDB_SuperClass)) > n) {
      colors <- colorRampPalette(ggsci.colors[[color.pal]][1:n])(length(unique(data$HMDB_SuperClass)))
    } else {
      colors <- ggsci.colors[[color.pal]][1:length(unique(data$HMDB_SuperClass))]
    }
    names(colors) <- unique(data$HMDB_SuperClass)
    data <- data %>% arrange(-Freq)
    data$HMDB_SuperClass <- factor(data$HMDB_SuperClass, levels = unique(data$HMDB_SuperClass))
    group.count <- length(unique(data$HMDB_SuperClass))
    legend.ncol <- min(group.count, ifelse(group.count > 12, 4, 3))
    p1 <- ggplot(data, aes(reorder(HMDB_Class, -Freq), Freq, fill = HMDB_SuperClass)) +
      geom_col() +
      facet_grid(. ~ HMDB_SuperClass, scales = "free", space = "free") +
      scale_fill_manual(values = colors) +
      theme_classic(base_size = font) +
      scale_x_discrete(expand = expansion(add = 0.8)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      coord_cartesian(clip = "off") +
      theme(
        strip.text = element_blank(), axis.text.x = element_text(angle = angle, hjust = 1), axis.text = element_text(color = "black"),
        plot.margin = unit(c(margin2, 5, 5, margin), "mm"),
        legend.background = element_blank(), legend.box.background = element_blank()
      ) +
      geom_text(aes(label = Freq), size = label.font, vjust = 0) +
      labs(x = "HMDB Category", y = "Number of compounds") +
      theme(legend.position = "top", legend.title.position = "top") +
      guides(fill = guide_legend(ncol = legend.ncol))
  } else {
    colors <- ggsci.colors[[color.pal]][1:length(unique(data$Categories))]
    names(colors) <- unique(data$Categories)
    p1 <- ggplot(data, aes(reorder(Name, Freq), Freq, fill = Categories)) +
      geom_col() +
      facet_grid(. ~ Categories, scales = "free", space = "free") +
      scale_fill_manual(values = colors) +
      theme_classic(base_size = font) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      theme(
        strip.text = element_blank(), axis.text.x = element_text(angle = angle, hjust = 1), axis.text = element_text(color = "black"),
        plot.margin = unit(c(margin2, 5, 5, margin), "mm")
      ) +
      geom_text(aes(label = Freq), size = label.font, vjust = 0) +
      theme(
        legend.background = element_blank(),
        legend.box.background = element_blank()
      ) +
      labs(x = "Lipid category", y = "Number of compounds") +
      theme(legend.position = "top", legend.title.position = "top") +
      guides(fill = guide_legend(ncol = 2))
  }
  return(p1)
}
diff.classification.plot.size <- function(data, min.width = 12, min.height = 7) {
  if (nrow(data) == 0) {
    return(list(width = min.width, height = min.height))
  }
  item.count <- length(unique(data$HMDB_Class))
  group.count <- length(unique(data$HMDB_SuperClass))
  legend.ncol <- min(group.count, ifelse(group.count > 12, 4, 3))
  legend.rows <- ceiling(group.count / legend.ncol)
  width <- max(min.width, min(15, min.width + max(0, item.count - 40) * 0.09 + max(0, group.count - 12) * 0.1))
  height <- max(min.height, min(7.5, min.height + max(0, legend.rows - 2) * 0.12))
  list(width = width, height = height)
}
diff.classy.plot2 <- function(data, font = 14, color.pal = "Set2", angle = 60, margin = 15, margin2 = 10, label.font = 4) {
  colors <- ggsci.colors[[color.pal]][1:length(unique(data$Pathway.level.1))]
  names(colors) <- unique(data$Pathway.level.1)
  p2 <- ggplot(
    data,
    aes(reorder(Pathway.level.2, Freq), Freq, fill = Pathway.level.1)
  ) +
    geom_col() +
    facet_grid(. ~ Pathway.level.1, scales = "free", space = "free") +
    scale_fill_manual(values = colors) +
    theme_classic(base_size = 14) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme(
      strip.text = element_blank(), axis.text.x = element_text(angle = angle, hjust = 1), axis.text = element_text(color = "black"),
      plot.margin = unit(c(margin2, 5, 5, margin), "mm")
    ) +
    geom_text(aes(label = Freq), vjust = 0, size = label.font) +
    theme(
      legend.background = element_blank(),
      legend.box.background = element_blank()
    ) +
    labs(x = "KEGG pathway", y = "Number of compounds", fill = "Pathway Level 1") +
    theme(legend.position = "top", legend.title.position = "top") +
    guides(fill = guide_legend(ncol = 2))
  return(p2)
}
ggsci.colors <- list(
  "Set1" = brewer.pal(9, "Set1"),
  "Set2" = brewer.pal(8, "Set2"),
  "Set3" = brewer.pal(12, "Set3"),
  "Pastel1" = brewer.pal(9, "Pastel1"),
  "Pastel2" = brewer.pal(8, "Pastel2"),
  "Paired" = brewer.pal(12, "Paired"),
  "Dark2" = brewer.pal(8, "Dark2"),
  "Accent" = brewer.pal(8, "Accent"),
  "npg" = pal_npg()(10),
  "aaas" = pal_aaas()(10),
  "nejm" = pal_nejm()(8),
  "lancet" = pal_lancet()(9),
  "jama" = pal_jama()(7),
  "jco" = pal_jco()(10),
  "ucscgb" = pal_ucscgb()(26),
  "d3_category10" = pal_d3("category10")(10),
  "d3_category20" = pal_d3("category20")(20),
  "d3_category20c" = pal_d3("category20c")(20),
  "igv" = pal_igv()(51),
  "uchicago" = pal_uchicago()(9),
  "rickandmorty" = pal_rickandmorty()(12),
  "simpsons" = pal_simpsons()(16)
)
enrich.test <- function(x) {
  a <- x[1]
  b <- x[2] - a
  c <- x[3] - a
  d <- x[4] - a - b - c
  p.value <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")$p.value
}
enrich.demo <- function(go.data, diff.list) {
  go.data.diff <- subset(go.data, `Metabolite name` %in% diff.list)
  mapping <- data.frame(table(go.data.diff$Terms))
  background <- data.frame(table(go.data$Terms))
  if (nrow(mapping) == 0) {
    return(data.frame(Terms = character(), mapping = integer(), `Fold enrich` = numeric(), p.value = numeric(), fdr = numeric(), Metabolites = character(), check.names = FALSE))
  }
  go.data.state <- merge(mapping, background, by = "Var1")
  colnames(go.data.state) <- c("Terms", "mapping", "background")
  go.data.state$`All mapping` <- length(unique(go.data.diff$`Metabolite name`))
  go.data.state$`All background` <- length(unique(go.data$`Metabolite name`))
  go.data.state$`Fold enrich` <- round(go.data.state$mapping / go.data.state$background / go.data.state$`All mapping` * go.data.state$`All background`, 2)
  go.data.state <- subset(go.data.state, mapping > 1)
  if (nrow(go.data.state) == 0) {
    return(data.frame(Terms = character(), mapping = integer(), `Fold enrich` = numeric(), p.value = numeric(), fdr = numeric(), Metabolites = character(), check.names = FALSE))
  }
  go.data.state$p.value <- signif(apply(go.data.state[2:5], 1, enrich.test), 3)
  go.data.state$fdr <- signif(p.adjust(go.data.state$p.value), 3)
  go.data.state <- subset(go.data.state, select = c("Terms", "mapping", "Fold enrich", "p.value", "fdr")) %>%
    merge(subset(go.data.diff, select = c("Terms", "Metabolite name")), by = "Terms") %>%
    group_by(Terms, mapping, `Fold enrich`, p.value, fdr) %>%
    summarize(Metabolites = paste(`Metabolite name`, collapse = " ; "))
  return(go.data.state)
}
gseaplot <- function(t) {
  temp.data <- data.frame(`Runing Enrichment Score` = t$running_score, reflist = t$reflist, inSet = t$inSet, Rank = rank(-t$reflist))
  if (t$nes < 0) {
    line.color <- "darkgreen"
  } else {
    line.color <- "red"
  }
  ggarrange(
    ggplot(temp.data) +
      geom_col(aes(Rank, reflist), color = "steelblue") +
      theme_classic() +
      theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()
      ) +
      annotate(
        geom = "text", color = line.color, label = sprintf("P value=%.0e; ES=%.3f; NES=%.3f", t$p.value, t$es, t$nes), x = nrow(temp.data) / 2,
        y = min(temp.data$reflist) / 2
      ) +
      labs(title = t$name, y = "Log2fc") +
      theme(plot.title = element_text(hjust = 0.5, size = 10)),
    ggplot(temp.data) +
      geom_line(aes(Rank, t$running_score), color = line.color) +
      geom_point(data = temp.data %>% filter(inSet == 1), aes(x = Rank, y = 0), shape = 3, size = 5) +
      theme_classic() +
      geom_hline(yintercept = 0) +
      labs(y = "Running Enrichment Score"),
    nrow = 2, heights = c(1, 2), align = "v"
  )
}
getpath <- function(obo, term) {
  p <- obo$parents[[term]]
  if (length(p) > 0) {
    out.data <- data.frame(source = p, target = term)
    out <- rbindlist(lapply(obo$parents[[term]], function(x) getpath(obo, x)))
    return(rbind(out.data, out))
  } else {
    return(NULL)
  }
}
obo.plot <- function(enrich.out, obo, col.1 = "blue", col.2 = "red") {
  enrich.out %>% slice_min(p.value, n = 5)
  network <- rbindlist(lapply(enrich.out$ID[enrich.out$p.value < 0.05], function(x) getpath(obo, x))) %>% distinct_all()
  node <- data.frame(ID = unique(c(network$source, network$target))) %>% merge(enrich.out, by = "ID", all.x = T)
  ig <- graph_from_data_frame(network, directed = T, vertices = node)
  ggraph(ig) + geom_edge_link(color = "black", arrow = arrow(length = unit(0.02, "npc"))) +
    scale_fill_gradientn(colours = c(col.1, col.2), limits = c(0, 5), oob = scales::squish) + theme_void() +
    geom_node_point(aes(size = `Fold enrich`, fill = -log10(p.value)), shape = 21) +
    geom_node_text(aes(label = stringr::str_wrap(Terms, width = 20)), size = 4, repel = T) +
    labs(title = enrich.out$compare[1]) + theme(plot.title = element_text(face = "bold"))
}

######################## KEGG dotplot JSON generation #########################
generate_kegg_dotplot_json <- function(data, output_path) {
  sig_vals <- sapply(data$p.value, function(p) {
    if (is.na(p)) {
      return("")
    }
    if (p < 0.01) {
      return("**")
    }
    if (p < 0.05) {
      return("*")
    }
    return("")
  })
  json_obj <- list(
    width = 800, height = 600, outerBorder = TRUE,
    margin = list(top = 40, right = 60, bottom = 60, left = 60),
    dataset = list(list(
      source = list(
        terms = as.character(data$Terms),
        FoldEnrich = round(as.numeric(data$`Fold enrich`), 4),
        pvalue = round(as.numeric(data$p.value), 4),
        mapping = as.integer(data$mapping),
        sig = as.character(sig_vals)
      )
    )),
    series = list(list(
      type = "bubble", datasetIndex = 0,
      encode = list(
        x = "FoldEnrich", y = "terms", size = "mapping", color = "pvalue",
        tooltip = list("terms", "FoldEnrich", "pvalue", "mapping", "sig"),
        annotation = "sig"
      ),
      border = list(show = TRUE)
    )),
    xAxis = list(list(
      type = "scaleLinear", name = "Fold enrich",
      nameTextStyle = list(color = "", fontFamily = "", fontSize = "16"),
      axisLabel = list(
        show = TRUE, color = "#000", fontStyle = "normal",
        fontWeight = "normal", fontFamily = "sans-serif", fontSize = "14", rotate = 0
      )
    )),
    yAxis = list(list(
      type = "scaleBand", name = "KEGG pathway",
      nameTextStyle = list(color = "", fontFamily = "", fontSize = "16"),
      axisLabel = list(
        show = TRUE, color = "#000", fontStyle = "normal",
        fontWeight = "normal", fontFamily = "sans-serif", fontSize = "14", rotate = 0
      )
    )),
    legend = list(list(show = FALSE, site = "right", orient = "vertical", top = 30, type = "circle", value = list())),
    sizeLegend = list(
      show = TRUE, title = "mapping", orient = "vertical",
      count = 4, right = 10, top = 200, bordered = TRUE, color = "#ffffff", fontSize = 12, site = "right"
    ),
    colorbar = list(
      show = TRUE, title = "p.value", width = 20, height = 150,
      site = "right", right = 10, top = 30, fontSize = 12, tickCount = 5, useScientificNotation = FALSE
    )
  )
  write_json(json_obj, path = output_path, pretty = TRUE, auto_unbox = TRUE)
}

######################## MSEA dotplot JSON generation ##########################
generate_msea_dotplot_json <- function(data, output_path) {
  neg_log10_pval <- round(-log10(data$p.value), 2)
  sig_vals <- sapply(data$p.value, function(p) {
    if (is.na(p)) {
      return("")
    }
    if (p < 0.01) {
      return("**")
    }
    if (p < 0.05) {
      return("*")
    }
    return("")
  })
  count_range <- range(data$Count)
  size_values <- as.integer(seq(count_range[1], count_range[2], length.out = 4))
  size_values <- unique(as.integer(round(size_values / 4) * 4))
  if (length(size_values) < 2) size_values <- as.integer(count_range)
  x_max <- ceiling(max(neg_log10_pval))
  json_obj <- list(
    width = 800, height = 600, outerBorder = TRUE,
    margin = list(top = 30, right = 60, bottom = 60, left = 420),
    dataset = list(list(
      source = list(
        terms = as.character(data$Terms),
        negLog10PValue = neg_log10_pval,
        NES = round(as.numeric(data$NES), 2),
        Count = as.integer(data$Count),
        sig = as.character(sig_vals)
      )
    )),
    series = list(list(
      type = "bubble", datasetIndex = 0,
      encode = list(
        x = "negLog10PValue", y = "terms", size = "Count", color = "NES",
        tooltip = list("terms", "negLog10PValue", "NES", "Count", "sig"),
        annotation = "sig"
      ),
      border = list(show = TRUE)
    )),
    xAxis = list(list(
      type = "scaleLinear", name = "-log10(p.value)", min = 0, max = x_max, tickCount = x_max + 1,
      nameTextStyle = list(color = "", fontFamily = "", fontSize = "16"),
      axisLabel = list(
        show = TRUE, color = "#000", fontStyle = "normal",
        fontWeight = "normal", fontFamily = "sans-serif", fontSize = "14", rotate = 0
      )
    )),
    yAxis = list(list(
      type = "scaleBand", name = "Terms",
      nameTextStyle = list(color = "", fontFamily = "", fontSize = "16"),
      axisLabel = list(
        show = TRUE, color = "#000", fontStyle = "normal",
        fontWeight = "normal", fontFamily = "sans-serif", fontSize = "14", rotate = 0
      )
    )),
    legend = list(list(show = FALSE, site = "right", orient = "vertical", top = 30, type = "circle", value = list())),
    sizeLegend = list(
      show = TRUE, title = "Count", orient = "vertical",
      count = 4, value = size_values, right = 10, top = 200, bordered = TRUE, color = "#ffffff", fontSize = 12, site = "right"
    ),
    colorbar = list(
      show = TRUE, title = "NES", width = 20, height = 150,
      colors = list("#10b0a4", "#ffffff", "#ef5337"), centerValue = 0, domain = list(-2, 2), reverse = FALSE,
      site = "right", right = 10, top = 30, fontSize = 12, tickCount = 5, useScientificNotation = FALSE
    )
  )
  write_json(json_obj, path = output_path, pretty = TRUE, auto_unbox = TRUE)
}

######################## Functional enrichment  ################################
set.seed(123)
compare.name <- sub(sub(diff_file, pattern = ".*of ", replacement = ""), pattern = ".xls$", replacement = "")
if (workflow == "Metabolomics") {
  print(db_dir1)
  if (db_dir1 != "-") {
    db_dir <- db_dir1
  } else {
    db_dir <- "/mnt/NexData_file/cosmicdarahub/Bioinformatics/Database/Metabolism/DB/HMDB_20250114"
  }
  db_dir <- gsub("/home/ztt/data/", "/mnt/NexData_file/cosmicdarahub/Bioinformatics/Database/", db_dir)
  db <- sub(".*\\/", "", db_dir)

  data.diff <- fread(diff_file)

  kegg_org_dir <- "resources/DB/KEGG_Org_20260428"
  kegg_org_file <- paste0(kegg_org_dir, "/", organism)
  cat("Reading local KEGG pathway list from:", kegg_org_file, "\n")
  kegg.data <- fread(paste0(db_dir, "/KEGG_compound_pathway_annotation.tsv")) %>%
    filter(`Pathway ID` %in% sub(fread(kegg_org_file, header = F)$V1, pattern = organism, replacement = "map"))
  diff.data <- separate_rows(data.diff, cpd_id, sep = "; ")

  if (db == "HMDB_20250114") {
    data.kegg <- merge(diff.data, kegg.data, by = c("cpd_id", "HMDB ID")) %>%
      distinct() %>%
      dplyr::select(`Pathway ID`, Terms, `Pathway level 1`, `Pathway level 2`, everything())
  } else {
    data.kegg <- merge(diff.data, kegg.data, by = c("cpd_id")) %>%
      distinct() %>%
      dplyr::select(`Pathway ID`, Terms, `Pathway level 1`, `Pathway level 2`, everything())
  }

  if (length(data.kegg$`Metabolite name`[data.kegg$regulation != "unchanged"]) > 0) {
    kegg.enrich <- enrich.demo(data.kegg %>% filter(`Pathway level 1` == "Metabolism") %>%
      distinct(Terms, `Metabolite name`), data.diff$`Metabolite name`[data.diff$regulation != "unchanged"])
  } else {
    kegg.enrich <- NULL
  }
  data.kegg.sub <- data.kegg %>%
    group_by(Terms) %>%
    filter(length(Terms) > 5, `Pathway level 1` == "Metabolism") %>%
    ungroup() %>%
    select(!cpd_id) %>%
    distinct()
  data.kegg.anno <- data.kegg %>%
    group_by(Terms) %>%
    summarise(Metabolites = paste(unique(`Metabolite name`), collapse = "; "))

  if (!is.null(kegg.enrich) && nrow(kegg.enrich) > 0) {
    if (nrow(kegg.enrich) > 30) {
      kegg.tmp <- kegg.enrich %>%
        as.data.frame() %>%
        slice_min(p.value, n = 30, with_ties = FALSE)
    } else {
      kegg.tmp <- kegg.enrich
    }
    xmax1 <- as.integer(max(kegg.tmp$`Fold enrich`)) + 1
    xmin1 <- as.integer(min(kegg.tmp$`Fold enrich`)) - 0.3
    enrich.plot1 <- ggplot(kegg.tmp, aes(`Fold enrich`, reorder(Terms, `Fold enrich`))) +
      geom_point(aes(size = mapping, fill = p.value), shape = 21) +
      theme_test(base_size = 14) +
      scale_size(range = c(3, 6)) +
      scale_fill_material("red", reverse = T) +
      xlim(xmin1 - 0.05 * xmin1, xmax1) +
      geom_text(aes(label = case_when(p.value < 0.01 ~ "**", p.value < 0.05 ~ "*")), hjust = 2, vjust = 0.8, size = 5) +
      labs(y = "KEGG pathway", x = "Fold enrich")
  } else {
    enrich.plot1 <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No terms enriched.", hjust = 0.5, vjust = 0.5, size = 6) +
      xlim(0, 1) +
      ylim(0, 1) +
      theme_void()
  }

  #---------MSEA enrichment -----------#
  set.seed(123)
  if (nrow(data.kegg.sub) > 0) {
    gsea.res <- lapply(split(data.kegg.sub, data.kegg.sub$`Pathway ID`), function(x) {
      t <- gsea(setNames(data.diff$log2fc, data.diff$`Metabolite name`), x$`Metabolite name`, method = "permutation", np = 500)
      t$name <- x$Terms[1]
      t$ID <- x$`Pathway ID`[1]
      t$N <- sum(t$inSet)
      return(t)
    })
    gsea.enrich <- rbindlist(lapply(gsea.res, function(x) {
      data.frame(
        Terms = x$name,
        NES = signif(x$nes, 3),
        p.value = x$p.value,
        Count = x$N
      )
    }), idcol = "Pathway ID") %>%
      merge(data.kegg.anno, by = "Terms")
  } else {
    gsea.res <- list()
    gsea.enrich <- data.frame(
      `Pathway ID` = character(),
      Terms = character(),
      NES = numeric(),
      p.value = numeric(),
      Count = integer(),
      Metabolites = character(),
      check.names = FALSE
    )
  }
  tmp.data <- gsea.enrich %>% filter(p.value < 0.05)
  if (!is.null(gsea.enrich) && nrow(tmp.data) > 0) {
    xmax2 <- as.integer(max(-log10(tmp.data$p.value))) + 1
    xmin2 <- 0
    enrich.plot2 <- ggplot(tmp.data, aes(-log10(p.value), reorder(Terms, -log10(p.value)), fill = NES)) +
      geom_point(aes(size = Count), shape = 21) +
      scale_fill_gradientn(colors = c("#00A087FF", "white", "#E64B35FF"), limits = c(-2, 2), oob = scales::squish) +
      theme_test(base_size = 14) +
      guides(color = "none") +
      xlim(xmin2, xmax2) +
      geom_text(aes(label = case_when(p.value < 0.01 ~ "**", p.value < 0.05 ~ "*")), hjust = 2, vjust = 0.8, size = 5) +
      labs(y = "Terms", x = "-log10(p.value)")
  } else {
    enrich.plot2 <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No terms enriched.", hjust = 0.5, vjust = 0.5, size = 6) +
      xlim(0, 1) +
      ylim(0, 1) +
      theme_void()
  }
  if (nrow(tmp.data) > 0) {
    res.plot <- lapply(gsea.res[sapply(gsea.res, function(x) x$p.value < 0.05)], gseaplot)
    ncol <- 4
    nrow <- ceiling(length(res.plot) / ncol)
    enrich.plot3 <- ggarrange(plotlist = res.plot, ncol = ncol, nrow = nrow)
  } else {
    res.plot <- NULL
    nrow <- 1
    enrich.plot3 <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No terms enriched.", hjust = 0.5, vjust = 0.5, size = 6) +
      xlim(0, 1) +
      ylim(0, 1) +
      theme_void()
  }

  # Diff compound classification
  classify.data1 <- data.diff %>%
    filter(regulation != "unchanged") %>%
    group_by(HMDB_SuperClass, HMDB_Class) %>%
    summarize(Freq = length(unique(`Metabolite name`)), Metabolites = paste(`Metabolite name`, collapse = " ; ")) %>%
    ungroup()
  if (nrow(classify.data1) > 0) {
    classify.plot.data1 <- classify.data1 %>% slice_max(Freq, n = 40)
    comclass.plot1 <- diff.classy.plot1(classify.plot.data1, font = 14, color.pal = "Set3", type = workflow, angle = 75, margin = 16, margin2 = 10, label.font = 4)
  } else {
    classify.plot.data1 <- classify.data1
    comclass.plot1 <- ggplot() +
      theme(panel.background = element_blank())
  }

  classify.data2 <- data.kegg %>%
    filter(regulation != "unchanged") %>%
    group_by(`Pathway level 1`, `Pathway level 2`) %>%
    summarize(Freq = length(unique(`Metabolite name`)), HMDBs = paste(`HMDB ID`, collapse = ";"), Metabolites = paste(`Metabolite name`, collapse = " ; ")) %>%
    ungroup()
  if (nrow(classify.data2) > 0) {
    colnames(classify.data2) <- gsub(" ", ".", colnames(classify.data2))
    comclass.plot2 <- diff.classy.plot2(classify.data2 %>% slice_max(Freq, n = 40), font = 14, color.pal = "Set2", angle = 60, margin = 18, margin2 = 14, label.font = 4)
  } else {
    comclass.plot2 <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No pathway mapping.", hjust = 0.5, vjust = 0.5, size = 6) +
      xlim(0, 1) +
      ylim(0, 1) +
      theme_void()
  }
}

if (workflow == "Lipidomics") {
  if (db_dir1 != "-") {
    db_dir <- db_dir1
  } else {
    db_dir <- "/mnt/NexData_file/cosmicdarahub/Bioinformatics/Database/Metabolism/DB/LION_20250116"
  }
  db_dir <- gsub("/home/ztt/data/", "/mnt/NexData_file/cosmicdarahub/Bioinformatics/Database/", db_dir)
  sample <- fread(sample_info)
  if (all(sample$Sample == "-") == T) {
    sample <- sample %>% mutate(Sample = ID)
  } else {
    sample <- sample
  }
  meta.data <- read_excel(paste(results_dir, "1-MS_identified_summary/Table 1.2 MS identified lipids detail_normalized.xlsx", sep = "/"))
  data.diff <- fread(diff_file)

  obo <- readRDS(paste0(db_dir, "/obo.RDS"))
  values <- list()
  #----------LIMAPS enrichment-----------#
  limap.enrich <- enrich.demo(data.diff %>% dplyr::rename(Terms = "Name"), data.diff$`Metabolite name`[data.diff$regulation != "unchanged"])

  #----------LION enrichment-------------#
  LION.annot <- rbindlist(lapply(split(meta.data, meta.data$`Metabolite name`), function(x) {
    id <- unique(getpath(obo, x$`LION ID`[1])$source)
    data.frame(ID = id, Name = obo$name[id])
  }), idcol = "Metabolite name", fill = T)
  LION.enrich <- enrich.demo(
    LION.annot %>% dplyr::rename(Terms = "Name"),
    data.diff$`Metabolite name`[data.diff$regulation != "unchanged"]
  ) %>%
    merge(LION.annot %>% distinct(Name, ID), by.x = "Terms", by.y = "Name")

  s <- unlist(str_split(compare.name, "-vs-"))

  if (nrow(limap.enrich) > 0) {
    xmax1 <- as.integer(max(limap.enrich$`Fold enrich`) + 1)
    xmin1 <- as.integer(min(limap.enrich$`Fold enrich`))
    enrich.plot1 <- ggplot(limap.enrich, aes(`Fold enrich`, reorder(Terms, `Fold enrich`))) +
      geom_point(aes(size = mapping, fill = p.value), shape = 21) +
      theme_test(base_size = 14) +
      scale_size(range = c(3, 6)) +
      scale_fill_material("red", reverse = T) +
      geom_text(aes(label = case_when(p.value < 0.01 ~ "**", p.value < 0.05 ~ "*")), hjust = 2, vjust = 0.7, size = 5) +
      labs(y = "LIMAPS Name", x = "Fold enrich") +
      xlim(xmin1 - 0.05 * xmin1, xmax1)
  } else {
    enrich.plot1 <- ggplot() +
      theme(panel.background = element_blank())
  }
  tmp <- LION.enrich %>% filter(p.value < 0.05)

  if (nrow(tmp) > 0) {
    xmax2 <- as.integer(max(LION.enrich$`Fold enrich`) + 1)
    xmin2 <- as.integer(min(LION.enrich$`Fold enrich`))
    enrich.plot2 <- ggplot(LION.enrich %>% filter(p.value < 0.05) %>% slice_min(p.value, n = 30), aes(`Fold enrich`, reorder(Terms, `Fold enrich`))) +
      geom_point(aes(size = mapping, fill = p.value), shape = 21) +
      theme_test(base_size = 14) +
      scale_size(range = c(3, 6)) +
      scale_fill_material("red", reverse = T) +
      geom_text(aes(label = case_when(p.value < 0.0001 ~ "***", p.value < 0.01 ~ "**", p.value < 0.05 ~ "*")), hjust = 1.5, vjust = 0.7, size = 5) +
      labs(y = "Lipid ontology", x = "Fold enrich") +
      xlim(xmin2 - 0.05 * xmin2, xmax2)
    enrich.plot3 <- obo.plot(LION.enrich, obo = obo)
  } else {
    enrich.plot2 <- ggplot() +
      theme(panel.background = element_blank())
    enrich.plot3 <- ggplot() +
      theme(panel.background = element_blank())
  }
  # Diff compound classification
  classify.data1 <- data.diff %>%
    filter(regulation != "unchanged") %>%
    ungroup() %>%
    group_by(Categories, Name) %>%
    summarize(Freq = length(unique(`Metabolite name`)), Metabolites = paste(`Metabolite name`, collapse = " ; "))
  comclass.plot1 <- diff.classy.plot1(classify.data1 %>% slice_max(Freq, n = 40), font = 14, color.pal = "Set3", type = workflow, angle = 60, margin = 10, margin2 = 10)
  classify.data2 <- data.diff %>%
    melt(measure.vars = sample$Sample[sample$Group %in% s], variable.name = "Sample") %>%
    merge(sample, by = "Sample") %>%
    filter(Group %in% s) %>%
    filter(regulation != "unchanged") %>%
    group_by(Name, Group) %>%
    summarize(mean = mean(value), ymin = mean_se(value)[[2]], ymax = mean_se(value)[[3]])
  comclass.plot2 <- classify.data2 %>% ggplot(aes(x = Name, y = mean, fill = Group)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    geom_errorbar(aes(x = Name, ymin = ymin, ymax = ymax), stat = "identity", position = position_dodge2(width = 0.4, padding = 0.5)) +
    theme_bw(base_size = 14) +
    scale_fill_manual(values = c("#71D0F5FF", "#197EC0FF")) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1), plot.margin = unit(c(5, 5, 5, 10), "mm")) +
    labs(x = "Category", y = "Lipid Abundance")
}

########################## Export Results #####################################
temp.dir <- file.path(results_dir, "4-Functional Annotation and Enrichment")
if (!dir.exists(temp.dir)) dir.create(path = temp.dir)


if (workflow == "Metabolomics") {
  dir.create(paste0(temp.dir, "/msea_ind"))

  fig4a.size <- diff.classification.plot.size(classify.plot.data1)
  ggsave(comclass.plot1, filename = paste0(temp.dir, "/Figure 4a. Diff Metabolite Classification by Category of ", compare.name, ".pdf"), width = fig4a.size$width, height = fig4a.size$height)
  ggsave(comclass.plot1, filename = paste0(temp.dir, "/Figure 4a. Diff Metabolite Classification by Category of ", compare.name, ".png"), width = fig4a.size$width, height = fig4a.size$height)

  # Figure 4a JSON (Metabolomics)
  if (nrow(classify.data1) > 0) {
    fig4a_data <- classify.data1 %>%
      slice_max(Freq, n = 40) %>%
      arrange(-Freq)
    fig4a_data$HMDB_SuperClass <- factor(fig4a_data$HMDB_SuperClass, levels = unique(fig4a_data$HMDB_SuperClass))
    fig4a_data <- fig4a_data[order(fig4a_data$HMDB_SuperClass, -fig4a_data$Freq), ]
    n_cat <- length(unique(fig4a_data$HMDB_SuperClass))
    n_max <- brewer.pal.info["Set3", "maxcolors"]
    if (n_cat > n_max) {
      fig4a_colors <- colorRampPalette(ggsci.colors[["Set3"]][1:n_max])(n_cat)
    } else {
      fig4a_colors <- ggsci.colors[["Set3"]][1:n_cat]
    }
    names(fig4a_colors) <- levels(fig4a_data$HMDB_SuperClass)
    fig4a_json <- list(
      width = 1200, height = 600,
      margin = list(top = 40, right = 360, bottom = 160, left = 100),
      dataset = list(list(source = list(
        subcategory = as.character(fig4a_data$HMDB_Class),
        count = as.integer(fig4a_data$Freq),
        category = as.character(fig4a_data$HMDB_SuperClass)
      ))),
      xAxis = list(list(type = "scaleBand", show = TRUE, silent = TRUE)),
      grid = list(list(width = "68%", height = "50%")),
      yAxis = list(list(type = "scaleLinear", show = TRUE, name = "Number of compounds", axisGrid = list(show = TRUE))),
      series = list(list(
        type = "hierarchicalBar", datasetIndex = 0,
        encode = list(x = "subcategory", y = "count", color = "category"),
        label = list(show = TRUE),
        border = list(show = FALSE),
        color = as.list(unname(fig4a_colors))
      )),
      legend = list(list(show = TRUE, title = "HMDB SuperClass", orient = "vertical", right = "14%"))
    )
    write_json(fig4a_json, path = paste0(temp.dir, "/Figure 4a. Diff Metabolite Classification by Category of ", compare.name, ".json"), pretty = TRUE, auto_unbox = TRUE)
  }

  ggsave(comclass.plot2, filename = paste0(temp.dir, "/Figure 4b. Diff Metabolite Classification by Pathway of ", compare.name, ".pdf"), width = 10, height = 7)
  ggsave(comclass.plot2, filename = paste0(temp.dir, "/Figure 4b. Diff Metabolite Classification by Pathway of ", compare.name, ".png"), width = 10, height = 7)

  # Figure 4b JSON (Metabolomics)
  if (nrow(classify.data2) > 0) {
    fig4b_data <- classify.data2 %>%
      ungroup() %>%
      slice_max(Freq, n = 40) %>%
      arrange(-Freq)
    fig4b_data <- fig4b_data[order(fig4b_data$Pathway.level.1, -fig4b_data$Freq), ]
    n_cat2 <- length(unique(fig4b_data$Pathway.level.1))
    fig4b_colors <- ggsci.colors[["Set2"]][1:n_cat2]
    names(fig4b_colors) <- unique(fig4b_data$Pathway.level.1)
    fig4b_json <- list(
      width = 1200, height = 600,
      margin = list(top = 40, right = 360, bottom = 160, left = 100),
      dataset = list(list(source = list(
        subcategory = as.character(fig4b_data$Pathway.level.2),
        count = as.integer(fig4b_data$Freq),
        category = as.character(fig4b_data$Pathway.level.1)
      ))),
      xAxis = list(list(type = "scaleBand", show = TRUE, silent = TRUE)),
      grid = list(list(width = "68%", height = "50%")),
      yAxis = list(list(type = "scaleLinear", show = TRUE, name = "Number of compounds", axisGrid = list(show = TRUE))),
      series = list(list(
        type = "hierarchicalBar", datasetIndex = 0,
        encode = list(x = "subcategory", y = "count", color = "category"),
        label = list(show = TRUE),
        border = list(show = FALSE),
        color = as.list(unname(fig4b_colors))
      )),
      legend = list(list(show = TRUE, title = "Pathway Level 1", orient = "vertical", right = "14%"))
    )
    write_json(fig4b_json, path = paste0(temp.dir, "/Figure 4b. Diff Metabolite Classification by Pathway of ", compare.name, ".json"), pretty = TRUE, auto_unbox = TRUE)
  }

  fwrite(data.kegg, file = paste0(temp.dir, "/Table 4a.KEGG Annotation of ", compare.name, ".xls"), sep = "\t")
  if (is.null(kegg.enrich)) {
    fwrite(
      data.table(
        Terms = character(),
        mapping = integer(),
        `Fold enrich` = numeric(),
        p.value = numeric(),
        fdr = numeric(),
        Metabolites = character()
      ),
      file = paste0(temp.dir, "/Table 4b.KEGG Enrichment of ", compare.name, ".xls"),
      sep = "\t"
    )
  } else {
    fwrite(kegg.enrich, file = paste0(temp.dir, "/Table 4b.KEGG Enrichment of ", compare.name, ".xls"), sep = "\t")
  }
  ggsave(enrich.plot1, filename = paste0(temp.dir, "/Figure 4c. KEGG enrichment of ", compare.name, ".pdf"), width = 8, height = 6)
  ggsave(enrich.plot1, filename = paste0(temp.dir, "/Figure 4c. KEGG enrichment of ", compare.name, ".png"), width = 8, height = 6)
  if (!is.null(kegg.enrich) && nrow(kegg.enrich) > 0) {
    kegg_json_data <- (if (nrow(kegg.enrich) > 30) {
      kegg.enrich %>%
        as.data.frame() %>%
        slice_min(p.value, n = 30, with_ties = FALSE)
    } else {
      kegg.enrich %>% as.data.frame()
    }) %>% arrange(`Fold enrich`)
    generate_kegg_dotplot_json(data = kegg_json_data, output_path = paste0(temp.dir, "/Figure 4c. KEGG enrichment of ", compare.name, ".json"))
  }
  fwrite(gsea.enrich, file = paste0(temp.dir, "/Table 4c.MSEA Enrichment of ", compare.name, ".xls"), sep = "\t")
  ggsave(enrich.plot2, filename = paste0(temp.dir, "/Figure 4d. MSEA enrichment1 of ", compare.name, ".pdf"), width = 8, height = 6)
  ggsave(enrich.plot2, filename = paste0(temp.dir, "/Figure 4d. MSEA enrichment1 of ", compare.name, ".png"), width = 8, height = 6)
  if (!is.null(gsea.enrich) && nrow(tmp.data) > 0) {
    msea_json_data <- tmp.data %>%
      as.data.frame() %>%
      arrange(-log10(p.value))
    generate_msea_dotplot_json(data = msea_json_data, output_path = paste0(temp.dir, "/Figure 4d. MSEA enrichment1 of ", compare.name, ".json"))
  }
  ggsave(enrich.plot3, filename = paste0(temp.dir, "/Figure 4e. MSEA enrichment2 of ", compare.name, ".pdf"), width = 12, height = 3 * nrow)
  ggsave(enrich.plot3, filename = paste0(temp.dir, "/Figure 4e. MSEA enrichment2 of ", compare.name, ".png"), width = 12, height = 3 * nrow, bg = "white")

  dir.create(paste0(temp.dir, "/msea_ind/", compare.name))
  saveRDS(gsea.res, file = paste0(temp.dir, "/MSEA res of ", compare.name, ".Rds"))
  if (!is.null(res.plot)) {
    lapply(names(res.plot), function(y) {
      p <- plot(res.plot[[y]])
      ggsave(paste0(temp.dir, "/msea_ind/", compare.name, "/MSEA Enrichment plot for ", y, ".pdf"), plot = p, width = 8, height = 8)
      ggsave(paste0(temp.dir, "/msea_ind/", compare.name, "/MSEA Enrichment plot for ", y, ".png"), plot = p, width = 8, height = 8)
    })
  } else {
    p <- ggplot() +
      theme(panel.background = element_blank())
    ggsave(paste0(temp.dir, "/msea_ind/", compare.name, "/MSEA Enrichment plot.pdf"), plot = p, width = 8, height = 8)
    ggsave(paste0(temp.dir, "/msea_ind/", compare.name, "/MSEA Enrichment plot.png"), plot = p, width = 8, height = 8)
  }
}
if (workflow == "Lipidomics") {
  fwrite(limap.enrich, file = paste0(temp.dir, "/Table 4a.LIPID MAPS Enrichment of ", compare.name, ".xls"), sep = "\t")
  fwrite(LION.enrich, file = paste0(temp.dir, "/Table 4b.LION Enrichment of ", compare.name, ".xls"), sep = "\t")

  ggsave(comclass.plot1, filename = paste0(temp.dir, "/Figure 4a. Diff Lipid Classification by Category of ", compare.name, ".pdf"), width = 10, height = 6)
  ggsave(comclass.plot1, filename = paste0(temp.dir, "/Figure 4a. Diff Lipid Classification by Category of ", compare.name, ".png"), width = 10, height = 6)

  # Figure 4a JSON (Lipidomics)
  if (nrow(classify.data1) > 0) {
    fig4a_data <- classify.data1 %>%
      slice_max(Freq, n = 40) %>%
      arrange(-Freq)
    fig4a_data <- fig4a_data[order(fig4a_data$Categories, -fig4a_data$Freq), ]
    n_cat <- length(unique(fig4a_data$Categories))
    fig4a_colors <- ggsci.colors[["Set3"]][1:n_cat]
    names(fig4a_colors) <- unique(fig4a_data$Categories)
    fig4a_json <- list(
      width = 1200, height = 600,
      margin = list(top = 40, right = 360, bottom = 160, left = 100),
      dataset = list(list(source = list(
        subcategory = as.character(fig4a_data$Name),
        count = as.integer(fig4a_data$Freq),
        category = as.character(fig4a_data$Categories)
      ))),
      xAxis = list(list(type = "scaleBand", show = TRUE, silent = TRUE)),
      grid = list(list(width = "68%", height = "50%")),
      yAxis = list(list(type = "scaleLinear", show = TRUE, name = "Number of compounds", axisGrid = list(show = TRUE))),
      series = list(list(
        type = "hierarchicalBar", datasetIndex = 0,
        encode = list(x = "subcategory", y = "count", color = "category"),
        label = list(show = TRUE),
        border = list(show = FALSE),
        color = as.list(unname(fig4a_colors))
      )),
      legend = list(list(show = TRUE, title = "Lipid Category", orient = "vertical", right = "14%"))
    )
    write_json(fig4a_json, path = paste0(temp.dir, "/Figure 4a. Diff Lipid Classification by Category of ", compare.name, ".json"), pretty = TRUE, auto_unbox = TRUE)
  }
  ggsave(comclass.plot2, filename = paste0(temp.dir, "/Figure 4b. Diff Lipid Abundance by Category of ", compare.name, ".pdf"), width = 10, height = 6)
  ggsave(comclass.plot2, filename = paste0(temp.dir, "/Figure 4b. Diff Lipid Abundance by Category of ", compare.name, ".png"), width = 10, height = 6)
  ggsave(enrich.plot1, filename = paste0(temp.dir, "/Figure 4c. LIPID MAPS enrichment of ", compare.name, ".pdf"), width = 8, height = 6)
  ggsave(enrich.plot1, filename = paste0(temp.dir, "/Figure 4c. LIPID MAPS enrichment of ", compare.name, ".png"), width = 8, height = 6)
  if (nrow(limap.enrich) > 0) {
    limap_json_data <- limap.enrich %>%
      as.data.frame() %>%
      arrange(`Fold enrich`)
    generate_dotplot_json(
      data = limap_json_data, terms_col = "Terms", x_col = "Fold enrich", color_col = "p.value", size_col = "mapping",
      x_axis_name = "Fold enrich", y_axis_name = "LIMAPS Name", colorbar_title = "P value",
      output_path = paste0(temp.dir, "/Figure 4c. LIPID MAPS enrichment of ", compare.name, ".json")
    )
  }
  ggsave(enrich.plot2, filename = paste0(temp.dir, "/Figure 4d. LION enrichment1 of ", compare.name, ".pdf"), width = 8, height = 6)
  ggsave(enrich.plot2, filename = paste0(temp.dir, "/Figure 4d. LION enrichment1 of ", compare.name, ".png"), width = 8, height = 6)
  if (nrow(tmp) > 0) {
    lion_json_data <- (LION.enrich %>% filter(p.value < 0.05) %>% slice_min(p.value, n = 30)) %>%
      as.data.frame() %>%
      arrange(`Fold enrich`)
    generate_dotplot_json(
      data = lion_json_data, terms_col = "Terms", x_col = "Fold enrich", color_col = "p.value", size_col = "mapping",
      x_axis_name = "Fold enrich", y_axis_name = "Lipid ontology", colorbar_title = "P value",
      output_path = paste0(temp.dir, "/Figure 4d. LION enrichment1 of ", compare.name, ".json")
    )
  }
  ggsave(enrich.plot3, filename = paste0(temp.dir, "/Figure 4e. LION enrichment2 of ", compare.name, ".pdf"), width = 11, height = 11)
  ggsave(enrich.plot3, filename = paste0(temp.dir, "/Figure 4e. LION enrichment2 of ", compare.name, ".png"), width = 11, height = 11, bg = "white")
}
