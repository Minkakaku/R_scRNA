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
# sample_info = "input/samples.csv"
# comparison = "input/comparison.txt"
# results_dir = "Results"

diff_file <- snakemake@input[["dep_result"]]
workflow <- "Lipidomics"
sample_info <- snakemake@input[["sample_info"]]
db_dir <- snakemake@params[["db_dir"]]
results_dir <- snakemake@params[["results_dir"]]
data_final_table <- if ("data_final_table" %in% names(snakemake@input)) {
  snakemake@input[["data_final_table"]]
} else {
  paste(results_dir, "1-MS_identified_summary/Table 1.2 MS identified lipids detail_normalized.xlsx", sep = "/")
}

############################ Functions defined #################################
diff.classy.plot1 <- function(data, font = 14, color.pal = "Set3", angle = 75, margin = 15, margin2 = 10, label.font = 4) {
  colors <- ggsci.colors[[color.pal]][1:length(unique(data$Categories))]
  names(colors) <- unique(data$Categories)
  p1 <- ggplot(data, aes(reorder(Name, Freq), Freq, fill = Categories)) +
    geom_col() +
    facet_grid(. ~ Categories, scales = "free", space = "free") +
    scale_fill_manual(values = colors) +
    theme_classic(base_size = font) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme(strip.text = element_blank(), axis.text.x = element_text(angle = angle, hjust = 1), axis.text = element_text(color = "black"), plot.margin = unit(c(margin2, 5, 5, margin), "mm")) +
    geom_text(aes(label = Freq), size = label.font, vjust = 0) +
    theme(legend.background = element_blank(), legend.box.background = element_blank()) +
    labs(x = "Lipid category", y = "Number of compounds") +
    theme(legend.position = "top", legend.title.position = "top") +
    guides(fill = guide_legend(ncol = 2))
  return(p1)
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

generate_dotplot_json <- function(data, terms_col, x_col, color_col, size_col, x_axis_name, y_axis_name, colorbar_title, output_path) {
  sig_vals <- sapply(data[[color_col]], function(p) {
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
  size_range <- range(data[[size_col]], na.rm = TRUE)
  size_values <- unique(as.integer(round(seq(size_range[1], size_range[2], length.out = 4))))
  if (length(size_values) < 2) {
    size_values <- as.integer(size_range)
  }
  json_obj <- list(
    width = 800, height = 600, outerBorder = TRUE,
    margin = list(top = 40, right = 60, bottom = 60, left = 300),
    dataset = list(list(
      source = list(
        terms = as.character(data[[terms_col]]),
        xValue = round(as.numeric(data[[x_col]]), 4),
        colorValue = round(as.numeric(data[[color_col]]), 4),
        sizeValue = as.integer(data[[size_col]]),
        sig = as.character(sig_vals)
      )
    )),
    series = list(list(
      type = "bubble", datasetIndex = 0,
      encode = list(
        x = "xValue", y = "terms", size = "sizeValue", color = "colorValue",
        tooltip = list("terms", "xValue", "colorValue", "sizeValue", "sig"),
        annotation = "sig"
      ),
      border = list(show = TRUE)
    )),
    xAxis = list(list(
      type = "scaleLinear", name = x_axis_name,
      nameTextStyle = list(color = "", fontFamily = "", fontSize = "16"),
      axisLabel = list(
        show = TRUE, color = "#000", fontStyle = "normal",
        fontWeight = "normal", fontFamily = "sans-serif", fontSize = "14", rotate = 0
      )
    )),
    yAxis = list(list(
      type = "scaleBand", name = y_axis_name,
      nameTextStyle = list(color = "", fontFamily = "", fontSize = "16"),
      axisLabel = list(
        show = TRUE, color = "#000", fontStyle = "normal",
        fontWeight = "normal", fontFamily = "sans-serif", fontSize = "14", rotate = 0
      )
    )),
    legend = list(list(show = FALSE, site = "right", orient = "vertical", top = 30, type = "circle", value = list())),
    sizeLegend = list(
      show = TRUE, title = size_col, orient = "vertical",
      count = 4, value = size_values, right = 10, top = 200, bordered = TRUE, color = "#ffffff", fontSize = 12, site = "right"
    ),
    colorbar = list(
      show = TRUE, title = colorbar_title, width = 20, height = 150,
      site = "right", right = 10, top = 30, fontSize = 12, tickCount = 5, useScientificNotation = FALSE
    )
  )
  write_json(json_obj, path = output_path, pretty = TRUE, auto_unbox = TRUE)
}

######################## Functional enrichment  ################################
set.seed(123)
compare.name <- sub(sub(diff_file, pattern = ".*of ", replacement = ""), pattern = ".xls$", replacement = "")

sample <- fread(sample_info)
if (all(sample$Sample == "-") == T) {
  sample <- sample %>% mutate(Sample = ID)
  sample <- sample
}
meta.data <- read_excel(data_final_table)
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
comclass.plot1.data <- classify.data1 %>% slice_max(Freq, n = 40)
comclass.plot1 <- diff.classy.plot1(comclass.plot1.data, font = 14, color.pal = "Set3", angle = 60, margin = 10, margin2 = 12, label.font = 4)
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
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1),
    axis.title.x = element_text(margin = margin(t = 2)),
    plot.margin = unit(c(5, 8, 8, 36), "mm")
  ) +
  labs(x = "Category", y = "Lipid Abundance")

########################## Export Results #####################################
temp.dir <- file.path(results_dir, "4-Functional Annotation and Enrichment")
if (!dir.exists(temp.dir)) dir.create(path = temp.dir)


fwrite(limap.enrich, file = paste0(temp.dir, "/Table 4a.LIPID MAPS Enrichment of ", compare.name, ".xls"), sep = "\t")
fwrite(LION.enrich, file = paste0(temp.dir, "/Table 4b.LION Enrichment of ", compare.name, ".xls"), sep = "\t")

ggsave(comclass.plot1, filename = paste0(temp.dir, "/Figure 4a. Diff Lipid Classification by Category of ", compare.name, ".pdf"), width = 16, height = 8)
ggsave(comclass.plot1, filename = paste0(temp.dir, "/Figure 4a. Diff Lipid Classification by Category of ", compare.name, ".png"), width = 16, height = 8)

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
ggsave(comclass.plot2, filename = paste0(temp.dir, "/Figure 4b. Diff Lipid Abundance by Category of ", compare.name, ".pdf"), width = 10, height = 7)
ggsave(comclass.plot2, filename = paste0(temp.dir, "/Figure 4b. Diff Lipid Abundance by Category of ", compare.name, ".png"), width = 10, height = 7)
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
