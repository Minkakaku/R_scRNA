library(data.table)
library(dplyr)

print("doneload")

annotation.files <- c(
  eggnog_mapper_annotation = "eggnog_mapper_annotation.xls",
  go_annotation = "go_annotation.xls",
  kegg_pathway_annotation = "kegg_pathway_annotation.xls",
  kegg_ko_annotation = "kegg_ko_annotation.xls",
  pfam_annotation = "pfam_annotation.xls",
  cog_annotation = "cog_annotation.xls",
  string_mapping = "string_mapping.txt"
)
functional.annotation.names <- setdiff(names(annotation.files), "string_mapping")

project.outputs <- vapply(
  names(annotation.files),
  function(name) snakemake@output[[name]],
  character(1)
)

stringdb.version <- sub("^v", "", as.character(snakemake@params[["stringdb_version"]]))
stringdb.dir <- snakemake@params[["stringdb_dir"]]
configured.tax.id <- as.character(snakemake@params[["tax_id"]])
if (length(configured.tax.id) == 0 || is.na(configured.tax.id)) {
  configured.tax.id <- ""
}
configured.tax.id <- trimws(configured.tax.id)
if (configured.tax.id == "-") {
  configured.tax.id <- ""
}
db.annotation.dir <- snakemake@params[["annotation_dir"]]
db.files <- setNames(file.path(db.annotation.dir, annotation.files), names(annotation.files))
db.meta.file <- file.path(db.annotation.dir, "annotation_cache_meta.txt")

is_annotation_complete <- function(paths) {
  version.line <- paste0("STRINGdb_version=", stringdb.version)
  tax.id.line <- paste0("tax_id=", configured.tax.id)
  all(file.exists(paths)) &&
    all(file.info(paths)$size > 0) &&
    file.exists(db.meta.file) &&
    version.line %in% readLines(db.meta.file, warn = FALSE) &&
    tax.id.line %in% readLines(db.meta.file, warn = FALSE)
}

is_file_group_complete <- function(paths) {
  all(file.exists(paths)) && all(file.info(paths)$size > 0)
}

is_nonempty_file <- function(path) {
  file.exists(path) && !is.na(file.info(path)$size) && file.info(path)$size > 0
}

is_valid_gzip <- function(path) {
  is_nonempty_file(path) &&
    identical(system2("gzip", c("-t", path), stdout = FALSE, stderr = FALSE), 0L)
}

has_fasta_record_in_gzip <- function(path) {
  if (!is_valid_gzip(path)) {
    return(FALSE)
  }
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  tryCatch(
    any(grepl("^>", readLines(con, n = 100, warn = FALSE))),
    error = function(e) FALSE
  )
}

decompress_string_sequence <- function(sequence.gz, sequence.fa) {
  if (!has_fasta_record_in_gzip(sequence.gz)) {
    return(FALSE)
  }
  unlink(sequence.fa)
  status <- system2("gzip", c("-dc", sequence.gz), stdout = sequence.fa, stderr = FALSE)
  if (!identical(status, 0L) && !identical(status, 0)) {
    unlink(sequence.fa)
    return(FALSE)
  }
  if (!is_nonempty_file(sequence.fa)) {
    unlink(sequence.fa)
    return(FALSE)
  }
  TRUE
}

copy_annotation_files <- function(src_paths, dst_paths) {
  for (name in names(dst_paths)) {
    dir.create(dirname(dst_paths[[name]]), recursive = TRUE, showWarnings = FALSE)
    ok <- file.copy(src_paths[[name]], dst_paths[[name]], overwrite = TRUE)
    if (!ok) {
      stop("Failed to copy annotation file: ", src_paths[[name]], " -> ", dst_paths[[name]])
    }
  }
}

if (is_annotation_complete(db.files)) {
  print(paste("Reuse annotation files from database:", db.annotation.dir))
  copy_annotation_files(db.files, project.outputs)
  quit(save = "no")
}

outdir <- snakemake@params[["outdir"]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setwd(outdir)

input.seqfile <- snakemake@params[["longest_protein"]]
species.list.file <- snakemake@params[["species_list"]]
cpu <- max(1, as.integer(snakemake@threads))
if (is.na(cpu)) {
  cpu <- 1
}

print(paste("Build annotation cache:", db.annotation.dir))
print(species.list.file)

run_command <- function(command, args) {
  print(paste(c(command, args), collapse = " "))
  status <- system2(command, args = args)
  if (!identical(status, 0L) && !identical(status, 0)) {
    stop("Command failed with exit status ", status, ": ", paste(c(command, args), collapse = " "))
  }
}

sequence_cache_file <- function(species.id, ext) {
  file.path(
    db.annotation.dir,
    paste0(species.id, ".protein.sequences.v", stringdb.version, ext)
  )
}

download_string_sequence <- function(species.id, dst.gz) {
  file.sequence <- paste0(
    "https://stringdb-downloads.org/download/protein.sequences.v",
    stringdb.version,
    "/",
    species.id,
    ".protein.sequences.v",
    stringdb.version,
    ".fa.gz"
  )
  tmp.gz <- paste0(dst.gz, ".tmp")
  unlink(tmp.gz)
  run_command("wget", c(
    "--connect-timeout=30",
    "--read-timeout=120",
    "--timeout=120",
    "--tries=10",
    "--waitretry=10",
    "--retry-connrefused",
    "-O", tmp.gz,
    file.sequence
  ))
  if (!is_valid_gzip(tmp.gz)) {
    unlink(tmp.gz)
    stop("Downloaded STRING protein sequence file is not a valid gzip file: ", file.sequence)
  }
  unlink(dst.gz)
  if (!file.rename(tmp.gz, dst.gz)) {
    unlink(tmp.gz)
    stop("Failed to move downloaded STRING protein sequence file to: ", dst.gz)
  }
}

prepare_string_sequence <- function(species.id) {
  db.sequence.gz <- sequence_cache_file(species.id, ".fa.gz")
  sequence.basename <- paste0(species.id, ".protein.sequences.v", stringdb.version)
  work.sequence.gz <- file.path(outdir, paste0(sequence.basename, ".fa.gz"))
  work.sequence.fa <- file.path(outdir, paste0(sequence.basename, ".fa"))

  dir.create(db.annotation.dir, recursive = TRUE, showWarnings = FALSE)

  if (decompress_string_sequence(work.sequence.gz, work.sequence.fa)) {
    return(list(fa = work.sequence.fa, gz = work.sequence.gz))
  }

  if (is_valid_gzip(db.sequence.gz)) {
    if (!file.copy(db.sequence.gz, work.sequence.gz, overwrite = TRUE)) {
      stop("Failed to copy STRING protein sequence cache to working directory: ",
           db.sequence.gz, " -> ", work.sequence.gz)
    }
    if (decompress_string_sequence(work.sequence.gz, work.sequence.fa)) {
      return(list(fa = work.sequence.fa, gz = work.sequence.gz))
    }
  }

  download_string_sequence(species.id, work.sequence.gz)
  if (!decompress_string_sequence(work.sequence.gz, work.sequence.fa)) {
    stop("STRING protein sequence FASTA is empty or invalid after decompression: ", work.sequence.fa)
  }
  list(fa = work.sequence.fa, gz = work.sequence.gz)
}

parse_emapper_annotations <- function(emapper.file) {
  if (!is_nonempty_file(emapper.file)) {
    stop("EggNOG-mapper annotation file is missing or empty: ", file.path(outdir, emapper.file))
  }

  data <- fread(cmd = paste("grep -v '^##'", shQuote(emapper.file)))
  required.columns <- c(
    "#query", "seed_ortholog", "Preferred_name", "Description", "COG_category",
    "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "PFAMs"
  )
  missing.columns <- setdiff(required.columns, colnames(data))
  if (length(missing.columns) > 0) {
    stop("EggNOG-mapper annotation file is missing columns: ", paste(missing.columns, collapse = ", "))
  }

  data %>%
    rename(gene_id = "#query") %>%
    select(gene_id, seed_ortholog, Preferred_name, Description, COG_category,
           GOs, EC, KEGG_ko, KEGG_Pathway, PFAMs)
}

empty.string.mapping <- data.frame(
  gene_id = character(),
  string_id = character(),
  perc_identity = numeric(),
  evalue = numeric(),
  bit_score = numeric()
)

write_empty_string_mapping <- function() {
  fwrite(empty.string.mapping, file = project.outputs[["string_mapping"]], sep = "\t", quote = FALSE)
}

run_string_mapping <- function(species.id) {
  sequence.files <- prepare_string_sequence(species.id)
  sequence.fa <- sequence.files$fa
  sequence.gz <- sequence.files$gz
  diamond.db <- "string_mapping.dmnd"
  diamond.out <- "string_mapping.diamond.tsv"

  run_command("diamond", c("makedb", "--in", sequence.fa, "--db", diamond.db))
  run_command(
    "diamond",
    c(
      "blastp",
      "--query", input.seqfile,
      "--db", diamond.db,
      "--outfmt", "6",
      "--out", diamond.out,
      "--top", "1",
      "--id", "80",
      "--threads", as.character(cpu)
    )
  )

  if (!is_nonempty_file(diamond.out)) {
    warning("DIAMOND did not produce STRING mappings. Writing an empty string_mapping.txt.")
    write_empty_string_mapping()
  } else {
    diamond.result <- fread(diamond.out, header = FALSE) %>%
      rename(
        gene_id = V1,
        string_id = V2,
        perc_identity = V3,
        alignment_length = V4,
        mismatches = V5,
        gap_opens = V6,
        q_start = V7,
        q_end = V8,
        s_start = V9,
        s_end = V10,
        evalue = V11,
        bit_score = V12
      ) %>%
      select(gene_id, string_id, perc_identity, evalue, bit_score) %>%
      arrange(-bit_score) %>%
      distinct(gene_id, .keep_all = TRUE)
    fwrite(diamond.result, file = project.outputs[["string_mapping"]], sep = "\t", quote = FALSE)
  }

  unlink(diamond.out)
  unlink(paste0(diamond.db, "*"))
  unlink(sequence.fa)
  list(species_id = species.id, sequence_gz = sequence.gz)
}

string.species.id <- NA_integer_
prepared.sequence.gz <- NA_character_

############################ EggNOG-mapper annotation ############################
if (is_file_group_complete(db.files[functional.annotation.names])) {
  print(paste("Reuse functional annotation files from database:", db.annotation.dir))
  copy_annotation_files(
    db.files[functional.annotation.names],
    project.outputs[functional.annotation.names]
  )
  data <- fread(project.outputs[["eggnog_mapper_annotation"]])
} else {
  emapper.file <- "diamond_out.emapper.annotations"
  data <- NULL
  if (is_nonempty_file(emapper.file)) {
    print(paste("Reuse existing EggNOG-mapper intermediate:", file.path(outdir, emapper.file)))
    data <- tryCatch(
      parse_emapper_annotations(emapper.file),
      error = function(e) {
        warning("Existing EggNOG-mapper intermediate is invalid and will be regenerated: ", conditionMessage(e))
        NULL
      }
    )
  }

  if (is.null(data)) {
    unlink("diamond_out*")
    run_command(
      "emapper.py",
      c(
        "-i", input.seqfile,
        "--itype", "proteins",
        "-o", "diamond_out",
        "--cpu", as.character(cpu),
        "--outfmt_short",
        "--override"
      )
    )
    data <- parse_emapper_annotations(emapper.file)
  }

  fwrite(data, file = project.outputs[["eggnog_mapper_annotation"]], sep = "\t", quote = FALSE)


  ############################ Functional annotation files ############################
  data.go <- subset(data, GOs != "-", select = c("gene_id", "GOs")) %>%
    tidyr::separate_rows(GOs, sep = ",") %>%
    distinct_all()
  data.kegg <- subset(data, KEGG_Pathway != "-", select = c("gene_id", "KEGG_Pathway")) %>%
    tidyr::separate_rows(KEGG_Pathway, sep = ",") %>%
    distinct_all() %>%
    filter(grepl("ko", KEGG_Pathway))
  data.ko <- subset(data, KEGG_ko != "-", select = c("gene_id", "KEGG_ko")) %>%
    tidyr::separate_rows(KEGG_ko, sep = ",") %>%
    distinct_all()
  data.pfam <- subset(data, PFAMs != "-", select = c("gene_id", "PFAMs")) %>%
    tidyr::separate_rows(PFAMs, sep = ",") %>%
    distinct_all()
  data.cog <- subset(data, COG_category != "-", select = c("gene_id", "COG_category")) %>%
    tidyr::separate_rows(COG_category, sep = "") %>%
    distinct_all() %>%
    filter(COG_category != "")

  write.table(data.go, file = project.outputs[["go_annotation"]], sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.kegg, file = project.outputs[["kegg_pathway_annotation"]], sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.ko, file = project.outputs[["kegg_ko_annotation"]], sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.pfam, file = project.outputs[["pfam_annotation"]], sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.cog, file = project.outputs[["cog_annotation"]], sep = "\t", quote = FALSE, row.names = FALSE)
  copy_annotation_files(project.outputs[functional.annotation.names], db.files[functional.annotation.names])
}


############################ STRING database mapping ############################
species.info <- fread(species.list.file)
species.info[["#taxon_id"]] <- as.character(species.info[["#taxon_id"]])
species.count <- data.frame(
  species.id = sub("\\..*", "", data$seed_ortholog),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(species.id), species.id != "-") %>%
  count(species.id, name = "Freq")

if (configured.tax.id != "") {
  species.id <- configured.tax.id
  if (!species.id %in% species.info[["#taxon_id"]]) {
    stop("Configured tax_id is not present in STRING species list: ", species.id)
  }
} else {
  species.id <- species.count %>%
    merge(species.info, by.x = "species.id", by.y = "#taxon_id") %>%
    filter(Freq > 200, STRING_type %in% c("core", "periphery")) %>%
    arrange(desc(Freq)) %>%
    slice(1) %>%
    pull(species.id) %>%
    as.character() %>%
    as.integer()
}

if (length(species.id) == 0 || is.na(species.id)) {
  warning("No compatible STRING species ID was inferred. Writing an empty string_mapping.txt.")
  write_empty_string_mapping()
} else {
  if (configured.tax.id != "") {
    print(paste("Run STRING mapping with configured tax_id:", species.id))
  } else {
    print(paste("Run STRING mapping with EggNOG-inferred taxid:", species.id))
  }
  string.mapping.result <- run_string_mapping(species.id)
  string.species.id <- string.mapping.result$species_id
  prepared.sequence.gz <- string.mapping.result$sequence_gz
}

dir.create(db.annotation.dir, recursive = TRUE, showWarnings = FALSE)
copy_annotation_files(project.outputs, db.files)
if (!is.na(string.species.id) && !is.na(prepared.sequence.gz) && is_valid_gzip(prepared.sequence.gz)) {
  db.sequence.gz <- sequence_cache_file(string.species.id, ".fa.gz")
  if (!identical(normalizePath(prepared.sequence.gz, mustWork = FALSE),
                 normalizePath(db.sequence.gz, mustWork = FALSE))) {
    ok <- file.copy(prepared.sequence.gz, db.sequence.gz, overwrite = TRUE)
    if (!ok) {
      stop("Failed to copy STRING protein sequence cache: ", prepared.sequence.gz, " -> ", db.sequence.gz)
    }
  }
}
writeLines(
  c(
    paste0("STRINGdb_version=", stringdb.version),
    paste0("tax_id=", configured.tax.id),
    paste0("STRINGdb_dir=", stringdb.dir)
  ),
  db.meta.file
)
print(paste("Annotation files have been cached in database:", db.annotation.dir))
unlink("diamond_out.emapper*")
