#!/usr/bin/env Rscript

# Compare the number of CDS features between each isolate and its 3D7 ortholog.
#
# Usage:
#   Rscript scripts/count_frameshifts.R \
#       --sample-table config/frameshift_samples.tsv \
#       --reference-gff path/to/PlasmoDB-3D7.gff \
#       --output-dir results/frameshifts
#
# The sample table must contain:
# isolate<TAB>orthologs<TAB>gff
#
# Paths in the sample table are relative to DATA_ROOT (environment variable)
# when set, otherwise relative to the repository root.

options(stringsAsFactors = FALSE)

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key)
    }
    key <- sub("^--", "", key)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      stop("Missing value for --", key)
    }
    out[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

if (is.null(args[["sample-table"]]) ||
    is.null(args[["reference-gff"]]) ||
    is.null(args[["output-dir"]])) {
  stop(
    "Usage: Rscript scripts/count_frameshifts.R ",
    "--sample-table config/frameshift_samples.tsv ",
    "--reference-gff path/to/PlasmoDB-3D7.gff ",
    "--output-dir results/frameshifts"
  )
}

# Repository root: scripts/count_frameshifts.R is one directory below it.
script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- script_args[grep("^--file=", script_args)]
if (length(file_arg)) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
  repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
} else {
  repo_root <- normalizePath(".", mustWork = FALSE)
}

data_root <- Sys.getenv("DATA_ROOT", unset = "")
resolve_path <- function(p) {
  if (grepl("^/", p)) {
    stop("Absolute paths are not allowed in the repository configuration: ", p)
  }
  if (nzchar(data_root)) normalizePath(file.path(data_root, p), mustWork = TRUE)
  else normalizePath(file.path(repo_root, p), mustWork = TRUE)
}

reference_gff <- resolve_path(args[["reference-gff"]])
sample_table_path <- resolve_path(args[["sample-table"]])
output_dir <- if (grepl("^/", args[["output-dir"]])) {
  stop("Absolute output paths are not allowed: ", args[["output-dir"]])
} else {
  file.path(repo_root, args[["output-dir"]])
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_gff <- function(path) {
  x <- read.delim(
    path,
    header = FALSE,
    comment.char = "#",
    sep = "\t",
    quote = "",
    fill = TRUE,
    stringsAsFactors = FALSE
  )
  if (ncol(x) < 9L) stop("GFF file has fewer than 9 columns: ", path)
  x
}

extract_id <- function(attributes, prefix) {
  x <- sub(paste0(".*(", prefix, "[^;:, ]*).*"), "\\1", attributes)
  x <- sub(";.*$", "", x)
  x
}

read_orthologs <- function(path, isolate) {
  x <- read.delim(
    path,
    header = FALSE,
    comment.char = "",
    quote = "\"",
    sep = "\t",
    fill = TRUE,
    stringsAsFactors = FALSE
  )
  if (ncol(x) < 2L) stop("Ortholog file needs at least two columns: ", path)
  x <- x[, 1:2]
  names(x) <- c("isolate_gene", "pf3d7_gene")
  x <- x[grepl("PF3D7_", x$pf3d7_gene, fixed = TRUE), , drop = FALSE]
  if (!nrow(x)) warning("No PF3D7 orthologs found for ", isolate)
  x
}

count_features <- function(gff, ids, feature_type) {
  if (length(ids) == 0L) return(setNames(integer(), character()))
  counts <- vapply(
    ids,
    function(id) sum(gff[[3]][grepl(id, gff$ID, fixed = TRUE)] == feature_type),
    integer(1)
  )
  names(counts) <- ids
  counts
}

# Support the original GFF ID conventions while avoiding hard-coded isolate IDs.
prepare_gff <- function(gff, isolate) {
  attrs <- gff[[9]]
  if (grepl("PF3D7", paste(head(attrs, 1000), collapse = " "))) {
    id <- sub(".*?(PF3D7[^;,: ]*).*", "\\1", attrs, perl = TRUE)
    id <- sub(";.*$", "", id)
    id <- sub("-.*$", "", id)
  } else {
    id <- sub(paste0(".*?(", isolate, "[^;,: ]*).*"), "\\1", attrs, perl = TRUE)
    id <- sub(";.*$", "", id)
    id <- sub(":.*$", "", id)
    id <- sub(",.*$", "", id)
  }
  gff$ID <- id
  gff
}

message("Loading reference annotation: ", reference_gff)
gff_3d7 <- prepare_gff(read_gff(reference_gff), "PF3D7")

sample_table <- read.delim(
  sample_table_path,
  header = TRUE,
  comment.char = "#",
  sep = "\t",
  quote = "",
  stringsAsFactors = FALSE
)

required_columns <- c("isolate", "orthologs", "gff")
if (!all(required_columns %in% names(sample_table))) {
  stop(
    "frameshift sample table must contain columns: ",
    paste(required_columns, collapse = ", ")
  )
}

summary_rows <- list()

for (row_i in seq_len(nrow(sample_table))) {
  isolate <- sample_table$isolate[[row_i]]
  ortholog_path <- resolve_path(sample_table$orthologs[[row_i]])
  isolate_gff_path <- resolve_path(sample_table$gff[[row_i]])

  message("Processing ", isolate)

  orthologs <- read_orthologs(ortholog_path, isolate)
  gff_isolate <- prepare_gff(read_gff(isolate_gff_path), isolate)

  ref_counts <- count_features(gff_3d7, orthologs$pf3d7_gene, "CDS")
  isolate_counts <- count_features(gff_isolate, orthologs$isolate_gene, "CDS")

  comparison <- data.frame(
    isolate = isolate,
    isolate_gene = orthologs$isolate_gene,
    pf3d7_gene = orthologs$pf3d7_gene,
    isolate_CDS = as.integer(isolate_counts[orthologs$isolate_gene]),
    pf3d7_CDS = as.integer(ref_counts[orthologs$pf3d7_gene]),
    stringsAsFactors = FALSE
  )
  comparison$CDS_difference <- comparison$isolate_CDS != comparison$pf3d7_CDS

  # Preserve the logic of the original analysis: when the isolate has zero CDS,
  # inspect pseudogenic_exon counts before classifying the difference.
  pseudo <- comparison[comparison$isolate_CDS == 0 & comparison$CDS_difference, , drop = FALSE]

  if (nrow(pseudo)) {
    pseudo_isolate_counts <- vapply(
      pseudo$isolate_gene,
      function(id) sum(gff_isolate[[3]][grepl(id, gff_isolate$ID, fixed = TRUE)] == "pseudogenic_exon"),
      integer(1)
    )
    pseudo_ref_counts <- vapply(
      pseudo$pf3d7_gene,
      function(id) sum(gff_3d7[[3]][grepl(id, gff_3d7$ID, fixed = TRUE)] == "CDS"),
      integer(1)
    )

    pseudo_out <- data.frame(
      isolate = isolate,
      isolate_gene = pseudo$isolate_gene,
      pf3d7_gene = pseudo$pf3d7_gene,
      isolate_pseudogenic_exon = as.integer(pseudo_isolate_counts),
      pf3d7_CDS = as.integer(pseudo_ref_counts),
      pseudogene_difference = pseudo_isolate_counts != pseudo_ref_counts,
      stringsAsFactors = FALSE
    )
  } else {
    pseudo_out <- data.frame(
      isolate = character(),
      isolate_gene = character(),
      pf3d7_gene = character(),
      isolate_pseudogenic_exon = integer(),
      pf3d7_CDS = integer(),
      pseudogene_difference = logical(),
      stringsAsFactors = FALSE
    )
  }

  comparison$core_call <- comparison$CDS_difference
  comparison$core_call[comparison$isolate_CDS == 0] <- FALSE

  nonpseudo <- comparison[comparison$core_call, , drop = FALSE]
  pseudo_frames <- pseudo_out[pseudo_out$pseudogene_difference, , drop = FALSE]

  nonpseudo_n <- nrow(nonpseudo)
  pseudo_n <- nrow(pseudo_frames)
  total_n <- nonpseudo_n + pseudo_n

  utils::write.table(
    comparison,
    file = file.path(output_dir, paste0(isolate, "_frameshift_all.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  utils::write.table(
    nonpseudo,
    file = file.path(output_dir, paste0(isolate, "_frameshift_non_pseudogene.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  utils::write.table(
    pseudo_frames,
    file = file.path(output_dir, paste0(isolate, "_frameshift_pseudogene.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  total_orthologs <- nrow(comparison)
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    isolate = isolate,
    orthologs_evaluated = total_orthologs,
    non_pseudogene_differences = nonpseudo_n,
    pseudogene_differences = pseudo_n,
    total_frameshift_candidates = total_n,
    percent_of_orthologs = if (total_orthologs) 100 * total_n / total_orthologs else NA_real_,
    stringsAsFactors = FALSE
  )
}

summary <- do.call(rbind, summary_rows)
utils::write.table(
  summary,
  file = file.path(output_dir, "frameshift_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

print(summary)
