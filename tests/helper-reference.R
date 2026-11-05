# tests/helper-reference.R — shared fixtures for the scorer tests (sourced by testthat::test_dir).
#
# Run the tests from the repository root:
#   Rscript -e 'testthat::test_dir("tests", stop_on_failure = FALSE)'
# PSPRS_RELEASE_R_DIR (default "R") selects the function files under test.

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(rbgen)
  library(dplyr); library(tibble); library(stringr); library(purrr)
})

example_dir <- normalizePath("../example", mustWork = TRUE)
bgen_dir    <- file.path(example_dir, "bgen")
release_rel <- Sys.getenv("PSPRS_RELEASE_R_DIR", "R")
release_dir <- normalizePath(if (grepl("^/", release_rel)) release_rel else file.path("..", release_rel), mustWork = TRUE)

pathway_cols_all  <- c("all_SNPs", "toy_pathway_A", "toy_pathway_B", "hsa99999: Toy/pathway, with punctuation")
pathway_cols_safe <- pathway_cols_all[1:3]

read_toy <- function() {
  read.delim(file.path(example_dir, "toy_weights.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
}

# toy table -> the GRanges the scorer takes (one chunk, seqnames "1", pathway columns logical)
toy_granges <- function(toy = read_toy(), seqlevels_order = NULL) {
  gr <- GRanges(seqnames = toy$chr, ranges = IRanges(start = toy$pos, end = toy$pos))
  mc <- S4Vectors::DataFrame(toy[, setdiff(colnames(toy), c("chr", "pos"))], check.names = FALSE)
  mcols(gr) <- mc
  if (!is.null(seqlevels_order)) seqlevels(gr) <- seqlevels_order
  gr
}

# source the release scorer into its own environment
release_env <- function() {
  e <- new.env(parent = globalenv())
  sys.source(file.path(release_dir, "genotype_dosage_fns.R"), envir = e)
  sys.source(file.path(release_dir, "psprs_score_bgen.R"), envir = e)
  e
}

# run the release scorer on a toy GRanges; returns the chunk summary matrix
run_release <- function(gr, out_dir, run_id = "test", blocksize = 50, effect_allele_col = "effect_allele",
                        pathway_columns = pathway_cols_all, env = release_env(), chunk_name = "1") {
  invisible(capture.output(
    path <- env$compute_PS_PGRS(PGRS_chunk = gr, PGRS_chunk_name = chunk_name, pathway_columns = pathway_columns,
                                genotype_dir = bgen_dir, genotype_filename_fn = function(chr) "example.16bits.bgen",
                                out_dir = out_dir, run_id = run_id, blocksize = blocksize,
                                effect_allele_col = effect_allele_col, verbose_TF = FALSE)))
  readRDS(path)
}

# collect warnings while evaluating expr; returns list(value, warnings)
with_warnings <- function(expr) {
  w <- character(0)
  value <- withCallingHandlers(expr, warning = function(cond) { w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning") })
  list(value = value, warnings = w)
}

# independent reference from example.gen: sum over PGRS rows and every bgen variant at the row's
# position whose allele set equals {A1, A2}, of membership x weight x effect-allele dosage.
# A missing genotype (all-zero triplet) makes that sample's score NA.
reference_scores <- function(toy = read_toy(), pathway_columns = pathway_cols_all) {
  g <- read.table(file.path(bgen_dir, "example.gen"), stringsAsFactors = FALSE)
  samples <- read.table(file.path(bgen_dir, "example.sample"), skip = 2, stringsAsFactors = FALSE)$V1
  p <- as.matrix(g[, -(1:6)])
  p0 <- p[, seq(1, ncol(p), 3)]; p1 <- p[, seq(2, ncol(p), 3)]; p2 <- p[, seq(3, ncol(p), 3)]
  d1 <- p1 + 2 * p2                             # dosage of allele1
  d1[(p0 + p1 + p2) == 0] <- NA
  out <- matrix(0, nrow = ncol(d1), ncol = length(pathway_columns),
                dimnames = list(samples, pathway_columns))
  for (r in seq_len(nrow(toy))) {
    j <- which(g$V4 == toy$pos[r] & ((g$V5 == toy$A1[r] & g$V6 == toy$A2[r]) | (g$V5 == toy$A2[r] & g$V6 == toy$A1[r])))
    for (jj in j) {
      d_ea <- if (toy$effect_allele[r] == g$V6[jj]) d1[jj, ] else 2 - d1[jj, ]
      m <- as.numeric(unlist(toy[r, pathway_columns]))
      out <- out + outer(d_ea, m * toy$wgt[r])
    }
  }
  out
}
