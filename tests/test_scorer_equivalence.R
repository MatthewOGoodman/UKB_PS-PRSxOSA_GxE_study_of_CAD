# T3 — the release scorer against the original scorer the paper ran, on the same example input.
# Needs the original code; skipped unless PSPRS_ORIGINAL_CODE_DIR names the directory holding
# 2.0_fns_PathSpec_PGRS_bgen.R and 2.0_fns_loadGenotypes.R.
#
# The original counted allele A2 of each PRS row. With effect_allele_col = "A2" the release
# must reproduce it exactly; with the default (effect allele = A1 on every row, as in the paper's
# weights) each matched row contributes 2*w - (original contribution), so
#   release = 2 * sum(w over matched rows, with multiplicity) - original.

orig_dir <- Sys.getenv("PSPRS_ORIGINAL_CODE_DIR")

original_env <- function() {
  e <- new.env(parent = globalenv())
  sys.source(file.path(orig_dir, "2.0_fns_loadGenotypes.R"), envir = e)
  sys.source(file.path(orig_dir, "2.0_fns_PathSpec_PGRS_bgen.R"), envir = e)
  e
}

run_original <- function(gr, blocksize = 50, pathway_columns = pathway_cols_safe) {
  work <- tempfile("orig_")
  dir.create(file.path(work, "tmp"), recursive = TRUE)
  old <- Sys.getenv(c("JOB_NAME", "JOB_ID", "HOSTNAME", "TMPDIR"), unset = NA)
  on.exit({
    for (v in names(old)) if (is.na(old[[v]])) Sys.unsetenv(v) else do.call(Sys.setenv, setNames(list(old[[v]]), v))
  }, add = TRUE)
  Sys.setenv(JOB_NAME = "eqv", JOB_ID = "1", HOSTNAME = "local", TMPDIR = file.path(work, "tmp"))
  # the original renames chromosome levels by position, so it needs all 22 autosomal levels in order
  seqlevels(gr) <- c(as.character(1:22))
  e <- original_env()
  invisible(capture.output(suppressWarnings(
    e$compute_PS_PGRS(PGRS_chunk = gr, PGRS_chunk_name = "1", pathway_columns = pathway_columns,
                      genotype_dir = bgen_dir, genotype_filename_fn = function(chr) "example.16bits.bgen",
                      out_dir = work, temp_dir = file.path(work, "tmp"), blocksize = blocksize,
                      subset_samples = NULL, existing_iters = NULL, verbose_TF = FALSE))))
  readRDS(file.path(work, "eqv_1", "PathSpecPGRS_chunk_summary_1.rds"))
}

# effect allele = A1 on every row, as in the paper's weights file
toy_eqv <- read_toy()
toy_eqv$effect_allele <- toy_eqv$A1

test_that("release with effect_allele_col = 'A2' is identical to the original", {
  skip_if(!nzchar(orig_dir), "PSPRS_ORIGINAL_CODE_DIR not set")
  gr <- toy_granges(toy_eqv)
  orig <- run_original(gr)
  rel  <- suppressWarnings(run_release(gr, file.path(tempfile(), "out"), effect_allele_col = "A2",
                                       pathway_columns = pathway_cols_safe))
  expect_identical(rel, orig)
})

test_that("release (effect allele) = 2 * matched weight - original", {
  skip_if(!nzchar(orig_dir), "PSPRS_ORIGINAL_CODE_DIR not set")
  gr <- toy_granges(toy_eqv)
  orig <- run_original(gr)
  rel  <- suppressWarnings(run_release(gr, file.path(tempfile(), "out"), pathway_columns = pathway_cols_safe))
  # matched weight per pathway: rows whose allele set equals a bgen variant's at that position,
  # counted once per matching bgen variant
  g <- read.table(file.path(bgen_dir, "example.gen"), stringsAsFactors = FALSE)[, 1:6]
  mult <- vapply(seq_len(nrow(toy_eqv)), function(r)
    sum(g$V4 == toy_eqv$pos[r] & ((g$V5 == toy_eqv$A1[r] & g$V6 == toy_eqv$A2[r]) |
                                  (g$V5 == toy_eqv$A2[r] & g$V6 == toy_eqv$A1[r]))), numeric(1))
  w_matched <- colSums(as.matrix(toy_eqv[, pathway_cols_safe]) * toy_eqv$wgt * mult)
  expected <- sweep(-orig, 2, 2 * w_matched, "+")
  expect_identical(is.na(rel), is.na(orig))
  expect_lt(max(abs(rel - expected), na.rm = TRUE), 1e-12)
})
