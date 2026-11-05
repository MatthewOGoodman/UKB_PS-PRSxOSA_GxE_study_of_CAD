# T8 and T6 — the distributable annotation export, and the annotation builder on one chromosome,
# against the paper's 2020 annotation files. Need the private 2020 files and the annotation
# sources; skipped unless PSPRS_PRIVATE_DATA_DIR (the directory holding the 2020 files and the
# PRS weights file) is set. T6 also needs PSPRS_GENCODE_GFF, PSPRS_GTEX_DIR, PSPRS_FANTOM5_BED,
# PSPRS_PCHIC_DIR and PSPRS_GENEHANCER_CSV, and network access (AnnotationHub liftover chain).

private_dir <- Sys.getenv("PSPRS_PRIVATE_DATA_DIR")
prs_name    <- "PGRS_CAD_Khera_LDPred_rho0.001.rds"
annot_file  <- function(dir, date, tag) file.path(dir, paste("gene_annot", date, tag, prs_name, sep = "_"))
rscript     <- file.path(R.home("bin"), "Rscript")

# run a repository script from the repository root in a child R, with extra environment
run_script <- function(script, env) {
  out <- withr::with_dir("..", withr::with_envvar(env,
    system2(rscript, script, stdout = TRUE, stderr = TRUE)))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) stop("script failed:\n", paste(tail(out, 30), collapse = "\n"))
  invisible(out)
}

pairs_of <- function(g) unique(paste(g$variant, g$ENSG, g$source))

test_that("the export drops GeneHancer and the weights and keeps every other row", {
  skip_if(!nzchar(private_dir), "PSPRS_PRIVATE_DATA_DIR not set")
  out_dir <- tempfile("export_")
  run_script("scripts/annotation/export_release_annotation.R",
             c(PSPRS_DATA_DIR = private_dir, PSPRS_EXPORT_DIR = out_dir))
  rel  <- readRDS(file.path(out_dir, "PSPRS_annotation_GRCh37_20201209.rds"))
  orig <- readRDS(annot_file(private_dir, "20201209", "combined_annotation_sources"))
  expect_length(rel, 7)
  expect_false(any(grepl("genehancer", c(names(rel), unlist(lapply(rel, function(g) unique(g$source)))), ignore.case = TRUE)))
  expect_false(any(vapply(rel, function(g) "wgt" %in% names(mcols(g)), logical(1))))
  for (nm in names(rel)) {
    o <- orig[[nm]]; mcols(o)$wgt <- NULL
    expect_identical(rel[[nm]], o, label = nm)
  }
})

test_that("the annotation builder reproduces the 2020 annotation on chromosome 22", {
  skip_if(!nzchar(private_dir), "PSPRS_PRIVATE_DATA_DIR not set")
  skip_if(!nzchar(Sys.getenv("PSPRS_GENEHANCER_CSV")), "annotation sources not set")
  work <- tempfile("annot_")
  dir.create(work)
  file.symlink(file.path(private_dir, "CoronaryArteryDisease_PRS_LDpred_rho0.001_v3.txt"),
               file.path(work, "CoronaryArteryDisease_PRS_LDpred_rho0.001_v3.txt"))
  run_script("scripts/annotation/run_annotation.R",
             c(PSPRS_DATA_DIR = work, PSPRS_CHROMOSOMES = "22", PSPRS_ANNOTATION_DATE = "test"))

  new_list <- readRDS(annot_file(work, "test", "combined_annotation_sources"))
  old_list <- readRDS(annot_file(private_dir, "20201209", "combined_annotation_sources"))
  expect_setequal(names(new_list), names(old_list))
  for (nm in names(old_list)) {
    o <- old_list[[nm]]; o <- o[as.character(seqnames(o)) == "chr22"]
    expect_setequal(pairs_of(new_list[[nm]]), pairs_of(o))
  }

  new_seq <- readRDS(annot_file(work, "test", "sequential.annotated.SNPs.to.ENSGs"))
  old_seq <- readRDS(annot_file(private_dir, "20201209", "sequential.annotated.SNPs.to.ENSGs"))
  old_seq <- old_seq[as.character(seqnames(old_seq)) == "chr22"]
  expect_setequal(pairs_of(new_seq), pairs_of(old_seq))
})
