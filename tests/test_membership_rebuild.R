# T5 and T7 — the cascade and the pathway membership against the paper's 2020 files, and the
# pathway builder on the current KEGG release.
#   T5 needs the private 2020 files: skipped unless PSPRS_PRIVATE_DATA_DIR is set.
#   T7 needs network access to KEGG: skipped unless PSPRS_NETWORK_TESTS=1 (the driver test also
#   needs PSPRS_PRIVATE_DATA_DIR and PSPRS_GENCODE_GFF).

suppressPackageStartupMessages(library(ggplot2))   # the cascade script plots its cumulative summary

private_dir <- Sys.getenv("PSPRS_PRIVATE_DATA_DIR")
network    <- identical(Sys.getenv("PSPRS_NETWORK_TESTS"), "1")
prs_suffix  <- "PGRS_CAD_Khera_LDPred_rho0.001.rds"
sources_tag <- "GENCODE_exons.FANTOM5.PCHiC.GTEx.genehancer.GENCODE_UTR.GENCODE_flanked.GENCODE_noncoding"
paper_runs  <- list(
  HIF1_VEGF = c(genes = "psPGRS_20201209_pathway_genes.rds",
                paths = paste0("psPGRS_20201209_PGRS_CAD_Khera_", sources_tag, "_annotated_PS.PGRS.rds")),
  NFkB_TNF  = c(genes = "psPGRS_20201209_PGRS_CAD_Khera_NFkB_TNF_pathway_genes.rds",
                paths = paste0("psPGRS_20201209_PGRS_CAD_Khera_NFkB_TNF_", sources_tag, "_annotated_PS.PGRS.rds")))

# evaluate only the annotate_PGRS_pathways definition from the pathway script
pathway_fn_env <- function() {
  e <- new.env(parent = globalenv())
  for (x in parse(file.path("..", "scripts", "pathways", "annotate_pathways.R")))
    if (is.call(x) && identical(x[[1]], as.name("=")) && identical(x[[2]], as.name("annotate_PGRS_pathways")))
      eval(x, e)
  stopifnot(exists("annotate_PGRS_pathways", envir = e, inherits = FALSE))
  e
}

test_that("the cascade reproduces the 2020 gene assignment", {
  skip_if(!nzchar(private_dir), "PSPRS_PRIVATE_DATA_DIR not set")
  e <- new.env(parent = globalenv())
  e$combined_annot_list <- readRDS(file.path(private_dir, paste("gene_annot", "20201209", "combined_annotation_sources", prs_suffix, sep = "_")))
  e$PGRS_variants_GR <- GRanges("chr1", IRanges(1, 1), variant = "scaffold")
  e$PGRS_data_path <- tempfile("cascade_"); dir.create(e$PGRS_data_path)
  e$annotation_date <- "test"; e$PGRS_out_name <- sub("[.]rds$", "", prs_suffix)
  invisible(capture.output(suppressMessages(
    sys.source(file.path("..", "scripts", "annotation", "assign_genes_cascade.R"), envir = e))))
  new <- readRDS(file.path(e$PGRS_data_path, paste("gene_annot", "test", "sequential.annotated.SNPs.to.ENSGs", prs_suffix, sep = "_")))
  old <- readRDS(file.path(private_dir, paste("gene_annot", "20201209", "sequential.annotated.SNPs.to.ENSGs", prs_suffix, sep = "_")))
  # compared by content: the 2020 object was serialised by older Bioconductor classes
  # (DataFrame, GenomeInfoDb::Seqinfo) that current versions name differently
  expect_identical(as.data.frame(unname(new)), as.data.frame(unname(old)))
  expect_identical(names(new), names(old))
})

test_that("pathway membership reproduces both 2020 pathway-annotated files", {
  skip_if(!nzchar(private_dir), "PSPRS_PRIVATE_DATA_DIR not set")
  e <- pathway_fn_env()
  seq_ann <- readRDS(file.path(private_dir, paste("gene_annot", "20201209", "sequential.annotated.SNPs.to.ENSGs", prs_suffix, sep = "_")))
  for (run in names(paper_runs)) {
    gene_lists <- readRDS(file.path(private_dir, paper_runs[[run]][["genes"]]))
    old <- readRDS(file.path(private_dir, paper_runs[[run]][["paths"]]))
    path_cols <- c(names(gene_lists), "all_SNPs")
    scaffold <- old
    mcols(scaffold) <- mcols(scaffold)[, setdiff(names(mcols(scaffold)), path_cols), drop = FALSE]
    res <- NULL
    invisible(capture.output(
      res <- e$annotate_PGRS_pathways(pathways_list = lapply(gene_lists, names), PGRS_multi_annotated = seq_ann,
                                      PGRS = scaffold, verbose_TF = FALSE)))
    new <- res$PGRS_paths
    new$all_SNPs <- TRUE
    for (p in path_cols) expect_identical(mcols(new)[[p]], mcols(old)[[p]], label = paste(run, p))
  }
})

test_that("the KEGG builder fetches a single pathway from the current release", {
  skip_if(!network, "PSPRS_NETWORK_TESTS not set")
  suppressPackageStartupMessages({ library(KEGGgraph); library(KEGGREST); library(org.Hs.eg.db) })
  e <- new.env(parent = globalenv())
  e$path_KEGG <- tempfile("kegg_")
  sys.source(file.path("..", "R", "kegg_fns.R"), envir = e)
  res <- NULL
  invisible(capture.output(res <- e$get_genes_in_KEGGREST_pathway("hsa04066")))
  expect_length(res$ENSG, 1)
  expect_match(names(res$ENSG), "^hsa04066: ")
  expect_gt(length(res$ENSG[[1]]), 50)
  expect_true(any(grepl("^KEGG_pathway_[0-9]{8}$", list.files(e$path_KEGG))))
})

test_that("the pathway driver builds membership for listed pathways and gene sets", {
  skip_if(!network || !nzchar(private_dir) || !nzchar(Sys.getenv("PSPRS_GENCODE_GFF")), "network, private data or GENCODE not set")
  work <- tempfile("paths_"); dir.create(work)
  for (f in c("CoronaryArteryDisease_PRS_LDpred_rho0.001_v3.txt",
              paste("gene_annot", "20201209", "sequential.annotated.SNPs.to.ENSGs", prs_suffix, sep = "_"),
              "gene_annot_20201209_all_ENSGs_annotated.txt"))
    file.symlink(file.path(private_dir, f), file.path(work, f))
  out <- withr::with_dir("..", withr::with_envvar(
    c(PSPRS_DATA_DIR = work, PSPRS_PATHWAYS_DIR = file.path(work, "pathways"), PSPRS_KEGG_DIR = file.path(work, "kgml"),
      PSPRS_ANNOTATION_DATE = "20201209", PSPRS_KEGG_PATHWAYS = "hsa04066,hsa04370", PSPRS_GENE_SETS = "HIF1_VEGF"),
    system2(file.path(R.home("bin"), "Rscript"), "scripts/pathways/run_pathways.R", stdout = TRUE, stderr = TRUE)))
  status <- attr(out, "status")
  expect_true(is.null(status) || status == 0, label = paste(tail(out, 20), collapse = "\n"))
  f <- list.files(work, pattern = "annotated_PS[.]PGRS[.]rds$", full.names = TRUE)
  expect_length(f, 1)
  m <- readRDS(f)
  cols <- names(mcols(m))
  expect_true(all(c("all_annotated_genes", "all_SNPs") %in% cols))
  expect_equal(sum(grepl("^hsa04066: |^hsa04370: ", cols)), 2)
  expect_true(any(grepl("^HIF1_only_genes[.]", cols)) && any(grepl("^HIF1_VEGF_genes[.]", cols)))
})
