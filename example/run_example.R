# example/run_example.R — score the toy PRS (example/toy_weights.tsv) on the bgen project's example
# genotypes with the same two scripts used for a full run: scripts/scoring/score_chunk.R, then
# scripts/scoring/combine_chunks.R.
#
# Run from the repository root, after setup/build_rbgen.sh and setup/setup_renv.R:
#   Rscript example/run_example.R
# Output: tests/_output/example/ (the toy membership file and the run directory)

suppressPackageStartupMessages({ library(GenomicRanges); library(S4Vectors) })

example_out <- file.path("tests", "_output", "example")
data_dir    <- file.path(example_out, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# 1. the toy PRS as a pathway-annotated GRanges (the shape the annotation and pathway builders
#    produce for a real PRS): one row per variant, wgt + alleles + logical pathway columns
toy <- read.delim(file.path("example", "toy_weights.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
membership <- GRanges(seqnames = toy$chr, ranges = IRanges(start = toy$pos, end = toy$pos))
mcols(membership) <- DataFrame(toy[, setdiff(colnames(toy), c("chr", "pos"))], check.names = FALSE)
saveRDS(membership, file.path(data_dir, "toy_membership.rds"))

# 2. the environment a scheduler task would receive
Sys.setenv(
  PSPRS_DATA_DIR              = data_dir,
  PSPRS_MEMBERSHIP_FILE       = "toy_membership.rds",
  PSPRS_GENOTYPE_DIR          = file.path("example", "bgen"),
  PSPRS_GENOTYPE_FILE_PATTERN = "example.16bits.bgen",   # one file; no "{chr}" placeholder needed
  PSPRS_CHR_LABEL             = "zero_pad",              # the example uses UK Biobank's "01"
  PSPRS_BLOCKSIZE             = "50",
  PSPRS_OUT_DIR               = file.path(example_out, "runs"),
  PSPRS_RUN_ID                = "example",
  PSPRS_TASK_ID               = "1")                      # the toy PRS is a single chunk

# 3. score the chunk, then combine (one chunk here; a real run has one task per chunk)
source(file.path("scripts", "scoring", "score_chunk.R"), local = new.env())
source(file.path("scripts", "scoring", "combine_chunks.R"), local = new.env())

scores <- readRDS(file.path(example_out, "runs", "example", "PathSpecPGRS_total_scores.rds"))
cat("\nPS-PRS for", nrow(scores), "example samples x", ncol(scores), "pathways\n")
print(round(head(scores), 4))
