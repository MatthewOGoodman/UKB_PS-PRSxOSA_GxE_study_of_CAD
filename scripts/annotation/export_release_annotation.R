# scripts/annotation/export_release_annotation.R — make the distributable copy of the per-source
# variant-to-gene annotation: the combined annotation list without the GeneHancer sources
# (licensed separately; not redistributable) and without the PRS weight column (the weights are
# published separately, e.g. PGS Catalog PGS000013).
#
# Run from the repository root:  Rscript scripts/annotation/export_release_annotation.R
# Required environment: PSPRS_DATA_DIR. Optional: PSPRS_ANNOTATION_DATE (default "20201209",
#   the annotation used for the paper), PSPRS_PRS_NAME, PSPRS_EXPORT_DIR (default
#   <PSPRS_DATA_DIR>/release).
# Output: <PSPRS_EXPORT_DIR>/PSPRS_annotation_GRCh37_<date>.rds and its .md5

suppressPackageStartupMessages({ library(GenomicRanges) })
source("config/paths.R")
psprs_require("PGRS_data_path")

annotation_date = Sys.getenv("PSPRS_ANNOTATION_DATE", "20201209")
PGRS_suffix = paste(Sys.getenv("PSPRS_PRS_NAME", "PGRS_CAD_Khera_LDPred_rho0.001"), "rds", sep = ".")

combined_annot_list = readRDS(file.path(PGRS_data_path,
  paste("gene_annot", annotation_date, "combined_annotation_sources", PGRS_suffix, sep = "_")))

release_list = combined_annot_list[!grepl("genehancer", names(combined_annot_list), ignore.case = TRUE)]
release_list = lapply(release_list, function(g) { mcols(g)$wgt = NULL; g })

stopifnot(
  "a GeneHancer source remains" = !any(grepl("genehancer", unlist(lapply(release_list, function(g) unique(g$source))), ignore.case = TRUE)),
  "a weight column remains"     = !any(vapply(release_list, function(g) "wgt" %in% names(mcols(g)), logical(1))))

out_dir = Sys.getenv("PSPRS_EXPORT_DIR", file.path(PGRS_data_path, "release"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file = file.path(out_dir, paste0("PSPRS_annotation_GRCh37_", annotation_date, ".rds"))
saveRDS(release_list, out_file)
writeLines(paste(unname(tools::md5sum(out_file)), basename(out_file)), paste0(out_file, ".md5"))
cat("wrote", out_file, "\nsources:", paste(names(release_list), collapse = ", "), "\n")
