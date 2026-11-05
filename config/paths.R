# config/paths.R — every location the release code reads or writes comes from an environment
# variable; nothing is hard-coded. Source it from the repository root:  source("config/paths.R")
#
# The R variable names are the ones the scripts already use, so the scripts read unchanged.
# Each script calls psprs_require() with the variables it needs; an unset variable stops the
# script with the name of the environment variable to set.
#
#   Environment variable     R variable(s)                     What it points at
#   PSPRS_DATA_DIR           Khera_path, PGRS_data_path        PRS weights file in; annotation and
#                                                              pathway-membership files out
#   PSPRS_PATHWAYS_DIR       pathways_path                     pathway gene lists out (KEGG, curated)
#   PSPRS_CURATED_CSV        curated_pathways_csv              curated pathway sheet
#                                                              (default: data/OSA_Pathways_20201117_freeze.csv)
#   PSPRS_GENCODE_GFF        path_gencode, GenCode.db          GENCODE GFF3 (GRCh37 lift), e.g.
#                                                              gencode.v34lift37.annotation.gff3.gz
#   PSPRS_GENEHANCER_CSV     path_genehancer, GeneHancer.file  GeneHancer csv (optional; licensed
#                                                              separately; when unset the source is skipped)
#   PSPRS_GTEX_DIR           path_GTEx                         GTEx eQTL directory (*.signif_variant_gene_pairs.txt.gz)
#   PSPRS_FANTOM5_BED        path_FANTOM5, FANTOM5.file        FANTOM5 enhancer_tss_associations.bed
#   PSPRS_PCHIC_DIR          path_PCHiC                        directory holding ActivePromoterEnhancerLinks.tsv
#   PSPRS_KEGG_DIR           path_KEGG                         KGML download cache
#   PSPRS_GENOTYPE_DIR       genotype_dir                      .bgen (+ .bgi) files
#   PSPRS_OUT_DIR            psprs_out_dir                     scoring output root
#   PSPRS_RUN_ID             psprs_run_id                      names the scoring run directory
#   PSPRS_TASK_ID            psprs_task_id                     chunk to score; falls back to
#                                                              SLURM_ARRAY_TASK_ID, then SGE_TASK_ID
#   PSPRS_ANALYSIS_DATA_DIR  UKB_data_path                     analysis data set (UK Biobank, approved access)

.psprs_getenv <- function(name, default = "") {
  v <- Sys.getenv(name, unset = "")
  if (nzchar(v)) v else default
}

.psprs_vars <- c(
  Khera_path = "PSPRS_DATA_DIR", PGRS_data_path = "PSPRS_DATA_DIR",
  pathways_path = "PSPRS_PATHWAYS_DIR", curated_pathways_csv = "PSPRS_CURATED_CSV",
  path_gencode = "PSPRS_GENCODE_GFF", GenCode.db = "PSPRS_GENCODE_GFF",
  path_genehancer = "PSPRS_GENEHANCER_CSV", GeneHancer.file = "PSPRS_GENEHANCER_CSV",
  path_GTEx = "PSPRS_GTEX_DIR",
  path_FANTOM5 = "PSPRS_FANTOM5_BED", FANTOM5.file = "PSPRS_FANTOM5_BED",
  path_PCHiC = "PSPRS_PCHIC_DIR", path_KEGG = "PSPRS_KEGG_DIR",
  genotype_dir = "PSPRS_GENOTYPE_DIR", psprs_out_dir = "PSPRS_OUT_DIR",
  psprs_run_id = "PSPRS_RUN_ID", psprs_task_id = "PSPRS_TASK_ID",
  UKB_data_path = "PSPRS_ANALYSIS_DATA_DIR")

Khera_path           <- .psprs_getenv("PSPRS_DATA_DIR")
PGRS_data_path       <- Khera_path
pathways_path        <- .psprs_getenv("PSPRS_PATHWAYS_DIR")
curated_pathways_csv <- .psprs_getenv("PSPRS_CURATED_CSV", file.path("data", "OSA_Pathways_20201117_freeze.csv"))

# Source files are given as files; the scripts join a directory and a file name, so both are set.
.gencode_gff    <- .psprs_getenv("PSPRS_GENCODE_GFF")
path_gencode    <- if (nzchar(.gencode_gff)) dirname(.gencode_gff) else ""
GenCode.db      <- if (nzchar(.gencode_gff)) basename(.gencode_gff) else ""
.genehancer_csv <- .psprs_getenv("PSPRS_GENEHANCER_CSV")
path_genehancer <- if (nzchar(.genehancer_csv)) dirname(.genehancer_csv) else ""
GeneHancer.file <- if (nzchar(.genehancer_csv)) basename(.genehancer_csv) else ""
.fantom5_bed    <- .psprs_getenv("PSPRS_FANTOM5_BED")
path_FANTOM5    <- if (nzchar(.fantom5_bed)) dirname(.fantom5_bed) else ""
FANTOM5.file    <- if (nzchar(.fantom5_bed)) basename(.fantom5_bed) else ""
path_GTEx       <- .psprs_getenv("PSPRS_GTEX_DIR")
path_PCHiC      <- .psprs_getenv("PSPRS_PCHIC_DIR")
path_KEGG       <- .psprs_getenv("PSPRS_KEGG_DIR")

genotype_dir  <- .psprs_getenv("PSPRS_GENOTYPE_DIR")
psprs_out_dir <- .psprs_getenv("PSPRS_OUT_DIR")
psprs_run_id  <- .psprs_getenv("PSPRS_RUN_ID")
psprs_task_id <- .psprs_getenv("PSPRS_TASK_ID",
                   .psprs_getenv("SLURM_ARRAY_TASK_ID", .psprs_getenv("SGE_TASK_ID")))

UKB_data_path <- .psprs_getenv("PSPRS_ANALYSIS_DATA_DIR")

# Stop unless every named R variable is set; the message names the environment variable.
psprs_require <- function(vars) {
  unknown <- setdiff(vars, names(.psprs_vars))
  if (length(unknown)) stop("psprs_require: unknown variable(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  caller <- parent.frame()
  empty <- vars[!vapply(vars, function(v) nzchar(get(v, envir = caller)), logical(1))]
  if (length(empty))
    stop("Set the environment variable(s): ", paste(unique(.psprs_vars[empty]), collapse = ", "),
         "  (needed for: ", paste(empty, collapse = ", "), ")", call. = FALSE)
  invisible(TRUE)
}
