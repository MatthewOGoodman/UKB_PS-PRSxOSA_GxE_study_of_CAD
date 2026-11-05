# scripts/annotation/run_annotation.R — annotate every PRS variant to genes from each source, then
# assign genes by the source-priority cascade.
#
# Run from the repository root:  Rscript scripts/annotation/run_annotation.R
#
# Required environment: PSPRS_DATA_DIR (holds the PRS weights file; outputs are written here),
#   PSPRS_GENCODE_GFF, PSPRS_GTEX_DIR, PSPRS_FANTOM5_BED, PSPRS_PCHIC_DIR.
# Optional: PSPRS_GENEHANCER_CSV (licensed separately; when unset the GeneHancer source is skipped),
#   PSPRS_PRS_FILE, PSPRS_PRS_FORMAT (khera | pgs_catalog), PSPRS_PRS_NAME, PSPRS_CHROMOSOMES
#   (comma-separated, e.g. "22", to annotate a subset), PSPRS_ANNOTATION_DATE (names the outputs;
#   default today).
# Outputs in PSPRS_DATA_DIR: gene_annot_<date>_<source>_<PRS name>.rds per source,
#   gene_annot_<date>_combined_annotation_sources_<PRS name>.rds, and
#   gene_annot_<date>_sequential.annotated.SNPs.to.ENSGs_<PRS name>.rds.

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(rtracklayer); library(GenomicFeatures)
  library(InteractionSet); library(data.table); library(dplyr); library(tidyr); library(readr)
  library(tibble); library(stringr); library(purrr); library(ggplot2)
})

source("config/paths.R")
source("R/genomic_utils.R")
source("R/load_prs_weights.R")
source("R/annotate_genes_fns.R")
psprs_require("PGRS_data_path")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PGRS options ----

# trained PGRS with GWAS variants (the paper: Khera et al. 2018 GPS_CAD; PGS Catalog PGS000013)
PGRS.trained_file = Sys.getenv("PSPRS_PRS_FILE", "CoronaryArteryDisease_PRS_LDpred_rho0.001_v3.txt")
PGRS_format = Sys.getenv("PSPRS_PRS_FORMAT", "khera")
PGRS_colnames_key = switch(PGRS_format,
  khera       = c(variant = "variant", chr = "chr", pos = "position_hg19", effect_weight = "effect_weight",
                  effect_allele = "effect_allele", A1 = "A1", A2 = "A2"),
  pgs_catalog = c(chr = "chr_name", pos = "chr_position", effect_weight = "effect_weight",
                  effect_allele = "effect_allele", A2 = "other_allele"),
  stop("PSPRS_PRS_FORMAT must be khera or pgs_catalog"))
# header lines before the column names (both formats carry a commented header)
PGRS_skip_lines = sum(grepl("^#", readLines(file.path(PGRS_data_path, PGRS.trained_file), n = 100)))


# annotations output:
###
PGRS_out_name = Sys.getenv("PSPRS_PRS_NAME", "PGRS_CAD_Khera_LDPred_rho0.001")
PGRS_suffix = paste(PGRS_out_name, "rds", sep = ".")

PGRS_variants_GR = load_split_save_PGRS(
  file_path = file.path(PGRS_data_path, PGRS.trained_file),
  out_name = PGRS_out_name,
  out_folder = PGRS_data_path,
  colnames_key = PGRS_colnames_key,
  skip_lines = PGRS_skip_lines,
  split_save_TF = FALSE,
  verbose_TF = FALSE,
  return_type = "GRanges")$GRanges

PGRS_chromosomes = Sys.getenv("PSPRS_CHROMOSOMES", "")
if(nzchar(PGRS_chromosomes)){
  keep_chr = paste0("chr", sub("^chr", "", trimws(strsplit(PGRS_chromosomes, ",")[[1]])))
  PGRS_variants_GR = PGRS_variants_GR[as.character(seqnames(PGRS_variants_GR)) %in% keep_chr]
  PGRS_variants_GR = keepSeqlevels(PGRS_variants_GR, unique(as.character(seqnames(PGRS_variants_GR))))
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# option whether to run annotation code from scratch (PSPRS_ANNOTATION_STEPS):
  # "all"     run annotate_variants_by_source.R to get raw annotations, then combine and cascade
  # "combine" re-combine existing per-source files (e.g. new GeneHancer settings), then cascade
  # "cascade" only the cascade, from an existing combined annotation file (e.g. the archived one)
annotation_steps = Sys.getenv("PSPRS_ANNOTATION_STEPS", "all")
if(!(annotation_steps %in% c("all", "combine", "cascade"))) stop("PSPRS_ANNOTATION_STEPS must be all, combine or cascade")
annotate_from_scratch_TF = annotation_steps == "all"
annotation_date_fixed = Sys.getenv("PSPRS_ANNOTATION_DATE") # if !annotate_from_scratch_TF, roll back to annotation files produced on this date

  if(annotate_from_scratch_TF){
    annotation_date = Sys.getenv("PSPRS_ANNOTATION_DATE", format(Sys.time(), "%Y%m%d"))
  } else {
    if(!nzchar(annotation_date_fixed)) stop("Set PSPRS_ANNOTATION_DATE to the date naming the existing annotation files")
    annotation_date = annotation_date_fixed
  }

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ANNOTATE GENES ----

if(annotate_from_scratch_TF){
  psprs_require(c("path_gencode", "path_GTEx", "path_FANTOM5", "path_PCHiC"))
  source( file.path("scripts", "annotation", "annotate_variants_by_source.R"))
  source( file.path("scripts", "annotation", "combine_annotation_sources.R"))
  source( file.path("scripts", "annotation", "assign_genes_cascade.R"))
} else if(annotation_steps == "combine"){
  source( file.path("scripts", "annotation", "combine_annotation_sources.R"))
  source( file.path("scripts", "annotation", "assign_genes_cascade.R"))
} else {
  source( file.path("scripts", "annotation", "assign_genes_cascade.R"))
}
