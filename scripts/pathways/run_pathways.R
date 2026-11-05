# scripts/pathways/run_pathways.R — fetch pathway gene lists and mark, for every PRS variant, the
# pathways whose genes it is assigned to. The output is the pathway-annotated PRS that
# scripts/scoring/score_chunk.R scores.
#
# Run from the repository root, after scripts/annotation/run_annotation.R:
#   Rscript scripts/pathways/run_pathways.R
#
# Required environment: PSPRS_DATA_DIR, PSPRS_PATHWAYS_DIR, PSPRS_KEGG_DIR, PSPRS_GENCODE_GFF,
#   PSPRS_ANNOTATION_DATE (the date that names the annotation outputs of run_annotation.R).
# Optional: PSPRS_KEGG_PATHWAYS (curated [default] | all | hsa04066,hsa04370,...),
#   PSPRS_GENE_SETS (single-gene sets: NFkB_TNF [default], HIF1_VEGF, or both comma-separated;
#   none to skip), PSPRS_CURATED_CSV, PSPRS_PRS_FILE, PSPRS_PRS_FORMAT, PSPRS_PRS_NAME.
# Outputs in PSPRS_DATA_DIR: psPGRS_<date>_<PRS>_<gene sets>_<sources>_annotated_PS.PGRS.rds and
#   psPGRS_<date>_<PRS>_<gene sets>_pathway_genes.rds; pathway gene lists in PSPRS_PATHWAYS_DIR.

suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomeInfoDb); library(rtracklayer); library(KEGGgraph)
  library(KEGGREST); library(org.Hs.eg.db); library(dplyr); library(tidyr); library(readr)
  library(tibble); library(stringr); library(purrr)
})

source("config/paths.R")
source("R/genomic_utils.R")
source("R/load_prs_weights.R")
source("R/annotate_genes_fns.R")
source("R/kegg_fns.R")
psprs_require(c("PGRS_data_path", "pathways_path", "path_KEGG", "path_gencode"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# SETTINGS

# curated genes and pathways spreadsheet ----
###

curated_pathways_name = "OSA_Pathways"
pathways_freeze_date = "20201117" # 20201117 adds the "ADDITIONAL PATHWAYS" for CALM1 based on the reverse lookup of pathways from gene

# choose pathways to annotate ----
###
kegg_mode = Sys.getenv("PSPRS_KEGG_PATHWAYS", "curated")
gene_sets_requested = trimws(strsplit(Sys.getenv("PSPRS_GENE_SETS", "NFkB_TNF"), ",")[[1]])
curated_pathways_name = "OSA_pathways"
annotation_date = Sys.getenv("PSPRS_ANNOTATION_DATE")
if(!nzchar(annotation_date)) stop("Set PSPRS_ANNOTATION_DATE to the date naming the run_annotation.R outputs")
dir.create(pathways_path, recursive = TRUE, showWarnings = FALSE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(kegg_mode == "curated") source(file.path("scripts", "pathways", "load_curated_pathways.R"))
source(file.path("scripts", "pathways", "fetch_kegg_gene_lists.R"))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# set PGRS input/output file paths ----
###

PGRS_name = "PGRS_CAD_Khera"
PGRS_out_name = Sys.getenv("PSPRS_PRS_NAME", "PGRS_CAD_Khera_LDPred_rho0.001")
# The sequentially annotated PGRS file with multiple rows of ENSGs per SNP
PGRS_gene_annot_infile = paste("gene_annot",  annotation_date, "sequential.annotated.SNPs.to.ENSGs", paste(PGRS_out_name, "rds", sep = "."), sep = "_")

# The PGRS with one row per SNP (the scaffold that receives the pathway columns)
PGRS.trained_file = Sys.getenv("PSPRS_PRS_FILE", "CoronaryArteryDisease_PRS_LDpred_rho0.001_v3.txt")
PGRS_format = Sys.getenv("PSPRS_PRS_FORMAT", "khera")
PGRS_colnames_key = switch(PGRS_format,
  khera       = c(variant = "variant", chr = "chr", pos = "position_hg19", effect_weight = "effect_weight",
                  effect_allele = "effect_allele", A1 = "A1", A2 = "A2"),
  pgs_catalog = c(chr = "chr_name", pos = "chr_position", effect_weight = "effect_weight",
                  effect_allele = "effect_allele", A2 = "other_allele"),
  stop("PSPRS_PRS_FORMAT must be khera or pgs_catalog"))
PGRS_MAF = load_split_save_PGRS(
  file_path = file.path(PGRS_data_path, PGRS.trained_file),
  colnames_key = PGRS_colnames_key,
  skip_lines = sum(grepl("^#", readLines(file.path(PGRS_data_path, PGRS.trained_file), n = 100))),
  split_save_TF = FALSE, verbose_TF = FALSE, return_type = "GRanges")$GRanges


# Filter annotation sources:
###
# Note this does not change the order previously set in assign_genes_cascade.R
select.annotation.sources = c(
  # unique(sequential_ann$source)
  GENCODE_exons = "GENCODE_coding_exon", FANTOM5 = "FANTOM5", PCHiC = "PCHiC", GTEx = "GTEx",
  genehancer = "genehancer_thrsh.elite", GENCODE_UTR = "GENCODE_coding_UTR", GENCODE_flanked = "GENCODE_flanked", GENCODE_noncoding = "GENCODE_noncoding"
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source(file.path("scripts", "pathways", "annotate_pathways.R"))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print(length(KEGG_allPaths_ENSG_list))
print(length(KEGG_slctPaths_ENSG_list))
print(length(GSPGRS_allPGRS_ENSG_list))
print(length(GSPGRS_slctPath_ENSG_list))

KEGG_pathwaysENSG_list = c(KEGG_slctPaths_ENSG_list, GSPGRS_slctPath_ENSG_list)

Khera_CAD_PGRS_MAF_pathgenes = annotate_PGRS_pathways(pathways_list = KEGG_pathwaysENSG_list,
                                                      PGRS_multi_annotated = sequential_ann,
                                                      PGRS = PGRS_MAF,
                                                      verbose_TF = FALSE)

Khera_CAD_PGRS_MAF_paths = Khera_CAD_PGRS_MAF_pathgenes$PGRS_paths
KEGGpathENSGs_in_PGRS    = Khera_CAD_PGRS_MAF_pathgenes$KEGGpathENSGs_in_PGRS

# add the "all_SNPs" pathway i.e. all SNPs in PGRS are included. (For comparison purposes).
# note this does not have a corresponding gene list as some SNPs are not assigned to genes.
Khera_CAD_PGRS_MAF_paths$all_SNPs = TRUE

# output files (named by the annotation sources actually present)
source.names = paste(names(select.annotation.sources)[select.annotation.sources %in% unique(sequential_ann$source)], collapse = ".")
PGRS_paths_annot_outfile =   paste( "psPGRS", annotation_date, PGRS_name, genes_name, source.names, "annotated_PS.PGRS.rds", sep = "_")
PGRS_pathway_genes_outfile = paste( "psPGRS", annotation_date, PGRS_name, genes_name, "pathway_genes.rds", sep = "_")

# SAVE pathway-annotated PGRS
Khera_CAD_PGRS_MAF_paths = sortSeqlevels(Khera_CAD_PGRS_MAF_paths)
Khera_CAD_PGRS_MAF_paths = sort(Khera_CAD_PGRS_MAF_paths, ignore.strand = TRUE)

saveRDS(Khera_CAD_PGRS_MAF_paths, file = file.path(PGRS_data_path, PGRS_paths_annot_outfile))
saveRDS(KEGGpathENSGs_in_PGRS, file = file.path(PGRS_data_path, PGRS_pathway_genes_outfile ))

rm(Khera_CAD_PGRS_MAF_paths, Khera_CAD_PGRS_MAF_pathgenes )
