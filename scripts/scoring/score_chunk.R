# scripts/scoring/score_chunk.R — score one chunk of the pathway-annotated PRS against .bgen files.
#
# One run of this script scores one chunk (a block of PRS variants within one chromosome). Run it
# once per chunk, e.g. as a scheduler array job with one task per chunk; the chunk number comes
# from PSPRS_TASK_ID (or SLURM_ARRAY_TASK_ID, or SGE_TASK_ID). Then run combine_chunks.R.
#
# Run from the repository root:  Rscript scripts/scoring/score_chunk.R
#
# Required environment: PSPRS_DATA_DIR, PSPRS_GENOTYPE_DIR, PSPRS_OUT_DIR, PSPRS_RUN_ID, and a task id.
# Optional environment (defaults are the settings used for the paper):
#   PSPRS_MEMBERSHIP_FILE         pathway-annotated PRS (.rds GRanges) in PSPRS_DATA_DIR
#   PSPRS_BLOCKSIZE               variants per genotype read                 (default 200)
#   PSPRS_CHUNK_SIZE              maximum variants per chunk                 (default 100000)
#   PSPRS_EFFECT_ALLELE_COL       column holding the weighted allele         (default effect_allele)
#   PSPRS_GENOTYPE_FILE_PATTERN   bgen file name; "{chr}" is replaced by the chromosome as
#                                 "chr1".."chr22"                            (default ukb_imp_{chr}_v3.bgen)
#   PSPRS_CHR_LABEL               chromosome names inside the bgen file:
#                                 zero_pad ("01".."22", UK Biobank; default), number ("1".."22"), chr ("chr1")

library(GenomicRanges)
library(GenomeInfoDb)
library(rbgen)
library(dplyr)
library(tibble)
library(stringr)
library(purrr)

source("config/paths.R")
source("R/genomic_utils.R")
source("R/genotype_dosage_fns.R")
source("R/psprs_score_bgen.R")
psprs_require(c("Khera_path", "genotype_dir", "psprs_out_dir", "psprs_run_id", "psprs_task_id"))

# set parameters
####
# set blocksize (n SNPs per iteration)
blocksize = as.integer(Sys.getenv("PSPRS_BLOCKSIZE", "200"))

# set chunksize
PGRS_chunk_col = "chunk"
chunk_size = as.integer(Sys.getenv("PSPRS_CHUNK_SIZE", "100000"))

effect_allele_col = Sys.getenv("PSPRS_EFFECT_ALLELE_COL", "effect_allele")

genotype_file_pattern = Sys.getenv("PSPRS_GENOTYPE_FILE_PATTERN", "ukb_imp_{chr}_v3.bgen")
genotype_filename_fn = function(chr_name) gsub("{chr}", chr_name, genotype_file_pattern, fixed = TRUE)

chr_label = Sys.getenv("PSPRS_CHR_LABEL", "zero_pad")
chr_label_fn = switch(chr_label,
  zero_pad = function(chr) { n = sub("^chr", "", chr); ifelse(grepl("^[0-9]$", n), paste0("0", n), n) },
  number   = function(chr) sub("^chr", "", chr),
  chr      = function(chr) chr,
  stop("PSPRS_CHR_LABEL must be zero_pad, number or chr"))

# grep for columns in pathway-annotated PGRS with logicals for SNP inclusion
pathway.column.grep = "hsa\\d{5}|gene|pathway|all_SNPs" # note this won't work for gene names unless you paste "gene" in front of each one

# set PGRS_directory
PGRS_directory = Khera_path
PGRS_file = Sys.getenv("PSPRS_MEMBERSHIP_FILE",
  "psPGRS_20201209_PGRS_CAD_Khera_NFkB_TNF_GENCODE_exons.FANTOM5.PCHiC.GTEx.genehancer.GENCODE_UTR.GENCODE_flanked.GENCODE_noncoding_annotated_PS.PGRS.rds")

# load pathway-annotated PGRS
####
PGRS = readRDS(file = file.path(PGRS_directory, PGRS_file))

# locate columns devoted to pathways:
# note: PGRS must already contain columns for pathway annotation via T/F inclusion/exclusion for each SNP
pathcols = colnames(mcols(PGRS))[ grepl( pathway.column.grep, colnames(mcols(PGRS)) ) ]
pathcols = pathcols[ ! pathcols %in% c("MAF_wgtd_effect", "chromosome")]


# set PGRS chunks
####

# using this function partitions PGRS by chunks that respect chromosomes (chunks do not overlap multiple chromosomes)
mcols(PGRS)[[PGRS_chunk_col]] = chunk_PGRS_within_chr( PGRS, chunk_size = chunk_size)$chunks
  print( "size of chunks for each task in array")
  print( sort(summary(as.factor(mcols(PGRS)[[PGRS_chunk_col]]))) )
total_chunks = length(unique( mcols(PGRS)[[PGRS_chunk_col]] ))

# run manifest: every task of one run must use the same input and settings
####
run_dir = file.path(psprs_out_dir, psprs_run_id)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
run_manifest = list(membership_file = PGRS_file,
                    membership_md5 = unname(tools::md5sum(file.path(PGRS_directory, PGRS_file))),
                    total_chunks = total_chunks, chunk_size = chunk_size, blocksize = blocksize,
                    effect_allele_col = effect_allele_col, pathcols = pathcols)
manifest_path = file.path(run_dir, "run_manifest.rds")
if(file.exists(manifest_path)){
  if(!identical(readRDS(manifest_path), run_manifest))
    stop("run ", psprs_run_id, " already exists with a different input or settings; use a new PSPRS_RUN_ID")
} else {
  saveRDS(run_manifest, paste0(manifest_path, ".", Sys.getpid(), ".partial"))
  file.rename(paste0(manifest_path, ".", Sys.getpid(), ".partial"), manifest_path)
}

# chunk for this task
####
PGRS_chunk_name = psprs_task_id
if(!(PGRS_chunk_name %in% as.character(seq_len(total_chunks)))){
  stop("task id ", PGRS_chunk_name, " is not a chunk number; this input has chunks 1..", total_chunks)
}

# subset to PGRS chunk from PGRS: chunk is amount of PGRS handled by one node, via iterations over multiple blocks
PGRS_chunk = PGRS[ grepl( paste0("^",PGRS_chunk_name,"$"), mcols(PGRS)[,PGRS_chunk_col]), ]

# remove full PGRS to save memory, leaving only PGRS_chunk
rm(PGRS)
gc(reset = TRUE)

##############
# Compute PGRS
####

compute_PS_PGRS(

  PGRS_chunk = PGRS_chunk,
  PGRS_chunk_name = PGRS_chunk_name,
  pathway_columns = pathcols,
  # note: PGRS_chunk must have columns with pathway annotation columns with T/F inclusion exclusion for each SNP

  genotype_dir = genotype_dir,
  genotype_filename_fn = genotype_filename_fn,
  # note: this function assumes that chr_name derived from the PGRS seqnames column is the basis for the filenames in the genotype_dir
  out_dir = psprs_out_dir,
  run_id = psprs_run_id,

  blocksize = blocksize,
  subset_samples = NULL,
  effect_allele_col = effect_allele_col,
  chr_label_fn = chr_label_fn,
  verbose_TF = TRUE

  )
