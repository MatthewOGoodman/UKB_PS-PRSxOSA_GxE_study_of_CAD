# Sourced by scripts/annotation/run_annotation.R after annotate_variants_by_source.R.
# GeneHancer elements are read and combined only when PSPRS_GENEHANCER_CSV is set.

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load Annotations Files and Clean Up / Unify Annotation Format  ----

PGRS_suffix = paste(PGRS_out_name, "rds", sep = ".")

###
# GenCode coding
###

PGRS_annot_GenCodeCoding = readRDS(
  file = file.path( PGRS_data_path, paste("gene_annot", annotation_date, "GENCODE_coding", PGRS_suffix, sep = "_"))
)
PGRS_annot_GenCodeCoding$ENSG = PGRS_annot_GenCodeCoding$ENSG_gencode

# "five_prime_UTR", "CDS", "start_codon", "stop_codon", "three_prime_UTR", "UTR", "stop_codon_redefined_as_selenocysteine"
summary(as.factor(PGRS_annot_GenCodeCoding$type)) 

PGRS_annot_GenCode_coding_exon = 
  PGRS_annot_GenCodeCoding[ PGRS_annot_GenCodeCoding$type %in% c("exon", "start_codon", "stop_codon") ]
PGRS_annot_GenCode_coding_exon$source = "GENCODE_coding_exon"

PGRS_annot_GenCode_coding_UTR = 
  PGRS_annot_GenCodeCoding[ PGRS_annot_GenCodeCoding$type  %in% c("five_prime_UTR", "three_prime_UTR", "UTR") ]
PGRS_annot_GenCode_coding_UTR$source = "GENCODE_coding_UTR"

###
# GenCode noncoding
###

PGRS_annot_GenCodeNoncoding = readRDS(
  file = file.path( PGRS_data_path, paste("gene_annot", annotation_date, "GENCODE_noncoding", PGRS_suffix, sep = "_"))
)
PGRS_annot_GenCodeNoncoding$ENSG = PGRS_annot_GenCodeNoncoding$ENSG_gencode
PGRS_annot_GenCodeNoncoding$source = "GENCODE_noncoding"

summary(as.factor(PGRS_annot_GenCodeNoncoding$type)) # type is all "gene"

###
# GenCode upstream downstream
###

PGRS_annot_GenCodeFlank = readRDS(
  file = file.path( PGRS_data_path, paste("gene_annot", annotation_date, "GENCODE_flanked_20kbUp_10kbDn", PGRS_suffix, sep = "_"))
)
PGRS_annot_GenCodeFlank$ENSG = PGRS_annot_GenCodeFlank$ENSG_gencode
PGRS_annot_GenCodeFlank$source = "GENCODE_flanked"

###
# GeneHancer
###
if(nzchar(GeneHancer.file)){

PGRS_annot_GeneHancer = readRDS(
  file = file.path( PGRS_data_path, paste("gene_annot", annotation_date, "GeneHancer", PGRS_suffix, sep = "_")) 
)
PGRS_annot_GeneHancer$ENSG = PGRS_annot_GeneHancer$ENSG_genehancer
PGRS_annot_GeneHancer$source = "genehancer"

### Limit to "more elite" enhancers and enhancer-gene associations
###

# unfortunately downloaded data does not indicate "elite" enhancers and enhancer-gene associations
# these are defined as enhancers and enhancer-gene associations with two ore more data sources supporting them
# Nevertheless it appears possible the genecards GeneHancer v4.4 may actually be limited to double-elite?

GeneHancer_elite.def = "threshold"
if(GeneHancer_elite.def == "threshold"){
  enhancer_score_thresh = 0.3; gene_score_thresh = 4
  PGRS_annot_GeneHancer_elite = PGRS_annot_GeneHancer[PGRS_annot_GeneHancer$enhancer_score > enhancer_score_thresh & PGRS_annot_GeneHancer$gene_score > gene_score_thresh, ]
  PGRS_annot_GeneHancer_elite$source = "genehancer_thrsh.elite"
  PGRS_annot_GeneHancer_non.elite = PGRS_annot_GeneHancer[PGRS_annot_GeneHancer$enhancer_score <= enhancer_score_thresh | PGRS_annot_GeneHancer$gene_score <= gene_score_thresh, ]
  PGRS_annot_GeneHancer_non.elite$source = "genehancer_thrsh.non.elite"
  length(PGRS_annot_GeneHancer_elite); length(PGRS_annot_GeneHancer_non.elite)
}

rm(PGRS_annot_GeneHancer)

# https://academic.oup.com/database/article/doi/10.1093/database/bax028/3737828
# "elite enhancers" generally appear to have enhancer scores of 0.5 and greater (see Fig. S5) or ~0.3 and greater (genecards examples)
# "elite gene-enhancer associations" generally have gene-enhancer association scores of ~3 or more (see genecards examples)
# the minimum eQTL enhancer score is ~5 and it roughly indicates an association p-value.
# Not all gene-enhancer association scores are commensurate  (see Fig. S9), 
# but the high-quality data (not just proximity based) has gene-enhancer scores of more than 5.
} # end GeneHancer


###
# GTEx
###

PGRS_annot_GTEx_art_hrt_wlbld_egenes = readRDS(
  file = file.path( PGRS_data_path, paste("gene_annot", annotation_date, "GTEx_art_hrt_wlbld_eQTLs", PGRS_suffix, sep = "_")) 
  # artery heart and whole blood
)
PGRS_annot_GTEx_art_hrt_wlbld_egenes$source = "GTEx"
PGRS_annot_GTEx_art_hrt_wlbld_egenes$ENSG = PGRS_annot_GTEx_art_hrt_wlbld_egenes$ENSG_GTEx


###
# FANTOM5
###

PGRS_annot_FANTOM5_17bloodcell_enhancer_genes = readRDS(
  file = file.path( PGRS_data_path, paste("gene_annot", annotation_date, "FANTOM5_17bloodcell_enhancergenes", PGRS_suffix, sep = "_"))
)
PGRS_annot_FANTOM5_17bloodcell_enhancer_genes = PGRS_annot_FANTOM5_17bloodcell_enhancer_genes[!is.na(PGRS_annot_FANTOM5_17bloodcell_enhancer_genes$ENSG_fantom)]
PGRS_annot_FANTOM5_17bloodcell_enhancer_genes$source = "FANTOM5"
PGRS_annot_FANTOM5_17bloodcell_enhancer_genes$ENSG = PGRS_annot_FANTOM5_17bloodcell_enhancer_genes$ENSG_fantom


###
# PCHiC
###

PGRS_annot_PCHiC_Javierre_enhancer_genes = readRDS(
  file = file.path( PGRS_data_path, paste("gene_annot", annotation_date, "PCHiC_Javierre_enhancergenes", PGRS_suffix, sep = "_"))
)
PGRS_annot_PCHiC_Javierre_enhancer_genes$source = "PCHiC"
PGRS_annot_PCHiC_Javierre_enhancer_genes$ENSG = PGRS_annot_PCHiC_Javierre_enhancer_genes$ENSG_PCHiC


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Combined Annotations List  ----
combined_annot_list = list()
####

combined_annot_list[["PGRS_annot_GENCODE_coding_exon"]] = 
  PGRS_annot_GenCode_coding_exon[ ,c("wgt", "effect_allele", "A1", "A2", "variant", "source", "ENSG")]

combined_annot_list[["PGRS_annot_GENCODE_coding_UTR"]] = 
  PGRS_annot_GenCode_coding_UTR[ ,c("wgt", "effect_allele", "A1", "A2", "variant", "source", "ENSG")]

combined_annot_list[["PGRS_annot_GENCODE_flanked_genes"]] = 
  PGRS_annot_GenCodeFlank [ PGRS_annot_GenCodeFlank$gene_type %in% c("protein_coding"),
                            c("wgt", "effect_allele", "A1", "A2", "variant", "source", "ENSG")]

combined_annot_list[["PGRS_annot_GENCODE_noncoding"]] = 
  PGRS_annot_GenCodeNoncoding [ , c("wgt", "effect_allele", "A1", "A2", "variant", "source", "ENSG")]

rm(PGRS_annot_GenCode_coding_exon, PGRS_annot_GenCode_coding_UTR, PGRS_annot_GenCodeFlank, PGRS_annot_GenCodeCoding, PGRS_annot_GenCodeNoncoding)

combined_annot_list[["PGRS_annot_GTEx"]] = 
  PGRS_annot_GTEx_art_hrt_wlbld_egenes [ ,c("wgt", "effect_allele", "A1", "A2", "variant", "source", "ENSG")]
rm(PGRS_annot_GTEx_art_hrt_wlbld_egenes)

combined_annot_list[["PGRS_annot_FANTOM5"]] = 
  PGRS_annot_FANTOM5_17bloodcell_enhancer_genes [ ,c("wgt", "effect_allele", "A1", "A2", "variant", "source", "ENSG")]
rm(PGRS_annot_FANTOM5_17bloodcell_enhancer_genes)

combined_annot_list[["PGRS_annot_PCHiC"]] = 
  PGRS_annot_PCHiC_Javierre_enhancer_genes [ ,c("wgt", "effect_allele", "A1", "A2", "variant", "source", "ENSG")]
rm(PGRS_annot_PCHiC_Javierre_enhancer_genes)

if(nzchar(GeneHancer.file)){
combined_annot_list[["PGRS_annot_GeneHancer_elite"]] = 
  PGRS_annot_GeneHancer_elite [ ,c("wgt", "effect_allele", "A1", "A2", "variant", "source", "ENSG")]
rm(PGRS_annot_GeneHancer_elite)

combined_annot_list[["PGRS_annot_GeneHancer_non.elite"]] = 
  PGRS_annot_GeneHancer_non.elite [ ,c("wgt", "effect_allele", "A1", "A2", "variant", "source", "ENSG")]
rm(PGRS_annot_GeneHancer_non.elite)
}

# GRangesList(combined_annot_list)

saveRDS(
  combined_annot_list,
  file = file.path( PGRS_data_path, paste("gene_annot", annotation_date, "combined_annotation_sources", PGRS_suffix, sep = "_")) 
)
