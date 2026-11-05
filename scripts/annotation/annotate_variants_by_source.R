# annotatePGRS_SNPs.to.Genes.R

# all genomic locations referenced to GRCh37 aka hg19

# Useful Bioconductor Resources:
# https://www.bioconductor.org/packages/devel/workflows/vignettes/annotation/inst/doc/Annotation_Resources.html
# Useful Paper for annotating genes:
# https://dx.doi.org/10.12688%2Ff1000research.13577.2

# Possible to add additional enhancer approaches:
# http://acgt.cs.tau.ac.il/focs/
# http://www.enhanceratlas.org/

# Sourced by scripts/annotation/run_annotation.R, which sets the paths (config/paths.R), the PRS
# (PGRS_variants_GR, PGRS_suffix) and annotation_date.

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# LOAD TRAINED PGRS ----
print( "loading PGRS with GWAS variants" )

if(!exists("PGRS_variants_GR")){
PGRS.trained_path = file.path(PGRS_data_path, PGRS.trained_file)

PGRS_variants_GR = load_split_save_PGRS(
  file_path = PGRS.trained_path, 
  out_name = PGRS_out_name, 
  out_folder = PGRS_data_path,
  split_save_TF = FALSE, # split and save by chr
  verbose_TF = FALSE,
  return_type = "GRanges")$GRanges
}

PGRS_suffix = paste(PGRS_out_name, "rds", sep = ".")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# LOAD LIBRARIES ----

library(GenomicFeatures)
# vignette("GenomicFeatures")
# ?txdbmaker::makeTxDbFromGFF


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# GET DATABASE OBJECTS ----

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GenCode consortium (Wellcome Trust Sanger Inst) ----

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GenCode TxDb database ----
####
print( "loading Gencode" )

# path_gencode and GenCode.db come from config/paths.R (PSPRS_GENCODE_GFF); the paper used
# GENCODE gencode.v34lift37 (GRCh37); the download is described in RUNBOOK.md

if(!exists("GenCode.txdb.hg19")){
  GenCode.txdb.hg19 <- txdbmaker::makeTxDbFromGFF(file.path(path_gencode, GenCode.db))
  GenCode.txdb.hg19 <- keepStandardChromosomes(GenCode.txdb.hg19)
}

# GenCode GenomicRanges database ----
####

# Use rtracklayer package to directly load GenCode gff3 file as GRanges
###
GenCode.txdb.hg19_GR <- rtracklayer::import(file.path(path_gencode, GenCode.db))
# discard superfluous columns
mcols(GenCode.txdb.hg19_GR) = 
  mcols(GenCode.txdb.hg19_GR)[, 
                              c( "source", "type", "score", "phase", 
                                 "ID", "Parent", # Note: ID can be ENSG or ENST, where exons are ENST ID's with suffix :1, :2 :3 etc
                                 "gene_id", "gene_type", "gene_name", "hgnc_id",                        
                                 "transcript_id",  "transcript_type", "transcript_name",
                                 "exon_id", "protein_id", "ccdsid" # exon_id is ENSE ID's
                              ) ]

# separate the strings in column gene_id
mcols(GenCode.txdb.hg19_GR) =
  separate( as_tibble(mcols(GenCode.txdb.hg19_GR)),
            col = "gene_id",
            sep = "\\.",
            into = c("ENSG_gencode", "ENSG_gencode_vrs")
            # GENCODE gene and transcript ids carry a version suffix after the dot.
            # See https://www.gencodegenes.org/pages/data_format.html
  )
mcols(GenCode.txdb.hg19_GR)$ENSG_gencode_vrs = NULL

# some summary stats, but e.g. the same exons appear in multiple transcript rows
# t(t(summary (mcols(GenCode.txdb.hg19_GR)$type)))
# t(t(summary ( as.factor( mcols(GenCode.txdb.hg19_GR)$gene_type))))
# table( mcols(GenCode.txdb.hg19_GR)$type, as.factor( mcols(GenCode.txdb.hg19_GR)$gene_type) )



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# ANNOTATE GENES ----


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GenCode genomic locations  ----
####

# Annotate: GenCode coding genes ----
###
print( "Annotate GenCode coding" )

GenCode.txdb.hg19_GR.coding = GenCode.txdb.hg19_GR[
  GenCode.txdb.hg19_GR$gene_type == "protein_coding" & 
    GenCode.txdb.hg19_GR$type %in% 
    c("exon", # exclude transcript and gene, to avoid involving intronic regions. Introns will be covered by the "flanked" genes below.
      "five_prime_UTR", "CDS", "start_codon", "stop_codon", "three_prime_UTR", "UTR", "stop_codon_redefined_as_selenocysteine")
  ]
summary(GenCode.txdb.hg19_GR.coding$type)
PGRS_anno_gencode_coding = get_PGRS_annotations_from_GRanges(PGRS_variants_GR = PGRS_variants_GR, annotations_GR = GenCode.txdb.hg19_GR.coding)
summary(PGRS_anno_gencode_coding$type)

# Hold out UTR region separately, because UTR's can regulate non-adjacent genes, even though UTR's are within exons
PGRS_anno_gencode_coding.df = data.frame(PGRS_anno_gencode_coding)
# reorder so first instance is UTR
PGRS_anno_gencode_coding.df = rbind(
  PGRS_anno_gencode_coding.df[PGRS_anno_gencode_coding.df$type %in% c( "five_prime_UTR", "three_prime_UTR", "UTR" ),],
  PGRS_anno_gencode_coding.df[PGRS_anno_gencode_coding.df$type %in% c( "exon", "CDS", "start_codon", "stop_codon", "stop_codon_redefined_as_selenocysteine" ),]
)
PGRS_anno_gencode_coding = makeGRangesFromDataFrame(PGRS_anno_gencode_coding.df, keep.extra.columns = TRUE)
# now when duplicates are removed UTR SNPs are assigned to UTR
dup_rows = duplicated(as.data.table(mcols(PGRS_anno_gencode_coding)[c("variant", "ENSG_gencode")]), sorted = FALSE)
PGRS_anno_gencode_coding = PGRS_anno_gencode_coding[!dup_rows] 
summary(PGRS_anno_gencode_coding$type)



length(PGRS_anno_gencode_coding); length(unique(PGRS_anno_gencode_coding$variant)); length(unique(PGRS_anno_gencode_coding$ENSG_gencode)); 

# save GenCode coding variants:
saveRDS(
  PGRS_anno_gencode_coding[ ,c("variant", "wgt", "effect_allele", "A1", "A2", 
                                  "ENSG_gencode", "gene_name", "hgnc_id",
                                  "type", "ID", "Parent", 
                                  "transcript_id", "transcript_type", "exon_id" )],
  file = file.path( PGRS_data_path, paste("gene_annot", annotation_date, "GENCODE_coding", PGRS_suffix, sep = "_") )
)


# Annotate: GenCode noncoding genes ----
###
print( "Annotate GenCode noncoding" )

GenCode.txdb.hg19.noncoding = GenCode.txdb.hg19_GR[GenCode.txdb.hg19_GR$gene_type != "protein_coding"]
# For  non-protein coding genes, not sure that you can safely exclude transcript and gene, and rely on exons alone.
# This will mean that all non-coding gene introns will be included, but this seems fine, 
# given that we are putting these non-coding genes last in the list of annoations,
# and KEGG  does not have very many of these noncoding genes/transcripts.
# NOTE: there is an argument for augmenting non-coding genes via flanking.

PGRS_anno_gencode_noncoding = get_PGRS_annotations_from_GRanges(PGRS_variants_GR = PGRS_variants_GR, annotations_GR = GenCode.txdb.hg19.noncoding)
# summary(as.factor(PGRS_anno_gencode_noncoding$gene_type))
# gene_types = unique(PGRS_anno_gencode_noncoding$gene_type)
# lapply(gene_types, function(x) head(PGRS_anno_gencode_noncoding[PGRS_anno_gencode_noncoding$gene_type %in% x]) )
# colnames(mcols(PGRS_anno_gencode_noncoding))

dup_rows = duplicated(as.data.table(mcols(PGRS_anno_gencode_noncoding)[c("variant", "ENSG_gencode")]), sorted = FALSE)
PGRS_anno_gencode_noncoding = PGRS_anno_gencode_noncoding[!dup_rows] 
length(PGRS_anno_gencode_noncoding); length(unique(PGRS_anno_gencode_noncoding$variant)); length(unique(PGRS_anno_gencode_noncoding$ENSG_gencode)); 

# save GenCode noncoding variants:
saveRDS(
  PGRS_anno_gencode_noncoding[ ,c("variant", "wgt", "effect_allele", "A1", "A2", 
                               "ENSG_gencode", "gene_name", "hgnc_id",
                               "type", "ID", "Parent", 
                               "transcript_id", "transcript_type", "exon_id" )],
  
  file = file.path( PGRS_data_path, paste("gene_annot", annotation_date, "GENCODE_noncoding", PGRS_suffix, sep = "_") )
)

rm(GenCode.txdb.hg19_GR.coding)
rm(GenCode.txdb.hg19.noncoding)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Flank GenCode gene coordinates ----
####
print("Annotate flanked GENCODE protein coding genes")


GenCode.txdb.hg19_GR.coding = GenCode.txdb.hg19_GR[
# subset to genes with standard UCSC chromosomal locations and protein_coding gene biotype
  seqnames(GenCode.txdb.hg19_GR) %in% c(paste0("chr", 1:22),"chrX", "chrY") &
    GenCode.txdb.hg19_GR$gene_type %in% c("protein_coding") &
  # subset to protein coding _genes_ only, excluding transcripts and exons, which are redundant for flanking:
    GenCode.txdb.hg19_GR$type %in% c("gene") 
  ]

# flank genes:
genesGenCode_5prime = trim(flank(GenCode.txdb.hg19_GR.coding, width = 20000, start = TRUE)) 
genesGenCode_3prime = trim(flank(GenCode.txdb.hg19_GR.coding, width = 10000, start = FALSE))
# recall: 5prime is start (left side on forward strand, right side on reverse), 3prime is stop.
genesGenCode_flanked = append(GenCode.txdb.hg19_GR.coding, genesGenCode_5prime)
genesGenCode_flanked = append(genesGenCode_flanked, genesGenCode_3prime)

rm(GenCode.txdb.hg19_GR.coding, genesGenCode_5prime, genesGenCode_3prime)


# Annotate: GenCode flanked genes ----
####

print( "Annotate GenCode flanked" )

PGRS_GenCodeFlankAnno = 
  get_PGRS_annotations_from_GRanges(PGRS_variants_GR = PGRS_variants_GR, annotations_GR = genesGenCode_flanked)


# only retain important columns
mcols(PGRS_GenCodeFlankAnno) = 
   mcols(PGRS_GenCodeFlankAnno)[
  c( "wgt", "effect_allele", "A1", "A2", "variant",
     "source", "type", "ID", "Parent",
     "ENSG_gencode", "gene_type", "gene_name", "hgnc_id",
     "transcript_id", "transcript_type", "transcript_name", 
     "exon_id", "protein_id", "ccdsid"  )
  ]

# save PGRS_GenCodeFlankAnno:
####

saveRDS(
  PGRS_GenCodeFlankAnno,
  file = file.path( PGRS_data_path, paste("gene_annot", annotation_date, "GENCODE_flanked_20kbUp_10kbDn", PGRS_suffix, sep = "_") )
)
rm(PGRS_GenCodeFlankAnno)
rm(genesGenCode_flanked)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Scripts below adapted from Ferrero et al, F1000 Research article ----
# Enrico Ferrero: Using regulatory genomics data to interpret the function of 
# disease variants and prioritise genes from expression studies
# https://dx.doi.org/10.12688%2Ff1000research.13577.2
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

####
# Annotate: GeneHancer genes ----
####


print( "Annotate GeneHancer genes" )

# GeneHancer is licensed separately and is not distributed with this repository. The section
# runs only when PSPRS_GENEHANCER_CSV names a GeneHancer csv; otherwise the source is skipped.
if(nzchar(GeneHancer.file)){

### Prepare GeneHancer database 
###

# GeneHancer: https://www.genecards.org/GeneHancer_version_4-4
# Note this GeneHancer version is based on data "converted to hg38" according to https://academic.oup.com/database/article/doi/10.1093/database/bax028/3737828
GeneHancer_data_wide = read_csv( file.path(path_genehancer, GeneHancer.file))
GeneHancer_data_wide = rename(GeneHancer_data_wide, enhancer_score = score)
hist(GeneHancer_data_wide$enhancer_score, xlim = c(0,2), breaks = seq(0,5,0.1))

GeneHancer_data_sep = separate(GeneHancer_data_wide, 
                               col = attributes,
                               into = c("genehancer_id", paste( c("connected_gene", "score"), rep(1:100, each = 2), sep = "." ) ),
                               sep = "\\;")

slct_cols = colnames(GeneHancer_data_sep)[! colnames(GeneHancer_data_sep) %in%
                                            c("chrom", "source", "feature name", "start", "end", "enhancer_score", "strand", "frame") ]

GeneHancer_data_enhancer.genes_list = apply(GeneHancer_data_sep[, slct_cols], 1, function(x) x[!is.na(x)])
Num_genes_per_enhancer = (sapply(GeneHancer_data_enhancer.genes_list, length)-1)/2
max( Num_genes_per_enhancer)
sum(Num_genes_per_enhancer == 0)
hist( Num_genes_per_enhancer, xlim = c(1,20), breaks = 0:100)


GeneHancer_data_enhancer.genes_list_dfs = lapply(
  GeneHancer_data_enhancer.genes_list,
  function(x) suppressWarnings(
    data.frame(
    genehancer_id =  x["genehancer_id"], 
    matrix(x[-1], byrow = TRUE, ncol = 2),
    stringsAsFactors = FALSE)) )

GeneHancer_data_enhancer.genes_df = bind_rows(GeneHancer_data_enhancer.genes_list_dfs)

nrow_GeneHancer_data = nrow(GeneHancer_data_wide)
GeneHancer_data_long = data.frame(
  GeneHancer_data_sep[ rep(1:nrow_GeneHancer_data, times = Num_genes_per_enhancer), 1:9],
  GeneHancer_data_enhancer.genes_df)


GeneHancer_data_long = 
  GeneHancer_data_long %>%
  separate(
    col = genehancer_id,
    into = c("genehancer_id_nm", "genehancer_id"),
    sep = "\\=") %>%
  separate(
    col = X1,
    into = c("connected_gene_nm", "connected_gene"),
    sep = "\\=") %>%
  separate(
    col = X2,
    into = c("score_nm", "gene_score"),
    sep = "\\=")

GeneHancer_data_long[,c("genehancer_id.1", "genehancer_id_nm", "connected_gene_nm", "score_nm")] = NULL

GeneHancer_data_GR = 
  makeGRangesFromDataFrame(df = GeneHancer_data_long,
                           keep.extra.columns=TRUE,
                           ignore.strand=FALSE,
                           seqinfo=NULL,
                           seqnames.field=c("chrom"),
                           start.field="start",
                           end.field=c("end"),
                           strand.field="strand",
                           starts.in.df.are.0based=FALSE)

chainfile.hg38.to.hg.19 = get_chain(hg_in = "hg38", hg_out = "hg19")
# str(chainfile.hg38.to.hg.19)
GeneHancer_data_hg.19 = unlist(liftOver(GeneHancer_data_GR, chainfile.hg38.to.hg.19))

# Obtain ENSG ID's for GeneHancer HGNC genes using GenCode database
###
# Direct lookup method
GeneHancer_data_hg.19$Gene = GeneHancer_data_hg.19$connected_gene
ENSG_grepl = (grepl( "^ENSG(\\d{11})$", GeneHancer_data_hg.19$Gene))
GeneHancer_data_hg.19$ENSG_genehancer = GeneHancer_data_hg.19$Gene
GenCode.txdb.hg19_GR <- rtracklayer::import(file.path(path_gencode, GenCode.db))
HGNC_vec = GeneHancer_data_hg.19$Gene [ !ENSG_grepl ]
ENSG_vec = HGNC_to_ENSG(HGNC_vec, GenCode.txdb.hg19_GR)
GeneHancer_data_hg.19$ENSG_genehancer [ !ENSG_grepl ] = ENSG_vec
head(GeneHancer_data_hg.19)

# Note: not all genes are mapped, unmapped genes appear to be mostly pseudogenes:
# GeneHancer_data_hg.19$Gene[ is.na(GeneHancer_data_hg.19$ENSG_genehancer)]


### Annotate PGRS with GeneHancer Data:
###

PGRS_anno_GeneHancer = 
  get_PGRS_annotations_from_GRanges(PGRS_variants_GR = PGRS_variants_GR, annotations_GR = GeneHancer_data_hg.19)

rm(GeneHancer_data_wide, 
   GeneHancer_data_sep,
   GeneHancer_data_enhancer.genes_list,
   GeneHancer_data_enhancer.genes_list_dfs,
   GeneHancer_data_enhancer.genes_df,
   GeneHancer_data_long,
   GeneHancer_data_GR
   )
rm(GeneHancer_data_hg.19)

### Get rid of NA's from missing genehancer ENSG's
PGRS_anno_GeneHancer = PGRS_anno_GeneHancer[!is.na(PGRS_anno_GeneHancer$ENSG_genehancer)]
length(unique(PGRS_anno_GeneHancer$variant))


# save data GeneHancer:
####
saveRDS(
  PGRS_anno_GeneHancer,
  file = file.path( PGRS_data_path, paste("gene_annot", annotation_date, "GeneHancer", PGRS_suffix, sep = "_") )
)

# sanity checks GeneHancer:
####
# length(PGRS_anno_GeneHancer)
# length(unique(PGRS_anno_GeneHancer$variant))
# length(PGRS_anno_GeneHancer) / length(unique(PGRS_anno_GeneHancer$variant)) # genes per variant
# length(unique(PGRS_anno_GeneHancer$variant)) / length(unique(PGRS_variants_GR$variant)) # percent of SNPs annotated

} # end GeneHancer (skipped when PSPRS_GENEHANCER_CSV is unset)




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

####
# Annotate: GTEx egenes ----
####

print( "Annotate GTEx egenes" )

# note that GTEx v7 is GRCh37.
# this may change with GTEx v8

all_gtex_tissues_raw = tibble(filename = dir(path_GTEx))  
all_gtex_tissues =  dplyr::filter( all_gtex_tissues_raw , grepl("v7.signif_variant_gene_pairs.txt.gz", filename))

all_gtex_tissues_tbl = 
  all_gtex_tissues %>% 
  tidyr::separate(filename, 
                  into = c("tissue", "file_type"),
                  sep = "\\.v7",
                  remove = FALSE #default:T remove original column
  )

# Select Tissue types: (Artery, Heart, Whole_Blood)
tissue_types = c(
  "Artery",
  "Heart",
  "Whole_Blood"
)

gtex_files_slct = dplyr::filter(all_gtex_tissues_tbl, grepl(paste(tissue_types, collapse = "|"), tissue))
# gtex_file_suffx = ".v7.signif_variant_gene_pairs.txt.gz"

# cycle through tissues and map eQTLs
####

gtex_tissues_list = list()
for(i in 1:nrow(gtex_files_slct) ){
  tissue_name = as.character(gtex_files_slct[i,"tissue"])
  file_name = as.character(gtex_files_slct[i,"filename"])
  
  PGRSvars_found_in_tissue_eQTLs_GR = tissue_eQTLs_GR_in_PGRS_variants = tissue_eQTLs_GR = NULL
  
  # load data 
  ####
  tissue_eQTLs = read_tsv( gzfile( file.path(path_GTEx, file_name) ))
  tissue_eQTLs$tissue = tissue_name
  
  tissue_eQTLs = 
    tissue_eQTLs %>% 
    tidyr::separate(
      variant_id, 
      into = c("chr", "pos", "A1_GTEx", "A2_GTEx", "hg_build"),
      sep = "\\_",
      remove = FALSE #default:T remove original column
    )
  tissue_eQTLs$start = tissue_eQTLs$end = tissue_eQTLs$pos 
  tissue_eQTLs$variant_GTEx = with(tissue_eQTLs, paste( paste0("chr", chr), pos, A1_GTEx, A2_GTEx, sep = ":"))
  
  # convert to GRanges 
  ####
  tissue_eQTLs_GR = makeGRangesFromDataFrame(tissue_eQTLs, keep.extra.columns = TRUE)
  seqlevelsStyle(tissue_eQTLs_GR) <- "UCSC"
  rm(tissue_eQTLs)
  
  # find overlaps with PGRS 
  ####
  
  PGRSvars_found_in_tissue_eQTLs_GR = 
    get_PGRS_annotations_from_GRanges(PGRS_variants_GR = PGRS_variants_GR, 
                                      annotations_GR = tissue_eQTLs_GR)
  
  PGRSvars_found_in_tissue_eQTLs_GR = as_tibble(PGRSvars_found_in_tissue_eQTLs_GR)
  
  PGRSvars_found_in_tissue_eQTLs_GR = 
    PGRSvars_found_in_tissue_eQTLs_GR %>% 
    tidyr::separate(
      gene_id, 
      into = c("ENSG_GTEx", "ENSG_GTEx_vrs"),
      sep = "\\."
    )

  gtex_tissues_list[[tissue_name]] = 
    PGRSvars_found_in_tissue_eQTLs_GR[,c("seqnames", "start", "end",
                                         "variant", "wgt", "effect_allele", "A1", "A2", 
                                         "variant_GTEx",  "A1_GTEx", "A2_GTEx",
                                         "tissue", "ENSG_GTEx", "ENSG_GTEx_vrs",
                                         # the following would include the GTEx tissue in the column name
                                         # paste("ENSG_GTEx", tissue_name, sep = "_"), paste("ENSG_GTEx", tissue_name, "vrs", sep = "_"),
                                         "tss_distance", "pval_beta")]
  
  rm(PGRSvars_found_in_tissue_eQTLs_GR, tissue_eQTLs_GR, tissue_eQTLs_GR_in_PGRS_variants)
  
}

# bind list to get one data frame for select tissues 
####
gtex_tissues_tbl = bind_rows(gtex_tissues_list)
rm(gtex_tissues_list)


# summaries on GTEx eQTL genes
print("dimension of gtex select tissues variant annotations for intron and intergenic")
print(dim(gtex_tissues_tbl))

# remove duplicates due to the same variants being intronic in multiple transcripts of the same gene 
####
GTEx_slct_tissues_eQTL_genes = distinct(gtex_tissues_tbl, 
                                        tissue, variant_GTEx, ENSG_GTEx,
                                        .keep_all = TRUE)
rm(gtex_tissues_tbl)

print("dimension of gtex select tissues variant annotations for intron and intergenic")
print(dim(GTEx_slct_tissues_eQTL_genes))
length(unique(GTEx_slct_tissues_eQTL_genes$variant_GTEx))

GTEx_slct_tissues_eQTL_genes = makeGRangesFromDataFrame(GTEx_slct_tissues_eQTL_genes, keep.extra.columns = TRUE)

saveRDS(
  GTEx_slct_tissues_eQTL_genes,
  file = file.path( PGRS_data_path, paste("gene_annot", annotation_date, "GTEx_art_hrt_wlbld_eQTLs", PGRS_suffix, sep = "_") )
  
)
# rm(GTEx_slct_tissues_eQTL_genes)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

####
# Annotate: FANTOM5 data ----
####

print( "Annotate FANTOM5 promoter/enhancer regions to genes" )

#   http://enhancer.binf.ku.dk/presets/
# Pre-defined enhancer sets and TSS-enhancer association tracks: Using statistical methods, we have defined several tracks used in Andersson et al.
# These include:
#   tracks where one facet is over-represented (similar, but not identical to so-called tissue-specific transcription)
#   tracks of enhancers used in a certain facet (but not necessarily over-represented)
#   tracks for ubiquitously used enhancers
# Most track sets are defined in the BED format unless specified otherwise and can be downloaded. 
# Each predicted enhancer is described in BED12 format with two blocks denoting the merged regions of transcription initiation on the minus and plus strands. 
# The thickStart and thickEnd columns denote the inferred mid position. 
# The score column gives the maximum pooled expression of TCs used to construct each bidirectional loci.

# Note: may also want to look at cell-type specific or miRNA data, if miRNA can be mapped to promoters and thus to genes.


# uncomment the following line to download the file
# download.file("http://enhancer.binf.ku.dk/presets/enhancer_tss_associations.bed", destfile = "enhancer_tss_associations.bed")
fantom <- read.delim(file.path(path_FANTOM5, FANTOM5.file), skip = 1, stringsAsFactors = FALSE)
head(fantom, 3)

# According to F1000 research: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5850119/
# Everything we need is in the fourth column, name: genomic location of the enhancer, gene identifiers, Pearson correlation coefficient and significance. 

# parse the "name" column: 
fantom = 
  separate( as_tibble(fantom),
            col = "name",
            into = c("chrom", "name", "symbol", "corr", "FDR"),
            sep = ";"
  )

# Now we can extract the genomic locations of the enhancers and the correlation values:
fantom = 
  fantom %>% 
  tidyr::separate(
    chrom, 
    into = c("chr","start", "end"),
    sep = "[:-]",
    remove = FALSE #default:T remove original column
  )

fantom$HGNC_fantom = fantom$symbol
fantom$corr <- sub("R:", "", fantom$corr)
fantom$fdr <- sub("FDR:", "", fantom$FDR)


# now select only enhancer-promoter pairs with corr >= 0.25 and fdr < 1e-5
fantom <- unique(subset(fantom, corr >= 0.25 & fdr < 1e-5, 
                        select = c("chr","start", "end", "HGNC_fantom", "corr", "fdr")))

fantom <- makeGRangesFromDataFrame(fantom, keep.extra.columns = TRUE)

###
# lookup table of GenCode HGNC to ENSG's
fantom$ENSG_fantom = HGNC_to_ENSG(fantom$HGNC_fantom, GenCode.txdb.hg19_GR)
# the MapIds function seems to get some wrong e.g. returning ENSG00000130035 = GALNT8 for KCNA6
# fantom$ENSG_fantom2 = HGNC_to_ENSG_AnnotationDb.MapIds(fantom$HGNC_fantom, GenCode.txdb.hg19, mVals = "first")
# sum(fantom$ENSG_fantom != fantom$ENSG_fantom2, na.rm = TRUE)
# bad = fantom$ENSG_fantom != fantom$ENSG_fantom2; bad[is.na(bad)] = FALSE 
# fantom[ bad,]



# We can now compute the overlap between SNPs and enhancers:
hits_f5 <- findOverlaps(query = PGRS_variants_GR, subject = fantom)
FANTOM5_17bloodcell_enhancer_genes = PGRS_variants_GR[queryHits(hits_f5)] 
fantom_annotation_hits = fantom[subjectHits(hits_f5)] 
mcols(FANTOM5_17bloodcell_enhancer_genes) <- cbind(mcols(FANTOM5_17bloodcell_enhancer_genes), mcols(fantom_annotation_hits))

# FANTOM5_17bloodcell_enhancer_genes <- as_tibble(FANTOM5_17bloodcell_enhancer_genes)
saveRDS(
  FANTOM5_17bloodcell_enhancer_genes,
  file = file.path( PGRS_data_path, paste("gene_annot", annotation_date, "FANTOM5_17bloodcell_enhancergenes", PGRS_suffix, sep = "_") )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

####
# Annotate: Promoter Capture Hi-C data. ----
####

print( "Annotate PC Hi-C genes" )

# https://www.cell.com/cell/fulltext/S0092-8674(16)31322-8
# https://osf.io/u8tzp/ (PCHiC data, CC BY 4.0)
# From DataS1-README.pdf: "In all cases, the data are mapped to the GRCh37 genomic assembly, with further details on each of these files given below.

# the thresholded promoter/open chromatin link file (from the Data S1 archive of Javierre et al. 2016):
PCHiC <- read.delim(file.path(path_PCHiC, "ActivePromoterEnhancerLinks.tsv"), stringsAsFactors = FALSE)
head(PCHiC, 3)

# We will use the InteractionSet package, 
# which is specifically designed for the representation of chromatin interaction data. 
# We start by creating a GInteractions object:

# BiocManager::install(pkgs = "InteractionSet")
# source("http://bioconductor.org/biocLite.R")
# biocLite("InteractionSet")

library(InteractionSet)
promoters <- GRanges(seqnames = PCHiC$baitChr, 
                     ranges = IRanges(start = PCHiC$baitSt, end = PCHiC$baitEnd))
enhancers <- GRanges(seqnames = PCHiC$oeChr, 
                     ranges = IRanges(start = PCHiC$oeSt, end = PCHiC$oeEnd))
PCHiC <- GInteractions(promoters, enhancers)

# we also have to map promoters to their nearest gene transcripts 
# so that we know which genes are regulated by which enhancers

# as above load txdb:
if(!exists("GenCode.txdb.hg19")){
  GenCode.txdb.hg19 <- txdbmaker::makeTxDbFromGFF(file.path(path_gencode, GenCode.db))
  GenCode.txdb.hg19 <- keepStandardChromosomes(GenCode.txdb.hg19)
}
# Add ENSG's to PCHiC object
# first get tsss position for associated gene:
tsss <- promoters(GenCode.txdb.hg19, upstream = 0, downstream = 1, columns = "gene_id")
# then get nearest match between PCHiC promoters and GenCode.txdb.hg19 tsss genes
hits_PCHiC <- nearest(promoters, tsss)
PCHiC$gene_id <- unlist(tsss[hits_PCHiC]$gene_id)

# Overlaps between PGRS and enhancers (the second region of the GInteractions object) :
hits <- findOverlaps(PGRS_variants_GR, PCHiC, use.region = "second") 
PCHiC_Javierre_enhancergenes = PGRS_variants_GR[queryHits(hits)] 
PCHiC_vars_in_PGRS_variants = PCHiC[subjectHits(hits)] 
mcols(PCHiC_Javierre_enhancergenes) <- cbind(mcols(PCHiC_Javierre_enhancergenes), mcols(PCHiC_vars_in_PGRS_variants))

mcols(PCHiC_Javierre_enhancergenes) = 
  separate( as_tibble(mcols(PCHiC_Javierre_enhancergenes)),
            col = "gene_id",
            into = c("ENSG_PCHiC", "ENSG_PCHiC_vrs")
  )

saveRDS(
  PCHiC_Javierre_enhancergenes,
  file = file.path( PGRS_data_path, paste("gene_annot", annotation_date, "PCHiC_Javierre_enhancergenes", PGRS_suffix, sep = "_") )
)



