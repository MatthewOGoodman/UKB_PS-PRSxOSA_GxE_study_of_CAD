# annotatePathways.R

# input: PGRS for prediction and PGRS sequentially annotated
# output: PGRS for prediction with new columns indicating SNP inclusion in each pathway (as TRUE/FALSE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# load data ----
####

# PGRS to add pathway columns onto (the PRS, one row per variant, loaded by run_pathways.R)
if(!exists("PGRS_MAF")) PGRS_MAF = readRDS(file = file.path(PGRS_data_path, PGRS_MAF_infile))

# SNP-to-gene annotations data
sequential_ann = readRDS(file = file.path(PGRS_data_path, PGRS_gene_annot_infile))
sequential_ann = sortSeqlevels(sequential_ann) ## seqlevels(sequential_ann) 
sequential_ann = sort(sequential_ann)
# limit to select.annotation.sources
sequential_ann = sequential_ann[sequential_ann$source %in% select.annotation.sources]


# get pathways lists ----
####

# PS-PGRS pathways
KEGG_allPaths_ENSG_list = readRDS(file = file.path(pathways_path, paste("KEGG_all_hsa_pathways_", annotation_date, "_ENSG_list.rds", sep = "")) )
KEGG_allPaths_HGNC_list = readRDS(file = file.path(pathways_path, paste("KEGG_all_hsa_pathways_", annotation_date, "_HGNC_list.rds", sep = "")) )

if(kegg_mode == "curated"){
KEGG_slctPaths_ENSG_list = readRDS(file = file.path(pathways_path, paste("KEGG_curated_pathways_", pathways_freeze_date, "_ENSG_list.rds", sep = "")) )
KEGG_slctPaths_HGNC_list = readRDS(file = file.path(pathways_path, paste("KEGG_curated_pathways_", pathways_freeze_date, "_HGNC_list.rds", sep = "")) )

# add curated gene-set, combined curated genes and pathways gene-set, and all genes, to KEGG pathways
###
curated_ENSGs   = readRDS(file = file.path(pathways_path, paste(paste("curated_ENSG", curated_pathways_name, pathways_freeze_date, sep = "_"), ".rds", sep = "")) )
KEGG_slctPaths_ENSG_list[["all_annotated_genes"]] = unique(sequential_ann$ENSG)
KEGG_slctPaths_ENSG_list[[curated_pathways_name]] = curated_ENSGs
# make a combined pathway out of all of the curated pathways plus the curated genes
slct_paths = paste( c(curated_KEGG_pathways, curated_pathways_name ), collapse = "|")
# make sure that combined_pathway below is not defined recursively accidentally:
combined_curated_KEGGpathways_ENSGgenes =  unique( unlist( KEGG_allPaths_ENSG_list[ grepl(slct_paths, names(KEGG_allPaths_ENSG_list) )] ))
KEGG_slctPaths_ENSG_list[["combined_pathway"]] = combined_curated_KEGGpathways_ENSGgenes
} else {
# "all" or a list of pathways: those pathways, plus all annotated genes
KEGG_slctPaths_ENSG_list = KEGG_allPaths_ENSG_list
KEGG_slctPaths_ENSG_list[["all_annotated_genes"]] = unique(sequential_ann$ENSG)
}

# GS-PGRS pathways
GSPGRS_allPGRS_ENSG_list = readRDS( file = file.path(pathways_path, paste("ENSG_annotPGRS_Ensembl_genes_", annotation_date, "_ENSG_list.rds", sep = "")) )
GSPGRS_slctPath_ENSG_list = if(length(genes_of_interest) > 0) readRDS( file = file.path(pathways_path, paste("ENSG_", genes_name, "_KEGGpath_Ensembl_genes_", annotation_date, "_ENSG_list.rds", sep = "")) ) else list()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# annotate PGRS with pathways ----

# Function that adds columns to Khera_CAD_PGRS: T/F for inclusion of SNPs based on genes in KEGG_pathwaysENSG_list
annotate_PGRS_pathways = function(pathways_list = KEGG_pathwaysENSG_list, 
                                  PGRS_multi_annotated = sequential_ann,
                                  PGRS = PGRS_MAF, 
                                  verbose_TF = FALSE){
  # pathway_name = "curated_genes"
  # pathways_list = KEGG_pathwaysENSG_list
  # PGRS = PGRS_MAF
  
  KEGGpathENSGs_in_PGRS = list()
  
  
  for(pathway_name in names(pathways_list)){
    pathway_ENSGs = SNPs_in_path = NULL
    rows_in_path = FALSE
    
    pathway_ENSGs = pathways_list[[pathway_name]] 
    pathway_ENSGs = unlist(unique(pathway_ENSGs))
    
    if(verbose_TF){
      print(pathway_name)
      print(pathway_ENSGs)
    }
    
    SNPs_in_path = PGRS_multi_annotated$variant[  PGRS_multi_annotated$ENSG %in% pathway_ENSGs ]
    SNPs_in_path = sort(unique(SNPs_in_path))
    
    ENSGs_in_path = pathway_ENSGs %in% unique(PGRS_multi_annotated$ENSG)
    names(ENSGs_in_path) = pathway_ENSGs
    KEGGpathENSGs_in_PGRS[[pathway_name]] = ENSGs_in_path
    
    mcols(PGRS)[pathway_name] = FALSE
    rows_in_path = mcols(PGRS)[["variant"]] %in% SNPs_in_path
    mcols(PGRS)[[pathway_name]] = rows_in_path
    
    # print checks:
    n_variants = sum(mcols(PGRS)[[pathway_name]])
    print(sprintf("There are %i variants in the pathway: %s", n_variants, pathway_name))
    if(verbose_TF){
      print(sprintf( "Sampling of 10 variants from pathway: %s ", pathway_name))
      print( PGRS[sample(which(rows_in_path),10),] )
    }
    
  }
  
  return(list(PGRS_paths = PGRS, KEGGpathENSGs_in_PGRS = KEGGpathENSGs_in_PGRS))
}
