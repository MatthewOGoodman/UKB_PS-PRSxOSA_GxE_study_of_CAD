# Save KEGG pathway gene lists


# Using KEGGREST package
###

# Get gene lists for curated pathways ---- 
###

# PSPRS_KEGG_PATHWAYS (set in run_pathways.R as kegg_mode):
#   "curated" - the curated KEGG pathways of the pathway sheet, plus all human pathways (needed
#               for combined_pathway and the gene-level sets); the method used for the paper
#   "all"     - all human pathways
#   "hsa04066,hsa04370" - the listed pathways only
# KEGG serves only its current release, so gene sets differ from the 2020-12-09 lists the paper used.

if(kegg_mode == "curated"){
# curated_KEGG_pathways loaded from load_curated_pathways.R
KEGG_curated.pathways_genes = get_genes_in_KEGGREST_pathway( slct_KEGG_pathways = curated_KEGG_pathways)

KEGG_curated.pathways_ENSG_list = KEGG_curated.pathways_genes$ENSG
KEGG_curated.pathways_HGNC_list = KEGG_curated.pathways_genes$HGNC

saveRDS(KEGG_curated.pathways_ENSG_list, file = file.path(pathways_path, paste("KEGG_curated_pathways_", pathways_freeze_date, "_ENSG_list.rds", sep = "")) )
saveRDS(KEGG_curated.pathways_HGNC_list, file = file.path(pathways_path, paste("KEGG_curated_pathways_", pathways_freeze_date, "_HGNC_list.rds", sep = "")) )
}




# Get gene lists for all pathways in database ----
###

# listDatabases()
# ?KEGGREST::keggList

if(kegg_mode %in% c("curated", "all")){
all_KEGG_pathways = KEGGREST::keggList(database = "pathway", organism = "hsa")
all_KEGG_pathways = gsub("path:", "", names(all_KEGG_pathways))
} else {
all_KEGG_pathways = trimws(strsplit(kegg_mode, ",")[[1]])
}

KEGG_all.pathways_genes = get_genes_in_KEGGREST_pathway( slct_KEGG_pathways = all_KEGG_pathways)

KEGG_all.pathways_ENSG_list = KEGG_all.pathways_genes$ENSG
KEGG_all.pathways_HGNC_list = KEGG_all.pathways_genes$HGNC

saveRDS(KEGG_all.pathways_ENSG_list, file = file.path(pathways_path, paste("KEGG_all_hsa_pathways_", annotation_date, "_ENSG_list.rds", sep = "")) )
saveRDS(KEGG_all.pathways_HGNC_list, file = file.path(pathways_path, paste("KEGG_all_hsa_pathways_", annotation_date, "_HGNC_list.rds", sep = "")) )

dir(pathways_path)

  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Get curated gene-level "pathways" ----
###
# Each PSPRS_GENE_SETS group splits two KEGG pathways into single-gene sets: genes in only the
# first, only the second, or both. The paper used both groups (two runs).
gene_set_groups = list(HIF1_VEGF = c("hsa04066", "hsa04370"), NFkB_TNF = c("hsa04064", "hsa04668"))
gene_set_groups = gene_set_groups[intersect(names(gene_set_groups), gene_sets_requested)]
kegg_ids_present = sub(":.*", "", names(KEGG_all.pathways_ENSG_list))
gene_set_groups = gene_set_groups[vapply(gene_set_groups, function(ids) all(ids %in% kegg_ids_present), logical(1))]
if(length(gene_set_groups) > 0){
# a pathway's genes, looked up by KEGG id (names change between KEGG releases); NULL if absent
kegg_get = function(lst, id) { i = grep(paste0("^", id, ":"), names(lst)); if (length(i)) lst[[i]] else NULL }

KEGG_pathwaysENSG_list = KEGG_all.pathways_ENSG_list
KEGG_pathwaysHGNC_list = KEGG_all.pathways_HGNC_list

HIF1_genes = kegg_get(KEGG_pathwaysENSG_list, "hsa04066")
names(HIF1_genes) = kegg_get(KEGG_pathwaysHGNC_list, "hsa04066")

VEGF_genes = kegg_get(KEGG_pathwaysENSG_list, "hsa04370")
names(VEGF_genes) = kegg_get(KEGG_pathwaysHGNC_list, "hsa04370")

NFkB_genes = kegg_get(KEGG_pathwaysENSG_list, "hsa04064")
names(NFkB_genes) = kegg_get(KEGG_pathwaysHGNC_list, "hsa04064")

TNF_genes = kegg_get(KEGG_pathwaysENSG_list, "hsa04668")
names(TNF_genes) = kegg_get(KEGG_pathwaysHGNC_list, "hsa04668")

length(HIF1_genes); length(VEGF_genes); length(NFkB_genes); length(TNF_genes)

HIF1_only_genes = setdiff(HIF1_genes, VEGF_genes)
names(HIF1_only_genes) = names(HIF1_genes)[match(HIF1_only_genes, HIF1_genes)]

VEGF_only_genes = setdiff(VEGF_genes, HIF1_genes)
names(VEGF_only_genes) = names(VEGF_genes)[match(VEGF_only_genes, VEGF_genes)]

HIF1_VEGF_genes = intersect(HIF1_genes, VEGF_genes)
names(HIF1_VEGF_genes) = names(HIF1_genes)[match(HIF1_VEGF_genes, HIF1_genes)]

###

NFkB_only_genes = setdiff(NFkB_genes, TNF_genes)
names(NFkB_only_genes) = names(NFkB_genes)[match(NFkB_only_genes, NFkB_genes)]

TNF_only_genes = setdiff(TNF_genes, NFkB_genes)
names(TNF_only_genes) = names(TNF_genes)[match(TNF_only_genes, TNF_genes)]

NFkB_TNF_genes = intersect(NFkB_genes, TNF_genes)
names(NFkB_TNF_genes) = names(NFkB_genes)[match(NFkB_TNF_genes, NFkB_genes)]

length(NFkB_only_genes); length(TNF_only_genes); length(NFkB_TNF_genes); 


genes_of_interest = c(
  if("HIF1_VEGF" %in% names(gene_set_groups)) list(
  HIF1_only_genes = as.list(HIF1_only_genes),
  VEGF_only_genes = as.list(VEGF_only_genes),
  HIF1_VEGF_genes = as.list(HIF1_VEGF_genes)),
  if("NFkB_TNF" %in% names(gene_set_groups)) list(
  NFkB_only_genes = as.list(NFkB_only_genes),
  TNF_only_genes  = as.list(TNF_only_genes),
  NFkB_TNF_genes  = as.list(NFkB_TNF_genes)),
  NULL
)
genes_of_interest = unlist(genes_of_interest, recursive = FALSE)
genes_name = paste ( unique(sapply(str_split(names(genes_of_interest), "\\_"), "[[", 1)), collapse = "_")
saveRDS( genes_of_interest, file = file.path(pathways_path, paste("ENSG_", genes_name, "_KEGGpath_Ensembl_genes_", annotation_date, "_ENSG_list.rds", sep = ""))  )
} else {
  genes_of_interest = list()
  genes_name = NULL
}


# Get all annotated genes ----
###
all_ENSGs = unlist( read.csv( file = file.path(PGRS_data_path, paste("gene_annot", annotation_date, "all_ENSGs_annotated.txt", sep = "_")), stringsAsFactors = FALSE ) )

# path_gencode and GenCode.db come from config/paths.R (PSPRS_GENCODE_GFF)
GenCode.txdb.hg19.track <- rtracklayer::import(file.path(path_gencode, GenCode.db))

all_HGNCs = geneID_translate_fn (  gene_vec = all_ENSGs, 
                                   gene_in_type = c("ensg_id"), 
                                   gene_out_type = c("gene_name"), 
                                   database = GenCode.txdb.hg19.track)
all_genes_ENSG = all_ENSGs
names(all_genes_ENSG) = all_HGNCs
all_genes_ENSG = as.list(all_genes_ENSG)

saveRDS( all_genes_ENSG, file = file.path(pathways_path, paste("ENSG_annotPGRS_Ensembl_genes_", annotation_date, "_ENSG_list.rds", sep = "")) )


