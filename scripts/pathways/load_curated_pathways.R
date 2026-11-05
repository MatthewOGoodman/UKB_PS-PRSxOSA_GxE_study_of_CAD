# loadCuratedPathways.R

# load curated genes and pathways spreadsheet ----
####

# the curated pathway sheet (default: data/OSA_Pathways_20201117_freeze.csv; config/paths.R)
pathways_mat = read_csv(curated_pathways_csv)
rows_all_missing = apply(pathways_mat, 1, function(x) all(is.na(x)))
pathways_mat = pathways_mat[!rows_all_missing, ]
dim(pathways_mat)

# get pre-curated genes ----
####
gene_idx = grep("Genes", colnames(pathways_mat))
gene_idxs = gene_idx:ncol(pathways_mat)
include = pathways_mat$IncludeTF == "TRUE"

curated_CAD_IH_genes =
  unique_from_columns(
    dataframe = pathways_mat, 
    rows = include, 
    columns = gene_idxs, 
    exclude = "?")

# prep curated genes
####

# Convert CAD_IH gene ID's from HGNC gene names to ENSG symbols
# TODO:  better function to convert gene symbols in genomic data manipulation section
library(org.Hs.eg.db)
curated_CAD_IH_genes_Entrez = sapply(mget(curated_CAD_IH_genes, org.Hs.egSYMBOL2EG, ifnotfound=NA), "[[",1)
curated_CAD_IH_genes_Entrez[is.na(curated_CAD_IH_genes_Entrez)]
curated_CAD_IH_genes_Entrez = curated_CAD_IH_genes_Entrez[!is.na(curated_CAD_IH_genes_Entrez)]
curated_CAD_IH_genes_ENSG = sapply(mget(curated_CAD_IH_genes_Entrez, org.Hs.egENSEMBL, ifnotfound=NA), "[[",1)
names(curated_CAD_IH_genes_ENSG) = names(curated_CAD_IH_genes_Entrez)
length(curated_CAD_IH_genes_ENSG) 

# save curated genes
####

saveRDS(curated_CAD_IH_genes_ENSG, file = file.path(pathways_path, paste(paste("curated_ENSG", curated_pathways_name, pathways_freeze_date, sep = "_"), ".rds", sep = "")))
dir(pathways_path)


# get curated KEGG pathways ----
####

curated_KEGG_pathways = 
  unique_from_columns(
    dataframe = pathways_mat, 
    rows = include, 
    columns = "KEGGmap", 
    exclude = "?")
curated_KEGG_pathways = sort(curated_KEGG_pathways)
# make sure that the e.g. "hsa04330" strings do not contain spaces entered in the spreadsheet by accident:
curated_KEGG_pathways = str_trim(curated_KEGG_pathways)





