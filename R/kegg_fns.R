# fns_KEGG.R

library(KEGGgraph)
# https://rdrr.io/bioc/KEGGgraph/f/inst/doc/KEGGgraph.pdf
# note that you can graph these networks: see docs

# get genes from KEGG pathway using KEGGgraph xml: can return graphs ----
####
get_genes_in_KEGGgraph_pathway = function(KEGGpathway, file_path_KEGG = path_KEGG, plot_TF = FALSE){

  # check there is a corresponding url
  getKGMLurl(pathwayid= KEGGpathway, organism = "hsa")
  
  KEGGfile = file.path(file_path_KEGG, paste0(KEGGpathway,".xml"))
  
  if(!file.exists(KEGGfile)){
    retrieveKGML(pathwayid = KEGGpathway, organism="hsa", destfile=KEGGfile, method="auto", quiet=FALSE)
  }
  # retrieveKGML reports a failed download only as a message, so check for the file
  if(!file.exists(KEGGfile) || file.size(KEGGfile) == 0){
    unlink(KEGGfile)
    stop("KGML download failed for ", KEGGpathway, " (", KEGGfile, ")")
  }
  
  KEGGgraph <- parseKGML2Graph(KEGGfile,expandGenes=TRUE) 
  # expand genes means translate paralogues i.e. "biological gene concepts" 
  # to specific transcribed genes/gene-products
  KEGGnodes = getKEGGnodeData(KEGGgraph)
  KEGGgenes = names(KEGGnodes)
  KEGG_EntrezIDs = translateKEGGID2GeneID(KEGGgenes)
  
  # translate to ENSG 
  library(org.Hs.eg.db)
  KEGG_ENSG = sapply(mget(KEGG_EntrezIDs, org.Hs.egENSEMBL, ifnotfound=NA), "[[",1)
  KEGG_ENSG[is.na(KEGG_ENSG)] = names(KEGG_ENSG)[is.na(KEGG_ENSG)]  # set back as KEGG names if not found
  KEGG_HGNC = sapply(mget(KEGG_EntrezIDs, org.Hs.egSYMBOL, ifnotfound=NA), "[[",1)
  KEGG_HGNC[is.na(KEGG_HGNC)] = names(KEGG_HGNC)[is.na(KEGG_HGNC)]  # set back as KEGG names if not found
  
  # graph node labels must be unique; two KEGG nodes can carry the same gene symbol
  nodes(KEGGgraph) = make.unique(KEGG_HGNC)
    if(plot_TF){
      plot(KEGGgraph, "neato",attrs=list(node=list(font=5, fillcolor="lightblue")))
  }

  return(list(graph = KEGGgraph, ENSG = KEGG_ENSG, HGNC = KEGG_HGNC))
  
  # ccan also look at pathway, not gene interaction network
  # (not needed, but contains additional information)
  # mapkPathway <- parseKGML(KEGGfile)
  # mapkPathway
  # nodes(mapkPathway)
  # note that "Type" includes orthologs and compounds as well as genes

}

# get genes from KEGG pathway: using KEGGREST database ----
library(KEGGREST)

# KGML cache directory for the current KEGG release, so that files downloaded from different
# KEGG releases are never mixed. keggInfo("kegg") reports the date the pathway database was last
# updated (a line such as "pathway  <entries>  YYYY/MM/DD"); that date names the subdirectory.
kegg_release_dir = function(base = path_KEGG){
  info = paste(KEGGREST::keggInfo("kegg"), collapse = "\n")
  hit = regmatches(info, regexec("pathway\\s+[0-9]+\\s+([0-9]{4})/([0-9]{2})/([0-9]{2})", info))[[1]]
  rel = if (length(hit) == 4) paste0("KEGG_pathway_", hit[2], hit[3], hit[4]) else paste0("KEGG_pathway_unknown_", format(Sys.Date(), "%Y%m%d"))
  d = file.path(base, rel)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# note this function is actually a wrapper for get_genes_in_KEGGgraph_pathway but uses KEGGrest to name the pathways
get_genes_in_KEGGREST_pathway = function( slct_KEGG_pathways = "all_KEGG_pathways", file_path_KEGG = kegg_release_dir()){
  if (identical(slct_KEGG_pathways, "all")) slct_KEGG_pathways = "all_KEGG_pathways"
  slct_TF = length(slct_KEGG_pathways) >= 1 & all(grepl("hsa", slct_KEGG_pathways))
  all_KEGG_TF = length(slct_KEGG_pathways) == 1 & slct_KEGG_pathways[1] == "all_KEGG_pathways"
  
  if(!(slct_TF|all_KEGG_TF)) warning("need specify KEGG pathways: e.g. hsa00001, or use character string: all_KEGG_pathways")
  if ( all_KEGG_TF ) {
    print( "using all human/hsa pathways" )
    # listDatabases()
    all_KEGG_pathways = KEGGREST::keggList(database = "pathway", organism = "hsa")
    KEGG_pathways = gsub("path:", "", names(all_KEGG_pathways))
  } else if(slct_TF){
    print( "using selected pathways" )
    KEGG_pathways = slct_KEGG_pathways
  }
  
  print("getting pathway names")
  # names come from the one list of human pathways (the KEGG API no longer lists a single pathway)
  hsa_pathway_names = KEGGREST::keggList(database = "pathway", organism = "hsa")
  names(hsa_pathway_names) = gsub("path:", "", names(hsa_pathway_names))
  unknown_pathways = setdiff(KEGG_pathways, names(hsa_pathway_names))
  if (length(unknown_pathways)) stop("not human KEGG pathways in the current release: ", paste(unknown_pathways, collapse = ", "))
  KEGG_path_longnames = hsa_pathway_names[KEGG_pathways]
  KEGG_path_shortnames = sapply(str_split(KEGG_path_longnames, "\\s-\\s"), "[[", 1)
  KEGG_path_shortnames = paste0(KEGG_pathways, ": ", KEGG_path_shortnames)
  names(KEGG_pathways) = KEGG_path_shortnames
  
  print("getting pathway genes")
  KEGG_pathways_list = lapply(KEGG_pathways, get_genes_in_KEGGgraph_pathway, file_path_KEGG = file_path_KEGG)
  
  KEGG_pathwaysHGNC_list = lapply(KEGG_pathways_list, '[[', "HGNC")
  names(KEGG_pathwaysHGNC_list) = KEGG_path_shortnames
  
  KEGG_pathwaysENSG_list = lapply(KEGG_pathways_list, '[[', "ENSG")
  names(KEGG_pathwaysENSG_list) = KEGG_path_shortnames
  
  return(list(ENSG = KEGG_pathwaysENSG_list, HGNC = KEGG_pathwaysHGNC_list))
}



# find pathways from genes
# query_gene = "HIF1A"
get_KEGG_genes_from_HGNC = function(query_gene, hsa_db = hsa){
    query_gene = paste0("(^|\\s)",query_gene,"(\\,|\\;)")
  idx = grep(query_gene, hsa_db)
  # hsa_db[idx]
  KEGG_gene = hsa_db[idx]
  print("gene info:")
  print(hsa_db[idx])
  return(KEGG_gene)
}
get_KEGG_pathways_from_KEGG_gene = function(KEGG_gene){
  res = keggLink("pathway", KEGG_gene)
  return(res)
}
