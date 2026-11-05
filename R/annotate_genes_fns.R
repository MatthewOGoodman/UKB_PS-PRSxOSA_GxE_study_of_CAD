# Annotate PGRS from GRanges


get_PGRS_annotations_from_GRanges = function(PGRS_variants_GR, annotations_GR){
  
  # Find PGRS variants in Annotation data:
  hits = findOverlaps(query = PGRS_variants_GR, subject = annotations_GR)
  
  # combine annotation columns as PGRS_annotated:
  PGRS_annotated = 
    PGRS_variants_GR[queryHits(hits)] 
  
  annotations_in_PGRS_variants = annotations_GR[subjectHits(hits)] 
  mcols(PGRS_annotated) = 
    cbind(mcols(PGRS_annotated), 
          mcols(annotations_in_PGRS_variants))
  
  return(PGRS_annotated)
}



# all((as_tibble(granges( annotations_GR[subjectHits(hits)] ))$start <= as_tibble(PGRS_variants_GR[queryHits(hits)])$start & 
       # as_tibble(granges( annotations_GR[subjectHits(hits)] ))$end >= as_tibble(PGRS_variants_GR[queryHits(hits)])$start ) )




# Find PGRS variants in TxDb data:
####

get_PGRS_annotations_from_txdb = function(PGRS_variants_GR, 
                                          annotations_TxDb, 
                                          
                                          AllVariants_TF = TRUE,
                                          ID_type = "gene",
                                          intgen_upstream = 1e6,
                                          intgen_downstream = 1e6,
                                          
                                          promo_upstream = 2000L,
                                          promo_downstream = 200L
                                          ){
  
  require(VariantAnnotation)
  
  if(AllVariants_TF){
    PGRS_anno_TxDb <- VariantAnnotation::locateVariants(query = PGRS_variants_GR, 
                                                        subject = annotations_TxDb, 
                                                        region = AllVariants(
                                                          promoter = PromoterVariants(
                                                            upstream   = promo_upstream,    # 2000L default
                                                            downstream = promo_downstream   # 200L default
                                                          ),
                                                          intergenic = IntergenicVariants( 
                                                            # for intergenic variants, 'GENEID' is NA
                                                            # In order to return a list of gene/tx names in PRECEDEID and FOLLOWID, 
                                                            # need to set upstream downstream and ID type ("gene" vs "tx")
                                                            upstream   = intgen_upstream,    # 1e6 default
                                                            downstream = intgen_downstream,  # 1e6 default
                                                            idType= ID_type) # "gene" default, "tx" alternative
                                                        ),
                                                        ignore.strand = FALSE #default
    )
  } else {
    PGRS_anno_TxDb <- VariantAnnotation::locateVariants(query = PGRS_variants_GR, 
                                                        subject = annotations_TxDb, 
                                                        region = CodingVariants(),
                                                        ignore.strand = FALSE #default
    )
  }



  
  # combine GenCode and PGRS metadata annotation columns:
  ####
  
  # expand PGRS metadata columns to new many-to-one setting:
  PGRS_metadata <- PGRS_variants_GR[PGRS_anno_TxDb$QUERYID]
  # Note: QUERYID is returned by locateVariants and is row location in original query GR
  
  # cbind results:
  mcols(PGRS_anno_TxDb) <- cbind(
    mcols(PGRS_metadata)[c("variant", "wgt", "effect_allele", "A1", "A2")],
    mcols(PGRS_anno_TxDb)
  )
  rm(PGRS_metadata)
  
  return (PGRS_anno_TxDb)
}


# Obtain ENSG ID's from HGNC ID's by hand
HGNC_to_ENSG = function(HGNC_vec, GenCode.txdb.hg19.track){
  
  # Create lookup table
  GenCode.txdb.hg19.lookup = distinct( as_tibble(mcols(GenCode.txdb.hg19.track))[ c("gene_id", "gene_name", "hgnc_id") ])
  GenCode.txdb.hg19.lookup = separate(GenCode.txdb.hg19.lookup,
                                      gene_id,
                                      into = c("ENSG", "ENSG_vrs"),
                                      sep = "\\.")
  
  # match gene names
  matches = match(x = HGNC_vec , table = GenCode.txdb.hg19.lookup$gene_name, nomatch = NA_character_, incomparables = NULL)
  ENSG_vec = GenCode.txdb.hg19.lookup$ENSG[matches]
  
  return(ENSG_vec)
}

# Obtain ENSG ID's from HGNC ID's using AnnotationDb function MapIds
HGNC_to_ENSG_AnnotationDb.mapIds = function(HGNC_vec, org.Hs.eg.db, mVals = "first"){
  ENSG_vec = 
    AnnotationDbi::mapIds(org.Hs.eg.db,
           keys= HGNC_vec ,
           keytype= "SYMBOL",
           column= "ENSEMBL",
           multiVals= mVals)
}


# Obtain ENSG ID's from HGNC ID's by hand
geneID_translate_fn = function(  gene_vec, 
                          gene_in_type = c("ensg_id", "gene_name", "hgnc_id"), 
                          gene_out_type = c("ensg_id", "gene_name", "hgnc_id"), 
                          database = GenCode.txdb.hg19.track
                          ){
  
  # Create lookup table
  GenCode.txdb.hg19.lookup = distinct( as_tibble(mcols(database))[ c("gene_id", "gene_name", "hgnc_id") ])
  GenCode.txdb.hg19.lookup = separate(GenCode.txdb.hg19.lookup,
                                      "gene_id",
                                      into = c("ensg_id", "ensg_vrs"),
                                      sep = "\\.")
  
  print(dim(GenCode.txdb.hg19.lookup))
  print(class(GenCode.txdb.hg19.lookup))
  
  # match gene names
  matches = match(x = gene_vec , table = GenCode.txdb.hg19.lookup[[gene_in_type]], nomatch = NA_character_, incomparables = NULL)
  out_vec = (GenCode.txdb.hg19.lookup[[gene_out_type]])[matches]
  
  return(out_vec)
}
