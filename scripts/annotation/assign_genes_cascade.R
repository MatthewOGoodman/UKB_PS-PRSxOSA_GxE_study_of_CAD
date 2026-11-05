#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Combine annotion sources sequentiallly orderd by perceived specificity
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###
# Final annotation dataframe constructed so that at each stage 
# the only SNPs left to annotate are the ones not annotated by the "more specific" prior methods.

# choose the sequential order: use names of "combined_annot_list"
ordered_annotations = c(
  "PGRS_annot_GENCODE_coding_exon", 
  "PGRS_annot_FANTOM5",
  "PGRS_annot_PCHiC", 
  "PGRS_annot_GTEx",
  "PGRS_annot_GeneHancer_elite", 
  "PGRS_annot_GENCODE_coding_UTR", 
  # "PGRS_annot_GeneHancer_non.elite", will be removed for PS-PGRS calculation. Will leave holes in flanked genes if used here
  "PGRS_annot_GENCODE_flanked_genes", 
  "PGRS_annot_GENCODE_noncoding")

PGRS_suffix = paste(PGRS_out_name, "rds", sep = ".")

# the PRS (one row per variant, loaded by run_annotation.R) is the scaffold that records which
# variants have been annotated
if(!exists("PGRS_MAF")){
  PGRS_MAF = PGRS_variants_GR
}
if(!exists("combined_annot_list")){
  combined_annot_list = readRDS( file = file.path( PGRS_data_path, paste("gene_annot", annotation_date, "combined_annotation_sources", PGRS_suffix, sep = "_" )) )
}
# sources absent from this build (GeneHancer when PSPRS_GENEHANCER_CSV is unset) drop out of the order
ordered_annotations = ordered_annotations[ordered_annotations %in% names(combined_annot_list)]

aggregated_ann = NULL
aggregated_var = NULL
PGRS_MAF$annotated_TF = FALSE
annotation_cumsum = c()
aggregated_ann_list_distinct = list()
for( ann_name in ordered_annotations ){
  short_name = paste(str_split(ann_name, pattern = "\\_")[[1]][-(1:2)], collapse = "_")
  annDB = combined_annot_list[[ann_name]]
  print(ann_name); print(length(unique(annDB$variant)))
  
  # aggregated_ann = append(aggregated_ann, annDB)
  # update list prior to updating variant
  aggregated_ann_list_distinct[[short_name]] = annDB[!annDB$variant %in% aggregated_var,]
  print("# new ann list variants"); print(length(unique(aggregated_ann_list_distinct[[short_name]]$variant)))
  print(length(setdiff(union(aggregated_var, annDB$variant), aggregated_var)))
  
  aggregated_var = union(aggregated_var, annDB$variant)
  PGRS_MAF$annotated_TF[PGRS_MAF$variant %in% aggregated_var] = TRUE

  annotation_cumsum[short_name] = sum(PGRS_MAF$annotated_TF)
  
}


aggregated_ann = unlist(GRangesList(aggregated_ann_list_distinct))
aggregated_ann = sortSeqlevels(aggregated_ann)
aggregated_ann = sort(aggregated_ann)

rm(aggregated_ann_list_distinct)

# checks
# lapply(aggregated_ann_list_distinct, function(x) names(mcols(x)))
# class(aggregated_ann)
summary( as.factor(aggregated_ann$source))


# plot cumulative annotation results:
####
annotation_percent = annotation_cumsum / length(PGRS_MAF) * 100

hist.cum.annot = ggplot(data =  as_tibble(data.frame("rows" = 1:length(annotation_percent), "database" = names(annotation_percent), "tot_vars" = annotation_percent)),
       aes( x = reorder(database, rows), y = tot_vars )) +
  theme(text = element_text(size=12.5),
        axis.text.x = element_text(angle=30, hjust=0.85)) +
  # aes(x = reorder(Var1, -Freq), y = Freq)) +
  geom_bar(stat = "identity") +
  xlab("Source Database of PGRS annotation") +
  ylab("Cumulative percent of PGRS variants annotated")

pdf(file.path(PGRS_data_path, paste("gene_annot", annotation_date, "sources_cumulative.pdf", sep = "_")),
    width=10, height=6)
print(hist.cum.annot)
dev.off()


annotations_used = gsub( "PGRS_annot_", "", ordered_annotations )
write.csv(annotations_used, file = file.path(PGRS_data_path, paste("gene_annot", annotation_date, "annotations_used.txt", sep = "_")) )
all_ENSGs = unique(mcols(aggregated_ann)$ENSG)
write.csv(all_ENSGs, file = file.path(PGRS_data_path, paste("gene_annot", annotation_date, "all_ENSGs_annotated.txt", sep = "_")) )

saveRDS(aggregated_ann, file = file.path( PGRS_data_path, paste("gene_annot", annotation_date, "sequential.annotated.SNPs.to.ENSGs", PGRS_suffix, sep = "_"))  )

rm(combined_annot_list, aggregated_var)




