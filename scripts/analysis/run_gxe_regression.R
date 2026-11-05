#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Scripts ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Run from the repository root:  Rscript scripts/analysis/run_gxe_regression.R
# Required environment: PSPRS_ANALYSIS_DATA_DIR (the analysis data set written by
# prepare_analysis_data.R; UK Biobank, approved access only).
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(mgcv) })
source("config/paths.R")
psprs_require("UKB_data_path")
source(file.path("scripts", "analysis", "curated_modules.R"))
source(file.path("R", "model_fit_fns.R"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# load results
UKB_pheno_PGRS = readRDS(file = file.path(UKB_data_path, "UKB_pheno_PGRS_38.psPGRS_HIF1.VEGF.NFkB.TNF.genes_analytic_20210407") )


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Pathways ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Check: core_genes, module_genes (heirarchical lists) and select.pathways.vec loaded from "4.0_curated_pathways_core.genes_module.genes.R"
####

select.pathways.vec
core_genes
module_genes

curated.gene.modules = c(names(core_genes), names(module_genes))
sum(!curated.gene.modules %in% colnames(UKB_pheno_PGRS))

# get pathway names from data columns
####

pathways.all      =  grepl("psPGRS_", colnames(UKB_pheno_PGRS) ) 
combined.pathways =  c( "psPGRS_all_SNPs",               # original PGRS
                        "psPGRS_all_annotated_genes",    # subset of PGRS annotated to genes
                        "psPGRS_combined_pathway",       # curated OSA/PGRS genes and their associated KEGG pathways as recored in the OSA pathways spreadsheet
                        "psPGRS_OSA_pathways" )          # curated OSA/PGRS genes as recored in the OSA pathways spreadsheet
pathways.KEGG     =  colnames(UKB_pheno_PGRS)[grepl("psPGRS", colnames(UKB_pheno_PGRS)) & grepl("hsa[[:digit:]]", colnames(UKB_pheno_PGRS) )]
pathways.KEGG     =  pathways.KEGG[!grepl("\\.1", pathways.KEGG)]

# get individual gene pathways
gene.pathways     =  colnames(UKB_pheno_PGRS)[ grepl("psPGRS_", colnames(UKB_pheno_PGRS) ) &  grepl("_genes_", colnames(UKB_pheno_PGRS) ) ]
gene.pathways     =  gene.pathways[!grepl("\\.[[:digit:]]|upstream|downstream|adjacent|core", gene.pathways)]
gene_PGRS_sums = colSums(UKB_pheno_PGRS[,gene.pathways])
valid_genes = !is.na(gene_PGRS_sums)
gene.pathways = gene.pathways[valid_genes]
# UKB_pheno_PGRS[,unique.gene.pathways[!valid_genes]]
# missing_genes = unique.gene.pathways[!valid_genes]
# missing_genes = separate(data.frame(genes = missing_genes), genes, sep = "_", into = c("a","b","c", "d", "gene"))$gene
# paste(sort(missing_genes), collapse = ", ")
genes = separate(data.frame(genes = gene.pathways), genes, sep = "_", into = c("a","b","c", "d", "gene"))
dup_genes = duplicated(genes$gene)
# length(dup_genes); nrow(genes); length(gene.pathways)
# gene.pathways[grepl("IL6", gene.pathways)]
unique.gene.pathways = gene.pathways[!dup_genes]
# sum(is.na(colSums(UKB_pheno_PGRS[,unique.gene.pathways])))

# break out core genes as singletons
core_gene_singletons = unique(unlist(core_genes))
names(core_gene_singletons) = core_gene_singletons
which(is.na(colSums(UKB_pheno_PGRS[,core_gene_singletons])))
core_gene_singletons = core_gene_singletons[!is.na(colSums(UKB_pheno_PGRS[,core_gene_singletons]))]
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Define: select.pathways ####
####

# select.pathways = c(select.pathways.vec, names(core_genes))

# select.pathways = "psPGRS_all_SNPs"

select.pathways = c(
  "psPGRS_hsa04066_HIF_1_signaling_pathway" ,
  "psPGRS_hsa04370_VEGF_signaling_pathway",
  "psPGRS_hsa04064_NF_kappa_B_signaling_pathway",
  "psPGRS_hsa04668_TNF_signaling_pathway",
  # "psPGRS_hsa04010_MAPK_signaling_pathway",
  
  "psPGRS_core_genes_HIF1",
  "psPGRS_core_genes_VEGF",
  "psPGRS_core_genes_NFkB",
  # "psPGRS_core_genes_NFkB_canonical",
  # "psPGRS_core_genes_NFkB_nonCanonical",
  "psPGRS_core_genes_TNF",
  # "psPGRS_core_genes_TNF_all_cells",
  # "psPGRS_core_genes_TNF_WBCs",
  "psPGRS_all_SNPs",
  # "psPGRS_all_annotated_genes",
  "psPGRS_OSA_pathways"
)

slct.modules = c(
  "VEGF_downstream_migration",
  "HIF1_downstream_angiogenesis",
  "HIF1_VEGF_shared_PI3K.AKT_genes",
  "HIF1_upstream_PI3K.AKT_HIF1a.translation",
  "HIF1_upstream_IH_HIF1a.degredation",
  "VEGF_downstream_survival",
  "VEGF_downstream_permeability",
  "glucose_inhibit_citrate"
)

slct.genes = c(
  "psPGRS_HIF1_only_genes_IL6R",
  "psPGRS_TNF_only_genes_IRF1",
  "psPGRS_HIF1_VEGF_genes_PIK3R3",
  "psPGRS_VEGF_only_genes_MAPKAPK3",
  "psPGRS_HIF1_VEGF_genes_PIK3R3")

slct.genes_neg.GxE = 
c(
"psPGRS_NFkB_only_genes_RELB",
"psPGRS_HIF1_VEGF_genes_PIK3R3", # also in TNF pathway
"psPGRS_VEGF_only_genes_MAPKAPK3",
"psPGRS_TNF_only_genes_IRF1", # co-activates NO2 with HIF1
"psPGRS_HIF1_VEGF_genes_VEGFA",
"psPGRS_NFkB_only_genes_TIRAP",
"psPGRS_HIF1_VEGF_genes_PLCG1", # also in NFkB pathway
"psPGRS_HIF1_only_genes_FLT1",
"psPGRS_NFkB_TNF_genes_CFLAR",
"psPGRS_HIF1_VEGF_genes_PLCG2", # also in NFkB pathway
"psPGRS_TNF_only_genes_CREB3L2"
)

slct.genes_pos.GxE = c(
  "psPGRS_HIF1_only_genes_ARNT",
  "psPGRS_HIF1_VEGF_genes_MAPK1", #also in TNF pathway
  "psPGRS_NFkB_only_genes_PIAS4",
  "psPGRS_HIF1_only_genes_PDHB",
  "psPGRS_HIF1_VEGF_genes_MAP2K1", #also in TNF pathway
  "psPGRS_HIF1_only_genes_IL6R",
  "psPGRS_NFkB_only_genes_NFKB2"
)

# slct.genes = slct.genes
slct.genes = c(slct.genes_neg.GxE, slct.genes_pos.GxE)
names(slct.genes) = slct.genes

# sum(!slct.genes %in% colnames(UKB_pheno_PGRS))
# colnames(UKB_pheno_PGRS)[grep("MAP2K1", colnames(UKB_pheno_PGRS))]
# colnames(UKB_pheno_PGRS)[grep("NFKB2", colnames(UKB_pheno_PGRS))]
# UKB_pheno_PGRS$psPGRS_NFkB_TNF_genes_TRAF2
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Choose pathways for analysis ####
####

# pathways.slct.string = "primary.analysis.HIF1.VEGF.NFkB.TNFa.allSNPs.curated"
# pathways = select.pathways

# pathways.slct.string = "secondary.analysis.HIF1.VEGF.NFkB.TNFa.KEGG.genes"
# pathways = unique.gene.pathways

pathways.slct.string = "secondary.analysis.select.genes"
pathways = slct.genes

# pathways.slct.string = "HIF1.VEGF.NFkB.TNF.core.genes.and.modules"
# pathways = curated.gene.modules

# pathways.slct.string = "HIF1.VEGF.NFkB.TNFa.core.genes"
# pathways = c(names(core_genes), core_gene_singletons)

# pathways.slct.string = "HIF1.VEGF.NFkB.TNFa.slct.gam.genes"
# pathways = c( core_gene_singletons, slct.genes)

# pathways.slct.string = "HIF1.VEGF.NFkB.TNFa.curated.paths"
# pathways = paste("psPGRS", names(curated_modules_reduce), sep = "_")

# pathways.slct.string = "primary.analysis.allSNPs.for.Table.1"
# pathways = select.pathways

####
print(pathways)
####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Set Models: ####
####

# covariates of the primary model (released as model_covariates_primary.R)
model.file = file.path("scripts", "analysis", "model_covariates_primary.R")


# OSA exposures:  ####
# (least to most specific)
####
Exposure_Vars = c("cmp_SA.snoring")
# Exposure_Vars = "cmp_snore.daysleep"
# Exposure_Vars = "cmp_sleep.apnea.w.snoring" # aka snore.SA
# Exposure_Vars = "cmp_OSA"   # T/F version
# Exposure_Vars = "OSA"       # 0/1 version
# Exposure_Vars = "cmp_OSA_noHF.noStr"
# Exposure_Vars = c( "cmp_snore.daysleep", "cmp_SA.snoring", "cmp_SA")
# Exposure_Vars = sleep.sr_insomnia==3

# Outcomes: ####
####
Outcome_Vars = c("cmp_CAD")     # T/F version
# Outcome_Vars = "CAD"       # 0/1 version
# Outcome_Vars = "cmp_AA.MI"
# Outcome_Vars = "cmp_Khera.MI"
# Outcome_Varss = "cmp_stroke"
# Outcome_Vars = "cmp_atherosclerosis"


# Subject subsets: ####
#### 

# primary analysis:
####
subsets = c( "(cmp_SA.snoring | cmp_SA.snoring_controls)" )

# controls exclude participants with both snoring and daytime sleepiness


# sensitivity analyses:
####
# subsets = c(
#   "(cmp_OSA.cases | cmp_OSA.controls)",         # see prepare_analysis_data.R for definitions
#   "(cmp_OSA.cases | !cmp_OSA.cases)",            # Sensitivity Analysis: All subjects
#   "(cmp_OSA.cases | cmp_OSA.controls) & white",  # Sensitivity: whites only
#   "(cmp_OSA.cases | cmp_OSA.controls) & sexM",   # Sensitivity: sex-specific
#   "(cmp_OSA.cases | cmp_OSA.controls) & !sexM",  # Sensitivity: sex-specific
#   NULL
# )

# Earlier subsets
# subsets = c( "cmp_OSA_noHF.noStr | cmp_OSA.controls_noHF.noStr" )  # ATS
# subsets = c( "(cmp_OSA.cases | cmp_OSA.controls)")          # see prepare_analysis_data.R for definitions
# subsets = c( "!(cmp_stroke | cmp_heart_failure)" )
# subsets = c( "(cmp_snr.day.cases | cmp_snr.day.controls)" )


# Set model frame
####
model_frame = expand.grid("Y" = Outcome_Vars, "E" = Exposure_Vars, "G" = pathways, "subset" = subsets)

# model.frame_subset = 
# (grepl("cmp_daysleep|cmp_snore.daysleep", model_frame$E) & grepl("cmp_daysleep", model_frame$subset)) | 
#   (grepl("cmp_SA.snoring", model_frame$E) & grepl("cmp_SA.snoring", model_frame$subset)) | 
#   grepl("TRUE", model_frame$subset)
# model_frame = model_frame[model.frame_subset,]

model_frame; dim(model_frame)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


# Compute results: ####
####

# options
testing_TF = FALSE
convert.E.to.numeric_TF = TRUE

glm.TF = TRUE
gam.TF = TRUE

glm_list = list()
gam_list = list()

set.seed(12345)

for( i in 1:nrow(model_frame)){
  
  model         = model_frame[i,]
  Outcome_Var   = as.character(model$Y)
  Exposure_Var  = as.character(model$E)
  pathway       = as.character(model$G)
  subset.char   = as.character(model$subset)
  
  return_model_TF   =  grepl("psPGRS_all_SNPs", pathway) # whether to return the model frame: only if running the full PGRS
  
  # need stratification variables to be factors in gam.fit
  Exposure_Var_fct = paste( Exposure_Var, "fct", sep = "_" )
  UKB_pheno_PGRS[[ Exposure_Var_fct ]] = factor( ( 1*UKB_pheno_PGRS[[Exposure_Var]]+1 ) , labels= c("no", "yes") )
  
  Exposure_Var_num = paste( Exposure_Var, "num", sep = "_" )
  UKB_pheno_PGRS[[ Exposure_Var_num ]] = ( 1*UKB_pheno_PGRS[[Exposure_Var]] )
  
  if(convert.E.to.numeric_TF){
    UKB_pheno_PGRS[[ Exposure_Var ]] = UKB_pheno_PGRS[[ Exposure_Var_num ]]
  }

  
  if(is.na(subset.char)){
    subset.vec = rep(TRUE, nrow(UKB_pheno_PGRS))
  } else if(is.logical(subset.char)){
    subset.vec = subset.char
  } else if(is.character(subset.char)){
    subset.vec = eval(parse(text = subset.char), UKB_pheno_PGRS) 
    subset.vec[is.na(subset.vec)] = FALSE
  }
  
  if(testing_TF){
    subset.vec[1:length(subset.vec)] = FALSE; subset.vec[1:10000] = TRUE
    model_frame = model_frame[1:2,]
  }
  
  print(model)
  
  
  print( paste("sourcing from:", model.file))
  source( model.file )
  source( file.path("scripts", "analysis", "model_formulas.R") )
  
  
  # for glm function:
  ####

  # glm
  ####
  if(glm.TF){
    rm(glm_fit)
    
    glm_fit = glm_function(
      formula         = formula_CAD_vs_SA.PGRS, 
      family          = binomial(link = "logit"), 
      subset.vec      = subset.vec, 
      data            = UKB_pheno_PGRS, 
      Q.cvts          = c( as.character(model$E), as.character(model$G), paste(model$E, model$G, sep = ":")),
      Q.list          = list( psPGRS.ctrl = c(0,1,0), psPGRS.case = c(0,1,1)),
      q.list          = list(0,0),
      CI.alpha        = 0.05,
      save_model_data = return_model_TF   
    )
    # glm_fit$coefficients
    # glm_fit$Q.mat
    
    glm_list[[i]] = glm_fit
    names(glm_list)[i] = paste("glm", model$Y, model$E, model$G, model$subset, sep = ",")  

  }
  
  # gam
  ####
  if(gam.TF){
    rm(gam_fit)

    gam_fit = gam_function(
      formula         = formula_CAD_vs_SA.PGRS_spline, 
      family          = binomial(link = "logit"), 
      subset.vec      = subset.vec, 
      data            = UKB_pheno_PGRS, 
      Q.cvts          = c( as.character(model$E), as.character(model$G), paste(model$E, model$G, sep = ":")),
      Q.list          = list( psPGRS.ctrl = c(0,1,0), psPGRS.case = c(0,1,1)),
      q.list          = list(0,0),
      CI.alpha        = 0.05,
      save_model_data = return_model_TF,
      gam_fit_fn = "bam" #"gam" or "bam" for mgcv::bam, supposed to be big data efficient
    )
    
    gam_list[[i]] = gam_fit
    names(gam_list)[i] = paste("gam", model$Y, model$E, model$G, model$subset, sep = ",") 

  }
  
} # end loop


# name of analysis:
###
options = paste( ifelse(gam.TF, "gam,", "glm,"), "Y=", paste(Outcome_Vars, collapse = ","), ",E=",  paste(Exposure_Vars, collapse = ","), ",subset=",  paste(subsets, collapse = " "), sep = "")
print(options)

filename = paste( format( Sys.Date(), "%Y%m%d" ), format(Sys.time(), "%H%M%S"), options, pathways.slct.string, basename(model.file), "RData" , sep = ".")
filepath = file.path( UKB_data_path, filename )

# save data
###
rm(UKB_data, UKB_pheno_PGRS)
save.image( file = filepath )
