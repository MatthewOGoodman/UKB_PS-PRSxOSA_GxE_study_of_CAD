#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# LOAD DATA ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Load scripts ####

# Run from the repository root:  Rscript scripts/analysis/prepare_analysis_data.R
# Required environment: PSPRS_ANALYSIS_DATA_DIR (UK Biobank phenotypes merged with the PS-PRS;
# available only under an approved UK Biobank application).
suppressPackageStartupMessages({ library(dplyr); library(stringr) })
source("config/paths.R")
psprs_require("UKB_data_path")
source(file.path("scripts", "analysis", "curated_modules.R"))

UKB_pheno_PGRS = readRDS(file = file.path(UKB_data_path, "UKB_pheno_PGRS_38.psPGRS_HIF1.VEGF.NFkB.TNF.genes_20210226") )

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Final data cleaning and variable definitions  ####
####


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Final Phenotype Definitions ####

# Covariates ####
####

# check composite covariates
####
colSums(  is.na( UKB_pheno_PGRS[,grepl("cmp_", colnames(UKB_pheno_PGRS))] ) )


# list the variable types, based on prefix
variable_classes = colnames(UKB_pheno_PGRS)[grep("_", colnames(UKB_pheno_PGRS))]
variable_classes = unlist(unique( lapply(stringr::str_split(variable_classes[-1], pattern = "_"), "[", 1) ))
variable_classes
# list the individual, unclassified variables
variable_singletons = colnames(UKB_pheno_PGRS)[!grepl("_", colnames(UKB_pheno_PGRS)) | grepl("ID_", colnames(UKB_pheno_PGRS))]
variable_singletons
# list the psPGRS
t(t(colnames(UKB_pheno_PGRS)[grepl("PGRS", colnames(UKB_pheno_PGRS), ignore.case = TRUE)]))
# list the phecodes
t(t(colnames(UKB_pheno_PGRS)[grepl("cmp.phecode", colnames(UKB_pheno_PGRS), ignore.case = TRUE)]))
# Note that cmp.phecode_'s are counts


# check some individual covariates
levels(as.factor(UKB_pheno_PGRS$HD_smoke_status_derived))
t(t(colnames(UKB_pheno_PGRS)[grepl("CAD", colnames(UKB_pheno_PGRS), ignore.case = TRUE)]))

# check for variables
####
colnames(UKB_pheno_PGRS) [grepl("cmp_sleep", colnames(UKB_pheno_PGRS))]
colnames(UKB_pheno_PGRS) [grepl("hyper", colnames(UKB_pheno_PGRS))]
colnames(UKB_pheno_PGRS) [grepl("snore", colnames(UKB_pheno_PGRS))]
colnames(UKB_pheno_PGRS) [grepl("cmp_CSA", colnames(UKB_pheno_PGRS))]



#Shorten names
####

UKB_pheno_PGRS =   dplyr::mutate(
  .data = UKB_pheno_PGRS,
  
  age.ge50   = Age0 - 50,
  sexM       = geno.PCs_Inferred_Gender == "M",
  sexM.age   = age.ge50 * sexM,
  sexM.BMI   = BMI0 * sexM,
  white      = geno.PCs_in_white_British_ancestry_subset,
  likelyT1dm = HD_probable_t1dm,
  likelyT2dm = HD_probable_t2dm,
  hypertn    = HD_hypertension,
  COPD       = cmp_COPD,
  FEV.FVC.z  = lung.fn_FEVFVC.ratio.zscore_20158,
  smoke      = as.factor(HD_smoke_status_derived)
)

# (cmp_stroke | cmp_heart_failure)
# 
# colnames(UKB_pheno_PGRS)[grepl("strok", colnames(UKB_pheno_PGRS))]
# with(UKB_pheno_PGRS, sum( AA_stroke_isch & !cmp_stroke_isch))
# with(UKB_pheno_PGRS, sum( cmp_stroke_isch & !AA_stroke_isch))
# 




# recompute sexM as sex_fct
####
UKB_pheno_PGRS$sex_fct      = factor( ( UKB_pheno_PGRS$sexM+1  ) , labels= c("f", "m") )

summary( with(UKB_pheno_PGRS, interaction(smoke,sex_fct)) )

UKB_pheno_PGRS$cmp_hypertension = with(UKB_pheno_PGRS, cmp.phecode_high.blood.pressure > 0 | card.sr_hypertension)
# with(UKB_pheno_PGRS, cor( cmp_hypertension , hypertn, use = "pairwise.complete"))
# with(UKB_pheno_PGRS, cor( card.sr_hypertension , hypertn, use = "pairwise.complete"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Final Sleep Exposure Definitions ####
####

# Earlier definitions from the phenotype pipeline:
###%%%%%%%%%%%%%%###
# HESIN exclusions
####
# (defined previously, based on evidence of CSA comorbidities prior to sleep apnea ICD10)
# cmp_CSA_exclusions_hesin =  (apnea_w.prior_ICD10_stroke | apnea_w.prior_ICD10_hrtfail | apnea_w.ICD10_brainstem ), # susan's definition

# Self-report exclusions: based on prevalent stroke, heart_failure or spinal/motor disease at time of sleep apnea self report
# cmp_CSA.exclusions_sr    =  (crbr.sr_stroke | card.sr_heart_failure | neur.sr.cmb_spinal_motor ), # susan's definition

# OSA CASES:####
###
# cmp_OSA       = (cmp_sleep.apnea.w.snoring_sr | cmp_sleep.apnea.w.snoring_hesin) & ! (cmp_CSA.exclusions_sr | cmp_CSA_exclusions_hesin) ,
# cmp_OSA_sr    = cmp_sleep.apnea.w.snoring_sr    & ! cmp_CSA.exclusions_sr,
# cmp_OSA_hesin = cmp_sleep.apnea.w.snoring_hesin & ! cmp_CSA_exclusions_hesin,
# cmp_OSA_noHF.noStr  = (cmp_sleep.apnea.w.snoring_sr | cmp_sleep.apnea.w.snoring_hesin) & ! (cmp_CSA.exclusions_sr | cmp_CSA_exclusions_hesin)  & ! (cmp_stroke | cmp_heart_failure),

# cmp_OSA_fct = as.factor(cmp_OSA),
# cmp_CAD_and_OSA_cases = cmp_CAD & cmp_OSA,

# Note: these are the definitions behind cmp_sleep.apnea.w.snoring
# cmp_sleep.apnea.w.snoring_hesin = (hesin_apnea_ICD10_TF & sleep.sr_snoring) , # ICD10 cases plus exclusions
# cmp_sleep.apnea.w.snoring_sr    = (brth.sr_sleep_apnea & sleep.sr_snoring),   # self-report cases plus exclusions
# cmp_sleep.apnea.w.snoring       = cmp_sleep.apnea.w.snoring_hesin | cmp_sleep.apnea.w.snoring_sr, # cases plus exclusions
###%%%%%%%%%%%%%%###

# Obstructive Sleep Apnea
###
# recall:   cmp_OSA is defined as (cmp_sleep.apnea.w.snoring_sr | cmp_sleep.apnea.w.snoring_hesin) & ! (cmp_CSA.exclusions_sr | cmp_CSA_exclusions_hesin) 
# Earlier definitions
# UKB_pheno_PGRS$cmp_OSA.controls  = cmp_OSA.controls  = with(UKB_pheno_PGRS, !cmp_OSA & !cmp_snore.daysleep & ! (cmp_stroke | cmp_heart_failure))
# UKB_pheno_PGRS$cmp_OSA.cases     = cmp_OSA.cases     = with(UKB_pheno_PGRS,  cmp_OSA & ! (cmp_stroke | cmp_heart_failure))
# 
# UKB_pheno_PGRS$cmp_snr.day.controls  = cmp_snr.day.controls  = with(UKB_pheno_PGRS, !cmp_snore.daysleep & ! (cmp_stroke | cmp_heart_failure))
# UKB_pheno_PGRS$cmp_snr.day.cases     = cmp_snr.day.cases     = with(UKB_pheno_PGRS,  cmp_snore.daysleep & ! (cmp_stroke | cmp_heart_failure))
# 
# UKB_pheno_PGRS$snore.SA         = with(UKB_pheno_PGRS, cmp_sleep.apnea.w.snoring) 

# Definitions used for the paper
UKB_pheno_PGRS =   dplyr::mutate(
  .data = UKB_pheno_PGRS,
  cmp_SA = (hesin_apnea_ICD10_TF | brth.sr_sleep_apnea),
  cmp_SA.snoring   =  cmp_sleep.apnea.w.snoring_sr | cmp_sleep.apnea.w.snoring_hesin,
  cmp_SA.snoring_controls  = !cmp_SA.snoring & !cmp_snore.daysleep, # moderate or severe daytime sleepiness from the phenotype pipeline
  
  sleep.sr_insomnia_often = sleep.sr_insomnia == 3,
  cmp_SA_AND_snore.AND.daysleep  = (hesin_apnea_ICD10_TF | brth.sr_sleep_apnea) & (sleep.sr_snoring & cmp_daysleep_mod.sev),
  cmp_SA_AND_snore.OR.daysleep   = (hesin_apnea_ICD10_TF | brth.sr_sleep_apnea) & (sleep.sr_snoring | cmp_daysleep_mod.sev),
  cmp_SA_AND_S.AND.D_controls  = !cmp_SA_AND_snore.AND.daysleep & !cmp_snore.daysleep,
  cmp_SA_AND_S.OR.D_controls  = !cmp_SA_AND_snore.OR.daysleep & !cmp_snore.daysleep
  )

# with(UKB_pheno_PGRS, all.equal(cmp_SA.snoring, cmp_sleep.apnea.w.snoring))
# with(UKB_pheno_PGRS, cor(hesin_apnea_ICD10_TF, brth.sr_sleep_apnea))
# with(UKB_pheno_PGRS, cor(cmp_SA, cmp_snore.daysleep))
# 
# colnames(UKB_pheno_PGRS) [grepl("apnea", colnames(UKB_pheno_PGRS))]
# with(UKB_pheno_PGRS, sum(cmp_SA))
# with(UKB_pheno_PGRS, sum(cmp_snore.daysleep))
# with(UKB_pheno_PGRS, sum(cmp_SA_AND_snore.OR.daysleep))
# with(UKB_pheno_PGRS, sum(cmp_SA.snoring))
# with(UKB_pheno_PGRS, sum(cmp_SA_AND_snore.AND.daysleep))
# 
# with(UKB_pheno_PGRS, cor(cmp_SA_AND_snore.AND.daysleep, cmp_SA_AND_snore.OR.daysleep))
# 
# with(UKB_pheno_PGRS, cor(cmp_snore.daysleep, hesin_apnea_ICD10_TF))
# with(UKB_pheno_PGRS, cor(cmp_snore.daysleep, brth.sr_sleep_apnea))
# with(UKB_pheno_PGRS, cor(cmp_snore.daysleep, hesin_apnea_ICD10_TF | brth.sr_sleep_apnea))
# 
# with(UKB_pheno_PGRS, cor(HD_icd10_apnea, hesin_apnea_ICD10_TF, use = "pairwise.complete"))
# with(UKB_pheno_PGRS, cor(HD_apnea_self, brth.sr_sleep_apnea, use = "pairwise.complete"))
# 
# sum(UKB_pheno_PGRS$HD_icd10_apnea, na.rm = TRUE)
# sum(UKB_pheno_PGRS$hesin_apnea_ICD10_TF, na.rm = TRUE)


# Pathways: ####
####

# calculate new curated module pathway psPGRSs using single-gene PS-PGRS
####

# combine pathways

curated_modules = c(core_genes, module_genes )
length(curated_modules)
unlist(curated_modules)[!unlist(curated_modules) %in% names(UKB_pheno_PGRS)]

psPGRS_columns = colnames(UKB_pheno_PGRS)[grepl("psPGRS_", colnames(UKB_pheno_PGRS))]
unavailable_genes = psPGRS_columns[ which(colSums(  is.na( UKB_pheno_PGRS[,psPGRS_columns] ) ) > 100000) ]
curated_modules_reduce = lapply(curated_modules, function(x) x[! x %in% unavailable_genes])


# calculate new psPGRS
new_curated_pathways = lapply(curated_modules_reduce, function(x) rowSums(UKB_pheno_PGRS[,x]) )
names(new_curated_pathways) = names(curated_modules_reduce) 
UKB_pheno_PGRS = data.frame( UKB_pheno_PGRS, as.data.frame(new_curated_pathways))
UKB_pheno_PGRS = UKB_pheno_PGRS[,!(grepl("^psPGRS", colnames(UKB_pheno_PGRS)) & grepl("\\.[[:digit:]]$", colnames(UKB_pheno_PGRS)) )]
# colnames(UKB_pheno_PGRS)[grepl("^psPGRS", colnames(UKB_pheno_PGRS)) & grepl("\\.[[:digit:]]$", colnames(UKB_pheno_PGRS)) ]
rm(new_curated_pathways)


# Scale the PGRSs: ####

scale_psPGRS_TF = TRUE
pathways.all = colnames(UKB_pheno_PGRS) [grepl("psPGRS", colnames(UKB_pheno_PGRS))]
if(scale_psPGRS_TF){ UKB_pheno_PGRS[,pathways.all] = data.frame(scale( UKB_pheno_PGRS[,pathways.all] )) }
# check PGRS standardization

# check and flip the PGRS's to be risk-increasing
####
PGRS.model = summary(lm(CAD ~ psPGRS_all_SNPs, data = UKB_pheno_PGRS))$coefficients
PGRS.reversed = PGRS.model["psPGRS_all_SNPs", "Estimate"] < 0
PGRS.reversed

if(PGRS.reversed){ UKB_pheno_PGRS[,pathways.all] = -1* UKB_pheno_PGRS[,pathways.all]}
PGRS.model = summary(lm(CAD ~ psPGRS_all_SNPs, data = UKB_pheno_PGRS))$coefficients
PGRS.reversed = PGRS.model["psPGRS_all_SNPs", "Estimate"] < 0
PGRS.reversed


# Save Results
# UKB_pheno_PGRS = saveRDS(UKB_pheno_PGRS, file = file.path(UKB_data_path, "UKB_pheno_PGRS_38.psPGRS_HIF1.VEGF.NFkB.TNF.genes_analytic_20210226") )
# UKB_pheno_PGRS = saveRDS(UKB_pheno_PGRS, file = file.path(UKB_data_path, "UKB_pheno_PGRS_38.psPGRS_HIF1.VEGF.NFkB.TNF.genes_analytic_20210315") )
UKB_pheno_PGRS = saveRDS(UKB_pheno_PGRS, file = file.path(UKB_data_path, "UKB_pheno_PGRS_38.psPGRS_HIF1.VEGF.NFkB.TNF.genes_analytic_20210407") ) # correct order with scaling after create new pathways


