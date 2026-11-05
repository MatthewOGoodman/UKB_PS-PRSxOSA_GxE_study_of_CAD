model_vars_long = c(
  
  # basic
  "sex_fct",
  "BMI0",
  "age.ge50",
  
  "sex_fct:age.ge50",
  "sex_fct:BMI0",
  
  "smoke",
  "smoke:sex_fct",
  
  # genetic
  "white",
  "geno.PCs_PC1", "geno.PCs_PC2", "geno.PCs_PC3", "geno.PCs_PC4", "geno.PCs_PC5",
  "BiLEVE",
  
  # CSA comorbidities 
  # "cmp_heart_failure",
  # "I(cmp_OSA*cmp_heart_failure)", # in case the exclusions cause heart_failure definition to be different between cases and control
  # "cmp_stroke",
  # "I(cmp_OSA*cmp_stroke)", # in case the exclusions cause stroke definition to be different between cases and control
  
  # CVD comorbidities
  # "likelyT1dm",
  # "likelyT1dm:sex_fct",
  # "likelyT2dm",
  # "likelyT2dm:sex_fct",
  # "cmp_hypertension",
  # "cmp_hypertension:sex_fct",
  
  # pulmonary comorbidities
  "cmp_COPD",
  "cmp_asthma",
  # "cmp_pulmn.fib",
  # "cmp_resp.failure",
  # "cmp_other.lung.collapse",
 
  NULL # so that you can turn variables on and off at will by commenting
)

model_vars_spline = c(
  
  # basic
  
  "sex_fct",
  # "sexM",
  # 's( age.ge50, k = 5)',
  # 's( BMI0, k = 5)',
  # 's( sexM.age, k = 5)',
  # 's( sexM.BMI, k = 5)',
  
  # tensor product basis for BMI and age, stratified by sex
  # 'te(age.ge50,BMI0,bs=c("tp","tp"),k=c(4,4),m=c(3,3),by=sex_fct)', #to 20200129
  'te(age.ge50,BMI0,bs=c("tp","tp"),k=c(7,7),m=c(2,2),by=sex_fct)', 
  # k is degrees of freedom, m is order of the penalty terms -- see s{mgcv} help file "e.g. 2 for normal cubic spline penalty with 2nd derivatives when using default t.p.r.s basis)"
  
  "smoke",
  "smoke:sex_fct",
  
  # genetic
  "white",
  "geno.PCs_PC1", "geno.PCs_PC2", "geno.PCs_PC3", "geno.PCs_PC4", "geno.PCs_PC5",
  "BiLEVE",
  
  # CSA comorbidities 
  # "cmp_heart_failure",
  # "I(cmp_OSA*cmp_heart_failure)", # in case the exclusions cause heart_failure definition to be different between cases and control
  # "cmp_stroke",
  # "I(cmp_OSA*cmp_stroke)", # in case the exclusions cause stroke definition to be different between cases and control
  
  # CVD comorbidities
  # "likelyT1dm",
  # "likelyT1dm:sex_fct",
  # "likelyT2dm",
  # "likelyT2dm:sex_fct",
  # "cmp_hypertension",
  # "cmp_hypertension:sex_fct",

  "cmp_COPD",
  "cmp_asthma",
  # "cmp_pulmn.fib",
  # "cmp_resp.failure", 
  # "cmp_other.lung.collapse",
  
  NULL
)

model_vars_long
model_vars_spline
