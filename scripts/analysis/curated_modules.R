select.pathways.vec = c(
  
  "psPGRS_all_SNPs",
  "psPGRS_all_annotated_genes",
  "psPGRS_OSA_pathways", # "formerly psPGRS_curated_CAD_IH_genes",
  # "psPGRS_combined_pathway",
  
  "psPGRS_hsa04066_HIF_1_signaling_pathway",
  "psPGRS_hsa04370_VEGF_signaling_pathway",
  "psPGRS_hsa04064_NF_kappa_B_signaling_pathway",
  "psPGRS_hsa04668_TNF_signaling_pathway",
  
  NULL
)

core_genes = list(
  # Core genes:
  core_genes_HIF1 = 
    c(
      "psPGRS_HIF1_only_genes_HIF1A",
      "psPGRS_HIF1_only_genes_ARNT",  # HIF1b
      NULL
    ),
  
  core_genes_VEGF = c(
    "psPGRS_HIF1_VEGF_genes_VEGFA",
    "psPGRS_HIF1_only_genes_FLT1",
    "psPGRS_VEGF_only_genes_KDR",
    NULL
  ),
  
  core_genes_TNF = c(  
    "psPGRS_NFkB_TNF_genes_TNF",
    "psPGRS_NFkB_TNF_genes_TNFRSF1A", # TNF receptor superfamily member 1A
    "psPGRS_TNF_only_genes_TNFRSF1B", # TNF receptor superfamily member 1B
    NULL
  ),
  
  core_genes_TNF_all_cells = c(  
    "psPGRS_NFkB_TNF_genes_TNF",
    "psPGRS_NFkB_TNF_genes_TNFRSF1A", # TNF receptor superfamily member 1A
    "psPGRS_TNF_only_genes_BAG4",   # now included below in NFkB upstream
    "psPGRS_NFkB_TNF_genes_TRADD",  # now included below in NFkB upstream
    "psPGRS_NFkB_TNF_genes_TRAF2",  # now included below in NFkB upstream
    "psPGRS_NFkB_TNF_genes_TRAF5",  # now included below in NFkB upstream
    "psPGRS_NFkB_TNF_genes_RIPK1",  # now included below in NFkB upstream
    NULL
  ),
  
  core_genes_TNF_WBCs = c(  
    "psPGRS_NFkB_TNF_genes_TNF",
    "psPGRS_NFkB_TNF_genes_LTA", # lymphotoxin alpha (TNF superfamily, member 1)
    "psPGRS_NFkB_only_genes_LTB", # lymphotoxin beta (TNF superfamily, member 3)
    "psPGRS_NFkB_only_genes_LTBR", # lymphotoxin beta receptor aka CD18
    "psPGRS_NFkB_only_genes_TNFSF14", #tumor necrosis factor ligand superfamily member 14
    # "TNFRSF14", # TNF receptor superfamily member 14 # Not available
    "psPGRS_TNF_only_genes_TNFRSF1B", # TNF receptor superfamily member 1B
    "psPGRS_NFkB_TNF_genes_TRAF2",  # now included below in NFkB upstream
    "psPGRS_NFkB_TNF_genes_TRAF3",  # now included below in NFkB upstream
    "psPGRS_NFkB_TNF_genes_TRAF5",  # now included below in NFkB upstream
    NULL
  ),
  
  core_genes_NFkB  = c(  
    "psPGRS_NFkB_TNF_genes_NFKBIA", # IkappaBalpha in KEGG diagram
    "psPGRS_NFkB_TNF_genes_NFKB1",  # p50
    "psPGRS_NFkB_TNF_genes_RELA",   # p65
    "psPGRS_NFkB_only_genes_NFKB2",  # p100/p52
    "psPGRS_NFkB_only_genes_RELB",   # RELB
    NULL
  ),
  
  core_genes_NFkB_canonical  = c(  
    "psPGRS_NFkB_TNF_genes_NFKBIA", # IkappaBalpha in KEGG diagram
    "psPGRS_NFkB_TNF_genes_NFKB1",  # p50
    "psPGRS_NFkB_TNF_genes_RELA",   # p65
    NULL
  ),
  
  core_genes_NFkB_nonCanonical  = c(  
    "psPGRS_NFkB_only_genes_NFKB2",  # p100/p52
    "psPGRS_NFkB_only_genes_RELB",   # RELB
    NULL
  )
)

names(core_genes) = paste("psPGRS", names(core_genes), sep = "_")
  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module_genes = list(
  
  HIF1_upstream_GF_RTK_HIF1a.translation = c(
    # GF in HIF1 diagram
    "psPGRS_HIF1_only_genes_EGF", 
    "psPGRS_HIF1_only_genes_IGF1", 
    "psPGRS_HIF1_only_genes_INS",  #insulin
    
    # RTK in HIF1 diagram
    "psPGRS_HIF1_only_genes_EGFR",
    "psPGRS_HIF1_only_genes_ERBB2", # aka HER2
    "psPGRS_HIF1_only_genes_IGF1R",
    "psPGRS_HIF1_only_genes_INSR",
    NULL
  ),
  
  HIF1_upstream_MTOR_HIF1a.translation = c(
    "psPGRS_HIF1_only_genes_MTOR", 
    "psPGRS_HIF1_only_genes_RPS6KB1",
    "psPGRS_HIF1_only_genes_RPS6KB2",
    "psPGRS_HIF1_only_genes_RPS6",
    "psPGRS_HIF1_only_genes_EIF4EBP1",
    "psPGRS_HIF1_only_genes_EIF4E", 
    "psPGRS_HIF1_only_genes_EIF4E1B", 
    "psPGRS_HIF1_only_genes_EIF4E2" 
  ),
  
  # HIF1_upstream_MAPK_HIF1a.translation = c(
  #   # identical to HIF1 VEGF shared MAPK pathway
  #   "psPGRS_HIF1_VEGF_genes_MAP2K1",
  #   "psPGRS_HIF1_VEGF_genes_MAP2K2",
  #   "psPGRS_HIF1_VEGF_genes_MAPK1",
  #   "psPGRS_HIF1_VEGF_genes_MAPK3",
  #   "psPGRS_HIF1_only_genes_MKNK1",
  #   "psPGRS_HIF1_only_genes_MKNK2",
  #   NULL
  # ),
  
  # HIF1_upstream_PI3K.AKT_HIF1a.translation = c(
  #   # identical to HIF1 VEGF shared PI3K/AKT pathway  
  #   "psPGRS_HIF1_VEGF_genes_PIK3CA",
  #   "psPGRS_HIF1_VEGF_genes_PIK3CB",
  #   "psPGRS_HIF1_VEGF_genes_PIK3CD",
  #   "psPGRS_HIF1_VEGF_genes_PIK3R1",
  #   "psPGRS_HIF1_VEGF_genes_PIK3R2",
  #   "psPGRS_HIF1_VEGF_genes_PIK3R3",
  #   "psPGRS_HIF1_VEGF_genes_AKT1",
  #   "psPGRS_HIF1_VEGF_genes_AKT2",
  #   "psPGRS_HIF1_VEGF_genes_AKT3",
  #   NULL
  # ),
 
  HIF1_upstream_IH_HIF1a.degredation = c(
    
    #NOX and IH downstream
    "psPGRS_HIF1_only_genes_CYBB",  # cytochrome b-245 beta chain (NOX in HIF1 diagram)
    
    "psPGRS_HIF1_VEGF_genes_PLCG1", # upstream of HIF1b (also in NFkB pathway)
    "psPGRS_HIF1_VEGF_genes_PLCG2", # upstream of HIF1b
    
    "psPGRS_HIF1_VEGF_genes_PRKCA", # classical protein kinase C alpha type (PKC upstream of PHD and mTOR in HIF1 diagram)
    
    #PHD
    "psPGRS_HIF1_only_genes_EGLN1", # egl-9 family hypoxia inducible factor 1 aka PHD in HIF1 diagram
    "psPGRS_HIF1_only_genes_EGLN2", # egl-9 family hypoxia inducible factor 2 aka PHD in HIF1 diagram
    "psPGRS_HIF1_only_genes_EGLN3", # egl-9 family hypoxia inducible factor 3 akNFKBIAa PHD in HIF1 diagram
    NULL
  ),
  
  HIF1_upstream_HIF1a.degredation_ubiquitin = c(
    "psPGRS_HIF1_only_genes_VHL",
    "psPGRS_HIF1_only_genes_RBX1",
    "psPGRS_HIF1_only_genes_ELOB",
    "psPGRS_HIF1_only_genes_ELOC",
    "psPGRS_HIF1_only_genes_CUL2",
    NULL
  ),
  

  HIF1_upstream_IH_HIF1b_activation = c(
    #NOX and IH downstream
    "psPGRS_HIF1_only_genes_CYBB",  # (NOX in HIF1 diagram) cytochrome b-245 beta chain 
    "psPGRS_HIF1_VEGF_genes_PLCG1", # (PLCgamma in HIF1 diagram) upstream of HIF1b (also in NFkB pathway)
    "psPGRS_HIF1_VEGF_genes_PLCG2", # (PLCgamma in HIF1 diagram) upstream of HIF1b
    
    # Hypoxia response element
    "psPGRS_HIF1_only_genes_CAMK2A", # CREB Binding Protein needed in conjunction with HIF1b in diagram (p300CBP in HIF1 diagram) 
    "psPGRS_HIF1_only_genes_CAMK2B", 
    "psPGRS_HIF1_only_genes_CAMK2D", 
    "psPGRS_HIF1_only_genes_CAMK2G",
    "psPGRS_HIF1_only_genes_CREBBP", 
    "psPGRS_HIF1_only_genes_EP300",
    "psPGRS_HIF1_only_genes_ARNT",
    NULL
  ),
  
  HIF1_downstream_RBC = c(
    "psPGRS_HIF1_only_genes_EPO", 
    "psPGRS_HIF1_only_genes_TF", 
    "psPGRS_HIF1_only_genes_TFRC",
    NULL
  ),
  
  HIF1_downstream_angiogenesis = c(
    "psPGRS_HIF1_VEGF_genes_VEGFA", # already covered in core VEGF pathway
    "psPGRS_HIF1_only_genes_FLT1",  # already covered in core VEGF pathway
    
    "psPGRS_HIF1_only_genes_EGF",
    "psPGRS_HIF1_only_genes_ANGPT1",
    "psPGRS_HIF1_only_genes_ANGPT2",
    "psPGRS_HIF1_only_genes_ANGPT4",
    
    "psPGRS_HIF1_only_genes_TIMP1",
    "psPGRS_HIF1_only_genes_SERPINE1",
    "psPGRS_HIF1_only_genes_TEK",
    NULL
    
  ),
  
  HIF1_downstream_vascular_tone = c(
    "psPGRS_HIF1_only_genes_EDN1",
    "psPGRS_HIF1_only_genes_NOS2",
    "psPGRS_HIF1_VEGF_genes_NOS3",
    "psPGRS_HIF1_only_genes_HMOX1",
    "psPGRS_HIF1_only_genes_NPPA",
    NULL
    
  ),
  
  HIF1_adjacent_genes_glucose_inhibit_citrate = c(
    "psPGRS_HIF1_only_genes_SLC2A1",
    "psPGRS_HIF1_only_genes_PDK1",
    "psPGRS_HIF1_only_genes_PDHA1",
    "psPGRS_HIF1_only_genes_PDHA2",
    "psPGRS_HIF1_only_genes_PDHB",

    NULL
    
  ),
  
  HIF1_adjacent_genes_glucose_promote_anaerobic = c(
    "psPGRS_HIF1_only_genes_HK2", 
    "psPGRS_HIF1_only_genes_HK3", 
    "psPGRS_HIF1_only_genes_HKDC1",
    "psPGRS_HIF1_only_genes_PFKL",
    # "psPGRS_HIF1_only_genes_PFKM", # not available in KEGG HIF1 pathway
    # "psPGRS_HIF1_only_genes_PFKP", # not available in KEGG HIF1 pathway
    "psPGRS_HIF1_only_genes_ALDOA",
    # "psPGRS_HIF1_only_genes_ALDOB", # not available in KEGG HIF1 pathway
    # "psPGRS_HIF1_only_genes_ALDOC", # not available in KEGG HIF1 pathway
    "psPGRS_HIF1_only_genes_GAPDH", 
    "psPGRS_HIF1_only_genes_ENO1", 
    "psPGRS_HIF1_only_genes_ENO2", 
    "psPGRS_HIF1_only_genes_ENO3", 
    "psPGRS_HIF1_only_genes_ENO4",
    "psPGRS_HIF1_only_genes_PGK1",    
    # "psPGRS_HIF1_only_genes_PGK2",  # not available in KEGG HIF1 pathway
    "psPGRS_HIF1_only_genes_PFKFB3",
    "psPGRS_HIF1_only_genes_LDHA",
    # "psPGRS_HIF1_only_genes_LDHB",    # not available in KEGG HIF1 pathway
    # "psPGRS_HIF1_only_genes_LDHC",    # not available in KEGG HIF1 pathway
    # "psPGRS_HIF1_only_genes_LDHAL6A", # not available in KEGG HIF1 pathway
    # "psPGRS_HIF1_only_genes_LDHAL6B", # not available in KEGG HIF1 pathway
    NULL
    
  ),
  
  HIF1_adjacent_genes_glucose_regulate_prolif = c(
    "psPGRS_HIF1_only_genes_BCL2", 
    "psPGRS_HIF1_only_genes_CDKN1A",
    "psPGRS_HIF1_only_genes_CDKN1B",
    NULL
  ),
  
  VEGF_downstream_proliferation = c(
    "psPGRS_HIF1_VEGF_genes_PLCG1", # upstream of HIF1b (also in NFkB and VEGF pathway)
    "psPGRS_HIF1_VEGF_genes_PLCG2", # upstream of HIF1b (also in NFkB and VEGF pathway)
    "psPGRS_HIF1_VEGF_genes_PRKCA", # classical protein kinase C alpha type (PKC upstream of PHD and mTOR in HIF1 diagram)
    "psPGRS_HIF1_VEGF_genes_PRKCB", 
    "psPGRS_HIF1_VEGF_genes_PRKCG", 
    "psPGRS_VEGF_only_genes_SPHK1",
    "psPGRS_VEGF_only_genes_SPHK2",
    "psPGRS_VEGF_only_genes_HRAS",
    "psPGRS_VEGF_only_genes_KRAS",
    "psPGRS_VEGF_only_genes_NRAS",
    "psPGRS_VEGF_only_genes_RAF1",
    "psPGRS_HIF1_VEGF_genes_MAP2K1", 
    "psPGRS_HIF1_VEGF_genes_MAP2K2",
    "psPGRS_HIF1_VEGF_genes_MAPK1", 
    "psPGRS_HIF1_VEGF_genes_MAPK3",
    NULL
  ),
  
  # VEGF_downstream_proliferation_no.MAPK = c(
  #   "psPGRS_HIF1_VEGF_genes_PLCG1", # upstream of HIF1b (also in NFkB and VEGF pathway)
  #   "psPGRS_HIF1_VEGF_genes_PLCG2", # upstream of HIF1b (also in NFkB and VEGF pathway)
  #   "psPGRS_HIF1_VEGF_genes_PRKCA", # classical protein kinase C alpha type (PKC upstream of PHD and mTOR in HIF1 diagram)
  #   "psPGRS_HIF1_VEGF_genes_PRKCB",
  #   "psPGRS_HIF1_VEGF_genes_PRKCG",
  #   "psPGRS_VEGF_only_genes_SPHK1",
  #   "psPGRS_VEGF_only_genes_SPHK2",
  #   "psPGRS_VEGF_only_genes_HRAS",
  #   "psPGRS_VEGF_only_genes_KRAS",
  #   "psPGRS_VEGF_only_genes_NRAS",
  #   "psPGRS_VEGF_only_genes_RAF1",
  #   # "psPGRS_HIF1_VEGF_genes_MAP2K1",
  #   # "psPGRS_HIF1_VEGF_genes_MAP2K2",
  #   # "psPGRS_HIF1_VEGF_genes_MAPK1",
  #   # "psPGRS_HIF1_VEGF_genes_MAPK3",
  #   NULL
  # ),
  
  VEGF_downstream_PGI2_production = c(
    "psPGRS_VEGF_only_genes_PPP3CA", 
    "psPGRS_VEGF_only_genes_PPP3CB",
    "psPGRS_VEGF_only_genes_PPP3CC",
    "psPGRS_VEGF_only_genes_PPP3R1",
    "psPGRS_VEGF_only_genes_PPP3R2",
    "psPGRS_VEGF_only_genes_NFATC2",
    "psPGRS_VEGF_only_genes_PTGS2",
    NULL
  ),
  
  VEGF_downstream_migration = c(
    "psPGRS_VEGF_only_genes_PTK2", # "protein tyrosine kinase 2" directly downstream from VEGF2 and upstream of all focal adhesion
    "psPGRS_VEGF_only_genes_PXN",  # "paxillin"                  directly downstream from VEGF2 and upstream of all focal adhesion
    
    "psPGRS_VEGF_only_genes_CDC42",   
    "psPGRS_VEGF_only_genes_MAPK14",
    "psPGRS_VEGF_only_genes_MAPK11",
    "psPGRS_VEGF_only_genes_MAPK13",
    "psPGRS_VEGF_only_genes_MAPK12",
    "psPGRS_VEGF_only_genes_MAPKAPK2",
    "psPGRS_VEGF_only_genes_MAPKAPK3", 
    NULL
  ),
  
  VEGF_downstream_permeability = c(
    
  "psPGRS_VEGF_only_genes_SRC",  
  "psPGRS_HIF1_VEGF_genes_PIK3CA",
  "psPGRS_HIF1_VEGF_genes_PIK3CB",
  "psPGRS_HIF1_VEGF_genes_PIK3CD",
  "psPGRS_HIF1_VEGF_genes_PIK3R1",
  "psPGRS_HIF1_VEGF_genes_PIK3R2",
  "psPGRS_HIF1_VEGF_genes_PIK3R3",
  "psPGRS_HIF1_VEGF_genes_AKT1",
  "psPGRS_HIF1_VEGF_genes_AKT2",
  "psPGRS_HIF1_VEGF_genes_AKT3",
  
  "psPGRS_VEGF_only_genes_RAC1",
  "psPGRS_VEGF_only_genes_RAC2",
  "psPGRS_VEGF_only_genes_RAC3",
  
  "psPGRS_HIF1_VEGF_genes_NOS3",
  NULL
  ),
  
  
  VEGF_downstream_survival = c(
    
    "psPGRS_VEGF_only_genes_SRC",  
    "psPGRS_HIF1_VEGF_genes_PIK3CA",
    "psPGRS_HIF1_VEGF_genes_PIK3CB",
    "psPGRS_HIF1_VEGF_genes_PIK3CD",
    "psPGRS_HIF1_VEGF_genes_PIK3R1",
    "psPGRS_HIF1_VEGF_genes_PIK3R2",
    "psPGRS_HIF1_VEGF_genes_PIK3R3",
    "psPGRS_HIF1_VEGF_genes_AKT1",
    "psPGRS_HIF1_VEGF_genes_AKT2",
    "psPGRS_HIF1_VEGF_genes_AKT3",
    
    "psPGRS_VEGF_only_genes_CASP9",
    "psPGRS_VEGF_only_genes_BAD",
    NULL
  ),
  
  # VEGF_downstream_migration_MAPK_genes = c(
  #   "psPGRS_VEGF_only_genes_MAPK14",
  #   "psPGRS_VEGF_only_genes_MAPK11",
  #   "psPGRS_VEGF_only_genes_MAPK13",
  #   "psPGRS_VEGF_only_genes_MAPK12",
  #   "psPGRS_VEGF_only_genes_MAPKAPK2",
  #   "psPGRS_VEGF_only_genes_MAPKAPK3",
  #   NULL
  # ),
  
  HIF1_VEGF_shared_MAPK_genes = c(
    # MAPK pathway in HIF1 and  VEGF migration
    "psPGRS_HIF1_VEGF_genes_MAP2K1",
    "psPGRS_HIF1_VEGF_genes_MAP2K2",

    "psPGRS_HIF1_VEGF_genes_MAPK1",
    "psPGRS_HIF1_VEGF_genes_MAPK3",
    "psPGRS_HIF1_only_genes_MKNK2",
    "psPGRS_HIF1_only_genes_MKNK1",

    NULL
  ),

  HIF1_VEGF_shared_PI3K.AKT_genes = c(
    # PI3K/AKT pathway
    "psPGRS_HIF1_VEGF_genes_PIK3CA",
    "psPGRS_HIF1_VEGF_genes_PIK3CB",
    "psPGRS_HIF1_VEGF_genes_PIK3CD",
    "psPGRS_HIF1_VEGF_genes_PIK3R1",
    "psPGRS_HIF1_VEGF_genes_PIK3R2",
    "psPGRS_HIF1_VEGF_genes_PIK3R3",
    "psPGRS_HIF1_VEGF_genes_AKT1",
    "psPGRS_HIF1_VEGF_genes_AKT2",
    "psPGRS_HIF1_VEGF_genes_AKT3",
    NULL
  ),
  
  HIF1_NFkB_upstream_genes = c(

    "psPGRS_HIF1_only_genes_IL6",    # upstream of STAT3
    "psPGRS_HIF1_only_genes_IL6R",
    "psPGRS_HIF1_only_genes_STAT3",  # co-regulates HIF1 activation along with NFkB in HIF1 diagram

    "psPGRS_HIF1_only_genes_IFNG",   # upstream of NFkB
    "psPGRS_HIF1_only_genes_IFNGR1",
    "psPGRS_HIF1_only_genes_IFNGR2",

    "psPGRS_HIF1_only_genes_TLR4",   # upstream of NFkB (also part of NFkB pathway)

    NULL
  ),
  
  NFkB_atypical_pathway = c(
    
    # Atypical pathways involve hypoxia as well as:
    
    # EGF:
    "psPGRS_HIF1_only_genes_EGF",   # upstream of MAPK, PI3K, and AKT, and thence HIF1A
    "psPGRS_HIF1_only_genes_EGFR",
    
    "psPGRS_NFkB_TNF_genes_NFKBIA", # NFKB inhibitor alpha
    "psPGRS_NFkB_TNF_genes_NFKB1", # nuclear factor kappa B subunit 1
    "psPGRS_NFkB_TNF_genes_RELA", # RELA proto-oncogene, NF-kB subunit
    
    "psPGRS_NFkB_only_genes_PLAU",
    "psPGRS_NFkB_only_genes_BCL2L1",
    
    # HER2:
    "psPGRS_HIF1_only_genes_ERBB2", # aka HER2
    
    "psPGRS_NFkB_only_genes_CSNK2A1", 
    "psPGRS_NFkB_only_genes_CSNK2A2",
    "psPGRS_NFkB_only_genes_CSNK2B",
    "psPGRS_NFkB_only_genes_CSNK2A3",
    "psPGRS_NFkB_only_genes_BCL2A1",

    # CNTF as listed in NFkB atypical path is not available. Corresponds to cytokine-cytokine interaction pathway IL6/12-like:
      # CLCF1
      # CNTF
      # CNTFR
      # LIFR
      # IL6ST

    # NGF as listed in NFkB atypical path is not available. Corresponds to cytokine-cytokine interaction pathway TNF-family:
      # NGF,
      # NGFR

    NULL

  ),

  NFkB_atypical_upstream = c(
    # IL1beta
    "psPGRS_NFkB_TNF_genes_IL1B",
    "psPGRS_NFkB_only_genes_IL1R1",
    "psPGRS_NFkB_only_genes_MYD88",
    "psPGRS_NFkB_only_genes_IRAK1",
    "psPGRS_NFkB_only_genes_IRAK4",
    "psPGRS_NFkB_only_genes_TRAF6",

    # TNFalpha (dup TNF core genes in all cells as listed above)
    "psPGRS_NFkB_TNF_genes_TNF",
    "psPGRS_NFkB_TNF_genes_TNFRSF1A", # TNF receptor superfamily member 1A
    "psPGRS_TNF_only_genes_BAG4",
    "psPGRS_NFkB_TNF_genes_TRADD",
    "psPGRS_NFkB_TNF_genes_TRAF2",
    "psPGRS_NFkB_TNF_genes_TRAF5",
    "psPGRS_NFkB_TNF_genes_RIPK1",

    # RANK
    "psPGRS_NFkB_only_genes_TNFSF11",   # aka RANKL
    "psPGRS_NFkB_only_genes_TNFRSF11A", # aka RANK
    # "psPGRS_NFkB_only_genes_TRAF6",   # dup above
    # "psPGRS_NFkB_only_genes_TRAF2",   # dup above

    NULL
  ),

  TNF_pathway_JNK = c(
    
    # already in core genes:
    # TNF, # tumor necrosis factor
    # TNFRSF1A, # TNF receptor superfamily member 1A
    # RIPK1, # receptor interacting serine/threonine kinase 1
    # TRADD, # TNFRSF1A associated via death domain
    # TRAF2, # TNF receptor associated factor 2
    # TRAF5, # TNF receptor associated factor 5
    
    "psPGRS_NFkB_TNF_genes_MAP3K7", # mitogen-activated protein kinase kinase kinase 7
    "psPGRS_NFkB_TNF_genes_TAB1", # TGF-beta activated kinase 1 (MAP3K7) binding protein 1
    "psPGRS_NFkB_TNF_genes_TAB2", # TGF-beta activated kinase 1 (MAP3K7) binding protein 2
    
    "psPGRS_TNF_only_genes_MAP2K4", # mitogen-activated protein kinase kinase 4
    "psPGRS_TNF_only_genes_MAP2K7", # mitogen-activated protein kinase kinase 7
    "psPGRS_TNF_only_genes_MAPK8", # mitogen-activated protein kinase 8
    "psPGRS_TNF_only_genes_MAPK9", # mitogen-activated protein kinase 9
    "psPGRS_TNF_only_genes_MAPK10", # mitogen-activated protein kinase 10
    "psPGRS_TNF_only_genes_FOS", # Fos proto-oncogene, AP-1 transcription factor subunit
    "psPGRS_TNF_only_genes_JUN", # Jun proto-oncogene, AP-1 transcription factor subunit
     NULL
  ),
  
  TNF_pathway_p38 = c(
    
    # already in core genes:
    # TNF, # tumor necrosis factor
    # TNFRSF1A, # TNF receptor superfamily member 1A
    # RIPK1, # receptor interacting serine/threonine kinase 1
    # TRADD, # TNFRSF1A associated via death domain
    # TRAF2, # TNF receptor associated factor 2
    # TRAF5, # TNF receptor associated factor 5
    
    "psPGRS_NFkB_TNF_genes_MAP3K7", # mitogen-activated protein kinase kinase kinase 7
    "psPGRS_NFkB_TNF_genes_TAB1", # TGF-beta activated kinase 1 (MAP3K7) binding protein 1
    "psPGRS_NFkB_TNF_genes_TAB2", # TGF-beta activated kinase 1 (MAP3K7) binding protein 2
    
    "psPGRS_TNF_only_genes_MAP2K3", # mitogen-activated protein kinase kinase 3
    "psPGRS_TNF_only_genes_MAP2K6", # mitogen-activated protein kinase kinase 6
    "psPGRS_TNF_only_genes_MAPK14", # mitogen-activated protein kinase 14
    "psPGRS_TNF_only_genes_RPS6KA4", # ribosomal protein S6 kinase A4
    "psPGRS_TNF_only_genes_RPS6KA5", # ribosomal protein S6 kinase A5
    "psPGRS_TNF_only_genes_CREB3",  # cAMP responsive element binding protein 3
    "psPGRS_TNF_only_genes_CREB1",  # cAMP responsive element binding protein 1
    "psPGRS_TNF_only_genes_ATF2",   # activating transcription factor 2
    "psPGRS_TNF_only_genes_ATF6B",  # activating transcription factor 6 beta
    "psPGRS_TNF_only_genes_CREB3L4", # cAMP responsive element binding protein 3 like 4
    "psPGRS_TNF_only_genes_ATF4",   # activating transcription factor 4
    "psPGRS_TNF_only_genes_CREB3L2", # cAMP responsive element binding protein 3 like 2
    "psPGRS_TNF_only_genes_CREB3L3", # cAMP responsive element binding protein 3 like 3
    "psPGRS_TNF_only_genes_CREB3L1", # cAMP responsive element binding protein 3 like 1
    "psPGRS_TNF_only_genes_CREB5",   # cAMP responsive element binding protein 5
    NULL
  ),
  
  NFKB_TNF_pathway= c(
    
    # already in core genes:
    "psPGRS_NFkB_TNF_genes_TNF", # tumor necrosis factor
    "psPGRS_NFkB_TNF_genes_TNFRSF1A", # TNF receptor superfamily member 1A
    "psPGRS_NFkB_TNF_genes_RIPK1", # receptor interacting serine/threonine kinase 1
    "psPGRS_NFkB_TNF_genes_TRADD", # TNFRSF1A associated via death domain
    "psPGRS_NFkB_TNF_genes_TRAF2", # TNF receptor associated factor 2
    "psPGRS_NFkB_TNF_genes_TRAF5", # TNF receptor associated factor 5
    
    "psPGRS_NFkB_TNF_genes_BIRC2",
    "psPGRS_NFkB_TNF_genes_BIRC3",
    
    "psPGRS_NFkB_TNF_genes_MAP3K14",  # mitogen-activated protein kinase kinase kinase 14
    "psPGRS_NFkB_TNF_genes_MAP3K7", # mitogen-activated protein kinase kinase kinase 7
    "psPGRS_NFkB_TNF_genes_TAB1", # TGF-beta activated kinase 1 (MAP3K7) binding protein 1
    "psPGRS_NFkB_TNF_genes_TAB2", # TGF-beta activated kinase 1 (MAP3K7) binding protein 2
    "psPGRS_NFkB_TNF_genes_TAB3", # TGF-beta activated kinase 1 (MAP3K7) binding protein 3
    
    "psPGRS_NFkB_TNF_genes_CHUK", # conserved helix-loop-helix ubiquitous kinase
    "psPGRS_NFkB_TNF_genes_IKBKB", # inhibitor of nuclear factor kappa B kinase subunit beta
    "psPGRS_NFkB_TNF_genes_IKBKG", # inhibitor of nuclear factor kappa B kinase subunit gamma
    "psPGRS_NFkB_TNF_genes_NFKBIA", # NFKB inhibitor alpha
    "psPGRS_NFkB_TNF_genes_NFKB1", # nuclear factor kappa B subunit 1
    "psPGRS_NFkB_TNF_genes_RELA", # RELA proto-oncogene, NF-kB subunit
    NULL
  ),
  
  NFKB_pathway_downstream_survival = c(
    
    "psPGRS_NFkB_TNF_genes_CFLAR", # CASP8 and FADD like apoptosis regulator [KO:K04724]
    "psPGRS_NFkB_only_genes_XIAP", # X-linked inhibitor of apoptosis [KO:K04725] [EC:2.3.2.27]
    "psPGRS_NFkB_only_genes_BCL2L1", # BCL2 like 1 [KO:K04570]
    "psPGRS_NFkB_only_genes_BCL2", # BCL2 apoptosis regulator [KO:K02161]
    # "GADD45A", # growth arrest and DNA damage inducible alpha [KO:K04402]
    "psPGRS_NFkB_only_genes_GADD45B", # growth arrest and DNA damage inducible beta [KO:K04402]
    # "GADD45G", # growth arrest and DNA damage inducible gamma [KO:K04402]
    "psPGRS_NFkB_TNF_genes_TRAF1", # TNF receptor associated factor 1 [KO:K03172]
    "psPGRS_NFkB_only_genes_BCL2A1", # BCL2 related protein A1 [KO:K02162]
    NULL
    
  ),

  NFKB_pathway_downstream_inflammation = c(
    "psPGRS_NFkB_only_genes_CXCL8", # C-X-C motif chemokine ligand 8 (IL8) [KO:K10030]
    "psPGRS_NFkB_TNF_genes_TNFAIP3", # TNF alpha induced protein 3 [KO:K11859] [EC:3.4.19.12]
    "psPGRS_NFkB_TNF_genes_PTGS2", # prostaglandin-endoperoxide synthase 2 [KO:K11987] [EC:1.14.99.1]
    "psPGRS_NFkB_only_genes_CCL4", # C-C motif chemokine ligand 4 [KO:K12964]
    "psPGRS_NFkB_only_genes_CCL4L2", # C-C motif chemokine ligand 4 like 2 [KO:K12964]
    "psPGRS_NFkB_only_genes_CCL4L1", # C-C motif chemokine ligand 4 like 1 [KO:K12964]
    "psPGRS_TNF_only_genes_CXCL1", # C-X-C motif chemokine ligand 1 [KO:K05505]
    "psPGRS_NFkB_TNF_genes_CXCL2", # C-X-C motif chemokine ligand 2 [KO:K05505]
    "psPGRS_TNF_only_genes_CXCL3", # C-X-C motif chemokine ligand 3 [KO:K05505]
    "psPGRS_NFkB_TNF_genes_VCAM1" # vascular cell adhesion molecule 1 [KO:K06527]
  ),

  TNF_pathway_apoptosis = c(
    "psPGRS_TNF_only_genes_FADD",
    "psPGRS_TNF_only_genes_ITCH",
    "psPGRS_NFkB_TNF_genes_CFLAR",
    "psPGRS_TNF_only_genes_CASP8",
    "psPGRS_TNF_only_genes_CASP10",
    "psPGRS_TNF_only_genes_CASP3",
    "psPGRS_TNF_only_genes_CASP7",
    NULL
  ),
  
  TNF_pathway_necroptosis = c(
    "psPGRS_NFkB_TNF_genes_RIPK1",
    "psPGRS_TNF_only_genes_RIPK3",
    "psPGRS_TNF_only_genes_MLKL",
    "psPGRS_TNF_only_genes_PGAM5",
    "psPGRS_TNF_only_genes_DNM1L",
    NULL
  ),
  
  TNF_downstream_immune_signaling = c(
    "psPGRS_TNF_only_genes_CCL2", # C-C motif chemokine ligand 2 [KO:K14624]
    "psPGRS_TNF_only_genes_CCL5", # C-C motif chemokine ligand 5 [KO:K12499]
    "psPGRS_TNF_only_genes_CCL20", # C-C motif chemokine ligand 20 [KO:K14625]
    "psPGRS_TNF_only_genes_CXCL1", # C-X-C motif chemokine ligand 1 [KO:K05505]
    "psPGRS_NFkB_TNF_genes_CXCL2", # C-X-C motif chemokine ligand 2 [KO:K05505]
    "psPGRS_TNF_only_genes_CXCL3", # C-X-C motif chemokine ligand 3 [KO:K05505]
    "psPGRS_TNF_only_genes_CXCL5", # C-X-C motif chemokine ligand 5 [KO:K05506]
    # "CXCL6", # C-X-C motif chemokine ligand 6 [KO:K05506]
    "psPGRS_TNF_only_genes_CXCL10", # C-X-C motif chemokine ligand 10 [KO:K12671]
    "psPGRS_TNF_only_genes_CX3CL1", # C-X3-C motif chemokine ligand 1 [KO:K05508]
    "psPGRS_TNF_only_genes_CSF1", # colony stimulating factor 1 [KO:K05453]
    "psPGRS_TNF_only_genes_CSF2", # colony stimulating factor 2 [KO:K05427]
    "psPGRS_TNF_only_genes_FAS", # Fas cell surface death receptor [KO:K04390]
    "psPGRS_TNF_only_genes_IL18R1", # interleukin 18 receptor 1 [KO:K05173]
    "psPGRS_TNF_only_genes_JAG1", # jagged canonical Notch ligand 1 [KO:K06052]
    "psPGRS_NFkB_TNF_genes_IL1B", # interleukin 1 beta [KO:K04519]
    "psPGRS_TNF_only_genes_IL6", # interleukin 6 [KO:K05405]
    "psPGRS_TNF_only_genes_IL15", # interleukin 15 [KO:K05433]
    "psPGRS_TNF_only_genes_LIF", # LIF interleukin 6 family cytokine [KO:K05419]
    "psPGRS_NFkB_TNF_genes_LTA", # lymphotoxin alpha [KO:K05468]
    "psPGRS_TNF_only_genes_BCL3", # BCL3 transcription coactivator [KO:K09258]
    "psPGRS_TNF_only_genes_SOCS3", # suppressor of cytokine signaling 3 [KO:K04696]
    "psPGRS_NFkB_TNF_genes_TNFAIP3", # TNF alpha induced protein 3 [KO:K11859] [EC:3.4.19.12]
    "psPGRS_VEGF_only_genes_RAF1", # TNF receptor associated factor 1 [KO:K03172]
    NULL
  ),
  
  TNF_downstream_intracellular_signaling = c(
    "psPGRS_TNF_only_genes_FOS",
    "psPGRS_TNF_only_genes_JUN",
    "psPGRS_TNF_only_genes_JUNB", # JunB proto-oncogene, AP-1 transcription factor subunit [KO:K09028]
    	"psPGRS_TNF_only_genes_MMP3", # matrix metallopeptidase 3 [KO:K01394] [EC:3.4.24.17]
    	"psPGRS_TNF_only_genes_MMP9", # matrix metallopeptidase 9 [KO:K01403] [EC:3.4.24.35]
    	"psPGRS_TNF_only_genes_MMP14", # matrix metallopeptidase 14 [KO:K07763] [EC:3.4.24.80]
    	"psPGRS_TNF_only_genes_EDN1", # endothelin 1 [KO:K16366]
    	"psPGRS_TNF_only_genes_VEGFC", # vascular endothelial growth factor C [KO:K05449]
    	# "VEGFD", # vascular endothelial growth factor D [KO:K05449]
      "psPGRS_TNF_only_genes_NOD2", # nucleotide binding oligomerization domain containing 2 [KO:K10165]
    	"psPGRS_NFkB_TNF_genes_ICAM1", # intercellular adhesion molecule 1 [KO:K06490]
    	"psPGRS_TNF_only_genes_SELE", # selectin E [KO:K06494]
    	"psPGRS_NFkB_TNF_genes_VCAM1", # vascular cell adhesion molecule 1 [KO:K06527]
    	"psPGRS_NFkB_TNF_genes_PTGS2", # prostaglandin-endoperoxide synthase 2 [KO:K11987] [EC:1.14.99.1]
    	NULL
  )
) #end of list

names(module_genes) = paste("psPGRS", names(module_genes), sep = "_")
length(module_genes)
module_gene_names = lapply( module_genes, function(x) stringr::str_split( x, "_", simplify = TRUE)[,5])
module_gene_names_out = lapply( module_gene_names, function(x) paste(sort(x), collapse = ", "))
print(module_gene_names_out)

