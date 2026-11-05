# loadPGRS ----

load_split_save_PGRS = function(file_path, 
                                out_name = "PGRS_CAD_Khera_LDPred_rho0.001", 
                                out_folder,
                                
                                # colnames_key: char vec is existing names in PGRS; names of char vec should not be user modified, needed below
                                colnames_key =   c(
                                  variant = "variant",
                                  chr = "chr",
                                  pos = "position_hg19",
                                  effect_weight = "effect_weight",
                                  effect_allele = "effect_allele",
                                  A1 = "A1",
                                  A2 = "A2"),
                                chr_prefix = "chr", # prefix, to convert e.g. "22" to "chr22", otherwise, put NULL for use original
                                
                                split_variant_str_TF = FALSE,
                                char_split_names = c("chr_str", "pos_str", "A1_str", "A2_str"),
                                char_sep_regexp = "\\:",
                                
                                skip_lines = 15, n_max_lines = Inf,
                                split_save_TF = TRUE, return_type = "NULL", 
                                verbose_TF = TRUE
                                # choose return_type %in% c("NULL", "tibble", "data.table", "GRanges")
){
  # for testing: file_path = Khera_file
  
  
  # One-time Setup for Reading data, needed to get "skip" value
  #####
  # ?readLines
  # readLines(file_path, n = 20)
  
  # Read PGRS Data
  #####
  PGRS_data_readr = read_tsv(file_path, skip = skip_lines, n_max = n_max_lines)
  # dim(PGRS_data_readr); colnames(PGRS_data_readr)
  
  # keep only the columns named in colnames_key, renamed to the internal names
  PGRS_data_readr = PGRS_data_readr[, colnames(PGRS_data_readr) %in% colnames_key]
  match_cols = match(colnames(PGRS_data_readr), table = colnames_key )
  colnames(PGRS_data_readr) = names(colnames_key)[match_cols]
  # files with only effect_allele + other_allele (e.g. the PGS Catalog format, other_allele mapped
  # to A2) have no A1: the effect allele is A1, so (A1, A2) = (effect, other)
  if(is.null(PGRS_data_readr$A1)) PGRS_data_readr$A1 = PGRS_data_readr$effect_allele

  if(split_variant_str_TF){
  PGRS_data = 
    PGRS_data_readr %>% 
    tidyr::separate(variant, 
                    into = char_split_names,
                    sep = char_sep_regexp,
                    remove = FALSE #default:T remove original column
    )
  } else {
    PGRS_data = PGRS_data_readr
  }
  rm(PGRS_data_readr)
  
  if(!is.null(chr_prefix)){
    # automatically convert chr
    variant_missing_chr_char = !grepl(chr_prefix, PGRS_data$chr)
    if( any(variant_missing_chr_char) ){
      PGRS_data$chr[variant_missing_chr_char] = paste(chr_prefix, PGRS_data$chr[variant_missing_chr_char], sep = "")
    }
    # PGRS_data = mutate(PGRS_data, position_hg19 = as.integer(pos_chr))
    PGRS_data = mutate(PGRS_data, variant = paste(chr, pos, A1, A2, sep = ":"))
  }
  if(is.null(PGRS_data$variant)){
    PGRS_data = mutate(PGRS_data, variant = paste(chr, pos, A1, A2, sep = ":"))
  }

  # convert to GRanges  
  PGRS_data_GR = GRanges(
    seqnames = Rle(
      PGRS_data[[ 'chr' ]], NULL
    ),
    ranges = IRanges(
      start = PGRS_data[[ "pos" ]], 
      end =   PGRS_data[[ "pos" ]], 
      names = PGRS_data[[ "variant" ]]
    ),
    wgt = PGRS_data[[ "effect_weight" ]],
    effect_allele = PGRS_data[[ "effect_allele" ]],
    A1 = PGRS_data[[ "A1" ]],
    A2 = PGRS_data[[ "A2" ]],
    variant = PGRS_data[[ "variant" ]]
  )
  
  print("PGRS GRanges:")
  print(PGRS_data_GR)
  
  if(verbose_TF){
    print("memory used to load PGRS \n and create GRanges copy:")
    print( format(utils::object.size(PGRS_data) + utils::object.size(PGRS_data_GR), units = "MB") )
    print("total current memory usage (MB):")
    print( sum(gc()[, 2]) )
  }


  if( split_save_TF ){
  
  PGRS_data_GR_split_chr = split(PGRS_data_GR, seqnames(PGRS_data_GR))
  unique(seqnames(PGRS_data_GR))
  # Note that the Khera PGRS does not include X or Y chromosomes
  
  for(idx in 1:length(PGRS_data_GR_split_chr) ){
    
    PGRS_data_chr = PGRS_data_GR_split_chr[[idx]]
    chr_name = names(PGRS_data_GR_split_chr)[idx]
    
    out_file = file.path(out_folder, paste0( paste(out_name, chr_name, sep = "_"), ".rds"))
    
    saveRDS(object = PGRS_data_chr, 
            compress = "gzip", #this is the default
            file = out_file)
    
    print(names(PGRS_data_GR_split_chr)[idx])
    print(PGRS_data_chr)
  }
  
  }
  
  return_list = list()
  # choose return_type %in% c("NULL", "tibble", "data.table", "GRanges")
  
  if(return_type == "NULL"){
    return(NULL)
  } else {
    
    if( "tibble" %in% return_type ){
      return_list$tibble = PGRS_data
    }
    if( "data.table" %in% return_type ){
      return_list$data.table = as.data.table(PGRS_data)
    }
    if( "GRanges" %in% return_type ){
      return_list$GRanges = PGRS_data_GR
    }
    
    return(return_list)
    
  }
  
}
