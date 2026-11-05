

# Pathway names become column names of the score matrix. The same cleaning is applied to the
# pathway names and to the joined variant table's column names, so that every pathway column is
# found whatever characters its name contains (":", "/", ",", "(", "-", spaces).
sanitize_pathway_names = function(x){
  x = make.names(x)
  x = str_replace_all(x, "\\.+", "\\_")
  x = str_replace_all(x, "\\_$", "")
  x
}

compute_PS_PGRS = function(    PGRS_chunk,
                               PGRS_chunk_name = "chr1",
                               pathway_columns,
                               # note: PGRS_chunk must have columns with pathway annotation columns with T/F inclusion exclusion for each SNP
                               genotype_dir,
                               genotype_filename_fn = function(chr_name) paste0("ukb_imp_", chr_name, "_v3.bgen"),
                               # note: this function assumes that chr_name derived from the PGRS seqnames column is the basis for the filenames in the genotype_dir
                               out_dir,
                               run_id,
                               blocksize = 1000,
                               subset_samples = NULL,
                               # column of PGRS_chunk holding the allele whose dosage the weight multiplies
                               effect_allele_col = "effect_allele",
                               # maps the PGRS chromosome names to the names used inside the bgen file
                               # (UK Biobank imputed v3 uses "01".."22")
                               chr_label_fn = function(chr) { n = sub("^chr", "", chr); ifelse(grepl("^[0-9]$", n), paste0("0", n), n) },
                               verbose_TF = FALSE){

  # set params ----
  ####
  # Block scores are checkpoints: each block is saved under the run directory as it completes, a
  # re-run of the same call skips the blocks already saved, and the chunk summary is written only
  # when every block is present. Block files hold per-sample scores, so the block directory is
  # removed once the summary is written.
  out_dir_job = file.path(out_dir, run_id)
  block_dir = file.path(out_dir_job, "blocks", PGRS_chunk_name)
  dir.create(block_dir, recursive = TRUE, showWarnings = FALSE)

  # input stamp: a block directory may only be resumed by the identical input
  stamp_tmp = tempfile(fileext = ".rds")
  saveRDS(list(as.data.frame(PGRS_chunk), pathway_columns, blocksize, effect_allele_col, subset_samples),
          stamp_tmp, compress = FALSE)
  stamp = unname(tools::md5sum(stamp_tmp))
  unlink(stamp_tmp)
  stamp_file = file.path(block_dir, "STAMP")
  if(file.exists(stamp_file)){
    if(!identical(readLines(stamp_file, warn = FALSE), stamp))
      stop("block directory ", block_dir, " belongs to a different input; use a new run_id or remove it")
  } else {
    writeLines(stamp, stamp_file)
  }

  # get chr name as UCSC style to correspond with ukb bgen data files
  # seqlevels(PGRS_chunk) = str_replace( seqlevels(PGRS_chunk), pattern = "^0", replacement = "") 
  
  seqlevelsStyle(PGRS_chunk) <- "UCSC"
  chr_name = unique(seqnames(PGRS_chunk))
  # warn if chromosome in chunk is not unique
  if( length(chr_name) > 1){
    stop("Chunks must respect chromosomes. It appears you are requesting data from different chromosomes")
  }
  
  
  # summary stats
  n_alleles_not_matched = 0
  n_alleles_counted_twice = 0
  
  # risk score params
  nSNP = length(PGRS_chunk)
  n_blocks = nSNP %/% blocksize
  n_last_block = nSNP %% blocksize
  # n_blocks * blocksize + n_last_block == nSNP

  # the join to the bgen variants within a block is by position, so two PGRS rows at one
  # position in the same block would each pick up the other's genotypes; checked for every block
  # before any genotypes are read
  block_of_row = (seq_len(nSNP) - 1) %/% blocksize
  dup_pos = duplicated(data.frame(block_of_row, as.character(seqnames(PGRS_chunk)), start(PGRS_chunk)))
  if(any(dup_pos)){
    stop("PGRS rows share a position within a block of chunk ", PGRS_chunk_name, ": ",
         paste(unique(paste(seqnames(PGRS_chunk), start(PGRS_chunk), sep = ":")[dup_pos]), collapse = ", "),
         ". Keep one row per position, or use a blocksize that separates them.")
  }
  absent_paths = !(sanitize_pathway_names(pathway_columns) %in% sanitize_pathway_names(names(mcols(PGRS_chunk))))
  if(any(absent_paths)){
    stop("pathway columns not in PGRS_chunk: ", paste(pathway_columns[absent_paths], collapse = ", "))
  }

  # re-set chromosome strings in PGRS to match variant locations in the bgen file
  # (rbgen does not use a standard seqlevelsStyle). Levels are renamed by name, after dropping
  # unused levels, so the result does not depend on the order of the levels.
  PGRS_chunk = keepSeqlevels(PGRS_chunk, unique(as.character(seqnames(PGRS_chunk))))
  PGRS_chunk = renameSeqlevels(PGRS_chunk, setNames(chr_label_fn(seqlevels(PGRS_chunk)), seqlevels(PGRS_chunk)))

  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # iterate over blocks to calculate PGRS ----
  ####
  
  iters = 0:n_blocks

  for(iter in iters){
    
    try({
      
      start = Sys.time()
      cat("\n\n")
      print("#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#")
      cat("\n\n")
      
      print(sprintf( "Iteration: %i", iter))
      
      block_filename = paste0(paste("PathSpecPGRS", "chunk", PGRS_chunk_name, "block", iter, sep = "_"), ".rds")

      if( file.exists( file.path(block_dir, block_filename)) ){
        print(sprintf("block %i already saved; skipping", iter))
        next # jump to next iteration if the output exists
      }

      # set PGRS block
      if(iter < n_blocks ){
        PGRS_block = PGRS_chunk[ (iter*blocksize + 1:blocksize ),] 
      } else if (iter == n_blocks & n_last_block != 0){ 
        PGRS_block = PGRS_chunk[ (iter*blocksize + 1:n_last_block ),]
      } else if (n_last_block == 0){
        break
      }
      PGRS_block = as_tibble(PGRS_block)
      colnames(PGRS_block)[1] = "chromosome"
      PGRS_block_ranges = as.data.frame(PGRS_block[,c("chromosome", "start", "end")])
      
      # load data
      ####
      chunk_file = file.path(genotype_dir, genotype_filename_fn(chr_name)) #ukb_imp_chr1_v3.bgen
      
      ukb3_bgen = bgen.load( filename = chunk_file, 
                             # index.filename = paste0(chunk_file,".bgi"), # not needed as default: sprintf( "%s.bgi", filename) 
                             ranges   = PGRS_block_ranges,
                             samples  = subset_samples
      )
      rm(PGRS_block_ranges)
      
      runtime = Sys.time() - start

      if(verbose_TF){
        print("time to load current SNP data block:")
        print(runtime)
        print("memory size of SNP data block:")
        print(format(utils::object.size(ukb3_bgen), units = "MB"))
        print("total current memory used:")
        print(gc())
        print("info on data loaded:")
        print(str(ukb3_bgen)) #note: this prints a NULL for some reason
      }
      
      
      # alleles need to be checked for matched pair, and flipped if reversed ----
      ####
      
      # to overcome duplicate rownames convert to tibble
      ukb3_bgen$variants = as_tibble(ukb3_bgen$variants)
      # (duplicate variant names/SNPids, due to multi-allelic)
      # rownames(ukb3_bgen$variants)[duplicated(rownames(ukb3_bgen$variants))]
      
      # left-join output (variants), as reference, with input (PGRS) 
      variants = left_join(ukb3_bgen$variants, PGRS_block, by = c("chromosome" = "chromosome", "position" = "start"))
      
      ukb3_bgen$variants = NULL; gc() # to save memory
      
      n_var_block = nrow(variants)
      # dim(variants)
      # view_dup_rows(variants, columns = c("chromosome", "position"))[,c("chromosome", "position", "rsid", "allele0", "allele1", "A1", "A2", "wgt")]
      # note that there are duplicate rows, when UKB data has multiple alleles for the given rsid or position
      # note that wgt is repeated across the duplicate rows, only one of which matches PGRS alleles. This will be corrected by setting the wgt and genotypes to 0 below.
      
      # which PGRS alleles are reversed?
      # note that same_order implies same alleles in same order, and reverse_order implies same alleles in reverse order
      same_alleles  = apply(variants[,c("allele0", "allele1", "A1", "A2")], 1, function(a){ w = a[1]; x = a[2]; y = a[3]; z = a[4]; setequal(c(w,x), c(y,z)) })
      same_order    = apply(variants[,c("allele0", "allele1", "A1", "A2")], 1, function(a){ w = a[1]; x = a[2]; y = a[3]; z = a[4]; all(c(w,x) == c(y,z)) }) 
      reverse_order = apply(variants[,c("allele0", "allele1", "A1", "A2")], 1, function(a){ w = a[1]; x = a[2]; y = a[3]; z = a[4]; all(c(w,x) == c(z,y)) })
      diff_alleles  = !same_order & !reverse_order

      variants[diff_alleles,"wgt"] = 0

      # the dosage counts allele1 unless the row is reversed; reverse the rows whose effect allele
      # is allele0, so that the dosage always counts the effect allele
      count_allele0 = as.character(variants$allele0) == as.character(variants[[effect_allele_col]]) & !diff_alleles
      count_allele0[is.na(count_allele0)] = FALSE

      # increment n_alleles_not_matched
      if(sum(same_alleles) < nrow(PGRS_block) ){
        warning("In iter ", iter, " there were ", nrow(PGRS_block)-sum(same_alleles), " alleles in the PGRS that were not matched" )
        n_alleles_not_matched = n_alleles_not_matched + nrow(PGRS_block)-sum(same_alleles)
      }
      # increment n_alleles_counted_twice: one PGRS row matched (in either orientation) by more
      # than one bgen variant, e.g. two bgen records at one position with the same alleles
      matched = same_order | reverse_order
      counted_twice = duplicated(variants[matched, c("chromosome", "position", "A1", "A2")])
      if( any(counted_twice) ){
        warning("In iter ", iter, " there were ", sum(counted_twice), " PGRS SNPs matched by more than one variant in the data; each match is counted" )
        n_alleles_counted_twice = n_alleles_counted_twice + sum(counted_twice)
      }
      
      # checks:
      ####
      # cbind(same_alleles, same_order, reverse_order, diff_alleles)
      # all( same_alleles ) # this will be false with multi-allelic SNPs in the data
      # all(variants$position == variants$end)
      # sum(same_alleles+diff_alleles) == nrow(variants)
      

      # reduce alleles to univariate 0,1,2 scale
      ####
      
      print("total memory used prior to reducing alleles:")
      print(gc())
      
      ukb3_bgen$data = t( reduce_bgen_allele_data(data = ukb3_bgen$data, reverse = count_allele0) )

      print("total memory used after reducing alleles:")
      print(gc())

      dim(ukb3_bgen$data); dim(variants)

      # fix problems with column names changing during left-join
      colnames(variants) = sanitize_pathway_names(colnames(variants))
      pathway_columns_block = sanitize_pathway_names(pathway_columns)
      if(!all(pathway_columns_block %in% colnames(variants))){
        stop("pathway columns not found after name cleaning: ",
             paste(pathway_columns[!(pathway_columns_block %in% colnames(variants))], collapse = ", "))
      }

      
      # get subset risk scores for each subset column in PGRS_subsets ----
      ####
      variants_subsets = as.matrix(1*variants[,pathway_columns_block])
      variants_wgt = variants$wgt
      scores = ukb3_bgen$data %*% (variants_subsets * variants_wgt)

      # save scores block (written to a temporary name first, so a block file is always complete)
      block_path = file.path(block_dir, block_filename)
      saveRDS(scores, file = paste0(block_path, ".partial"))
      file.rename(paste0(block_path, ".partial"), block_path)
      
      rm(scores)
      rm(ukb3_bgen, variants, variants_subsets, variants_wgt)
      
      if(verbose_TF){
        print(sprintf("Time to complete one iteration for block %i: %s", iter,  format(Sys.time() - start, digits = 2) ))
        print(sprintf("Memory usage at end of iteration for block %i: %s MB", iter, round(sum(gc()[, 2]), 2) ))
      } else {
        gc()
      }
      
    }) # end try
  }  # end iters over blocks
  
  # chunk summary: written only when every block is present ----
  ####
  expected_iters = if(n_last_block == 0) seq_len(n_blocks) - 1 else 0:n_blocks
  block_files = file.path(block_dir, paste0(paste("PathSpecPGRS", "chunk", PGRS_chunk_name, "block", expected_iters, sep = "_"), ".rds"))
  missing_iters = expected_iters[!file.exists(block_files)]
  if( length(missing_iters) > 0 ){
    stop("chunk ", PGRS_chunk_name, ": block(s) ", paste(missing_iters, collapse = ", "),
         " missing (see the errors above). The completed blocks are saved in ", block_dir,
         "; re-run the same call to compute only the missing blocks.")
  }
  
  chunk_filename = paste0(paste("PathSpecPGRS", "chunk_summary", PGRS_chunk_name, sep = "_"), ".rds")
  summary_path = file.path(out_dir_job, chunk_filename)
  chunk_scores = Reduce(`+`, lapply(block_files, readRDS))
  saveRDS(chunk_scores, file = paste0(summary_path, ".partial"))
  file.rename(paste0(summary_path, ".partial"), summary_path)
  unlink(block_dir, recursive = TRUE)
  
  invisible(summary_path)
}
