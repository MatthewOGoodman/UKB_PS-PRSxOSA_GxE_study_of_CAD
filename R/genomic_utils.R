# LIFTOVER FUNCTIONS -----

# References
# http://genome.ucsc.edu/cgi-bin/hgLiftOver
# http://hgdownload.soe.ucsc.edu/downloads.html#human

# Examples
# ?rtracklayer::liftOver

# Get chain file from annotation hub -----

get_chain = function(hg_in = hg_input, hg_out = hg_output){
  # function to acquire chain file from AnnotationHub for liftover
  # AnnotationHub is loaded here, not at file top, so that sourcing this file for the scorer does
  # not require it. Only the GeneHancer liftover calls get_chain.
  library(AnnotationHub)
  ahub <- AnnotationHub()
  ahub.chain <- subset(ahub, rdataclass == "ChainFile" & species == "Homo sapiens")
  querychains = query(ahub.chain, c(hg_in, hg_out))
  # the UCSC liftOver chain, e.g. hg38ToHg19.over.chain.gz; other chain records (e.g. an rRNA-only
  # chain) also mention both builds, so the file-name pattern is matched exactly
  chain_title = querychains$title[grep(paste0("^", hg_in, "To", hg_out, "\\.over\\.chain"), querychains$title , ignore.case = TRUE)]
  chain_ahub_name = names(querychains)[ querychains$title %in% chain_title  ]
  if(length(chain_ahub_name) != 1){
    stop("expected one AnnotationHub chain for ", hg_in, " to ", hg_out, ", found: ", paste(chain_title, collapse = ", "))
  }
  chain = ahub[[chain_ahub_name]]
}

# UTILITY FUNCTIONS ----
###

# sub-group a grouping variable
###
chunk_PGRS_within_chr = function(PGRS, chunk_size = 100){
  # the purpose of chunk_chrs is to split a parallel calculation into chunks (smaller groups than the chr level)
  # chrs is a vector with repeats containing the chr for each variant being processed
  # chunk_size is the maximum number of variants in each chunk. 
  # chunks respect chromosomes, so there will be smaller chunks at the end of each chromosome.
  # PGRS must be sorted by chromosome: chunk ids follow the order of the variants.

  # setup
  ####
  # drop chromosome levels with no variants; the per-chromosome start/end arithmetic below
  # assumes every level is populated
  PGRS = GenomeInfoDb::keepSeqlevels(PGRS, unique(as.character(seqnames(PGRS))))
  n = length(PGRS)
  mcols(PGRS)$chromosome = seqnames(PGRS)
  chrs = tibble(chr = as.factor(seqnames(PGRS)) )
  
  PGRS_list <- split(PGRS, seqnames(PGRS))
  n_chr = length(PGRS_list)
  
  # get start and end of each chromosome
  ####
  chr_lengths = tibble( .rows = n_chr )
  chr_lengths$lengths = sapply(PGRS_list, length)
  chr_lengths$start_pos = c(1, 1+cumsum(chr_lengths$lengths[-n_chr]) )
  chr_lengths$end_pos = cumsum(chr_lengths$lengths)
  
  # get start and end of each chunk
  ####
  # method 1: handy function to convert array (rows of tbl data frame ) to list elements
  chr_lengths_list = array_tree(array = chr_lengths, margin = 1)
  chunk_start_pos_in_chr = map(chr_lengths_list, ~ seq(1, .x['lengths'], chunk_size) )
  
  # method 2: more straightforward but less flexible approach
  # chr_lengths_list = as.list(chr_lengths$lengths)
  # chunk_start_pos_in_chr = map(chr_lengths_list, ~ seq(1, .x, chunk_size) )
  
  prior_length = as.list(chr_lengths$start_pos - 1)
  
  chunk_starts = unlist(map2( prior_length, chunk_start_pos_in_chr, ~.x + .y))
  chunk_ends = c( (chunk_starts - 1)[-1], chr_lengths$end_pos[n_chr] )
  chunk_sizes = chunk_ends - chunk_starts + 1
  n_chunks = length(chunk_starts)
  chrs$chunks = rep(1:n_chunks, times = chunk_sizes)
  # length(chunks); length(chrs$chr)
  
  # sanity checks:
  PGRS$chunks = chrs$chunks
  chr_per_chunk = tapply(PGRS$chromosome, PGRS$chunks, function(x) unique(x))
  
  print(sprintf("there is a maximum of %s chromosome(s) per chunk", max( sapply(chr_per_chunk, length)) ))
  print("here is the assignment of chromosomes to chunks:")
  print(chr_per_chunk)
  
  if( max( sapply(chr_per_chunk, length)) >1){ stop("please make sure that chunks respect chromosomes, allow only one chromosome per chunk") }
  
  return(chrs)
}


# duplicated rows
####
get_dup_rows = function(data, columns = NULL, index_TF = FALSE){
  if(is.null(columns)) {columns = 1:ncol(data)}
  data_col = data[,columns]
  logical = duplicated(data_col) | duplicated(data_col[nrow(data_col):1, ])[nrow(data_col):1]
  if(index_TF){
    index = which(logical)
    return(index)
  } else {
    return(logical)
  }
}

view_dup_rows = function(data, columns = NULL){
  if(is.null(columns)){columns = 1:ncol(data)}
  data[get_dup_rows(data, columns),]
}

# get unique elements from selected rows and columns of matrix, excluding NA and 'exclude'
####
unique_from_columns = function(dataframe, rows, columns, exclude = "?"){
  data = c( t(dataframe[rows,columns]) )
  data = unique(data)
  data = data[!is.na(data)]
  data = data[! data %in% exclude]
}
