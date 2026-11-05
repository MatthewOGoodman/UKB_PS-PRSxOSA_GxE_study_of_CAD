# fns_loadGenotypes.R

####
# Reduce genotype probabilities (p(AA),p(AB),p(BB)) to matrix of expected genotypes on the (0,1,2) additive scale ----
####

reduce_bgen_allele_data = function(data, reverse=NULL){
  # input data is in rbgen array format: SNPs x Samples x Genotypes, where
  #   Genotypes is the slices of the array, corresponding to p(AA), p(AB), and p(BB)
  # output has SNPs x Samples (p x n) dimensions
  
  #reverse the alleles e.g. to match PGRS:
  if( !is.null(reverse) & length(reverse > 0) ){
    data[reverse , , 1:3] = data[reverse , , 3:1]
  }
  #combine the levels to create mean effect-allele count
  # drop = FALSE keeps the SNPs x Samples shape when a block holds a single SNP (or a single sample)
  dn = dimnames(data)[1:2]
  data = 0 + data[,,2, drop = FALSE] + 2*data[,,3, drop = FALSE]
  dim(data) = dim(data)[1:2]
  dimnames(data) = dn
  data
}
