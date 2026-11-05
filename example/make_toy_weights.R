# example/make_toy_weights.R — build example/toy_weights.tsv, a toy PRS with pathway membership
# for the bgen project's example genotypes (example/bgen/example.gen; 199 synthetic variants on
# chromosome "01", 500 synthetic samples).
#
# Run from the repository root:  Rscript example/make_toy_weights.R
# The output is committed; re-running with the same seed reproduces it.
#
# The table exercises the scorer's allele handling:
#   - one row per position: position 10000 carries two A/G variants in the bgen file, and the
#     single PRS row there is matched (and counted) by both, with a warning
#   - rows whose A1/A2 order is the reverse of the bgen allele0/allele1 order
#   - rows whose alleles match neither bgen allele (weight set to 0 by the scorer)
#   - an effect allele that is A1 on most rows and A2 on some
#   - pathway columns, one of them with ":", "/" and "," in its name

set.seed(2026)

gen <- read.table("example/bgen/example.gen", stringsAsFactors = FALSE)[, 1:6]
colnames(gen) <- c("chr", "snpid", "rsid", "pos", "allele0", "allele1")
gen <- gen[!duplicated(gen$pos), ]
n <- nrow(gen)

toy <- data.frame(chr = "1", pos = gen$pos, A1 = gen$allele0, A2 = gen$allele1, stringsAsFactors = FALSE)

swap <- sample(n, 10)
toy[swap, c("A1", "A2")] <- toy[swap, c("A2", "A1")]

mismatch <- sample(setdiff(seq_len(n), swap), 5)
toy$A1[mismatch] <- "C"
toy$A2[mismatch] <- "T"
stopifnot(all(!(toy$A1[mismatch] %in% c(gen$allele0[mismatch], gen$allele1[mismatch]) &
                toy$A2[mismatch] %in% c(gen$allele0[mismatch], gen$allele1[mismatch]))))

toy$effect_allele <- toy$A1
ea_a2 <- sample(seq_len(n), 20)
toy$effect_allele[ea_a2] <- toy$A2[ea_a2]

toy$wgt <- round(rnorm(n, 0, 0.1), 6)

toy$all_SNPs <- TRUE
toy$toy_pathway_A <- runif(n) < 0.5
toy$toy_pathway_B <- seq_len(n) <= 50
toy[["hsa99999: Toy/pathway, with punctuation"]] <- runif(n) < 0.3

write.table(toy, "example/toy_weights.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote example/toy_weights.tsv:", nrow(toy), "rows\n")
