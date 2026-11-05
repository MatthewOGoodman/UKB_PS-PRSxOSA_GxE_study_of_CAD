# T2 — the release scorer on the bgen project's example genotypes, against an independent
# reference computed from the example's text-format copy (example.gen).

toy <- read_toy()
tol <- 1e-5 * sum(abs(toy$wgt))   # 16-bit genotype probabilities: error scales with the summed weights

test_that("scores match the reference computed from example.gen", {
  res <- with_warnings(run_release(toy_granges(toy), file.path(tempfile(), "out")))
  ref <- reference_scores(toy)
  expect_equal(dim(res$value), dim(ref))
  expect_identical(rownames(res$value), rownames(ref))
  expect_identical(is.na(unname(res$value)), is.na(unname(ref)))
  expect_lt(max(abs(unname(res$value) - unname(ref)), na.rm = TRUE), tol)
  expect_true(any(grepl("matched by more than one variant", res$warnings)))
})

test_that("pathway names with ':', '/' and ',' are scored", {
  res <- with_warnings(run_release(toy_granges(toy), file.path(tempfile(), "out")))$value
  expect_true(any(grepl("Toy_pathway_with_punctuation", colnames(res))))
})

test_that("two PGRS rows at one position in a block stop before any genotype is read", {
  dup <- rbind(toy[toy$pos == 10000, ], toy[toy$pos == 10000, ])
  dup$A1[2] <- "C"
  expect_error(suppressWarnings(run_release(toy_granges(dup), file.path(tempfile(), "out"))),
               "share a position")
})

test_that("chromosome levels in a non-standard order are mapped by name", {
  gr <- toy_granges(toy, seqlevels_order = c("10", "2", "1"))
  res <- suppressWarnings(run_release(gr, file.path(tempfile(), "out")))
  expect_lt(max(abs(unname(res) - unname(reference_scores(toy))), na.rm = TRUE), tol)
})

test_that("a final block holding a single variant is scored", {
  res <- suppressWarnings(run_release(toy_granges(toy), file.path(tempfile(), "out"), blocksize = nrow(toy) - 1))
  expect_lt(max(abs(unname(res) - unname(reference_scores(toy))), na.rm = TRUE), tol)
})
