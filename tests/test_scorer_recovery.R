# T4 — block checkpoints: a failed block stops the chunk without a summary, a re-run of the same
# call completes it, a block directory is never reused by a different input, and two runs in one
# session do not share blocks.

toy <- read_toy()

# a release environment whose bgen.load fails for the block containing position `fail_pos`
failing_env <- function(fail_pos) {
  e <- release_env()
  e$bgen.load <- function(filename, ranges, ...) {
    if (any(ranges$start == fail_pos)) stop("simulated read failure")
    rbgen::bgen.load(filename = filename, ranges = ranges, ...)
  }
  e
}

test_that("a failed block stops the chunk, keeps the other blocks, and a re-run completes it", {
  out <- file.path(tempfile(), "out")
  fail_pos <- toy$pos[120]                      # row 120 is in block 2 at blocksize 50
  expect_error(suppressWarnings(run_release(toy_granges(toy), out, env = failing_env(fail_pos))),
               "block\\(s\\) 2 missing")
  block_dir <- file.path(out, "test", "blocks", "1")
  expect_setequal(list.files(block_dir, pattern = "rds$"),
                  sprintf("PathSpecPGRS_chunk_1_block_%d.rds", c(0, 1, 3)))
  expect_false(file.exists(file.path(out, "test", "PathSpecPGRS_chunk_summary_1.rds")))

  resumed <- suppressWarnings(run_release(toy_granges(toy), out))
  clean   <- suppressWarnings(run_release(toy_granges(toy), file.path(tempfile(), "out")))
  expect_identical(resumed, clean)
  expect_false(dir.exists(block_dir))
})

test_that("a block directory is not reused by a different input", {
  out <- file.path(tempfile(), "out")
  expect_error(suppressWarnings(run_release(toy_granges(toy), out, env = failing_env(toy$pos[120]))))
  changed <- toy
  changed$wgt <- changed$wgt * 2
  expect_error(suppressWarnings(run_release(toy_granges(changed), out)), "different input")
})

test_that("two runs in one session do not share blocks", {
  out <- file.path(tempfile(), "out")
  changed <- toy
  changed$wgt <- -changed$wgt
  a <- suppressWarnings(run_release(toy_granges(toy), out, run_id = "A"))
  b <- suppressWarnings(run_release(toy_granges(changed), out, run_id = "B"))
  expect_equal(b, -a)
})
