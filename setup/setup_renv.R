#!/usr/bin/env Rscript
# setup/setup_renv.R — (re)build this project's renv library from ./DESCRIPTION.
#
# THE DECLARATION IS ./DESCRIPTION. renv's "explicit" snapshot type reads it, so the library is a
# function of a reviewable list rather than of whatever happens to be installed. Imports are what
# the scorer needs; Suggests are the annotation and pathway builders plus the tests.
#
# Five steps, run as ONE R session:
#   1. reset    renv::deactivate(clean = TRUE) removes any earlier renv state, so a re-run
#               converges on the declared set instead of accumulating on top of it.
#   2. init     renv::init(bare = TRUE): an empty project library.
#   3. hydrate  link every declared package the running R already has in its own library
#               (for a conda R, the environment's prebuilt binaries), so nothing is recompiled.
#   4. install  the rest: CRAN, Bioconductor (BIOC_PKGS below), and rbgen from the source package
#               setup/build_rbgen.sh builds.
#   5. snapshot renv::snapshot(type = "explicit") writes renv.lock.
#
# RUN from the project root, with the R that should own the library, under --vanilla:
#   bash setup/build_rbgen.sh
#   Rscript --vanilla setup/setup_renv.R
# --vanilla is required: step 1 deletes the project library, so renv must load from the running
# R's own library, and a user or project .Rprofile must not run first.
# A conda R must run inside its activated environment, because R compiles with the environment's
# own prefixed compilers:
#   conda run -p <env prefix> Rscript --vanilla setup/setup_renv.R
#
# Environment (optional):
#   RENV_HYDRATE_SOURCE  library to hydrate from (default: the running R's own library, .Library)
#   RBGEN_TARBALL        rbgen source package (default: the newest .rbgen_build/rbgen_*_cxx14.tar.gz)

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!file.exists("DESCRIPTION") || !file.exists(file.path("setup", "setup_renv.R")))
  stop("Run setup/setup_renv.R from the project root (it needs ./DESCRIPTION). getwd() = ", getwd(),
       call. = FALSE)
if (!requireNamespace("renv", quietly = TRUE))
  stop("renv is not installed in the running R. Install it first: install.packages(\"renv\")",
       call. = FALSE)

# --- declared dependencies ---------------------------------------------------
dcf   <- read.dcf("DESCRIPTION")
field <- function(f) if (f %in% colnames(dcf))
  setdiff(trimws(gsub("[(].*?[)]", "", strsplit(dcf[1, f], ",")[[1]])), c("", "R")) else character(0)
imports  <- field("Imports")
suggests <- field("Suggests")
deps     <- sort(unique(c(imports, suggests)))

# Routing only: where a declared package comes from. rbgen has its own source (step 4).
BIOC_PKGS <- c("GenomicRanges", "GenomeInfoDb", "IRanges", "S4Vectors", "rtracklayer",
               "GenomicFeatures", "txdbmaker", "VariantAnnotation", "InteractionSet", "AnnotationHub",
               "KEGGREST", "KEGGgraph", "org.Hs.eg.db", "AnnotationDbi")
LOCAL_PKG <- "rbgen"
spec_for  <- function(p) if (p %in% BIOC_PKGS) paste0("bioc::", p) else p

hydrate_src <- Sys.getenv("RENV_HYDRATE_SOURCE", .Library)

rbgen_tgz <- Sys.getenv("RBGEN_TARBALL", "")
if (!nzchar(rbgen_tgz)) {
  cand <- list.files(".rbgen_build", pattern = "^rbgen_.*_cxx14[.]tar[.]gz$", full.names = TRUE)
  rbgen_tgz <- if (length(cand)) cand[order(file.mtime(cand), decreasing = TRUE)][1] else ""
}
if (LOCAL_PKG %in% deps && (!nzchar(rbgen_tgz) || !file.exists(rbgen_tgz)))
  stop("No rbgen source package found. Run `bash setup/build_rbgen.sh` first, or set RBGEN_TARBALL.",
       call. = FALSE)

message("== R        : ", R.version.string, " (", R.home(), ")")
message("== declared : ", length(imports), " Imports + ", length(suggests), " Suggests")
message("== hydrate  : ", hydrate_src)

# --- 1. reset ------------------------------------------------------------------
message("== [1/5] reset")
if (dir.exists("renv") || file.exists("renv.lock")) renv::deactivate(clean = TRUE)

# --- 2. init -------------------------------------------------------------------
message("== [2/5] init (bare)")
renv::init(bare = TRUE, restart = FALSE)
renv::settings$snapshot.type("explicit", persist = TRUE)
renv::settings$snapshot.dev(TRUE, persist = TRUE)

# --- 3. hydrate ----------------------------------------------------------------
avail <- rownames(installed.packages(lib.loc = hydrate_src))
keep  <- setdiff(intersect(deps, avail), LOCAL_PKG)
message(sprintf("== [3/5] hydrate: %d of %d declared are prebuilt in %s", length(keep), length(deps), hydrate_src))
if (length(keep))
  renv::hydrate(packages = keep, library = renv::paths$library(), sources = hydrate_src,
                update = FALSE, prompt = FALSE)

# --- 4. install the rest ---------------------------------------------------------
# C17 COMPILER. Rgraphviz (required by KEGGgraph) requests the C17 standard, and its bundled
# graphviz sources do not compile under the newer C default. An R whose build configuration leaves
# CC17 empty (the conda-forge R does) refuses such a package. When that is the case, a Makevars file
# under renv/ defines CC17 as R's own C compiler in gnu17 mode. The compiler is written literally:
# R sets CC to $(CC17) while compiling, so a CC17 that referred to $(CC) would be circular.
r_bin <- file.path(R.home("bin"), "R")
r_cfg <- function(v) trimws(suppressWarnings(system2(r_bin, c("CMD", "config", v), stdout = TRUE))[1])
if (!nzchar(r_cfg("CC17"))) {
  if (nzchar(Sys.getenv("R_MAKEVARS_USER"))) {
    message("   note: this R defines no CC17 and R_MAKEVARS_USER is already set (",
            Sys.getenv("R_MAKEVARS_USER"), "); it must define CC17 for Rgraphviz to build")
  } else {
    mk <- file.path(getwd(), "renv", "Makevars.c17")
    writeLines(paste0("CC17 = ", r_cfg("CC"), " -std=gnu17"), mk)
    Sys.setenv(R_MAKEVARS_USER = mk)
    message("   CC17 undefined in this R; using ", mk, ": ", readLines(mk))
  }
}

have    <- rownames(installed.packages(lib.loc = renv::paths$library()))
missing <- setdiff(deps, c(have, LOCAL_PKG))
message(sprintf("== [4/5] install: %d from CRAN/Bioconductor, plus rbgen from %s", length(missing), rbgen_tgz))
if (length(missing)) renv::install(unname(vapply(missing, spec_for, "")), prompt = FALSE)
renv::install(normalizePath(rbgen_tgz), prompt = FALSE)

# --- 5. snapshot -------------------------------------------------------------------
message("== [5/5] snapshot (explicit, with Suggests)")
renv::snapshot(type = "explicit", dev = TRUE, prompt = FALSE, force = TRUE)
if (!file.exists("renv.lock")) stop("renv::snapshot() wrote no renv.lock; see the messages above.", call. = FALSE)

still <- setdiff(deps, rownames(installed.packages(lib.loc = renv::paths$library())))
if (length(still)) stop("Declared but not installed: ", paste(still, collapse = ", "), call. = FALSE)
message("\n== done: every declared package is in ", renv::paths$library())
