#!/usr/bin/env bash
# setup/build_rbgen.sh — build the rbgen R source package from a pinned bgen release, and unpack
# the bgen example files that the scorer example uses.
#
# rbgen (the R reader for .bgen files) is not on CRAN or Bioconductor. It is built from the bgen
# repository by the repository's own `waf build_rbgen` step. This script pins that repository to
# release v1.1.7 by checksum, runs the step, and records the one change the package needs to
# compile under current R (see "C++ standard" below). Compiling happens afterwards, in
# setup/setup_renv.R, with the R that owns the project library.
#
# Usage, from the project root:
#   bash setup/build_rbgen.sh          # default: work in ./.rbgen_build, examples to ./example/bgen
#
# Environment (all optional):
#   PYTHON             python 3 interpreter for bgen's waf build          (default: python3)
#   RBGEN_BUILD_DIR    work directory                                       (default: ./.rbgen_build)
#   RBGEN_EXAMPLE_DIR  destination for example.16bits.bgen/.bgi, example.gen, example.sample
#                                                                          (default: ./example/bgen)
#
# Requires: bash, curl, tar, gzip, python 3, and sha256sum or shasum. waf's configure step also
# looks for a C/C++ compiler and zlib.
#
# Output: $RBGEN_BUILD_DIR/rbgen_v1.1.7_cxx14.tar.gz, and the example files. Re-running reuses the
# verified download and rebuilds everything else from it.
set -euo pipefail

BGEN_TAG="v1.1.7"
BGEN_URL="https://enkre.net/cgi-bin/code/bgen/tarball/${BGEN_TAG}/bgen.tgz"
# The host generates this tarball on request; it is byte-identical across downloads, so the
# compressed file itself is pinned.
BGEN_SHA256="22a5e9254b71a08f233779350472f03c6f32a60cb6e07750305a9115454d7470"

PYTHON="${PYTHON:-python3}"
WORK="${RBGEN_BUILD_DIR:-$PWD/.rbgen_build}"
EXAMPLE_DIR="${RBGEN_EXAMPLE_DIR:-$PWD/example/bgen}"

# Linux ships sha256sum; macOS ships shasum.
if command -v sha256sum >/dev/null 2>&1; then
  SHA256=(sha256sum)
elif command -v shasum >/dev/null 2>&1; then
  SHA256=(shasum -a 256)
else
  echo "build_rbgen: no sha256sum or shasum on PATH" >&2
  exit 1
fi

for tool in curl tar gzip "$PYTHON"; do
  command -v "$tool" >/dev/null 2>&1 || { echo "build_rbgen: '$tool' not found on PATH" >&2; exit 1; }
done

mkdir -p "$WORK"
archive="$WORK/bgen_${BGEN_TAG}.tgz"

# --- 1. download and verify ------------------------------------------------
have=""
if [[ -f "$archive" ]]; then have="$("${SHA256[@]}" "$archive" | cut -d' ' -f1)"; fi
if [[ "$have" != "$BGEN_SHA256" ]]; then
  echo "== downloading bgen ${BGEN_TAG}"
  curl -fsSL -o "$archive.partial" "$BGEN_URL"
  mv "$archive.partial" "$archive"
fi
got="$("${SHA256[@]}" "$archive" | cut -d' ' -f1)"
if [[ "$got" != "$BGEN_SHA256" ]]; then
  echo "build_rbgen: checksum mismatch for $archive" >&2
  echo "  expected $BGEN_SHA256" >&2
  echo "  got      $got" >&2
  exit 1
fi
echo "== verified $archive"

# --- 2. unpack a fresh source tree -----------------------------------------
# The archive's top-level directory is named after the download ("bgen.tgz"), so it is located
# by the waf script it contains rather than by name.
rm -rf "$WORK/src" "$WORK/tmp" "$WORK/pkg"
mkdir -p "$WORK/src" "$WORK/tmp" "$WORK/pkg"
tar -xzf "$archive" -C "$WORK/src"
src_dir="$(dirname "$(find "$WORK/src" -maxdepth 2 -name waf -type f | head -n 1)")"
[[ -f "$src_dir/wscript" ]] || { echo "build_rbgen: no waf build in the unpacked archive" >&2; exit 1; }

# --- 3. bgen's own rbgen packaging step ------------------------------------
# build_rbgen writes the package tarball into python's tempfile directory, which follows TMPDIR;
# pointing TMPDIR into the work directory makes the output location known.
echo "== waf configure && waf build_rbgen (log: $WORK/waf.log)"
(
  cd "$src_dir"
  TMPDIR="$WORK/tmp" "$PYTHON" waf configure
  TMPDIR="$WORK/tmp" "$PYTHON" waf build_rbgen
) > "$WORK/waf.log" 2>&1 || { echo "build_rbgen: waf failed; see $WORK/waf.log" >&2; exit 1; }
pkg_tgz="$(find "$WORK/tmp" -name "rbgen_${BGEN_TAG}.tgz" -type f | head -n 1)"
[[ -n "$pkg_tgz" ]] || { echo "build_rbgen: waf produced no rbgen_${BGEN_TAG}.tgz; see $WORK/waf.log" >&2; exit 1; }

# --- 4. C++ standard -------------------------------------------------------
# R >= 4.3 compiles packages as C++17 unless the package asks otherwise. bgen's bundled Boost 1.55
# and its db/ layer use std::auto_ptr and std::unary_function, which C++17 removed, so the package
# is pinned to C++14, a standard R still supports. This line is the only change to bgen's code.
tar -xzf "$pkg_tgz" -C "$WORK/pkg"
makevars="$WORK/pkg/rbgen_${BGEN_TAG}/src/Makevars"
[[ -f "$makevars" ]] || { echo "build_rbgen: $makevars not found" >&2; exit 1; }
{ printf 'CXX_STD = CXX14\n'; cat "$makevars"; } > "$makevars.new"
mv "$makevars.new" "$makevars"
out_tgz="$WORK/rbgen_${BGEN_TAG}_cxx14.tar.gz"
tar -czf "$out_tgz" -C "$WORK/pkg" "rbgen_${BGEN_TAG}"

# --- 5. example files --------------------------------------------------------
mkdir -p "$EXAMPLE_DIR"
for f in example.16bits.bgen example.16bits.bgen.bgi example.gen example.sample; do
  cp "$src_dir/example/$f" "$EXAMPLE_DIR/$f"
done

echo "== rbgen source package : $out_tgz"
echo "== example files        : $EXAMPLE_DIR"
