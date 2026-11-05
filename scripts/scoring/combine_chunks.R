# scripts/scoring/combine_chunks.R — sum the chunk summaries of one scoring run into the final
# samples x pathways score matrix.
#
# Run from the repository root after every chunk of the run has been scored:
#   Rscript scripts/scoring/combine_chunks.R
# Required environment: PSPRS_OUT_DIR, PSPRS_RUN_ID.
# Output: <PSPRS_OUT_DIR>/<PSPRS_RUN_ID>/PathSpecPGRS_total_scores.rds

source("config/paths.R")
psprs_require(c("psprs_out_dir", "psprs_run_id"))

out_dir_job = file.path(psprs_out_dir, psprs_run_id)
run_manifest = readRDS(file.path(out_dir_job, "run_manifest.rds"))

# every chunk of the run must have its summary; a chunk summary exists only once all of that
# chunk's blocks were scored
chunk_files = file.path(out_dir_job, paste0("PathSpecPGRS_chunk_summary_", seq_len(run_manifest$total_chunks), ".rds"))
missing_chunks = seq_len(run_manifest$total_chunks)[!file.exists(chunk_files)]
if(length(missing_chunks) > 0){
  stop("run ", psprs_run_id, ": chunk summaries missing for chunk(s) ", paste(missing_chunks, collapse = ", "),
       "; score those chunks first")
}
print(head(chunk_files))

for( iter in 1:length(chunk_files)){


  if(iter == 1){
    final_scores = readRDS(file = chunk_files[iter])
  } else {
    chunk_scores = readRDS(file = chunk_files[iter])
    final_scores = final_scores + chunk_scores

    rm(chunk_scores)
  }

  print(sprintf("loading iteration %s", iter))

}

print("dim(final_scores)")
print(dim(final_scores))

PGRS_filename = "PathSpecPGRS_total_scores.rds"
print(PGRS_filename)

saveRDS(final_scores, file = file.path(out_dir_job, PGRS_filename ) )
