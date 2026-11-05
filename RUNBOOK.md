# RUNBOOK

Ordered procedures for this repository. Run every command from the repository root.

## Setting file locations

The scripts contain no file paths. Before you run a step, tell it where your files are by
setting environment variables in the shell you run it from:

```sh
export PSPRS_DATA_DIR=/path/to/prs_and_outputs
export PSPRS_GENOTYPE_DIR=/path/to/ukb/imputed_v3
Rscript scripts/scoring/score_chunk.R
```

On a cluster, put the `export` lines in the job script, or pass them to the scheduler (SLURM:
`sbatch --export=ALL,PSPRS_DATA_DIR=...`).

`config/paths.R` reads these variables into the names the scripts use (for example
`PSPRS_GENOTYPE_DIR` becomes `genotype_dir`). If a step needs a variable you have not set, it
stops with "Set the environment variable(s): …" and names the variable to set. The table lists
every variable and the steps that use it.

| Variable | Used by | Meaning |
|---|---|---|
| `PSPRS_DATA_DIR` | §3–§5 | PRS weights file in; annotation and pathway-membership files out |
| `PSPRS_GENCODE_GFF` | §3, §4 | GENCODE GFF3 (GRCh37 lift) |
| `PSPRS_GTEX_DIR` | §3 | directory of GTEx `*.signif_variant_gene_pairs.txt.gz` files |
| `PSPRS_FANTOM5_BED` | §3 | FANTOM5 `enhancer_tss_associations.bed` |
| `PSPRS_PCHIC_DIR` | §3 | directory holding `ActivePromoterEnhancerLinks.tsv` |
| `PSPRS_GENEHANCER_CSV` | §3 | GeneHancer csv (optional; your own licensed copy) |
| `PSPRS_PATHWAYS_DIR`, `PSPRS_KEGG_DIR` | §4 | pathway gene lists out; KGML download cache |
| `PSPRS_ANNOTATION_DATE` | §3, §4 | date naming the annotation outputs (§3 defaults to today; §4 requires it) |
| `PSPRS_GENOTYPE_DIR`, `PSPRS_OUT_DIR`, `PSPRS_RUN_ID`, `PSPRS_TASK_ID` | §5, §6 | `.bgen` files; scoring output root; run name; chunk number |
| `PSPRS_ANALYSIS_DATA_DIR` | §7 | UK Biobank analysis data set |

Optional settings are listed at the top of each script.

## §1 Set up R and rbgen

1. Use an R ≥ 4.3 with a C/C++ compiler, `make` and zlib, plus python 3, curl and the `renv` package.
2. `bash setup/build_rbgen.sh`
   - Downloads bgen v1.1.7 from `https://enkre.net/cgi-bin/code/bgen/tarball/v1.1.7/bgen.tgz` and checks its SHA-256.
   - Builds the rbgen source package with bgen's own `waf build_rbgen` step.
   - Sets the package to compile as C++14, which current R needs for bgen's bundled Boost.
   - Unpacks the example genotypes to `example/bgen/`.
   - Set `PYTHON` if `python3` is not the interpreter to use.
3. `Rscript --vanilla setup/setup_renv.R`
   - Builds the project library from `./DESCRIPTION`, reusing packages the running R already has.
   - Installs rbgen from `.rbgen_build/`.
   - With a conda R, run it inside the environment: `conda run -p <env> Rscript --vanilla setup/setup_renv.R`. R compiles with the environment's own compilers.
   - If your R defines no C17 compiler, the script supplies one (needed by Rgraphviz, used by KEGGgraph).
4. Check: `Rscript -e 'testthat::test_dir("tests", stop_on_failure = FALSE)'`. The example tests must pass.

## §2 Run the example

`Rscript example/run_example.R`. This scores `example/toy_weights.tsv` on the bgen example genotypes (chromosome `"01"`, 500 samples). It uses §5 and §6's scripts, and writes to `tests/_output/example/`.

## §3 Annotate PRS variants to genes

1. Get the sources (the versions used for the paper). Each file is available from its original host and, except GeneHancer, from the project's Zenodo deposit ([doi:10.5281/zenodo.23990992](https://doi.org/10.5281/zenodo.23990992)), which holds the copies used for the paper.
   - GENCODE: `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.annotation.gff3.gz`, or `gencode.v34lift37.annotation.gff3.gz` on Zenodo → `PSPRS_GENCODE_GFF`.
   - GTEx v7 eQTL: `https://storage.googleapis.com/adult-gtex/bulk-qtl/v7/single-tissue-cis-qtl/GTEx_Analysis_v7_eQTL.tar.gz` (unpack it), or the six `*.v7.signif_variant_gene_pairs.txt.gz` files on Zenodo. Point `PSPRS_GTEX_DIR` at the directory holding them. The builder uses the artery, heart and whole-blood tissues.
   - FANTOM5: `http://enhancer.binf.ku.dk/presets/enhancer_tss_associations.bed` (offline as of November 2026), or `enhancer_tss_associations.bed` on Zenodo → `PSPRS_FANTOM5_BED`.
   - Promoter capture Hi-C: download "Data S1" from the article page (Javierre et al. 2016, *Cell*, https://doi.org/10.1016/j.cell.2016.09.037) and extract `ActivePromoterEnhancerLinks.tsv` from it, or download `ActivePromoterEnhancerLinks.tsv` from Zenodo. Point `PSPRS_PCHIC_DIR` at the directory holding it.
   - GeneHancer (optional): request it from LifeMap Sciences (GeneCards). The paper used `GeneHancer_version_4-4.csv` → `PSPRS_GENEHANCER_CSV`. Leave the variable unset to skip the source.
   - PRS weights: put the weights file in `PSPRS_DATA_DIR`.
     - The paper's file, `CoronaryArteryDisease_PRS_LDpred_rho0.001_v3.txt`, is the default.
     - For the PGS Catalog copy, download `https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000013/ScoringFiles/PGS000013.txt.gz`, then set `PSPRS_PRS_FILE=PGS000013.txt.gz` and `PSPRS_PRS_FORMAT=pgs_catalog`. PGS000013 holds the same variants, alleles and weights.
2. `Rscript scripts/annotation/run_annotation.R`
   - The GeneHancer liftover (GRCh38 → GRCh37) fetches the UCSC chain through AnnotationHub, so it needs network access.
   - `PSPRS_CHROMOSOMES=22` annotates one chromosome, for a quick check.
3. Outputs in `PSPRS_DATA_DIR`:
   - one `gene_annot_<date>_<source>_<PRS>.rds` per source;
   - `gene_annot_<date>_combined_annotation_sources_<PRS>.rds`;
   - `gene_annot_<date>_sequential.annotated.SNPs.to.ENSGs_<PRS>.rds`, the source-priority cascade, with each variant assigned by the first source that annotates it.
4. To skip step 2 for the paper's variants, use the paper's per-source annotation from Zenodo (`PSPRS_annotation_GRCh37_20201209.rds`). It covers every source except GeneHancer.
   - Copy it to `PSPRS_DATA_DIR` as `gene_annot_20201209_combined_annotation_sources_PGRS_CAD_Khera_LDPred_rho0.001.rds`.
   - Run `PSPRS_ANNOTATION_STEPS=cascade PSPRS_ANNOTATION_DATE=20201209 Rscript scripts/annotation/run_annotation.R`.
   - `PSPRS_ANNOTATION_STEPS=combine` instead re-combines existing per-source files.

## §4 Pathway gene lists and membership

1. Set `PSPRS_ANNOTATION_DATE` to the date of the §3 outputs.
2. `PSPRS_KEGG_PATHWAYS` chooses the pathways:
   - `curated` (default, the paper's method): the KEGG pathways listed in the curated sheet `data/OSA_Pathways_20201117_freeze.csv` (override with `PSPRS_CURATED_CSV`), the curated gene set, their union (`combined_pathway`), and all annotated genes.
   - `all`: every human KEGG pathway.
   - `hsa04066,hsa04370`: the listed pathways.
3. `PSPRS_GENE_SETS` (`NFkB_TNF` default; `HIF1_VEGF`; both; `none`) adds single-gene sets that split two pathways into genes found only in the first, only in the second, or in both. The paper made two runs, one per group.
4. `Rscript scripts/pathways/run_pathways.R`. KEGG serves only its current release, so gene sets differ from those used for the paper.
5. Output: `psPGRS_<date>_<PRS>_<gene sets>_<sources>_annotated_PS.PGRS.rds`. It holds one row per PRS variant: weight, alleles, and one TRUE/FALSE column per pathway, plus `all_SNPs`.

## §5 Score genotypes

Scoring runs as one job per chunk (at most `PSPRS_CHUNK_SIZE` variants, default 100000, never spanning chromosomes). Each chunk is read from the `.bgen` file in blocks of `PSPRS_BLOCKSIZE` variants (default 200).

1. Settings:
   - `PSPRS_DATA_DIR` and `PSPRS_MEMBERSHIP_FILE` (the §4 output).
   - `PSPRS_GENOTYPE_DIR` with the `.bgen` and `.bgi` files.
   - `PSPRS_GENOTYPE_FILE_PATTERN` (default `ukb_imp_{chr}_v3.bgen`; `{chr}` becomes `chr1`..`chr22`).
   - `PSPRS_CHR_LABEL`: the chromosome names inside the files. `zero_pad` gives `"01"` (UK Biobank, the default); `number` gives `"1"`.
   - `PSPRS_EFFECT_ALLELE_COL` (default `effect_allele`).
   - `PSPRS_OUT_DIR` and `PSPRS_RUN_ID`.
2. Find the number of chunks. The first task prints it, and writes `<PSPRS_OUT_DIR>/<PSPRS_RUN_ID>/run_manifest.rds` with `total_chunks`.
3. Run one task per chunk:
   - SLURM: `sbatch --array=1-<total_chunks> --wrap 'Rscript scripts/scoring/score_chunk.R'` (the task id comes from `SLURM_ARRAY_TASK_ID`).
   - SGE: the task id comes from `SGE_TASK_ID`.
   - Without a scheduler: `for k in $(seq 1 <total_chunks>); do PSPRS_TASK_ID=$k Rscript scripts/scoring/score_chunk.R; done`.
4. Resources: memory scales with block size × sample size, and with the number of pathway columns. Scoring every KEGG pathway needs more memory per task than the paper's 43 curated pathways.
5. Recovery: each completed block is saved under `<out>/<run>/blocks/<chunk>/`. A chunk with a failed block stops and lists the missing blocks; rerunning the same task scores only those. A run directory refuses a different input, so use a new `PSPRS_RUN_ID` for a new membership file or new settings.

## §6 Combine chunks

`Rscript scripts/scoring/combine_chunks.R` (same `PSPRS_OUT_DIR` and `PSPRS_RUN_ID`). It stops if any chunk summary is missing. Output: `PathSpecPGRS_total_scores.rds`, a samples × pathways matrix of PS-PRS.

## §7 Analysis record (UK Biobank)

These scripts need UK Biobank data and run only under an approved application.

1. `scripts/analysis/prepare_analysis_data.R`:
   - derives the covariates and sleep apnea exposures;
   - builds the curated-module scores from the single-gene scores (`scripts/analysis/curated_modules.R`);
   - standardises the scores, and aligns their sign so that the genome-wide score is risk-increasing for CAD.
2. `scripts/analysis/run_gxe_regression.R`:
   - fits CAD ~ covariates + E + G + E × G by logistic regression (`glm`) and by a generalised additive model (`mgcv::bam`), using `R/model_fit_fns.R`, `scripts/analysis/model_covariates_primary.R` and `scripts/analysis/model_formulas.R`;
   - the primary analysis compares snoring sleep apnea with controls.
