# Pathway-specific polygenic risk scores (PS-PRS) for coronary artery disease

Code accompanying:

> Goodman MO, Cade BE, Shah NA, Huang T, Dashti HS, Saxena R, Rutter MK, Libby P, Sofer T, Redline S.
> **Pathway-Specific Polygenic Risk Scores Identify Obstructive Sleep Apnea-Related Pathways
> Differentially Moderating Genetic Susceptibility to Coronary Artery Disease.**
> *Circ Genom Precis Med* 2022;15(5):e003535. doi:[10.1161/CIRCGEN.121.003535](https://doi.org/10.1161/CIRCGEN.121.003535)

A pathway-specific PRS sums a genome-wide PRS over only the variants assigned to the genes of one
pathway. The paper split the CAD PRS of Khera et al. 2018 (PGS Catalog
[PGS000013](https://www.pgscatalog.org/score/PGS000013/)) into KEGG and curated pathways, scored
UK Biobank participants, and tested pathway × obstructive sleep apnea interactions on CAD.

This repository contains:

- **a reusable scorer** for imputed genotypes in `.bgen` files (UK Biobank imputed v3 layout by
  default), with a runnable example on public synthetic genotypes;
- **the builders** that assign PRS variants to genes (GENCODE, GTEx, FANTOM5, promoter capture
  Hi-C, GeneHancer) and mark each variant's pathways (KEGG, curated);
- **the analysis record**: the scripts that prepared the UK Biobank data and fitted the primary
  GxE models. They need UK Biobank data and are included for transparency.

## Repository layout

```
setup/               build rbgen and the project R library (renv)
config/paths.R       every input and output location, from environment variables
R/                   functions: scorer, dosage, PRS loading, gene annotation, KEGG, model fitting
scripts/annotation/  variant -> gene annotation by source, and the source-priority cascade
scripts/pathways/    pathway gene lists and per-variant pathway membership
scripts/scoring/     score one chunk of the PRS (array job), then combine the chunks
scripts/analysis/    analysis record: data preparation and the primary GxE models
data/                curated pathway sheet used for the paper
example/             toy PRS for the bgen example genotypes, and the example run
tests/               scorer, annotation and membership tests
RUNBOOK.md           step-by-step procedures
```

## Quick start (example data)

Requires R (tested with 4.5) with a C/C++ compiler, python 3, curl and the `renv` package.

```sh
bash setup/build_rbgen.sh            # fetch bgen v1.1.7, build rbgen, unpack the example genotypes
Rscript --vanilla setup/setup_renv.R # project library (a conda R: run both through `conda run -p <env>`)
Rscript example/run_example.R        # score the toy PRS on the example genotypes
```

The example scores 500 synthetic samples on four toy pathways with the same scripts used for a
full run. [RUNBOOK.md](RUNBOOK.md) covers scoring your own genotypes, building the annotation and
pathway membership, and the analysis record.

## Effect allele

The scorer multiplies each weight by the dosage of the allele named in the PRS column
`effect_allele_col` (default `effect_allele`), so weights coded on either allele score correctly.

The scores used for the paper counted the other PRS allele (A2), and so were sign-reversed. The
analysis step aligned the sign before modelling (`scripts/analysis/prepare_analysis_data.R`
flips the scores when the genome-wide score is not risk-increasing for CAD), and standardised
them, so the published results are unaffected. With `effect_allele_col = "A2"` the scorer
reproduces the paper's raw scores exactly (`tests/test_scorer_equivalence.R`).

## Data sources and licences

No individual-level data are included. The annotation uses public datasets with their own terms.
The Zenodo deposit (next section) holds copies of the files marked "yes", as downloaded for the
paper, so the annotation can be rebuilt if an original host is unavailable.

| Source | Used for | Original download | On Zenodo | Licence / terms |
|---|---|---|---|---|
| PGS Catalog PGS000013 (Khera et al. 2018) | PRS weights | [PGS000013](https://www.pgscatalog.org/score/PGS000013/) | no | Academic research use; see the PGS Catalog entry |
| GENCODE gencode.v34lift37 (GRCh37) | genes, exons, UTRs, flanks | [EBI FTP](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.annotation.gff3.gz) | yes | Open access |
| GTEx v7 single-tissue eQTL (artery, heart, whole blood) | eQTL gene links | [GTEx Portal](https://storage.googleapis.com/adult-gtex/bulk-qtl/v7/single-tissue-cis-qtl/GTEx_Analysis_v7_eQTL.tar.gz) | yes (6 tissues used) | GTEx Portal terms (below) |
| FANTOM5 enhancer–TSS associations (Andersson et al. 2014) | enhancer gene links | `http://enhancer.binf.ku.dk/presets/enhancer_tss_associations.bed` (offline as of November 2026) | yes | CC BY 4.0 |
| Promoter capture Hi-C (Javierre et al. 2016) | enhancer gene links | `ActivePromoterEnhancerLinks.tsv` in Data S1 of [the article](https://doi.org/10.1016/j.cell.2016.09.037) | yes | CC BY 4.0 |
| GeneHancer (Fishilevich et al. 2017), version 4.4 | enhancer gene links | LifeMap Sciences (GeneCards) | no | Licensed by LifeMap Sciences; this repository does not include it. To annotate with GeneHancer, as the paper did, obtain the GeneHancer csv from LifeMap Sciences and set `PSPRS_GENEHANCER_CSV` to its path. If you leave it unset, the builder skips GeneHancer and assigns variants from the other sources |
| KEGG pathways (KEGGREST) | pathway gene lists | KEGG REST API, fetched by the builder | no | Academic use. The builder fetches the current KEGG release, which differs from the 2020-12-09 lists the paper used (gene sets, and some pathway and gene names, which also name the score columns); those lists are available from the authors on request |
| UK Biobank | genotypes and phenotypes | UK Biobank | no | Approved application only |

**GTEx.** The data used for the analyses described in the paper were obtained from the GTEx
Portal (release v7; local copy dated 04/10/2019). Data redistributed with this project may not
reflect the most current or accurate data available from the GTEx Portal, and the GTEx Project
does not endorse this project.

## Annotation data on Zenodo

Zenodo: [doi:10.5281/zenodo.23990992](https://doi.org/10.5281/zenodo.23990992) (version 1, the files this code expects; all versions: [10.5281/zenodo.23990991](https://doi.org/10.5281/zenodo.23990991)). The deposit holds:

- `PSPRS_annotation_GRCh37_20201209.rds`: the paper's per-source variant → gene annotation of
  the PGS000013 variants, for every source except GeneHancer and without the PRS weights. It is
  the `.rds` list of GRanges that `scripts/annotation/assign_genes_cascade.R` reads
  (`scripts/annotation/export_release_annotation.R` produces it).
- Copies of the GENCODE, GTEx (six tissues), FANTOM5 and promoter capture Hi-C source files
  marked "yes" above, with checksums and a README giving each file's source, licence and
  citation.

RUNBOOK §3 shows how to use either the original downloads or the Zenodo copies.

## Tests

```sh
Rscript -e 'testthat::test_dir("tests", stop_on_failure = FALSE)'
```

The example-data tests always run. Tests against the original code, the paper's 2020 files, or
the network are skipped unless their environment variables are set (see each test file).

## Licence

Code: MIT ([LICENSE](LICENSE)), copyright Mass General Brigham. Data files keep their sources' licences (above). The bgen library
and rbgen are distributed by their authors under the Boost Software License; this repository
downloads them and does not redistribute them.
