# Perturb-ME

Analysis code for *"Perturb-ME: Scalable mechanism discovery from phenotype-enriched
genome-wide screens"* (Wang, Gu, Frangieh et al.).

Perturb-ME combines genome-scale CRISPR screening, marker-based phenotypic enrichment by FACS,
and multimodal single-cell profiling by Perturb-CITE-seq. Applied to surface MHC-I (HLA-A/B/C)
in the A375 melanoma line, it profiles HLA-low and HLA-high cells with matched RNA,
surface-protein and sgRNA measurements, fits a regularised linear model of perturbation
effects, and decomposes that model into co-functional regulatory modules and gene programs.

## Data

This repository contains **code only**. The data live in two places:

- **Raw and processed sequencing data** — NCBI Gene Expression Omnibus, accession in the
  paper's Data Availability statement. Single-cell gene expression, CITE-seq, cell hashing and
  CROP-seq guide-capture libraries from 20 10x channels, plus bulk sgRNA enrichment
  sequencing from genomic DNA.
- **Analysis-ready objects** — figshare,
  [doi.org/10.6084/m9.figshare.33264225](https://doi.org/10.6084/m9.figshare.33264225): the
  two AnnData objects (353,091 cells and the 277,152 sgRNA-singlet subset), the EM-refined
  ElasticNet coefficient matrix, the target-module and gene-program partitions, MAGeCK
  outputs, and the Co-Scientist prompt and response.

## Layout

```
src/            67 analysis modules, organised by stage
  RNA/          cell and gene filtering, QC
  Hash/         hashtag demultiplexing (Cumulus)
  CITE/         antibody-derived tag normalisation
  CROP/         guide calling and per-target coverage
    linear_model/   lm_target_EM.py         the published ElasticNet + EM model
                    cluster_by_kmeans*.py   target modules and gene programs
  gDNA/         bulk sgRNA enrichment, MAGeCK inputs
  clustering/   condition and expression clustering, cell-cycle scoring
  manuscript_figures/   one module per published panel, plus source-data builders
scripts/        batch drivers and the data-deposit builders
notebook/       5 notebooks that produced published panels
```

Figure panels map to `src/manuscript_figures/` by name — `1C_cluster_conditions.py`,
`2E_cluster_by_kmeans_genome_wide_features.py`, and so on. The Fig. 2b,c module and
gene-program partition is produced by `notebook/20260519_Figure2B_Rerun_BetaCorr.ipynb`.

## Environment

```bash
bash setup.sh          # creates the `perturbme` conda environment (Python 3.10)
conda activate perturbme
```

Core dependencies: scanpy, anndata, scikit-learn, scipy, pandas, harmonypy, scvi-tools.
MAGeCK v0.5.9.2 is used for the bulk screen analysis and is invoked separately.

## ⚠️ Paths are institution-specific

This code was written to run on a single institutional HPC cluster and **contains hardcoded
absolute paths** (`/gstore/...`, `/gnet/...`, and in the older linear-model modules
`/ahg/regevdata/...`). They appear as module-level constants, usually named `REPO`, `DATA` or
similar near the top of each file.

To run this elsewhere, point those constants at your own copy of the data from GEO/figshare.
Nothing else about the code is site-specific. We have deliberately left the paths as they were
when the analysis was run rather than rewriting them after the fact.

## Citation

Please cite the paper.

## License

MIT — see [LICENSE](LICENSE).
