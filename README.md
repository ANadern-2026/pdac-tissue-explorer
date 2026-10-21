# PDAC Tissue Explorer

[![Code DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23264879.svg)](https://doi.org/10.5281/zenodo.23264879) [![Data DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23264998.svg)](https://doi.org/10.5281/zenodo.23264998)

An interactive web app for exploring where any human gene is expressed in pancreatic ductal adenocarcinoma (PDAC) tissue: in tumour versus matched normal tissue, versus healthy pancreas, in which cell types, across tissue states (healthy donor, adjacent normal, primary tumour, liver metastasis, treated tumour) and in published fibroblast and tumour-cell subtypes. It combines three public datasets and needs no programming: type a gene, read the summary, download the figures and the numbers.

**Live app:** <https://pdac.biointeract.space/> · **Code:** [doi:10.5281/zenodo.23264879](https://doi.org/10.5281/zenodo.23264879) · **Data:** [doi:10.5281/zenodo.23264998](https://doi.org/10.5281/zenodo.23264998) · **Author:** Ali Nadernezhad, Leibniz-Institut für Polymerforschung Dresden ([ORCID 0000-0002-3761-1601](https://orcid.org/0000-0002-3761-1601))

## What you can do

- **Explore a gene.** A plain-language summary, then sample-level views: paired tumour vs adjacent-normal protein, RNA and phosphosites (CPTAC; phosphosites also relative to the protein level), PDAC vs healthy pancreas (TCGA/GTEx), expression by cell type with transcript shares, expression across tissue states per cell type, fibroblast and tumour-cell subtypes, the gene on the atlas map, its clinical context (survival, stage and the classical/basal-like tumour subtype in TCGA-PAAD and CPTAC) and the genes it is co-expressed with within each cell type across tumours. Search by symbol, alias or full name; links to the gene's pages at NCBI Gene, UniProt, the Human Protein Atlas, Open Targets and DepMap.
- **Compare genes.** Up to 40 genes at once as a cell-type dot heatmap and a tumour-vs-normal forest plot.
- **Download.** Every figure as PDF/PNG and the underlying per-sample values as CSV.

## Data sources

| Dataset | Content used | Reference |
|---|---|---|
| CPTAC PDAC (LinkedOmics) | Gene-level proteome and RNA, site-level phosphoproteome, clinical data (stage, overall survival); 135 PDAC tumours, 64 protein/phosphosite and 19 RNA tumour–normal pairs | Cao L. et al., *Cell* 2021;184:5031 |
| TCGA Pan-Cancer Clinical Data Resource | Overall survival and progression-free interval, stage, age for the TCGA-PAAD tumours | Liu J. et al., *Cell* 2018;173:400 |
| TCGA-PAAD + GTEx pancreas (recount3) | Bulk RNA: 147 PDAC (ductal histology), 4 adjacent normal, 330 healthy pancreas | Wilks C. et al., *Genome Biol* 2021;22:323; TCGA, *Cancer Cell* 2017; GTEx, *Science* 2020 |
| PDAC single-cell atlas | ~620,000 cells, 234 samples from 12 studies, with published cell-type, tumour-cell and fibroblast subclusters | Loveless I.M. et al., *Clin Cancer Res* 2025;31:756 (Zenodo 10.5281/zenodo.14199536) |

Please cite this resource together with the original data sources (see `CITATION.cff`):

> Nadernezhad A. PDAC Tissue Explorer: an interactive resource for gene expression in pancreatic ductal adenocarcinoma tissue. Zenodo. doi:10.5281/zenodo.23264879 (software); processed data: doi:10.5281/zenodo.23264998.

## Repository layout

```
app/                 Shiny app (app.R, R/ modules, www/ styles and About page)
app/data/            prepared data (not in Git; get_data.R or data-raw/): gene index + on-disk store (*.sqlite)
get_data.R           downloads the prepared data files from Zenodo (checksum-verified) and builds the data store
data-raw/            data preparation, run in order:
  00_scatlas_readiness.R   filtering of the published atlas object (cycling cells, ambient-RNA outliers)
  01_scatlas_summaries.R   per-sample x cell-type summaries, atlas map subsample
  01b_sc_coexpression.R    co-expression within each cell type across tumour samples (from the 01 summaries)
  02_cptac.R               CPTAC download and processing
  02b_cptac_phospho.R      CPTAC phosphosites (site level)
  03_recount3.R            TCGA-PAAD / GTEx via recount3
  03b_clinical.R           survival, stage and classical/basal-like (Moffitt) score of the CPTAC and TCGA tumours
  04_gene_index.R          searchable gene index, symbol reconciliation across datasets
  04b_gene_xrefs.R         Ensembl IDs and gene types for the links to other resources
  05_build_store.R         on-disk data store: one SQLite file per dataset, one compressed record per gene
  logs/                    provenance logs (input checksums, counts, session info)
tests/               testthat tests: pipeline on a miniature atlas, data store, app statistics;
                     regression-v01.R compares every output of the data store with the v0.1 in-memory data
deploy/connect-cloud/           legacy Connect Cloud address: relocation notice (deploy.R), later a forwarding page (retire.R)
deploy/docker/Dockerfile        minimal container for local use
deploy/server/                  self-hosting on a server: Docker Compose + Caddy (HTTPS), see its README
```

## Run it locally

Requires R ≥ 4.3 with `shiny`, `bslib`, `ggplot2`, `plotly`, `DT`, `Matrix`, `DBI`, `RSQLite`, `survival`, `markdown`, `patchwork`, `scales`.

```r
source("get_data.R")    # downloads the prepared data (~245 MB) from Zenodo into app/data/ and builds the data store
shiny::runApp("app")
```

The app reads each gene's values on demand from the data store (`app/data/*.sqlite`), so it runs in well under 1 GB of memory.

## Rebuild the data

```r
source("data-raw/02_cptac.R")        # downloads from LinkedOmics
source("data-raw/02b_cptac_phospho.R")  # phosphosites, from LinkedOmics
source("data-raw/03_recount3.R")     # needs Bioconductor recount3 on the first run
source("data-raw/03b_clinical.R")    # clinical data: LinkedOmics, TCGA-CDR (GDC) and recount3 sample metadata; needs readxl
# single-cell: download scAtlas.rds from Zenodo, then
source("data-raw/00_scatlas_readiness.R")   # ~ several hours; produces sc_ref_base.rds
source("data-raw/01_scatlas_summaries.R")   # ~ 1 h, ~50 GB RAM; set the paths at the top
source("data-raw/04_gene_index.R")
source("data-raw/04b_gene_xrefs.R")         # identifiers for the links to other resources
source("data-raw/01b_sc_coexpression.R")    # ~15 min; co-expression partners per cell type
source("data-raw/05_build_store.R")         # converts the tables into the app's on-disk data store
testthat::test_dir("tests/testthat")
```

## Methods, limitations and interpretation

See the app's **About** page (`app/www/about.md`). In short: tumour–normal comparisons are paired within patients (Wilcoxon signed-rank, bootstrap CI); single-cell summaries are per sample (≥ 20 cells of a type) with the library as replication unit; TCGA vs GTEx differences are descriptive (no cross-study batch correction); survival associations are exploratory (median-split Kaplan–Meier and continuous Cox models, no optimised cut-points); co-expression is a study-adjusted Spearman correlation across tumour samples within one cell type; p-values are unadjusted (phosphosites: false discovery rate across the sites of the selected gene). Abundance does not show activity or function, adjacent "normal" pancreas is not healthy tissue, and in the atlas tissue state is partly confounded with study.

## Licence

Code: MIT (`LICENSE`). Derived data tables: CC BY 4.0, consistent with the source datasets' terms; please cite the sources.

## Acknowledgements and AI assistance

This resource builds entirely on the generosity of the CPTAC, TCGA, GTEx and recount3 teams and the authors of the twelve single-cell studies and of the integrated atlas. Code was drafted with the assistance of an AI model (Anthropic Claude); all data processing was executed, checked and curated by the author.
