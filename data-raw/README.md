# Data preparation

Run from the repository root, in order. Each script writes a provenance log to `data-raw/logs/` (input checksums, counts, session info).

| Script | Input | Output (`app/data/`) | Where it runs |
|---|---|---|---|
| `00_scatlas_readiness.R` | `scAtlas.rds` from Zenodo 10.5281/zenodo.14199536 | `sc_ref_base.rds` (filtered atlas) | workstation, several hours, high memory |
| `01_scatlas_summaries.R` | `sc_ref_base.rds` + subcluster label files from github.com/PDAC-MULTIOMIC/PDAC_Atlas | `sc_groups.rds`, `sc_expr.rds`, `sc_umap.rds` | workstation, ~1 h, ~50 GB RAM |
| `01b_sc_coexpression.R` | outputs of 01 | `sc_coexp.rds` (most co-expressed genes per gene and cell type across tumour samples) | any machine, ~15 min, ~4 GB RAM |
| `02_cptac.R` | downloads from LinkedOmics | `cptac.rds` | any machine, minutes |
| `02b_cptac_phospho.R` | downloads from LinkedOmics (site-level phosphoproteome) | `cptac_phospho.rds` | any machine, ~2 min |
| `03_recount3.R` | downloads via recount3 | `bulk_tcga_gtex.rds` | any machine, minutes |
| `03b_clinical.R` | `cptac.rds`, `bulk_tcga_gtex.rds`; downloads the CPTAC clinical table (LinkedOmics), the TCGA Pan-Cancer Clinical Data Resource (GDC) and the recount3 TCGA-PAAD sample metadata | `clinical.rds` (survival, stage, age, sex and Moffitt classical/basal-like score per tumour) | any machine, a minute; needs `readxl` |
| `04_gene_index.R` | outputs of 01–03 | `genes.rds` | any machine, a minute |
| `04b_gene_xrefs.R` | `genes.rds` | `gene_xrefs.rds` (Ensembl ID and gene type for the links to other resources) | any machine, a minute |
| `05_build_store.R` | outputs of 01–04, 01b, 02b and 03b | `cptac.sqlite`, `bulk.sqlite`, `sc.sqlite`, `phospho.sqlite`, `clinical.sqlite`, `coexp.sqlite` (the data store the app reads) | any machine, ~7 min, ~2 GB RAM |

The app itself needs `genes.rds`, `cptac.sqlite` and `bulk.sqlite`; `sc.sqlite`, `phospho.sqlite`, `clinical.sqlite`, `coexp.sqlite` and `gene_xrefs.rds` are optional (the corresponding views or direct links are left out without them). Each `.sqlite` file holds one compressed record per gene plus the small tables (samples, cell groups, atlas-map coordinates), and `05` checks every stored record against the source tables. `get_data.R` (repository root) downloads the outputs of 01–04 from Zenodo and runs `05` (the outputs of 01b, 02b, 03b and 04b are added to Zenodo with the v0.2 data release).

`helpers.R` contains the logging helpers. Tests: `testthat::test_dir("tests/testthat")` runs `01` on a miniature atlas and checks the app's statistics.
