# Data preparation

Run from the repository root, in order. Each script writes a provenance log to `data-raw/logs/` (input checksums, counts, session info).

| Script | Input | Output (`app/data/`) | Where it runs |
|---|---|---|---|
| `00_scatlas_readiness.R` | `scAtlas.rds` from Zenodo 10.5281/zenodo.14199536 | `sc_ref_base.rds` (filtered atlas) | workstation, several hours, high memory |
| `01_scatlas_summaries.R` | `sc_ref_base.rds` + subcluster label files from github.com/PDAC-MULTIOMIC/PDAC_Atlas | `sc_groups.rds`, `sc_expr.rds`, `sc_umap.rds` | workstation, ~1 h, ~50 GB RAM |
| `02_cptac.R` | downloads from LinkedOmics | `cptac.rds` | any machine, minutes |
| `03_recount3.R` | downloads via recount3 | `bulk_tcga_gtex.rds` | any machine, minutes |
| `04_gene_index.R` | outputs of 01–03 | `genes.rds` | any machine, a minute |

`helpers.R` contains the logging helpers. Tests: `testthat::test_dir("tests/testthat")` runs `01` on a miniature atlas and checks the app's statistics.
