# ---- Download the prepared data of the PDAC Tissue Explorer from Zenodo and build the app's data store --------
# Dataset: Nadernezhad A. PDAC Tissue Explorer - processed data (v0.1.0). Zenodo. doi:10.5281/zenodo.23264998
# Run from the repository root:  source("get_data.R")  (or: Rscript get_data.R)   then   shiny::runApp("app")
#  1. downloads the six prepared tables (~245 MB) into app/data/ and verifies each file's MD5 checksum
#     (files already present and intact are skipped);
#  2. converts them into the on-disk store the app reads (cptac.sqlite, bulk.sqlite, sc.sqlite; data-raw/05_build_store.R),
#     checking every gene against the downloaded tables. Takes 5-15 minutes and ~2 GB RAM; skipped when up to date.
# Needs the R packages Matrix, DBI and RSQLite.

ZENODO_RECORD <- "23264998"
FILES <- c(                                    # file name = MD5 checksum (as published on Zenodo)
  "genes.rds"          = "136b0b189bee803bec7730dcbafee554",
  "cptac.rds"          = "2f46c794435d49dfe78005f20774bf4f",
  "bulk_tcga_gtex.rds" = "468c88a4c3d785e25ed476c7f98977ff",
  "sc_groups.rds"      = "38d9fc2793eb35d12ae125f04e586b9c",
  "sc_expr.rds"        = "8dfaba06afcdf84e0a3da5e3e3a82d6e",
  "sc_umap.rds"        = "3012ec7a4296887661c0b6ace1eef642")
OUT_DIR <- file.path("app", "data")

if (!dir.exists("app") || !file.exists(file.path("data-raw", "05_build_store.R")))
  stop("Run this script from the repository root (the folder that contains app/ and data-raw/).")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
old_timeout <- options(timeout = max(1200, getOption("timeout")))

intact <- function(path, md5) file.exists(path) && unname(tools::md5sum(path)) == md5
for (f in names(FILES)) {
  dest <- file.path(OUT_DIR, f)
  if (intact(dest, FILES[[f]])) { message(f, ": already present"); next }
  url <- sprintf("https://zenodo.org/records/%s/files/%s?download=1", ZENODO_RECORD, f)
  message(f, ": downloading ...")
  tmp <- paste0(dest, ".part")
  utils::download.file(url, tmp, mode = "wb", quiet = TRUE)
  if (!intact(tmp, FILES[[f]])) { unlink(tmp); stop(f, ": checksum mismatch after download; please try again.") }
  file.rename(tmp, dest)
  message(f, ": ok")
}
options(old_timeout)

message("Building the data store (skipped for layers that are already up to date) ...")
old_opt <- options(pdacx.build_store.autorun = FALSE)
source(file.path("data-raw", "05_build_store.R"))
options(old_opt)
build_store(OUT_DIR)
message("All data files are in ", normalizePath(OUT_DIR), ". Start the app with shiny::runApp(\"app\").")
