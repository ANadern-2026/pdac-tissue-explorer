# ---- Download the prepared data of the PDAC Tissue Explorer from Zenodo --------------------------------------
# Dataset: Nadernezhad A. PDAC Tissue Explorer - processed data (v0.3.0). Zenodo. doi:10.5281/zenodo.23977456
# Run from the repository root:  source("get_data.R")  (or: Rscript get_data.R)   then   shiny::runApp("app")
# Downloads the app's data store (~640 MB: genes.rds, gene_xrefs.rds and one .sqlite file per data layer) into app/data/
# and verifies each file's MD5 checksum. Files already present and intact are skipped. A file that is present but differs
# from the published version (e.g. rebuilt locally with data-raw/) is kept, with a message; set
# options(pdacx.get_data.overwrite = TRUE) to replace it.
# Needs the R packages DBI and RSQLite.
#
# Rebuild from tables instead (the v0.1 route): options(pdacx.get_data.from_tables = TRUE) downloads the six v0.1 tables
# (~245 MB; carried over unchanged into this deposit) and converts them into cptac.sqlite, bulk.sqlite and sc.sqlite with
# data-raw/05_build_store.R (5-15 minutes, ~2 GB RAM, also needs Matrix). The other layers are then downloaded as above.

DATA_RECORD <- "23977456"                   # Zenodo record of data v0.3.0
STORE_FILES <- c(                                # file name = MD5 checksum (as published on Zenodo)
  "genes.rds"       = "136b0b189bee803bec7730dcbafee554",
  "gene_xrefs.rds"  = "f69c7fbca9200f01afc22d5012aacf85",
  "cptac.sqlite"    = "ecfae858a474bd9447379b057fea176a",
  "bulk.sqlite"     = "6e6ed930fc80c88e6879699e7a2d4771",
  "sc.sqlite"       = "40db9da64c023ebbc26cfdfd78ee1e74",
  "phospho.sqlite"  = "23bec8c165bb4d1ac6d81d660aeb5cc8",
  "clinical.sqlite" = "b905eef43b487a98db08759ceac56d12",
  "coexp.sqlite"    = "add1577c37caf2d0999eb96dcf7f6509",
  "models.sqlite"   = "4141d400a0ce96b92a83db831dfce78b")
TABLE_FILES <- c(                                # v0.1 tables (doi:10.5281/zenodo.23264998), input of 05_build_store.R
  "cptac.rds"          = "2f46c794435d49dfe78005f20774bf4f",
  "bulk_tcga_gtex.rds" = "468c88a4c3d785e25ed476c7f98977ff",
  "sc_groups.rds"      = "38d9fc2793eb35d12ae125f04e586b9c",
  "sc_expr.rds"        = "8dfaba06afcdf84e0a3da5e3e3a82d6e",
  "sc_umap.rds"        = "3012ec7a4296887661c0b6ace1eef642")
FROM_TABLES <- c("cptac.sqlite", "bulk.sqlite", "sc.sqlite")   # store files that the table route builds locally
OUT_DIR <- file.path("app", "data")

if (!dir.exists("app") || !file.exists(file.path("app", "R", "01_data.R")))
  stop("Run this script from the repository root (the folder that contains app/ and data-raw/).")
if (grepl("@", DATA_RECORD, fixed = TRUE)) stop("DATA_RECORD is not set in get_data.R.")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
from_tables <- isTRUE(getOption("pdacx.get_data.from_tables", FALSE))
overwrite   <- isTRUE(getOption("pdacx.get_data.overwrite", FALSE))
old_timeout <- options(timeout = max(3600, getOption("timeout")))

md5_of <- function(path) unname(tools::md5sum(path))
fetch <- function(f, md5) {
  dest <- file.path(OUT_DIR, f)
  if (file.exists(dest)) {
    if (md5_of(dest) == md5) { message(f, ": already present"); return(invisible("present")) }
    if (!overwrite) {
      message(f, ": present but differs from the published file (local build?) - kept. ",
              "Use options(pdacx.get_data.overwrite = TRUE) to replace it."); return(invisible("kept")) }
  }
  url <- sprintf("https://zenodo.org/records/%s/files/%s?download=1", DATA_RECORD, f)
  message(f, ": downloading ...")
  tmp <- paste0(dest, ".part")
  for (try in 1:3) {
    ok <- tryCatch({ utils::download.file(url, tmp, mode = "wb", quiet = TRUE); file.exists(tmp) && md5_of(tmp) == md5 },
                   error = function(e) { message("  ", conditionMessage(e)); FALSE })
    if (ok) break
    unlink(tmp); if (try < 3) message("  retrying ...")
  }
  if (!ok) stop(f, ": download failed or checksum mismatch; please try again.")
  if (file.exists(dest)) unlink(dest)
  file.rename(tmp, dest)
  message(f, ": ok")
  invisible("downloaded")
}

if (from_tables) {
  for (f in names(TABLE_FILES)) fetch(f, TABLE_FILES[[f]])
  for (f in setdiff(names(STORE_FILES), FROM_TABLES)) fetch(f, STORE_FILES[[f]])
  message("Building cptac.sqlite, bulk.sqlite and sc.sqlite from the tables (skipped when up to date) ...")
  old_opt <- options(pdacx.build_store.autorun = FALSE)
  source(file.path("data-raw", "05_build_store.R"))
  options(old_opt)
  build_store(OUT_DIR, layers = c("cptac", "bulk", "sc"))
} else {
  for (f in names(STORE_FILES)) fetch(f, STORE_FILES[[f]])
}
options(old_timeout)

# every store file must open with this version of the app (format and layer checks of store_open())
suppressPackageStartupMessages({ library(DBI); library(RSQLite) })
for (f in c("00_config.R", "01_data.R")) source(file.path("app", "R", f))
for (f in grep("\\.sqlite$", names(STORE_FILES), value = TRUE)) {
  path <- file.path(OUT_DIR, f); if (!file.exists(path)) next
  st <- store_open(path); store_close(st)
}
message("All data files are in ", normalizePath(OUT_DIR), ". Start the app with shiny::runApp(\"app\").")
