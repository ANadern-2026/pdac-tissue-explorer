# 06a_fetch_models.R: download the model-system source files for 06_models.R and write an inventory.
# Run from the repository root:  source("data-raw/06a_fetch_models.R")
# Downloads (skipped when already present with the right checksum):
#   1. DepMap (Broad), newest "DepMap Public" release (DepMap's site now asks browsers to pass a verification check, so the
#      file index must be saved by hand first: open https://depmap.org/portal/api/download/files in the browser and save it as
#      data-raw/source/models/depmap/depmap_file_index.csv): Model.csv, CRISPRGeneEffect.csv, the protein-coding TPM expression
#      matrix and the model/profile mapping files; plus the newest harmonised public proteomics files (Broad data only, CC BY 4.0).
#   2. GDC open-access RNA-seq gene expression (STAR counts, same pipeline as TCGA) for
#      ORGANOID-PANCREATIC (Tiriac et al., Cancer Discov 2018) and HCMI-CMDC models of the pancreas.
# Writes data-raw/source/models/{depmap,gdc}/ and the log data-raw/logs/06a_fetch_models_<date>.txt (inventory of files, checksums, headers).
# Expect ~1.5-2 GB in total and 10-40 min depending on the connection.
source(file.path("data-raw", "helpers.R"))
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install the jsonlite package: install.packages(\"jsonlite\")")
options(timeout = max(3600, getOption("timeout")))
log_file <- start_log("06a_fetch_models")
SRC <- file.path("data-raw", "source", "models")
D_DM <- file.path(SRC, "depmap"); D_GDC <- file.path(SRC, "gdc")
dir.create(D_DM, recursive = TRUE, showWarnings = FALSE); dir.create(D_GDC, recursive = TRUE, showWarnings = FALSE)

is_html <- function(f) { x <- tryCatch(readChar(f, 200L, useBytes = TRUE), error = function(e) ""); grepl("^\\s*(<!DOCTYPE html|<html)", x, ignore.case = TRUE) }
fetch <- function(url, dest, md5 = NA_character_) {
  if (file.exists(dest) && file.size(dest) > 0 && (is.na(md5) || !nzchar(md5) || identical(unname(tools::md5sum(dest)), md5)))
    return(invisible(TRUE))
  tmp <- paste0(dest, ".part"); unlink(tmp)
  ok <- tryCatch({ utils::download.file(url, tmp, mode = "wb", quiet = TRUE); file.size(tmp) > 0 }, error = function(e) { logc(log_file, "  error: ", conditionMessage(e)); FALSE },
                 warning = function(w) { logc(log_file, "  warning: ", conditionMessage(w)); FALSE })
  if (!isTRUE(ok)) { unlink(tmp); logc(log_file, "FAILED ", basename(dest)); return(invisible(FALSE)) }
  if (is_html(tmp)) { logc(log_file, "BLOCKED ", basename(dest), " (the site returned a web page, e.g. a verification check, instead of the file)"); unlink(tmp); return(invisible(FALSE)) }
  got <- unname(tools::md5sum(tmp))
  if (!is.na(md5) && nzchar(md5) && !identical(got, md5)) { logc(log_file, "MD5 MISMATCH ", basename(dest), " (expected ", md5, ", got ", got, ")"); unlink(tmp); return(invisible(FALSE)) }
  file.rename(tmp, dest); invisible(TRUE)
}
peek_csv <- function(f, n_fields = 6) {        # header fields, first row start, number of lines (without reading the whole file into R)
  con <- file(f, "r"); on.exit(close(con))
  h <- readLines(con, n = 2L, warn = FALSE)
  n <- 2L; repeat { x <- readLines(con, n = 20000L, warn = FALSE); if (!length(x)) break; n <- n + length(x) }
  hf <- strsplit(h[1], ",", fixed = TRUE)[[1]]; r1 <- strsplit(h[2], ",", fixed = TRUE)[[1]]
  sprintf("%d lines, %d header fields; header: %s ... %s | row 1: %s", n, length(hf), paste(utils::head(hf, n_fields), collapse = " | "),
          paste(utils::tail(hf, 2), collapse = " | "), paste(utils::head(r1, n_fields), collapse = " | "))
}

## ---- 1. DepMap -------------------------------------------------------------------------------------------
logc(log_file, "== DepMap")
manual <- character(0)                        # DepMap files to download by hand (if the site blocks scripts)
depmap_part <- function() {
idx_file <- file.path(D_DM, "depmap_file_index.csv")
unlink(list.files(D_DM, pattern = "^depmap_file_index_.*\\.csv$", full.names = TRUE))   # verification pages saved by the first version of this script
if (!file.exists(idx_file) || is_html(idx_file))
  stop("Save the DepMap file index first: open https://depmap.org/portal/api/download/files in your browser and save it as ", idx_file)
logc(log_file, "file index: ", idx_file, " (", format(file.mtime(idx_file), "%Y-%m-%d %H:%M"), "), MD5 ", unname(tools::md5sum(idx_file)))
idx <- utils::read.csv(idx_file, stringsAsFactors = FALSE, check.names = FALSE)
logc(log_file, "index columns: ", paste(names(idx), collapse = ", "), "; ", nrow(idx), " files")
col <- function(...) { for (pat in c(...)) { k <- grep(pat, names(idx), ignore.case = TRUE, value = TRUE); if (length(k)) return(k[1]) }; stop("no index column like ", paste(c(...), collapse = " / ")) }
c_rel <- col("^release$", "^release.?name$", "release"); c_file <- col("^file.?name$", "file.?name"); c_url <- col("^url$", "url", "download")
c_md5 <- grep("md5", names(idx), ignore.case = TRUE, value = TRUE)[1]; c_date <- grep("date", names(idx), ignore.case = TRUE, value = TRUE)[1]
rel_order <- function(r) { d <- if (!is.na(c_date)) as.Date(idx[[c_date]][match(r, idx[[c_rel]])]) else as.Date(NA)
  if (all(is.na(d))) order(r, decreasing = TRUE) else order(d, decreasing = TRUE) }
rels <- unique(idx[[c_rel]])
logc(log_file, "releases (newest first): ", paste(utils::head(rels[rel_order(rels)], 12), collapse = "; "))
public <- rels[grepl("^DepMap Public", rels)]; public <- public[rel_order(public)]
if (!length(public)) stop("no 'DepMap Public' release in the index")
REL <- public[1]; logc(log_file, "using release: ", REL)
in_rel <- idx[idx[[c_rel]] == REL, , drop = FALSE]
logc(log_file, "files in ", REL, ": ", paste(sort(in_rel[[c_file]]), collapse = ", "))
want <- c("^Model\\.csv$", "^CRISPRGeneEffect\\.csv$", "^OmicsProfiles\\.csv$", "^OmicsDefaultModelProfiles\\.csv$", "^ModelCondition\\.csv$",
          "^OmicsExpression.*(ProteinCoding|HumanProteinCodingGenes).*TPM.*\\.csv$|^OmicsExpression.*TPM.*(ProteinCoding|HumanProteinCodingGenes).*\\.csv$")
pick <- unique(unlist(lapply(want, function(p) grep(p, in_rel[[c_file]], value = TRUE))))
pick <- pick[!grepl("BatchCorrected|Stranded|Profile.*Expression|Transcript", pick) | pick %in% c("OmicsProfiles.csv", "OmicsDefaultModelProfiles.csv")]
# proteomics: newest release whose name mentions proteomics (Broad harmonised files only; Sanger data are not redistributable)
prot_rels <- rels[grepl("proteom", rels, ignore.case = TRUE) & !grepl("sanger", rels, ignore.case = TRUE)]; prot_rels <- prot_rels[rel_order(prot_rels)]
prot <- if (length(prot_rels)) idx[idx[[c_rel]] == prot_rels[1], , drop = FALSE] else idx[0, , drop = FALSE]
if (length(prot_rels)) logc(log_file, "proteomics release: ", prot_rels[1], " - files: ", paste(prot[[c_file]], collapse = ", "))
dl <- rbind(in_rel[in_rel[[c_file]] %in% pick, , drop = FALSE], prot[!grepl("sanger", prot[[c_file]], ignore.case = TRUE), , drop = FALSE])
logc(log_file, "download hosts: ", paste(unique(sub("^(https?://[^/]+).*$", "\\1", dl[[c_url]])), collapse = ", "))
for (i in seq_len(nrow(dl))) {
  f <- dl[[c_file]][i]; sub <- if (dl[[c_rel]][i] == REL) "" else "proteomics"
  dd <- if (nzchar(sub)) file.path(D_DM, sub) else D_DM; dir.create(dd, showWarnings = FALSE)
  dest <- file.path(dd, f); md5 <- if (!is.na(c_md5)) dl[[c_md5]][i] else NA_character_
  if (file.exists(dest) && !is_html(dest) && (is.na(md5) || !nzchar(md5) || identical(unname(tools::md5sum(dest)), md5))) { logc(log_file, "present ", if (nzchar(sub)) file.path(sub, f) else f); next }
  logc(log_file, "downloading ", if (nzchar(sub)) file.path(sub, f) else f, " ...")
  if (!fetch(dl[[c_url]][i], dest, md5)) manual <<- c(manual, sprintf("%s  (release '%s') -> save as %s", f, dl[[c_rel]][i], dest))
}
writeLines(c(REL, if (length(prot_rels)) prot_rels[1]), file.path(D_DM, "RELEASE.txt"))
}
tryCatch(depmap_part(), error = function(e) { logc(log_file, "DepMap part stopped: ", conditionMessage(e)); manual <<- c(manual, conditionMessage(e)) })

## ---- 2. GDC organoids --------------------------------------------------------------------------------------
logc(log_file, "== GDC")
gdc_query <- function(filters) {
  fields <- c("file_id", "file_name", "md5sum", "file_size", "access", "data_type", "analysis.workflow_type", "cases.case_id", "cases.submitter_id",
              "cases.project.project_id", "cases.primary_site", "cases.disease_type", "cases.diagnoses.primary_diagnosis",
              "cases.samples.submitter_id", "cases.samples.sample_type", "cases.samples.tissue_type", "cases.samples.tumor_descriptor",
              "cases.samples.specimen_type", "cases.samples.preservation_method")
  url <- paste0("https://api.gdc.cancer.gov/files?size=5000&format=JSON&fields=", paste(fields, collapse = ","),
                "&filters=", utils::URLencode(jsonlite::toJSON(filters, auto_unbox = TRUE), reserved = TRUE))
  tmp <- tempfile(fileext = ".json"); utils::download.file(url, tmp, mode = "wb", quiet = TRUE)
  h <- jsonlite::fromJSON(tmp, simplifyVector = FALSE)$data$hits
  g <- function(x, ...) { for (k in list(...)) { if (is.null(x)) return(NA_character_); x <- if (is.numeric(k)) (if (length(x) >= k) x[[k]] else NULL) else x[[k]] }
                          if (is.null(x)) NA_character_ else paste(unlist(x), collapse = ";") }
  do.call(rbind, lapply(h, function(x) data.frame(
    file_id = x$file_id, file_name = x$file_name, md5sum = x$md5sum, file_size = as.numeric(x$file_size), access = x$access, workflow = g(x, "analysis", "workflow_type"),
    project = g(x, "cases", 1, "project", "project_id"), case = g(x, "cases", 1, "submitter_id"), primary_site = g(x, "cases", 1, "primary_site"),
    disease_type = g(x, "cases", 1, "disease_type"), diagnosis = g(x, "cases", 1, "diagnoses", 1, "primary_diagnosis"),
    n_samples = length(x$cases[[1]]$samples), sample = g(x, "cases", 1, "samples", 1, "submitter_id"), sample_type = g(x, "cases", 1, "samples", 1, "sample_type"),
    tissue_type = g(x, "cases", 1, "samples", 1, "tissue_type"), tumor_descriptor = g(x, "cases", 1, "samples", 1, "tumor_descriptor"),
    specimen_type = g(x, "cases", 1, "samples", 1, "specimen_type"), preservation = g(x, "cases", 1, "samples", 1, "preservation_method"), stringsAsFactors = FALSE)))
}
f_and <- function(...) list(op = "and", content = list(...)); f_in <- function(field, value) list(op = "in", content = list(field = field, value = as.list(value)))
expr <- list(f_in("data_type", "Gene Expression Quantification"), f_in("access", "open"))
q <- list(
  ORGANOID = do.call(f_and, c(list(f_in("cases.project.project_id", "ORGANOID-PANCREATIC")), expr)),
  HCMI = do.call(f_and, c(list(f_in("cases.project.project_id", "HCMI-CMDC"), f_in("cases.primary_site", c("Pancreas", "pancreas"))), expr)))
meta <- do.call(rbind, lapply(names(q), function(k) { m <- gdc_query(q[[k]]); if (is.null(m)) { logc(log_file, k, ": no files found"); return(NULL) }; m$set <- k; m }))
if (is.null(meta) || !nrow(meta)) stop("GDC returned no expression files")
utils::write.table(meta, file.path(D_GDC, "gdc_expression_files.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
for (k in unique(meta$set)) { m <- meta[meta$set == k, ]
  logc(log_file, k, ": ", nrow(m), " files, ", length(unique(m$case)), " cases, ", round(sum(m$file_size) / 1e6), " MB; workflows: ", paste(unique(m$workflow), collapse = ", "))
  tb <- table(paste(m$sample_type, m$tissue_type, m$tumor_descriptor, m$specimen_type, sep = " / ")); for (n in names(tb)) logc(log_file, "   ", tb[[n]], " x ", n)
  logc(log_file, "   diagnoses: ", paste(names(sort(table(m$diagnosis), decreasing = TRUE)), collapse = "; "))
  logc(log_file, "   files with >1 sample: ", sum(m$n_samples > 1)) }
for (i in seq_len(nrow(meta))) {
  dd <- file.path(D_GDC, meta$set[i]); dir.create(dd, showWarnings = FALSE)
  for (k in 1:4) {                            # GDC drops connections now and then ("Failure when receiving data from the peer"): retry
    if (fetch(paste0("https://api.gdc.cancer.gov/data/", meta$file_id[i]), file.path(dd, meta$file_name[i]), meta$md5sum[i])) break
    Sys.sleep(5 * k)
  }
  if (i %% 20 == 0) message(i, " of ", nrow(meta), " GDC files")
}

## ---- 3. inventory ---------------------------------------------------------------------------------------------
logc(log_file, "== inventory")
fs <- list.files(SRC, recursive = TRUE, full.names = TRUE); fs <- fs[!grepl("\\.part$", fs)]
inv <- data.frame(file = sub(paste0("^", SRC, "/"), "", fs), bytes = file.size(fs), md5 = unname(tools::md5sum(fs)))
for (i in which(!grepl("^gdc/(ORGANOID|HCMI)/", inv$file))) logc(log_file, inv$file[i], "  ", inv$bytes[i], " bytes  MD5 ", inv$md5[i])
g <- inv[grepl("^gdc/(ORGANOID|HCMI)/", inv$file), ]
logc(log_file, "GDC expression files on disk: ", nrow(g), " (", round(sum(g$bytes) / 1e6), " MB); missing: ", sum(!meta$file_name %in% basename(g$file)))
for (f in fs[grepl("\\.csv$", fs) & grepl("/depmap/", fs)]) logc(log_file, "peek ", basename(f), ": ", tryCatch(peek_csv(f), error = function(e) conditionMessage(e)))
g1 <- fs[grepl("^.*/gdc/(ORGANOID|HCMI)/", fs)][1]
if (!is.na(g1)) { logc(log_file, "first lines of ", basename(g1), ":"); for (l in readLines(g1, n = 8L)) logc(log_file, "   ", substr(l, 1, 200)) }
if (length(manual)) {
  logc(log_file, "== to do by hand (DepMap blocked the script): download these on https://depmap.org/portal/data_page/?tab=allData, then run this script again")
  for (m in manual) logc(log_file, "   ", m)
}
end_log(log_file)
message("Done. Log: ", log_file, if (length(manual)) "  -- some DepMap files need a manual download, see the end of the log" else "")
