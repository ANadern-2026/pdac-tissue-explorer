# =====================================================================================================
# 02_cptac.R
# CPTAC pancreatic ductal adenocarcinoma (Cao et al., Cell 2021) via LinkedOmics.
# Downloads gene-level proteome (median-normalised log2 intensity) and RNA (RSEM upper-quartile log2) for
# tumours and adjacent normal tissue plus the clinical table, keeps PDAC cases (adenosquamous excluded),
# and writes app/data/cptac.rds.
# Run from the repository root:  source("data-raw/02_cptac.R")      Requires base R only.
# =====================================================================================================

SRC_URL <- "https://linkedomics.org/data_download/CPTAC-PDAC/"
FILES <- c(clinical = "clinical_table_140.tsv",
           protein_tumor = "proteomics_gene_level_MD_abundance_tumor.cct",
           protein_normal = "proteomics_gene_level_MD_abundance_normal.cct",
           rna_tumor = "mRNA_RSEM_UQ_log2_Tumor.cct",
           rna_normal = "mRNA_RSEM_UQ_log2_Normal.cct")
SRC_DIR <- file.path("data-raw", "source", "cptac"); OUT <- file.path("app", "data", "cptac.rds")
dir.create(SRC_DIR, recursive = TRUE, showWarnings = FALSE); dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
source(file.path("data-raw", "helpers.R"))
log_file <- start_log("02_cptac")

## ---- download (skipped for files already present) -------------------------------------------------------
for (f in FILES) {
  dest <- file.path(SRC_DIR, f)
  if (!file.exists(dest)) { logc(log_file, "downloading ", f); utils::download.file(paste0(SRC_URL, f), dest, mode = "wb", quiet = TRUE) }
  logc(log_file, "source ", f, " MD5 ", unname(tools::md5sum(dest)), " (", file.size(dest), " bytes)")
}

read_cct <- function(path) {
  d <- read.delim(path, check.names = FALSE, row.names = NULL, quote = "", comment.char = "", stringsAsFactors = FALSE,
                  na.strings = c("", "NA", "NaN", "N/A", "null", "na"))
  names(d) <- trimws(gsub('^"|"$', "", names(d)))
  genes <- trimws(as.character(d[[1]])); cols <- grep("^C3[NL]-[0-9]+$", names(d))
  x <- as.matrix(d[, cols]); storage.mode(x) <- "double"; dimnames(x) <- list(genes, names(d)[cols])
  stopifnot(!anyDuplicated(genes), !anyDuplicated(colnames(x)))
  x
}
mean_semicolon <- function(v) vapply(strsplit(as.character(v), ";"), function(z) {
  z <- suppressWarnings(as.numeric(trimws(z))); if (any(is.finite(z))) mean(z, na.rm = TRUE) else NA_real_ }, numeric(1))

## ---- clinical ---------------------------------------------------------------------------------------------
cl <- read.delim(file.path(SRC_DIR, FILES[["clinical"]]), check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA", "na", "N/A"))
cl$case_id <- trimws(cl$case_id)
pdac <- cl$case_id[which(trimws(cl$histology_diagnosis) == "PDAC" & tolower(trimws(cl$tumor_included_for_the_study)) == "yes")]
norm_ok <- cl$case_id[which(tolower(trimws(cl$normal_included_for_the_study)) == "yes")]
stage <- sub("^Stage (IV|III|II|I)[AB]?$", "\\1", trimws(cl$tumor_stage_pathological)); stage[!stage %in% c("I", "II", "III", "IV")] <- NA
clinical <- data.frame(case_id = cl$case_id, histology = cl$histology_diagnosis, stage = stage,
                       stage_detail = cl$tumor_stage_pathological, sex = cl$sex, age = suppressWarnings(as.numeric(cl$age)),
                       stromal_fraction = mean_semicolon(cl$Stromal_fraction), neoplastic_cellularity = mean_semicolon(cl$Neoplastic_cellularity),
                       stringsAsFactors = FALSE)
clinical <- clinical[clinical$case_id %in% pdac, ]; rownames(clinical) <- NULL

## ---- matrices: PDAC tumours; adjacent normals only from PDAC cases flagged for the study -----------------
m <- lapply(FILES[-1], function(f) read_cct(file.path(SRC_DIR, f)))
keep_t <- function(x) x[, colnames(x) %in% pdac, drop = FALSE]
keep_n <- function(x) x[, colnames(x) %in% intersect(pdac, norm_ok), drop = FALSE]
out <- list(protein_tumor = keep_t(m$protein_tumor), protein_normal = keep_n(m$protein_normal),
            rna_tumor = keep_t(m$rna_tumor), rna_normal = keep_n(m$rna_normal))
out <- lapply(out, function(x) round(x, 3))
for (k in names(out)) logc(log_file, k, ": ", nrow(out[[k]]), " genes x ", ncol(out[[k]]), " samples")
logc(log_file, "protein pairs: ", length(intersect(colnames(out$protein_tumor), colnames(out$protein_normal))),
     " | RNA pairs: ", length(intersect(colnames(out$rna_tumor), colnames(out$rna_normal))))
out$clinical <- clinical
out$source <- list(citation = "Cao L. et al. Proteogenomic characterization of pancreatic ductal adenocarcinoma. Cell 2021;184:5031-5052.",
                   url = SRC_URL, units = c(protein = "log2 median-normalised MS intensity (TMT, gene level)", rna = "log2 RSEM upper-quartile"),
                   files = FILES, md5 = vapply(FILES, function(f) unname(tools::md5sum(file.path(SRC_DIR, f))), character(1)), built = format(Sys.time(), tz = "UTC"))
saveRDS(out, OUT, compress = "xz")
logc(log_file, "wrote ", OUT, " (", round(file.size(OUT) / 1e6, 1), " MB)")
end_log(log_file)
