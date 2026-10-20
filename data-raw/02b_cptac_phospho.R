# =====================================================================================================
# 02b_cptac_phospho.R
# CPTAC PDAC phosphoproteome (Cao et al., Cell 2021) via LinkedOmics: site-level, median-normalised log2
# intensities of tumours and adjacent normal tissue. Same sample selection as 02_cptac.R (PDAC cases,
# adenosquamous excluded; normals only from PDAC cases flagged for the study). Writes app/data/cptac_phospho.rds.
# Sites are identified as in the source: RefSeq protein accession + residue + position (e.g. NP_000413.1_S10).
# -Inf entries in the source (log2 of zero) are treated as not quantified.
# Run from the repository root:  source("data-raw/02b_cptac_phospho.R")      Requires base R only.
# =====================================================================================================

SRC_URL <- "https://linkedomics.org/data_download/CPTAC-PDAC/"
FILES <- c(clinical = "clinical_table_140.tsv",
           phospho_tumor = "phosphoproteomics_site_level_MD_abundance_tumor.cct",
           phospho_normal = "phosphoproteomics_site_level_MD_abundance_normal.cct")
SRC_DIR <- file.path("data-raw", "source", "cptac"); OUT <- file.path("app", "data", "cptac_phospho.rds")
dir.create(SRC_DIR, recursive = TRUE, showWarnings = FALSE); dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
source(file.path("data-raw", "helpers.R"))
log_file <- start_log("02b_cptac_phospho")
old_timeout <- options(timeout = max(1200, getOption("timeout")))

## ---- download (skipped for files already present) -------------------------------------------------------
for (f in FILES) {
  dest <- file.path(SRC_DIR, f)
  if (!file.exists(dest)) { logc(log_file, "downloading ", f); utils::download.file(paste0(SRC_URL, f), dest, mode = "wb", quiet = TRUE) }
  logc(log_file, "source ", f, " MD5 ", unname(tools::md5sum(dest)), " (", file.size(dest), " bytes)")
}
options(old_timeout)

## ---- clinical: same selection as 02_cptac.R ----------------------------------------------------------------
cl <- read.delim(file.path(SRC_DIR, FILES[["clinical"]]), check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA", "na", "N/A"))
cl$case_id <- trimws(cl$case_id)
pdac <- cl$case_id[which(trimws(cl$histology_diagnosis) == "PDAC" & tolower(trimws(cl$tumor_included_for_the_study)) == "yes")]
norm_ok <- cl$case_id[which(tolower(trimws(cl$normal_included_for_the_study)) == "yes")]

## ---- site tables --------------------------------------------------------------------------------------------
read_sites <- function(path) {
  d <- read.delim(path, check.names = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE, na.strings = c("", "NA", "NaN", "N/A", "null", "na"))
  names(d) <- trimws(gsub('^"|"$', "", names(d)))
  stopifnot(identical(names(d)[1:3], c("Index", "Gene", "Peptide")))
  cols <- grep("^C3[NL]-[0-9]+$", names(d))
  x <- as.matrix(d[, cols]); storage.mode(x) <- "double"; dimnames(x) <- list(trimws(d$Index), names(d)[cols])
  stopifnot(!anyDuplicated(rownames(x)), !anyDuplicated(colnames(x)))
  list(x = x, gene = trimws(d$Gene), peptide = d$Peptide)
}
tu <- read_sites(file.path(SRC_DIR, FILES[["phospho_tumor"]])); no <- read_sites(file.path(SRC_DIR, FILES[["phospho_normal"]]))
stopifnot(identical(rownames(tu$x), rownames(no$x)), identical(tu$gene, no$gene))
for (k in c("tumour", "normal")) { m <- if (k == "tumour") tu$x else no$x; logc(log_file, k, ": ", sum(is.infinite(m)), " -Inf values (in ", sum(rowSums(is.infinite(m)) > 0), " sites) set to NA") }
tu$x[is.infinite(tu$x)] <- NA; no$x[is.infinite(no$x)] <- NA

x_t <- round(tu$x[, colnames(tu$x) %in% pdac, drop = FALSE], 3)
x_n <- round(no$x[, colnames(no$x) %in% intersect(pdac, norm_ok), drop = FALSE], 3)

## ---- site annotation ------------------------------------------------------------------------------------------
id <- rownames(x_t)
re <- "^((?:NP|XP)_[0-9]+\\.[0-9]+)_([STY])([0-9]+)$"
if (!all(grepl(re, id, perl = TRUE))) stop("unexpected site identifiers, e.g. ", id[!grepl(re, id, perl = TRUE)][1])
pep <- vapply(strsplit(gsub(";+$", "", tu$peptide), ";", fixed = TRUE), function(p) paste(unique(trimws(p[nzchar(p)])), collapse = ";"), character(1))
sites <- data.frame(site_id = id, gene = tu$gene, refseq = sub(re, "\\1", id, perl = TRUE), residue = sub(re, "\\2", id, perl = TRUE),
                    position = as.integer(sub(re, "\\3", id, perl = TRUE)), peptide = pep, stringsAsFactors = FALSE)
if (anyNA(sites$gene) || any(!nzchar(sites$gene))) stop("sites without gene symbol")

keep <- rowSums(is.finite(x_t)) + rowSums(is.finite(x_n)) > 0               # quantified in at least one kept sample
sites <- sites[keep, ]; x_t <- x_t[keep, , drop = FALSE]; x_n <- x_n[keep, , drop = FALSE]
o <- order(sites$gene, sites$refseq, sites$position, sites$residue); sites <- sites[o, ]; x_t <- x_t[o, , drop = FALSE]; x_n <- x_n[o, , drop = FALSE]
rownames(sites) <- NULL

pairs <- intersect(colnames(x_t), colnames(x_n))
np <- rowSums(is.finite(x_t[, pairs, drop = FALSE]) & is.finite(x_n[, pairs, drop = FALSE]))
logc(log_file, "phospho_tumor: ", nrow(x_t), " sites x ", ncol(x_t), " samples | phospho_normal: ", ncol(x_n), " samples | pairs: ", length(pairs))
logc(log_file, "genes: ", length(unique(sites$gene)), " | residues: ", paste(names(t <- table(sites$residue)), t, sep = " = ", collapse = ", "),
     " | sites quantified in >= 10 pairs: ", sum(np >= 10), " | median pairs per site: ", stats::median(np))
cp <- file.path("app", "data", "cptac.rds")
if (file.exists(cp)) {
  pc <- readRDS(cp)
  same <- setequal(colnames(x_t), colnames(pc$protein_tumor)) && setequal(colnames(x_n), colnames(pc$protein_normal))
  logc(log_file, "samples identical to the protein data in cptac.rds: ", same)
  if (!same) warning("phospho samples differ from the protein samples in cptac.rds")
}

out <- list(sites = sites, phospho_tumor = x_t, phospho_normal = x_n,
            source = list(citation = "Cao L. et al. Proteogenomic characterization of pancreatic ductal adenocarcinoma. Cell 2021;184:5031-5052.",
                          url = SRC_URL, units = c(phospho = "log2 median-normalised MS intensity (TMT, phosphosite level)"),
                          files = FILES, md5 = vapply(FILES, function(f) unname(tools::md5sum(file.path(SRC_DIR, f))), character(1)),
                          built = format(Sys.time(), tz = "UTC")))
saveRDS(out, OUT, compress = "xz")
logc(log_file, "wrote ", OUT, " (", round(file.size(OUT) / 1e6, 1), " MB)")
end_log(log_file)
