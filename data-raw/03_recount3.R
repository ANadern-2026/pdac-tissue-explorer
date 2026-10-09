# =====================================================================================================
# 03_recount3.R
# TCGA-PAAD and GTEx pancreas, uniformly processed by recount3 (gene sums, GENCODE v26).
# Converts counts to TPM, keeps TCGA-PAAD primary tumours of ductal histology ("PDAC"), the TCGA
# adjacent-normal samples and GTEx pancreas, computes composition signature scores, and writes
# app/data/bulk_tcga_gtex.rds.
# Run from the repository root:  source("data-raw/03_recount3.R")
# Requires recount3 + SummarizedExperiment (Bioconductor) for the first run only; afterwards the
# processed caches in data-raw/source/recount3/ are reused.
# =====================================================================================================

SRC_DIR <- file.path("data-raw", "source", "recount3"); OUT <- file.path("app", "data", "bulk_tcga_gtex.rds")
dir.create(SRC_DIR, recursive = TRUE, showWarnings = FALSE); dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
source(file.path("data-raw", "helpers.R"))
log_file <- start_log("03_recount3")

SIGNATURES <- list(fibroblast = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRB", "FAP"),
                   endothelial = c("PECAM1", "VWF", "CDH5", "CLDN5"),
                   epithelial = c("KRT19", "EPCAM", "KRT8", "KRT18"),
                   acinar = c("PRSS1", "CPA1", "CELA3A", "CTRB1"),
                   immune = c("PTPRC"))

## ---- fetch: recount3 gene sums -> log2(TPM+1), one row per gene symbol (highest-expressed row kept) ------
fetch <- function(project, source, tag) {
  f <- file.path(SRC_DIR, paste0(tag, "_log2tpm.rds"))
  if (file.exists(f)) { logc(log_file, "using cached ", f, " MD5 ", unname(tools::md5sum(f))); return(readRDS(f)) }
  for (pk in c("recount3", "SummarizedExperiment")) if (!requireNamespace(pk, quietly = TRUE)) stop("Install ", pk, " (Bioconductor) for the first run.")
  hp <- recount3::available_projects()
  pi <- hp[hp$project == project & hp$file_source == source & hp$project_type == "data_sources", ]
  if (nrow(pi) != 1) stop("recount3 project not found: ", project, " / ", source)
  rse <- recount3::create_rse(pi)
  cnt <- recount3::transform_counts(rse); len <- SummarizedExperiment::rowRanges(rse)$bp_length
  rpk <- cnt / (len / 1000); tpm <- sweep(rpk, 2, colSums(rpk), "/") * 1e6
  sym <- as.character(SummarizedExperiment::rowRanges(rse)$gene_name); ok <- !is.na(sym) & sym != ""
  tpm <- tpm[ok, ]; sym <- sym[ok]; o <- order(sym, -rowMeans(tpm)); tpm <- tpm[o, ]; sym <- sym[o]
  d <- !duplicated(sym); tpm <- tpm[d, ]; rownames(tpm) <- sym[d]
  res <- list(log2tpm = log2(tpm + 1), meta = as.data.frame(SummarizedExperiment::colData(rse)), project = project, source = source)
  saveRDS(res, f); logc(log_file, "downloaded and cached ", f); res
}
tc <- fetch("PAAD", "tcga", "TCGA_PAAD"); gt <- fetch("PANCREAS", "gtex", "GTEx_pancreas")

## ---- cohorts ------------------------------------------------------------------------------------------------
cd <- tc$meta; gd <- gt$meta
st <- cd[["tcga.gdc_cases.samples.sample_type"]]; hist <- cd[["tcga.cgc_case_histological_diagnosis"]]
is_pdac <- grepl("ductal", hist, ignore.case = TRUE) & st == "Primary Tumor"
is_adjn <- st == "Solid Tissue Normal"
gt_ok <- gd[["gtex.smafrze"]] %in% c("RNASEQ", "USE ME")
stage_raw <- tolower(trimws(as.character(cd[["tcga.gdc_cases.diagnoses.tumor_stage"]])))
stage <- rep(NA_character_, length(stage_raw))
stage[grepl("stage i(a|b)?$", stage_raw)] <- "I"; stage[grepl("stage ii(a|b)?$", stage_raw)] <- "II"
stage[grepl("stage iii", stage_raw)] <- "III"; stage[grepl("stage iv", stage_raw)] <- "IV"

common <- intersect(rownames(tc$log2tpm), rownames(gt$log2tpm))
X <- cbind(gt$log2tpm[common, gt_ok], tc$log2tpm[common, is_adjn], tc$log2tpm[common, is_pdac])
samples <- data.frame(sample = colnames(X),
                      cohort = c(rep("Healthy pancreas (GTEx)", sum(gt_ok)), rep("Adjacent normal (TCGA)", sum(is_adjn)), rep("PDAC (TCGA)", sum(is_pdac))),
                      stage = c(rep(NA, sum(gt_ok) + sum(is_adjn)), stage[is_pdac]), stringsAsFactors = FALSE)
for (k in names(SIGNATURES)) samples[[paste0("score_", k)]] <- round(colMeans(X[intersect(SIGNATURES[[k]], rownames(X)), , drop = FALSE]), 3)

## ---- genes: keep those expressed (log2 TPM+1 >= 1, i.e. TPM >= 1) in >= 10 % of samples of any cohort ---
expressed <- Reduce(`|`, lapply(split(seq_len(ncol(X)), samples$cohort), function(j) rowMeans(X[, j, drop = FALSE] >= 1) >= 0.10))
X <- X[expressed, ]
logc(log_file, "cohorts: ", paste(names(table(samples$cohort)), as.integer(table(samples$cohort)), sep = " = ", collapse = "; "))
logc(log_file, "genes kept (TPM >= 1 in >= 10% of any cohort): ", nrow(X), " of ", length(common))
out <- list(log2tpm_x1000 = matrix(as.integer(round(X * 1000)), nrow = nrow(X), dimnames = dimnames(X)), samples = samples,
            signatures = SIGNATURES,
            source = list(citation = c("Wilks C. et al. recount3: summaries and queries for large-scale RNA-seq expression and splicing. Genome Biol 2021;22:323.",
                                       "The Cancer Genome Atlas Research Network. Integrated genomic characterization of pancreatic ductal adenocarcinoma. Cancer Cell 2017;32:185-203.",
                                       "GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science 2020;369:1318-1330."),
                          units = "log2(TPM + 1); stored x1000 as integer", annotation = "GENCODE v26 gene sums (recount3)", built = format(Sys.time(), tz = "UTC")))
saveRDS(out, OUT, compress = "xz")
logc(log_file, "wrote ", OUT, " (", round(file.size(OUT) / 1e6, 1), " MB)")
end_log(log_file)
