# =====================================================================================================
# 03b_clinical.R
# Clinical context for the tumours of the two bulk cohorts: survival, pathological stage, age, sex and a
# classical vs basal-like transcriptional score.
#   CPTAC-PDAC (Cao et al., Cell 2021): LinkedOmics clinical table (same file as 02_cptac.R): overall survival from
#     follow_up_days + vital_status; pathological stage.
#   TCGA-PAAD: survival endpoints and stage from the TCGA Pan-Cancer Clinical Data Resource (TCGA-CDR; Liu et al.,
#     Cell 2018; overall survival OS and progression-free interval PFI, both recommended for PAAD), linked to the
#     recount3 samples of 03_recount3.R through the recount3 sample metadata (sample ID -> TCGA barcode).
#   Moffitt score (Moffitt et al., Nat Genet 2015): per tumour, mean z-score of the 25 basal-like genes minus mean
#     z-score of the 25 classical genes (z across the tumours of the cohort; RNA, log2 scale). Higher = more basal-like.
# Inputs : app/data/cptac.rds, app/data/bulk_tcga_gtex.rds (02, 03) and the source files below (downloaded if absent)
# Output : app/data/clinical.rds
# Run from the repository root:  source("data-raw/03b_clinical.R")     Requires readxl; a minute.
# =====================================================================================================

SRC_CPTAC <- file.path("data-raw", "source", "cptac"); SRC_TCGA <- file.path("data-raw", "source", "tcga")
IN_DIR <- Sys.getenv("PDACX_DATA_DIR", file.path("app", "data")); OUT <- file.path(IN_DIR, "clinical.rds")
dir.create(SRC_CPTAC, recursive = TRUE, showWarnings = FALSE); dir.create(SRC_TCGA, recursive = TRUE, showWarnings = FALSE)
source(file.path("data-raw", "helpers.R")); log_file <- start_log("03b_clinical")
if (!requireNamespace("readxl", quietly = TRUE)) stop("Install the readxl package: install.packages(\"readxl\")")

## ---- sources (downloaded once; checksums logged) ------------------------------------------------------------
RECOUNT3_MIRRORS <- c("https://recount-opendata.s3.amazonaws.com/recount3/release/human", "http://duffel.rail.bio/recount3/human",
                      "https://data.idies.jhu.edu/recount3/data/human")
SOURCES <- list(
  cptac_clinical = list(file = file.path(SRC_CPTAC, "clinical_table_140.tsv"), urls = "https://linkedomics.org/data_download/CPTAC-PDAC/clinical_table_140.tsv"),
  tcga_cdr = list(file = file.path(SRC_TCGA, "TCGA-CDR-SupplementalTableS1.xlsx"), urls = "https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81"),
  recount3_tcga = list(file = file.path(SRC_TCGA, "tcga.tcga.PAAD.MD.gz"), urls = paste0(RECOUNT3_MIRRORS, "/data_sources/tcga/metadata/AD/PAAD/tcga.tcga.PAAD.MD.gz")))
old_timeout <- options(timeout = max(600, getOption("timeout")))
for (s in SOURCES) {
  if (!file.exists(s$file) || file.size(s$file) == 0) {
    for (u in s$urls) {
      ok <- tryCatch({ utils::download.file(u, s$file, mode = "wb", quiet = TRUE); file.size(s$file) > 0 }, error = function(e) FALSE, warning = function(w) FALSE)
      if (isTRUE(ok)) { logc(log_file, "downloaded ", basename(s$file), " from ", u); break }
      unlink(s$file)
    }
    if (!file.exists(s$file)) stop("could not download ", basename(s$file))
  }
  logc(log_file, "source ", basename(s$file), " MD5 ", unname(tools::md5sum(s$file)), " (", file.size(s$file), " bytes)")
}
options(old_timeout)
md5 <- vapply(SOURCES, function(s) unname(tools::md5sum(s$file)), character(1))

cptac <- readRDS(file.path(IN_DIR, "cptac.rds")); bulk <- readRDS(file.path(IN_DIR, "bulk_tcga_gtex.rds"))
if (!identical(unname(md5[["cptac_clinical"]]), unname(cptac$source$md5[["clinical"]])))
  stop("clinical_table_140.tsv differs from the file used for cptac.rds (02_cptac.R); rerun 02 or restore the file")

stage_main <- function(x) {                                    # "Stage IIB" -> "II"; anything else -> NA
  s <- sub("^Stage (IV|III|II|I)[ABC]?$", "\\1", trimws(x)); s[!s %in% c("I", "II", "III", "IV")] <- NA; s
}
num <- function(x) suppressWarnings(as.numeric(x))

## ---- Moffitt classical / basal-like genes (Moffitt et al. 2015, Fig. 1; current HGNC symbols) -----------------
MOFFITT <- list(
  basal = c("VGLL1", "UCA1", "S100A2", "LY6D", "SPRR3", "SPRR1B", "LEMD1", "KRT15", "CTSV", "DHRS9", "AREG", "CST6", "SERPINB3",
            "KRT6C", "KRT6A", "SERPINB4", "FAM83A", "SCEL", "FGFBP1", "KRT7", "KRT17", "GPR87", "TNS4", "SLC2A1", "ANXA8L1"),
  classical = c("BTNL8", "FAM3D", "PRR15L", "AGR3", "CTSE", "LOC400573", "LYZ", "TFF2", "TFF1", "ANXA10", "LGALS4", "PLA2G10", "CEACAM6",
                "VSIG2", "TSPAN8", "ST6GALNAC1", "AGR2", "TFF3", "CYP3A7", "MYO1A", "CLRN3", "KRT20", "CDH17", "SPINK4", "REG4"))
ALIASES <- c(CTSV = "CTSL2", PRR15L = "ATAD4", ANXA8L1 = "ANXA8L2")      # symbols used in the 2015 paper
moffitt_score <- function(X, cohort) {                          # X: genes x tumours (log2 scale)
  resolve <- function(g) { a <- unname(ALIASES[g]); ifelse(g %in% rownames(X), g, ifelse(!is.na(a) & a %in% rownames(X), a, NA)) }
  pick <- function(g) { h <- stats::na.omit(resolve(g)); h[apply(X[h, , drop = FALSE], 1, function(v) isTRUE(stats::sd(v, na.rm = TRUE) > 0))] }
  z <- function(g) colMeans(t(scale(t(X[g, , drop = FALSE]))), na.rm = TRUE)
  b <- pick(MOFFITT$basal); cl <- pick(MOFFITT$classical)
  all <- c(MOFFITT$basal, MOFFITT$classical); miss <- all[is.na(resolve(all)) | !resolve(all) %in% c(b, cl)]
  logc(log_file, cohort, ": Moffitt score from ", length(b), " of 25 basal-like and ", length(cl), " of 25 classical genes",
       if (length(miss)) paste0(" (not available or constant: ", paste(miss, collapse = ", "), ")") else "")
  list(score = round(z(b) - z(cl), 3), genes = list(basal = b, classical = cl))
}

## ---- CPTAC ----------------------------------------------------------------------------------------------------
cl <- read.delim(SOURCES$cptac_clinical$file, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA", "na", "N/A"))
cl$case_id <- trimws(cl$case_id)
tum <- colnames(cptac$rna_tumor)
cl <- cl[match(cptac$clinical$case_id, cl$case_id), ]
stopifnot(identical(cl$case_id, cptac$clinical$case_id), setequal(tum, cl$case_id) || all(tum %in% cl$case_id))
vs <- trimws(cl$vital_status); fu <- num(cl$follow_up_days)
ms_c <- moffitt_score(cptac$rna_tumor, "CPTAC RNA")
cp <- data.frame(case_id = cl$case_id, age = num(cl$age), sex = cl$sex, stage = stage_main(cl$tumor_stage_pathological),
                 stage_detail = trimws(cl$tumor_stage_pathological), residual_tumor = sub(":.*$", "", trimws(cl$residual_tumor)),
                 os_time = ifelse(vs %in% c("Deceased", "Living"), fu, NA), os_event = ifelse(vs == "Deceased", 1L, ifelse(vs == "Living", 0L, NA_integer_)),
                 moffitt = unname(ms_c$score[cl$case_id]), stringsAsFactors = FALSE)
cp$os_event[is.na(cp$os_time)] <- NA_integer_
logc(log_file, "CPTAC: ", nrow(cp), " tumours; overall survival for ", sum(!is.na(cp$os_time)), " (", sum(cp$os_event, na.rm = TRUE), " deaths; median follow-up of survivors ",
     round(stats::median(cp$os_time[cp$os_event %in% 0L]) / (365.25 / 12), 1), " months); stage ",
     paste(names(t <- table(cp$stage, useNA = "ifany")), t, sep = " = ", collapse = ", "))

## ---- TCGA-PAAD ------------------------------------------------------------------------------------------------
md <- read.delim(gzfile(SOURCES$recount3_tcga$file), check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
names(md) <- sub("^tcga\\.", "", names(md))
smp <- bulk$samples[bulk$samples$cohort == "PDAC (TCGA)", ]
m <- md[match(smp$sample, md$external_id), ]
if (anyNA(m$external_id)) stop(sum(is.na(m$external_id)), " TCGA samples of bulk_tcga_gtex.rds are not in the recount3 metadata")
barcode <- m$tcga_barcode; patient <- substr(barcode, 1, 12)
stopifnot(!anyDuplicated(patient), all(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-01", barcode)))
cdr <- as.data.frame(readxl::read_excel(SOURCES$tcga_cdr$file, sheet = "TCGA-CDR", col_types = "text", .name_repair = "minimal", na = c("", "NA", "#N/A", "[Not Available]", "[Not Applicable]", "[Unknown]", "[Discrepancy]", "[Not Evaluated]")))
cdr <- cdr[cdr$type == "PAAD", ]
if (anyDuplicated(cdr$bcr_patient_barcode)) stop("TCGA-CDR: duplicated patient barcodes")
redacted <- cdr$bcr_patient_barcode[!is.na(cdr$Redaction)]
k <- match(patient, cdr$bcr_patient_barcode); k[patient %in% redacted] <- NA
logc(log_file, "TCGA: ", nrow(smp), " PDAC samples (", length(unique(patient)), " patients); in TCGA-CDR: ", sum(!is.na(k)),
     if (length(redacted)) paste0("; redacted cases excluded: ", sum(patient %in% redacted)) else "")
ms_t <- moffitt_score(bulk$log2tpm_x1000[, smp$sample] / 1000, "TCGA RNA")
tc <- data.frame(sample = smp$sample, patient = patient, barcode = barcode, age = num(cdr$age_at_initial_pathologic_diagnosis[k]),
                 sex = tools::toTitleCase(tolower(cdr$gender[k])), stage = stage_main(cdr$ajcc_pathologic_tumor_stage[k]),
                 stage_detail = cdr$ajcc_pathologic_tumor_stage[k],
                 os_time = num(cdr$OS.time[k]), os_event = as.integer(num(cdr$OS[k])), pfi_time = num(cdr$PFI.time[k]), pfi_event = as.integer(num(cdr$PFI[k])),
                 moffitt = unname(ms_t$score[smp$sample]), stringsAsFactors = FALSE)
for (ep in c("os", "pfi")) { bad <- is.na(tc[[paste0(ep, "_time")]]) | is.na(tc[[paste0(ep, "_event")]]); tc[bad, paste0(ep, c("_time", "_event"))] <- NA }
logc(log_file, "TCGA: OS for ", sum(!is.na(tc$os_time)), " (", sum(tc$os_event, na.rm = TRUE), " deaths), PFI for ", sum(!is.na(tc$pfi_time)), " (",
     sum(tc$pfi_event, na.rm = TRUE), " events); median follow-up of survivors ", round(stats::median(tc$os_time[tc$os_event %in% 0L]) / (365.25 / 12), 1),
     " months; stage ", paste(names(t <- table(tc$stage, useNA = "ifany")), t, sep = " = ", collapse = ", "))

## ---- write ----------------------------------------------------------------------------------------------------
out <- list(cptac = cp, tcga = tc, moffitt = MOFFITT, moffitt_used = list(cptac = ms_c$genes, tcga = ms_t$genes),
            source = list(citation = c("Cao L. et al. Proteogenomic characterization of pancreatic ductal adenocarcinoma. Cell 2021;184:5031-5052.",
                                       "Liu J. et al. An integrated TCGA pan-cancer clinical data resource to drive high-quality survival outcome analytics. Cell 2018;173:400-416.",
                                       "Moffitt R.A. et al. Virtual microdissection identifies distinct tumor- and stroma-specific subtypes of pancreatic ductal adenocarcinoma. Nat Genet 2015;47:1168-1178."),
                          units = c(time = "days", moffitt = "mean z-score of basal-like genes minus mean z-score of classical genes (RNA, z across the tumours of the cohort)"),
                          files = vapply(SOURCES, function(s) basename(s$file), character(1)), md5 = md5,
                          input_md5 = unname(tools::md5sum(file.path(IN_DIR, c("cptac.rds", "bulk_tcga_gtex.rds")))), built = format(Sys.time(), tz = "UTC")))
saveRDS(out, OUT, compress = "xz")
logc(log_file, "wrote ", OUT, " (", round(file.size(OUT) / 1e3), " kB)")
end_log(log_file)
