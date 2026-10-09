# =====================================================================================================
# 01_scatlas_summaries.R
# Builds the single-cell data layer of the PDAC tissue explorer from the Loveless et al. (2025) atlas.
#
# Input  : the filtered atlas object produced by 00_scatlas_readiness.R (sc_ref_base.rds, Seurat v4),
#          plus the study's per-cell subcluster labels (DuctalTumor_Subcluster.txt.gz,
#          Fibroblast_Subcluster.txt.gz from https://github.com/PDAC-MULTIOMIC/PDAC_Atlas).
# Output : app/data/sc_groups.rds        one row per (sample x fine cell group) with metadata and n_cells
#          app/data/sc_expr.rds          genes x groups: mean CP10k (integer, x100) and fraction expressing
#                                        (integer, x10000) -> any coarser grouping is an exact weighted mean
#          app/data/sc_umap.rds          stratified subsample of cells: UMAP coordinates, metadata, and a
#                                        sparse cells x genes log1p(CP10k) matrix for feature plots
#          data-raw/logs/01_scatlas_summaries_<date>.txt   provenance (input checksums, counts, sessionInfo)
#
# Run from the repository root:  source("data-raw/01_scatlas_summaries.R")
# Requires: Seurat/SeuratObject (to read the object), Matrix. Peak memory ~45-50 GB (dominated by readRDS).
# =====================================================================================================

## ---- paths: edit to your local copies (default: data-raw/source/scatlas/, excluded from Git) -------
ATLAS_RDS     <- "data-raw/source/scatlas/sc_ref_base.rds"
DUCTAL_LABELS <- "data-raw/source/scatlas/DuctalTumor_Subcluster.txt.gz"   # column 'seurat_clusters': Ductal 1..9
FIBRO_LABELS  <- "data-raw/source/scatlas/Fibroblast_Subcluster.txt.gz"    # column 'Coarse': myFb 1..5, iFb 1..3
OUT_DIR       <- file.path("app", "data")
LOG_DIR       <- file.path("data-raw", "logs")

## ---- settings ------------------------------------------------------------------------------------------
MIN_FEATURES    <- 200     # light per-cell QC: cells with fewer detected genes are dropped
MIN_GENE_CELLS  <- 50      # genes detected in fewer cells overall are not shipped
UMAP_CELLS      <- 40000   # size of the stratified subsample shipped for UMAP feature plots
UMAP_MIN_TYPE   <- 1500    # every coarse cell type gets at least this many cells (or all it has)
set.seed(20261009)

suppressPackageStartupMessages({ library(Matrix); library(SeuratObject) })
stopifnot(file.exists(ATLAS_RDS), file.exists(DUCTAL_LABELS), file.exists(FIBRO_LABELS))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE); dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(LOG_DIR, paste0("01_scatlas_summaries_", format(Sys.Date()), ".txt"))
logc <- function(...) { msg <- paste0(...); cat(msg, "\n"); cat(msg, "\n", file = log_file, append = TRUE) }
cat("", file = log_file)
logc("01_scatlas_summaries.R started ", format(Sys.time(), tz = "UTC"), " UTC")
logc("label file MD5: ductal ", unname(tools::md5sum(DUCTAL_LABELS)), " | fibroblast ", unname(tools::md5sum(FIBRO_LABELS)))
logc("atlas file size (bytes): ", file.size(ATLAS_RDS))

## ---- 1. read the object and keep only what is needed --------------------------------------------------
logc("reading atlas object (this is the slow, memory-heavy step) ...")
obj <- readRDS(ATLAS_RDS)
assay <- obj[["RNA"]]
dat <- tryCatch(SeuratObject::LayerData(obj, assay = "RNA", layer = "data"), error = function(e) slot(assay, "data"))
meta <- obj@meta.data
umap <- tryCatch(obj@reductions$umap@cell.embeddings, error = function(e) NULL)
rm(obj, assay); invisible(gc())
stopifnot(inherits(dat, "dgCMatrix"), identical(colnames(dat), rownames(meta)))
logc("object: ", nrow(dat), " genes x ", ncol(dat), " cells; UMAP present: ", !is.null(umap))

# sanity: 'data' must be log1p(CP10k) (checked on 2,000 random cells)
chk <- sample(ncol(dat), 2000)
cs <- Matrix::colSums(expm1(dat[, chk]))
logc("normalisation check, median column sum of expm1(data) = ", round(median(cs), 1), " (expect ~10000)")
if (abs(median(cs) - 1e4) > 100) stop("'data' layer is not log1p(CP10k); adapt the script before continuing.")

## ---- 2. cell metadata, labels and QC ------------------------------------------------------------------
read_lab <- function(path, col) {
  d <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!col %in% names(d)) stop("column '", col, "' not found in ", path)
  setNames(trimws(as.character(d[[col]])), rownames(d))
}
duct <- read_lab(DUCTAL_LABELS, "seurat_clusters"); fib <- read_lab(FIBRO_LABELS, "Coarse")
cells <- rownames(meta)
m <- data.frame(cell = cells,
                sample = as.character(meta$Name),
                study = as.character(meta$GSE.SRA..Study.),
                disease_state = as.character(meta$DiseaseState),
                met_site = as.character(meta$If.metastatic..location),
                treatment = as.character(meta$Treatment),
                cell_type = as.character(meta$Clusters),
                n_features = meta$nFeature_RNA, stringsAsFactors = FALSE)
m$met_site[!(m$disease_state == "Metastatic lesion") | m$met_site %in% c("2", "NA", "")] <- NA   # stray '2' level in the source
m$treatment[m$disease_state == "Donor"] <- "Donor (untreated)"
m$treatment[grepl("^Treatment na.{1,2}ve$", m$treatment)] <- "Treatment-naive"   # locale-proof match of "naïve"
m$treatment[m$treatment == "Treatment"] <- "Treated"
m$subcluster <- NA_character_
m$subcluster[m$cell_type == "DUCTAL"] <- ifelse(cells[m$cell_type == "DUCTAL"] %in% names(duct), duct[cells[m$cell_type == "DUCTAL"]], "Ductal (unlabelled)")
m$subcluster[m$cell_type == "FIBROBLASTS"] <- ifelse(cells[m$cell_type == "FIBROBLASTS"] %in% names(fib), fib[cells[m$cell_type == "FIBROBLASTS"]], "Fibroblast (unlabelled)")
m$fine_group <- ifelse(is.na(m$subcluster), m$cell_type, m$subcluster)
logc("ductal labels matched: ", sum(m$cell_type == "DUCTAL" & m$subcluster != "Ductal (unlabelled)"),
     " | fibroblast labels matched: ", sum(m$cell_type == "FIBROBLASTS" & m$subcluster != "Fibroblast (unlabelled)"))

keep <- m$n_features >= MIN_FEATURES
logc("cells kept with >= ", MIN_FEATURES, " detected genes: ", sum(keep), " of ", length(keep))
if (sum(keep) < 0.5 * length(keep)) stop("QC would remove more than half of the cells; check MIN_FEATURES and nFeature_RNA.")
dat <- dat[, keep]; m <- m[keep, ]; if (!is.null(umap)) umap <- umap[m$cell, , drop = FALSE]
invisible(gc())

# genes detected in enough cells
nnz_gene <- tabulate(dat@i + 1L, nbins = nrow(dat))
gkeep <- nnz_gene >= MIN_GENE_CELLS
logc("genes detected in >= ", MIN_GENE_CELLS, " cells: ", sum(gkeep), " of ", nrow(dat))
dat <- dat[gkeep, ]; invisible(gc())

## ---- 3. UMAP subsample (before dat is modified in place) ----------------------------------------------
if (!is.null(umap)) {
  ct <- m$cell_type; tab <- table(ct)
  quota <- pmin(tab, pmax(UMAP_MIN_TYPE, round(UMAP_CELLS * tab / sum(tab))))
  pick <- unlist(lapply(names(quota), function(k) { i <- which(ct == k); if (length(i) <= quota[[k]]) i else sample(i, quota[[k]]) }))
  pick <- sort(pick)
  X <- t(dat[, pick]); X@x <- round(X@x, 2)                       # cells x genes, gene = column (fast lookup)
  X <- as(X, "CsparseMatrix")
  umap_obj <- list(coords = data.frame(UMAP_1 = umap[pick, 1], UMAP_2 = umap[pick, 2]),
                   meta = m[pick, c("sample", "study", "disease_state", "met_site", "treatment", "cell_type", "fine_group")],
                   expr = X, note = "stratified subsample by coarse cell type; values log1p(CP10k) rounded to 2 decimals")
  rownames(umap_obj$meta) <- NULL
  saveRDS(umap_obj, file.path(OUT_DIR, "sc_umap.rds"), compress = "xz")
  logc("UMAP subsample saved: ", length(pick), " cells; per type: ", paste(names(quota), as.integer(quota), sep = "=", collapse = ", "))
  rm(X, umap_obj); invisible(gc())
}

## ---- 4. per (sample x fine group): mean CP10k and fraction expressing ---------------------------------
key <- paste(m$sample, m$fine_group, sep = "\r")
g <- factor(key); ng <- nlevels(g)
ind <- sparseMatrix(i = seq_along(g), j = as.integer(g), x = 1, dims = c(length(g), ng))   # cells x groups
n_cells <- as.numeric(Matrix::colSums(ind))

logc("computing mean CP10k per group (", ng, " groups) ...")
dat@x <- expm1(dat@x)                                  # in place: log1p(CP10k) -> CP10k
mean_cp10k <- as.matrix(dat %*% ind); mean_cp10k <- sweep(mean_cp10k, 2, n_cells, "/")
logc("computing fraction expressing ...")
dat@x[] <- 1
frac <- as.matrix(dat %*% ind); frac <- sweep(frac, 2, n_cells, "/")
genes <- rownames(dat); rm(dat, ind); invisible(gc())

parts <- do.call(rbind, strsplit(levels(g), "\r", fixed = TRUE))
first <- !duplicated(key); smeta <- m[first, c("sample", "study", "disease_state", "met_site", "treatment")]
smeta <- smeta[!duplicated(smeta$sample), ]
groups <- data.frame(group_id = seq_len(ng), sample = parts[, 1], fine_group = parts[, 2], n_cells = as.integer(n_cells), stringsAsFactors = FALSE)
groups$cell_type <- m$cell_type[match(levels(g), key)]
groups <- cbind(groups, smeta[match(groups$sample, smeta$sample), c("study", "disease_state", "met_site", "treatment")])
rownames(groups) <- NULL

expr <- list(genes = genes,
             mean_cp10k_x100 = matrix(as.integer(round(mean_cp10k * 100)), nrow = length(genes), dimnames = list(genes, NULL)),
             frac_x10000 = matrix(as.integer(round(frac * 10000)), nrow = length(genes), dimnames = list(genes, NULL)),
             note = "columns follow sc_groups$group_id; divide by 100 / 10000; average expression for display = log1p(mean CP10k)")
saveRDS(groups, file.path(OUT_DIR, "sc_groups.rds"), compress = "xz")
saveRDS(expr, file.path(OUT_DIR, "sc_expr.rds"), compress = "xz")
logc("saved sc_groups.rds (", nrow(groups), " groups, ", length(unique(groups$sample)), " samples) and sc_expr.rds (", length(genes), " genes)")
logc("cells per coarse type: ", paste(names(table(m$cell_type)), as.integer(table(m$cell_type)), sep = "=", collapse = ", "))
logc("samples per disease state: ", paste(names(table(smeta$disease_state)), as.integer(table(smeta$disease_state)), sep = "=", collapse = ", "))
for (f in list.files(OUT_DIR, pattern = "^sc_", full.names = TRUE)) logc("output ", basename(f), " ", round(file.size(f) / 1e6, 1), " MB, MD5 ", unname(tools::md5sum(f)))
logc("finished ", format(Sys.time(), tz = "UTC"), " UTC")
capture.output(sessionInfo(), file = log_file, append = TRUE)
