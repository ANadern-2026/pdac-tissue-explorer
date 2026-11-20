# =====================================================================================================
# 00_scatlas_readiness_lowmem.R
# Low-memory variant of 00_scatlas_readiness.R (the original is left untouched).
#
# Same filtering logic, same thresholds, same output files as the original:
#   A) remove cells whose coarse label is a cycling cluster;
#   B) decontX per sample (constrained by the coarse labels) to score ambient-RNA contamination;
#   C) remove "ambient outliers": contamination > 0.30 OR above the 98th percentile within each
#      sample x coarse label;
#   D-G) keep/drop lists, gene whitelists (side files only), the filtered object and a summary.
#
# What is different (memory only, not results):
#   * The 55 GB Seurat object is never subset() as a whole. The original copied the full object in
#     step B (sc_noCycl) and again in step F (sc_ref_base), so peak memory was ~2x the object.
#   * Step 0 takes the object apart: metadata, UMAP (and other reductions' embeddings), the RNA
#     'counts' matrix and the RNA 'data' matrix. Everything else (scale.data, other assays such as
#     SCT/integrated, graphs) is not needed by 01 or 07 and is released.
#   * The RNA 'data' layer (log1p CP10k) IS needed downstream: 01_scatlas_summaries.R and
#     07_cell_populations.R read it and 01 stops if it is not log1p(CP10k). It is therefore not
#     thrown away but taken out of memory for step B (decontX only needs counts) and restored for
#     step F, in one of two ways (DATA_MODE):
#       "auto"      (default) checks on 2,000 random cells that 'data' == log1p(counts / colSums * 1e4).
#                   If so, it is released now and recomputed exactly for the kept cells in step F
#                   (no disk needed). If not, it falls back to "disk".
#       "recompute" as above, but stops if the check fails.
#       "disk"      writes 'data' uncompressed to WORK_DIR, releases it, reads it back in step F.
#                   Needs free disk space about the size of the 'data' matrix.
#   * decontX runs on the per-sample counts matrix directly (no Seurat -> SingleCellExperiment
#     conversion); it is the same decontX call on the same counts and labels (decontX uses its own
#     fixed seed, 12345). Results are checkpointed per sample, so an interrupted run resumes.
#   * Step F rebuilds a slim Seurat object (RNA counts + data, full meta.data incl. decontX_contam,
#     reductions' embeddings, identities) from the kept cells only, after the full matrices are gone.
#
# Usage in RStudio, with the atlas already loaded as `scatlas` (otherwise it is read from IN_RDS):
#   1. Make sure nothing else holds the object: no copies under other names, and run e.g. `gc()`
#      once at the console first (the console's .Last.value can still point at the object if the
#      last thing you ran was the readRDS() call).
#   2. setwd() to the folder you ran the original from, then source("data-raw/00_scatlas_readiness_lowmem.R").
#   With FREE_INPUT_OBJECT <- TRUE the script removes `scatlas` from your workspace after step 0;
#   that is what frees the memory.
# Requires Seurat/SeuratObject, celda, Matrix, dplyr, stringr, data.table, jsonlite
# (SingleCellExperiment is no longer needed).
# =====================================================================================================

suppressPackageStartupMessages({
  library(SeuratObject); library(celda); library(Matrix)
  library(dplyr); library(stringr); library(data.table); library(jsonlite)
})
set.seed(42)

## ---- paths and settings ---------------------------------------------------------------------------------
IN_RDS            <- "data/scAtlas.rds"                 # used only if `scatlas` is not already in the workspace
readiness_dir     <- file.path("results", "01_readiness")
WORK_DIR          <- file.path(readiness_dir, "lowmem_work")   # checkpoints and (DATA_MODE "disk") the stashed data layer
DATA_MODE         <- "auto"      # "auto", "recompute" or "disk" (see header)
FREE_INPUT_OBJECT <- TRUE        # remove `scatlas` from the workspace after step 0 (needed to free memory)
KEEP_COUNTS       <- TRUE        # TRUE = sc_ref_base.rds holds counts + data, like the original.
                                 # FALSE = data only: 01 and 07 do not read counts, and the file is about
                                 # half the size, so their readRDS() peak drops too. Not what the original wrote.
SCALE_FACTOR      <- 1e4

abs_contam_cut <- 0.30   # absolute contamination cut
tail_q         <- 0.98   # top 2 % per (sample, label)
labels_with_Ig <- c("B CELLS", "PLASMA")

dir.create(readiness_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(WORK_DIR, "decontx_by_sample"), recursive = TRUE, showWarnings = FALSE)
stopifnot(DATA_MODE %in% c("auto", "recompute", "disk"))

mem_used <- function() { g <- gc(); paste0(round(sum(g[, 2]) / 1024, 1), " GB in use") }
gb <- function(x) format(object.size(x), units = "GB", standard = "legacy", digits = 3)
get_layer <- function(obj, layer) {
  m <- tryCatch(SeuratObject::LayerData(obj, assay = "RNA", layer = layer), error = function(e) NULL)
  if (is.null(m)) m <- slot(obj[["RNA"]], layer)
  if (!inherits(m, "dgCMatrix")) m <- as(m, "CsparseMatrix")
  m
}
# log1p(counts / colSums * scale) column block by column block, filling one preallocated vector
# (avoids several nnz-sized temporaries at once). The result shares @i and @p with `cnt` in memory.
lognorm <- function(cnt, scale = SCALE_FACTOR, block = 20000L) {
  lib <- Matrix::colSums(cnt); p <- cnt@p; nc <- ncol(cnt)
  newx <- numeric(length(cnt@x))
  for (st in seq(1L, nc, by = block)) {
    en <- min(nc, st + block - 1L)
    r0 <- p[st] + 1; r1 <- p[en + 1]
    if (r1 >= r0) {
      r <- r0:r1
      newx[r] <- log1p(cnt@x[r] / rep.int(lib[st:en], diff(p[st:(en + 1)])) * scale)
    }
  }
  out <- cnt; out@x <- newx
  out
}

## ---- 0) take the object apart and release what is not needed --------------------------------------------
if (!exists("scatlas", envir = globalenv(), inherits = FALSE)) {
  message("`scatlas` not in the workspace, reading ", IN_RDS, " ...")
  assign("scatlas", readRDS(IN_RDS), envir = globalenv())
}
scatlas <- get("scatlas", envir = globalenv())
message("Input object: ", ncol(scatlas), " cells; assays: ", paste(names(scatlas@assays), collapse = ", "),
        "; reductions: ", paste(names(scatlas@reductions), collapse = ", "), " | ", mem_used())

meta <- scatlas@meta.data
broad_col  <- if ("Clusters" %in% colnames(meta)) "Clusters" else stop("Broad label column 'Clusters' not found.")
sample_col <- if ("Name" %in% colnames(meta)) "Name" else if ("GSE.SRA..Study." %in% colnames(meta)) "GSE.SRA..Study." else stop("No suitable sample column found.")

idents       <- tryCatch(Idents(scatlas), error = function(e) NULL)
project_name <- tryCatch(scatlas@project.name, error = function(e) "SeuratProject")
reds <- lapply(scatlas@reductions, function(r) list(
  emb = r@cell.embeddings, key = r@key,
  loadings = tryCatch(r@feature.loadings, error = function(e) NULL),
  stdev = tryCatch(r@stdev, error = function(e) numeric(0))))

cnt <- get_layer(scatlas, "counts")
dat <- get_layer(scatlas, "data")
stopifnot(identical(colnames(cnt), rownames(meta)), identical(dim(cnt), dim(dat)),
          identical(colnames(cnt), colnames(dat)), identical(rownames(cnt), rownames(dat)))
message("counts: ", gb(cnt), " | data: ", gb(dat))

# can 'data' be recomputed exactly from 'counts'?
use_recompute <- FALSE
if (DATA_MODE %in% c("auto", "recompute")) {
  chk <- sort(sample(ncol(cnt), min(2000L, ncol(cnt))))
  ref <- dat[, chk]; rec <- lognorm(cnt[, chk])
  same_pattern <- identical(ref@i, rec@i) && identical(ref@p, rec@p)
  max_diff <- if (same_pattern) max(abs(ref@x - rec@x)) else Inf
  message("data vs log1p(counts/colSums*", SCALE_FACTOR, ") on ", length(chk), " cells: same sparsity = ",
          same_pattern, ", max |diff| = ", signif(max_diff, 3))
  use_recompute <- same_pattern && max_diff < 1e-6
  if (!use_recompute && DATA_MODE == "recompute") stop("'data' cannot be recomputed from 'counts'; use DATA_MODE <- \"disk\".")
  rm(ref, rec)
}
data_stash <- file.path(WORK_DIR, "data_layer_full.rds")
if (!use_recompute) {
  message("'data' will be stashed to disk: ", data_stash, " (uncompressed, ~", gb(dat), ") ...")
  saveRDS(dat, data_stash, compress = FALSE)
}
fwrite(data.frame(data_mode = if (use_recompute) "recompute" else "disk"), file.path(WORK_DIR, "data_mode_used.csv"))

rm(dat)
if (FREE_INPUT_OBJECT) {
  rm("scatlas", envir = globalenv())
  message("`scatlas` removed from the workspace.")
}
if (exists("scatlas", inherits = FALSE)) rm(scatlas)
invisible(gc()); invisible(gc())
message("After step 0: ", mem_used(), " (should be roughly the size of 'counts' plus metadata)")

## ---- A) drop cycling labels ---------------------------------------------------------------------------
cycling_idx <- grepl("(?i)\\bCYCLING\\b|(?i)^CYCLING|(?i)CYCLING\\.", meta[[broad_col]])
cycling_labels_present <- sort(unique(meta[[broad_col]][cycling_idx]))
cycling_cells <- rownames(meta)[cycling_idx]
fwrite(data.frame(label = cycling_labels_present), file.path(readiness_dir, "cycling_labels_found.csv"))
message("Cycling labels detected: ", paste(cycling_labels_present, collapse = ", "), " | cells: ", length(cycling_cells))

## ---- B) decontX per sample (non-cycling cells), on the counts matrix, checkpointed ----------------------
nc_cells <- setdiff(rownames(meta), cycling_cells)            # same order as colnames(sc_noCycl) in the original
meta_nc  <- meta[nc_cells, , drop = FALSE]
nc_col   <- match(nc_cells, colnames(cnt))
samples  <- unique(meta_nc[[sample_col]])
contam_vec <- setNames(numeric(length(nc_cells)), nc_cells)
for (k in seq_along(samples)) {
  s <- samples[k]
  idx <- which(meta_nc[[sample_col]] == s)
  if (length(idx) < 50) { contam_vec[idx] <- NA_real_; next }          # too small to estimate robustly
  ck <- file.path(WORK_DIR, "decontx_by_sample", sprintf("sample_%04d.rds", k))
  if (file.exists(ck)) {
    r <- readRDS(ck)
    if (identical(r$sample, s) && identical(r$cells, nc_cells[idx])) { contam_vec[idx] <- r$contam; next }
    message("Checkpoint ", basename(ck), " does not match sample ", s, "; recomputing.")
  }
  t0 <- Sys.time()
  cnt_s <- cnt[, nc_col[idx], drop = FALSE]
  z_s   <- meta_nc[idx, broad_col, drop = TRUE]
  res   <- celda::decontX(cnt_s, z = z_s, estimateDelta = TRUE, verbose = TRUE)
  cs    <- as.numeric(res$contamination)
  stopifnot(length(cs) == length(idx))
  contam_vec[idx] <- cs
  saveRDS(list(sample = s, cells = nc_cells[idx], contam = cs), ck)
  rm(cnt_s, res); invisible(gc())
  message(sprintf("[%d/%d] %s: %d cells, %.1f min | %s", k, length(samples), s, length(idx),
                  as.numeric(difftime(Sys.time(), t0, units = "mins")), mem_used()))
}
saveRDS(contam_vec, file.path(readiness_dir, "contamination_scores.rds"))
meta$decontX_contam <- NA_real_
meta[names(contam_vec), "decontX_contam"] <- contam_vec

## ---- C) flag & remove severe ambient outliers ------------------------------------------------------------
df_contam <- data.frame(cell = names(contam_vec), sample_id = meta_nc[[sample_col]], broad_label = meta_nc[[broad_col]],
                        contam = contam_vec, stringsAsFactors = FALSE)
thr_tbl <- df_contam %>% group_by(sample_id, broad_label) %>%
  summarise(n = dplyr::n(), p98 = suppressWarnings(as.numeric(quantile(contam, tail_q, na.rm = TRUE))),
            thr_union = pmin(abs_contam_cut, p98), .groups = "drop")
df_contam <- df_contam %>% left_join(thr_tbl, by = c("sample_id", "broad_label")) %>%
  mutate(decontX_outlier = ifelse(is.na(contam) | is.na(thr_union), FALSE, contam > thr_union))
fwrite(df_contam, file.path(readiness_dir, "decontx_per_cell_nonCycl.csv"))
fwrite(thr_tbl,   file.path(readiness_dir, "decontx_thresholds_by_sample_label.csv"))
ambient_cells <- df_contam$cell[df_contam$decontX_outlier %in% TRUE]
message("Ambient-outlier cells to remove: ", length(ambient_cells))

## ---- D) keep/drop lists --------------------------------------------------------------------------------------
drop_df <- rbind(
  data.frame(cell = cycling_cells, reason = "cycling_label", contam = NA_real_,
             sample_id = meta[cycling_cells, sample_col, drop = TRUE], broad_label = meta[cycling_cells, broad_col, drop = TRUE]),
  df_contam %>% filter(decontX_outlier) %>% transmute(cell, reason = "ambient_outlier", contam, sample_id, broad_label))
drop_df <- drop_df[!duplicated(drop_df$cell), ]
keep_cells <- setdiff(rownames(meta), drop_df$cell)
fwrite(drop_df, file.path(readiness_dir, "cells_to_drop.csv")); writeLines(keep_cells, file.path(readiness_dir, "cells_keep.txt"))

## ---- E) gene whitelists (side files; the object keeps all genes) ----------------------------------------------
genes_all <- rownames(cnt)
is_mt <- grepl("^MT-", genes_all); is_rps <- grepl("^RPS\\d", genes_all); is_rpl <- grepl("^RPL\\d", genes_all); is_hb <- grepl("^HB[ABDEGQMZ]", genes_all)
is_ig <- grepl("^IGH[VDJC]", genes_all) | grepl("^IGK", genes_all) | grepl("^IGL", genes_all) |
  genes_all %in% c("IGKC", "IGLC1", "IGLC2", "IGLC3", "IGLC6", "IGHM", "IGHD", "IGHE", "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4")
exclude_global <- is_mt | is_rps | is_rpl | is_hb
fwrite(data.frame(gene = genes_all[!exclude_global]), file.path(readiness_dir, "whitelist_withIg.csv"))
fwrite(data.frame(gene = genes_all[!exclude_global & !is_ig]), file.path(readiness_dir, "whitelist_noIg.csv"))

## ---- G) summary (written before F so it is not lost if F runs out of memory) -----------------------------------
meta_keep <- meta[keep_cells, , drop = FALSE]
fwrite(as.data.frame(sort(table(meta[[broad_col]]), decreasing = TRUE)), file.path(readiness_dir, "label_sizes_before.csv"))
fwrite(as.data.frame(sort(table(meta_keep[[broad_col]]), decreasing = TRUE)), file.path(readiness_dir, "label_sizes_after.csv"))
write_json(list(input_cells = nrow(meta), cycling_removed = length(cycling_cells), ambient_removed = length(ambient_cells),
                kept = length(keep_cells), abs_cut = abs_contam_cut, tail_quantile = tail_q),
           file.path(readiness_dir, "00_readiness_summary.json"), auto_unbox = TRUE, pretty = TRUE)

## ---- F) filtered object, rebuilt from the kept cells only -----------------------------------------------------------
keep_col <- match(keep_cells, colnames(cnt))
cnt_keep <- cnt[, keep_col, drop = FALSE]
rm(cnt); invisible(gc())
message("Kept counts: ", gb(cnt_keep), " | ", mem_used())

if (use_recompute) {
  dat_keep <- lognorm(cnt_keep)
} else {
  tmp_cnt <- file.path(WORK_DIR, "counts_keep_tmp.rds")
  saveRDS(cnt_keep, tmp_cnt, compress = FALSE); rm(cnt_keep); invisible(gc())
  dat <- readRDS(data_stash); dat_keep <- dat[, keep_col, drop = FALSE]; rm(dat); invisible(gc())
  cnt_keep <- readRDS(tmp_cnt)
}
stopifnot(identical(colnames(dat_keep), keep_cells), identical(rownames(dat_keep), rownames(cnt_keep)))
message("Kept data: ", gb(dat_keep), " | ", mem_used())

op <- options(Seurat.object.assay.version = "v3")    # classic Assay, as in the original Seurat v4 object
if (KEEP_COUNTS) {
  assay <- CreateAssayObject(counts = cnt_keep)
  slot(assay, "data") <- dat_keep
} else {
  assay <- CreateAssayObject(data = dat_keep)
}
rm(cnt_keep, dat_keep); invisible(gc())
stopifnot(identical(rownames(assay), genes_all), identical(colnames(assay), keep_cells))
sc_ref_base <- CreateSeuratObject(counts = assay, project = project_name, assay = "RNA")
rm(assay); options(op)
sc_ref_base@meta.data <- meta_keep                    # original metadata, incl. nCount/nFeature_RNA and decontX_contam
if (!is.null(idents)) Idents(sc_ref_base) <- droplevels(idents[keep_cells])
for (nm in names(reds)) {
  r <- reds[[nm]]
  emb <- r$emb[keep_cells, , drop = FALSE]
  load <- if (!is.null(r$loadings) && length(r$loadings)) r$loadings else new("matrix")
  sc_ref_base[[nm]] <- tryCatch(
    CreateDimReducObject(embeddings = emb, loadings = load, stdev = r$stdev, key = r$key, assay = "RNA"),
    error = function(e) CreateDimReducObject(embeddings = emb, key = r$key, assay = "RNA"))
}
rm(reds); invisible(gc())

# checks equivalent to what 01 / 07 rely on
chk_dat <- get_layer(sc_ref_base, "data")
chk <- sample(ncol(chk_dat), min(2000L, ncol(chk_dat)))
cs <- Matrix::colSums(expm1(chk_dat[, chk])); rm(chk_dat)
stopifnot(abs(median(cs) - SCALE_FACTOR) < 100, identical(colnames(sc_ref_base), keep_cells),
          "decontX_contam" %in% colnames(sc_ref_base@meta.data))
message("Output check: median column sum of expm1(data) = ", round(median(cs), 1),
        "; UMAP present: ", "umap" %in% names(sc_ref_base@reductions))

saveRDS(sc_ref_base, file = file.path(readiness_dir, "sc_ref_base.rds"))
message("Saved ", file.path(readiness_dir, "sc_ref_base.rds"), " | kept cells: ", length(keep_cells))
if (!use_recompute) message("You can delete the stashed files in ", WORK_DIR, " (data_layer_full.rds, counts_keep_tmp.rds).")
message("Per-sample decontX checkpoints are in ", file.path(WORK_DIR, "decontx_by_sample"), "; delete them once you are happy with the result.")
