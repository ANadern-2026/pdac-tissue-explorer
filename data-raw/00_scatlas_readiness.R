# =====================================================================================================
# 00_scatlas_readiness.R
# Preprocessing of the Loveless et al. (2025) PDAC single-cell atlas (Zenodo 10.5281/zenodo.14199536,
# scAtlas.rds) that produced sc_ref_base.rds, the input of 01_scatlas_summaries.R.
#
#   A) removes cells whose coarse label is a cycling cluster;
#   B) runs decontX per sample (constrained by the coarse labels) to score ambient-RNA contamination;
#   C) removes "ambient outliers": contamination > 0.30 OR above the 98th percentile within each
#      sample x coarse label;
#   D-G) writes keep/drop lists, gene whitelists (side files only) and the filtered object.
# Counts and normalised values of retained cells are NOT altered (decontX scores only).
#
# This is the script used for the published data layer; it is included for transparency and
# reproducibility. Logic is unchanged from the original run; only the input/output paths were made
# editable at the top. Requires Seurat, SingleCellExperiment, celda, Matrix, dplyr, data.table, jsonlite.
# =====================================================================================================

suppressPackageStartupMessages({
  library(Seurat); library(SingleCellExperiment); library(celda); library(Matrix)
  library(dplyr); library(stringr); library(data.table); library(jsonlite)
})
set.seed(42)

## ---- paths -----------------------------------------------------------------------------------------------
IN_RDS        <- "data/scAtlas.rds"                 # the Zenodo object, decompressed
readiness_dir <- file.path("results", "01_readiness")
if (!dir.exists(readiness_dir)) dir.create(readiness_dir, recursive = TRUE)

scatlas <- readRDS(IN_RDS)
meta <- scatlas@meta.data
broad_col  <- if ("Clusters" %in% colnames(meta)) "Clusters" else stop("Broad label column 'Clusters' not found.")
sample_col <- if ("Name" %in% colnames(meta)) "Name" else if ("GSE.SRA..Study." %in% colnames(meta)) "GSE.SRA..Study." else stop("No suitable sample column found.")

abs_contam_cut <- 0.30   # absolute contamination cut
tail_q         <- 0.98   # top 2 % per (sample, label)
labels_with_Ig <- c("B CELLS", "PLASMA")

## ---- A) drop cycling labels ---------------------------------------------------------------------------
cycling_idx <- grepl("(?i)\\bCYCLING\\b|(?i)^CYCLING|(?i)CYCLING\\.", meta[[broad_col]])
cycling_labels_present <- sort(unique(meta[[broad_col]][cycling_idx]))
cycling_cells <- rownames(meta)[cycling_idx]
fwrite(data.frame(label = cycling_labels_present), file.path(readiness_dir, "cycling_labels_found.csv"))
message("Cycling labels detected: ", paste(cycling_labels_present, collapse = ", "), " | cells: ", length(cycling_cells))

## ---- B) decontX per sample (non-cycling cells) ---------------------------------------------------------
DefaultAssay(scatlas) <- "RNA"
sc_noCycl <- subset(scatlas, cells = setdiff(colnames(scatlas), cycling_cells))
meta_nc   <- sc_noCycl@meta.data
samples <- unique(meta_nc[[sample_col]])
contam_vec <- setNames(numeric(ncol(sc_noCycl)), colnames(sc_noCycl))
for (s in samples) {
  idx <- which(meta_nc[[sample_col]] == s)
  if (length(idx) < 50) { contam_vec[idx] <- NA_real_; next }          # too small to estimate robustly
  sc_s <- subset(sc_noCycl, cells = colnames(sc_noCycl)[idx])
  sce  <- as.SingleCellExperiment(sc_s, assay = "RNA")
  sce$coarse_lab <- sc_s@meta.data[, broad_col, drop = TRUE]
  sce_dc <- celda::decontX(sce, z = sce$coarse_lab, assayName = "counts", estimateDelta = TRUE, verbose = TRUE)
  contam_vec[idx] <- colData(sce_dc)$decontX_contamination
  rm(sc_s, sce, sce_dc); gc()
}
saveRDS(contam_vec, file.path(readiness_dir, "contamination_scores.rds"))
scatlas@meta.data$decontX_contam <- NA_real_
scatlas@meta.data[names(contam_vec), "decontX_contam"] <- contam_vec

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
             sample_id = scatlas@meta.data[cycling_cells, sample_col, drop = TRUE], broad_label = scatlas@meta.data[cycling_cells, broad_col, drop = TRUE]),
  df_contam %>% filter(decontX_outlier) %>% transmute(cell, reason = "ambient_outlier", contam, sample_id, broad_label))
drop_df <- drop_df[!duplicated(drop_df$cell), ]
keep_cells <- setdiff(colnames(scatlas), drop_df$cell)
fwrite(drop_df, file.path(readiness_dir, "cells_to_drop.csv")); writeLines(keep_cells, file.path(readiness_dir, "cells_keep.txt"))

## ---- E) gene whitelists (side files; the object keeps all genes) ----------------------------------------------
genes_all <- rownames(scatlas@assays$RNA)
is_mt <- grepl("^MT-", genes_all); is_rps <- grepl("^RPS\\d", genes_all); is_rpl <- grepl("^RPL\\d", genes_all); is_hb <- grepl("^HB[ABDEGQMZ]", genes_all)
is_ig <- grepl("^IGH[VDJC]", genes_all) | grepl("^IGK", genes_all) | grepl("^IGL", genes_all) |
  genes_all %in% c("IGKC", "IGLC1", "IGLC2", "IGLC3", "IGLC6", "IGHM", "IGHD", "IGHE", "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4")
exclude_global <- is_mt | is_rps | is_rpl | is_hb
fwrite(data.frame(gene = genes_all[!exclude_global]), file.path(readiness_dir, "whitelist_withIg.csv"))
fwrite(data.frame(gene = genes_all[!exclude_global & !is_ig]), file.path(readiness_dir, "whitelist_noIg.csv"))

## ---- F) filtered object ---------------------------------------------------------------------------------------
sc_ref_base <- subset(scatlas, cells = keep_cells)
sc_ref_base@meta.data$decontX_contam <- scatlas@meta.data[rownames(sc_ref_base@meta.data), "decontX_contam", drop = TRUE]
saveRDS(sc_ref_base, file = file.path(readiness_dir, "sc_ref_base.rds"))

## ---- G) summary ---------------------------------------------------------------------------------------------------
fwrite(as.data.frame(sort(table(scatlas@meta.data[[broad_col]]), decreasing = TRUE)), file.path(readiness_dir, "label_sizes_before.csv"))
fwrite(as.data.frame(sort(table(sc_ref_base@meta.data[[broad_col]]), decreasing = TRUE)), file.path(readiness_dir, "label_sizes_after.csv"))
write_json(list(input_cells = ncol(scatlas), cycling_removed = length(cycling_cells), ambient_removed = length(ambient_cells),
                kept = length(keep_cells), abs_cut = abs_contam_cut, tail_quantile = tail_q),
           file.path(readiness_dir, "00_readiness_summary.json"), auto_unbox = TRUE, pretty = TRUE)
message("Kept cells: ", length(keep_cells))
