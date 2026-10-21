# =====================================================================================================
# 01b_sc_coexpression.R
# Co-expression within each cell type across tumour samples, from the single-cell summaries written by 01.
# For every cell type: one value per sample = log1p of the mean normalised expression (CP10k) of all cells of that
# type in the sample; samples = treatment-naive primary tumours with >= 20 cells of the type, from studies that
# contribute >= 3 such samples. Each gene is centred within study (removes study/batch offsets), and genes are
# correlated by Spearman's rho of the centred values. Only genes detected in >= 2% of the type's cells in at least
# half of the samples are analysed. For each analysed gene the 50 most positively and 20 most negatively correlated
# genes are kept, with p-values (n - 1 - number of studies degrees of freedom) and Benjamini-Hochberg FDR across all
# genes analysed in that cell type. Functions and settings are shared with the app (app/R/00_config.R, 02_stats.R).
# Input : app/data/sc_groups.rds, app/data/sc_expr.rds (outputs of 01, also on Zenodo)
# Output: app/data/sc_coexp.rds
# Run from the repository root:  source("data-raw/01b_sc_coexpression.R")   Requires base R; ~15 min, ~4 GB RAM.
# =====================================================================================================

source(file.path("data-raw", "helpers.R")); log_file <- start_log("01b_sc_coexpression")
for (f in c("00_config.R", "01_data.R", "02_stats.R")) source(file.path("app", "R", f))
IN_DIR <- Sys.getenv("PDACX_DATA_DIR", file.path("app", "data"))      # (00_config.R's DATA_DIR is relative to app/)

src <- file.path(IN_DIR, c("sc_groups.rds", "sc_expr.rds"))
for (f in src) logc(log_file, "source ", basename(f), " MD5 ", unname(tools::md5sum(f)))
groups <- prepare_sc_groups(readRDS(src[1])); e <- readRDS(src[2])
stopifnot(identical(groups$group_id, seq_len(nrow(groups))), ncol(e$mean_cp10k_x100) == nrow(groups))
logc(log_file, "settings: state = ", COEXP_STATE, "; min cells per sample = ", MIN_CELLS, "; min samples per study = ", COEXP_MIN_STUDY,
     "; min samples = ", COEXP_MIN_SAMPLES, "; min median fraction expressing = ", COEXP_MIN_FRAC, "; partners kept = ", COEXP_TOP_POS, " + ", COEXP_TOP_NEG)

res <- list(); info <- list()
for (cc in CELLTYPE_ORDER) {
  t0 <- Sys.time()
  u <- coexp_units(groups, cc)
  if (length(u$samples) < COEXP_MIN_SAMPLES) { logc(log_file, cc, ": ", length(u$samples), " samples, skipped"); next }
  cp <- coexp_top_pairs(e$mean_cp10k_x100, e$frac_x10000, u)
  res[[cc]] <- cp$pairs
  info[[cc]] <- data.frame(compartment = cc, n_samples = length(u$samples), n_studies = u$n_studies, df = cp$df, n_genes = cp$n_genes,
                           studies = paste(sort(unique(u$study)), collapse = ", "))
  logc(log_file, cc, ": ", length(u$samples), " samples from ", u$n_studies, " studies, ", cp$n_genes, " genes analysed; ",
       sum(cp$pairs$r > 0 & cp$pairs$q < 0.05), " stored positive pairs at FDR < 5%; ", round(as.numeric(Sys.time() - t0, units = "secs")), " s")
  rm(cp); invisible(gc())
}

tab <- do.call(rbind, unname(Map(function(cc, d) data.frame(compartment = match(cc, CELLTYPE_ORDER), d), names(res), res)))
rownames(tab) <- NULL
OUT <- file.path(IN_DIR, "sc_coexp.rds")
saveRDS(list(pairs = tab, genes = e$genes, compartments = CELLTYPE_ORDER, info = do.call(rbind, unname(info)),
             settings = list(state = COEXP_STATE, min_cells = MIN_CELLS, min_study = COEXP_MIN_STUDY, min_samples = COEXP_MIN_SAMPLES,
                             min_frac = COEXP_MIN_FRAC, top_pos = COEXP_TOP_POS, top_neg = COEXP_TOP_NEG),
             note = "gene, partner = row indices into `genes` (the single-cell gene symbols); compartment = index into `compartments`; r = Spearman's rho of study-centred values (3 d.p.); q = BH FDR across the genes analysed in that cell type",
             source_md5 = unname(tools::md5sum(src)), built = format(Sys.time(), tz = "UTC")), OUT, compress = "xz")
logc(log_file, "pairs stored: ", nrow(tab), " | wrote ", OUT, " (", round(file.size(OUT) / 1e6, 1), " MB)")
end_log(log_file)
