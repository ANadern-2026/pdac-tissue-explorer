# ---- figures: every builder returns a ggplot (on screen via plotly with tooltips; downloads via ggsave) -----
suppressPackageStartupMessages({ library(ggplot2) })

theme_px <- function(base = 12) {
  theme_minimal(base_size = base) +
    theme(text = element_text(colour = INK[["primary"]]), axis.text = element_text(colour = INK[["secondary"]]),
          axis.title = element_text(colour = INK[["secondary"]]), panel.grid.major = element_line(colour = INK[["grid"]], linewidth = 0.3),
          panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = base + 1),
          plot.subtitle = element_text(colour = INK[["secondary"]]), plot.caption = element_text(colour = INK[["muted"]], size = base - 3, hjust = 0),
          strip.text = element_text(face = "bold", colour = INK[["primary"]], hjust = 0), legend.position = "bottom",
          plot.background = element_rect(fill = "white", colour = NA))
}
empty_plot <- function(msg) ggplot() + annotate("text", x = 0, y = 0, label = msg, colour = INK[["secondary"]], size = 4.5) + theme_void()
src_caption <- function(...) paste0(..., "  |  ", APP_TITLE, " v", APP_VERSION)

# CPTAC: one line per patient, adjacent normal -> tumour
# (also used for one phosphosite: p$assay "phospho" or "phospho_rel" and p$site = the site label)
plot_cptac_paired <- function(p, symbol, assay = "protein") {
  if (is.null(p) || nrow(p$pairs) < 3) return(empty_plot(paste0(symbol, ": not measured in enough matched patients (", assay, ")")))
  lab <- switch(p$assay, protein = "Protein abundance (log2, median-normalised)", rna = "RNA (log2 RSEM upper-quartile)",
                phospho = "Phosphosite (log2, median-normalised)", phospho_rel = "Phosphosite minus protein (log2)")
  ttl <- switch(p$assay, protein = paste0(symbol, " protein: tumour vs matched adjacent normal"), rna = paste0(symbol, " RNA: tumour vs matched adjacent normal"),
                paste0(symbol, " ", p$site, if (p$assay == "phospho_rel") ": relative to protein" else ": phosphosite"))
  long <- rbind(data.frame(p$pairs[, "case_id", drop = FALSE], tissue = "Adjacent normal", value = p$pairs$normal),
                data.frame(p$pairs[, "case_id", drop = FALSE], tissue = "Tumour", value = p$pairs$tumour))
  long$tissue <- factor(long$tissue, levels = c("Adjacent normal", "Tumour"))
  long$direction <- rep(ifelse(p$pairs$diff > 0, "higher in tumour", "lower in tumour"), 2)
  long$text <- paste0("Patient ", long$case_id, "<br>", long$tissue, ": ", fmt_num(long$value), "<br>change: ", fmt_signed(rep(p$pairs$diff, 2)))
  ggplot(long, aes(tissue, value, group = case_id)) +
    geom_line(aes(colour = direction), linewidth = 0.4, alpha = 0.55) +
    geom_point(aes(text = text), shape = 21, fill = INK[["secondary"]], colour = "white", size = 2.4, stroke = 0.5) +
    scale_colour_manual(values = c("higher in tumour" = "#eb6834", "lower in tumour" = "#2a78d6"), guide = "none") +
    labs(x = NULL, y = lab, title = ttl,
         subtitle = paste0(nrow(p$pairs), " patients with both samples measured; orange line = higher in tumour, blue = lower"),
         caption = src_caption("Data: Cao et al., Cell 2021 (CPTAC) via LinkedOmics")) + theme_px()
}

# phosphosites of one gene: median paired change per site with 95% bootstrap CI; colour = direction if FDR < 5%
plot_phospho_sites <- function(tab, symbol, relative = FALSE, has_sites = !is.null(tab)) {
  rel <- relative
  if (!has_sites) return(empty_plot(paste0(symbol, ": no phosphosites quantified in the CPTAC data")))
  if (is.null(tab)) return(empty_plot(paste0(symbol, ": protein not quantified, so values relative to protein are not available")))
  tt <- tab[!is.na(tab$p), ]
  if (!nrow(tt)) return(empty_plot(paste0(symbol, ": no phosphosite quantified in at least ", PHOSPHO_MIN_PAIRS, " matched patients")))
  shown <- if (nrow(tt) > PHOSPHO_MAX_PLOT) tt[order(tt$p, -tt$n_pairs)[seq_len(PHOSPHO_MAX_PLOT)], ] else tt
  shown <- shown[order(match(shown$label, tab$label)), ]                 # protein order: isoform, then position
  shown$class <- ifelse(shown$q < 0.05, ifelse(shown$median_diff > 0, "higher in tumour", "lower in tumour"), "no clear difference")
  shown$text <- paste0(symbol, " ", shown$label, "<br>median change: ", fmt_signed(shown$median_diff), " (95% CI ", fmt_num(shown$ci_low), " to ", fmt_num(shown$ci_high), ")",
                       "<br>", shown$n_up, " of ", shown$n_pairs, " patients higher in tumour<br>p = ", fmt_p(shown$p), ", FDR = ", fmt_p(shown$q),
                       ifelse(shown$source_symbol != symbol, paste0("<br>listed in CPTAC under ", shown$source_symbol), ""))
  ggplot(shown, aes(median_diff, label, text = text)) + geom_vline(xintercept = 0, colour = INK[["muted"]], linewidth = 0.3) +
    geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, colour = INK[["secondary"]], linewidth = 0.4, orientation = "y") +
    geom_point(aes(fill = class), shape = 21, size = 2.6, colour = INK[["primary"]], stroke = 0.3) +
    scale_fill_manual(values = c("higher in tumour" = "#eb6834", "lower in tumour" = "#2a78d6", "no clear difference" = "white"), guide = "none") +
    scale_y_discrete(limits = rev(shown$label)) +
    labs(x = if (rel) "Median change relative to protein, tumour minus normal (log2; bars = 95% bootstrap CI)" else "Median change, tumour minus normal (log2; bars = 95% bootstrap CI)",
         y = NULL, title = paste0(symbol, " phosphosites: tumour vs matched adjacent normal", if (rel) " (relative to protein level)" else ""),
         subtitle = paste0(nrow(tt), " sites with >= ", PHOSPHO_MIN_PAIRS, " patient pairs",
                           if (nrow(tt) > nrow(shown)) paste0(" (", nrow(shown), " with smallest p shown; all in the table)") else "",
                           "; orange = higher, blue = lower in tumour (FDR < 5%), white = no clear difference"),
         caption = src_caption("Data: Cao et al., Cell 2021 (CPTAC) via LinkedOmics; sites numbered on the RefSeq protein")) + theme_px()
}

plot_bulk <- function(b, symbol) {
  if (is.null(b)) return(empty_plot(paste0(symbol, " is not in the TCGA/GTEx expression table")))
  b$text <- paste0(b$cohort, "<br>", b$sample, "<br>log2(TPM+1): ", fmt_num(b$value))
  n <- table(b$cohort); lv <- paste0(levels(b$cohort), "\n(n = ", as.integer(n), ")"); b$x <- factor(lv[as.integer(b$cohort)], levels = lv)
  ggplot(b, aes(x, value)) +
    geom_jitter(aes(colour = cohort, text = text), width = 0.18, height = 0, size = 1.3, alpha = 0.55) +
    geom_boxplot(fill = NA, outlier.shape = NA, width = 0.45, colour = INK[["primary"]], linewidth = 0.35) +
    scale_colour_manual(values = BULK_COLS, guide = "none") +
    labs(x = NULL, y = "RNA, log2(TPM + 1)", title = paste0(symbol, " RNA: PDAC vs healthy and adjacent pancreas"),
         subtitle = "Bulk tissue; healthy pancreas is mostly acinar tissue, adjacent normal is fibrotic, stroma-rich pancreas",
         caption = src_caption("Data: TCGA-PAAD and GTEx, uniformly processed by recount3 (no cross-study batch correction)")) + theme_px()
}

plot_celltypes <- function(cs, symbol, state) {
  if (is.null(cs) || !nrow(cs$per_sample)) return(empty_plot(paste0("No single-cell data for ", symbol, " in this context")))
  ps <- cs$per_sample; sm <- cs$summary
  ps$text <- paste0(ps$compartment, "<br>sample ", ps$sample, " (", ps$n_cells, " cells)<br>average expression: ", fmt_num(ps$avg_expr),
                    "<br>cells expressing: ", fmt_pct(ps$frac))
  sm$text <- paste0(sm$compartment, "<br>median over ", sm$n_samples, " samples: ", fmt_num(sm$median_avg_expr),
                    "<br>median % expressing: ", fmt_pct(sm$median_frac), "<br>share of all ", symbol, " transcripts: ", fmt_pct(sm$share))
  ggplot(ps, aes(avg_expr, compartment)) +
    geom_jitter(aes(text = text), height = 0.18, width = 0, colour = ONE_SERIES, size = 1.2, alpha = 0.45) +
    geom_point(data = sm, aes(median_avg_expr, compartment, text = text), shape = 23, size = 3.2, fill = "white", colour = INK[["primary"]], stroke = 0.8) +
    scale_y_discrete(limits = rev(CELLTYPE_ORDER), drop = TRUE) +
    labs(x = "Expression per sample (log1p CP10k)", y = NULL,
         title = paste0(symbol, " by cell type: ", tolower(state), " samples"),
         subtitle = paste0("One dot per sample (>= ", MIN_CELLS, " cells of that type); white diamond = median"),
         caption = src_caption("Data: Loveless et al., Clin Cancer Res 2025 single-cell atlas (12 studies)")) + theme_px()
}

plot_share <- function(cs, symbol) {
  if (is.null(cs)) return(empty_plot("No data"))
  sm <- cs$summary[!is.na(cs$summary$share), ]
  sm$text <- paste0(sm$compartment, "<br>", fmt_pct(sm$share), " of ", symbol, " transcripts<br>", fmt_pct(sm$cell_share), " of all cells")
  ggplot(sm, aes(share, compartment, text = text)) +
    geom_col(fill = ONE_SERIES, width = 0.7) +
    geom_point(aes(cell_share), shape = 23, size = 2.2, fill = "white", colour = INK[["primary"]]) +
    scale_y_discrete(limits = rev(CELLTYPE_ORDER), drop = TRUE) + scale_x_continuous(labels = function(x) paste0(round(100 * x), "%")) +
    labs(x = "Share of transcripts (bar) and of cells (diamond)", y = NULL,
         title = "Which cells carry the transcripts", caption = src_caption("Pooled over all cells of the selected context")) + theme_px()
}

plot_frac <- function(cs, symbol) {
  if (is.null(cs)) return(empty_plot("No data"))
  sm <- cs$summary[!is.na(cs$summary$median_frac), ]
  sm$text <- paste0(sm$compartment, "<br>median ", fmt_pct(sm$median_frac), " of cells express ", symbol)
  ggplot(sm, aes(median_frac, compartment, text = text)) + geom_col(fill = ONE_SERIES, width = 0.7) +
    scale_y_discrete(limits = rev(CELLTYPE_ORDER), drop = TRUE) + scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, NA)) +
    labs(x = "Cells with detected transcript (median)", y = NULL, title = "How many cells express it",
         caption = src_caption("Single-cell detection underestimates low-expressed genes (dropout)")) + theme_px()
}

plot_states <- function(v, symbol, compartments) {
  if (is.null(v)) return(empty_plot("No single-cell data for this gene"))
  ps <- sc_aggregate(v[!is.na(v$state) & v$compartment %in% compartments, ], c("sample", "state", "compartment"))
  ps <- ps[ps$n_cells >= MIN_CELLS, ]
  if (!nrow(ps)) return(empty_plot("No samples with enough cells"))
  ps$state <- factor(ps$state, levels = STATE_LEVELS); ps$compartment <- factor(ps$compartment, levels = compartments)
  ps$text <- paste0(ps$compartment, " | ", ps$state, "<br>sample ", ps$sample, " (", ps$n_cells, " cells)<br>average expression: ", fmt_num(ps$avg_expr),
                    "<br>cells expressing: ", fmt_pct(ps$frac))
  ggplot(ps, aes(state, avg_expr)) +
    geom_jitter(aes(colour = state, text = text), width = 0.2, height = 0, size = 1.3, alpha = 0.6) +
    geom_boxplot(fill = NA, outlier.shape = NA, width = 0.5, colour = INK[["primary"]], linewidth = 0.3) +
    facet_wrap(~compartment, scales = "free_y", ncol = min(3, length(compartments))) +
    scale_colour_manual(values = STATE_COLS, guide = "none", drop = FALSE) +
    scale_x_discrete(labels = STATE_SHORT) +
    labs(x = NULL, y = "Expression per sample (log1p CP10k)", title = paste0(symbol, " across tissue states, by cell type"),
         subtitle = paste0("One dot per sample with >= ", MIN_CELLS, " cells of that type; Liver met = treatment-naive liver metastases; Treated = primary tumours after therapy"),
         caption = src_caption("Data: Loveless et al., Clin Cancer Res 2025 single-cell atlas")) + theme_px() +
    theme(axis.text.x = element_text(size = 9, angle = 30, hjust = 1))
}

plot_fibro <- function(v, symbol) {
  if (is.null(v)) return(empty_plot("No single-cell data for this gene"))
  ps <- sc_aggregate(v[v$state %in% "Primary tumour" & !is.na(v$fibro_subtype), ], c("sample", "fibro_subtype"))
  ps <- ps[ps$n_cells >= MIN_CELLS, ]
  if (!nrow(ps)) return(empty_plot("No fibroblast samples with enough cells"))
  ps$class <- ifelse(grepl("^myFb", ps$fibro_subtype), "Myofibroblasts (myFb)", "Inflammatory fibroblasts (iFb)")
  ps$label <- factor(ifelse(ps$fibro_subtype %in% names(FIBRO_NOTE), paste0(ps$fibro_subtype, "\n(", FIBRO_NOTE[as.character(ps$fibro_subtype)], ")"), as.character(ps$fibro_subtype)),
                     levels = ifelse(FIBRO_ORDER %in% names(FIBRO_NOTE), paste0(FIBRO_ORDER, "\n(", FIBRO_NOTE[FIBRO_ORDER], ")"), FIBRO_ORDER))
  ps$text <- paste0(ps$fibro_subtype, " | sample ", ps$sample, " (", ps$n_cells, " cells)<br>average expression: ", fmt_num(ps$avg_expr), "<br>cells expressing: ", fmt_pct(ps$frac))
  ggplot(ps, aes(label, avg_expr)) +
    geom_jitter(aes(colour = class, text = text), width = 0.2, height = 0, size = 1.3, alpha = 0.6) +
    geom_boxplot(fill = NA, outlier.shape = NA, width = 0.5, colour = INK[["primary"]], linewidth = 0.3) +
    scale_colour_manual(values = c("Myofibroblasts (myFb)" = "#2a78d6", "Inflammatory fibroblasts (iFb)" = "#eb6834"), guide = "none") +
    labs(x = NULL, y = "Expression per sample (log1p CP10k)", title = paste0(symbol, " in fibroblast subtypes (treatment-naive primary tumours)"),
         subtitle = "Blue = myofibroblast subtypes (myFb), orange = inflammatory fibroblast subtypes (iFb); one dot per sample",
         caption = src_caption("Subtype labels as published with the atlas (Loveless et al. 2025)")) + theme_px()
}

plot_ductal <- function(v, symbol) {
  if (is.null(v)) return(empty_plot("No single-cell data for this gene"))
  ps <- sc_aggregate(v[v$state %in% c("Primary tumour", "Liver metastasis") & !is.na(v$ductal_group), ], c("sample", "state", "ductal_group"))
  ps <- ps[ps$n_cells >= MIN_CELLS & ps$ductal_group != "Secretory", ]
  if (!nrow(ps)) return(empty_plot("No tumour samples with enough cells"))
  ps$state <- factor(ps$state, levels = c("Primary tumour", "Liver metastasis"))
  ps$text <- paste0(ps$ductal_group, " | ", ps$state, "<br>sample ", ps$sample, " (", ps$n_cells, " cells)<br>average expression: ", fmt_num(ps$avg_expr), "<br>cells expressing: ", fmt_pct(ps$frac))
  ggplot(ps, aes(ductal_group, avg_expr)) +
    geom_jitter(aes(colour = state, text = text), width = 0.2, height = 0, size = 1.3, alpha = 0.6) +
    geom_boxplot(fill = NA, outlier.shape = NA, width = 0.5, colour = INK[["primary"]], linewidth = 0.3) +
    facet_wrap(~state) + scale_colour_manual(values = STATE_COLS, guide = "none") +
    labs(x = NULL, y = "Expression per sample (log1p CP10k)", title = paste0(symbol, " in tumour-cell states (treatment-naive)"),
         subtitle = "Groups of the atlas's ductal tumour subclusters: classical-like 1/2/3/9, hybrid 5/6, basal-like 4/7",
         caption = src_caption("Data: Loveless et al., Clin Cancer Res 2025 single-cell atlas")) + theme_px()
}

# multi-gene dot heatmap: colour = z-score of the per-cell-type median across cell types, size = % expressing
plot_multi_dot <- function(tab) {
  if (is.null(tab) || !nrow(tab)) return(empty_plot("No single-cell data for these genes"))
  tab$text <- paste0(tab$gene, " | ", tab$compartment, "<br>median average expression: ", fmt_num(tab$median_avg_expr), "<br>median % expressing: ", fmt_pct(tab$median_frac))
  ggplot(tab, aes(compartment, gene, text = text)) +
    geom_point(aes(size = pmax(median_frac, 0.01), fill = z), shape = 21, colour = INK[["secondary"]], stroke = 0.25) +
    scale_fill_gradient2(low = DIV_COLS[["low"]], mid = DIV_COLS[["mid"]], high = DIV_COLS[["high"]], midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish,
                         name = "Relative expression\n(z-score across cell types)") +
    scale_size_area(max_size = 7, name = "Cells expressing", labels = function(x) paste0(round(100 * x), "%"), breaks = c(0.1, 0.25, 0.5, 0.75)) +
    scale_x_discrete(limits = CELLTYPE_ORDER, drop = TRUE) + scale_y_discrete(limits = rev(unique(tab$gene))) +
    labs(x = NULL, y = NULL, title = "Genes by cell type (treatment-naive primary tumours)",
         subtitle = "Colour = expression relative to the gene's own average across cell types; dot size = % of cells expressing",
         caption = src_caption("Data: Loveless et al., Clin Cancer Res 2025 single-cell atlas")) + theme_px() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.box = "vertical")
}

plot_multi_forest <- function(tab) {
  if (is.null(tab) || !nrow(tab)) return(empty_plot("No tumour vs normal data for these genes"))
  tab$text <- paste0(tab$gene, " | ", tab$measure, "<br>median change: ", fmt_signed(tab$diff), "<br>n = ", tab$n, "<br>p = ", fmt_p(tab$p))
  tab$sig <- ifelse(!is.na(tab$p) & tab$p < 0.05, "p < 0.05", "p >= 0.05")
  ggplot(tab, aes(diff, gene, text = text)) + geom_vline(xintercept = 0, colour = INK[["muted"]], linewidth = 0.3) +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, colour = INK[["secondary"]], linewidth = 0.4, orientation = "y") +
    geom_point(aes(shape = sig), size = 2.4, colour = ONE_SERIES, fill = "white") + scale_shape_manual(values = c("p < 0.05" = 16, "p >= 0.05" = 21), guide = "none") +
    facet_wrap(~measure, scales = "free_x") + scale_y_discrete(limits = rev(unique(tab$gene))) +
    labs(x = "Median change, tumour minus normal (log2 units; bars = 95% bootstrap CI)", y = NULL, title = "Tumour vs normal tissue",
         subtitle = "Filled dot = p < 0.05, open dot = p >= 0.05 (unadjusted)",
         caption = src_caption("CPTAC: matched pairs; TCGA vs GTEx: PDAC vs healthy pancreas (unpaired, interval not shown)")) + theme_px()
}

plot_umap_static <- function(u, symbol, values) {     # values: umap_values() for the gene, NULL if not in the atlas
  if (is.null(u)) return(empty_plot("Atlas map not available"))
  df <- data.frame(u$coords, value = if (is.null(values)) 0 else values)
  df <- df[order(df$value), ]
  ggplot(df, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = value), size = 0.25) +
    scale_colour_gradientn(colours = c(INK[["zero"]], SEQ_RAMP[-1]), name = "log1p(CP10k)") + coord_equal() +
    labs(title = paste0(symbol, " on the atlas map"), caption = src_caption("Subsample of atlas cells; coordinates from the published atlas")) + theme_px() +
    theme(axis.text = element_blank(), panel.grid = element_blank())
}

# static export: wrap long titles/subtitles to the page width
for_export <- function(p, width_in = 8) {
  w <- function(x, n) if (is.null(x)) NULL else paste(strwrap(x, n), collapse = "\n")
  p$labels$subtitle <- w(p$labels$subtitle, round(width_in * 10.5)); p$labels$title <- w(p$labels$title, round(width_in * 9))
  p$labels$caption <- w(p$labels$caption, round(width_in * 14)); p
}
