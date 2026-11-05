# ---- model systems (models.sqlite, data-raw/06): summaries, figures and text --------------------------------------
# Design (user decisions 2026-11-04): cell lines and organoids serve different purposes and are shown separately; values are
# compared only within one technology (DepMap lines with DepMap lines of all lineages; GDC organoids with GDC organoids and their
# source tumours). Patient tissue appears only as background (where patients' cells rank a gene), never on one axis with models,
# because a cell-type average from single-cell data and a monoculture differ in technology and composition. Organoid availability
# (biobanks) is not used to suggest models.
MODEL_GROUPS <- c("Cell lines (DepMap)", "Organoids (DepMap)", "Organoids (HCMI)", "Organoids (Tiriac 2018)", "Normal organoids (Tiriac 2018)")
MODEL_COLS   <- stats::setNames(c("#eb6834", "#2a78d6", "#1baf7a", "#4a3aa7", "#e87ba4"), MODEL_GROUPS)   # validated order (dataviz validator, light mode)
ORGANOID_GROUPS <- MODEL_GROUPS[2:5]
ORIGIN_COLS  <- c("Primary tumour" = "#2a78d6", "Metastasis" = "#eb6834", "Not recorded" = "#b5b4ae")           # validated pair + neutral
OTHER_FILL   <- "#b5b4ae"
SIM_TUMOUR   <- "Ductal: tumour subclusters"
DEP_STRONG   <- -0.5     # gene effect below this: the knockout clearly slows growth (DepMap convention: -1 = median of common essential genes)

nth <- function(n) if (is.na(n)) "n/a" else paste0(n, if (n %% 100 %in% 11:13) "th" else c("th", "st", "nd", "rd", "th", "th", "th", "th", "th", "th")[n %% 10 + 1])
model_group <- function(m) factor(ifelse(m$type == "Normal organoid", "Normal organoids (Tiriac 2018)",
                                  ifelse(m$source == "DepMap", ifelse(m$type == "Cell line", "Cell lines (DepMap)", "Organoids (DepMap)"),
                                         ifelse(m$source == "HCMI", "Organoids (HCMI)", "Organoids (Tiriac 2018)"))), levels = MODEL_GROUPS)
model_caption <- function(...) src_caption("DepMap Public 26Q3 (Broad, CC BY 4.0); GDC ORGANOID-PANCREATIC (Tiriac et al. 2018) and HCMI", ...)
subtype_rel_label <- function(z, ref = z) { q <- stats::quantile(ref, c(1 / 3, 2 / 3), na.rm = TRUE)
  ifelse(is.na(z), NA, ifelse(z <= q[[1]], "more classical", ifelse(z >= q[[2]], "more basal-like", "intermediate"))) }
fmt_supplier <- function(s, c) ifelse(is.na(s) | !nzchar(s) | is.na(c) | !nzchar(c), "", paste0(s, " ", c))   # only with a catalogue number (a purchasable line)

# median (white diamond) and interquartile range (bar) per row of a horizontal strip plot (y = factor column "y")
iqr_marks <- function(x, v) {
  s <- do.call(rbind, lapply(split(x[[v]], x$y, drop = TRUE), function(z) data.frame(q1 = stats::quantile(z, 0.25), med = stats::median(z), q3 = stats::quantile(z, 0.75), n = length(z))))
  s$y <- factor(rownames(s), levels = levels(x$y)); s$text <- paste0("median ", fmt_num(s$med), "<br>middle half: ", fmt_num(s$q1), " to ", fmt_num(s$q3))
  list(geom_errorbar(data = s, aes(xmin = q1, xmax = q3, y = y), inherit.aes = FALSE, width = 0.35, colour = INK[["primary"]], linewidth = 0.5, orientation = "y"),
       geom_point(data = s, aes(med, y, text = text), inherit.aes = FALSE, shape = 23, size = 3.4, fill = "white", colour = INK[["primary"]], stroke = 0.9))
}
strip_rows <- function(x, grp, levels) { lab <- paste0(levels, "\n(n = ", as.integer(table(factor(grp, levels = levels))), ")")
  factor(lab[match(grp, levels)], levels = rev(lab[table(factor(grp, levels = levels)) > 0])) }
strip_points <- function(sizes, alphas) list(
  geom_jitter(aes(text = text, fill = grp, alpha = grp, size = grp), shape = 21, colour = "white", stroke = 0.25, height = 0.22, width = 0),
  scale_alpha_manual(values = alphas, guide = "none"), scale_size_manual(values = sizes, guide = "none"))

# ---- cell lines (DepMap: one pipeline for every lineage) ----
line_expr_stats <- function(lines) {
  if (is.null(lines)) return(NULL)
  x <- lines[!is.na(lines$log2tpm) & lines$type %in% "Cell Line", ]; a <- x$log2tpm[x$pdac]; b <- x$log2tpm[!x$pdac]
  if (length(a) < 3) return(NULL)
  list(n_pdac = length(a), n_other = length(b), med_pdac = stats::median(a), med_other = stats::median(b),
       p = tryCatch(suppressWarnings(stats::wilcox.test(a, b)$p.value), error = function(e) NA_real_), pct_other_below = mean(b < stats::median(a)))
}
plot_line_expression <- function(mv, symbol) {
  s <- line_expr_stats(mv$lines)
  if (is.null(s)) return(empty_plot(paste0(symbol, ": not in the DepMap RNA data")))
  x <- mv$lines[!is.na(mv$lines$log2tpm) & mv$lines$type %in% "Cell Line", ]
  lv <- c("PDAC cell lines", "Other cancer cell lines"); x$grp <- factor(ifelse(x$pdac, lv[1], lv[2]), levels = lv); x$y <- strip_rows(x, as.character(x$grp), lv)
  x$text <- paste0(x$name, " (", x$lineage, ")<br>log2(TPM+1): ", fmt_num(x$log2tpm))
  set.seed(7L)
  ggplot(x, aes(log2tpm, y)) + strip_points(c(2.6, 1.3), c(0.95, 0.4)) + iqr_marks(x, "log2tpm") +
    scale_fill_manual(values = c("PDAC cell lines" = MODEL_COLS[[1]], "Other cancer cell lines" = OTHER_FILL), guide = "none") +
    labs(x = "RNA, log2(TPM + 1)", y = NULL, title = paste0(symbol, " RNA in PDAC cell lines and other cancer cell lines"),
         subtitle = paste0("Same sequencing and processing for every line. Median ", fmt_num(s$med_pdac), " in PDAC lines vs ", fmt_num(s$med_other),
                           " in other cancer lines (p = ", fmt_p(s$p), "); white diamond = median, bar = middle half"),
         caption = src_caption("Data: DepMap Public 26Q3, OmicsExpressionTPMLogp1HumanProteinCodingGenes (Broad, CC BY 4.0)")) + theme_px()
}

dependency_stats <- function(cr) {
  if (is.null(cr)) return(NULL)
  x <- cr[!is.na(cr$effect), ]; a <- x$effect[x$pdac]; b <- x$effect[!x$pdac]
  if (length(a) < 3) return(NULL)
  p <- tryCatch(suppressWarnings(stats::wilcox.test(a, b)$p.value), error = function(e) NA_real_)
  list(n_pdac = length(a), n_other = length(b), med_pdac = stats::median(a), med_other = stats::median(b), p = p,
       n_dep = sum(a < DEP_STRONG), frac_dep_other = mean(b < DEP_STRONG))
}
plot_dependency <- function(cr, symbol) {
  s <- dependency_stats(cr)
  if (is.null(s)) return(empty_plot(paste0(symbol, ": not in the DepMap CRISPR screens")))
  x <- cr[!is.na(cr$effect), ]
  lv <- c("PDAC cell lines", "PDAC organoids (DepMap)", "Other cancer lines")
  x$grp <- factor(ifelse(!x$pdac, lv[3], ifelse(x$type %in% "Organoid", lv[2], lv[1])), levels = lv); x$y <- strip_rows(x, as.character(x$grp), lv)
  x$text <- paste0(x$name, " (", x$lineage, ")<br>gene effect: ", fmt_num(x$effect))
  set.seed(7L)
  ggplot(x, aes(effect, y)) +
    geom_vline(xintercept = 0, colour = INK[["secondary"]], linewidth = 0.4) + geom_vline(xintercept = -1, colour = "#e34948", linewidth = 0.5, linetype = "dashed") +
    strip_points(c(2.6, 2.6, 1.3), c(0.95, 0.95, 0.4)) + iqr_marks(x, "effect") +
    scale_fill_manual(values = c("PDAC cell lines" = MODEL_COLS[[1]], "PDAC organoids (DepMap)" = MODEL_COLS[[2]], "Other cancer lines" = OTHER_FILL), guide = "none") +
    labs(x = "Effect of knocking the gene out on growth (0 = none; red dashed line -1 = a typical essential gene)", y = NULL,
         title = paste0(symbol, ": do PDAC cells need it to grow?"),
         subtitle = paste0("Median effect in PDAC models ", fmt_signed(s$med_pdac), " vs ", fmt_signed(s$med_other), " in other cancer lines (p = ", fmt_p(s$p), "); ",
                           s$n_dep, " of ", s$n_pdac, " PDAC models below ", DEP_STRONG, ". White diamond = median, bar = middle half"),
         caption = src_caption("Data: DepMap Public 26Q3 CRISPRGeneEffect (Chronos), Broad Institute, CC BY 4.0")) + theme_px()
}

plot_line_subtype <- function(mv, symbol) {
  if (is.null(mv$rna)) return(empty_plot(paste0(symbol, ": not in the model RNA data")))
  r <- mv$rna[mv$rna$type == "Cell line" & !is.na(mv$rna$subtype_rel) & !is.na(mv$rna$log2tpm), ]
  if (nrow(r) < 5) return(empty_plot("Too few PDAC cell lines"))
  rr <- spearman(r$subtype_rel, r$log2tpm)
  r$text <- paste0(r$name, "<br>log2(TPM+1): ", fmt_num(r$log2tpm), "<br>subtype score (vs other PDAC lines): ", fmt_num(r$subtype_rel))
  ggplot(r, aes(subtype_rel, log2tpm)) + geom_vline(xintercept = 0, colour = INK[["grid"]], linewidth = 0.4) +
    geom_point(aes(text = text), shape = 21, fill = MODEL_COLS[[1]], colour = "white", size = 3.4, stroke = 0.4, alpha = 0.95) +
    labs(x = "More classical  <-  0  ->  more basal-like (relative to the other PDAC cell lines)", y = "RNA, log2(TPM + 1)",
         title = paste0(symbol, " along the classical / basal-like axis of PDAC cell lines"),
         subtitle = paste0("Each dot is one PDAC cell line; rho = ", fmt_num(rr$rho), " (p = ", fmt_p(rr$p), "). Score: mean z-score of the Moffitt basal-like minus classical genes across PDAC cell lines"),
         caption = src_caption("Data: DepMap Public 26Q3; genes: Moffitt et al., Nat Genet 2015")) + theme_px()
}

plot_model_protein <- function(pr, symbol) {
  if (is.null(pr) || sum(!is.na(pr$value[pr$pdac])) < 3) return(empty_plot(paste0(symbol, ": not quantified in enough PDAC lines by CCLE proteomics")))
  x <- pr[!is.na(pr$value), ]; lv <- c("PDAC cell lines", "Other cancer lines")
  x$grp <- factor(ifelse(x$pdac, lv[1], lv[2]), levels = lv); x$y <- strip_rows(x, as.character(x$grp), lv)
  x$text <- paste0(x$name, " (", x$lineage, ")<br>relative abundance: ", fmt_num(x$value))
  set.seed(7L)
  ggplot(x, aes(value, y)) + strip_points(c(2.6, 1.3), c(0.95, 0.4)) + iqr_marks(x, "value") +
    scale_fill_manual(values = c("PDAC cell lines" = MODEL_COLS[[1]], "Other cancer lines" = OTHER_FILL), guide = "none") +
    labs(x = "Relative protein abundance (log2, centred across lines)", y = NULL, title = paste0(symbol, " protein in cancer cell lines"),
         subtitle = "Mass-spectrometry proteomics of the Cancer Cell Line Encyclopedia; relative between lines. White diamond = median, bar = middle half",
         caption = src_caption("Data: Nusinow et al., Cell 2020, harmonised by DepMap (26Q1), CC BY 4.0")) + theme_px()
}

# ---- organoids ----
plot_organoid_expression <- function(mv, symbol) {
  if (is.null(mv$rna)) return(empty_plot(paste0(symbol, ": not in the model RNA data")))
  r <- mv$rna[mv$rna$type != "Cell line" & !is.na(mv$rna$log2tpm), ]; if (!nrow(r)) return(empty_plot(paste0(symbol, ": not measured in the organoids")))
  r$grp <- factor(as.character(model_group(r)), levels = ORGANOID_GROUPS); r$y <- strip_rows(r, as.character(r$grp), ORGANOID_GROUPS)
  r$text <- paste0(r$name, " (", r$grp, ")<br>log2(TPM+1): ", fmt_num(r$log2tpm))
  set.seed(7L)
  ggplot(r, aes(log2tpm, y)) + strip_points(stats::setNames(rep(2.6, 4), ORGANOID_GROUPS), stats::setNames(rep(0.9, 4), ORGANOID_GROUPS)) + iqr_marks(r, "log2tpm") +
    scale_fill_manual(values = MODEL_COLS[ORGANOID_GROUPS], guide = "none") +
    labs(x = "RNA, log2(TPM + 1)", y = NULL, title = paste0(symbol, " RNA in patient-derived organoids"),
         subtitle = "Compare organoids within a row: DepMap organoids and the GDC collections (HCMI, Tiriac) were sequenced and processed differently. Normal organoids: from non-tumour pancreas",
         caption = model_caption()) + theme_px()
}
pair_stats <- function(pp) {
  if (is.null(pp)) return(NULL)
  ok <- is.finite(pp$organoid_log2tpm) & is.finite(pp$tumour_log2tpm); if (sum(ok) < 10) return(NULL)
  r <- spearman(pp$tumour_log2tpm[ok], pp$organoid_log2tpm[ok])
  list(n = sum(ok), rho = r$rho, p = r$p, med_org = stats::median(pp$organoid_log2tpm[ok]), med_tum = stats::median(pp$tumour_log2tpm[ok]))
}
plot_hcmi_gene_pairs <- function(mv, symbol) {
  s <- pair_stats(mv$pairs); if (is.null(s)) return(empty_plot(paste0(symbol, ": not in the HCMI organoid / tumour pairs")))
  pp <- mv$pairs; pp$text <- paste0(pp$case, "<br>source tumour: ", fmt_num(pp$tumour_log2tpm), "<br>organoid: ", fmt_num(pp$organoid_log2tpm))
  ggplot(pp, aes(tumour_log2tpm, organoid_log2tpm)) +
    geom_point(aes(text = text), shape = 21, fill = MODEL_COLS[[3]], colour = "white", size = 3.4, stroke = 0.4, alpha = 0.95) +
    labs(x = "Source tumour (tissue), log2(TPM + 1)", y = "Organoid from the same patient, log2(TPM + 1)",
         title = paste0(symbol, ": do organoids keep the differences between patients?"),
         subtitle = paste0(s$n, " HCMI patients, organoid and own tumour sequenced and processed alike; rho = ", fmt_num(s$rho), " (p = ", fmt_p(s$p), "). ",
                           "Tumour tissue also contains stroma and immune cells, so genes of those cells are not expected to agree"),
         caption = src_caption("Data: Human Cancer Models Initiative (HCMI), GDC RNA-seq (STAR)")) + theme_px()
}

# ---- background: patients ----
plot_tissue_background <- function(mv, symbol, n_universe) {
  t <- mv$tissue; if (is.null(t)) return(empty_plot(paste0(symbol, ": no rank in the single-cell atlas (only protein-coding genes measured in every dataset)")))
  t <- t[t$group %in% CELLTYPE_ORDER & !is.na(t$pct), ]; t$group <- factor(t$group, levels = intersect(CELLTYPE_ORDER, unique(t$group)))
  sm <- stats::aggregate(pct ~ group, t, stats::median)
  t$grp <- ifelse(t$group == SIM_TUMOUR, "Tumour cells", "Other cells"); t$text <- paste0(t$group, "<br>sample ", t$sample, " (", t$n_cells, " cells)<br>rank: ", round(t$pct))
  sm$text <- paste0(sm$group, ": median rank ", round(sm$pct))
  set.seed(7L)
  ggplot(t, aes(pct, group)) +
    geom_jitter(aes(text = text, fill = grp), shape = 21, colour = "white", stroke = 0.25, height = 0.2, width = 0, size = 2, alpha = 0.85) +
    geom_point(data = sm, aes(pct, group, text = text), shape = 23, size = 3.4, fill = "white", colour = INK[["primary"]], stroke = 0.9) +
    scale_fill_manual(values = c("Tumour cells" = "#3d3c39", "Other cells" = OTHER_FILL), guide = "none") +
    scale_y_discrete(limits = rev(levels(t$group))) + scale_x_continuous(limits = c(-3, 103), breaks = seq(0, 100, 25), expand = c(0, 0)) +
    labs(x = "Rank among the sample's genes (100 = highest-expressed, 0 = not detected)", y = NULL,
         title = paste0("Background: where patients' cells rank ", symbol),
         subtitle = paste0("Single-cell atlas, treatment-naive primary tumours; one dot per tumour and cell type (rank among ", format(n_universe, big.mark = ","),
                           " protein-coding genes). Not comparable to the model values above: different technology and sample composition"),
         caption = src_caption("Data: Loveless et al., Clin Cancer Res 2025 single-cell atlas")) + theme_px()
}

# ---- tables and plain-language answers ----
model_table <- function(mv) {
  if (is.null(mv$rna)) return(NULL)
  r <- mv$rna; g <- as.character(model_group(r))
  out <- data.frame(Model = r$name, Group = g, From = ifelse(is.na(r$origin), "", r$origin), `log2(TPM+1)` = round(r$log2tpm, 2), check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(mv$crispr)) out$`Knockout effect` <- round(mv$crispr$effect[match(r$model, mv$crispr$model)], 2)
  if (!is.null(mv$protein)) out$`Protein (CCLE)` <- round(mv$protein$value[match(r$model, mv$protein$model)], 2)
  out$`Subtype (within its group)` <- r$subtype_rel; out$Supplier <- ifelse(r$type == "Cell line", fmt_supplier(r$supplier, r$catalog), "")
  out$RRID <- ifelse(is.na(r$rrid), "", r$rrid); out$`DepMap ID` <- ifelse(is.na(r$depmap_id), "", r$depmap_id)
  out[order(match(out$Group, MODEL_GROUPS), -out$`log2(TPM+1)`), ]
}

model_answers <- function(mv, symbol) {
  out <- list()
  s <- line_expr_stats(mv$lines)
  if (!is.null(s) && !is.null(mv$rna)) {
    r <- mv$rna[mv$rna$type == "Cell line" & !is.na(mv$rna$log2tpm), ]; r <- utils::head(r[order(-r$log2tpm), ], 5)
    lead <- if (s$med_pdac - s$med_other > 1 && isTRUE(s$p < 0.01)) "Higher than in most other cancer lines" else
            if (s$med_other - s$med_pdac > 1 && isTRUE(s$p < 0.01)) "Lower than in most other cancer lines" else
            if (s$med_pdac < 1) "Low in most PDAC cell lines" else "Similar to other cancer lines"
    out$lines <- list(q = "PDAC cell lines: highest expression", a = lead,
                      d = c(paste0(r$name, " (", fmt_num(r$log2tpm), if_supplier(r), ")")),
                      foot = paste0("log2(TPM+1); median ", fmt_num(s$med_pdac), " in ", s$n_pdac, " PDAC lines vs ", fmt_num(s$med_other), " in other cancer lines"))
  } else out$lines <- list(q = "PDAC cell lines: highest expression", a = "Not measured", d = paste0(symbol, " is not in the DepMap RNA data."))
  k <- dependency_stats(mv$crispr)
  out$dep <- if (is.null(k)) list(q = "Do PDAC cells need it to grow?", a = "Not tested", d = paste0(symbol, " is not in the DepMap knockout screens (or in too few PDAC models)."))
    else list(q = "Do PDAC cells need it to grow?",
              a = if (k$med_other < -0.8) "Yes, like most cancer cells" else
                  if (k$med_pdac < DEP_STRONG && k$med_pdac < k$med_other - 0.2 && isTRUE(k$p < 0.01)) "Yes, more than other cancers" else
                  if (k$n_dep / k$n_pdac >= 0.2) "Some PDAC models do" else "No clear effect",
              d = paste0("Knocking it out slowed growth clearly (effect below ", DEP_STRONG, ") in ", k$n_dep, " of ", k$n_pdac, " PDAC models and in ", fmt_pct(k$frac_dep_other), " of other cancer lines."))
  ps <- pair_stats(mv$pairs)
  out$org <- if (is.null(ps)) list(q = "Do organoids keep the patient's level?", a = "Not measured", d = paste0(symbol, " is not in the HCMI organoid / tumour pairs."))
    else list(q = "Do organoids keep the patient's level?",
              a = if (max(ps$med_org, ps$med_tum) < 1) "Too low to tell" else if (isTRUE(ps$rho >= 0.5 && ps$p < 0.05)) "Yes, largely" else
                  if (isTRUE(ps$rho >= 0.3 && ps$p < 0.05)) "Partly" else "No clear relation",
              d = paste0("Across ", ps$n, " HCMI patients, organoid vs own source tumour: rho = ", fmt_num(ps$rho),
                         " (p = ", fmt_p(ps$p), "). Patients with more ", symbol, " in the tumour ", if (isTRUE(ps$rho >= 0.3 && ps$p < 0.05)) "tend to" else "do not clearly", " have organoids with more."))
  out
}
if_supplier <- function(r) { s <- fmt_supplier(r$supplier, r$catalog); ifelse(nzchar(s), paste0("; ", s), "") }
answer_cards <- function(ans) div(class = "q-cards", lapply(ans, function(x) div(class = "q-card",
  div(class = "q-title", x$q), if (!is.null(x$a)) div(class = "q-verdict", x$a),
  if (length(x$d) > 1) tags$ol(class = "q-list", lapply(x$d, tags$li)) else p(class = "q-detail", x$d),
  if (!is.null(x$foot)) p(class = "q-foot", x$foot))))

# ---- overview page ----
model_overview <- function(info, groups = MODEL_GROUPS) { m <- info$models; m$group <- model_group(m); m[m$group %in% groups, ] }
ref_note <- function(info) paste0("Reference: ", info$source$references[["primary"]], " (single-cell RNA)")
overview_answers <- function(info) {
  m <- model_overview(info); cl <- m[m$type == "Cell line", ]; p <- info$hcmi_pairs
  list(
    list(q = "PDAC cell lines", a = paste0(nrow(cl), " lines with RNA"),
         d = paste0(sum(cl$has_crispr), " with genome-wide knockout screens, ", sum(cl$has_protein), " with proteomics; ", sum(nzchar(fmt_supplier(cl$supplier, cl$catalog)) & !is.na(cl$catalog)),
                    " with a supplier catalogue number in DepMap (shown in the tables).")),
    if (!is.null(p) && nrow(p)) list(q = "Do organoids keep their own patient's features?", a = paste0("Yes, ", sum(p$rank_own == 1), " of ", nrow(p)),
         d = paste0(sum(p$rank_own == 1), " of ", nrow(p), " HCMI organoids resemble their own source tumour more than any of the other ", nrow(p) - 1, " tumours (same sequencing pipeline).")),
    list(q = "Patient-derived organoids", a = paste0(sum(m$type == "Organoid"), " tumour organoids"),
         d = paste0("From three collections (DepMap ", sum(m$source == "DepMap" & m$type == "Organoid"), ", HCMI ", sum(m$source == "HCMI"), ", Tiriac et al. 2018 ",
                    sum(m$source == "Tiriac 2018" & m$type == "Organoid"), "), plus ", sum(m$type == "Normal organoid"), " organoids from non-tumour pancreas.")))
}
plot_overview_lines <- function(info) {
  m <- model_overview(info, MODEL_GROUPS[1]); m$from <- factor(ifelse(is.na(m$origin), "Not recorded", m$origin), levels = names(ORIGIN_COLS))
  m$text <- paste0(m$name, " (", m$from, ")<br>resemblance: ", fmt_num(m$sim_tumour), "<br>subtype vs other PDAC lines: ", fmt_num(m$subtype_rel),
                   ifelse(nzchar(fmt_supplier(m$supplier, m$catalog)), paste0("<br>", fmt_supplier(m$supplier, m$catalog)), ""))
  ggplot(m, aes(subtype_rel, sim_tumour)) + geom_vline(xintercept = 0, colour = INK[["grid"]], linewidth = 0.4) +
    geom_point(aes(text = text, fill = from), shape = 21, colour = "white", size = 3.6, stroke = 0.4, alpha = 0.95) +
    scale_fill_manual(values = ORIGIN_COLS, name = NULL, drop = TRUE) +
    labs(x = "More classical  <-  0  ->  more basal-like (relative to the other PDAC cell lines)", y = "Resemblance to patients' tumour cells",
         title = "PDAC cell lines", subtitle = paste0("Blue: derived from a primary tumour; orange: from a metastasis. Resemblance = Spearman rho with the reference over the 2,000 genes that best separate atlas cell types\n",
                                                     ref_note(info)), caption = model_caption()) + theme_px()
}
plot_overview_organoids <- function(info) {
  m <- model_overview(info, ORGANOID_GROUPS); m <- m[m$type != "Normal organoid", ]; m$group <- droplevels(m$group)
  m$text <- paste0(m$name, " (", m$group, ")<br>resemblance: ", fmt_num(m$sim_tumour), "<br>subtype vs other organoids of the collection: ", fmt_num(m$subtype_rel))
  ggplot(m, aes(subtype_rel, sim_tumour)) + geom_vline(xintercept = 0, colour = INK[["grid"]], linewidth = 0.4) +
    geom_point(aes(text = text, fill = group), shape = 21, colour = "white", size = 3.4, stroke = 0.4, alpha = 0.92) +
    facet_wrap(~group, nrow = 1) + scale_fill_manual(values = MODEL_COLS, guide = "none") +
    labs(x = "More classical  <-  0  ->  more basal-like (relative to the other organoids of the same collection)", y = "Resemblance to patients' tumour cells",
         title = "Patient-derived organoids, by collection",
         subtitle = paste0("Collections were sequenced and processed differently: compare organoids within a panel. ", ref_note(info)), caption = model_caption()) + theme_px()
}
plot_similarity_groups <- function(info, groups = MODEL_GROUPS) {
  m <- model_overview(info, groups); if (!nrow(m)) return(empty_plot("No models"))
  S <- info$similarity[m$model, , drop = FALSE]
  tab <- do.call(rbind, lapply(split(seq_len(nrow(m)), m$group, drop = TRUE), function(i) data.frame(group = m$group[i[1]], cell = colnames(S), r = apply(S[i, , drop = FALSE], 2, stats::median), n = length(i))))
  tab$cell <- factor(tab$cell, levels = intersect(CELLTYPE_ORDER, colnames(S))); tab$group <- factor(tab$group, levels = MODEL_GROUPS)
  tab$text <- paste0(tab$group, " (n = ", tab$n, ")<br>", tab$cell, ": median rho ", fmt_num(tab$r))
  ggplot(tab, aes(cell, group, fill = r, text = text)) + geom_tile(colour = "white", linewidth = 1.5) +
    geom_text(aes(label = formatC(r, format = "f", digits = 2)), size = 3, colour = ifelse(tab$r > stats::median(tab$r), "white", INK[["primary"]])) +
    scale_fill_gradientn(colours = SEQ_RAMP, name = "Median rho") + scale_y_discrete(limits = rev(levels(droplevels(tab$group)))) +
    labs(x = NULL, y = NULL, title = "Which patient cell type do the models resemble?",
         subtitle = "Median Spearman correlation of each model group with the average profile of each atlas cell type (treatment-naive primary tumours)",
         caption = model_caption("; atlas: Loveless et al. 2025")) + theme_px() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), panel.grid = element_blank())
}
plot_hcmi_pairs <- function(info) {
  p <- info$hcmi_pairs; if (is.null(p) || !nrow(p) || is.null(p$r_own)) return(empty_plot("No organoid - source tumour pairs"))
  p$own_first <- factor(ifelse(p$rank_own == 1, "Own tumour is the closest", "Another tumour is closer"), levels = c("Own tumour is the closest", "Another tumour is closer"))
  p$text <- paste0(p$case, "<br>own tumour: rho ", fmt_num(p$r_own), "<br>other tumours (median): ", fmt_num(p$r_other), "<br>own tumour ranks ", p$rank_own, " of ", nrow(p))
  ggplot(p, aes(r_other, r_own)) + geom_abline(slope = 1, intercept = 0, colour = INK[["muted"]], linewidth = 0.3, linetype = "dashed") +
    geom_point(aes(text = text, fill = own_first, shape = own_first), colour = "white", size = 3.4, stroke = 0.4, alpha = 0.92) +
    scale_fill_manual(values = c("Own tumour is the closest" = MODEL_COLS[[3]], "Another tumour is closer" = OTHER_FILL), name = NULL) +
    scale_shape_manual(values = c("Own tumour is the closest" = 21, "Another tumour is closer" = 22), name = NULL) +
    labs(x = "Correlation with the other patients' tumours (median)", y = "Correlation with its own tumour",
         title = "HCMI organoids keep features of their own patient's tumour",
         subtitle = paste0("Green circles: own tumour closest (", sum(p$rank_own == 1), " of ", nrow(p), " organoids); grey squares: another tumour is closer; dashed line = equal. ",
                           "Genes centred within organoids and within tumours, so the stroma of tissue samples does not dominate"),
         caption = src_caption("Data: HCMI organoids and their source tumours, GDC RNA-seq (STAR, same pipeline)")) + theme_px()
}
lines_table <- function(info) {
  m <- model_overview(info, MODEL_GROUPS[1])
  out <- data.frame(`Cell line` = m$name, From = ifelse(is.na(m$origin), "", paste0(m$origin, ifelse(is.na(m$collection_site), "", paste0(" (", gsub("_", " ", m$collection_site), ")")))),
                    Subtype = subtype_rel_label(m$subtype_rel), `Resemblance to tumour cells` = formatC(m$sim_tumour, format = "f", digits = 2),
                    `Knockout screen` = ifelse(m$has_crispr, "yes", ""), Proteomics = ifelse(m$has_protein, "yes", ""), Supplier = fmt_supplier(m$supplier, m$catalog),
                    RRID = ifelse(is.na(m$rrid), "", m$rrid), check.names = FALSE, stringsAsFactors = FALSE)
  out[order(-m$sim_tumour), ]
}
organoids_table <- function(info) {
  m <- model_overview(info, ORGANOID_GROUPS)
  out <- data.frame(Organoid = m$name, Collection = as.character(m$group), From = ifelse(is.na(m$origin), "", m$origin),
                    `Subtype (within collection)` = ifelse(m$type == "Normal organoid", "", unname(unlist(lapply(seq_len(nrow(m)), function(i) subtype_rel_label(m$subtype_rel[i], m$subtype_rel[m$group == m$group[i]]))))),
                    `Resemblance to tumour cells` = formatC(m$sim_tumour, format = "f", digits = 2),
                    `Own tumour (HCMI)` = if (is.null(m$own_tumour_rank)) "" else ifelse(is.na(m$own_tumour_rank), "", ifelse(m$own_tumour_rank == 1, "closest", paste0(m$own_tumour_rank, ". closest of ", m$n_source_tumours))),
                    `Knockout screen` = ifelse(m$has_crispr, "yes", ""), check.names = FALSE, stringsAsFactors = FALSE)
  out[order(match(out$Collection, MODEL_GROUPS), -m$sim_tumour), ]
}
fidelity_table <- function(info) {      # full numbers for the CSV download
  m <- model_overview(info); for (k in c("own_tumour_rank", "sim_livermet")) if (is.null(m[[k]])) m[[k]] <- NA
  data.frame(model = m$name, group = as.character(m$group), origin = m$origin, collection_site = m$collection_site, sex = m$sex, age = m$age,
             subtype_score_within_group = m$subtype_rel, resemblance_primary_tumour_cells = m$sim_tumour, resemblance_liver_metastasis_tumour_cells = m$sim_livermet,
             closest_atlas_cell_type = m$closest, resemblance_fibroblasts = info$similarity[m$model, "Fibroblasts"], own_tumour_rank_hcmi = m$own_tumour_rank, knockout_screen = m$has_crispr, proteomics = m$has_protein,
             supplier = m$supplier, catalog = m$catalog, depmap_id = m$depmap_id, hcmi_id = m$hcmi_id, rrid = m$rrid, stringsAsFactors = FALSE)
}

# ---- compare genes (each column one technology; medians of log2(TPM + 1)) ----
MULTI_MODEL_COLS <- c("PDAC cell lines", "Other cancer cell lines", ORGANOID_GROUPS)
multi_models_table <- function(d, genes) {
  rows <- lapply(genes, function(g) { mv <- model_values(d, g); if (is.null(mv$lines) && is.null(mv$rna)) return(NULL)
    v <- list()
    if (!is.null(mv$lines)) { l <- mv$lines[!is.na(mv$lines$log2tpm), ]
      v[["PDAC cell lines"]] <- l$log2tpm[l$pdac & l$type %in% "Cell Line"]; v[["Other cancer cell lines"]] <- l$log2tpm[!l$pdac & l$type %in% "Cell Line"] }
    if (!is.null(mv$rna)) { r <- mv$rna[!is.na(mv$rna$log2tpm), ]; gg <- as.character(model_group(r)); for (k in ORGANOID_GROUPS) v[[k]] <- r$log2tpm[gg == k] }
    v <- v[vapply(v, length, integer(1)) > 0]
    data.frame(gene = g, group = names(v), median = vapply(v, stats::median, numeric(1)), n = vapply(v, length, integer(1))) })
  tab <- do.call(rbind, rows); if (is.null(tab)) return(NULL)
  tab$group <- factor(tab$group, levels = MULTI_MODEL_COLS); tab$gene <- factor(tab$gene, levels = genes); droplevels(tab)
}
multi_dependency_table <- function(d, genes) {
  rows <- lapply(genes, function(g) { s <- dependency_stats(model_values(d, g)$crispr); if (is.null(s)) return(NULL)
    data.frame(gene = g, med_pdac = s$med_pdac, med_other = s$med_other, p = s$p, n_dep = s$n_dep, n_pdac = s$n_pdac) })
  tab <- do.call(rbind, rows); if (is.null(tab)) return(NULL); tab$gene <- factor(tab$gene, levels = genes); tab
}
plot_multi_models <- function(tab) {
  if (is.null(tab) || !nrow(tab)) return(empty_plot("No model data for these genes"))
  tab$text <- paste0(tab$gene, " | ", tab$group, "<br>median log2(TPM+1) ", fmt_num(tab$median), " (n = ", tab$n, ")")
  hi <- max(tab$median, na.rm = TRUE)
  ggplot(tab, aes(group, gene, fill = median, text = text)) + geom_tile(colour = "white", linewidth = 1.2) +
    geom_text(aes(label = formatC(median, format = "f", digits = 1)), size = 2.9, colour = ifelse(tab$median > 0.6 * hi, "white", INK[["primary"]])) +
    geom_vline(xintercept = 2.5, colour = INK[["primary"]], linewidth = 0.5) +
    scale_fill_gradientn(colours = SEQ_RAMP, name = "Median\nlog2(TPM+1)") + scale_y_discrete(limits = rev(levels(tab$gene))) +
    labs(x = NULL, y = NULL, title = "Genes in cell lines (left) and patient-derived organoids (right)",
         subtitle = "Median RNA, log2(TPM + 1). Left of the line: DepMap, one pipeline for all lines. Right: organoid collections; DepMap and GDC collections were processed differently",
         caption = model_caption()) + theme_px() + theme(axis.text.x = element_text(angle = 40, hjust = 1), panel.grid = element_blank())
}
plot_multi_dependency <- function(tab) {
  if (is.null(tab) || !nrow(tab)) return(empty_plot("None of these genes is in the DepMap CRISPR screens"))
  long <- rbind(data.frame(gene = tab$gene, who = "PDAC models", effect = tab$med_pdac, p = tab$p, n = tab$n_pdac, dep = tab$n_dep),
                data.frame(gene = tab$gene, who = "Other cancer lines", effect = tab$med_other, p = tab$p, n = NA, dep = NA))
  long$who <- factor(long$who, levels = c("PDAC models", "Other cancer lines"))
  long$text <- paste0(long$gene, " | ", long$who, "<br>median effect ", fmt_signed(long$effect), ifelse(long$who == "PDAC models", paste0("<br>", long$dep, " of ", long$n, " below ", DEP_STRONG, "<br>p = ", fmt_p(long$p)), ""))
  ggplot(long, aes(effect, gene)) + geom_vline(xintercept = 0, colour = INK[["secondary"]], linewidth = 0.4) + geom_vline(xintercept = -1, colour = "#e34948", linewidth = 0.5, linetype = "dashed") +
    geom_line(aes(group = gene), colour = INK[["grid"]], linewidth = 1.2) +
    geom_point(aes(fill = who, shape = who, text = text), colour = "white", size = 3.2, stroke = 0.4) +
    scale_fill_manual(values = c("PDAC models" = MODEL_COLS[[1]], "Other cancer lines" = OTHER_FILL), name = NULL) +
    scale_shape_manual(values = c("PDAC models" = 21, "Other cancer lines" = 22), name = NULL) + scale_y_discrete(limits = rev(levels(tab$gene))) +
    labs(x = "Median effect of knocking the gene out (0 = none; red dashed line -1 = a typical essential gene)", y = NULL, title = "Dependency in PDAC models vs other cancer lines",
         subtitle = "Orange circle: PDAC models; grey square: other cancer lines",
         caption = src_caption("Data: DepMap Public 26Q3 CRISPRGeneEffect (Chronos), CC BY 4.0")) + theme_px()
}

# one sentence for the Summary tab (NULL if the gene is not in the model data)
txt_models <- function(mv, symbol) {
  if (is.null(mv)) return(NULL)
  a <- model_answers(mv, symbol)
  paste0("Model systems: in PDAC cell lines ", tolower(a$lines$a), "; knockout effect: ", tolower(a$dep$a), "; organoids keep the patient's level: ", tolower(a$org$a),
         " (see the 'Model systems' tab).")
}
