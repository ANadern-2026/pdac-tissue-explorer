# ---- plain-language summaries ---------------------------------------------------------------------------------
# A p-value is NA when every value is identical (e.g. a gene not expressed in any sample): reported as no difference.
fold_txt <- function(log2diff) { f <- 2^abs(log2diff); if (f < 1.15) "" else paste0(" (about ", if (f >= 10) round(f) else formatC(f, format = "f", digits = 1), "-fold ", if (log2diff > 0) "higher" else "lower", ")") }
prose <- function(x) unname(CELLTYPE_PROSE[as.character(x)])

txt_cptac <- function(p, s, symbol, assay) {
  what <- if (assay == "protein") "Protein" else "RNA"
  if (is.null(p)) return(paste0(what, " (CPTAC): ", symbol, " is not in the ", tolower(what), " data."))
  if (!isTRUE(s$ok)) return(paste0(what, " (CPTAC): measured in only ", s$n, " matched tumour/normal pairs, too few to compare."))
  dir <- if (!isTRUE(s$p < 0.05)) "no clear difference between tumour and matched adjacent normal tissue" else if (s$median_diff > 0) "higher in tumour than in matched adjacent normal tissue" else "lower in tumour than in matched adjacent normal tissue"
  cov <- if (assay == "protein" && p$n_tumours_measured < 0.8 * p$n_tumours_total)
    paste0(" Note: quantified in only ", p$n_tumours_measured, " of ", p$n_tumours_total, " tumours, typically a sign of low abundance; the comparison covers that subset.") else ""
  paste0(what, " (CPTAC, ", s$n, " patients): ", dir, "; ", s$n_up, " of ", s$n, " patients higher in tumour, median change ",
         fmt_signed(s$median_diff), " log2 units", if (isTRUE(s$p < 0.05)) fold_txt(s$median_diff) else "", " (paired Wilcoxon p = ", fmt_p(s$p), ").", cov)
}

# tab: phospho_table() of site abundance (with or without CI); rel: the same relative to protein (NULL if no protein)
txt_phospho <- function(tab, rel, symbol, n_patients) {
  if (is.null(tab)) return(paste0("Phosphosites (CPTAC): no ", symbol, " phosphosites were quantified."))
  tt <- tab[!is.na(tab$p), ]
  head <- paste0("Phosphosites (CPTAC, ", n_patients, " patients): ", nrow(tab), " site", if (nrow(tab) != 1) "s", " quantified")
  if (!nrow(tt)) return(paste0(head, ", none in at least ", PHOSPHO_MIN_PAIRS, " matched tumour/normal pairs, too few to compare."))
  cnt <- function(t) c(up = sum(t$q < 0.05 & t$median_diff > 0, na.rm = TRUE), down = sum(t$q < 0.05 & t$median_diff < 0, na.rm = TRUE))
  a <- cnt(tt); big <- tt[tt$q < 0.05, ]; big <- big[which.max(abs(big$median_diff)), ]
  out <- paste0(head, "; ", nrow(tt), " compared in at least ", PHOSPHO_MIN_PAIRS, " matched pairs: ", a[["up"]], " higher and ", a[["down"]],
                " lower in tumour (paired Wilcoxon, false discovery rate < 5% within the gene)",
                if (nrow(big)) paste0("; largest change at ", big$label, " (median ", fmt_signed(big$median_diff), " log2 units)") else "", ".")
  if (is.null(rel)) return(paste0(out, " The protein itself was not quantified, so these changes cannot be separated from changes in protein abundance."))
  r <- cnt(rel[!is.na(rel$p), ])
  paste0(out, " Relative to the protein level (phosphorylation beyond the change in protein abundance): ", r[["up"]], " higher and ", r[["down"]], " lower.")
}

txt_bulk <- function(b, symbol) {
  if (is.null(b)) return(paste0("Bulk RNA (TCGA/GTEx): ", symbol, " is not in this table (not expressed at TPM >= 1 in at least 10% of any cohort)."))
  pd <- b$value[b$cohort == "PDAC (TCGA)"]; gx <- b$value[b$cohort == "Healthy pancreas (GTEx)"]; an <- b$value[b$cohort == "Adjacent normal (TCGA)"]
  t <- two_group(pd, gx)
  dir <- if (!isTRUE(t$p < 0.05) || abs(t$diff) < 0.2) "similar in PDAC and healthy pancreas" else if (t$diff > 0) "higher in PDAC than in healthy pancreas" else "lower in PDAC than in healthy pancreas"
  paste0("Bulk RNA (TCGA/GTEx): ", dir, ": median ", fmt_num(t$median_x), " vs ", fmt_num(t$median_y), " log2(TPM+1)",
         if (isTRUE(t$p < 0.05)) fold_txt(t$diff) else "", " (", t$n_x, " vs ", t$n_y, " samples). The ", length(an),
         " TCGA adjacent-normal samples have a median of ", fmt_num(stats::median(an)),
         "; differences against healthy pancreas partly reflect the much higher stromal content of tumours.")
}

txt_sc <- function(cs, symbol) {
  if (is.null(cs)) return(paste0("Single-cell (atlas): ", symbol, " is not in the single-cell data."))
  sm <- cs$summary[!is.na(cs$summary$share), ]
  if (max(sm$median_frac, na.rm = TRUE) < 0.05)
    return(paste0("Single-cell (atlas): rarely detected; fewer than 5% of cells express it in every cell type. Droplet single-cell RNA sequencing often misses low-abundance transcripts, so absence here is not proof of absence."))
  top <- sm[order(-sm$share), ][1:min(3, nrow(sm)), ]
  fr <- sm[which.max(sm$median_frac), ]
  paste0("Single-cell (atlas, treatment-naive primary tumours): most ", symbol, " transcripts are in ",
         paste0(prose(top$compartment), " (", fmt_pct(top$share), ")", collapse = ", "),
         ". The highest fraction of expressing cells is in ", prose(fr$compartment), " (median ", fmt_pct(fr$median_frac), " of cells).")
}

txt_state <- function(v, symbol) {
  if (is.null(v)) return(NULL)
  cs <- celltype_summary(v); if (is.null(cs)) return(NULL)
  topc <- as.character(cs$summary$compartment[which.max(cs$summary$share)])
  ps <- sc_aggregate(v[!is.na(v$state) & v$compartment == topc, ], c("sample", "state")); ps <- ps[ps$n_cells >= MIN_CELLS, ]
  a <- two_group(ps$avg_expr[ps$state == "Primary tumour"], ps$avg_expr[ps$state == "Adjacent normal"])
  m <- two_group(ps$avg_expr[ps$state == "Liver metastasis"], ps$avg_expr[ps$state == "Primary tumour"])
  word <- function(t, a, b) if (!isTRUE(t$ok)) paste0(a, " vs ", b, ": too few samples") else
    paste0(a, " vs ", b, ": ", if (!isTRUE(t$p < 0.05)) "no clear difference" else if (t$diff > 0) "higher" else "lower", " (median ", fmt_num(t$median_x), " vs ", fmt_num(t$median_y), ", p = ", fmt_p(t$p), ", ", t$n_x, " vs ", t$n_y, " samples)")
  paste0("Within ", prose(topc), " (the main ", symbol, "-carrying cells): ", word(a, "primary tumour", "adjacent normal"), "; ", word(m, "liver metastasis", "primary tumour"), ".")
}

# ---- links out ------------------------------------------------------------------------------------------------
# One entry per external resource: label, hint (tooltip) and URL. Direct pages use the Entrez or Ensembl ID when
# known (gene_xrefs.rds), otherwise a search for the symbol. Protein resources are left out for genes known to be
# non-coding.
gene_links <- function(gr) {
  sym <- sub(" \\(.*\\)$", "", gr$symbol)                          # duplicated symbols carry a " (Entrez ...)" suffix
  q <- utils::URLencode(sym, reserved = TRUE)
  ens <- if ("ensembl" %in% names(gr)) gr$ensembl else NA_character_
  type <- if ("gene_type" %in% names(gr)) gr$gene_type else NA_character_
  coding <- is.na(type) || type == "protein-coding"
  url <- function(k, id = NA_character_) { L <- LINKS_OUT[[k]]
    if (!is.na(id) && !is.null(L$id)) sprintf(L$id, utils::URLencode(id, reserved = TRUE)) else sprintf(L$search, q) }
  keep <- c(ncbi = TRUE, uniprot = coding, hpa = coding, opentargets = TRUE, depmap = TRUE)
  ids <- c(ncbi = gr$entrez, uniprot = NA, hpa = ens, opentargets = ens, depmap = NA)
  out <- lapply(names(LINKS_OUT)[keep[names(LINKS_OUT)]], function(k) list(key = k, label = LINKS_OUT[[k]]$label, hint = LINKS_OUT[[k]]$hint, url = url(k, ids[[k]])))
  stats::setNames(out, vapply(out, `[[`, "", "key"))
}

# ---- multi-gene tables ----------------------------------------------------------------------------------------
multi_dot_table <- function(d, genes) {
  rows <- lapply(genes, function(g) { v <- sc_values(d, g); if (is.null(v)) return(NULL); cs <- celltype_summary(v); if (is.null(cs)) return(NULL)
    s <- cs$summary; s$gene <- g; s })
  tab <- do.call(rbind, rows); if (is.null(tab)) return(NULL)
  tab <- tab[!is.na(tab$median_avg_expr), ]
  tab$z <- ave(tab$median_avg_expr, tab$gene, FUN = function(x) if (stats::sd(x) > 0) (x - mean(x)) / stats::sd(x) else 0 * x)
  tab$gene <- factor(tab$gene, levels = genes); tab
}

multi_forest_table <- function(d, genes) {
  rows <- list()
  for (g in genes) {
    for (a in c("protein", "rna")) { p <- cptac_pairs(d, g, a); s <- paired_stats(p)
      if (isTRUE(s$ok)) rows[[length(rows) + 1]] <- data.frame(gene = g, measure = if (a == "protein") "Protein, CPTAC pairs" else "RNA, CPTAC pairs", n = s$n, diff = s$median_diff, lo = s$ci_low, hi = s$ci_high, p = s$p) }
    b <- bulk_values(d, g)
    if (!is.null(b)) { t <- two_group(b$value[b$cohort == "PDAC (TCGA)"], b$value[b$cohort == "Healthy pancreas (GTEx)"])
      rows[[length(rows) + 1]] <- data.frame(gene = g, measure = "RNA, PDAC vs healthy (TCGA/GTEx)", n = t$n_x, diff = t$diff, lo = NA, hi = NA, p = t$p) }
  }
  tab <- do.call(rbind, rows); if (is.null(tab)) return(NULL)
  tab$gene <- factor(tab$gene, levels = genes)
  tab$measure <- factor(tab$measure, levels = c("Protein, CPTAC pairs", "RNA, CPTAC pairs", "RNA, PDAC vs healthy (TCGA/GTEx)"))
  tab
}

parse_gene_list <- function(txt, d) {
  toks <- unique(toupper(trimws(unlist(strsplit(txt, "[,;\\s]+", perl = TRUE))))); toks <- toks[nzchar(toks)]
  up <- toupper(d$genes$symbol); hit <- character(0); miss <- character(0)
  for (t in toks) {
    i <- match(t, up)
    if (is.na(i)) { j <- which(vapply(strsplit(toupper(d$genes$aliases), ", "), function(a) t %in% a, logical(1))); i <- if (length(j)) j[1] else NA }
    if (is.na(i)) miss <- c(miss, t) else hit <- c(hit, d$genes$symbol[i])
  }
  list(genes = unique(hit), missing = miss)
}

# ---- clinical context and co-expression -----------------------------------------------------------------------
txt_survival <- function(tab, symbol) {
  if (is.null(tab)) return(NULL)
  u <- tab[tab$model == "unadjusted" & !is.na(tab$hr), ]
  if (!nrow(u)) return(paste0("Survival: not enough data for ", symbol, "."))
  each <- paste0(u$label, " ", fmt_num(u$hr), " (95% CI ", fmt_num(u$ci_low), "-", fmt_num(u$ci_high), ", p = ", fmt_p(u$p), ")")
  sig <- u[u$p < 0.05, ]
  verdict <- if (!nrow(sig)) "no dataset shows a clear association (p < 0.05)" else
    paste0("p < 0.05 in ", paste0(sig$label, " (higher values with ", ifelse(sig$hr > 1, "shorter", "longer"), " survival)", collapse = "; "))
  paste0("Survival (exploratory; Cox model, hazard ratio per standard deviation, > 1 = higher values with shorter survival): ", paste(each, collapse = "; "),
         ". ", toupper(substr(verdict, 1, 1)), substr(verdict, 2, nchar(verdict)), ". Not corrected for testing many genes; see the Clinical context tab, including models adjusted for age and stage.")
}

txt_clin_assoc <- function(a, symbol, in_score = FALSE) {
  if (is.null(a) || !nrow(a)) return(NULL)
  part <- function(kind, hi, lo) {
    r <- a[[paste0("rho_", kind)]]; p <- a[[paste0("p_", kind)]]; ok <- !is.na(r)
    if (!any(ok)) return("too few tumours")
    each <- paste0(fmt_num(r[ok]), " (", a$dataset[ok], ", p = ", fmt_p(p[ok]), ")")
    s <- ok & p < 0.05; dir <- function(x) paste("higher in", if (x > 0) hi else lo)
    v <- if (!any(s)) "no clear association" else if (length(unique(sign(r[s]))) > 1) "inconsistent between datasets" else
      if (sum(s) >= 2 || sum(ok) == 1) dir(r[s][1]) else paste0("p < 0.05 only in ", a$dataset[s], " (", dir(r[s]), ")")
    paste0(v, "; rho = ", paste(each, collapse = ", "))
  }
  paste0("Tumour subtype and stage (Spearman correlation across tumours): with the basal-like minus classical score: ", part("moffitt", "basal-like tumours", "classical tumours"),
         if (in_score) paste0(" (", symbol, " is itself one of the score's genes)") else "", ". With stage: ", part("stage", "higher-stage tumours", "lower-stage tumours"), ".")
}

txt_coexp <- function(px, symbol, compartment, info) {
  if (is.null(px) || is.null(compartment)) return(NULL)
  x <- px[px$compartment == compartment & px$r > 0 & px$q < 0.05, ]
  head <- paste0("Co-expression (single-cell, ", prose(compartment), " of ", info$n_samples, " treatment-naive tumours, study-adjusted): ")
  if (!nrow(x)) return(paste0(head, "no gene correlates with ", symbol, " at FDR < 5%."))
  x <- utils::head(x, 5)
  paste0(head, "most correlated across tumours are ", paste0(x$partner, " (rho = ", fmt_num(x$r), ")", collapse = ", "),
         ". Co-expression across patients suggests shared regulation or a shared cell state, not a physical or functional interaction.")
}

txt_geneset <- function(x, h, ha, genes, dirs, dataset, endpoint, mode, n_all) {
  if (is.null(x) || !nrow(x)) return(paste0(dataset, ": no patients with values for all of ", paste(genes, collapse = ", "), " after the filters."))
  unit <- switch(mode, graded = "each additional gene in its chosen direction", all = "all genes in their chosen direction (vs fewer)", score = "each standard deviation of the combined score")
  eff <- function(h) if (is.na(h$hr)) "not estimable (too few events or no variation)" else
    paste0("hazard ratio ", fmt_num(h$hr), " (95% CI ", fmt_num(h$lo), "-", fmt_num(h$hi), ", p = ", fmt_p(h$p), ")")
  paste0(dataset, ", ", tolower(SURV_ENDPOINTS[[endpoint]]), ": ", nrow(x), " patients (", sum(x$event), " events)",
         if (nrow(x) < n_all) paste0(" of ", n_all, " with values for all genes") else "", ". For ", unit, ": ", eff(h),
         "; adjusted for age and stage: ", eff(ha), ". Hazard ratio > 1 = shorter survival. Exploratory; see the notes below.")
}
