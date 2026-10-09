# ---- statistics computed on the fly --------------------------------------------------------------------
boot_median_ci <- function(x, n = BOOT_N, seed = 1L) {
  if (length(x) < 3) return(c(NA_real_, NA_real_))
  set.seed(seed); b <- replicate(n, stats::median(sample(x, replace = TRUE)))
  unname(stats::quantile(b, c(0.025, 0.975)))
}

# paired tumour vs normal (CPTAC)
paired_stats <- function(p) {
  if (is.null(p) || nrow(p$pairs) < 5) return(list(n = if (is.null(p)) 0L else nrow(p$pairs), ok = FALSE))
  d <- p$pairs$diff; nz <- d[d != 0]
  exact <- length(nz) < 50 && !anyDuplicated(abs(nz))
  pv <- suppressWarnings(stats::wilcox.test(p$pairs$tumour, p$pairs$normal, paired = TRUE, exact = exact)$p.value)
  ci <- boot_median_ci(d)
  list(ok = TRUE, n = length(d), median_diff = stats::median(d), ci_low = ci[1], ci_high = ci[2], p = pv,
       n_up = sum(d > 0), n_down = sum(d < 0), median_normal = stats::median(p$pairs$normal), median_tumour = stats::median(p$pairs$tumour))
}

# unpaired two-group comparison on per-sample values
two_group <- function(x, y, min_n = 4L) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < min_n || length(y) < min_n) return(list(ok = FALSE, n_x = length(x), n_y = length(y)))
  list(ok = TRUE, n_x = length(x), n_y = length(y), median_x = stats::median(x), median_y = stats::median(y),
       diff = stats::median(x) - stats::median(y), p = suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)$p.value))
}

fmt_p <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 1e-4, formatC(p, format = "e", digits = 1), formatC(p, format = "g", digits = 2)))
fmt_num <- function(x, digits = 2) ifelse(is.na(x), "n/a", formatC(x, format = "f", digits = digits))
fmt_signed <- function(x, digits = 2) ifelse(is.na(x), "n/a", paste0(ifelse(x > 0, "+", ""), formatC(x, format = "f", digits = digits)))
fmt_pct <- function(x) ifelse(is.na(x), "n/a", paste0(round(100 * x), "%"))

# cell-type summary for one context (state): per-sample values -> medians, fraction expressing, transcript share
celltype_summary <- function(v, state = "Primary tumour", min_cells = MIN_CELLS) {
  v <- v[v$state %in% state, ]
  if (!nrow(v)) return(NULL)
  ps <- sc_aggregate(v, c("sample", "compartment")); ps <- ps[ps$n_cells >= min_cells, ]
  pooled <- sc_aggregate(transform(v, all = "all"), c("compartment"))
  pooled$transcripts <- pooled$mean_cp10k * pooled$n_cells
  pooled$share <- pooled$transcripts / sum(pooled$transcripts)
  pooled$cell_share <- pooled$n_cells / sum(pooled$n_cells)
  med <- stats::aggregate(cbind(avg_expr, frac) ~ compartment, ps, stats::median)
  ns <- stats::aggregate(n_cells ~ compartment, ps, length); names(ns)[2] <- "n_samples"
  out <- merge(merge(med, ns, by = "compartment", all = TRUE), pooled[, c("compartment", "n_cells", "share", "cell_share")], by = "compartment", all = TRUE)
  names(out)[names(out) == "avg_expr"] <- "median_avg_expr"; names(out)[names(out) == "frac"] <- "median_frac"
  out$compartment <- factor(out$compartment, levels = CELLTYPE_ORDER)
  list(per_sample = ps, summary = out[order(out$compartment), ])
}
