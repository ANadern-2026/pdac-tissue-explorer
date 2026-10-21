# ---- statistics computed on the fly --------------------------------------------------------------------
# Same resamples and result as replicate(n, median(sample(x, replace = TRUE))) after set.seed(seed) (v0.1), but
# vectorised: all resamples drawn at once (identical random stream), medians from column-wise sorted values.
boot_median_ci <- function(x, n = BOOT_N, seed = 1L) {
  if (length(x) < 3) return(c(NA_real_, NA_real_))
  set.seed(seed); k <- length(x)
  m <- matrix(sample(x, k * n, replace = TRUE), k)
  m <- matrix(m[order(col(m), m)], k)
  b <- (m[(k + 1L) %/% 2L, ] + m[k %/% 2L + 1L, ]) / 2
  unname(stats::quantile(b, c(0.025, 0.975)))
}

# paired tumour vs normal (CPTAC)
paired_stats <- function(p, min_n = 5L, ci = TRUE) {
  if (is.null(p) || nrow(p$pairs) < min_n) return(list(n = if (is.null(p)) 0L else nrow(p$pairs), ok = FALSE))
  d <- p$pairs$diff; nz <- d[d != 0]
  exact <- length(nz) < 50 && !anyDuplicated(abs(nz))
  pv <- suppressWarnings(stats::wilcox.test(p$pairs$tumour, p$pairs$normal, paired = TRUE, exact = exact)$p.value)
  ci <- if (ci) boot_median_ci(d) else c(NA_real_, NA_real_)
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

# Phosphosites of one gene: one row per site with the paired tumour vs normal comparison (sites quantified in at least
# min_pairs patients) and Benjamini-Hochberg adjusted p-values across the compared sites of this gene.
# relative = TRUE compares site minus protein level, i.e. phosphorylation beyond the change in protein abundance.
phospho_table <- function(ph, relative = FALSE, min_pairs = PHOSPHO_MIN_PAIRS, ci = TRUE) {
  if (is.null(ph)) return(NULL)
  rows <- lapply(seq_len(nrow(ph$sites)), function(k) {
    pr <- phospho_pairs(ph, k, relative); s <- paired_stats(list(pairs = pr), min_n = min_pairs, ci = ci)
    data.frame(n_pairs = nrow(pr), median_normal = if (s$ok) s$median_normal else NA_real_, median_tumour = if (s$ok) s$median_tumour else NA_real_,
               median_diff = if (s$ok) s$median_diff else NA_real_, ci_low = if (s$ok) s$ci_low else NA_real_, ci_high = if (s$ok) s$ci_high else NA_real_,
               p = if (s$ok) s$p else NA_real_, n_up = if (s$ok) s$n_up else NA_integer_, n_down = if (s$ok) s$n_down else NA_integer_)
  })
  tab <- data.frame(ph$sites[, c("label", "site_id", "residue", "position", "refseq", "source_symbol", "peptide")],
                    n_tumours = rowSums(is.finite(ph$tumor)), do.call(rbind, rows), stringsAsFactors = FALSE)
  tab$q <- stats::p.adjust(tab$p, "BH")      # NA p-values (sites not compared) are not counted
  tab$values <- if (relative) "relative to protein" else "site abundance"
  tab
}

# ---- co-expression within one cell type across samples (Co-expression tab; partners precomputed by data-raw/01b) ----
# One value per sample: log1p of the mean normalised expression of all cells of the cell type in that sample (as in the
# cell-type view). Samples: treatment-naive primary tumours (COEXP_STATE) with >= MIN_CELLS cells of the type, from studies
# contributing >= COEXP_MIN_STUDY such samples. Study effects are removed by centring each gene within study, and the
# correlation is Spearman's rho of the centred values; p-values use n - 1 - (number of studies) degrees of freedom.
coexp_units <- function(groups, compartment, state = COEXP_STATE, min_cells = MIN_CELLS, min_study = COEXP_MIN_STUDY) {
  rows <- which(groups$compartment %in% compartment & groups$state %in% state)
  n <- tapply(groups$n_cells[rows], groups$sample[rows], sum)
  st <- tapply(groups$study[rows], groups$sample[rows], function(s) s[1])
  ok <- sort(names(n)[n >= min_cells]); ns <- table(st[ok]); ok <- ok[st[ok] %in% names(ns)[ns >= min_study]]
  rows <- rows[groups$sample[rows] %in% ok]
  # weights: rows (fine groups) -> samples, so that values %*% W is the cell-weighted mean per sample
  W <- matrix(0, length(rows), length(ok)); W[cbind(seq_along(rows), match(groups$sample[rows], ok))] <- groups$n_cells[rows]
  W <- sweep(W, 2, colSums(W), "/")
  list(compartment = as.character(compartment), state = state, samples = ok, study = unname(st[ok]), n_cells = unname(n[ok]),
       n_studies = length(unique(st[ok])), rows = rows, W = W)
}
# per-sample mean CP10k (or fraction) of one or more genes: x = values over all fine groups (vector or genes x groups matrix)
coexp_mean <- function(x, u) { if (is.null(dim(x))) x <- matrix(x, 1); x[, u$rows, drop = FALSE] %*% u$W }
# centre each row within study, then rank (Spearman); returns a matrix of the same shape
coexp_ranks <- function(X, study) {
  for (s in unique(study)) { j <- which(study == s); X[, j] <- X[, j, drop = FALSE] - rowMeans(X[, j, drop = FALSE]) }
  X <- round(X, 10)                                   # equal values stay tied irrespective of floating-point noise
  t(apply(X, 1, rank, ties.method = "average")) |> matrix(nrow = nrow(X))
}
coexp_df <- function(n, n_studies) n - 1L - n_studies
coexp_p <- function(r, df) { r <- pmin(pmax(r, -1), 1); t <- r * sqrt(df / pmax(1 - r^2, 1e-300)); 2 * stats::pt(-abs(t), df) }
# rho between two genes in one cell type (both: values over all fine groups, e.g. sc_values()$mean_cp10k)
coexp_pair <- function(x, y, u) {
  R <- coexp_ranks(log1p(rbind(coexp_mean(x, u), coexp_mean(y, u))), u$study)
  r <- if (stats::sd(R[1, ]) > 0 && stats::sd(R[2, ]) > 0) stats::cor(R[1, ], R[2, ]) else NA_real_
  df <- coexp_df(length(u$samples), u$n_studies)
  list(r = r, p = if (is.na(r)) NA_real_ else coexp_p(r, df), n = length(u$samples), n_studies = u$n_studies, df = df)
}

# ---- survival (Clinical context tab) -------------------------------------------------------------------------------
# Exploratory associations of one gene with survival, as decided for the app: (1) Kaplan-Meier curves for a median split
# (no search for an optimal cut-point) with a log-rank test, and (2) Cox proportional-hazards models of the continuous value,
# per standard deviation within the cohort, unadjusted and adjusted for age and stage (I, II, III-IV).
SURV_ENDPOINTS <- c(os = "Overall survival", pfi = "Progression-free interval")
surv_setups <- function(cv) {     # dataset x endpoint combinations available for these clinical values
  if (is.null(cv)) return(NULL)
  s <- data.frame(dataset = rep(CLIN_DATASETS, each = 2), endpoint = rep(names(SURV_ENDPOINTS), 3), stringsAsFactors = FALSE)
  s <- s[!(s$dataset != "TCGA-PAAD RNA" & s$endpoint == "pfi") & s$dataset %in% unique(as.character(cv$dataset)), ]
  s$label <- paste0(s$dataset, ": ", tolower(SURV_ENDPOINTS[s$endpoint])); s$key <- paste(s$dataset, s$endpoint, sep = "|"); rownames(s) <- NULL; s
}
surv_data <- function(cv, dataset, endpoint) {
  x <- cv[cv$dataset == dataset, ]
  data.frame(patient = x$patient, value = x$value, time = x[[paste0(endpoint, "_time")]] / DAYS_PER_MONTH, event = x[[paste0(endpoint, "_event")]],
             age = x$age, stage = x$stage, stringsAsFactors = FALSE)[!is.na(x[[paste0(endpoint, "_time")]]) & !is.na(x[[paste0(endpoint, "_event")]]), ]
}
stage3 <- function(s) factor(ifelse(s %in% c("III", "IV"), "III-IV", s), levels = c("I", "II", "III-IV"))
cox_hr <- function(x, adjusted = FALSE) {
  na <- list(hr = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, n = nrow(x), events = sum(x$event))
  if (adjusted) { x <- x[!is.na(x$age) & !is.na(x$stage), ]; x$stage <- droplevels(stage3(x$stage)); na$n <- nrow(x); na$events <- sum(x$event) }
  if (sum(x$event) < SURV_MIN_EVENTS || !isTRUE(stats::sd(x$value) > 0)) return(na)
  x$z <- (x$value - mean(x$value)) / stats::sd(x$value)
  f <- if (adjusted && nlevels(x$stage) > 1) survival::Surv(time, event) ~ z + age + stage else if (adjusted) survival::Surv(time, event) ~ z + age else survival::Surv(time, event) ~ z
  fit <- tryCatch(suppressWarnings(survival::coxph(f, data = x)), error = function(e) NULL)
  if (is.null(fit) || !is.finite(stats::coef(fit)[["z"]])) return(na)
  b <- stats::coef(fit)[["z"]]; se <- sqrt(stats::vcov(fit)["z", "z"])
  list(hr = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = 2 * stats::pnorm(-abs(b / se)), n = nrow(x), events = sum(x$event))
}
# median split: "high" = above the median, "low" = at or below; NULL if a group is too small (many identical values)
km_split <- function(x) {
  g <- factor(ifelse(x$value > stats::median(x$value), "High", "Low"), levels = c("High", "Low"))
  if (min(table(g)) < SURV_MIN_GROUP || sum(x$event) < SURV_MIN_EVENTS) return(NULL)
  fit <- survival::survfit(survival::Surv(time, event) ~ g, data = data.frame(x, g = g))
  lr <- survival::survdiff(survival::Surv(time, event) ~ g, data = data.frame(x, g = g))
  p <- stats::pchisq(lr$chisq, df = 1, lower.tail = FALSE)
  med <- summary(fit)$table[, "median"]; names(med) <- sub("^g=", "", names(med))
  steps <- do.call(rbind, lapply(levels(g), function(k) {
    s <- fit[paste0("g=", k)]
    data.frame(group = k, time = c(0, s$time), surv = c(1, s$surv), n_risk = c(s$n, s$n.risk), n_event = c(0, s$n.event), n_censor = c(0, s$n.censor))
  }))
  ticks <- seq(0, max(x$time), by = 12)
  at_risk <- do.call(rbind, lapply(levels(g), function(k) data.frame(group = k, time = ticks, n = vapply(ticks, function(t) sum(x$time[g == k] >= t), integer(1)))))
  list(steps = steps, at_risk = at_risk, p = p, median = med, n = c(table(g)), events = c(tapply(x$event, g, sum)), cut = stats::median(x$value), group = g)
}
surv_table <- function(cv) {      # all dataset x endpoint combinations: unadjusted and adjusted HR per SD
  s <- surv_setups(cv); if (is.null(s) || !nrow(s)) return(NULL)
  do.call(rbind, lapply(seq_len(nrow(s)), function(i) {
    x <- surv_data(cv, s$dataset[i], s$endpoint[i])
    do.call(rbind, lapply(c(FALSE, TRUE), function(adj) { h <- cox_hr(x, adj)
      data.frame(dataset = s$dataset[i], endpoint = s$endpoint[i], label = s$label[i], model = if (adj) "adjusted for age and stage" else "unadjusted",
                 n = h$n, events = h$events, hr = h$hr, ci_low = h$lo, ci_high = h$hi, p = h$p, stringsAsFactors = FALSE) }))
  }))
}

# ---- stage and classical/basal-like score ----
# Spearman correlation with stage as an ordered variable (I < II < III < IV) and with the Moffitt basal-minus-classical score
spearman <- function(x, y, min_n = 10L) {
  ok <- is.finite(x) & is.finite(y); n <- sum(ok)
  if (n < min_n || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) return(list(rho = NA_real_, p = NA_real_, n = n))
  t <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  list(rho = unname(t$estimate), p = t$p.value, n = n)
}
STAGE_NUM <- c(I = 1, II = 2, III = 3, IV = 4)
clin_assoc <- function(cv) {      # per dataset: correlation with stage and with the Moffitt score
  if (is.null(cv)) return(NULL)
  do.call(rbind, lapply(split(cv, cv$dataset, drop = TRUE), function(x) {
    st <- spearman(x$value, unname(STAGE_NUM[x$stage])); mo <- spearman(x$value, x$moffitt)
    data.frame(dataset = x$dataset[1], n_stage = st$n, rho_stage = st$rho, p_stage = st$p, n_moffitt = mo$n, rho_moffitt = mo$rho, p_moffitt = mo$p)
  }))
}

# All genes of one cell type (data-raw/01b): mean_x100 / frac_x10000 = the integer genes x fine-groups matrices of sc_expr.rds.
# Genes detected in >= min_frac of the cells in at least half of the samples are correlated with each other (same arithmetic
# as coexp_pair()); for each gene the top_pos most positive and top_neg most negative partners are returned with rho
# (3 decimals) and the Benjamini-Hochberg FDR across all analysed genes. gene/partner = row indices of the input matrices.
coexp_top_pairs <- function(mean_x100, frac_x10000, u, min_frac = COEXP_MIN_FRAC, top_pos = COEXP_TOP_POS, top_neg = COEXP_TOP_NEG, block = 1500L) {
  fr <- (frac_x10000[, u$rows, drop = FALSE] / 1e4) %*% u$W                  # sc_values() divides first, as here
  keep <- unname(which(apply(fr, 1, stats::median) >= min_frac)); rm(fr)
  R <- coexp_ranks(log1p((mean_x100[keep, u$rows, drop = FALSE] / 100) %*% u$W), u$study)
  ok <- apply(R, 1, stats::sd) > 0; keep <- keep[ok]; R <- R[ok, , drop = FALSE]
  Z <- R - rowMeans(R); Z <- Z / sqrt(rowSums(Z^2))                         # unit-norm rows: Z %*% t(Z) = Pearson correlation of ranks
  m <- length(keep); df <- coexp_df(length(u$samples), u$n_studies); out <- vector("list", m)
  for (b in split(seq_len(m), ceiling(seq_len(m) / block))) {
    C <- tcrossprod(Z[b, , drop = FALSE], Z)
    for (k in seq_along(b)) {
      i <- b[k]; r <- C[k, ]; r[i] <- NA
      q <- stats::p.adjust(coexp_p(r, df), "BH")                               # the gene itself (NA) is not counted
      o <- order(-r, na.last = NA); pos <- utils::head(o, top_pos)
      neg <- utils::head(rev(o), top_neg); neg <- setdiff(neg[r[neg] < 0], pos)
      j <- c(pos, neg)
      out[[i]] <- data.frame(gene = keep[i], partner = keep[j], r = round(r[j], 3), q = signif(q[j], 3))
    }
  }
  list(pairs = if (m) do.call(rbind, out) else NULL, n_genes = m, df = df)
}
