# Survival, stage/subtype association and single-cell co-expression functions (S3), on synthetic data with known answers.
root <- normalizePath(file.path(test_path(), "..", ".."))
for (f in sort(list.files(file.path(root, "app", "R"), full.names = TRUE))) source(f)
suppressPackageStartupMessages(library(survival))

set.seed(5)
n <- 120
sx <- data.frame(patient = paste0("P", 1:n), value = rnorm(n, 8, 2), age = round(runif(n, 40, 85)), stage = sample(c("I", "II", "III", "IV"), n, TRUE, c(.2, .5, .2, .1)))
sx$time <- rexp(n, 0.03 * exp(0.4 * (sx$value - 8) / 2)); sx$event <- as.integer(runif(n) < 0.7)

test_that("Cox hazard ratio per SD equals coxph on the standardised value", {
  h <- cox_hr(sx); fit <- coxph(Surv(time, event) ~ scale(value)[, 1], sx)
  expect_equal(h$hr, unname(exp(coef(fit)))); expect_equal(h$p, unname(summary(fit)$coefficients[, "Pr(>|z|)"]))
  expect_equal(c(h$lo, h$hi), unname(exp(confint(fit))[1, ]), tolerance = 1e-3)
  expect_gt(h$hr, 1)
  ha <- cox_hr(sx, adjusted = TRUE)
  fa <- coxph(Surv(time, event) ~ scale(value)[, 1] + age + stage3(stage), sx)
  expect_equal(ha$hr, unname(exp(coef(fa))[1]))
  expect_true(is.na(cox_hr(transform(sx, value = 1))$hr))                     # no variation
  expect_true(is.na(cox_hr(transform(sx, event = c(rep(1L, 5), rep(0L, n - 5))))$hr))   # too few events
})

test_that("median split: groups, log-rank test and numbers at risk", {
  k <- km_split(sx)
  g <- factor(ifelse(sx$value > median(sx$value), "High", "Low"), levels = c("High", "Low"))
  expect_identical(k$group, g); expect_equal(unname(k$n), as.integer(table(g)))
  expect_equal(k$p, pchisq(survdiff(Surv(time, event) ~ g, sx)$chisq, 1, lower.tail = FALSE))
  expect_equal(k$at_risk$n[k$at_risk$group == "High" & k$at_risk$time == 12], sum(sx$time[g == "High"] >= 12))
  expect_equal(k$steps$surv[k$steps$time == 0], c(1, 1))
  expect_null(km_split(transform(sx, value = c(rep(0, n - 5), 1:5))))          # mostly identical values: no split
})

test_that("Spearman associations with stage and the subtype score", {
  cv <- data.frame(dataset = factor("CPTAC RNA", levels = CLIN_DATASETS), patient = sx$patient, value = sx$value, stage = sx$stage,
                   moffitt = sx$value + rnorm(n))
  a <- clin_assoc(cv)
  expect_equal(a$rho_moffitt, unname(cor(cv$value, cv$moffitt, method = "spearman")))
  expect_equal(a$rho_stage, unname(cor(cv$value, STAGE_NUM[cv$stage], method = "spearman")))
  expect_match(txt_clin_assoc(a, "G"), "higher in basal-like tumours")
  expect_true(is.na(spearman(1:5, 1:5)$rho))                                   # fewer than 10 values
})

test_that("survival summary text reports every dataset and the direction", {
  tab <- data.frame(label = c("TCGA-PAAD RNA: overall survival", "CPTAC RNA: overall survival"), model = "unadjusted", hr = c(1.5, 0.9),
                    ci_low = c(1.1, 0.7), ci_high = c(2, 1.2), p = c(0.01, 0.4))
  t <- txt_survival(tab, "G")
  expect_match(t, "TCGA-PAAD RNA: overall survival 1.50"); expect_match(t, "shorter survival\\)"); expect_match(t, "CPTAC RNA: overall survival 0.90")
  expect_match(txt_survival(tab[2, ], "G"), "No dataset shows a clear association")
})

# single-cell co-expression: 3 studies with strong study offsets; two genes share a per-sample factor
mk_sc <- function() {
  set.seed(9)
  smp <- paste0("s", 1:30); study <- rep(c("A", "B", "C"), each = 10)
  g <- expand.grid(sample = smp, fine_group = c("myFb 1", "iFb 2"), stringsAsFactors = FALSE)
  g$group_id <- seq_len(nrow(g)); g$n_cells <- sample(15:60, nrow(g), TRUE); g$cell_type <- "FIBROBLASTS"
  g$study <- study[match(g$sample, smp)]; g$disease_state <- "Primary tumor"; g$treatment <- "Treatment-naive"; g$met_site <- NA
  g$n_cells[g$sample == "s30"] <- 5                                            # s30: too few cells
  g <- prepare_sc_groups(g)
  f <- rnorm(30); off <- c(A = 0, B = 3, C = -2)[study]
  base <- function(sig) pmax(0, exp(1 + sig + off + rnorm(30, 0, 0.3)) - 1)
  genes <- c("G1", "G2", "G3", "RARE")
  m <- rbind(G1 = base(f), G2 = base(f), G3 = base(-f), RARE = 0)[, match(g$sample, smp)] * runif(nrow(g), 0.9, 1.1)
  m["RARE", 1:3] <- 50
  mean_x100 <- matrix(as.integer(round(m * 100)), nrow(m), dimnames = list(genes, NULL))
  frac <- matrix(as.integer(ifelse(mean_x100 > 0, 3000L, 0L)), nrow(m), dimnames = list(genes, NULL))
  list(groups = g, mean = mean_x100, frac = frac)
}

test_that("co-expression: units, study adjustment and precomputed partners agree with the pairwise statistic", {
  sc <- mk_sc(); u <- coexp_units(sc$groups, "Fibroblasts")
  expect_identical(u$samples, sort(setdiff(paste0("s", 1:30), "s30"))); expect_equal(u$n_studies, 3)
  expect_equal(c(coexp_mean(sc$mean[1, ] / 100, u))[1], with(sc$groups[sc$groups$sample == u$samples[1], ], sum(sc$mean[1, group_id] / 100 * n_cells) / sum(n_cells)))
  p12 <- coexp_pair(sc$mean["G1", ] / 100, sc$mean["G2", ] / 100, u); p13 <- coexp_pair(sc$mean["G1", ] / 100, sc$mean["G3", ] / 100, u)
  expect_gt(p12$r, 0.8); expect_lt(p13$r, -0.8); expect_equal(p12$df, 29 - 1 - 3)
  # study offsets do not change the study-adjusted correlation
  x <- log1p(c(coexp_mean(sc$mean["G1", ] / 100, u))); y <- log1p(c(coexp_mean(sc$mean["G2", ] / 100, u)))
  cen <- function(v) v - ave(v, u$study)
  expect_equal(p12$r, cor(rank(cen(x)), rank(cen(y))))
  tp <- coexp_top_pairs(sc$mean, sc$frac, u, top_pos = 2, top_neg = 2)
  expect_equal(tp$n_genes, 3)                                                  # RARE: detected in too few samples
  r12 <- tp$pairs$r[tp$pairs$gene == 1 & tp$pairs$partner == 2]; expect_equal(r12, round(p12$r, 3))
  expect_equal(tp$pairs$r[tp$pairs$gene == 1 & tp$pairs$partner == 3], round(p13$r, 3))
  expect_false(any(tp$pairs$gene == tp$pairs$partner))
  expect_equal(tp$pairs$q[tp$pairs$gene == 1 & tp$pairs$partner == 2], signif(p.adjust(coexp_p(c(p12$r, p13$r), tp$df), "BH")[1], 3))
})

test_that("data-driven cut-off: maxstat estimate, corrected p-value and transfer to another cohort", {
  suppressPackageStartupMessages(library(maxstat))
  m <- cut_maxstat(sx); ref <- maxstat.test(Surv(time, event) ~ value, sx, smethod = "LogRank", pmethod = "Lau92", minprop = 0.2, maxprop = 0.8)
  expect_equal(m$cut, unname(ref$estimate)); expect_equal(m$p, min(1, unname(ref$p.value)))
  expect_true(m$q >= 0.19 && m$q <= 0.81)
  k <- km_split(sx, m$cut); expect_equal(unname(k$n[["High"]]), sum(sx$value > m$cut)); expect_equal(k$q, mean(sx$value <= m$cut))
  # the corrected p is larger than the naive log-rank p at the chosen cut-off
  expect_gt(m$p, k$p)
  other <- transform(sx, value = value + 5)
  tr <- transfer_split(other, m$q); expect_equal(tr$cut, unname(quantile(other$value, m$q))); expect_equal(tr$q, m$q)
  expect_null(cut_maxstat(transform(sx, event = 0L)))
})

test_that("Kaplan-Meier for several groups and the spline hazard curve", {
  g <- factor(cut(sx$value, quantile(sx$value, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("a", "b", "c")))
  k <- km_groups(sx, g)
  expect_equal(k$p, pchisq(survdiff(Surv(time, event) ~ g, sx)$chisq, 2, lower.tail = FALSE)); expect_identical(k$levels, c("a", "b", "c"))
  expect_null(km_groups(sx, factor(c(rep("a", n - 5), rep("b", 5)))))      # a group below SURV_MIN_GROUP
  sp <- surv_spline(sx)
  fit <- coxph(Surv(time, event) ~ splines::ns(value, df = 3), sx)
  expect_equal(sp$p_any, pchisq(2 * diff(fit$loglik), 3, lower.tail = FALSE))
  i <- which.min(abs(sp$curve$value - median(sx$value))); expect_lt(abs(log(sp$curve$hr[i])), 0.05)   # HR ~1 at the median
  expect_true(all(sp$curve$lo <= sp$curve$hr & sp$curve$hr <= sp$curve$hi))
})

test_that("gene sets: calls, counts, merged groups, filters and the Cox trend", {
  set.seed(3)
  M <- cbind(A = rnorm(60), B = rnorm(60), C = rnorm(60))
  cl <- geneset_calls(M, c("high", "low", "high"))
  expect_identical(unname(cl$calls[, "A"]), unname(M[, "A"] > median(M[, "A"]))); expect_identical(unname(cl$calls[, "B"]), unname(M[, "B"] <= median(M[, "B"])))
  expect_equal(unname(cl$count), unname(rowSums(cl$calls)))
  expect_equal(unname(cl$score), unname(drop(scale(M) %*% c(1, -1, 1)) / 3))
  ct <- geneset_calls(M, c("high", "low", "high"), "tertile")
  expect_equal(sum(ct$calls[, "A"]), sum(M[, "A"] > quantile(M[, "A"], 2 / 3)))
  gr <- geneset_groups(c(rep(0L, 3), rep(1L, 20), rep(2L, 25), rep(3L, 4)), 3)
  expect_identical(levels(gr), c("0-1 of 3", "2-3 of 3")); expect_equal(as.integer(table(gr)), c(23L, 29L))
  expect_identical(levels(geneset_groups(c(rep(0L, 12), rep(1L, 12), rep(2L, 12)), 2)), c("0 of 2", "1 of 2", "2 of 2"))
  gv <- list(M = M, clin = data.frame(patient = paste0("p", 1:60), age = round(runif(60, 40, 80)), sex = rep(c("Female", "Male"), 30),
                                      stage = rep(c("I", "II", "III"), 20), residual = NA_character_, os_time = rexp(60, 1 / 600), os_event = rbinom(60, 1, 0.7),
                                      pfi_time = NA, pfi_event = NA))
  x <- geneset_frame(gv, c("high", "low", "high"), "os")
  expect_equal(nrow(x), 60); expect_equal(x$value_count, unname(cl$count)); expect_equal(x$time, gv$clin$os_time / DAYS_PER_MONTH)
  xf <- geneset_frame(gv, c("high", "low", "high"), "os", filters = list(stage = "III", sex = c("Female", "Male")))
  expect_true(all(xf$stage == "III")); expect_equal(xf$value_count, unname(cl$count[gv$clin$stage == "III"]))   # calls made before filtering
  h <- geneset_cox(x, "graded"); fit <- coxph(Surv(time, event) ~ value_count, x)
  expect_equal(h$hr, unname(exp(coef(fit))))
  hs <- geneset_cox(x, "score"); expect_equal(hs$hr, unname(exp(coef(coxph(Surv(time, event) ~ scale(score)[, 1], x)))))
  expect_null(geneset_frame(list(missing = "X"), "high", "os"))
})
