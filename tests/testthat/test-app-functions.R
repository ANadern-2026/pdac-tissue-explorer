# Checks the app's statistics, aggregation and text on a small synthetic dataset with known answers.
root <- normalizePath(file.path(test_path(), "..", ".."))
for (f in sort(list.files(file.path(root, "app", "R"), full.names = TRUE))) source(f)

test_that("paired statistics recover a known shift", {
  set.seed(2); n <- 30; normal <- rnorm(n, 10); tumour <- normal + 1 + rnorm(n, 0, 0.2)
  s <- paired_stats(list(pairs = data.frame(case_id = paste0("P", 1:n), normal = normal, tumour = tumour, diff = tumour - normal)))
  expect_true(s$ok); expect_equal(s$n, n); expect_equal(s$median_diff, median(tumour - normal))
  expect_lt(s$p, 1e-6); expect_true(s$ci_low < s$median_diff && s$median_diff < s$ci_high); expect_equal(s$n_up, sum(tumour > normal))
  expect_false(paired_stats(list(pairs = data.frame(case_id = "a", normal = 1, tumour = 2, diff = 1)))$ok)
})

test_that("weighted aggregation of fine groups is exact", {
  v <- data.frame(sample = c("s1", "s1", "s2"), compartment = c("Fibroblasts", "Fibroblasts", "Fibroblasts"),
                  n_cells = c(10, 30, 5), mean_cp10k = c(2, 6, 1), frac = c(0.5, 0.1, 1))
  a <- sc_aggregate(v, c("sample", "compartment"))
  expect_equal(a$n_cells[a$sample == "s1"], 40)
  expect_equal(a$mean_cp10k[a$sample == "s1"], (10 * 2 + 30 * 6) / 40)
  expect_equal(a$frac[a$sample == "s1"], (10 * 0.5 + 30 * 0.1) / 40)
  expect_equal(a$avg_expr, log1p(a$mean_cp10k))
})

test_that("cell-type summary: transcript shares sum to one and respect the cell minimum", {
  g <- expand.grid(sample = paste0("s", 1:6), cell_type = c("FIBROBLASTS", "TNK"), stringsAsFactors = FALSE)
  g$fine_group <- ifelse(g$cell_type == "FIBROBLASTS", "myFb 1", "TNK"); g$n_cells <- c(rep(50, 6), rep(10, 5), 40)
  g$disease_state <- "Primary tumor"; g$treatment <- "Treatment-naive"; g$met_site <- NA; g$study <- "GSE1"
  g <- prepare_sc_groups(g); g$mean_cp10k <- ifelse(g$cell_type == "FIBROBLASTS", 3, 1); g$frac <- 0.4
  cs <- celltype_summary(g, "Primary tumour")
  expect_equal(sum(cs$summary$share), 1)
  expect_equal(cs$summary$n_samples[cs$summary$compartment == "T and NK cells"], 1)   # only one sample has >= 20 T/NK cells
  expect_equal(cs$summary$share[cs$summary$compartment == "Fibroblasts"], (300 * 3) / (300 * 3 + 90 * 1))
})

test_that("tissue states are assigned as documented", {
  g <- data.frame(sample = 1:6, cell_type = "TNK", fine_group = "TNK", n_cells = 30, study = "x",
                  disease_state = c("Donor", "Adjacent normal", "Primary tumor", "Primary tumor", "Metastatic lesion", "Metastatic lesion"),
                  treatment = c("Donor (untreated)", "Treatment-naive", "Treatment-naive", "Treated", "Treatment-naive", "Treatment-naive"),
                  met_site = c(NA, NA, NA, NA, "Liver", "Lung"))
  s <- as.character(prepare_sc_groups(g)$state)
  expect_identical(s, c("Healthy donor", "Adjacent normal", "Primary tumour", "Primary tumour (treated)", "Liver metastasis", NA))
})

test_that("plain-language text handles missing data gracefully", {
  expect_match(txt_cptac(NULL, NULL, "XYZ", "protein"), "not in the protein data")
  expect_match(txt_sc(NULL, "XYZ"), "not in the single-cell data")
  expect_match(fold_txt(log2(186.3)), "about 186-fold higher")
  expect_match(fold_txt(-1), "about 2.0-fold lower")
  expect_identical(fold_txt(0.1), "")
  zero <- list(pairs = data.frame(case_id = paste0("P", 1:8), normal = 0, tumour = 0, diff = 0), n_tumours_measured = 8, n_tumours_total = 8)
  expect_match(txt_cptac(zero, paired_stats(zero), "XYZ", "rna"), "no clear difference.*p = n/a")   # all values identical: p is NA
})
