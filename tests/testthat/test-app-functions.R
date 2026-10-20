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

test_that("vectorised bootstrap interval is identical to the v0.1 replicate() version", {
  ref <- function(x, n = BOOT_N, seed = 1L) { if (length(x) < 3) return(c(NA_real_, NA_real_))
    set.seed(seed); b <- replicate(n, stats::median(sample(x, replace = TRUE))); unname(stats::quantile(b, c(0.025, 0.975))) }
  set.seed(9)
  for (i in 1:60) { k <- sample(3:70, 1); x <- round(rnorm(k), sample(1:4, 1)); if (i %% 4 == 0) x[1:3] <- x[4]; expect_identical(boot_median_ci(x), ref(x)) }
  expect_identical(boot_median_ci(c(1, 2)), c(NA_real_, NA_real_))
})

test_that("phosphosite statistics: minimum pairs, relative to protein, FDR within the gene", {
  set.seed(4); cs <- paste0("C3L-", 1:20)
  tu <- rbind(S1 = rnorm(20, 20) + 2, T5 = rnorm(20, 18), Y9 = c(rnorm(5, 15), rep(NA, 15)))
  no <- rbind(S1 = rnorm(20, 20), T5 = tu["T5", ] + rep(c(0.1, -0.1), 10), Y9 = rnorm(20, 15)); colnames(tu) <- colnames(no) <- cs   # T5: no change
  ph <- list(sites = data.frame(label = rownames(tu), site_id = paste0("NP_1.1_", rownames(tu)), residue = c("S", "T", "Y"), position = c(1L, 5L, 9L),
                                refseq = "NP_1.1", source_symbol = "G", peptide = "x"),
             tumor = tu, normal = no, has_protein = TRUE, prot_tumor = rep(2, 20), prot_normal = rep(0, 20))
  a <- phospho_table(ph); r <- phospho_table(ph, relative = TRUE)
  expect_identical(a$n_pairs, c(20L, 20L, 5L)); expect_true(is.na(a$p[3]))                      # Y9: 5 pairs < PHOSPHO_MIN_PAIRS
  expect_lt(a$q[1], 0.05); expect_equal(r$median_diff, a$median_diff - 2)                         # protein up by 2 in every tumour
  expect_identical(a$q[1:2], p.adjust(a$p[1:2], "BH"))
  expect_identical(phospho_table(ph, ci = FALSE)$p, a$p); expect_true(all(is.na(phospho_table(ph, ci = FALSE)$ci_low)))
  txt <- txt_phospho(a, r, "G", 20)
  expect_match(txt, "3 sites quantified; 2 compared in at least 10 matched pairs: 1 higher and 0 lower in tumour")
  expect_match(txt, "largest change at S1"); expect_match(txt, "Relative to the protein level")
  expect_match(txt_phospho(a, NULL, "G", 20), "protein itself was not quantified")
  expect_match(txt_phospho(NULL, NULL, "G", 20), "no G phosphosites")
  expect_s3_class(plot_phospho_sites(a, "G"), "ggplot"); expect_s3_class(plot_phospho_sites(NULL, "G", has_sites = FALSE), "ggplot")
})

test_that("links out use identifiers when known and searches otherwise", {
  gr <- data.frame(symbol = "MSLN", entrez = "10232", stringsAsFactors = FALSE)
  l <- gene_links(gr)
  expect_identical(names(l), c("ncbi", "uniprot", "hpa", "opentargets", "depmap"))
  expect_identical(l$ncbi$url, "https://www.ncbi.nlm.nih.gov/gene/10232"); expect_match(l$hpa$url, "proteinatlas.org/search/MSLN$")
  x <- data.frame(symbol = c("MSLN", "MIR1"), entrez = c("10232", "406904"), ensembl = c("ENSG00000102854", NA), gene_type = c("protein-coding", "ncRNA"))
  g <- add_xrefs(data.frame(symbol = c("MSLN", "MIR1", "ABC"), entrez = c("10232", "999", NA), stringsAsFactors = FALSE), x)
  expect_identical(g$ensembl, c("ENSG00000102854", NA, NA)); expect_identical(g$gene_type, c("protein-coding", NA, NA))   # MIR1: Entrez disagrees
  l <- gene_links(g[1, ])
  expect_identical(l$hpa$url, "https://www.proteinatlas.org/ENSG00000102854"); expect_identical(l$opentargets$url, "https://platform.opentargets.org/target/ENSG00000102854")
  nc <- gene_links(data.frame(symbol = "MIR1 (Entrez 406904)", entrez = NA, ensembl = NA, gene_type = "ncRNA"))
  expect_identical(names(nc), c("ncbi", "opentargets", "depmap")); expect_match(nc$ncbi$url, "term=MIR1%5Bsym%5D"); expect_match(nc$depmap$url, "/gene/MIR1\\?")
})
