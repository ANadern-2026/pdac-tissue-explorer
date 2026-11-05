# Model systems (data-raw/06, app/R/07_models.R): store round trip on a tiny synthetic layer, plain-language answers, labels.
root <- normalizePath(file.path(test_path(), "..", ".."))
suppressPackageStartupMessages({ library(ggplot2) })
for (f in sort(list.files(file.path(root, "app", "R"), full.names = TRUE))) source(f)
local({ old <- options(pdacx.build_store.autorun = FALSE); on.exit(options(old)); source(file.path(root, "data-raw", "05_build_store.R"), local = FALSE) })   # option restored afterwards

toy_models <- function() {
  set.seed(3)
  models <- data.frame(model = c("ACH-1", "ACH-2", "ACH-3", "ORG-1", "N-1"), name = c("L1", "L2", "O-DM", "O-T", "N-T"),
                       source = c("DepMap", "DepMap", "DepMap", "Tiriac 2018", "Tiriac 2018"), type = c("Cell line", "Cell line", "Cell line", "Organoid", "Normal organoid"),
                       tumour = "Adenocarcinoma", origin = "Primary tumour", collection_site = "pancreas", sex = NA, age = NA, depmap_id = c("ACH-1", "ACH-2", "ACH-3", NA, NA),
                       hcmi_id = NA, rrid = NA, has_rna = TRUE, moffitt = c(30, 15, -20, -35, 0), sim_tumour = c(0.8, 0.82, 0.86, 0.84, 0.78), closest = SIM_TUMOUR,
                       has_crispr = c(TRUE, TRUE, TRUE, FALSE, FALSE), has_protein = c(TRUE, FALSE, FALSE, FALSE, FALSE), supplier = c("ATCC", "Academic lab", NA, NA, NA),
                       catalog = c("CRL-1", NA, NA, NA, NA), subtype_rel = c(0.8, -0.5, 0.1, -0.2, NA), sim_livermet = 0.8, family = NA, stringsAsFactors = FALSE)
  units <- data.frame(unit = c(paste0(SIM_TUMOUR, "|s", 1:3), paste0("Fibroblasts|s", 1:3), "TCGA|t1"), group = c(rep(SIM_TUMOUR, 3), rep("Fibroblasts", 3), "PDAC bulk tumours (TCGA)"),
                      sample = c("s1", "s2", "s3", "s1", "s2", "s3", "t1"), study = "x", n_cells = 50L, moffitt = c(-10, 5, 20, NA, NA, NA, 0), stringsAsFactors = FALSE)
  g <- c("KRT17", "COL1A1", "KRAS"); lines <- data.frame(model = c(models$model[1:3], paste0("X", 1:6)), name = c("L1", "L2", "O-DM", paste0("other", 1:6)),
                                                        lineage = c(rep("Pancreas", 3), rep("Lung", 6)), code = c(rep("PAAD", 3), rep("LUAD", 6)),
                                                        type = c("Cell Line", "Cell Line", "Organoid", rep("Cell Line", 6)), pdac = c(rep(TRUE, 3), rep(FALSE, 6)))
  x <- list(models = models, rna_log2tpm_x1000 = matrix(c(9000L, 8000L, 7000L, 6000L, 0L, 1000L, 2000L, 1500L, 1200L, 900L, 5000L, 5100L, 5200L, 5300L, 4000L), 3, 5, byrow = TRUE, dimnames = list(g, models$model)),
       universe = g[1:2], pct_models_x10 = matrix(c(950L, 900L, 300L, 250L, 0L, 400L, 420L, 380L, 390L, 410L), 2, 5, byrow = TRUE, dimnames = list(g[1:2], models$model)),
       tissue_units = units, pct_tissue_x10 = matrix(c(990L, 985L, 995L, 0L, 0L, 10L, 980L, 100L, 120L, 90L, 1000L, 995L, 990L, 1000L), 2, 7, byrow = TRUE, dimnames = list(g[1:2], units$unit)),
       similarity = matrix(c(0.8, 0.5, 0.82, 0.5, 0.86, 0.45, 0.84, 0.44, 0.78, 0.43), 5, 2, byrow = TRUE, dimnames = list(models$model, c(SIM_TUMOUR, "Fibroblasts"))),
       crispr_x1000 = matrix(c(-1800L, -1500L, -1200L, -400L, -500L, -600L, -300L, -450L, -550L,
                                                          10L, -20L, 30L, 0L, 5L, -10L, 15L, -5L, 20L), 2, 9, byrow = TRUE, dimnames = list(c("KRAS", "KRT17"), lines$model)),
       crispr_lines = lines, protein_x1000 = matrix(c(500L, NA, -100L), 1, 3, dimnames = list("KRAS", c("ACH-1", "X1", "X2"))),
       protein_lines = data.frame(model = c("ACH-1", "X1", "X2"), name = c("L1", "o1", "o2"), lineage = c("Pancreas", "Lung", "Lung"), pdac = c(TRUE, FALSE, FALSE)),
       moffitt_genes = list(basal = "KRT17", classical = "GATA6"), source = list(depmap_release = "test", references = c(primary = "atlas tumour cells, test")),
       rna_lines = data.frame(model = c("ACH-1", "ACH-2", "ACH-3", paste0("X", 1:6)), name = c("L1", "L2", "O-DM", paste0("other", 1:6)), lineage = c(rep("Pancreas", 3), rep("Lung", 6)),
                              type = rep("Cell Line", 9), pdac = c(rep(TRUE, 3), rep(FALSE, 6))),
       rna_lines_x1000 = matrix(c(9000L, 8000L, 7000L, 1000L, 1200L, 900L, 1100L, 800L, 1000L, 5000L, 5100L, 5200L, 5000L, 4900L, 5300L, 5100L, 5000L, 4800L), 2, 9, byrow = TRUE,
                                dimnames = list(c("KRT17", "KRAS"), c("ACH-1", "ACH-2", "ACH-3", paste0("X", 1:6)))),
       pair_organoid_x1000 = matrix(c(1:12 * 500L), 1, 12, dimnames = list("KRT17", paste0("H", 1:12))), pair_tumour_x1000 = matrix(c(1:12 * 400L), 1, 12, dimnames = list("KRT17", paste0("T", 1:12))))
  x$hcmi_pairs <- data.frame(model = paste0("H", 1:12), case = paste0("C", 1:12), tumour = paste0("T", 1:12), r_own = 0.3, r_other = 0, rank_own = 1L)
  x
}

test_that("models layer: every record survives the store round trip", {
  dir <- tempfile(); dir.create(dir); x <- toy_models(); saveRDS(x, file.path(dir, "models.rds"))
  build_models(dir, function(...) invisible(NULL))
  d <- list(has_models = TRUE, models = store_open(file.path(dir, STORE_LAYERS[["models"]])))
  d$model_info <- d$models$info[c("models", "tissue_units", "crispr_lines", "protein_lines", "rna_lines", "similarity", "hcmi_pairs", "moffitt_genes", "n_universe", "source")]
  on.exit({ store_close(d$models); store_cache_clear(); unlink(dir, recursive = TRUE) })
  v <- model_values(d, "KRT17")
  expect_equal(v$rna$log2tpm, unname(x$rna_log2tpm_x1000["KRT17", ] / 1000))
  expect_equal(v$rna$pct, unname(x$pct_models_x10["KRT17", ] / 10))
  expect_equal(v$tissue$pct, unname(x$pct_tissue_x10["KRT17", ] / 10))
  expect_equal(v$crispr$effect, unname(x$crispr_x1000["KRT17", ] / 1000))
  expect_null(v$protein)
  k <- model_values(d, "KRAS"); expect_true(all(is.na(k$rna$pct))); expect_null(k$tissue); expect_equal(k$protein$value, c(0.5, NA, -0.1))
  expect_null(model_values(d, "NOTAGENE")); expect_equal(d$model_info$n_universe, 2L)

  expect_equal(v$lines$log2tpm, c(9, 8, 7, 1, 1.2, 0.9, 1.1, 0.8, 1))
  expect_equal(v$pairs$organoid_log2tpm, 1:12 * 0.5); expect_equal(v$pairs$tumour_log2tpm, 1:12 * 0.4)
  a <- model_answers(v, "KRT17")
  expect_equal(a$lines$a, "Similar to other cancer lines")            # 8 vs 1 but 3 vs 6 lines: p = 0.024, not < 0.01 required for "higher"
  expect_equal(a$lines$d[1], "L1 (9.00; ATCC CRL-1)"); expect_equal(a$lines$d[2], "L2 (8.00)")   # no supplier without a catalogue number
  expect_equal(a$dep$a, "No clear effect")
  expect_equal(a$org$a, "Yes, largely")                               # perfect rank agreement across 12 pairs
  ak <- model_answers(k, "KRAS"); expect_equal(ak$org$a, "Not measured"); expect_equal(ak$dep$a, "Some PDAC models do")   # 3 vs 6 lines: p cannot reach < 0.01
  expect_equal(ak$dep$d, "Knocking it out slowed growth clearly (effect below -0.5) in 3 of 3 PDAC models and in 33% of other cancer lines.")
  s <- dependency_stats(k$crispr); expect_equal(s$med_pdac, -1.5); expect_equal(s$med_other, -0.475); expect_equal(s$n_dep, 3L)
  ls <- line_expr_stats(v$lines); expect_equal(ls$n_pdac, 3L); expect_equal(ls$n_other, 6L)
  expect_match(txt_models(v, "KRT17"), "^Model systems: in PDAC cell lines similar to other cancer lines; knockout effect: no clear effect")
  # figures and tables build without errors
  for (p in list(plot_line_expression(v, "KRT17"), plot_dependency(k$crispr, "KRAS"), plot_line_subtype(v, "KRT17"), plot_model_protein(k$protein, "KRAS"),
                 plot_organoid_expression(v, "KRT17"), plot_hcmi_gene_pairs(v, "KRT17"), plot_tissue_background(v, "KRT17", 2L),
                 plot_overview_lines(d$model_info), plot_overview_organoids(d$model_info), plot_similarity_groups(d$model_info), plot_hcmi_pairs(d$model_info)))
    expect_s3_class(ggplot_build(p), "ggplot_built")
  expect_equal(nrow(model_table(v)), 5L); expect_equal(nrow(lines_table(d$model_info)), 3L); expect_equal(nrow(organoids_table(d$model_info)), 2L)
  expect_equal(nrow(fidelity_table(d$model_info)), 5L)
})

test_that("labels", {
  expect_equal(subtype_rel_label(c(-1, 0, 1, NA), c(-1, -0.5, 0, 0.5, 1)), c("more classical", "intermediate", "more basal-like", NA))
  expect_equal(fmt_supplier(c("ATCC", "Academic lab", NA), c("CRL-1", NA, "X")), c("ATCC CRL-1", "", ""))
  m <- data.frame(source = c("DepMap", "DepMap", "HCMI", "Tiriac 2018", "Tiriac 2018"), type = c("Cell line", "Organoid", "Organoid", "Organoid", "Normal organoid"))
  expect_equal(as.character(model_group(m)), MODEL_GROUPS)
  expect_equal(nth(c(1)), "1st"); expect_equal(nth(12), "12th"); expect_equal(nth(23), "23rd"); expect_equal(nth(NA), "n/a")
})
