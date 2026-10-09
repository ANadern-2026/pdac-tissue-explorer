# Runs data-raw/01_scatlas_summaries.R on a miniature atlas and checks its group summaries against a
# direct per-cell recomputation. Needs SeuratObject + Matrix; skipped otherwise.
test_that("single-cell summaries match per-cell recomputation", {
  skip_if_not_installed("SeuratObject")
  suppressPackageStartupMessages({ library(Matrix); library(SeuratObject) })
  root <- normalizePath(file.path(test_path(), "..", ".."))
  source(file.path(root, "tests", "fixtures", "make_mini_atlas.R"))
  tmp <- withr_tempdir <- tempfile("mini"); dir.create(tmp); make_mini_atlas(tmp)
  out <- file.path(tmp, "out"); dir.create(out)
  script <- readLines(file.path(root, "data-raw", "01_scatlas_summaries.R"))
  script <- sub('^ATLAS_RDS .*', sprintf('ATLAS_RDS <- "%s"', file.path(tmp, "sc_ref_base.rds")), script)
  script <- sub('^DUCTAL_LABELS .*', sprintf('DUCTAL_LABELS <- "%s"', file.path(tmp, "DuctalTumor_Subcluster.txt.gz")), script)
  script <- sub('^FIBRO_LABELS .*', sprintf('FIBRO_LABELS <- "%s"', file.path(tmp, "Fibroblast_Subcluster.txt.gz")), script)
  script <- sub('^OUT_DIR .*', sprintf('OUT_DIR <- "%s"', out), script); script <- sub('^LOG_DIR .*', sprintf('LOG_DIR <- "%s"', out), script)
  script <- sub('^MIN_FEATURES .*', "MIN_FEATURES <- 10", script); script <- sub('^MIN_GENE_CELLS .*', "MIN_GENE_CELLS <- 5", script)
  script <- sub('^UMAP_CELLS .*', "UMAP_CELLS <- 1000", script); script <- sub('^UMAP_MIN_TYPE .*', "UMAP_MIN_TYPE <- 100", script)
  f <- file.path(tmp, "run01.R"); writeLines(script, f)
  expect_no_error(suppressMessages(capture.output(source(f, local = new.env()))))

  g <- readRDS(file.path(out, "sc_groups.rds")); e <- readRDS(file.path(out, "sc_expr.rds")); u <- readRDS(file.path(out, "sc_umap.rds"))
  o <- readRDS(file.path(tmp, "sc_ref_base.rds")); d <- LayerData(o, layer = "data"); md <- o@meta.data; keep <- md$nFeature_RNA >= 10
  dl <- read.delim(file.path(tmp, "DuctalTumor_Subcluster.txt.gz")); fl <- read.delim(file.path(tmp, "Fibroblast_Subcluster.txt.gz"))
  fine <- md$Clusters; cn <- rownames(md)
  fine[md$Clusters == "DUCTAL"] <- ifelse(cn[md$Clusters == "DUCTAL"] %in% rownames(dl), dl[cn[md$Clusters == "DUCTAL"], "seurat_clusters"], "Ductal (unlabelled)")
  fine[md$Clusters == "FIBROBLASTS"] <- fl[cn[md$Clusters == "FIBROBLASTS"], "Coarse"]
  set.seed(1)
  for (k in sample(nrow(g), 20)) {
    cells <- which(keep & md$Name == g$sample[k] & fine == g$fine_group[k])
    expect_equal(length(cells), g$n_cells[k])
    for (gene in c("MSLN", "KRT19", "G50")) {
      v <- d[gene, cells]
      expect_equal(unname(e$mean_cp10k_x100[gene, g$group_id[k]] / 100), mean(expm1(v)), tolerance = 0.006, scale = 1)
      expect_equal(unname(e$frac_x10000[gene, g$group_id[k]] / 1e4), mean(v > 0), tolerance = 1e-4, scale = 1)
    }
  }
  expect_setequal(unique(g$treatment), c("Treatment-naive", "Donor (untreated)", "Treated"))
  expect_true(all(is.na(g$met_site[g$disease_state != "Metastatic lesion"])))      # stray "2" level removed
  expect_identical(colnames(u$expr), e$genes)
  expect_equal(nrow(u$coords), nrow(u$expr))
})
