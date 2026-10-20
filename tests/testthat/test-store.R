# Builds the on-disk store (data-raw/05_build_store.R) from small synthetic v0.1-style tables and checks that the
# app's accessors return exactly what direct indexing of the in-memory matrices returned in v0.1.
root <- normalizePath(file.path(test_path(), "..", ".."))
for (f in sort(list.files(file.path(root, "app", "R"), full.names = TRUE))) source(f)
suppressPackageStartupMessages(library(Matrix))
old_opt <- options(pdacx.build_store.autorun = FALSE, pdacx.root = root)
source(file.path(root, "data-raw", "05_build_store.R"))
options(old_opt)

make_v01_tables <- function(dir) {
  set.seed(11)
  mk <- function(genes, cases, na = 0.2, digits = 3) { m <- matrix(round(rnorm(length(genes) * length(cases), 20, 3), digits), length(genes), dimnames = list(genes, cases)); m[runif(length(m)) < na] <- NA; m }
  pg <- c("MSLN", "KRT19", "COL1A1", paste0("P", 1:27)); rg <- c("MSLN", "KRT19", "COL1A1", "GENE.X", paste0("R", 1:36))
  tc <- paste0("C3L-", sprintf("%05d", 1:30)); nc <- sample(tc, 12)
  cptac <- list(protein_tumor = mk(pg, tc), protein_normal = mk(pg, nc), rna_tumor = mk(rg, tc, 0), rna_normal = mk(rg, nc[1:6], 0),
                clinical = data.frame(case_id = tc, stage = "II"), source = list(note = "synthetic"))
  cptac$rna_tumor["GENE.X", ] <- 1.123456789                                      # not exactly representable at 3 decimals
  bs <- c(paste0("GTEX-", 1:8), paste0("TCGA-", 1:10))
  bulk <- list(log2tpm_x1000 = matrix(sample(0:15000, 22 * 18, TRUE), 22, dimnames = list(c("MSLN", "KRT19", paste0("B", 1:20)), bs)),
               samples = data.frame(sample = bs, cohort = c(rep("Healthy pancreas (GTEx)", 8), rep("Adjacent normal (TCGA)", 2), rep("PDAC (TCGA)", 8))),
               signatures = list(), source = list(note = "synthetic"))
  groups <- expand.grid(sample = paste0("s", 1:8), fine_group = c("myFb 1", "iFb 2", "Ductal 4", "TNK"), stringsAsFactors = FALSE)
  groups$group_id <- seq_len(nrow(groups)); groups$n_cells <- sample(5:80, nrow(groups), TRUE)
  groups$cell_type <- c(FIBROBLASTS = "FIBROBLASTS", DUCT = "DUCTAL", TNK = "TNK")[ifelse(grepl("Fb", groups$fine_group), "FIBROBLASTS", ifelse(groups$fine_group == "TNK", "TNK", "DUCT"))]
  groups$study <- "GSE1"; groups$disease_state <- ifelse(groups$sample %in% c("s7", "s8"), "Adjacent normal", "Primary tumor")
  groups$met_site <- NA; groups$treatment <- "Treatment-naive"; groups <- groups[, c("group_id", "sample", "fine_group", "n_cells", "cell_type", "study", "disease_state", "met_site", "treatment")]
  sg <- c("MSLN", "KRT19", "COL1A1", paste0("S", 1:17)); ng <- nrow(groups)
  mean_x <- matrix(sample(c(0L, 0L, 1:900000), length(sg) * ng, TRUE), length(sg), dimnames = list(sg, NULL))
  expr <- list(genes = sg, mean_cp10k_x100 = mean_x, frac_x10000 = matrix(sample(0:10000, length(sg) * ng, TRUE), length(sg), dimnames = list(sg, NULL)), note = "synthetic")
  X <- rsparsematrix(500, length(sg), 0.1, rand.x = function(n) round(runif(n, 0, 9), 2)); colnames(X) <- sg
  X@x[1:5] <- 0                                                                  # explicit zeros, as in the real atlas
  X[, "S17"] <- 0; X <- as(X, "CsparseMatrix")                                     # one gene not detected in any subsampled cell
  umap <- list(coords = data.frame(UMAP_1 = rnorm(500), UMAP_2 = rnorm(500)),
               meta = data.frame(sample = "s1", study = "GSE1", disease_state = "Primary tumor", met_site = NA, treatment = "Treatment-naive",
                                 cell_type = sample(c("DUCTAL", "TNK"), 500, TRUE), fine_group = "TNK"), expr = X, note = "synthetic")
  all <- unique(c(pg, rg, rownames(bulk$log2tpm_x1000), sg, "NOWHERE"))
  genes <- data.frame(symbol = all, sym_cptac_protein = ifelse(all %in% pg, all, NA), sym_cptac_rna = ifelse(all %in% rg, all, NA),
                      sym_tcga_gtex = ifelse(all %in% rownames(bulk$log2tpm_x1000), all, NA), sym_sc = ifelse(all %in% sg, all, NA),
                      name = NA_character_, aliases = "", stringsAsFactors = FALSE)
  genes <- genes[genes$symbol != "NOWHERE", ]
  saveRDS(cptac, file.path(dir, "cptac.rds")); saveRDS(bulk, file.path(dir, "bulk_tcga_gtex.rds")); saveRDS(groups, file.path(dir, "sc_groups.rds"))
  saveRDS(expr, file.path(dir, "sc_expr.rds")); saveRDS(umap, file.path(dir, "sc_umap.rds")); saveRDS(list(genes = genes), file.path(dir, "genes.rds"))
  list(cptac = cptac, bulk = bulk, groups = groups, expr = expr, umap = umap, genes = genes)
}

tmp <- tempfile("store"); dir.create(tmp); src <- make_v01_tables(tmp)
suppressMessages(build_store(tmp, log = FALSE))
D <- load_app_data(tmp)

test_that("store files are built and the large tables are not held in memory", {
  expect_true(all(file.exists(file.path(tmp, STORE_LAYERS[c("cptac", "bulk", "sc")]))))
  expect_true(is.environment(D$sc) && is.environment(D$cptac) && is.environment(D$bulk))
  expect_null(D$sc_mean); expect_null(D$umap$expr)
  expect_identical(D$sc_groups, prepare_sc_groups(src$groups))
  expect_identical(D$umap$coords, src$umap$coords)
})

test_that("CPTAC pairs are identical to direct matrix indexing (v0.1)", {
  for (g in c("MSLN", "KRT19", "P5", "R10", "GENE.X")) for (a in c("protein", "rna")) {
    s <- gene_row(D, g)[[paste0("sym_cptac_", a)]]
    if (is.na(s)) { expect_null(cptac_pairs(D, g, a)); next }
    t <- src$cptac[[paste0(a, "_tumor")]][s, ]; n <- src$cptac[[paste0(a, "_normal")]][s, ]
    ids <- intersect(names(t), names(n)); ok <- ids[is.finite(t[ids]) & is.finite(n[ids])]
    p <- cptac_pairs(D, g, a)
    expect_identical(p$pairs, data.frame(case_id = ok, normal = unname(n[ok]), tumour = unname(t[ok]), diff = unname(t[ok] - n[ok])))
    expect_identical(p$tumour_values, t[is.finite(t)]); expect_identical(p$n_tumours_total, length(t))
  }
})

test_that("bulk, single-cell and atlas-map values are identical to v0.1", {
  smp <- src$bulk$samples; smp$cohort <- factor(smp$cohort, levels = BULK_LEVELS)
  expect_identical(bulk_values(D, "KRT19"), data.frame(smp, value = src$bulk$log2tpm_x1000["KRT19", ] / 1000, check.names = FALSE))
  expect_null(bulk_values(D, "S3"))
  g <- prepare_sc_groups(src$groups)
  for (s in c("MSLN", "S4", "S17")) {
    expect_identical(sc_values(D, s), data.frame(g, mean_cp10k = src$expr$mean_cp10k_x100[s, ] / 100, frac = src$expr$frac_x10000[s, ] / 1e4))
    expect_identical(umap_values(D, s), as.numeric(src$umap$expr[, s]))
  }
  expect_null(sc_values(D, "P5")); expect_null(umap_values(D, NA_character_))
  expect_identical(multi_dot_table(D, c("MSLN", "KRT19", "S2"))$gene, factor(rep(c("MSLN", "KRT19", "S2"), each = 3), levels = c("MSLN", "KRT19", "S2")))
})

test_that("codecs are exact and the record cache stays within its budget", {
  x <- list(a = c(1.5, NA, 28.674), b = 1:3, c = c(1 / 3, 2), d = c(NaN, 1))
  expect_identical(store_unpack(pack_rds(x), "rds-gzip"), x)
  expect_s3_class(store_scale(x)$a, "pdacx_scaled"); expect_false(inherits(store_scale(x)$c, "pdacx_scaled"))
  expect_error(pack_sparse16(c(1L, 3L), c(0.5, 1 / 3), 5L), "cannot represent")
  expect_identical(store_unpack(pack_sparse16(c(2L, 5L), c(0.25, 9.05), 6L), "sparse16-bzip2"), c(0, 0.25, 0, 0, 9.05, 0))
  store_cache_clear(); old <- STORE_CACHE_MB; STORE_CACHE_MB <<- 0.01            # ~10 kB: room for a couple of records
  for (s in c("MSLN", "KRT19", "COL1A1", "S1", "S2")) umap_values(D, s)
  expect_lte(sum(.store_cache$bytes), 0.01 * 2^20 + max(.store_cache$bytes)); expect_true(length(.store_cache$keys) >= 1)
  expect_identical(umap_values(D, "MSLN"), as.numeric(src$umap$expr[, "MSLN"]))  # evicted records are read again correctly
  STORE_CACHE_MB <<- old; store_cache_clear()
})

test_that("the optional phosphosite layer is skipped without its source, then built and read exactly", {
  expect_false(D$has_phospho); expect_null(phospho_sites(D, "MSLN"))
  set.seed(12)
  sites <- data.frame(site_id = c("NP_000001.1_S5", "NP_000001.1_T12", "NP_000001.1_S100", "NP_000002.1_S10", "NP_000002.1_Y3", "NP_000003.1_S7", "NP_000009.1_S1"),
                      gene = c("MSLN", "MSLN", "MSLN", "KRT19", "KRT19", "K19OLD", "NOGENE"), stringsAsFactors = FALSE)
  sites$refseq <- sub("_[STY][0-9]+$", "", sites$site_id); sites$residue <- substr(sub("^.*_", "", sites$site_id), 1, 1)
  sites$position <- as.integer(sub("^.*_[STY]", "", sites$site_id)); sites$peptide <- "PEPsTIDE"
  o <- order(sites$gene, sites$refseq, sites$position); sites <- sites[o, ]; rownames(sites) <- NULL
  tc <- colnames(src$cptac$protein_tumor); nc <- colnames(src$cptac$protein_normal)
  mk <- function(cols) { m <- matrix(round(rnorm(nrow(sites) * length(cols), 20, 2), 3), nrow(sites), dimnames = list(sites$site_id, cols)); m[runif(length(m)) < 0.25] <- NA; m }
  ph_src <- list(sites = sites, phospho_tumor = mk(rev(tc)), phospho_normal = mk(nc), source = list(note = "synthetic"))   # sample order differs from protein
  saveRDS(ph_src, file.path(tmp, "cptac_phospho.rds"))
  g <- readRDS(file.path(tmp, "genes.rds")); g$genes$aliases[g$genes$symbol == "KRT19"] <- "K19OLD, CK19"; saveRDS(g, file.path(tmp, "genes.rds"))
  msg <- capture_messages(build_store(tmp, layers = "phospho", log = FALSE))
  expect_true(any(grepl("not linked: NOGENE", msg))); expect_true(any(grepl("combined: K19OLD \\+ KRT19", msg)))
  D2 <- load_app_data(tmp); on.exit(close_app_data(D2))
  expect_true(D2$has_phospho)
  ph <- phospho_sites(D2, "KRT19")
  expect_identical(ph$sites$site_id, c("NP_000002.1_Y3", "NP_000002.1_S10", "NP_000003.1_S7"))           # isoform with most sites first, then position
  expect_identical(ph$sites$label, c("Y3 (NP_000002.1)", "S10 (NP_000002.1)", "S7 (NP_000003.1)"))
  expect_identical(ph$sites$source_symbol, c("KRT19", "KRT19", "K19OLD"))
  expect_identical(unname(ph$tumor), unname(ph_src$phospho_tumor[ph$sites$site_id, ])); expect_identical(colnames(ph$tumor), rev(tc))
  expect_identical(unname(ph$normal), unname(ph_src$phospho_normal[ph$sites$site_id, ]))
  pt <- src$cptac$protein_tumor["KRT19", ]; pn <- src$cptac$protein_normal["KRT19", ]
  expect_identical(ph$prot_tumor, unname(pt[rev(tc)])); expect_identical(ph$prot_normal, unname(pn[nc]))
  t <- ph_src$phospho_tumor["NP_000002.1_S10", ] - pt[names(ph_src$phospho_tumor["NP_000002.1_S10", ])]
  n <- ph_src$phospho_normal["NP_000002.1_S10", ] - pn[nc]
  ids <- intersect(names(t), names(n)); ok <- ids[is.finite(t[ids]) & is.finite(n[ids])]
  expect_identical(phospho_pairs(ph, 2, relative = TRUE), data.frame(case_id = ok, normal = unname(n[ok]), tumour = unname(t[ok]), diff = unname(t[ok] - n[ok])))
  expect_identical(phospho_sites(D2, "MSLN")$sites$label, c("S5", "T12", "S100"))
  expect_null(phospho_sites(D2, "COL1A1"))
  tab <- phospho_table(ph, min_pairs = 3)
  expect_identical(tab$n_pairs, vapply(1:3, function(k) nrow(phospho_pairs(ph, k)), integer(1)))
  expect_identical(tab$q, p.adjust(tab$p, "BH"))
  # the builder keeps an existing phospho.sqlite when its source is absent (server installs without the pipeline tables)
  file.remove(file.path(tmp, "cptac_phospho.rds"))
  msg <- capture_messages(build_store(tmp, layers = "phospho", log = FALSE)); expect_true(any(grepl("existing phospho.sqlite kept", msg)))
})

test_that("the builder skips up-to-date layers and the app refuses stores of another format", {
  msg <- capture_messages(build_store(tmp, log = FALSE)); expect_true(any(grepl("up to date", msg)))
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(tmp, "bulk.sqlite"))
  DBI::dbExecute(con, "UPDATE info SET value = ? WHERE key = 'format'", params = list(blob::as_blob(list(pack_rds(99L)))))
  DBI::dbDisconnect(con)
  expect_error(store_open(file.path(tmp, "bulk.sqlite")), "store format 99")
})

close_app_data(D)
unlink(tmp, recursive = TRUE)
