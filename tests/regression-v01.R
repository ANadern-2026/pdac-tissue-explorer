# Regression check for the on-disk store (v0.2 data layer) against the v0.1 in-memory data layer.
# Needs BOTH the v0.1 .rds files and the store files (built by data-raw/05_build_store.R) in app/data/.
# Compares, for every gene of the search index, the values the app works with (CPTAC pairs, TCGA/GTEx values,
# single-cell group values, atlas-map values) and, for a random subset plus the example genes, all derived
# statistics, summary sentences, multi-gene tables and the data behind every figure.
# Run from the repository root:  Rscript tests/regression-v01.R      (~50 min, ~3 GB RAM; not part of test_dir())
#                                 Rscript tests/regression-v01.R --derived-only    (skips the all-gene loop, ~5 min)
suppressPackageStartupMessages({ library(Matrix); library(ggplot2) })
DIR <- Sys.getenv("PDACX_DATA_DIR", file.path("app", "data"))
for (f in sort(list.files(file.path("app", "R"), pattern = "\\.R$", full.names = TRUE))) source(f)

# ---- v0.1 reference implementation (verbatim logic of app/R/01_data.R in v0.1.0) ----
v01_load <- function(dir) {
  f <- function(x) file.path(dir, x)
  d <- list(genes = readRDS(f("genes.rds"))$genes, cptac = readRDS(f("cptac.rds")), bulk = readRDS(f("bulk_tcga_gtex.rds")))
  d$has_sc <- TRUE; d$sc_groups <- prepare_sc_groups(readRDS(f("sc_groups.rds")))
  e <- readRDS(f("sc_expr.rds")); d$sc_mean <- e$mean_cp10k_x100; d$sc_frac <- e$frac_x10000
  d$umap <- readRDS(f("sc_umap.rds")); d$umap$meta$compartment <- compartment_of(d$umap$meta$cell_type, d$umap$meta$fine_group)
  d$bulk$samples$cohort <- factor(d$bulk$samples$cohort, levels = BULK_LEVELS)
  d
}
v01_cptac_pairs <- function(d, symbol, assay) {
  gr <- gene_row(d, symbol); s <- gr[[if (assay == "protein") "sym_cptac_protein" else "sym_cptac_rna"]]
  if (is.na(s)) return(NULL)
  t <- d$cptac[[paste0(assay, "_tumor")]][s, ]; n <- d$cptac[[paste0(assay, "_normal")]][s, ]
  ids <- intersect(names(t), names(n)); ok <- ids[is.finite(t[ids]) & is.finite(n[ids])]
  tumours <- t[is.finite(t)]
  list(pairs = data.frame(case_id = ok, normal = unname(n[ok]), tumour = unname(t[ok]), diff = unname(t[ok] - n[ok])),
       n_tumours_measured = length(tumours), n_tumours_total = length(t), tumour_values = tumours, assay = assay, symbol_used = s)
}
v01_bulk_values <- function(d, symbol) {
  s <- gene_row(d, symbol)$sym_tcga_gtex; if (is.na(s)) return(NULL)
  data.frame(d$bulk$samples, value = d$bulk$log2tpm_x1000[s, ] / 1000, check.names = FALSE)
}
v01_sc_values <- function(d, symbol) {
  s <- gene_row(d, symbol)$sym_sc; if (is.na(s)) return(NULL)
  i <- match(s, rownames(d$sc_mean))
  data.frame(d$sc_groups, mean_cp10k = d$sc_mean[i, ] / 100, frac = d$sc_frac[i, ] / 1e4)
}
v01_umap <- function(d, symbol) { s <- gene_row(d, symbol)$sym_sc; if (is.na(s)) NULL else as.numeric(d$umap$expr[, s]) }
# v0.1 multi-gene tables called the v0.1 accessors; run them against the old data with the old accessors in scope
v01_scope <- function(fun, d) { environment(fun) <- list2env(list(cptac_pairs = v01_cptac_pairs, bulk_values = v01_bulk_values, sc_values = v01_sc_values), parent = environment(fun)); fun }

old <- v01_load(DIR); new <- load_app_data(DIR)
t0 <- Sys.time(); bad <- character(0)
chk <- function(ok, what) if (!isTRUE(ok)) { bad <<- c(bad, what); message("DIFFERS: ", what) }
stopifnot(identical(old$genes, new$genes), identical(old$sc_groups, new$sc_groups), identical(old$bulk$samples, new$bulk_samples),
          identical(old$umap$coords, new$umap$coords), identical(old$umap$meta, new$umap$meta))

# 1. primary values, every gene
genes <- new$genes$symbol
for (k in if ("--derived-only" %in% commandArgs(TRUE)) integer(0) else seq_along(genes)) {
  g <- genes[k]
  chk(identical(v01_cptac_pairs(old, g, "protein"), cptac_pairs(new, g, "protein")), paste(g, "cptac protein"))
  chk(identical(v01_cptac_pairs(old, g, "rna"), cptac_pairs(new, g, "rna")), paste(g, "cptac rna"))
  chk(identical(v01_bulk_values(old, g), bulk_values(new, g)), paste(g, "bulk"))
  chk(identical(v01_sc_values(old, g), sc_values(new, g)), paste(g, "sc"))
  chk(identical(v01_umap(old, g), umap_values(new, gene_row(new, g)$sym_sc)), paste(g, "umap"))
  if (k %% 5000 == 0) message(k, " / ", length(genes), " genes compared")
}
if (!"--derived-only" %in% commandArgs(TRUE)) message("primary values: ", length(genes), " genes compared, ", length(bad), " differences")

# 2. derived statistics, text and figure data on a subset
set.seed(7)
sub <- unique(c("MSLN", "KRT17", "FAP", "GATA6", "CEACAM5", "COL11A1", "CD274", "KRT19", sample(genes[new$genes$n_datasets == 4], 150), sample(genes, 50)))
sub <- intersect(sub, genes)
fig_data <- function(p) { set.seed(1); if (inherits(p, "ggplot")) ggplot_build(p)$data else p }   # same seed: jittered points are random
for (g in sub) {
  po <- v01_cptac_pairs(old, g, "protein"); pn <- cptac_pairs(new, g, "protein"); ro <- v01_cptac_pairs(old, g, "rna"); rn <- cptac_pairs(new, g, "rna")
  bo <- v01_bulk_values(old, g); bn <- bulk_values(new, g); so <- v01_sc_values(old, g); sn <- sc_values(new, g)
  chk(identical(paired_stats(po), paired_stats(pn)) && identical(paired_stats(ro), paired_stats(rn)), paste(g, "paired_stats"))
  for (ctx in c("Primary tumour", "Liver metastasis")) chk(identical(if (is.null(so)) NULL else celltype_summary(so, ctx), if (is.null(sn)) NULL else celltype_summary(sn, ctx)), paste(g, "celltype_summary", ctx))
  cso <- if (is.null(so)) NULL else celltype_summary(so); csn <- if (is.null(sn)) NULL else celltype_summary(sn)
  chk(identical(list(txt_cptac(po, paired_stats(po), g, "protein"), txt_cptac(ro, paired_stats(ro), g, "rna"), txt_bulk(bo, g), txt_sc(cso, g), txt_state(so, g)),
                list(txt_cptac(pn, paired_stats(pn), g, "protein"), txt_cptac(rn, paired_stats(rn), g, "rna"), txt_bulk(bn, g), txt_sc(csn, g), txt_state(sn, g))), paste(g, "summary text"))
  figs_o <- list(plot_cptac_paired(po, g, "protein"), plot_cptac_paired(ro, g, "RNA"), plot_bulk(bo, g), plot_celltypes(cso, g, "Primary tumour"), plot_share(cso, g),
                 plot_frac(cso, g), plot_states(so, g, CELLTYPE_ORDER), plot_fibro(so, g), plot_ductal(so, g))
  figs_n <- list(plot_cptac_paired(pn, g, "protein"), plot_cptac_paired(rn, g, "RNA"), plot_bulk(bn, g), plot_celltypes(csn, g, "Primary tumour"), plot_share(csn, g),
                 plot_frac(csn, g), plot_states(sn, g, CELLTYPE_ORDER), plot_fibro(sn, g), plot_ductal(sn, g))
  for (i in seq_along(figs_o)) chk(identical(fig_data(figs_o[[i]]), fig_data(figs_n[[i]])), paste(g, "figure", i))
  uo <- old$umap; uo_val <- v01_umap(old, g)
  chk(identical(fig_data(plot_umap_static(uo, g, uo_val)), fig_data(plot_umap_static(new$umap, g, umap_values(new, gene_row(new, g)$sym_sc)))), paste(g, "atlas map figure"))
}
message("derived outputs: ", length(sub), " genes compared")

# 3. multi-gene tables (example lists + a random list of 40)
lists <- c(lapply(c("COL1A1, DCN, LUM, FAP, ACTA2, PDGFRB, CXCL12", "GATA6, TFF1, LGALS4, CLDN18, KRT17, KRT5, S100A2, VGLL1",
                    "CD274, PDCD1LG2, PDCD1, CTLA4, LAG3, HAVCR2, TIGIT"), function(x) parse_gene_list(x, new)$genes), list(sample(genes, 40)))
for (gl in lists) {
  chk(identical(v01_scope(multi_dot_table, old)(old, gl), multi_dot_table(new, gl)), paste("multi_dot_table:", paste(utils::head(gl, 3), collapse = ",")))
  chk(identical(v01_scope(multi_forest_table, old)(old, gl), multi_forest_table(new, gl)), paste("multi_forest_table:", paste(utils::head(gl, 3), collapse = ",")))
}
close_app_data(new)
message(if (length(bad)) paste("FAILED:", length(bad), "differences") else "PASSED: the store reproduces every v0.1 output that was compared",
        " (", round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min)")
if (length(bad)) quit(status = 1)
