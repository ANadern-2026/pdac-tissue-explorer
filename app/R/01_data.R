# ---- data loading and per-gene accessors ------------------------------------------------------------------
suppressPackageStartupMessages(library(Matrix))
load_app_data <- function(dir = DATA_DIR) {
  f <- function(x) file.path(dir, x)
  need <- c("genes.rds", "cptac.rds", "bulk_tcga_gtex.rds")
  miss <- need[!file.exists(f(need))]
  if (length(miss)) stop("Missing data files in '", dir, "': ", paste(miss, collapse = ", "), ". Run the scripts in data-raw/ first.")
  d <- list(genes = readRDS(f("genes.rds"))$genes, cptac = readRDS(f("cptac.rds")), bulk = readRDS(f("bulk_tcga_gtex.rds")))
  d$has_sc <- all(file.exists(f(c("sc_groups.rds", "sc_expr.rds"))))
  if (d$has_sc) {
    d$sc_groups <- prepare_sc_groups(readRDS(f("sc_groups.rds")))
    e <- readRDS(f("sc_expr.rds")); d$sc_mean <- e$mean_cp10k_x100; d$sc_frac <- e$frac_x10000
    d$sc_note <- e$note
  }
  d$umap <- if (file.exists(f("sc_umap.rds"))) readRDS(f("sc_umap.rds")) else NULL
  if (!is.null(d$umap)) d$umap$meta$compartment <- compartment_of(d$umap$meta$cell_type, d$umap$meta$fine_group)
  d$bulk$samples$cohort <- factor(d$bulk$samples$cohort, levels = BULK_LEVELS)
  d$choices <- setNames(d$genes$symbol, ifelse(is.na(d$genes$name), d$genes$symbol, paste0(d$genes$symbol, " - ", d$genes$name)))
  d
}

compartment_of <- function(cell_type, fine_group) {
  key <- ifelse(cell_type == "DUCTAL", ifelse(fine_group %in% names(DUCTAL_GROUP), "TUMOUR_DUCTAL", "OTHER_DUCTAL"), cell_type)
  factor(unname(CELLTYPE_LABELS[key]), levels = CELLTYPE_ORDER)
}

prepare_sc_groups <- function(g) {
  g$compartment <- compartment_of(g$cell_type, g$fine_group)
  g$state <- with(g, ifelse(disease_state == "Donor", "Healthy donor",
                     ifelse(disease_state == "Adjacent normal", ifelse(treatment == "Treated", NA, "Adjacent normal"),
                     ifelse(disease_state == "Primary tumor", ifelse(treatment == "Treated", "Primary tumour (treated)", ifelse(treatment == "Treatment-naive", "Primary tumour", NA)),
                     ifelse(disease_state == "Metastatic lesion" & met_site %in% "Liver" & treatment != "Treated", "Liver metastasis", NA)))))
  g$state <- factor(g$state, levels = STATE_LEVELS)
  g$ductal_group <- factor(unname(DUCTAL_GROUP[g$fine_group]), levels = DUCTAL_GROUP_ORDER)
  g$fibro_subtype <- factor(ifelse(g$fine_group %in% FIBRO_ORDER, g$fine_group, NA), levels = FIBRO_ORDER)
  g$fibro_class <- factor(ifelse(grepl("^myFb", g$fine_group), "Myofibroblasts (myFb)", ifelse(grepl("^iFb", g$fine_group), "Inflammatory fibroblasts (iFb)", NA)),
                          levels = c("Myofibroblasts (myFb)", "Inflammatory fibroblasts (iFb)"))
  g$study_label <- unname(STUDY_CITATION[g$study])
  g
}

gene_row <- function(d, symbol) d$genes[match(symbol, d$genes$symbol), , drop = FALSE]

# CPTAC paired values for one gene and assay ("protein" | "rna")
cptac_pairs <- function(d, symbol, assay = c("protein", "rna")) {
  assay <- match.arg(assay); gr <- gene_row(d, symbol)
  s <- gr[[if (assay == "protein") "sym_cptac_protein" else "sym_cptac_rna"]]
  if (is.na(s)) return(NULL)
  t <- d$cptac[[paste0(assay, "_tumor")]][s, ]; n <- d$cptac[[paste0(assay, "_normal")]][s, ]
  ids <- intersect(names(t), names(n)); ok <- ids[is.finite(t[ids]) & is.finite(n[ids])]
  tumours <- t[is.finite(t)]
  list(pairs = data.frame(case_id = ok, normal = unname(n[ok]), tumour = unname(t[ok]), diff = unname(t[ok] - n[ok])),
       n_tumours_measured = length(tumours), n_tumours_total = length(t), tumour_values = tumours, assay = assay, symbol_used = s)
}

bulk_values <- function(d, symbol) {
  s <- gene_row(d, symbol)$sym_tcga_gtex
  if (is.na(s)) return(NULL)
  data.frame(d$bulk$samples, value = d$bulk$log2tpm_x1000[s, ] / 1000, check.names = FALSE)
}

# single-cell: one row per (sample x fine group) for a gene
sc_values <- function(d, symbol) {
  if (!isTRUE(d$has_sc)) return(NULL)
  s <- gene_row(d, symbol)$sym_sc
  if (is.na(s)) return(NULL)
  i <- match(s, rownames(d$sc_mean))
  data.frame(d$sc_groups, mean_cp10k = d$sc_mean[i, ] / 100, frac = d$sc_frac[i, ] / 1e4)
}

# exact weighted aggregation of fine groups to any coarser grouping (per sample)
sc_aggregate <- function(v, by) {
  v$w_mean <- v$mean_cp10k * v$n_cells; v$w_frac <- v$frac * v$n_cells
  a <- stats::aggregate(cbind(n_cells, w_mean, w_frac) ~ ., data = v[, c(by, "n_cells", "w_mean", "w_frac")], FUN = sum)
  a$mean_cp10k <- a$w_mean / a$n_cells; a$frac <- a$w_frac / a$n_cells; a$avg_expr <- log1p(a$mean_cp10k)
  a[, c(by, "n_cells", "mean_cp10k", "avg_expr", "frac")]
}
