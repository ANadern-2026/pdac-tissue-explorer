# ---- data loading and per-gene accessors ------------------------------------------------------------------
# Since v0.2 the large tables live on disk: one SQLite file per data layer (cptac.sqlite, bulk.sqlite, sc.sqlite),
# built by data-raw/05_build_store.R. Each file holds
#   info(key, value)  small R objects kept in memory (sample tables, column names, UMAP coordinates, notes, ...)
#   <table>(id, data) one compressed record per gene (id = the gene symbol used in that dataset)
# Only the gene index and the small tables are loaded at start-up; per-gene values are read on demand and kept
# in a small least-recently-used cache, so memory use no longer grows with the size of the data.
suppressPackageStartupMessages({ library(DBI); library(RSQLite) })

STORE_FORMAT  <- 1L      # bump when the on-disk layout changes (05_build_store.R writes the same constant)
STORE_INFO_CODEC <- "rds-gzip"
STORE_LAYERS  <- c(cptac = "cptac.sqlite", bulk = "bulk.sqlite", sc = "sc.sqlite", phospho = "phospho.sqlite")   # sc, phospho: optional
STORE_CACHE_MB <- 64     # memory budget of the per-gene record cache (shared by all sessions)

# ---- record codecs (encoders are in data-raw/05_build_store.R) ----
# "rds-<gzip|bzip2|xz>": a serialised R object; numeric vectors with few decimals may be stored as scaled integers
#                        (class pdacx_scaled, exact: the encoder only does this when decoding returns identical values)
# "sparse16-<...>":      a sparse numeric column: int32 length, int32 nnz, uint16 index gaps, uint16 values x100
store_unscale <- function(x) {
  if (inherits(x, "pdacx_scaled")) return(x$v / x$scale)
  if (is.list(x) && !is.data.frame(x)) { a <- attributes(x); x <- lapply(x, store_unscale); attributes(x) <- a }
  x
}
store_unpack <- function(b, codec) {
  kind <- sub("-.*$", "", codec); type <- sub("^.*-", "", codec)
  b <- memDecompress(b, type)
  switch(kind,
    "rds" = store_unscale(unserialize(b)),
    "sparse16" = {
      h <- readBin(b[1:8], "integer", 2L, size = 4L, endian = "little"); n <- h[1]; k <- h[2]
      v <- numeric(n)
      if (k > 0L) {
        i <- cumsum(readBin(b[8L + seq_len(2L * k)], "integer", k, size = 2L, signed = FALSE, endian = "little"))
        v[i] <- readBin(b[8L + 2L * k + seq_len(2L * k)], "integer", k, size = 2L, signed = FALSE, endian = "little") / 100
      }
      v
    },
    stop("Unknown record codec '", codec, "'"))
}

# ---- least-recently-used record cache ----
.store_cache <- new.env(parent = emptyenv())
.store_cache$keys <- character(0); .store_cache$bytes <- numeric(0)     # keys in order of last use, oldest first
store_cache_clear <- function() {
  rm(list = setdiff(ls(.store_cache, all.names = TRUE), c("keys", "bytes")), envir = .store_cache)
  .store_cache$keys <- character(0); .store_cache$bytes <- numeric(0)
}

store_open <- function(path, layer = sub("\\.sqlite$", "", basename(path))) {
  con <- DBI::dbConnect(RSQLite::SQLite(), path, flags = RSQLite::SQLITE_RO)
  inf <- DBI::dbGetQuery(con, "SELECT key, value FROM info")
  info <- stats::setNames(lapply(inf$value, store_unpack, codec = STORE_INFO_CODEC), inf$key)
  if (!identical(info$format, STORE_FORMAT)) {
    DBI::dbDisconnect(con)
    stop("'", path, "' has store format ", format(info$format), "; this app version needs format ", STORE_FORMAT,
         ". Rebuild it with data-raw/05_build_store.R (or run get_data.R).")
  }
  st <- new.env(parent = emptyenv())
  st$con <- con; st$path <- normalizePath(path); st$layer <- layer; st$info <- info; st$tables <- info$tables
  st
}

store_close <- function(st) if (is.environment(st) && DBI::dbIsValid(st$con)) DBI::dbDisconnect(st$con)

# one record (a vector or list) for one gene; NULL if the gene is not in that table
store_get <- function(st, table, id) {
  if (!table %in% names(st$tables)) stop("Table '", table, "' is not in store '", st$layer, "'")
  if (length(id) != 1L || is.na(id)) return(NULL)
  key <- paste(st$path, table, id, sep = "\r")
  hit <- .store_cache[[key]]
  if (!is.null(hit)) {                                   # move to the most-recent end
    j <- match(key, .store_cache$keys); b <- .store_cache$bytes[j]
    .store_cache$keys <- c(.store_cache$keys[-j], key); .store_cache$bytes <- c(.store_cache$bytes[-j], b)
    return(hit)
  }
  res <- DBI::dbGetQuery(st$con, paste0("SELECT data FROM \"", table, "\" WHERE id = ?"), params = list(id))
  if (!nrow(res)) return(NULL)
  val <- store_unpack(res$data[[1]], st$tables[[table]])
  assign(key, val, envir = .store_cache)
  .store_cache$keys <- c(.store_cache$keys, key); .store_cache$bytes <- c(.store_cache$bytes, as.numeric(utils::object.size(val)))
  over <- cumsum(rev(.store_cache$bytes)) > STORE_CACHE_MB * 2^20          # evict oldest records beyond the budget
  if (any(over[-1])) {                                                       # (the newest record always stays)
    keep <- rev(!over); keep[length(keep)] <- TRUE
    rm(list = .store_cache$keys[!keep], envir = .store_cache)
    .store_cache$keys <- .store_cache$keys[keep]; .store_cache$bytes <- .store_cache$bytes[keep]
  }
  val
}

store_ids <- function(st, table) DBI::dbGetQuery(st$con, paste0("SELECT id FROM \"", table, "\""))$id

# ---- app data ----
load_app_data <- function(dir = DATA_DIR) {
  f <- function(x) file.path(dir, x)
  need <- c("genes.rds", STORE_LAYERS[c("cptac", "bulk")])
  miss <- need[!file.exists(f(need))]
  if (length(miss)) {
    old <- file.exists(f("sc_expr.rds")) || file.exists(f("cptac.rds"))
    stop("Missing data files in '", dir, "': ", paste(miss, collapse = ", "), ". ",
         if (old) "The v0.1 .rds files are present: build the on-disk store with source(\"data-raw/05_build_store.R\") from the repository root."
         else "Run get_data.R from the repository root (or the scripts in data-raw/).")
  }
  d <- list(genes = readRDS(f("genes.rds"))$genes, dir = normalizePath(dir))
  if (file.exists(f("gene_xrefs.rds"))) d$genes <- add_xrefs(d$genes, readRDS(f("gene_xrefs.rds"))$xrefs)   # optional (links out)
  d$cptac <- store_open(f(STORE_LAYERS[["cptac"]]))
  d$bulk  <- store_open(f(STORE_LAYERS[["bulk"]]))
  d$bulk_samples <- d$bulk$info$samples
  d$bulk_samples$cohort <- factor(d$bulk_samples$cohort, levels = BULK_LEVELS)
  d$has_sc <- file.exists(f(STORE_LAYERS[["sc"]]))
  d$umap <- NULL
  if (d$has_sc) {
    d$sc <- store_open(f(STORE_LAYERS[["sc"]]))
    d$sc_groups <- prepare_sc_groups(d$sc$info$groups)
    d$sc_note <- d$sc$info$note
    if ("umap" %in% names(d$sc$tables)) {
      d$umap <- list(coords = d$sc$info$umap_coords, meta = d$sc$info$umap_meta, note = d$sc$info$umap_note)
      d$umap$meta$compartment <- compartment_of(d$umap$meta$cell_type, d$umap$meta$fine_group)
    }
  }
  d$has_phospho <- file.exists(f(STORE_LAYERS[["phospho"]]))
  if (d$has_phospho) { d$phospho <- store_open(f(STORE_LAYERS[["phospho"]])); d$phospho_map <- d$phospho$info$gene_map }
  d$choices <- stats::setNames(d$genes$symbol, ifelse(is.na(d$genes$name), d$genes$symbol, paste0(d$genes$symbol, " - ", d$genes$name)))
  d
}

close_app_data <- function(d) { for (k in c("cptac", "bulk", "sc", "phospho")) store_close(d[[k]]); store_cache_clear(); invisible(NULL) }

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

# Ensembl ID and gene type (data-raw/04b_gene_xrefs.R) joined by symbol; used only where the Entrez ID agrees
add_xrefs <- function(genes, x) {
  i <- match(genes$symbol, x$symbol); ok <- !is.na(i) & (x$entrez[i] %in% genes$entrez | (is.na(x$entrez[i]) & is.na(genes$entrez)))
  genes$ensembl <- ifelse(ok, x$ensembl[i], NA_character_); genes$gene_type <- ifelse(ok, x$gene_type[i], NA_character_)
  genes
}

gene_row <- function(d, symbol) d$genes[match(symbol, d$genes$symbol), , drop = FALSE]

# CPTAC paired values for one gene and assay ("protein" | "rna")
cptac_pairs <- function(d, symbol, assay = c("protein", "rna")) {
  assay <- match.arg(assay); gr <- gene_row(d, symbol)
  s <- gr[[if (assay == "protein") "sym_cptac_protein" else "sym_cptac_rna"]]
  if (is.na(s)) return(NULL)
  rec <- store_get(d$cptac, assay, s); if (is.null(rec)) return(NULL)
  cols <- d$cptac$info$colnames
  t <- stats::setNames(rec$tumor, cols[[paste0(assay, "_tumor")]]); n <- stats::setNames(rec$normal, cols[[paste0(assay, "_normal")]])
  ids <- intersect(names(t), names(n)); ok <- ids[is.finite(t[ids]) & is.finite(n[ids])]
  tumours <- t[is.finite(t)]
  list(pairs = data.frame(case_id = ok, normal = unname(n[ok]), tumour = unname(t[ok]), diff = unname(t[ok] - n[ok])),
       n_tumours_measured = length(tumours), n_tumours_total = length(t), tumour_values = tumours, assay = assay, symbol_used = s)
}

# CPTAC phosphosites of one gene: sites x samples matrices (log2) for tumours and adjacent normals, the gene's protein
# level aligned to the same samples (for values relative to protein), and one annotation row per site; NULL if none
phospho_sites <- function(d, symbol) {
  if (!isTRUE(d$has_phospho)) return(NULL)
  ps <- d$phospho_map$sym_cptac_phospho[d$phospho_map$symbol %in% symbol]
  recs <- lapply(ps, function(s) store_get(d$phospho, "site", s)); have <- !vapply(recs, is.null, logical(1))
  ps <- ps[have]; recs <- recs[have]
  if (!length(recs)) return(NULL)
  cols <- d$phospho$info$colnames
  mat <- function(part, cn) { m <- do.call(rbind, lapply(recs, function(r) matrix(r[[part]], nrow = length(r$site)))); dimnames(m) <- list(NULL, cn); m }
  site_id <- unlist(lapply(recs, `[[`, "site"))
  re <- "^((?:NP|XP)_[0-9]+\\.[0-9]+)_([STY])([0-9]+)$"
  sites <- data.frame(site_id = site_id, refseq = sub(re, "\\1", site_id, perl = TRUE), residue = sub(re, "\\2", site_id, perl = TRUE),
                      position = as.integer(sub(re, "\\3", site_id, perl = TRUE)), peptide = unlist(lapply(recs, `[[`, "peptide")),
                      source_symbol = rep(ps, vapply(recs, function(r) length(r$site), integer(1))),
                      stringsAsFactors = FALSE)
  tu <- mat("tumor", cols$phospho_tumor); no <- mat("normal", cols$phospho_normal)
  ref_n <- table(sites$refseq); main <- names(ref_n)[order(-ref_n, names(ref_n))]       # isoform with most sites first
  o <- order(match(sites$refseq, main), sites$position, sites$residue); sites <- sites[o, ]; tu <- tu[o, , drop = FALSE]; no <- no[o, , drop = FALSE]
  sites$label <- paste0(sites$residue, sites$position, if (length(main) > 1) paste0(" (", sites$refseq, ")") else "")
  rownames(sites) <- NULL; rownames(tu) <- rownames(no) <- sites$label
  ps_prot <- gene_row(d, symbol)$sym_cptac_protein
  prot <- if (is.na(ps_prot)) NULL else store_get(d$cptac, "protein", ps_prot)
  pcols <- d$cptac$info$colnames
  align <- function(v, cn, target) if (is.null(v)) rep(NA_real_, length(target)) else unname(v[match(target, cn)])
  list(sites = sites, tumor = tu, normal = no, has_protein = !is.null(prot),
       prot_tumor = align(prot$tumor, pcols$protein_tumor, colnames(tu)), prot_normal = align(prot$normal, pcols$protein_normal, colnames(no)))
}

# matched tumour/normal values of site k (row index), as in cptac_pairs()$pairs; relative = site minus protein level
phospho_pairs <- function(ph, k, relative = FALSE) {
  t <- ph$tumor[k, ]; n <- ph$normal[k, ]
  if (relative) { t <- t - ph$prot_tumor; n <- n - ph$prot_normal }
  ids <- intersect(names(t), names(n)); ok <- ids[is.finite(t[ids]) & is.finite(n[ids])]
  data.frame(case_id = ok, normal = unname(n[ok]), tumour = unname(t[ok]), diff = unname(t[ok] - n[ok]))
}

bulk_values <- function(d, symbol) {
  s <- gene_row(d, symbol)$sym_tcga_gtex
  if (is.na(s)) return(NULL)
  v <- store_get(d$bulk, "log2tpm_x1000", s); if (is.null(v)) return(NULL)
  data.frame(d$bulk_samples, value = stats::setNames(v, d$bulk$info$colnames) / 1000, check.names = FALSE)
}

# single-cell: one row per (sample x fine group) for a gene
sc_values <- function(d, symbol) {
  if (!isTRUE(d$has_sc)) return(NULL)
  s <- gene_row(d, symbol)$sym_sc
  if (is.na(s)) return(NULL)
  rec <- store_get(d$sc, "expr", s); if (is.null(rec)) return(NULL)
  data.frame(d$sc_groups, mean_cp10k = rec$mean / 100, frac = rec$frac / 1e4)
}

# atlas map: log1p(CP10k) of one gene in every subsampled cell (order of d$umap$coords); NULL if unavailable
umap_values <- function(d, sym_sc) {
  if (is.null(d$umap) || length(sym_sc) != 1L || is.na(sym_sc)) return(NULL)
  store_get(d$sc, "umap", sym_sc)
}

# exact weighted aggregation of fine groups to any coarser grouping (per sample)
sc_aggregate <- function(v, by) {
  v$w_mean <- v$mean_cp10k * v$n_cells; v$w_frac <- v$frac * v$n_cells
  a <- stats::aggregate(cbind(n_cells, w_mean, w_frac) ~ ., data = v[, c(by, "n_cells", "w_mean", "w_frac")], FUN = sum)
  a$mean_cp10k <- a$w_mean / a$n_cells; a$frac <- a$w_frac / a$n_cells; a$avg_expr <- log1p(a$mean_cp10k)
  a[, c(by, "n_cells", "mean_cp10k", "avg_expr", "frac")]
}
