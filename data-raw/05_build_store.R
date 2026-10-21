# =====================================================================================================
# 05_build_store.R
# Converts the prepared tables written by 01-03 (cptac.rds, bulk_tcga_gtex.rds, sc_groups.rds, sc_expr.rds,
# sc_umap.rds) into the app's on-disk store: one SQLite file per data layer with one compressed record per
# gene, so the app reads a gene's values on demand instead of holding every matrix in memory.
#   cptac.sqlite  tables protein, rna       record = list(tumor, normal)   (log2 values, NA = not quantified)
#   bulk.sqlite   table  log2tpm_x1000      record = integer vector over samples
#   sc.sqlite     table  expr               record = list(mean = mean CP10k x100, frac = fraction x10000) over groups
#                 table  umap               record = log1p(CP10k) of the atlas-map subsample (sparse, 2 decimals)
#   phospho.sqlite table  site               record = list(site, peptide, tumor, normal): all phosphosites of one gene
#                                            (values flattened site-major-by-column, i.e. as.vector() of the sites x samples matrix)
#                                            optional layer (02b); info$gene_map links index symbols to phospho symbols
#   clinical.sqlite info only                 per-tumour clinical tables (survival, stage, Moffitt score) of CPTAC and TCGA;
#                                            optional layer (03b)
#   coexp.sqlite  table  partners           record = list(compartment, partner, r_x1000, q_mlog10_x100): the most co-expressed
#                                            genes of one gene in each cell type (indices into info$genes / info$compartments);
#                                            optional layer (01b)
# Small tables (sample metadata, column names, single-cell groups, UMAP coordinates) go into each file's
# `info` table. Every record is decoded again after writing and compared with the source (identical()).
# Run from the repository root after 01-04 (and 01b, 02b, 03b):  source("data-raw/05_build_store.R")
# get_data.R calls build_store() itself. Peak memory ~2 GB (single-cell layer); ~6 min.
# =====================================================================================================
suppressPackageStartupMessages({ library(Matrix); library(DBI); library(RSQLite) })

# reader + codec definitions shared with the app (STORE_FORMAT, store_unpack, store_open, store_get, ...)
if (!exists("store_unpack", mode = "function"))
  for (f in c("00_config.R", "01_data.R")) source(file.path(getOption("pdacx.root", "."), "app", "R", f))

# ---- encoders (decoders: store_unpack in app/R/01_data.R) ----
# exact scaled-integer form for double vectors with few decimals (e.g. CPTAC log2 values with 3 decimals)
store_scale <- function(x) {
  if (is.list(x) && !is.data.frame(x)) { a <- attributes(x); x <- lapply(x, store_scale); attributes(x) <- a; return(x) }
  if (!is.double(x) || !is.null(attributes(x)) || !length(x) || any(is.nan(x) | is.infinite(x))) return(x)
  for (sc in c(1, 10, 100, 1000, 10000)) {
    v <- round(x * sc); if (any(abs(v) > .Machine$integer.max, na.rm = TRUE)) return(x)
    v <- as.integer(v); if (identical(v / sc, x)) return(structure(list(v = v, scale = sc), class = "pdacx_scaled"))
  }
  x
}
pack_rds <- function(x, type = "gzip") memCompress(serialize(store_scale(x), NULL, xdr = FALSE), type)
pack_sparse16 <- function(i, x, n, type = "bzip2") {     # i: 1-based row indices (increasing), x: values with 2 decimals
  keep <- x != 0; i <- i[keep]; x <- x[keep]
  x100 <- round(x * 100); gaps <- diff(c(0L, i))
  if (n > 2147483647 || any(gaps > 65535) || any(x100 < 0 | x100 > 65535) || !identical(x100 / 100, x))
    stop("sparse16 codec cannot represent this column exactly")
  memCompress(c(writeBin(as.integer(c(n, length(i))), raw(), size = 4L, endian = "little"),
                writeBin(as.integer(gaps), raw(), size = 2L, endian = "little"),
                writeBin(as.integer(x100), raw(), size = 2L, endian = "little")), type)
}
codec_type <- function(codec) sub("^.*-", "", codec)

store_create <- function(path, layer, tables, info) {
  part <- paste0(path, ".part"); unlink(part)
  con <- DBI::dbConnect(RSQLite::SQLite(), part)
  DBI::dbExecute(con, "PRAGMA page_size = 8192"); DBI::dbExecute(con, "PRAGMA journal_mode = OFF"); DBI::dbExecute(con, "PRAGMA synchronous = OFF")
  DBI::dbExecute(con, "CREATE TABLE info (key TEXT PRIMARY KEY, value BLOB NOT NULL)")
  for (t in names(tables)) DBI::dbExecute(con, paste0("CREATE TABLE \"", t, "\" (id TEXT PRIMARY KEY, data BLOB NOT NULL)"))
  info <- c(list(format = STORE_FORMAT, layer = layer, tables = tables, built = format(Sys.time(), tz = "UTC", usetz = TRUE),
                 built_with = paste0(R.version.string, "; RSQLite ", packageVersion("RSQLite"))), info)
  s <- list(con = con, path = path, part = part); store_add_info(s, info); s
}
store_add_info <- function(s, info) DBI::dbExecute(s$con, "INSERT INTO info (key, value) VALUES (?, ?)",
                                                 params = list(names(info), blob::as_blob(lapply(info, pack_rds, type = codec_type(STORE_INFO_CODEC)))))
store_write <- function(s, table, ids, blobs) {
  DBI::dbWithTransaction(s$con, DBI::dbExecute(s$con, paste0("INSERT INTO \"", table, "\" (id, data) VALUES (?, ?)"),
                                               params = list(ids, blob::as_blob(blobs))))
}
store_finish <- function(s) { DBI::dbExecute(s$con, "VACUUM"); DBI::dbDisconnect(s$con); invisible(s$part) }
store_publish <- function(s) {     # only verified files get their final name (an interrupted build leaves only *.part)
  unlink(s$path); if (!file.rename(s$part, s$path)) stop("could not replace ", s$path, " (is the app still running and holding it open?)")
  invisible(s$path)
}
chunks <- function(n, size = 1000L) split(seq_len(n), ceiling(seq_len(n) / size))
check_ids <- function(ids, what) { if (anyNA(ids) || anyDuplicated(ids) || any(!nzchar(ids))) stop(what, ": row names must be unique and non-empty") }

# verify every record against ref(k), the source record of ids[k]; returns the number of genes checked
verify_table <- function(path, table, ids, ref) {
  st <- store_open(path); on.exit(store_close(st))
  stored <- store_ids(st, table); if (!setequal(stored, ids)) stop(basename(path), "/", table, ": stored ids differ from the source")
  for (ch in chunks(length(ids), 2000L)) {
    store_cache_clear()
    for (k in ch) if (!identical(store_get(st, table, ids[k]), ref(k))) stop(basename(path), "/", table, ": record for '", ids[k], "' differs from the source")
  }
  store_cache_clear(); length(ids)
}

# ---- layers ----
build_cptac <- function(dir, say) {
  src <- file.path(dir, "cptac.rds"); x <- readRDS(src); out <- file.path(dir, STORE_LAYERS[["cptac"]])
  tabs <- c(protein = "rds-gzip", rna = "rds-gzip")
  for (a in names(tabs)) {
    check_ids(rownames(x[[paste0(a, "_tumor")]]), a)
    if (!identical(rownames(x[[paste0(a, "_tumor")]]), rownames(x[[paste0(a, "_normal")]]))) stop("cptac ", a, ": tumour and normal rows differ")
  }
  cols <- lapply(x[c("protein_tumor", "protein_normal", "rna_tumor", "rna_normal")], colnames)
  s <- store_create(out, "cptac", tabs, list(colnames = cols, clinical = x$clinical, source = x$source, source_md5 = unname(tools::md5sum(src))))
  for (a in names(tabs)) {
    tu <- x[[paste0(a, "_tumor")]]; no <- x[[paste0(a, "_normal")]]; ids <- rownames(tu)
    for (k in chunks(length(ids))) store_write(s, a, ids[k], lapply(k, function(r) pack_rds(list(tumor = unname(tu[r, ]), normal = unname(no[r, ])), codec_type(tabs[[a]]))))
  }
  store_finish(s)
  for (a in names(tabs)) {
    tu <- x[[paste0(a, "_tumor")]]; no <- x[[paste0(a, "_normal")]]
    n <- verify_table(s$part, a, rownames(tu), function(k) list(tumor = unname(tu[k, ]), normal = unname(no[k, ])))
    say("cptac.sqlite/", a, ": ", n, " genes verified")
  }
  store_publish(s)
  out
}

build_bulk <- function(dir, say) {
  src <- file.path(dir, "bulk_tcga_gtex.rds"); x <- readRDS(src); out <- file.path(dir, STORE_LAYERS[["bulk"]])
  m <- x$log2tpm_x1000; ids <- rownames(m); check_ids(ids, "bulk")
  s <- store_create(out, "bulk", c(log2tpm_x1000 = "rds-gzip"),
                    list(colnames = colnames(m), samples = x$samples, signatures = x$signatures, source = x$source, source_md5 = unname(tools::md5sum(src))))
  for (k in chunks(length(ids))) store_write(s, "log2tpm_x1000", ids[k], lapply(k, function(r) pack_rds(unname(m[r, ]), "gzip")))
  store_finish(s)
  say("bulk.sqlite/log2tpm_x1000: ", verify_table(s$part, "log2tpm_x1000", ids, function(k) unname(m[k, ])), " genes verified")
  store_publish(s)
  out
}

build_sc <- function(dir, say) {      # one large matrix in memory at a time (peak ~1 GB)
  out <- file.path(dir, STORE_LAYERS[["sc"]])
  src <- stats::setNames(file.path(dir, c("sc_groups.rds", "sc_expr.rds", "sc_umap.rds")), c("groups", "expr", "umap"))
  has_umap <- file.exists(src[["umap"]])
  groups <- readRDS(src[["groups"]])
  read_expr <- function() {
    e <- readRDS(src[["expr"]]); check_ids(e$genes, "sc")
    if (!identical(rownames(e$mean_cp10k_x100), e$genes) || !identical(rownames(e$frac_x10000), e$genes) || ncol(e$mean_cp10k_x100) != nrow(groups))
      stop("sc_expr.rds does not match sc_groups.rds")
    e
  }
  e <- read_expr(); ids <- e$genes
  s <- store_create(out, "sc", c(expr = "rds-bzip2", if (has_umap) c(umap = "sparse16-bzip2")),
                    list(groups = groups, genes = ids, note = e$note, source_md5 = unname(tools::md5sum(src[file.exists(src)]))))
  for (k in chunks(length(ids))) store_write(s, "expr", ids[k], lapply(k, function(r) pack_rds(list(mean = unname(e$mean_cp10k_x100[r, ]), frac = unname(e$frac_x10000[r, ])), "bzip2")))
  rm(e); invisible(gc())
  if (has_umap) {
    u <- readRDS(src[["umap"]]); X <- u$expr; n <- nrow(X)
    if (!identical(colnames(X), ids)) stop("sc_umap.rds genes differ from sc_expr.rds")
    if (nrow(u$coords) != n || nrow(u$meta) != n) stop("sc_umap.rds: coordinates, metadata and expression differ in cell number")
    store_add_info(s, list(umap_coords = u$coords, umap_meta = u$meta, umap_note = u$note)); rm(u)
    col <- function(j) { a <- X@p[j]; b <- X@p[j + 1L]; if (b > a) list(i = X@i[(a + 1L):b] + 1L, x = X@x[(a + 1L):b]) else list(i = integer(0), x = numeric(0)) }
    for (k in chunks(length(ids))) store_write(s, "umap", ids[k], lapply(k, function(j) { cj <- col(j); pack_sparse16(cj$i, cj$x, n, "bzip2") }))
  }
  store_finish(s)
  if (has_umap) {
    say("sc.sqlite/umap: ", verify_table(s$part, "umap", ids, function(k) { cj <- col(k); v <- numeric(n); v[cj$i] <- cj$x; v }), " genes verified")
    rm(X); invisible(gc())
  }
  e <- read_expr()
  say("sc.sqlite/expr: ", verify_table(s$part, "expr", ids, function(k) list(mean = unname(e$mean_cp10k_x100[k, ]), frac = unname(e$frac_x10000[k, ]))), " genes verified")
  store_publish(s)
  out
}

# phosphosite symbols -> gene index: the CPTAC RNA/protein symbol of the same release first, then the official
# symbol, then an unambiguous alias. Returns data.frame(symbol = index symbol, sym_cptac_phospho); an index gene can have
# several phospho symbols (e.g. a retired symbol for another protein isoform), whose sites the app shows together.
phospho_gene_map <- function(sym, genes) {
  i <- match(sym, genes$sym_cptac_rna)
  j <- is.na(i); i[j] <- match(sym[j], genes$sym_cptac_protein)
  j <- is.na(i); i[j] <- match(sym[j], genes$symbol)
  j <- is.na(i)
  if (any(j)) {
    al <- strsplit(genes$aliases, ", ", fixed = TRUE); a <- unlist(al); ix <- rep(seq_along(al), lengths(al))
    k <- match(sym[j], a); k[sym[j] %in% a[duplicated(a)]] <- NA; i[j] <- ix[k]
  }
  m <- data.frame(symbol = genes$symbol[i], sym_cptac_phospho = sym, stringsAsFactors = FALSE)
  m <- m[!is.na(m$symbol), , drop = FALSE]; rownames(m) <- NULL; m
}

build_phospho <- function(dir, say) {
  src <- file.path(dir, "cptac_phospho.rds"); x <- readRDS(src); out <- file.path(dir, STORE_LAYERS[["phospho"]])
  sites <- x$sites; tu <- x$phospho_tumor; no <- x$phospho_normal
  check_ids(sites$site_id, "phospho sites")
  if (!identical(rownames(tu), sites$site_id) || !identical(rownames(no), sites$site_id)) stop("cptac_phospho.rds: site rows differ between tables")
  genes <- readRDS(file.path(dir, "genes.rds"))$genes
  rows <- split(seq_len(nrow(sites)), factor(sites$gene, levels = unique(sites$gene))); ids <- names(rows)
  map <- phospho_gene_map(ids, genes)
  dup <- map[map$symbol %in% map$symbol[duplicated(map$symbol)], ]
  say("phospho: ", length(ids), " gene symbols with sites, ", nrow(map), " linked to ", length(unique(map$symbol)), " index genes",
      if (length(un <- setdiff(ids, map$sym_cptac_phospho))) paste0("; not linked: ", paste(un, collapse = ", ")) else "",
      if (nrow(dup)) paste0("; combined: ", paste(tapply(dup$sym_cptac_phospho, dup$symbol, paste, collapse = " + "), collapse = ", ")) else "")
  rec <- function(k) { r <- rows[[k]]
    list(site = sites$site_id[r], peptide = sites$peptide[r], tumor = as.vector(unname(tu[r, , drop = FALSE])), normal = as.vector(unname(no[r, , drop = FALSE]))) }
  s <- store_create(out, "phospho", c(site = "rds-gzip"),
                    list(colnames = list(phospho_tumor = colnames(tu), phospho_normal = colnames(no)), gene_map = map, n_sites = nrow(sites),
                         source = x$source, source_md5 = unname(tools::md5sum(src)), genes_md5 = unname(tools::md5sum(file.path(dir, "genes.rds")))))
  for (k in chunks(length(ids))) store_write(s, "site", ids[k], lapply(k, function(j) pack_rds(rec(j), "gzip")))
  store_finish(s)
  say("phospho.sqlite/site: ", verify_table(s$part, "site", ids, rec), " genes (", nrow(sites), " sites) verified")
  store_publish(s)
  out
}

build_clinical <- function(dir, say) {          # small per-patient tables only: everything goes into info
  src <- file.path(dir, "clinical.rds"); x <- readRDS(src); out <- file.path(dir, STORE_LAYERS[["clinical"]])
  keys <- c("cptac", "tcga", "moffitt", "moffitt_used", "source")
  s <- store_create(out, "clinical", list(), c(x[keys], list(source_md5 = unname(tools::md5sum(src)))))
  store_finish(s)
  st <- store_open(s$part); ok <- vapply(keys, function(k) identical(st$info[[k]], x[[k]]), logical(1)); store_close(st)
  if (!all(ok)) stop("clinical.sqlite: info differs from the source for ", paste(keys[!ok], collapse = ", "))
  say("clinical.sqlite: ", nrow(x$cptac), " CPTAC and ", nrow(x$tcga), " TCGA tumours verified")
  store_publish(s)
  out
}

# co-expression partners (01b): one record per single-cell gene symbol with every cell type in which it was analysed
coexp_record <- function(tab) {
  o <- order(tab$compartment, -tab$r, tab$partner); tab <- tab[o, ]
  list(compartment = as.integer(tab$compartment), partner = as.integer(tab$partner), r_x1000 = as.integer(round(tab$r * 1000)),
       q_mlog10_x100 = as.integer(pmin(round(-log10(pmax(tab$q, 1e-300)) * 100), 30000L)))
}
build_coexp <- function(dir, say) {
  src <- file.path(dir, "sc_coexp.rds"); x <- readRDS(src); out <- file.path(dir, STORE_LAYERS[["coexp"]])
  rows <- split(seq_len(nrow(x$pairs)), x$pairs$gene); ids <- x$genes[as.integer(names(rows))]; check_ids(ids, "coexp")
  rec <- function(k) coexp_record(x$pairs[rows[[k]], c("compartment", "partner", "r", "q")])
  s <- store_create(out, "coexp", c(partners = "rds-gzip"),
                    list(genes = x$genes, compartments = x$compartments, cell_types = x$info, settings = x$settings, note = x$note,
                         source_md5 = unname(tools::md5sum(src)), built_from = x$built))
  for (k in chunks(length(ids))) store_write(s, "partners", ids[k], lapply(k, function(j) pack_rds(rec(j), "gzip")))
  store_finish(s)
  say("coexp.sqlite/partners: ", verify_table(s$part, "partners", ids, rec), " genes (", nrow(x$pairs), " pairs) verified")
  store_publish(s)
  out
}

# every gene of the search index must have a record under its dataset symbol
check_index <- function(dir, say) {
  g <- readRDS(file.path(dir, "genes.rds"))$genes
  need <- list(list("cptac", "protein", "sym_cptac_protein"), list("cptac", "rna", "sym_cptac_rna"),
               list("bulk", "log2tpm_x1000", "sym_tcga_gtex"), list("sc", "expr", "sym_sc"), list("sc", "umap", "sym_sc"))
  for (n in need) {
    path <- file.path(dir, STORE_LAYERS[[n[[1]]]]); if (!file.exists(path)) next
    st <- store_open(path); if (!n[[2]] %in% names(st$tables)) { store_close(st); next }
    miss <- setdiff(stats::na.omit(g[[n[[3]]]]), store_ids(st, n[[2]])); store_close(st)
    if (length(miss)) stop("gene index: ", length(miss), " ", n[[3]], " symbols have no record in ", basename(path), "/", n[[2]], " (e.g. ", miss[1], ")")
  }
  path <- file.path(dir, STORE_LAYERS[["phospho"]])
  if (file.exists(path)) {
    st <- store_open(path); map <- st$info$gene_map; ids <- store_ids(st, "site"); store_close(st)
    if (!all(map$symbol %in% g$symbol)) stop("phospho.sqlite: gene map refers to symbols missing from genes.rds; rebuild the phospho layer")
    if (!all(map$sym_cptac_phospho %in% ids)) stop("phospho.sqlite: gene map refers to genes without a record")
  }
  path <- file.path(dir, STORE_LAYERS[["coexp"]])
  if (file.exists(path)) {
    st <- store_open(path); sc_genes <- st$info$genes; store_close(st)
    if (!all(stats::na.omit(g$sym_sc) %in% sc_genes)) stop("coexp.sqlite: built from other single-cell genes than genes.rds; rebuild 01b and the coexp layer")
  }
  say("gene index: all ", nrow(g), " genes resolve to stored records")
}

# Build (or rebuild) the store in `dir`. layers: any of "cptac", "bulk", "sc", "phospho", "clinical", "coexp". Layers whose store
# file is newer than its source .rds files are skipped unless force = TRUE. The single-cell, phosphosite, clinical and
# co-expression layers are optional: without their source files they are skipped (an existing store file is kept).
build_store <- function(dir = Sys.getenv("PDACX_DATA_DIR", file.path("app", "data")), layers = c("cptac", "bulk", "sc", "phospho", "clinical", "coexp"), force = FALSE, log = TRUE) {
  log_file <- NULL
  if (log && file.exists(file.path("data-raw", "helpers.R")) && file.access("data-raw", 2) == 0) {
    source(file.path("data-raw", "helpers.R"), local = TRUE); log_file <- start_log("05_build_store")
  }
  say <- function(...) if (is.null(log_file)) message(...) else logc(log_file, ...)
  srcs <- list(cptac = "cptac.rds", bulk = "bulk_tcga_gtex.rds", sc = c("sc_groups.rds", "sc_expr.rds", "sc_umap.rds"), phospho = c("cptac_phospho.rds", "genes.rds"),
               clinical = "clinical.rds", coexp = "sc_coexp.rds")
  optional <- c("phospho", "clinical", "coexp")      # skipped without their source files (an existing store file is kept)
  for (l in layers) {
    out <- file.path(dir, STORE_LAYERS[[l]]); src <- file.path(dir, srcs[[l]]); src <- src[file.exists(src)]
    if (l == "sc" && !all(file.exists(file.path(dir, c("sc_groups.rds", "sc_expr.rds"))))) { say("sc: source files not found, skipped"); next }
    if (l %in% optional && !all(file.exists(file.path(dir, srcs[[l]])))) {
      say(l, ": source files not found, skipped", if (file.exists(out)) paste0(" (existing ", basename(out), " kept)") else ""); next }
    if (!length(src)) stop(l, ": source file(s) not found in ", dir)
    if (!force && file.exists(out) && all(file.mtime(out) > file.mtime(src))) {
      ok <- tryCatch({ st <- store_open(out); store_close(st); TRUE }, error = function(e) FALSE)
      if (ok) { say(basename(out), ": up to date, skipped"); next }
    }
    t0 <- Sys.time(); say(l, ": building ", basename(out), " from ", paste(basename(src), collapse = ", "), " ...")
    switch(l, cptac = build_cptac(dir, say), bulk = build_bulk(dir, say), sc = build_sc(dir, say), phospho = build_phospho(dir, say),
           clinical = build_clinical(dir, say), coexp = build_coexp(dir, say))
    invisible(gc())
    say(basename(out), ": ", round(file.size(out) / 2^20, 1), " MB, MD5 ", unname(tools::md5sum(out)), ", ", round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min")
  }
  if (file.exists(file.path(dir, "genes.rds"))) check_index(dir, say)
  if (!is.null(log_file)) end_log(log_file)
  invisible(file.path(dir, STORE_LAYERS))
}

# sourcing this file builds the store; get_data.R sets options(pdacx.build_store.autorun = FALSE) and calls build_store() itself
if (isTRUE(getOption("pdacx.build_store.autorun", TRUE))) build_store()
