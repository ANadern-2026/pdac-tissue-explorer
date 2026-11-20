# =====================================================================================================
# 07_cell_populations.R
# Builds the two-level cell-population layer behind the "Cell populations" page, one compartment at a time
# (phase 1: myeloid; later phases add compartments to COMPARTMENTS and populations.csv).
#
#   Level 1 (lineages)      PUBLISHED labels where they exist (level1 = "published"), not re-derived here; for myeloid cells the atlas authors' subclustering
#                           (Loveless et al. 2025; Myeloid_Subcluster.txt.gz: Macrophage, Monocyte, Dendritic,
#                           Neutrophil) plus the atlas's coarse MAST type. Each label is checked against the cited
#                           marker set of its lineage (data-raw/reference/populations.csv) and the check is shipped
#                           with the data, so the app can show how well the published labels agree with the markers.
#   Level 2 (sub-lineages)  ASSIGNED HERE. Within each lineage: variable genes, PCA, Harmony by sample, kNN graph,
#                           Leiden clusters. Every cluster is named after the literature signature (populations.csv,
#                           label_method = "signature") with the highest mean module score, if that score is
#                           >= MIN_SCORE and exceeds the runner-up by >= MIN_MARGIN; otherwise the cluster is
#                           "unresolved". Programmes with label_method = "score_only" (neutrophil TAN states) are
#                           scored and shown, never used as labels.
#   Level 1 by signature    For compartments without an established lineage nomenclature (fibroblasts; level1 = "signature"):
#                           the whole compartment is clustered and every cluster named after the level-1 signature
#                           with the highest score (same rules as level 2). Candidates are the level-1 labels and the
#                           'boundary' populations (another cell type inside the atlas type, e.g. mural or Schwann cells:
#                           shown, but kept out of the compartment's make-up). The atlas authors' labels, if any, are
#                           kept per cell only for comparison (atlas_compare).
#   Compartment-wide        level-1 populations with label_method "score_only" (e.g. dissociation stress) are scored
#   programmes              on every cell, reported with every cluster ("strongest programme") and can veto a label:
#                           cfg$veto = c(<label> = <programme>) leaves a cluster unresolved if the programme scores
#                           higher than the label's own signature (currently used by no compartment; see decision FIB_04).
# Every population therefore carries its label source (published / signature), its signature genes and the
# references that define it. 07 refuses to run if a population cites a reference missing from references.csv.
#
# Input : data-raw/source/scatlas/sc_ref_base.rds   filtered atlas object from 00_scatlas_readiness.R (step A only)
#         data-raw/source/scatlas/Myeloid_Subcluster.txt.gz   https://github.com/PDAC-MULTIOMIC/PDAC_Atlas
#         data-raw/reference/populations.csv, data-raw/reference/references.csv (decisions.csv is checked here and
#         attached to the store by 05_build_store.R)
# Output: data-raw/work/pop_<id>_base.rds      step-A cache per compartment: its cells only (log1p CP10k, metadata)
#         app/data/pop_<id>.rds                 everything the app needs (05_build_store.R turns it into pop.sqlite)
#         data-raw/work/pop_<id>_cell_labels.tsv.gz   level-1/level-2 label and cluster of every cell
#         data-raw/logs/07_cell_populations_<date>.txt   provenance (checksums, counts, cluster calls, sessionInfo)
#
# Two steps. A is heavy (~50 GB RAM, as 01) and runs once; B reads the small cache (~3-6 GB RAM, ~20-40 min).
# After editing populations.csv (signatures, thresholds), re-run with RUN_STEP_A <- FALSE.
# Run from the repository root:  source("data-raw/07_cell_populations.R")   (PDACX_COMPARTMENT=<id> selects the compartment)
# Requires: Matrix, irlba, uwot, igraph (Leiden), harmony; SeuratObject for step A.
# =====================================================================================================

## ---- paths and settings --------------------------------------------------------------------------------
ATLAS_RDS   <- "data-raw/source/scatlas/sc_ref_base.rds"
POP_FILE    <- "data-raw/reference/populations.csv"
REF_FILE    <- "data-raw/reference/references.csv"
DEC_FILE    <- "data-raw/reference/decisions.csv"     # why each label is (or is not) used; attached by 05_build_store.R, checked here
WORK_DIR    <- file.path("data-raw", "work")
OUT_DIR     <- file.path("app", "data")
LOG_DIR     <- file.path("data-raw", "logs")

# One entry per compartment. Phases: add a compartment here and its populations to populations.csv, then run 07
# for it; nothing else changes. Level 1 comes from published labels (label_file, label_map) and/or from the atlas's
# coarse types (type_map); cells of the compartment's atlas types without either become `unlabelled`.
#   label_format "space": space-delimited file, barcode first, label `label_from_end` fields from the end
#   label_format "tab"  : tab-delimited with row names = barcodes, label in column `label_col` (as the ductal and
#                         fibroblast files read by 01_scatlas_summaries.R)
COMPARTMENTS <- list(
  myeloid = list(name = "Myeloid", atlas_types = c("MYELOID", "MAST"), level1 = "published",
                 label_file = "data-raw/source/scatlas/Myeloid_Subcluster.txt.gz", label_format = "space", label_from_end = 3L,
                 label_map = c(Macrophage = "MYE_MAC", Monocyte = "MYE_MONO", Dendritic = "MYE_DC", Neutrophil = "MYE_NEU"),
                 type_map = c(MAST = "MYE_MAST"), unlabelled = "MYE_UNLABELLED"),
  # fibroblasts: no lineage-based human nomenclature, so level 1 is assigned here by signature (decision FIB_01); the atlas's
  # myFb 1-5 / iFb 1-3 labels (column Coarse) are kept only for comparison
  fibroblasts = list(name = "Fibroblasts", atlas_types = "FIBROBLASTS", level1 = "signature",
                     label_file = "data-raw/source/scatlas/Fibroblast_Subcluster.txt.gz", label_format = "tab", label_col = "Coarse",
                     unlabelled = NA_character_)   # no veto: the published iCAF signature does not follow the stress genes in this atlas (decision FIB_04)
  # Template for a later phase (populations for it must be added to populations.csv first), e.g. T and NK cells,
  # which have no published sub-labels in the atlas, so level 1 would come from signatures as well:
  # tnk = list(name = "T and NK cells", atlas_types = "TNK", label_file = NULL, type_map = c(TNK = "TNK_ALL"), unlabelled = NA)
)
COMPARTMENT_ID <- Sys.getenv("PDACX_COMPARTMENT", "myeloid")      # which compartment step B builds
cfg <- COMPARTMENTS[[COMPARTMENT_ID]]; if (is.null(cfg)) stop("unknown compartment '", COMPARTMENT_ID, "'")
COMPARTMENT <- cfg$name; UNLABELLED <- cfg$unlabelled
SIG_L1 <- identical(cfg$level1, "signature")                       # level 1 assigned here (else: published labels)
VETO   <- if (is.null(cfg$veto)) character(0) else cfg$veto
cache_file  <- function(id) file.path(WORK_DIR, paste0("pop_", id, "_base.rds"))
STEP_A_TODO <- names(COMPARTMENTS)[!file.exists(vapply(names(COMPARTMENTS), cache_file, ""))]   # compartments not cached yet
RUN_STEP_A  <- length(STEP_A_TODO) > 0                             # one atlas read caches all of them

MIN_FEATURES   <- 200      # as 01_scatlas_summaries.R
MIN_GENE_CELLS <- 50       # as 01: same gene universe as the rest of the single-cell layer
N_HVG       <- 2000
N_PC        <- 30
HARMONY_BY  <- "sample"
N_NEIGHBORS <- 30
LEIDEN_RES  <- 1.0
MIN_CLUSTER <- 50          # smaller clusters stay unresolved
MIN_CLUSTER_SAMPLES <- 3L  # clusters from fewer samples stay unresolved (one patient's cells, e.g. FIB_MYO_c10 of 2026-11-13: 57 cells, 1 sample)
MIN_SCORE   <- 0.15        # cluster mean module score (log1p CP10k units) needed to take a signature's name ...
                           # (0.10 in the first run; raised after weak calls near 0.1 proved to be dominated by other programmes)
MIN_MARGIN  <- 0.05        # ... and its lead over the runner-up signature
LOW_GENES_FRAC <- 0.5      # clusters whose cells detect (median) fewer than this share of the genes typical of their own study stay unresolved
                           # (relative to the study, not the compartment: single-nucleus studies detect fewer genes per cell
                           # without being of low quality; first fibroblast run 2026-11-13 lost GSE202051 this way)
MIN_CELLS   <- 20L         # per-sample summaries: as the app's MIN_CELLS (stored for all groups; the app filters)
MARKER_MIN_SAMPLES <- 5L   # markers of a population only from samples with >= MIN_CELLS of it and of its comparison group
MARKER_TOP  <- 50L
MAP_PER_GROUP <- 2500L     # cells per level-2 group shipped for the maps (all cells of smaller groups)
set.seed(20261106)

suppressPackageStartupMessages({ library(Matrix) })
for (p in c("irlba", "uwot", "igraph", "harmony")) if (!requireNamespace(p, quietly = TRUE)) stop("package '", p, "' is required")
dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(LOG_DIR, paste0("07_cell_populations_", format(Sys.Date()), ".txt"))
logc <- function(...) { msg <- paste0(...); cat(msg, "\n"); cat(msg, "\n", file = log_file, append = TRUE) }
cat("", file = log_file)
logc("07_cell_populations.R started ", format(Sys.time(), tz = "UTC"), " UTC; compartment ", COMPARTMENT)

## ---- 0. taxonomy: every population must cite references that exist ---------------------------------------
pops <- read.csv(POP_FILE, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
refs <- read.csv(REF_FILE, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
pops <- pops[pops$compartment == COMPARTMENT, ]
split_ids <- function(x) { x <- unlist(strsplit(x, "[;,]")); x <- trimws(sub("\\s*\\(.*$", "", trimws(x))); x[nzchar(x)] }
cited <- unique(c(split_ids(pops$definition_refs), unlist(regmatches(pops$signature_ref, gregexpr("[a-z]+[0-9]{4}", pops$signature_ref)))))
missing_refs <- setdiff(cited, refs$ref_id)
if (length(missing_refs)) stop("populations.csv cites references missing from references.csv: ", paste(missing_refs, collapse = ", "))
no_ref <- pops$population_id[!nzchar(trimws(pops$definition_refs))]
if (length(no_ref)) stop("populations without a defining reference: ", paste(no_ref, collapse = ", "))
if (!nrow(pops)) stop("populations.csv has no populations for compartment '", COMPARTMENT, "'")
bad_m <- setdiff(pops$label_method, c("published", "signature", "score_only", "boundary"))
if (length(bad_m)) stop("populations.csv: unknown label_method ", paste(bad_m, collapse = ", "))
if (SIG_L1 && any(pops$level == 1 & pops$label_method == "published")) stop(COMPARTMENT, ": level 1 is assigned by signature; no level-1 population may be 'published'")
if (!SIG_L1 && any(pops$label_method == "boundary")) stop(COMPARTMENT, ": 'boundary' populations need level1 = \"signature\"")
if (length(VETO) && !all(c(names(VETO), VETO) %in% pops$population_id)) stop("cfg$veto names populations missing from populations.csv")
need_sig <- pops$level >= 1
if (any(!nzchar(pops$signature_genes[need_sig] ))) stop("every level-1/2 population needs signature_genes (published labels: for the marker check)")
sigs <- setNames(lapply(strsplit(pops$signature_genes, ";", fixed = TRUE), trimws), pops$population_id)
sigs <- sigs[lengths(sigs) > 0 & nzchar(vapply(sigs, `[`, "", 1))]
logc("taxonomy: ", nrow(pops), " populations, ", length(cited), " references cited; MD5 populations.csv ",
     unname(tools::md5sum(POP_FILE)), ", references.csv ", unname(tools::md5sum(REF_FILE)))
if (file.exists(DEC_FILE)) {                     # the decision record must name populations and references that exist (05 attaches it)
  decs <- read.csv(DEC_FILE, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8", colClasses = "character")
  decs <- decs[decs$compartment %in% c("all", COMPARTMENT), ]
  bad <- setdiff(split_ids(decs$applies_to), pops$population_id)
  if (length(bad)) stop("decisions.csv names populations missing from populations.csv: ", paste(bad, collapse = ", "))
  bad <- setdiff(split_ids(decs$refs), refs$ref_id)
  if (length(bad)) stop("decisions.csv cites references missing from references.csv: ", paste(bad, collapse = ", "))
  logc("decisions: ", nrow(decs), " for this compartment (incl. 'all'); MD5 decisions.csv ", unname(tools::md5sum(DEC_FILE)))
}
if (any(grepl("pending", pops$status))) logc("NOTE populations with markers pending verification: ",
                                             paste(pops$population_id[grepl("pending", pops$status)], collapse = ", "))

## ---- A. heavy step: compartment cells from the atlas (once) ---------------------------------------------
read_labels <- function(cf) {
  if (identical(cf$label_format, "tab")) {
    d <- read.delim(cf$label_file, stringsAsFactors = FALSE, check.names = FALSE)
    return(setNames(trimws(as.character(d[[cf$label_col]])), rownames(d)))
  }
  con <- gzfile(cf$label_file); on.exit(close(con))   # space-delimited; some values contain spaces, so the label is
  l <- readLines(con, encoding = "UTF-8")              # taken counting from the end of the line
  p <- strsplit(l[-1], " ", fixed = TRUE); k <- cf$label_from_end - 1L
  setNames(vapply(p, function(z) z[length(z) - k], ""), vapply(p, `[`, "", 1L))
}

# the atlas authors' own label of every cell (NA without one): kept for comparison in every compartment
atlas_labels <- function(cf, cells, atlas_type) {
  if (is.null(cf$label_file)) return(rep(NA_character_, length(cells)))
  lab <- read_labels(cf); v <- unname(lab[cells]); n_pub <- sum(!atlas_type %in% names(cf$type_map))
  logc(cf$name, ": atlas labels found for ", sum(!is.na(v)), " of ", length(cells), " cells; ", length(lab), " labels in file, MD5 ", unname(tools::md5sum(cf$label_file)))
  if (sum(!is.na(v)) < 0.8 * n_pub) stop("fewer than 80% of the cells carry an atlas label; check the barcode format of ", cf$label_file)
  v
}
# level-1 label of every cell of compartment `cf`: published label > coarse-type rule > unlabelled (signature mode: step B)
level1 <- function(cf, cells, atlas_type) {
  if (identical(cf$level1, "signature")) return(rep(NA_character_, length(cells)))
  l1 <- rep(NA_character_, length(cells))
  if (length(cf$type_map)) { w <- atlas_type %in% names(cf$type_map); l1[w] <- unname(cf$type_map[atlas_type[w]]) }
  if (!is.null(cf$label_file)) {
    lab <- read_labels(cf); w <- is.na(l1) & cells %in% names(lab)
    v <- unname(cf$label_map[lab[cells[w]]])
    if (anyNA(v)) stop("label values without a population id in ", cf$label_file, ": ", paste(setdiff(unique(lab[cells[w]]), names(cf$label_map)), collapse = ", "))
    l1[w] <- v
    n_pub <- sum(!atlas_type %in% names(cf$type_map))
    logc(cf$name, ": published labels matched for ", sum(w), " of ", n_pub, " cells (", round(100 * sum(w) / n_pub, 1), "%); ",
         length(lab), " labels in file, MD5 ", unname(tools::md5sum(cf$label_file)))
    if (sum(w) < 0.8 * n_pub) stop("fewer than 80% of the cells carry a published label; check the barcode format of ", cf$label_file)
  }
  l1[is.na(l1)] <- cf$unlabelled
  l1
}

if (RUN_STEP_A) {
  suppressPackageStartupMessages(library(SeuratObject))
  stopifnot(file.exists(ATLAS_RDS))
  logc("step A: reading atlas object (slow, memory-heavy) ...")
  obj <- readRDS(ATLAS_RDS)
  dat <- tryCatch(SeuratObject::LayerData(obj, assay = "RNA", layer = "data"), error = function(e) slot(obj[["RNA"]], "data"))
  meta <- obj@meta.data; rm(obj); invisible(gc())
  stopifnot(inherits(dat, "dgCMatrix"), identical(colnames(dat), rownames(meta)))
  keep <- meta$nFeature_RNA >= MIN_FEATURES
  dat <- dat[, keep]; meta <- meta[keep, ]; invisible(gc())
  nnz <- tabulate(dat@i + 1L, nbins = nrow(dat)); gkeep <- nnz >= MIN_GENE_CELLS
  logc("atlas: ", ncol(dat), " cells after QC; genes detected in >= ", MIN_GENE_CELLS, " cells: ", sum(gkeep))
  for (id in STEP_A_TODO) {
    cf <- COMPARTMENTS[[id]]
    inc <- as.character(meta$Clusters) %in% cf$atlas_types
    X <- dat[gkeep, inc]; m <- meta[inc, ]
    cm <- data.frame(cell = colnames(X), sample = as.character(m$Name), study = as.character(m$GSE.SRA..Study.),
                     disease_state = as.character(m$DiseaseState), met_site = as.character(m$If.metastatic..location),
                     treatment = as.character(m$Treatment), cell_type = as.character(m$Clusters), stringsAsFactors = FALSE)
    cm$l1 <- level1(cf, cm$cell, cm$cell_type); cm$atlas_label <- atlas_labels(cf, cm$cell, cm$cell_type)
    # same harmonisation of metadata values as 01_scatlas_summaries.R
    cm$met_site[!(cm$disease_state == "Metastatic lesion") | cm$met_site %in% c("2", "NA", "")] <- NA
    cm$treatment[cm$disease_state == "Donor"] <- "Donor (untreated)"
    cm$treatment[grepl("^Treatment na.{1,2}ve$", cm$treatment)] <- "Treatment-naive"
    cm$treatment[cm$treatment == "Treatment"] <- "Treated"
    saveRDS(list(X = X, meta = cm, note = "log1p(CP10k) of the compartment's cells; genes as in 01"), cache_file(id))
    logc("step A: ", id, " cached, ", ncol(X), " cells x ", nrow(X), " genes")
    rm(X, m, cm); invisible(gc())
  }
  rm(dat, meta); invisible(gc())
} else logc("step A skipped: using cache ", cache_file(COMPARTMENT_ID))
if (!COMPARTMENT_ID %in% STEP_A_TODO) logc("cache ", cache_file(COMPARTMENT_ID), " (step A not repeated for this compartment)")

## ---- B. light step ------------------------------------------------------------------------------------------
base <- readRDS(cache_file(COMPARTMENT_ID)); X <- base$X; cm <- base$meta; rm(base)
if (is.null(cm$atlas_label)) cm$atlas_label <- NA_character_      # caches written before atlas labels were kept
ROOT <- pops$population_id[pops$level == 0][1]
NGENES <- diff(X@p)                                                # detected genes per cell
NGENES_REL <- NGENES / ave(NGENES, cm$study, FUN = stats::median)  # ... relative to the median cell of the same study
if (!SIG_L1) logc("level-1 cells: ", paste(names(table(cm$l1)), as.integer(table(cm$l1)), sep = "=", collapse = ", "))

# module score of Tirosh et al. 2016: signature mean minus the mean of expression-matched control genes
module_scores <- function(X, sigs, nbin = 24L, nctrl = 100L, seed = 1L) {
  set.seed(seed); avg <- Matrix::rowMeans(X)
  bin <- cut(rank(avg, ties.method = "random") / length(avg), nbin, labels = FALSE)
  out <- matrix(NA_real_, ncol(X), length(sigs), dimnames = list(colnames(X), names(sigs)))
  used <- list()
  for (s in names(sigs)) {
    gi <- match(sigs[[s]], rownames(X)); gi <- gi[!is.na(gi)]; used[[s]] <- rownames(X)[gi]
    if (length(gi) < 3) next
    ctrl <- unique(unlist(lapply(gi, function(i) { pool <- setdiff(which(bin == bin[i]), gi); pool[sample.int(length(pool), min(nctrl, length(pool)))] })))
    out[, s] <- Matrix::colMeans(X[gi, , drop = FALSE]) - Matrix::colMeans(X[ctrl, , drop = FALSE])
  }
  attr(out, "genes_used") <- used; out
}

hvg <- function(X, n = N_HVG, min_frac = 0.01, exclude = "^(MT-|RPL|RPS|MRPL|MRPS|HB[ABDGQEMZ]|IG[HKL][VCJDG]|JCHAIN)") {
  nnz <- tabulate(X@i + 1L, nbins = nrow(X)); mu <- Matrix::rowMeans(X); v <- Matrix::rowMeans(X^2) - mu^2
  ok <- nnz >= min_frac * ncol(X) & !grepl(exclude, rownames(X)) & v > 0
  fit <- stats::loess(log10(v[ok]) ~ log10(mu[ok]), span = 0.3)
  res <- log10(v[ok]) - stats::predict(fit)
  rownames(X)[ok][order(res, decreasing = TRUE)][seq_len(min(n, sum(ok)))]
}

run_harmony <- function(pcs, meta, by) {
  f <- if (exists("RunHarmony", envir = asNamespace("harmony"))) harmony::RunHarmony else harmony::HarmonyMatrix
  e <- tryCatch(f(pcs, meta, by, verbose = FALSE), error = function(e1)
    tryCatch(f(data_mat = pcs, meta_data = meta, vars_use = by, do_pca = FALSE, verbose = FALSE), error = function(e2) stop(e2)))
  as.matrix(e)
}

leiden <- function(nn_idx, res) {
  k <- ncol(nn_idx); el <- cbind(rep(seq_len(nrow(nn_idx)), k - 1L), as.vector(nn_idx[, -1, drop = FALSE]))
  g <- igraph::simplify(igraph::graph_from_edgelist(el, directed = FALSE))
  arg <- if ("resolution" %in% names(formals(igraph::cluster_leiden))) "resolution" else "resolution_parameter"
  cl <- do.call(igraph::cluster_leiden, c(list(graph = g, objective_function = "modularity", n_iterations = 10L), setNames(list(res), arg)))
  igraph::membership(cl)
}

# map and clusters get their own seed per level (from the tag), so they do not depend on the signatures scored before
# them (module scores draw control genes at random); clusters are therefore the same whatever the taxonomy says
embed <- function(Xs, meta, tag) {          # HVG -> PCA -> Harmony -> UMAP + kNN
  set.seed(20261106L + sum(utf8ToInt(tag)))
  genes <- hvg(Xs)
  pc <- irlba::prcomp_irlba(Matrix::t(Xs[genes, ]), n = N_PC, center = TRUE, scale. = TRUE)$x
  rownames(pc) <- colnames(Xs)
  h <- if (length(unique(meta[[HARMONY_BY]])) > 1) run_harmony(pc, meta, HARMONY_BY) else pc
  u <- uwot::umap(h, n_neighbors = N_NEIGHBORS, min_dist = 0.3, ret_nn = TRUE, n_threads = max(1L, parallel::detectCores() - 1L))
  logc(tag, ": ", ncol(Xs), " cells, ", length(genes), " variable genes, ", N_PC, " PCs, Harmony by ", HARMONY_BY)
  list(umap = u$embedding, nn = u$nn$euclidean$idx, hvg = genes)
}

# name each cluster after the candidate signature with the highest mean score: S = scores of the candidate labels (cells of
# `idx`), P = score-only programmes (notes; vetoes), R = TRUE adds the 'resembles another level-1 population' note
call_clusters <- function(lin, cl, idx, S, P, R = TRUE) {
  ngenes <- NGENES[idx]
  l1_own <- if (lin %in% colnames(sc1)) lin else NA_character_
  do.call(rbind, lapply(sort(unique(cl)), function(k) {
    w <- cl == k; mu <- colMeans(S[w, , drop = FALSE]); o <- order(mu, decreasing = TRUE)
    best <- names(mu)[o[1]]; second <- if (length(o) > 1) mu[o[2]] else NA_real_
    med_genes <- stats::median(ngenes[w]); low_q <- stats::median(NGENES_REL[idx[w]]) < LOW_GENES_FRAC
    # why a cluster stays unresolved (first reason that applies); level-1 resemblance is reported, never relabelled
    resembles <- FALSE
    if (R) { l1m <- colMeans(sc1[idx[w], , drop = FALSE]); other <- setdiff(names(l1m), c(l1_own, UNLABELLED))
             oth <- other[which.max(l1m[other])]
             resembles <- !is.na(l1_own) && length(oth) && l1m[oth] >= MIN_SCORE && l1m[oth] - l1m[l1_own] >= MIN_MARGIN }
    pm <- if (ncol(P)) colMeans(P[w, , drop = FALSE]) else numeric(0)
    v_id <- unname(VETO[best]); vetoed <- length(v_id) == 1 && !is.na(v_id) && v_id %in% names(pm) && pm[[v_id]] > mu[o[1]]
    n_smp <- length(unique(cm$sample[idx][w]))
    reason <- if (sum(w) < MIN_CLUSTER) "small cluster" else if (n_smp < MIN_CLUSTER_SAMPLES) paste0("cells from only ", n_smp, " sample", if (n_smp > 1) "s") else if (low_q) "few genes detected per cell for its study (possible low-quality cells)" else
      if (mu[o[1]] < MIN_SCORE) "no signature reaches the minimum score" else
      if (!is.na(second) && mu[o[1]] - second < MIN_MARGIN) "two signatures score almost equally" else
      if (vetoed) paste0(pops$label[match(v_id, pops$population_id)], " scores higher than its signature (possible technical effect)") else ""
    ok <- !nzchar(reason)
    # notes shown with the cluster whether or not it is labelled: it looks like another lineage by the level-1 markers,
    # or a score-only programme (e.g. heat shock, metallothionein) is stronger than the label's own signature
    strong <- pm[pm >= MIN_SCORE & pm > mu[o[1]]]
    note <- c(if (resembles) paste0("resembles ", pops$label[match(oth, pops$population_id)], " by the level-1 markers"),
              if (length(strong)) paste0("strongest programme: ", pops$label[match(names(strong), pops$population_id)], sprintf(" (score %.2f)", strong)))
    note <- paste(note, collapse = "; ")
    data.frame(lineage = lin, cluster = paste0(lin, "_c", k), n_cells = sum(w), n_samples = length(unique(cm$sample[idx][w])),
               median_genes = med_genes, call = if (ok) best else paste0(lin, "_UNRESOLVED"), reason = reason, note = note, best = best,
               best_score = unname(mu[o[1]]), runner_up = if (length(o) > 1) names(mu)[o[2]] else NA_character_,
               margin = unname(mu[o[1]] - second), t(round(mu, 4)), check.names = FALSE, stringsAsFactors = FALSE)
  }))
}
log_calls <- function(tab) for (r in seq_len(nrow(tab))) logc(sprintf("  %-14s n=%6d samples=%3d genes=%5.0f best=%-13s %6.3f runner-up=%-13s margin=%6.3f -> %s%s%s",
  tab$cluster[r], tab$n_cells[r], tab$n_samples[r], tab$median_genes[r], tab$best[r], tab$best_score[r], tab$runner_up[r], tab$margin[r], tab$call[r],
  if (nzchar(tab$reason[r])) paste0(" [", tab$reason[r], "]") else "", if (nzchar(tab$note[r])) paste0(" {", tab$note[r], "}") else ""))

# ---- B0. compartment-wide score-only programmes (level 1, label_method score_only): scored on every cell ----
prog0 <- intersect(pops$population_id[pops$level == 1 & pops$label_method == "score_only"], names(sigs))
sc0 <- module_scores(X, sigs[prog0])
if (length(prog0)) logc("compartment-wide programmes: ", paste(prog0, collapse = ", "))

# ---- B1. level-1 map; published labels: marker check; signature mode: clusters named here ----
l1_ids <- pops$population_id[pops$level == 1 & pops$label_method != "score_only"]
sc1 <- module_scores(X, sigs[intersect(l1_ids, names(sigs))])
e1 <- embed(X, cm, "compartment map")
cm$umap1_x <- e1$umap[, 1]; cm$umap1_y <- e1$umap[, 2]
cluster_calls <- list(); sc2_all <- list(); cm$cluster1 <- NA_character_
if (SIG_L1) {
  cl1 <- leiden(e1$nn, LEIDEN_RES); cm$cluster1 <- paste0(ROOT, "_c", cl1)
  tab1 <- call_clusters(ROOT, cl1, seq_len(ncol(X)), sc1, sc0, R = FALSE)
  cm$l1 <- tab1$call[match(cm$cluster1, tab1$cluster)]; cluster_calls[[ROOT]] <- tab1
  logc(ROOT, ": ", nrow(tab1), " clusters -> ", paste(names(table(tab1$call)), as.integer(table(tab1$call)), sep = " x", collapse = ", ")); log_calls(tab1)
  logc("level-1 cells: ", paste(names(table(cm$l1)), as.integer(table(cm$l1)), sep = "=", collapse = ", "))
}
top1 <- colnames(sc1)[max.col(replace(sc1, is.na(sc1), -Inf), ties.method = "first")]
if (SIG_L1) {                                   # no published level-1 labels to check; the atlas labels are compared below
  check <- as.data.frame.matrix(table(population = cm$l1, top_marker_set = top1))
  check_mean <- aggregate(sc1, by = list(population = cm$l1), FUN = mean); agree <- stats::setNames(numeric(0), character(0))
} else {
  check <- as.data.frame.matrix(table(published = cm$l1, top_marker_set = top1))
  check_mean <- aggregate(sc1, by = list(published = cm$l1), FUN = mean)
  agree <- vapply(rownames(check), function(r) if (r %in% colnames(check)) check[r, r] / sum(check[r, ]) else NA_real_, numeric(1))
  logc("published level-1 labels vs their own marker set (share of cells scoring highest on it): ",
       paste(names(agree), round(100 * agree), sep = "=", collapse = "%, "), "%")
}

# ---- B2. level 2 within each lineage ----
cm$cluster <- NA_character_; cm$l2 <- cm$l1; cm$umap2_x <- NA_real_; cm$umap2_y <- NA_real_
l2_tab <- pops[pops$level == 2, ]
for (lin in unique(l2_tab$parent_id)) {
  kids <- l2_tab[l2_tab$parent_id == lin, ]; label_kids <- kids$population_id[kids$label_method == "signature"]
  idx <- which(cm$l1 == lin); if (length(idx) < 2 * MIN_CLUSTER) { logc(lin, ": too few cells (", length(idx), "), skipped"); next }
  Xs <- X[, idx]; sc2 <- module_scores(Xs, sigs[intersect(kids$population_id, names(sigs))]); sc2_all[[lin]] <- sc2
  for (s in colnames(sc2)) { miss <- setdiff(sigs[[s]], attr(sc2, "genes_used")[[s]]); if (length(miss)) logc("  ", s, ": signature genes not in the data: ", paste(miss, collapse = ", ")) }
  if (!length(label_kids)) { logc(lin, ": score-only programmes (", paste(colnames(sc2), collapse = ", "), "), no level-2 labels"); next }
  e2 <- embed(Xs, cm[idx, ], lin)
  cl <- leiden(e2$nn, LEIDEN_RES); cm$cluster[idx] <- paste0(lin, "_c", cl); cm$umap2_x[idx] <- e2$umap[, 1]; cm$umap2_y[idx] <- e2$umap[, 2]
  P <- cbind(sc2[, setdiff(colnames(sc2), label_kids), drop = FALSE], sc0[idx, , drop = FALSE])
  tab <- call_clusters(lin, cl, idx, sc2[, label_kids, drop = FALSE], P)
  cluster_calls[[lin]] <- tab
  cm$l2[idx] <- tab$call[match(cm$cluster[idx], tab$cluster)]
  logc(lin, ": ", nrow(tab), " clusters -> ", paste(names(table(tab$call)), as.integer(table(tab$call)), sep = " x", collapse = ", ")); log_calls(tab)
}
cluster_calls <- do.call(rbind, lapply(cluster_calls, function(d) d[, c("lineage", "cluster", "n_cells", "n_samples", "median_genes", "call", "reason", "note", "best", "best_score", "runner_up", "margin")]))
cluster_scores <- lapply(names(sc2_all), function(lin) { cl <- cm$cluster[cm$l1 == lin]
  if (all(is.na(cl))) return(NULL); aggregate(sc2_all[[lin]], by = list(cluster = cl), FUN = mean) })
names(cluster_scores) <- names(sc2_all)
if (SIG_L1) cluster_scores[[ROOT]] <- aggregate(cbind(sc1, sc0), by = list(cluster = cm$cluster1), FUN = mean)
# the atlas authors' labels of the cells of every population (both levels): the comparison shown on the cards
atlas_compare <- if (all(is.na(cm$atlas_label))) NULL else {
  ac <- rbind(data.frame(level = 1L, population_id = cm$l1, atlas_label = cm$atlas_label), data.frame(level = 2L, population_id = cm$l2, atlas_label = cm$atlas_label)[!is.na(cm$cluster), ])
  ac$atlas_label[is.na(ac$atlas_label)] <- "(no atlas label)"
  ac <- stats::aggregate(list(n_cells = rep(1L, nrow(ac))), ac[, c("level", "population_id", "atlas_label")], sum); ac[order(ac$level, ac$population_id, -ac$n_cells), ] }
if (!is.null(atlas_compare)) { a1 <- atlas_compare[atlas_compare$level == 1, ]
  for (pid in unique(a1$population_id)) { r <- a1[a1$population_id == pid, ]
    logc("  atlas labels of ", pid, ": ", paste0(r$atlas_label, " ", round(100 * r$n_cells / sum(r$n_cells)), "%", collapse = ", ")) } }
logc("level-2 cells: ", paste(names(table(cm$l2)), as.integer(table(cm$l2)), sep = "=", collapse = ", "))

# ---- B3. per (sample x population) summaries, as sc_expr: mean CP10k x100 and fraction expressing x10000 ----
summarise_groups <- function(level) {
  pop <- if (level == 1) cm$l1 else cm$l2
  if (level == 2) { has2 <- cm$l1 %in% l2_tab$parent_id[l2_tab$label_method == "signature"]; pop[!has2] <- NA }
  w <- which(!is.na(pop)); key <- factor(paste(cm$sample[w], pop[w], sep = "\r")); ng <- nlevels(key)
  ind <- sparseMatrix(i = seq_along(key), j = as.integer(key), x = 1, dims = c(length(key), ng))
  n <- as.numeric(Matrix::colSums(ind)); Xw <- X[, w]
  mean_cp10k <- sweep(as.matrix(expm1(Xw) %*% ind), 2, n, "/")
  Xb <- Xw; Xb@x[] <- 1; frac <- sweep(as.matrix(Xb %*% ind), 2, n, "/")
  parts <- do.call(rbind, strsplit(levels(key), "\r", fixed = TRUE))
  g <- data.frame(level = level, sample = parts[, 1], population_id = parts[, 2], n_cells = as.integer(n), stringsAsFactors = FALSE)
  list(groups = g, mean = mean_cp10k, frac = frac)
}
s1 <- summarise_groups(1); s2 <- summarise_groups(2)
groups <- rbind(s1$groups, s2$groups); groups$group_id <- seq_len(nrow(groups))
smeta <- cm[!duplicated(cm$sample), c("sample", "study", "disease_state", "met_site", "treatment")]
groups <- cbind(groups, smeta[match(groups$sample, smeta$sample), -1])
mean_all <- cbind(s1$mean, s2$mean); frac_all <- cbind(s1$frac, s2$frac); rownames(mean_all) <- rownames(frac_all) <- rownames(X)
expr <- list(genes = rownames(X),
             mean_cp10k_x100 = matrix(as.integer(round(mean_all * 100)), nrow(X), dimnames = list(rownames(X), NULL)),
             frac_x10000 = matrix(as.integer(round(frac_all * 10000)), nrow(X), dimnames = list(rownames(X), NULL)),
             note = "columns follow groups$group_id; divide by 100 / 10000; average expression for display = log1p(mean CP10k)")
logc("groups: ", sum(groups$level == 1), " sample x lineage, ", sum(groups$level == 2), " sample x sub-lineage; with >= ", MIN_CELLS,
     " cells: ", sum(groups$n_cells >= MIN_CELLS))

# ---- B4. markers: each population vs the rest of its parent, paired within samples ----
markers_for <- function(pid, parent_groups, level) {
  gi <- groups[groups$level == level & groups$population_id == pid, ]
  both <- intersect(gi$sample, parent_groups$sample)
  gi <- gi[match(both, gi$sample), ]; j <- match(both, parent_groups$sample)
  pg <- list(n_cells = parent_groups$n_cells[j], mean = parent_groups$mean[, j, drop = FALSE])
  n_rest <- pg$n_cells - gi$n_cells; ok <- gi$n_cells >= MIN_CELLS & n_rest >= MIN_CELLS
  if (sum(ok) < MARKER_MIN_SAMPLES) return(NULL)
  gi <- gi[ok, ]; pg <- list(n_cells = pg$n_cells[ok], mean = pg$mean[, ok, drop = FALSE]); n_rest <- n_rest[ok]
  mp <- mean_all[, gi$group_id, drop = FALSE]; mparent <- pg$mean
  mrest <- sweep(sweep(mparent, 2, pg$n_cells, "*") - sweep(mp, 2, gi$n_cells, "*"), 2, n_rest, "/"); mrest[mrest < 0] <- 0
  lfc <- log1p(mp) - log1p(mrest); fr <- frac_all[, gi$group_id, drop = FALSE]
  med <- apply(lfc, 1, stats::median); pos <- rowSums(lfc > 0); nn <- ncol(lfc); medfrac <- apply(fr, 1, stats::median)
  p <- stats::pbinom(pos - 1, nn, 0.5, lower.tail = FALSE)           # one-sided sign test across samples
  d <- data.frame(population_id = pid, gene = rownames(X), median_log_ratio = med, samples_up = pos, samples = nn,
                  median_frac = medfrac, p = p, stringsAsFactors = FALSE)
  d$q <- stats::p.adjust(d$p, "BH"); d <- d[d$median_frac >= 0.10 & d$q < 0.05 & d$median_log_ratio > 0, ]
  head(d[order(-d$median_log_ratio), ], MARKER_TOP)
}
parent_table <- function(level, parent) {    # per sample: cells and mean CP10k of the parent (all its cells)
  if (level == 1) { G <- groups[groups$level == 1, ] } else G <- groups[groups$level == 1 & groups$population_id == parent, ]
  smp <- unique(G$sample)
  n <- vapply(smp, function(s) sum(G$n_cells[G$sample == s]), numeric(1))
  m <- vapply(smp, function(s) { r <- G[G$sample == s, ]; drop(mean_all[, r$group_id, drop = FALSE] %*% r$n_cells) / sum(r$n_cells) }, numeric(nrow(X)))
  list(sample = smp, n_cells = n, mean = m)
}
markers <- list()
P1 <- parent_table(1, NULL)
for (pid in setdiff(unique(groups$population_id[groups$level == 1]), UNLABELLED)) markers[[pid]] <- markers_for(pid, P1, 1)
for (par in unique(cluster_calls$lineage)) { Pp <- parent_table(2, par)
  for (pid in unique(groups$population_id[groups$level == 2 & groups$population_id %in% cm$l2[cm$l1 == par]])) markers[[pid]] <- markers_for(pid, Pp, 2) }
markers <- do.call(rbind, markers); rownames(markers) <- NULL
logc("markers: ", if (is.null(markers)) 0 else nrow(markers), " population-gene pairs (q < 0.05, detected in >= 10% of cells)")

# ---- B5. map subsample: coordinates, labels, signature scores and expression for gene overlays ----
grp <- ifelse(is.na(cm$cluster), cm$l1, cm$l2)
pick <- sort(unlist(lapply(split(seq_len(nrow(cm)), grp), function(i) if (length(i) <= MAP_PER_GROUP) i else sample(i, MAP_PER_GROUP))))
scores <- matrix(NA_real_, nrow(cm), 0); rownames(scores) <- cm$cell
scores <- cbind(scores, round(sc1, 3), round(sc0, 3))
for (lin in names(sc2_all)) { s <- sc2_all[[lin]]; full <- matrix(NA_real_, nrow(cm), ncol(s), dimnames = list(cm$cell, colnames(s)))
  full[rownames(s), ] <- round(s, 3); scores <- cbind(scores, full) }
Xm <- Matrix::t(X[, pick]); Xm@x <- round(Xm@x, 2); Xm <- as(Xm, "CsparseMatrix")
map <- list(meta = cm[pick, c("sample", "study", "disease_state", "met_site", "treatment", "l1", "l2", "cluster", "umap1_x", "umap1_y", "umap2_x", "umap2_y")],
            scores = scores[pick, , drop = FALSE], expr = Xm,
            note = paste0("up to ", MAP_PER_GROUP, " cells per population; expression log1p(CP10k), 2 decimals; scores: module scores (Tirosh et al. 2016)"))
rownames(map$meta) <- NULL; rownames(map$scores) <- NULL
logc("map subsample: ", length(pick), " cells")

## ---- save ----------------------------------------------------------------------------------------------------
out <- list(compartment = COMPARTMENT, level1_method = if (SIG_L1) "signature" else "published", atlas_compare = atlas_compare,
            taxonomy = pops, references = refs, groups = groups, expr = expr,
            label_check = list(counts = check, mean_scores = check_mean, agreement = agree),
            clusters = list(calls = cluster_calls, scores = cluster_scores), markers = markers, map = map,
            params = list(MIN_FEATURES = MIN_FEATURES, N_HVG = N_HVG, N_PC = N_PC, HARMONY_BY = HARMONY_BY, N_NEIGHBORS = N_NEIGHBORS,
                          LEIDEN_RES = LEIDEN_RES, MIN_CLUSTER = MIN_CLUSTER, MIN_CLUSTER_SAMPLES = MIN_CLUSTER_SAMPLES, MIN_SCORE = MIN_SCORE, MIN_MARGIN = MIN_MARGIN, LOW_GENES_FRAC = LOW_GENES_FRAC,
                          MIN_CELLS = MIN_CELLS, MARKER_MIN_SAMPLES = MARKER_MIN_SAMPLES, VETO = VETO),
            built = format(Sys.time(), tz = "UTC", usetz = TRUE))
lab_out <- file.path(WORK_DIR, paste0("pop_", COMPARTMENT_ID, "_cell_labels.tsv.gz"))   # every cell's labels, for reuse and checking
write.table(cm[, c("cell", "sample", "atlas_label", "cluster1", "l1", "cluster", "l2")], gzfile(lab_out), sep = "\t", quote = FALSE, row.names = FALSE)
logc("cell labels ", lab_out, " MD5 ", unname(tools::md5sum(lab_out)))
f <- file.path(OUT_DIR, paste0("pop_", COMPARTMENT_ID, ".rds")); saveRDS(out, f, compress = "xz")
logc("output ", basename(f), " ", round(file.size(f) / 1e6, 1), " MB, MD5 ", unname(tools::md5sum(f)))
logc("finished ", format(Sys.time(), tz = "UTC"), " UTC")
capture.output(sessionInfo(), file = log_file, append = TRUE)
