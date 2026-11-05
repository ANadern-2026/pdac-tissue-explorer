# =====================================================================================================
# 06_models.R
# Model systems of PDAC next to patient tissue: cell lines and organoids from DepMap (Broad) and organoids from the
# GDC (ORGANOID-PANCREATIC, Tiriac et al. 2018; HCMI-CMDC, pancreas).
#  - RNA: log2(TPM + 1) per model (compared only within one technology: DepMap lines with DepMap lines, including all
#    other cancer lineages; GDC organoids with GDC organoids and their source tumours), and each gene's within-sample
#    percentile; for patient tissue (single-cell atlas, cell-weighted mean per treatment-naive primary tumour and cell type)
#    the percentile is shown as background only, never on one axis with the models (different technology and composition).
#  - CRISPR knockout effect (DepMap Chronos gene effect) in all screened lines (PDAC and every other lineage).
#  - Protein: DepMap harmonised mass-spectrometry proteomics of the CCLE (Nusinow et al. 2020; Broad data only).
#  - Per model: classical vs basal-like score relative to the other models of the same type and collection (Moffitt genes,
#    mean z-score of basal-like minus classical genes), similarity to atlas tumour cells of treatment-naive primary tumours
#    and of liver metastases (Spearman over the genes that best separate cell types), and for HCMI organoids with a sequenced
#    source tumour: correlation with their own tumour vs the other source tumours (same GDC pipeline); per gene, the
#    organoid and source-tumour values of these pairs.
# Inputs : data-raw/source/models/ (06a_fetch_models.R), app/data/genes.rds, gene_xrefs.rds, sc_groups.rds, sc_expr.rds,
#          bulk_tcga_gtex.rds, clinical.rds (Moffitt gene lists)
# Output : app/data/models.rds (then 05_build_store.R builds models.sqlite)
# Run from the repository root:  source("data-raw/06_models.R")   ~5-10 min, ~8 GB RAM.
# =====================================================================================================
source(file.path("data-raw", "helpers.R")); log_file <- start_log("06_models")
for (f in c("00_config.R", "01_data.R", "02_stats.R")) source(file.path("app", "R", f))
IN_DIR <- Sys.getenv("PDACX_DATA_DIR", file.path("app", "data"))
SRC <- file.path("data-raw", "source", "models"); D_DM <- file.path(SRC, "depmap"); D_GDC <- file.path(SRC, "gdc")

MODEL_CODES <- c(PAAD = "Pancreatic adenocarcinoma", PAASC = "Adenosquamous carcinoma of the pancreas")   # DepMap OncotreeCode
HCMI_PDAC_DX <- c("Infiltrating duct carcinoma, NOS", "Adenocarcinoma, NOS", "Adenocarcinoma, metastatic, NOS")
FID_GENES <- 2000L          # genes used for the similarity scores (most variable between cell types / between source tumours)

src_files <- c(file.path(D_DM, c("Model.csv", "CRISPRGeneEffect.csv", "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")),
               file.path(D_DM, "proteomics", c("harmonized_MS_CCLE_Gygi.csv", "uniprot_hugo_entrez_id_mapping.csv")),
               file.path(D_GDC, "gdc_expression_files.tsv"), file.path(IN_DIR, c("genes.rds", "gene_xrefs.rds", "sc_groups.rds", "sc_expr.rds", "bulk_tcga_gtex.rds", "clinical.rds")))
miss <- src_files[!file.exists(src_files)]; if (length(miss)) stop("missing inputs (run data-raw/06a_fetch_models.R first?): ", paste(miss, collapse = ", "))
src_md5 <- stats::setNames(unname(tools::md5sum(src_files)), sub("^.*/(source/models/)?", "", src_files))
for (k in names(src_md5)) logc(log_file, "source ", k, " MD5 ", src_md5[[k]])
release <- readLines(file.path(D_DM, "RELEASE.txt")); logc(log_file, "DepMap releases: ", paste(release, collapse = "; "))

read_wide <- function(f) {     # DepMap wide CSV -> data.frame (data.table::fread when available: much faster)
  if (requireNamespace("data.table", quietly = TRUE)) as.data.frame(data.table::fread(f, showProgress = FALSE), stringsAsFactors = FALSE)
  else utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
}
pct_cols <- function(X) {      # within-column percentile x10 (0..1000), ties -> lowest rank, so undetected genes are 0
  n <- nrow(X); out <- matrix(NA_integer_, n, ncol(X), dimnames = dimnames(X))
  for (j in seq_len(ncol(X))) { v <- X[, j]; ok <- !is.na(v); m <- sum(ok)
    if (m > 1) out[ok, j] <- as.integer(round(1000 * (rank(v[ok], ties.method = "min") - 1) / (m - 1))) }
  out
}

## ---- gene index --------------------------------------------------------------------------------------------
genes <- readRDS(file.path(IN_DIR, "genes.rds"))$genes; xr <- readRDS(file.path(IN_DIR, "gene_xrefs.rds"))$xrefs
stopifnot(identical(xr$symbol, genes$symbol))
idx_entrez <- function(eg) match(as.character(eg), genes$entrez)
idx_ensembl <- function(ens) match(sub("\\..*$", "", ens), xr$ensembl)

## ---- DepMap models ---------------------------------------------------------------------------------------------
mdl <- utils::read.csv(file.path(D_DM, "Model.csv"), stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA"))
dm_pdac <- mdl[mdl$OncotreeCode %in% names(MODEL_CODES), ]
logc(log_file, "DepMap PDAC models (", paste(names(MODEL_CODES), collapse = ", "), "): ", nrow(dm_pdac), " (", sum(dm_pdac$ModelType == "Cell Line", na.rm = TRUE),
     " cell lines, ", sum(dm_pdac$ModelType == "Organoid", na.rm = TRUE), " organoids)")

ex <- read_wide(file.path(D_DM, "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"))
meta_cols <- intersect(c("", "V1", "SequencingID", "ModelConditionID", "ModelID", "IsDefaultEntryForMC", "IsDefaultEntryForModel"), names(ex))
ex <- ex[ex$IsDefaultEntryForModel %in% c("Yes", "TRUE", TRUE) & ex$ModelID %in% mdl$ModelID, , drop = FALSE]
stopifnot(!anyDuplicated(ex$ModelID))
gcols <- setdiff(names(ex), meta_cols); eg <- sub("^.*\\((\\d+)\\)$", "\\1", gcols); gi <- idx_entrez(eg)
keep <- !is.na(gi) & !duplicated(gi)
R_all <- t(as.matrix(ex[, gcols[keep]])); dimnames(R_all) <- list(genes$symbol[gi[keep]], ex$ModelID); storage.mode(R_all) <- "double"; rm(ex)
rlines <- data.frame(model = colnames(R_all), name = mdl$CellLineName[match(colnames(R_all), mdl$ModelID)], lineage = mdl$OncotreeLineage[match(colnames(R_all), mdl$ModelID)],
                     type = mdl$ModelType[match(colnames(R_all), mdl$ModelID)], pdac = colnames(R_all) %in% dm_pdac$ModelID, stringsAsFactors = FALSE)
R_dm <- R_all[, colnames(R_all) %in% dm_pdac$ModelID, drop = FALSE]
logc(log_file, "DepMap RNA: ", ncol(R_all), " models of all lineages (", ncol(R_dm), " PDAC) x ", length(gcols), " genes; ", sum(keep), " mapped to the gene index (Entrez)")

## ---- GDC organoids ----------------------------------------------------------------------------------------------
gm <- utils::read.delim(file.path(D_GDC, "gdc_expression_files.tsv"), stringsAsFactors = FALSE, na.strings = c("", "NA"))
gm$path <- file.path(D_GDC, gm$set, gm$file_name); stopifnot(all(file.exists(gm$path)))
is_model <- grepl("Next Generation Cancer Model", gm$sample_type)
# ORGANOID-PANCREATIC: tumour and normal organoids; one per case and tissue type
org <- gm[gm$set == "ORGANOID" & is_model, ]; org <- org[order(org$case, org$tissue_type, org$sample), ]
org <- org[!duplicated(paste(org$case, org$tissue_type)), ]
# HCMI: PDAC organoids (one per case, 'Next Generation' before 'Expanded') and their source tumours (solid tissue)
hc <- gm[gm$set == "HCMI", ]
hc_dx <- tapply(hc$diagnosis, hc$case, function(x) x[!is.na(x)][1])
hm <- hc[grepl("Next Generation Cancer Model", hc$sample_type) & hc$tissue_type %in% "Tumor" & hc_dx[hc$case] %in% HCMI_PDAC_DX, ]
hm <- hm[order(hm$case, grepl("^Expanded", hm$sample_type), hm$sample), ]; hm <- hm[!duplicated(hm$case), ]
ht <- hc[hc$sample_type %in% c("Primary Tumor", "Metastatic") & hc$case %in% hm$case & !grepl("FFPE", hc$sample_type), ]
ht <- ht[order(ht$case, ht$sample), ]; ht <- ht[!duplicated(ht$case), ]
logc(log_file, "GDC: ORGANOID-PANCREATIC ", sum(org$tissue_type == "Tumor"), " tumour + ", sum(org$tissue_type == "Normal"), " normal organoids; HCMI PDAC organoids ",
     nrow(hm), " (", nrow(ht), " with a sequenced source tumour); excluded HCMI pancreas models with other diagnoses: ",
     paste(names(table(hc_dx[setdiff(unique(hc$case[grepl("Model", hc$sample_type)]), hm$case)])), collapse = "; "))
read_star <- function(f) {
  x <- utils::read.delim(f, comment.char = "#", stringsAsFactors = FALSE)
  x <- x[grepl("^ENSG", x$gene_id) & x$gene_type == "protein_coding", ]
  structure(log2(x$tpm_unstranded + 1), names = x$gene_id, symbol = x$gene_name)
}
gdc_mat <- function(m) {
  v <- lapply(m$path, read_star); ids <- names(v[[1]]); stopifnot(all(vapply(v, function(z) identical(names(z), ids), logical(1))))
  X <- do.call(cbind, v)
  gi <- idx_ensembl(ids)                                          # Ensembl ID first (unique Entrez-Ensembl pairs only, 04b) ...
  gs <- match(attr(v[[1]], "symbol"), genes$symbol); fill <- is.na(gi) & !is.na(gs) & !gs %in% gi; gi[fill] <- gs[fill]   # ... then the symbol
  keep <- !is.na(gi) & !duplicated(gi)
  X <- X[keep, , drop = FALSE]; dimnames(X) <- list(genes$symbol[gi[keep]], m$sample); X
}
gall <- rbind(org, hm, ht); R_gdc <- gdc_mat(gall)
logc(log_file, "GDC RNA: ", ncol(R_gdc), " samples x ", nrow(R_gdc), " protein-coding genes mapped to the gene index (Ensembl)")

## ---- model table -----------------------------------------------------------------------------------------------
num <- function(x) suppressWarnings(as.numeric(x))
site <- function(x) ifelse(grepl("^Primary", x, ignore.case = TRUE), "Primary tumour", ifelse(grepl("Metasta", x, ignore.case = TRUE), "Metastasis", NA))
dm_tab <- data.frame(model = dm_pdac$ModelID, name = dm_pdac$CellLineName, source = "DepMap",
                     type = ifelse(dm_pdac$ModelType %in% "Organoid", "Organoid", "Cell line"), tumour = ifelse(dm_pdac$OncotreeCode == "PAASC", "Adenosquamous", "Adenocarcinoma"),
                     origin = site(dm_pdac$PrimaryOrMetastasis), collection_site = dm_pdac$SampleCollectionSite, sex = dm_pdac$Sex, age = num(dm_pdac$Age),
                     depmap_id = dm_pdac$ModelID, hcmi_id = dm_pdac$HCMIID, rrid = dm_pdac$RRID, supplier = dm_pdac$SourceType, catalog = dm_pdac$CatalogNumber, stringsAsFactors = FALSE)
# (if, not ifelse(), for values that depend on the single `type`: ifelse() would recycle its first element over all rows)
gdc_tab <- function(m, src, type) data.frame(model = m$sample, name = if (src == "HCMI") m$case else paste0("Organoid case ", m$case, if (type == "Normal organoid") " (normal)" else ""),
                                             source = src, type = type, tumour = if (type == "Normal organoid") NA_character_ else "Adenocarcinoma",
                                             origin = if (type == "Normal organoid") NA_character_ else site(m$tumor_descriptor), collection_site = NA_character_, sex = NA_character_, age = NA_real_,
                                             depmap_id = NA_character_, hcmi_id = if (src == "HCMI") m$case else NA_character_, rrid = NA_character_, supplier = NA_character_, catalog = NA_character_,
                                             stringsAsFactors = FALSE)
models <- rbind(dm_tab, gdc_tab(org[org$tissue_type == "Tumor", ], "Tiriac 2018", "Organoid"), gdc_tab(org[org$tissue_type == "Normal", ], "Tiriac 2018", "Normal organoid"),
                gdc_tab(hm, "HCMI", "Organoid"))
# HCMI models also profiled by DepMap are shown once (DepMap profile: it also has CRISPR data); their GDC profile is kept for the pairing
dup <- models$source == "HCMI" & models$hcmi_id %in% stats::na.omit(dm_tab$hcmi_id)
logc(log_file, "HCMI models also in DepMap (shown once, as the DepMap model): ", sum(dup))
models <- models[!dup, ]
models$has_rna <- models$model %in% c(colnames(R_dm), colnames(R_gdc))
models <- models[models$has_rna, ]

## ---- universe and percentiles --------------------------------------------------------------------------------------
sc_groups <- prepare_sc_groups(readRDS(file.path(IN_DIR, "sc_groups.rds"))); e <- readRDS(file.path(IN_DIR, "sc_expr.rds"))
bulk <- readRDS(file.path(IN_DIR, "bulk_tcga_gtex.rds"))
pc <- genes$symbol[xr$gene_type %in% "protein-coding"]
in_sc <- genes$symbol[!is.na(genes$sym_sc) & genes$sym_sc %in% e$genes]; in_tcga <- genes$symbol[!is.na(genes$sym_tcga_gtex) & genes$sym_tcga_gtex %in% rownames(bulk$log2tpm_x1000)]
U <- Reduce(intersect, list(pc, rownames(R_dm), rownames(R_gdc), in_sc)); U <- genes$symbol[genes$symbol %in% U]
U_tcga <- intersect(U, in_tcga)       # TCGA (recount3, GENCODE v26) lacks some genes: its percentiles are over the universe genes it has
logc(log_file, "universe: ", length(U), " protein-coding genes measured in DepMap, GDC and the atlas (", length(U_tcga), " of them also in TCGA)")
# models: one matrix over the union of mapped genes (NA where a source lacks the gene)
allg <- genes$symbol[genes$symbol %in% union(rownames(R_dm), rownames(R_gdc))]
M <- matrix(NA_real_, length(allg), nrow(models), dimnames = list(allg, models$model))
j <- models$model[models$source == "DepMap"]; M[rownames(R_dm), j] <- R_dm[, j]
j <- models$model[models$source != "DepMap"]; M[rownames(R_gdc), j] <- R_gdc[, j]
P_models <- pct_cols(M[U, , drop = FALSE])
# tissue units: atlas (treatment-naive primary tumours, cell-weighted mean per sample and cell type, >= MIN_CELLS cells) + TCGA PDAC
units <- list(); Ts <- list()
scg <- genes$sym_sc[match(U, genes$symbol)]; scr <- match(scg, e$genes)
for (cc in CELLTYPE_ORDER) {
  u <- coexp_units(sc_groups, cc, state = "Primary tumour", min_cells = MIN_CELLS, min_study = 1L)
  if (!length(u$samples)) next
  X <- (e$mean_cp10k_x100[scr, u$rows, drop = FALSE] / 100) %*% u$W; rownames(X) <- U
  Ts[[cc]] <- X; units[[cc]] <- data.frame(unit = paste(cc, u$samples, sep = "|"), group = cc, sample = u$samples, study = u$study, n_cells = u$n_cells, stringsAsFactors = FALSE)
}
T_sc <- do.call(cbind, Ts); u_sc <- do.call(rbind, units); dimnames(T_sc) <- list(U, u_sc$unit)
tc <- bulk$samples$sample[bulk$samples$cohort == "PDAC (TCGA)"]
T_tcga <- matrix(NA_real_, length(U), length(tc), dimnames = list(U, paste0("TCGA|", tc)))
T_tcga[U_tcga, ] <- bulk$log2tpm_x1000[genes$sym_tcga_gtex[match(U_tcga, genes$symbol)], tc, drop = FALSE] / 1000
tissue_units <- rbind(u_sc, data.frame(unit = colnames(T_tcga), group = "PDAC bulk tumours (TCGA)", sample = tc, study = "TCGA-PAAD", n_cells = NA_integer_))
P_tissue <- cbind(pct_cols(T_sc), pct_cols(T_tcga))
logc(log_file, "tissue units: ", paste(names(table(tissue_units$group)), table(tissue_units$group), sep = " = ", collapse = "; "))

## ---- classical / basal-like score --------------------------------------------------------------------------------------
moff <- readRDS(file.path(IN_DIR, "clinical.rds"))$moffitt
mb <- intersect(moff$basal, U); mc <- intersect(moff$classical, U)
logc(log_file, "Moffitt score from ", length(mb), " basal-like and ", length(mc), " classical genes in the universe (missing: ",
     paste(setdiff(c(moff$basal, moff$classical), U), collapse = ", "), ")")
moffitt_pct <- function(P) round((colMeans(P[mb, , drop = FALSE], na.rm = TRUE) - colMeans(P[mc, , drop = FALSE], na.rm = TRUE)) / 10, 1)
models$moffitt <- unname(moffitt_pct(P_models)[models$model])
tissue_units$moffitt <- unname(moffitt_pct(P_tissue)[tissue_units$unit])
# relative subtype: within each model type and collection (same technology), mean z-score of basal-like minus classical genes
models$family <- ifelse(models$type == "Normal organoid", NA, paste(models$source, models$type))
models$subtype_rel <- NA_real_
for (fam in stats::na.omit(unique(models$family))) {
  j <- models$model[models$family %in% fam]; if (length(j) < 5) next
  X <- M[intersect(c(mb, mc), rownames(M)), j, drop = FALSE]; X <- X[apply(X, 1, function(v) isTRUE(stats::sd(v, na.rm = TRUE) > 0)), , drop = FALSE]
  Z <- t(scale(t(X))); models$subtype_rel[match(j, models$model)] <- round(colMeans(Z[intersect(mb, rownames(Z)), , drop = FALSE]) - colMeans(Z[intersect(mc, rownames(Z)), , drop = FALSE]), 2)
}
logc(log_file, "relative subtype score per family: ", paste(names(table(models$family[!is.na(models$subtype_rel)])), table(models$family[!is.na(models$subtype_rel)]), sep = " = ", collapse = "; "))

## ---- similarity to atlas cell types ---------------------------------------------------------------------------------------
prof <- sapply(names(Ts), function(cc) rowMeans(log1p(Ts[[cc]])))                    # mean log1p(CP10k) per cell type
fid_genes <- names(sort(apply(prof, 1, stats::var), decreasing = TRUE))[seq_len(min(FID_GENES, nrow(prof)))]
Mf <- M[fid_genes, models$model, drop = FALSE]
sim <- stats::cor(Mf, prof[fid_genes, ], method = "spearman", use = "pairwise.complete.obs")
models$sim_tumour <- round(unname(sim[models$model, "Ductal: tumour subclusters"]), 3)
ulm <- coexp_units(sc_groups, "Ductal: tumour subclusters", state = "Liver metastasis", min_cells = MIN_CELLS, min_study = 1L)
prof_lm <- rowMeans(log1p((e$mean_cp10k_x100[scr, ulm$rows, drop = FALSE] / 100) %*% ulm$W)); names(prof_lm) <- U
models$sim_livermet <- round(unname(stats::cor(Mf, prof_lm[fid_genes], method = "spearman", use = "pairwise.complete.obs")[models$model, 1]), 3)
logc(log_file, "liver-metastasis tumour-cell reference: ", length(ulm$samples), " samples from ", ulm$n_studies, " studies")
models$closest <- colnames(sim)[apply(sim, 1, which.max)][match(models$model, rownames(sim))]
logc(log_file, "similarity to atlas cell types over ", length(fid_genes), " genes; closest cell type: ",
     paste(names(table(models$closest)), table(models$closest), sep = " = ", collapse = "; "))
# HCMI: own source tumour vs other source tumours (GDC profiles of both)
pairs <- data.frame(model = hm$sample, case = hm$case, tumour = ht$sample[match(hm$case, ht$case)], stringsAsFactors = FALSE); pairs <- pairs[!is.na(pairs$tumour), ]
if (nrow(pairs) >= 5) {
  # each gene centred within the organoids and within the tumours (removes the organoid-vs-tissue offset, e.g. stroma), genes most
  # variable between the organoids: the correlation then reflects what is specific to a patient (tested 2026-11-04: own tumour
  # closest for 31 of 42 this way vs 3 of 42 with uncentred values over genes variable between tumours, which mostly track stroma)
  Mo <- R_gdc[U, pairs$model, drop = FALSE]; Tt <- R_gdc[U, pairs$tumour, drop = FALSE]
  vg <- names(sort(apply(Mo, 1, stats::var), decreasing = TRUE))[seq_len(min(FID_GENES, length(U)))]
  C <- stats::cor(Mo[vg, ] - rowMeans(Mo[vg, ]), Tt[vg, ] - rowMeans(Tt[vg, ]), method = "spearman")
  pairs$r_own <- round(diag(C), 3); pairs$r_other <- round(vapply(seq_len(nrow(pairs)), function(i) stats::median(C[i, -i]), numeric(1)), 3)
  pairs$rank_own <- vapply(seq_len(nrow(pairs)), function(i) sum(C[i, ] >= C[i, i]), integer(1))     # 1 = own tumour is the closest
  logc(log_file, "HCMI pairs: ", nrow(pairs), "; own tumour closest for ", sum(pairs$rank_own == 1), "; median r own ", stats::median(pairs$r_own), " vs other ", stats::median(pairs$r_other))
  shown <- ifelse(pairs$case %in% dm_tab$hcmi_id, dm_tab$model[match(pairs$case, dm_tab$hcmi_id)], pairs$model)
  models$r_own_tumour <- pairs$r_own[match(models$model, shown)]; models$r_other_tumours <- pairs$r_other[match(models$model, shown)]
  models$own_tumour_rank <- pairs$rank_own[match(models$model, shown)]; models$n_source_tumours <- ifelse(is.na(models$own_tumour_rank), NA_integer_, nrow(pairs))
  pair_organoid <- R_gdc[, pairs$model, drop = FALSE]; pair_tumour <- R_gdc[, pairs$tumour, drop = FALSE]
}

## ---- CRISPR ------------------------------------------------------------------------------------------------------------------
cr <- read_wide(file.path(D_DM, "CRISPRGeneEffect.csv")); idc <- names(cr)[1]
gcols <- names(cr)[-1]; gi <- idx_entrez(sub("^.*\\((\\d+)\\)$", "\\1", gcols)); keep <- !is.na(gi) & !duplicated(gi)
CR <- t(as.matrix(cr[, gcols[keep]])); dimnames(CR) <- list(genes$symbol[gi[keep]], cr[[idc]])
lines <- data.frame(model = colnames(CR), name = mdl$CellLineName[match(colnames(CR), mdl$ModelID)], lineage = mdl$OncotreeLineage[match(colnames(CR), mdl$ModelID)],
                    code = mdl$OncotreeCode[match(colnames(CR), mdl$ModelID)], type = mdl$ModelType[match(colnames(CR), mdl$ModelID)], stringsAsFactors = FALSE)
lines$pdac <- lines$code %in% names(MODEL_CODES)
logc(log_file, "CRISPR: ", ncol(CR), " lines (", sum(lines$pdac), " PDAC) x ", nrow(CR), " genes mapped (of ", length(gcols), ")")
models$has_crispr <- models$model %in% colnames(CR)

## ---- proteomics (CCLE, Gygi) ------------------------------------------------------------------------------------------------
pr <- read_wide(file.path(D_DM, "proteomics", "harmonized_MS_CCLE_Gygi.csv")); idp <- names(pr)[1]
um <- utils::read.csv(file.path(D_DM, "proteomics", "uniprot_hugo_entrez_id_mapping.csv"), stringsAsFactors = FALSE)
pcols <- names(pr)[-1]; acc <- sub("-\\d+$", "", pcols); canon <- !grepl("-\\d+$", pcols)
eg <- um$EntrezID[match(acc, um$UniprotID)]; gi <- idx_entrez(ifelse(is.na(eg), NA, sub("\\.0$", "", format(eg, scientific = FALSE, trim = TRUE))))
PX <- as.matrix(pr[, pcols]); n_ok <- colSums(!is.na(PX))
o <- order(is.na(gi), !canon, -n_ok); pick <- o[!is.na(gi[o]) & !duplicated(gi[o])]
PR <- t(PX[, pick, drop = FALSE]); dimnames(PR) <- list(genes$symbol[gi[pick]], pr[[idp]])
plines <- data.frame(model = colnames(PR), name = mdl$CellLineName[match(colnames(PR), mdl$ModelID)], lineage = mdl$OncotreeLineage[match(colnames(PR), mdl$ModelID)], stringsAsFactors = FALSE)
plines$pdac <- mdl$OncotreeCode[match(colnames(PR), mdl$ModelID)] %in% names(MODEL_CODES)
logc(log_file, "proteomics: ", ncol(PR), " lines (", sum(plines$pdac), " PDAC) x ", nrow(PR), " genes (one protein entry per gene: canonical isoform, most lines quantified)")
models$has_protein <- models$model %in% colnames(PR)

## ---- save ------------------------------------------------------------------------------------------------------------------------
i1000 <- function(X) { Y <- round(X * 1000); storage.mode(Y) <- "integer"; Y }
out <- list(models = models, rna_log2tpm_x1000 = i1000(M), rna_lines_x1000 = i1000(R_all), rna_lines = rlines,
            pair_organoid_x1000 = i1000(pair_organoid), pair_tumour_x1000 = i1000(pair_tumour), universe = U, pct_models_x10 = P_models, tissue_units = tissue_units, pct_tissue_x10 = P_tissue,
            similarity = round(sim, 3), similarity_genes = fid_genes, hcmi_pairs = pairs,
            crispr_x1000 = i1000(CR), crispr_lines = lines, protein_x1000 = i1000(PR), protein_lines = plines,
            moffitt_genes = list(basal = mb, classical = mc),
            source = list(depmap_release = release, source_md5 = src_md5, built = format(Sys.time(), tz = "UTC"),
                          citation = c("DepMap, Broad (2026): DepMap Public 26Q3. https://depmap.org (CC BY 4.0).",
                                       "Nusinow D.P. et al. Quantitative proteomics of the Cancer Cell Line Encyclopedia. Cell 2020;180:387-402.",
                                       "Tiriac H. et al. Organoid profiling identifies common responders to chemotherapy in pancreatic cancer. Cancer Discov 2018;8:1112-1129 (GDC ORGANOID-PANCREATIC).",
                                       "Human Cancer Models Initiative (HCMI), https://ocg.cancer.gov/programs/HCMI; dbGaP phs001486 (GDC HCMI-CMDC)."),
                          references = c(primary = paste0("atlas tumour cells (ductal tumour subclusters), treatment-naive primary tumours, ", sum(tissue_units$group == "Ductal: tumour subclusters"), " samples"),
                                         liver_met = paste0("atlas tumour cells, liver metastases, ", length(ulm$samples), " samples")),
                          units = c(rna = "log2(TPM + 1)", pct = "within-sample percentile among the universe genes x10 (0 = lowest / not detected)",
                                    crispr = "Chronos gene effect (0 = no effect, -1 = median of common essential genes)", protein = "relative abundance (harmonised, log2 scale)",
                                    moffitt = "mean percentile of basal-like minus classical genes (percentage points)")))
stopifnot(!anyDuplicated(out$models$model), identical(colnames(out$pct_models_x10), out$models$model), identical(colnames(out$pct_tissue_x10), out$tissue_units$unit))
saveRDS(out, file.path(IN_DIR, "models.rds"))
logc(log_file, "wrote ", file.path(IN_DIR, "models.rds"), " (", round(file.size(file.path(IN_DIR, "models.rds")) / 1e6, 1), " MB): ", nrow(models), " models (",
     paste(names(table(paste(models$source, models$type))), table(paste(models$source, models$type)), sep = " = ", collapse = "; "), ")")
end_log(log_file)
