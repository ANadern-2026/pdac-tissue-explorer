# =====================================================================================================
# 04_gene_index.R
# Builds the searchable gene index used by the app: one row per gene (Entrez-level where possible) with
# the current official symbol, full name, aliases, and the symbol under which it appears in each dataset
# (CPTAC protein, CPTAC RNA, TCGA/GTEx, single-cell atlas). Datasets use different annotation releases,
# so symbols are reconciled through Entrez IDs (official symbol first, unambiguous alias second).
# Run after 01-03:  source("data-raw/04_gene_index.R")      Requires AnnotationDbi + org.Hs.eg.db.
# =====================================================================================================

DATA_DIR <- Sys.getenv("PDACX_DATA_DIR", file.path("app", "data"))
source(file.path("data-raw", "helpers.R")); log_file <- start_log("04_gene_index")
suppressPackageStartupMessages({ library(AnnotationDbi); library(org.Hs.eg.db) })

cptac <- readRDS(file.path(DATA_DIR, "cptac.rds")); bulk <- readRDS(file.path(DATA_DIR, "bulk_tcga_gtex.rds"))
sc_file <- file.path(DATA_DIR, "sc_expr.rds")
sets <- list(cptac_protein = rownames(cptac$protein_tumor), cptac_rna = rownames(cptac$rna_tumor),
             tcga_gtex = rownames(bulk$log2tpm_x1000), sc = if (file.exists(sc_file)) readRDS(sc_file)$genes else character(0))
if (!length(sets$sc)) logc(log_file, "WARNING: ", sc_file, " not found - single-cell genes not indexed yet")

sym2eg <- AnnotationDbi::toTable(org.Hs.egSYMBOL); alias2eg <- AnnotationDbi::toTable(org.Hs.egALIAS2EG)
uniq_alias <- alias2eg[!alias2eg$alias_symbol %in% alias2eg$alias_symbol[duplicated(alias2eg$alias_symbol)], ]
to_eg <- function(s) {
  eg <- sym2eg$gene_id[match(s, sym2eg$symbol)]
  miss <- is.na(eg); eg[miss] <- uniq_alias$gene_id[match(s[miss], uniq_alias$alias_symbol)]
  eg
}
all_sym <- sort(unique(unlist(sets)))
idx <- data.frame(dataset_symbol = all_sym, entrez = to_eg(all_sym), stringsAsFactors = FALSE)
idx$symbol <- sym2eg$symbol[match(idx$entrez, sym2eg$gene_id)]
idx$symbol[is.na(idx$symbol)] <- idx$dataset_symbol[is.na(idx$symbol)]
idx$key <- ifelse(is.na(idx$entrez), paste0("sym:", idx$dataset_symbol), paste0("eg:", idx$entrez))

genes <- do.call(rbind, lapply(split(idx, idx$key), function(d) {
  pick <- function(set) { s <- intersect(d$dataset_symbol, sets[[set]]); if (!length(s)) NA_character_ else if (d$symbol[1] %in% s) d$symbol[1] else s[1] }
  data.frame(key = d$key[1], symbol = d$symbol[1], entrez = d$entrez[1],
             sym_cptac_protein = pick("cptac_protein"), sym_cptac_rna = pick("cptac_rna"), sym_tcga_gtex = pick("tcga_gtex"), sym_sc = pick("sc"),
             stringsAsFactors = FALSE)
}))
genes$name <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = ifelse(is.na(genes$entrez), "0", genes$entrez), keytype = "ENTREZID", column = "GENENAME", multiVals = "first")
genes$name[is.na(genes$entrez)] <- NA
al <- split(alias2eg$alias_symbol, alias2eg$gene_id)
genes$aliases <- vapply(seq_len(nrow(genes)), function(i) {
  a <- unique(c(if (!is.na(genes$entrez[i])) al[[genes$entrez[i]]], idx$dataset_symbol[idx$key == genes$key[i]]))
  paste(setdiff(a, genes$symbol[i]), collapse = ", ") }, character(1))
genes$n_datasets <- rowSums(!is.na(genes[, c("sym_cptac_protein", "sym_cptac_rna", "sym_tcga_gtex", "sym_sc")]))
genes <- genes[genes$n_datasets > 0, ]
# a symbol shown to users must be unique
dup <- duplicated(genes$symbol) | duplicated(genes$symbol, fromLast = TRUE)
genes$symbol[dup] <- paste0(genes$symbol[dup], " (", ifelse(is.na(genes$entrez[dup]), genes$sym_tcga_gtex[dup], paste0("Entrez ", genes$entrez[dup])), ")")
genes <- genes[order(-genes$n_datasets, genes$symbol), ]; rownames(genes) <- NULL

logc(log_file, "indexed genes: ", nrow(genes), " | in all four datasets: ", sum(genes$n_datasets == 4),
     " | reconciled symbol differences: ", sum(apply(genes[, 4:7], 1, function(r) length(unique(na.omit(r))) > 1)))
logc(log_file, "per dataset: ", paste(names(sets), lengths(sets), sep = " = ", collapse = "; "))
saveRDS(list(genes = genes, annotation = paste0("org.Hs.eg.db ", as.character(packageVersion("org.Hs.eg.db")))), file.path(DATA_DIR, "genes.rds"), compress = "xz")
logc(log_file, "wrote ", file.path(DATA_DIR, "genes.rds"))
end_log(log_file)
