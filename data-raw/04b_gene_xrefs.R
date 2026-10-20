# =====================================================================================================
# 04b_gene_xrefs.R
# Cross-references for the links out on the app's gene page: Ensembl gene ID and gene type for every gene of
# the search index, via its Entrez ID (org.Hs.eg.db). An Ensembl ID is kept only when the Entrez -> Ensembl
# mapping is unique; otherwise the app links to a search for the symbol. genes.rds itself is not changed.
# Run after 04:  source("data-raw/04b_gene_xrefs.R")      Requires AnnotationDbi + org.Hs.eg.db.
# =====================================================================================================

DATA_DIR <- Sys.getenv("PDACX_DATA_DIR", file.path("app", "data"))
source(file.path("data-raw", "helpers.R")); log_file <- start_log("04b_gene_xrefs")
suppressPackageStartupMessages({ library(AnnotationDbi); library(org.Hs.eg.db) })

gfile <- file.path(DATA_DIR, "genes.rds"); genes <- readRDS(gfile)$genes
eg <- sort(unique(stats::na.omit(genes$entrez)))
m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = eg, keytype = "ENTREZID", columns = c("ENSEMBL", "GENETYPE")))
ens <- unique(m[!is.na(m$ENSEMBL), c("ENTREZID", "ENSEMBL")])
n_ens <- table(ens$ENTREZID); ens1 <- ens[ens$ENTREZID %in% names(n_ens)[n_ens == 1], ]
typ <- unique(m[!is.na(m$GENETYPE), c("ENTREZID", "GENETYPE")]); typ <- typ[!duplicated(typ$ENTREZID), ]

xrefs <- data.frame(symbol = genes$symbol, entrez = genes$entrez,
                    ensembl = ens1$ENSEMBL[match(genes$entrez, ens1$ENTREZID)],
                    gene_type = typ$GENETYPE[match(genes$entrez, typ$ENTREZID)], stringsAsFactors = FALSE)
logc(log_file, "genes: ", nrow(xrefs), " | with Entrez ID: ", sum(!is.na(xrefs$entrez)), " | unique Ensembl ID: ", sum(!is.na(xrefs$ensembl)),
     " | several Ensembl IDs (search link): ", sum(xrefs$entrez %in% names(n_ens)[n_ens > 1]))
logc(log_file, "gene types: ", paste(names(t <- sort(table(xrefs$gene_type, useNA = "ifany"), decreasing = TRUE)), t, sep = " = ", collapse = "; "))

out <- file.path(DATA_DIR, "gene_xrefs.rds")
saveRDS(list(xrefs = xrefs, annotation = paste0("org.Hs.eg.db ", as.character(packageVersion("org.Hs.eg.db"))),
             genes_md5 = unname(tools::md5sum(gfile)), built = format(Sys.time(), tz = "UTC")), out, compress = "xz")
logc(log_file, "wrote ", out, " (", round(file.size(out) / 1e6, 2), " MB)")
end_log(log_file)
