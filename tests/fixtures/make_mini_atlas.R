# Test fixture: a tiny Seurat object with the atlas's metadata layout, label files and UMAP, used to test
# data-raw/01_scatlas_summaries.R without the 28 GB atlas.  Usage: make_mini_atlas(out_dir)
make_mini_atlas <- function(out_dir) {
suppressPackageStartupMessages({ library(Matrix); library(SeuratObject) })
options(Seurat.object.assay.version = "v3")
set.seed(3); ng <- 400; nc <- 3000
counts <- rsparsematrix(ng, nc, density = 0.08, rand.x = function(n) rpois(n, 3) + 1)
rownames(counts) <- c("MSLN", "KRT19", "COL1A1", paste0("G", 4:ng)); colnames(counts) <- paste0("CELL", seq_len(nc), "-1_1")
counts@x <- round(counts@x)
lib <- colSums(counts); dat <- counts; dat@x <- log1p(dat@x / rep(lib, diff(dat@p)) * 1e4)
obj <- CreateSeuratObject(counts = counts)
obj <- SetAssayData(obj, layer = "data", new.data = dat)
ct <- sample(c("DUCTAL", "FIBROBLASTS", "TNK", "MYELOID", "ENDOTHELIAL"), nc, TRUE, prob = c(.4, .2, .2, .1, .1))
smp <- sample(paste0("S", 1:12), nc, TRUE)
ds <- c(S1 = "Primary tumor", S2 = "Primary tumor", S3 = "Primary tumor", S4 = "Primary tumor", S5 = "Primary tumor", S6 = "Primary tumor",
        S7 = "Adjacent normal", S8 = "Adjacent normal", S9 = "Metastatic lesion", S10 = "Metastatic lesion", S11 = "Donor", S12 = "Donor")
md <- data.frame(Name = smp, GSE.SRA..Study. = ifelse(smp %in% c("S1", "S2", "S9"), "GSE1", "GSE2"), DiseaseState = factor(ds[smp]),
                 If.metastatic..location = ifelse(ds[smp] == "Metastatic lesion", "Liver", ifelse(smp == "S3", "2", NA)),
                 Treatment = ifelse(ds[smp] == "Donor", "N_A", ifelse(smp == "S4", "Treatment", "Treatment na\u00efve")), Clusters = ct,
                 row.names = colnames(obj), stringsAsFactors = FALSE)
obj <- AddMetaData(obj, md)
obj@meta.data$nFeature_RNA[1:20] <- 5L
emb <- cbind(UMAP_1 = rnorm(nc), UMAP_2 = rnorm(nc)); rownames(emb) <- colnames(obj)
obj[["umap"]] <- CreateDimReducObject(embeddings = emb, key = "UMAP_", assay = "RNA")
saveRDS(obj, file.path(out_dir, "sc_ref_base.rds"))
cells <- colnames(obj)
dl <- data.frame(seurat_clusters = sample(paste("Ductal", 1:9), sum(ct == "DUCTAL"), TRUE), Name = smp[ct == "DUCTAL"], row.names = cells[ct == "DUCTAL"])
dl <- dl[-(1:30), ]                                  # some ductal cells unlabelled, as in the real atlas
fl <- data.frame(Coarse = sample(c(paste("myFb", 1:5), paste("iFb", 1:3)), sum(ct == "FIBROBLASTS"), TRUE), row.names = cells[ct == "FIBROBLASTS"])
write.table(dl, gzfile(file.path(out_dir, "DuctalTumor_Subcluster.txt.gz")), sep = "\t", quote = FALSE)
write.table(fl, gzfile(file.path(out_dir, "Fibroblast_Subcluster.txt.gz")), sep = "\t", quote = FALSE)
invisible(out_dir)
}
