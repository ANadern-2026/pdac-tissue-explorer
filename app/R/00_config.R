# ---- configuration: labels, palettes, constants -------------------------------------------------------
APP_TITLE   <- "PDAC Tissue Explorer"
APP_VERSION <- "0.1.0"
MIN_CELLS   <- 20L          # a sample contributes to a cell-type summary only with >= this many cells of that type
BOOT_N      <- 1000L        # bootstrap resamples for confidence intervals of median changes

# validated categorical palette (dataviz reference instance, light mode); fixed order, never cycled
STATE_LEVELS <- c("Healthy donor", "Adjacent normal", "Primary tumour", "Liver metastasis", "Primary tumour (treated)")
STATE_COLS   <- c("Healthy donor" = "#1baf7a", "Adjacent normal" = "#2a78d6", "Primary tumour" = "#eb6834",
                  "Liver metastasis" = "#4a3aa7", "Primary tumour (treated)" = "#e87ba4")
STATE_SHORT  <- c("Healthy donor" = "Donor", "Adjacent normal" = "Adj. normal", "Primary tumour" = "Primary", "Liver metastasis" = "Liver met", "Primary tumour (treated)" = "Treated")
BULK_LEVELS  <- c("Healthy pancreas (GTEx)", "Adjacent normal (TCGA)", "PDAC (TCGA)")
BULK_COLS    <- c("Healthy pancreas (GTEx)" = "#1baf7a", "Adjacent normal (TCGA)" = "#2a78d6", "PDAC (TCGA)" = "#eb6834")
TISSUE_COLS  <- c("Adjacent normal" = "#2a78d6", "Tumour" = "#eb6834")
SEQ_RAMP     <- c("#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b")   # single-hue blue, light -> dark
DIV_COLS     <- c(low = "#2a78d6", mid = "#f0efec", high = "#e34948")                 # blue <-> red, grey midpoint
ONE_SERIES   <- "#2a78d6"
INK <- c(primary = "#0b0b0b", secondary = "#52514e", muted = "#8a8984", grid = "#e9e8e4", surface = "#fcfcfb", zero = "#d9d8d4")

CELLTYPE_LABELS <- c(FIBROBLASTS = "Fibroblasts", PERICYTES = "Pericytes", ENDOTHELIAL = "Endothelial cells",
                     TUMOUR_DUCTAL = "Ductal: tumour subclusters", OTHER_DUCTAL = "Ductal: other",
                     ACINAR = "Acinar cells", ENDOCRINE = "Endocrine cells", MYELOID = "Myeloid cells", MAST = "Mast cells",
                     TNK = "T and NK cells", `B CELLS` = "B cells", PLASMA = "Plasma cells")
CELLTYPE_ORDER <- unname(CELLTYPE_LABELS)
CELLTYPE_PROSE <- setNames(c("fibroblasts", "pericytes", "endothelial cells", "ductal cells of the tumour subclusters", "other ductal cells", "acinar cells",
                             "endocrine cells", "myeloid cells", "mast cells", "T and NK cells", "B cells", "plasma cells"), CELLTYPE_ORDER)
FIBRO_ORDER  <- c(paste("myFb", 1:5), paste("iFb", 1:3))
FIBRO_NOTE   <- c("myFb 5" = "CXCL10+ myofibroblasts")
DUCTAL_GROUP <- c("Ductal 1" = "Classical-like", "Ductal 2" = "Classical-like", "Ductal 3" = "Classical-like", "Ductal 9" = "Classical-like",
                  "Ductal 4" = "Basal-like", "Ductal 7" = "Basal-like", "Ductal 5" = "Hybrid", "Ductal 6" = "Hybrid", "Ductal 8" = "Secretory")
DUCTAL_GROUP_ORDER <- c("Classical-like", "Hybrid", "Basal-like", "Secretory")

STUDY_CITATION <- c(
  GSE154778 = "Lin et al., Genome Med 2020", GSE155698 = "Steele et al., Nat Cancer 2020", GSE156405 = "Lee et al., Clin Cancer Res 2021",
  GSE158356 = "Kemp et al., Life Sci Alliance 2021", GSE194247 = "Unpublished at atlas release", GSE202051 = "Shiau et al., Front Oncol 2022",
  GSE205013 = "Werba et al., Nat Commun 2023", GSE211644 = "Schalck et al., Cancer Discov 2022", GSE229413 = "Carpenter et al., Cancer Discov 2023",
  EGAS00001002543 = "Chan-Seng-Yue et al., Nat Genet 2020", phs001840.v1.p1 = "Elyada et al., Cancer Discov 2019", PRJCA001063 = "Peng et al., Cell Res 2019")

DATA_DIR <- Sys.getenv("PDACX_DATA_DIR", "data")
REPO_URL <- "https://github.com/ANadern-2026/pdac-tissue-explorer"
APP_URL  <- "https://pdac.biointeract.space/"   # self-hosted (deploy/server/)
AUTHOR   <- "Ali Nadernezhad (Leibniz-Institut für Polymerforschung Dresden)"
ORCID    <- "https://orcid.org/0000-0002-3761-1601"
