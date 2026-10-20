# ---- configuration: labels, palettes, constants -------------------------------------------------------
APP_TITLE   <- "PDAC Tissue Explorer"
APP_VERSION <- "0.1.0"
MIN_CELLS   <- 20L          # a sample contributes to a cell-type summary only with >= this many cells of that type
BOOT_N      <- 1000L        # bootstrap resamples for confidence intervals of median changes
PHOSPHO_MIN_PAIRS <- 10L    # a phosphosite is compared only if quantified in tumour and normal tissue of >= this many patients
PHOSPHO_MAX_PLOT  <- 50L    # the site overview plot shows at most this many sites (those with the smallest p-values)

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

# links out on the gene page (gene_links() in 04_text.R): direct pages when the identifier is known, otherwise a search
LINKS_OUT <- list(
  ncbi = list(label = "NCBI Gene", hint = "Gene summary, function and literature (NCBI)",
              id = "https://www.ncbi.nlm.nih.gov/gene/%s", search = "https://www.ncbi.nlm.nih.gov/gene/?term=%s%%5Bsym%%5D+AND+human%%5Borgn%%5D"),
  uniprot = list(label = "UniProt", hint = "Protein function, domains, sites and isoforms (UniProt, reviewed human entries)",
                 search = "https://www.uniprot.org/uniprotkb?query=gene_exact%%3A%s+AND+organism_id%%3A9606+AND+reviewed%%3Atrue"),
  hpa = list(label = "Human Protein Atlas", hint = "Protein localisation in tissues, cells and cancers (antibody staining)",
             id = "https://www.proteinatlas.org/%s", search = "https://www.proteinatlas.org/search/%s"),
  opentargets = list(label = "Open Targets", hint = "Disease associations, genetics, drugs and tractability",
                     id = "https://platform.opentargets.org/target/%s", search = "https://platform.opentargets.org/search?q=%s"),
  depmap = list(label = "DepMap", hint = "Cancer cell-line dependency (CRISPR) and expression",
                search = "https://depmap.org/portal/gene/%s?tab=overview"))

DATA_DIR <- Sys.getenv("PDACX_DATA_DIR", "data")
REPO_URL <- "https://github.com/ANadern-2026/pdac-tissue-explorer"
APP_URL  <- "https://pdac.biointeract.space/"   # self-hosted (deploy/server/)
AUTHOR   <- "Ali Nadernezhad (Leibniz-Institut für Polymerforschung Dresden)"
ORCID    <- "https://orcid.org/0000-0002-3761-1601"
# Relocation notice: only when a file MOVED_TO sits next to app.R (written by deploy/connect-cloud/deploy.R for the
# old Connect Cloud address). Line 1 = new address; line 2 (optional) = extra sentence. Absent everywhere else.
MOVED_TO <- if (file.exists("MOVED_TO")) trimws(readLines("MOVED_TO", warn = FALSE)) else character(0)
