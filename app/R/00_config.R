# ---- configuration: labels, palettes, constants -------------------------------------------------------
APP_TITLE   <- "PDAC Tissue Explorer"
APP_VERSION <- "0.3.0"
MIN_CELLS   <- 20L          # a sample contributes to a cell-type summary only with >= this many cells of that type
BOOT_N      <- 1000L        # bootstrap resamples for confidence intervals of median changes
PHOSPHO_MIN_PAIRS <- 10L    # a phosphosite is compared only if quantified in tumour and normal tissue of >= this many patients
PHOSPHO_MAX_PLOT  <- 50L    # the site overview plot shows at most this many sites (those with the smallest p-values)
COEXP_STATE       <- "Primary tumour"   # co-expression across samples: treatment-naive primary tumours only (data-raw/01b)
COEXP_MIN_STUDY   <- 3L     # ... from studies contributing at least this many samples of the cell type (study effects are removed)
COEXP_MIN_SAMPLES <- 20L    # a cell type is analysed only with at least this many such samples
COEXP_MIN_FRAC    <- 0.02   # a gene is analysed in a cell type only if detected in >= 2% of its cells in at least half of the samples
COEXP_TOP_POS     <- 50L    # partners stored per gene and cell type: the most positively ...
COEXP_TOP_NEG     <- 20L    # ... and the most negatively correlated genes
COEXP_SHOW        <- 25L    # partners shown in the co-expression plot (all stored partners are in the table)
SURV_MIN_GROUP    <- 10L    # Kaplan-Meier median split only if each group has at least this many patients
SURV_MIN_EVENTS   <- 10L    # survival models only with at least this many events
DAYS_PER_MONTH    <- 365.25 / 12

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

# ---- page menu icons: small two-colour drawings on a tinted tile, one colour per page (user 2026-11-13: larger,
# colourful, not plain symbols). Colours come from CSS (#nav a[data-value=...]: --c, --c2, --t); classes p/s/r/w below.
NAV_ICON_SVG <- c(
  gene  = '<path class="p" d="M7 2.5c0 5 10 4.5 10 9.5s-10 4.5-10 9.5" fill="none" stroke-width="2.3" stroke-linecap="round"/><path class="s" d="M17 2.5c0 5-10 4.5-10 9.5s10 4.5 10 9.5" fill="none" stroke-width="2.3" stroke-linecap="round"/><path class="r" d="M8.3 5.5h7.4M9.6 9h4.8M9.6 15h4.8M8.3 18.5h7.4" stroke-width="1.7" stroke-linecap="round"/>',
  multi = '<circle class="pf" cx="6" cy="6" r="2.6"/><circle class="sf" cx="12" cy="6" r="1.6"/><circle class="pf" cx="18" cy="6" r="2.1"/><circle class="sf" cx="6" cy="12" r="1.9"/><circle class="pf" cx="12" cy="12" r="2.8"/><circle class="sf" cx="18" cy="12" r="1.3"/><circle class="pf" cx="6" cy="18" r="1.4"/><circle class="sf" cx="12" cy="18" r="2.4"/><circle class="pf" cx="18" cy="18" r="2.7"/>',
  pops  = '<circle class="pf" cx="8.6" cy="9.4" r="5.4" opacity="0.9"/><circle class="sf" cx="16" cy="8.4" r="4.3" opacity="0.9"/><circle class="pf2" cx="13.4" cy="16.3" r="4.9" opacity="0.85"/><circle class="w" cx="8.3" cy="9.2" r="1.5"/><circle class="w" cx="16.2" cy="8.2" r="1.2"/><circle class="w" cx="13.4" cy="16.6" r="1.4"/>',
  sets  = '<path class="r" d="M4 3.5v16.5h16.5" fill="none" stroke-width="1.8" stroke-linecap="round"/><path class="p" d="M4.2 5h4v3.2h3.4v3.4h3v3.6h4.6" fill="none" stroke-width="2.3" stroke-linejoin="round" stroke-linecap="round"/><path class="s" d="M4.2 5h3v1.6h3.6v1.8h3.2v1.6h5" fill="none" stroke-width="2.3" stroke-linejoin="round" stroke-linecap="round"/>',
  models = '<path class="pf2" d="M7.2 14.5h9.6l2.6 4.7c.5 1-.2 2.3-1.4 2.3H6c-1.2 0-1.9-1.3-1.4-2.3z"/><path class="p" d="M9.5 2.8h5M10.3 3v6.3L4.6 19.2c-.6 1.1.2 2.3 1.4 2.3h12c1.2 0 2-1.2 1.4-2.3L13.7 9.3V3" fill="none" stroke-width="1.9" stroke-linejoin="round" stroke-linecap="round"/><circle class="sf" cx="10.3" cy="17.6" r="1.1"/><circle class="sf" cx="13.8" cy="16.2" r="0.8"/>',
  about = '<circle class="pf" cx="12" cy="12" r="9.5"/><circle class="w" cx="12" cy="7.4" r="1.5"/><path class="wl" d="M12 11v6.3" stroke-width="2.5" stroke-linecap="round"/>',
  usage = '<rect class="sf" x="4" y="12" width="3.6" height="8" rx="1"/><rect class="pf" x="10.2" y="6.5" width="3.6" height="13.5" rx="1"/><rect class="sf" x="16.4" y="9.5" width="3.6" height="10.5" rx="1"/>')
nav_icon <- function(page) htmltools::span(class = "nav-ico", htmltools::HTML(paste0('<svg viewBox="0 0 24 24" aria-hidden="true">', NAV_ICON_SVG[[page]], '</svg>')))
