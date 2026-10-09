# ---- Deploy the PDAC Tissue Explorer to Posit Connect Cloud ------------------------------------------------
# Free plan: public apps, 4 GB RAM, 1 CPU; uploads from R up to 1 GiB (this app: ~240 MB).
# The data files are not in Git, so the app is published from R (not from the GitHub repository).
# Run from the repository root (the folder that contains app/ and deploy/), e.g. in RStudio.

if (!requireNamespace("rsconnect", quietly = TRUE)) install.packages("rsconnect")
if (!exists("connectCloudUser", envir = asNamespace("rsconnect")))
  stop("This rsconnect version cannot publish to Connect Cloud; run install.packages('rsconnect') to update it.")

# 1. the six prepared data files must be present (built by data-raw/ or downloaded from Zenodo)
needed  <- c("bulk_tcga_gtex.rds", "cptac.rds", "genes.rds", "sc_expr.rds", "sc_groups.rds", "sc_umap.rds")
missing <- needed[!file.exists(file.path("app", "data", needed))]
if (length(missing)) stop("Missing in app/data/: ", paste(missing, collapse = ", "))

# 2. the app must start locally (also confirms that all packages are installed)
#    shiny::runApp("app")   # optional check; close the app window before continuing

# 3. one-time: link this computer to your Connect Cloud account (opens the browser to authorise)
if (!nrow(rsconnect::accounts()) || !any(grepl("posit\\.cloud", rsconnect::accounts()$server)))
  rsconnect::connectCloudUser()

# 4. publish exactly the app files (app.R, R/, www/, data/); re-run this line to publish an update
app_files <- list.files("app", recursive = TRUE)
app_files <- app_files[!grepl("(^|/)(\\.|rsconnect/)", app_files)]
rsconnect::deployApp(appDir = "app", appFiles = app_files, appTitle = "PDAC Tissue Explorer", launch.browser = TRUE)

# After the first publish, in the Connect Cloud web page of the app:
#   - Settings > Sharing: make sure access is "Public".
#   - Settings > URL / custom name: choose "pdac-tissue-explorer" if available.
#   - Copy the final public URL into app/R/00_config.R (APP_URL), README.md and CITATION.cff, then publish again.
