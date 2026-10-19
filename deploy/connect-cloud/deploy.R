# ---- Deploy the PDAC Tissue Explorer to Posit Connect Cloud (legacy address; see also retire.R) ------------
# Free plan: public apps, 4 GB RAM, 1 CPU; uploads from R up to 1 GiB (this app: ~280 MB).
# The data files are not in Git, so the app is published from R (not from the GitHub repository).
# Run from the repository root (the folder that contains app/ and deploy/), e.g. in RStudio.

if (!requireNamespace("rsconnect", quietly = TRUE)) install.packages("rsconnect")
if (!exists("connectCloudUser", envir = asNamespace("rsconnect")))
  stop("This rsconnect version cannot publish to Connect Cloud; run install.packages('rsconnect') to update it.")

# 1. the gene index and the data store must be present (get_data.R, or data-raw/ incl. 05_build_store.R)
needed  <- c("genes.rds", "cptac.sqlite", "bulk.sqlite", "sc.sqlite")
missing <- needed[!file.exists(file.path("app", "data", needed))]
if (length(missing)) stop("Missing in app/data/: ", paste(missing, collapse = ", "))

# 2. the app must start locally (also confirms that all packages are installed)
#    shiny::runApp("app")   # optional check; close the app window before continuing

# 3. one-time: link this computer to your Connect Cloud account (opens the browser to authorise)
if (!nrow(rsconnect::accounts()) || !any(grepl("posit\\.cloud", rsconnect::accounts()$server)))
  rsconnect::connectCloudUser()

# 4. relocation notice: the app's permanent home is now https://pdac.biointeract.space/ (deploy/server/).
#    With MOVED_TO set, this copy shows a banner and a dialog pointing there. Set MOVED_TO <- "" to publish without it.
MOVED_TO   <- "https://pdac.biointeract.space/"
MOVED_NOTE <- "This old address will stop working in November 2026."
if (nzchar(MOVED_TO)) writeLines(c(MOVED_TO, MOVED_NOTE), file.path("app", "MOVED_TO"))   # removed again after publishing

# 5. publish exactly the app files (app.R, R/, www/, data store, MOVED_TO); re-run to publish an update
app_files <- list.files("app", recursive = TRUE)
app_files <- app_files[!grepl("(^|/)(\\.|rsconnect/)", app_files)]
app_files <- app_files[!grepl("^data/", app_files) | app_files %in% file.path("data", needed)]   # not the source .rds tables
rsconnect::deployApp(appDir = "app", appFiles = app_files, appTitle = "PDAC Tissue Explorer", launch.browser = TRUE)
if (nzchar(MOVED_TO)) unlink(file.path("app", "MOVED_TO"))

# After the first publish, in the Connect Cloud web page of the app:
#   - Settings > Sharing: make sure access is "Public".
#   - Settings > URL / custom name: choose "pdac-tissue-explorer" if available.
#   - Copy the final public URL into app/R/00_config.R (APP_URL), README.md and CITATION.cff, then publish again.
