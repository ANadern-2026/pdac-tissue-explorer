# ---- Retire the Connect Cloud address (run after the relocation period, e.g. in November 2026) ----------------
# Replaces the full app at the SAME Connect Cloud address with a tiny page that forwards visitors to
# https://pdac.biointeract.space/ . It uploads no data and needs very little memory, so it uses almost no credits.
# Run from the repository root, on the computer that published the app (it reuses app/rsconnect/, the record of
# that publication, so the update lands on the existing address instead of creating a new one).

if (!dir.exists(file.path("app", "rsconnect")))
  stop("No publication record in app/rsconnect/. Run this on the computer that published the app with deploy.R.")

tmp <- file.path(tempdir(), "pdac-redirect")
unlink(tmp, recursive = TRUE); dir.create(tmp)
file.copy(file.path("deploy", "connect-cloud", "redirect", "app.R"), tmp)
file.copy(file.path("app", "rsconnect"), tmp, recursive = TRUE)     # same target as the full app

rsconnect::deployApp(appDir = tmp, appFiles = "app.R", appTitle = "PDAC Tissue Explorer", forceUpdate = TRUE, launch.browser = TRUE)

# Afterwards, in Connect Cloud (the app's Settings > Compute): set memory and CPU to the smallest option.
