# ---- usage counting: anonymous visits and in-app usage -------------------------------------------------------
# Active only when the environment variable PDAC_USAGE_DB names a database file (set on the server in compose.yaml),
# so local runs in RStudio record nothing. No cookies, no IP addresses, no user-agent strings are stored:
#  - a visit is a random id created when a browser opens a live session; only browsers that run the page's script are counted;
#  - per visit: start/end time, the referring site (host name only), a phone/desktop class from the window width,
#    and the link tag if the address had ?ref=... (e.g. https://<domain>/?ref=linkedin in a post);
#  - events: pages and gene tabs opened, genes looked up, gene lists compared, gene sets, downloads (button id only).
# The private stats page needs a secret key in the address (?stats=<PDAC_STATS_KEY>) and PDAC_STATS_PASSWORD; any other
# address shows the normal site with no hint of the page. Sessions that unlock it are not counted.
# Every write is wrapped in tryCatch: counting must never break the app.

USAGE_DB <- Sys.getenv("PDAC_USAGE_DB", "")
STATS_PASSWORD <- Sys.getenv("PDAC_STATS_PASSWORD", "")
STATS_KEY <- Sys.getenv("PDAC_STATS_KEY", "")
USAGE_ON <- nzchar(USAGE_DB)
STATS_ON <- USAGE_ON && nzchar(STATS_PASSWORD) && nchar(STATS_KEY) >= 16     # short keys could be guessed: page stays off
# same-length comparison that does not stop at the first differing character
secret_equal <- function(x, secret) {
  x <- if (is.null(x) || !length(x) || is.na(x[1])) "" else enc2utf8(as.character(x[1])); a <- utf8ToInt(x); b <- utf8ToInt(enc2utf8(secret))
  if (length(a) != length(b) || !length(b)) return(FALSE)
  sum(bitwXor(a, b)) == 0
}
USAGE_TZ <- "Europe/Berlin"
USAGE_MAX_VALUE <- 300L      # characters kept per event value (gene lists are cut)

.usage <- new.env()
usage_con <- function() {
  if (!USAGE_ON) return(NULL)
  if (!is.null(.usage$con) && DBI::dbIsValid(.usage$con)) return(.usage$con)
  con <- tryCatch({
    c0 <- DBI::dbConnect(RSQLite::SQLite(), USAGE_DB)
    DBI::dbExecute(c0, "PRAGMA journal_mode = WAL"); DBI::dbExecute(c0, "PRAGMA busy_timeout = 2000")
    DBI::dbExecute(c0, "CREATE TABLE IF NOT EXISTS visits (sid TEXT PRIMARY KEY, start INTEGER NOT NULL, last INTEGER, end INTEGER, referrer TEXT, device TEXT, tag TEXT, excluded INTEGER NOT NULL DEFAULT 0)")
    DBI::dbExecute(c0, "CREATE TABLE IF NOT EXISTS events (sid TEXT NOT NULL, t INTEGER NOT NULL, kind TEXT NOT NULL, value TEXT)")
    DBI::dbExecute(c0, "CREATE INDEX IF NOT EXISTS events_sid ON events (sid)")
    DBI::dbExecute(c0, "CREATE INDEX IF NOT EXISTS visits_start ON visits (start)")
    c0
  }, error = function(e) { message("usage counting off: ", conditionMessage(e)); NULL })
  .usage$con <- con
  con
}
usage_close <- function() if (!is.null(.usage$con)) { try(DBI::dbDisconnect(.usage$con), silent = TRUE); .usage$con <- NULL }
usage_now <- function() as.integer(Sys.time())
usage_exec <- function(sql, params) {
  con <- usage_con(); if (is.null(con)) return(invisible(FALSE))
  tryCatch({ DBI::dbExecute(con, sql, params = params); invisible(TRUE) },
           error = function(e) { message("usage write failed: ", conditionMessage(e)); invisible(FALSE) })
}

# referrer -> host name only ("" = typed/bookmarked or hidden by the browser; own domain = internal)
referrer_host <- function(ref, own = NULL) {
  if (is.null(ref) || !length(ref) || is.na(ref[1]) || !nzchar(ref[1])) return("(direct)")
  h <- tolower(sub("^[a-z][a-z0-9+.-]*://([^/:?#]+).*$", "\\1", ref[1], ignore.case = TRUE))
  if (identical(h, tolower(ref[1])) || !grepl("^[a-z0-9.-]+$", h)) return("(other)")
  h <- sub("^www\\.", "", h)
  if (!is.null(own) && nzchar(own) && h == sub("^www\\.", "", tolower(own))) return("(direct)")
  # group the many hosts of large sites
  if (grepl("(^|\\.)google\\.[a-z.]+$", h)) return("google")
  if (h %in% c("lnkd.in") || grepl("(^|\\.)linkedin\\.com$", h)) return("linkedin.com")
  if (h %in% c("t.co", "x.com", "twitter.com")) return("x.com")
  if (grepl("(^|\\.)bing\\.com$", h)) return("bing.com")
  if (grepl("(^|\\.)bsky\\.(app|social)$", h)) return("bsky.app")
  h
}
usage_device <- function(w) { w <- suppressWarnings(as.numeric(w)); if (!length(w) || is.na(w)) NA_character_ else if (w < 768) "phone" else if (w < 1100) "tablet / small window" else "desktop" }

# one tracker per Shiny session; the visit row is written only when the browser reports in (real browser with scripts)
usage_tracker <- function(session) {
  # Shiny's session token is random and unique (drawn from Shiny's private seed); sample() here would repeat ids, because the bootstrap code sets fixed seeds
  sid <- substr(gsub("[^0-9a-f]", "", tolower(session$token)), 1, 20)
  if (nchar(sid) < 8) sid <- paste0(format(Sys.time(), "%H%M%OS6"), "-", Sys.getpid())
  st <- new.env(); st$on <- FALSE; st$excluded <- FALSE
  list(
    start = function(client) {
      if (!USAGE_ON || st$on || st$excluded) return(invisible())
      if (isTRUE(client$nocount)) { st$excluded <- TRUE; return(invisible()) }
      own <- session$clientData$url_hostname
      tag <- if (length(client$tag) && nzchar(client$tag)) substr(gsub("[^A-Za-z0-9_.-]", "", as.character(client$tag)), 1, 40) else NA_character_
      now <- usage_now()
      st$on <- isTRUE(usage_exec("INSERT OR IGNORE INTO visits (sid, start, last, referrer, device, tag) VALUES (?, ?, ?, ?, ?, ?)",
                                 list(sid, now, now, referrer_host(client$ref, own), usage_device(client$w), tag)))
    },
    log = function(kind, value = NA_character_) {
      if (!st$on || st$excluded) return(invisible())
      v <- if (is.null(value) || !length(value)) NA_character_ else substr(paste(as.character(value), collapse = ","), 1, USAGE_MAX_VALUE)
      key <- paste(kind, v); if (kind != "download" && identical(key, st$last)) return(invisible()); st$last <- key   # repeats of the same event
      now <- usage_now()
      usage_exec("INSERT INTO events (sid, t, kind, value) VALUES (?, ?, ?, ?)", list(sid, now, kind, v))
      usage_exec("UPDATE visits SET last = ? WHERE sid = ?", list(now, sid))
    },
    exclude = function() {        # the owner unlocked the stats page: drop this visit
      st$excluded <- TRUE
      if (st$on) { usage_exec("DELETE FROM events WHERE sid = ?", list(sid)); usage_exec("DELETE FROM visits WHERE sid = ?", list(sid)) }
      st$on <- FALSE
    },
    end = function() if (st$on && !st$excluded) usage_exec("UPDATE visits SET end = ? WHERE sid = ?", list(usage_now(), sid)),
    active = function() st$on
  )
}

# client side: report referrer/width once connected; count download clicks; "don't count this browser" flag (owner only)
usage_js <- function() tags$script(HTML("
(function() {
  function nocount() { try { return window.localStorage.getItem('pdac_nocount') === '1'; } catch (e) { return false; } }
  var sent = false;
  function report() {
    if (sent || !(window.Shiny && Shiny.setInputValue)) return; sent = true;
    var q = new URLSearchParams(window.location.search);
    Shiny.setInputValue('usage_client', { ref: document.referrer || '', w: window.innerWidth, nocount: nocount(),
      tag: q.has('stats') ? '' : (q.get('ref') || q.get('utm_source') || '') }, { priority: 'event' });
  }
  $(document).on('shiny:connected', report);
  $(function() { if (window.Shiny && Shiny.shinyapp && Shiny.shinyapp.isConnected && Shiny.shinyapp.isConnected()) report(); });
  $(document).on('click', 'a.shiny-download-link', function() { Shiny.setInputValue('usage_dl', this.id, { priority: 'event' }); });
  Shiny.addCustomMessageHandler('usage_nocount', function(on) {
    try { if (on) window.localStorage.setItem('pdac_nocount', '1'); else window.localStorage.removeItem('pdac_nocount'); } catch (e) {}
  });
})();"))

# ---- reading the counts (stats page) -------------------------------------------------------------------------
usage_read <- function(days = NULL) {
  con <- usage_con(); if (is.null(con)) return(NULL)
  from <- if (is.null(days)) 0L else usage_now() - as.integer(days) * 86400L
  v <- DBI::dbGetQuery(con, "SELECT * FROM visits WHERE start >= ? AND excluded = 0", params = list(from))
  e <- DBI::dbGetQuery(con, "SELECT e.* FROM events e JOIN visits v ON v.sid = e.sid WHERE v.start >= ? AND v.excluded = 0", params = list(from))
  v$dur <- pmax(0, ifelse(is.na(v$end), v$last, v$end) - v$start)
  list(visits = v, events = e, from = from)
}

USAGE_DL_NAMES <- c(cptac = "Tumour vs normal (CPTAC)", bulk = "Tumour vs normal (TCGA/GTEx)", phospho = "Phosphosites", surv = "Survival",
                    stage = "Stage", moffitt = "Classical/basal", ct = "Cell types", states = "Tissue states", fibro = "Fibroblast subtypes",
                    ductal = "Tumour-cell states", coexp = "Co-expression", umap = "Atlas map", sets = "Gene sets & survival",
                    mdot = "Compare: cell types", mforest = "Compare: tumour vs normal")

usage_kpis <- function(u) {
  v <- u$visits; e <- u$events
  n <- nrow(v); used <- unique(e$sid[e$kind %in% c("gene", "compare", "geneset", "download")])
  list(visits = n,
       days = if (n) length(unique(as.Date(as.POSIXct(v$start, origin = "1970-01-01", tz = USAGE_TZ)))) else 0L,
       median_min = if (n) stats::median(v$dur) / 60 else NA_real_,
       engaged = if (n) mean(v$sid %in% used) else NA_real_,
       genes = length(unique(e$value[e$kind == "gene"])),
       downloads = sum(e$kind == "download"))
}

usage_daily <- function(u, days) {
  tz <- USAGE_TZ
  end <- as.Date(Sys.time(), tz = tz)
  start <- if (is.null(days)) { if (nrow(u$visits)) min(as.Date(as.POSIXct(u$visits$start, origin = "1970-01-01", tz = tz), tz = tz)) else end } else end - as.integer(days) + 1L
  d <- seq(start, end, by = "day")
  vd <- as.Date(as.POSIXct(u$visits$start, origin = "1970-01-01", tz = tz), tz = tz)
  data.frame(date = d, visits = as.integer(table(factor(as.character(vd), levels = as.character(d)))))
}

usage_count <- function(x, label, n = NULL, total = NULL) {
  if (!length(x)) { out <- data.frame(character(0), integer(0)); names(out) <- c(label, "Count"); return(out) }
  t <- sort(table(x), decreasing = TRUE)
  out <- data.frame(names(t), as.integer(t), check.names = FALSE, stringsAsFactors = FALSE); names(out) <- c(label, "Count")
  if (!is.null(total) && total > 0) out$Share <- sprintf("%.0f%%", 100 * out$Count / total)
  if (!is.null(n)) out <- utils::head(out, n)
  out
}

usage_tables <- function(u) {
  v <- u$visits; e <- u$events; n <- nrow(v)
  first <- function(kind) { x <- e[e$kind == kind, , drop = FALSE]; unique(x[, c("sid", "value")]) }  # once per visit
  dl <- e$value[e$kind == "download"]
  dl_type <- ifelse(grepl("_csv$", dl), "Data (CSV)", ifelse(grepl("_png$", dl), "Figure (PNG)", ifelse(grepl("_pdf$", dl), "Figure (PDF)", "Table copy/CSV")))
  dl_what <- sub("^dl_", "", sub("_(pdf|png|csv)$", "", dl))
  dl_what <- ifelse(dl_what %in% names(USAGE_DL_NAMES), USAGE_DL_NAMES[dl_what], dl_what)
  list(
    referrers = usage_count(v$referrer, "Came from", total = n),
    tags = usage_count(v$tag[!is.na(v$tag)], "Link tag (?ref=)", total = n),
    devices = usage_count(ifelse(is.na(v$device), "unknown", v$device), "Screen", total = n),
    pages = usage_count(unname(c(gene = "Explore a gene", multi = "Compare genes", sets = "Gene sets & survival", about = "About")[first("page")$value]) |> (\(x) x[!is.na(x)])(),
                        "Page opened (visits)", total = n),
    tabs = usage_count(first("tab")$value, "Gene tab opened (visits)", total = n),
    genes = usage_count(first("gene")$value, "Gene", n = 100),
    downloads = if (length(dl)) usage_count(paste0(dl_what, " - ", dl_type), "Download") else usage_count(character(0), "Download"),
    sets = usage_count(e$value[e$kind == "geneset"], "Gene set", n = 50),
    compare = usage_count(e$value[e$kind == "compare"], "Gene list compared", n = 50))
}

plot_usage_daily <- function(d) {
  ggplot(d, aes(date, visits, text = paste0(format(date, "%a %d %b %Y"), "<br>", visits, " visit", ifelse(visits == 1, "", "s")))) +
    geom_col(fill = ONE_SERIES, width = 0.8) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08)), breaks = function(l) unique(floor(pretty(c(0, l[2]))))) +
    labs(x = NULL, y = "Visits", title = "Visits per day", subtitle = paste("Days in", USAGE_TZ, "time; your own unlocked visits are not counted")) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.grid.major.y = element_line(colour = INK[["grid"]]),
          axis.text = element_text(colour = INK[["secondary"]]), plot.title = element_text(face = "bold"))
}
