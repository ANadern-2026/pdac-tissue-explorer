# Usage counting (app/R/05_usage.R): off without PDAC_USAGE_DB; referrer grouping; visits, events and exclusions in a temporary database.
root <- normalizePath(file.path(test_path(), "..", ".."))
load_usage <- function(db = "", pw = "") {
  old <- Sys.getenv(c("PDAC_USAGE_DB", "PDAC_STATS_PASSWORD"), unset = NA)
  Sys.setenv(PDAC_USAGE_DB = db, PDAC_STATS_PASSWORD = pw)
  on.exit(for (k in names(old)) if (is.na(old[[k]])) Sys.unsetenv(k) else do.call(Sys.setenv, stats::setNames(list(old[[k]]), k)))
  e <- new.env()
  for (f in sort(list.files(file.path(root, "app", "R"), full.names = TRUE))) sys.source(f, envir = e)
  e
}
fake_session <- function(token) list(token = token, clientData = list(url_hostname = "pdac.example.org"))

test_that("counting is off without PDAC_USAGE_DB (local runs)", {
  e <- load_usage()
  expect_false(e$USAGE_ON); expect_false(e$STATS_ON); expect_null(e$usage_con())
  tr <- e$usage_tracker(fake_session("abc123def4567890"))
  tr$start(list(ref = "https://www.linkedin.com/", w = 1400)); tr$log("gene", "KRT17")
  expect_false(tr$active())
  expect_false(load_usage(pw = "x")$STATS_ON)            # stats page needs both variables
})

test_that("referrers are reduced to grouped host names", {
  e <- load_usage()
  expect_equal(e$referrer_host(""), "(direct)"); expect_equal(e$referrer_host(NULL), "(direct)")
  expect_equal(e$referrer_host("https://www.linkedin.com/feed/update/xyz?trk=1"), "linkedin.com")
  expect_equal(e$referrer_host("https://lnkd.in/abc"), "linkedin.com")
  expect_equal(e$referrer_host("https://www.google.de/"), "google")
  expect_equal(e$referrer_host("https://t.co/abc"), "x.com")
  expect_equal(e$referrer_host("https://pdac.example.org/?gene=X", own = "pdac.example.org"), "(direct)")
  expect_equal(e$referrer_host("https://Scholar.Example.edu:8443/path"), "scholar.example.edu")
  expect_equal(e$referrer_host("not a url"), "(other)")
  expect_equal(e$usage_device(390), "phone"); expect_equal(e$usage_device(900), "tablet / small window")
  expect_equal(e$usage_device(1440), "desktop"); expect_true(is.na(e$usage_device(NULL)))
})

test_that("visits and events are stored, repeats skipped, exclusions removed", {
  db <- tempfile(fileext = ".sqlite"); e <- load_usage(db, "pw")
  expect_true(e$USAGE_ON); expect_false(e$STATS_ON)        # password without a key: stats page off
  t1 <- e$usage_tracker(fake_session("aaaaaaaaaaaaaaaaaaaaaaaa"))
  t1$log("gene", "EARLY")                                   # before the browser reports in: not counted
  t1$start(list(ref = "https://www.linkedin.com/", w = 1440, tag = "post<1>"))
  t1$log("gene", "KRT17"); t1$log("gene", "KRT17"); t1$log("tab", "Cell types"); t1$log("gene", "KRT17")
  t1$log("download", "dl_ct_csv"); t1$log("download", "dl_ct_csv"); t1$log("geneset", c("A", "B")); t1$end()
  t2 <- e$usage_tracker(fake_session("bbbbbbbbbbbbbbbbbbbbbbbb")); t2$start(list(ref = "", w = 390)); t2$log("page", "about")
  t3 <- e$usage_tracker(fake_session("cccccccccccccccccccccccc")); t3$start(list(ref = "", w = 1440)); t3$log("gene", "MSLN"); t3$exclude()
  t4 <- e$usage_tracker(fake_session("dddddddddddddddddddddddd")); t4$start(list(ref = "", w = 1440, nocount = TRUE)); t4$log("gene", "X")
  u <- e$usage_read()
  expect_equal(sort(u$visits$referrer), c("(direct)", "linkedin.com"))
  expect_equal(u$visits$tag[u$visits$referrer == "linkedin.com"], "post1")
  ev <- u$events[order(u$events$t, seq_len(nrow(u$events))), ]
  expect_equal(ev$value[ev$kind == "gene"], c("KRT17", "KRT17"))      # immediate repeat skipped, later lookup kept
  expect_equal(sum(ev$kind == "download"), 2L)                         # downloads always counted
  expect_equal(ev$value[ev$kind == "geneset"], "A,B")
  expect_false(any(c("MSLN", "X", "EARLY") %in% ev$value))
  k <- e$usage_kpis(u)
  expect_equal(k$visits, 2L); expect_equal(k$engaged, 0.5); expect_equal(k$genes, 1L); expect_equal(k$downloads, 2L)
  tb <- e$usage_tables(u)
  expect_equal(tb$genes$Count[tb$genes$Gene == "KRT17"], 1L)          # once per visit
  expect_equal(tb$downloads$Download, "Cell types - Data (CSV)")
  expect_equal(tb$pages$`Page opened (visits)`, "About")
  expect_equal(sum(e$usage_daily(u, 7)$visits), 2L); expect_equal(nrow(e$usage_daily(u, 7)), 7L)
  e$usage_close(); unlink(paste0(db, c("", "-wal", "-shm")))
})

test_that("stats page needs a long key and a password; secrets compared exactly", {
  e <- load_usage(tempfile(fileext = ".sqlite"), "pw")
  expect_false(e$STATS_ON)                                  # no key
  old <- Sys.getenv("PDAC_STATS_KEY", unset = NA); on.exit(if (is.na(old)) Sys.unsetenv("PDAC_STATS_KEY") else Sys.setenv(PDAC_STATS_KEY = old))
  Sys.setenv(PDAC_STATS_KEY = "short"); expect_false(load_usage(tempfile(), "pw")$STATS_ON)
  Sys.setenv(PDAC_STATS_KEY = "0123456789abcdef0123"); e <- load_usage(tempfile(), "pw"); expect_true(e$STATS_ON)
  expect_true(e$secret_equal("0123456789abcdef0123", e$STATS_KEY))
  expect_false(e$secret_equal("0123456789abcdef012", e$STATS_KEY)); expect_false(e$secret_equal("0123456789abcdef0124", e$STATS_KEY))
  expect_false(e$secret_equal(NULL, e$STATS_KEY)); expect_false(e$secret_equal("", e$STATS_KEY)); expect_false(e$secret_equal("", ""))
  expect_false(e$secret_equal(NA, e$STATS_KEY))
})

test_that("a session excluded before the browser reports in is never written", {
  db <- tempfile(fileext = ".sqlite"); e <- load_usage(db, "pw")
  t1 <- e$usage_tracker(fake_session("eeeeeeeeeeeeeeeeeeeeeeee")); t1$exclude(); t1$start(list(ref = "", w = 1440)); t1$log("gene", "KRT17")
  expect_equal(nrow(e$usage_read()$visits), 0L)
  e$usage_close(); unlink(paste0(db, c("", "-wal", "-shm")))
})
