# ---- private stats page: open the app with ?stats=<PDAC_STATS_KEY> and enter PDAC_STATS_PASSWORD --------------
# Without the exact key the app behaves as if the page did not exist (no prompt, no message). With it, a password prompt;
# on success a "Usage" page is added to the navigation of this browser tab only, and the visit is removed from the counts.
STATS_PERIODS <- c("Last 7 days" = "7", "Last 30 days" = "30", "Last 90 days" = "90", "All time" = "all")
STATS_MAX_TRIES <- 5L

stats_tile <- function(value, label, sub = NULL) div(class = "stat-tile", div(class = "stat-value", value), div(class = "stat-label", label),
                                                    if (!is.null(sub)) div(class = "stat-sub", sub))

stats_page_ui <- function(host) nav_panel("Usage", value = "usage",
  div(class = "settings-strip", div(class = "strip-label", "Usage of this app (private page)"),
      div(class = "strip-controls",
          radioButtons("stats_period", NULL, choices = STATS_PERIODS, selected = "30", inline = TRUE),
          actionButton("stats_refresh", "Refresh", class = "btn-sm btn-outline-secondary"),
          uiOutput("stats_nocount_ui", inline = TRUE))),
  uiOutput("stats_tiles"),
  plotlyOutput("p_stats_daily", height = "320px"),
  layout_columns(col_widths = breakpoints(sm = 12, lg = c(4, 4, 4)),
    card(card_header("Where visitors came from"), DTOutput("t_stats_ref", fill = FALSE),
         note("Host name of the referring page as the browser reports it. '(direct)': typed, bookmarked, opened from an app or e-mail, or hidden by the browser. ",
              "Tag your own links with ?ref=..., e.g. ", code(paste0("https://", host, "/?ref=linkedin")), ", to see them under 'Link tag'.")),
    card(card_header("Pages and gene tabs"), DTOutput("t_stats_pages", fill = FALSE), DTOutput("t_stats_tabs", fill = FALSE),
         note("Share = proportion of visits that opened the page or tab at least once (the start page and the Summary tab are open by default and not counted).")),
    card(card_header("Link tags and screens"), DTOutput("t_stats_tags", fill = FALSE), DTOutput("t_stats_dev", fill = FALSE))),
  layout_columns(col_widths = breakpoints(sm = 12, lg = c(4, 4, 4)),
    card(card_header("Genes looked up"), DTOutput("t_stats_genes", fill = FALSE),
         note("Counted once per visit; the random start gene is not counted.")),
    card(card_header("Downloads"), DTOutput("t_stats_dl", fill = FALSE)),
    card(card_header("Gene sets and gene lists"), DTOutput("t_stats_sets", fill = FALSE), DTOutput("t_stats_compare", fill = FALSE),
         note("Gene sets chosen on the survival page (the pre-filled example is not counted) and lists run in Compare genes."))),
  note("What is stored: a random visit id, start and end time, the referring host name, a screen-size class, an optional ?ref= tag, and the pages, tabs, genes and downloads used. ",
       "No cookies, IP addresses or browser details. Only browsers that run the page's script open a live session, so most bots are not counted; visits shorter than a few seconds still are."))

stats_server <- function(input, output, session, tr) {
  tries <- 0L; unlocked <- reactiveVal(FALSE)
  observeEvent(session$clientData$url_search, {
    q <- parseQueryString(session$clientData$url_search)
    if (!secret_equal(q$stats, STATS_KEY)) return()              # wrong or missing key: normal site
    tr$exclude()                                                  # the admin's visit is not counted, even before the password
    showModal(modalDialog(title = "Usage statistics", easyClose = TRUE,
      passwordInput("stats_pw", "Password", width = "100%"), uiOutput("stats_pw_msg"),
      footer = tagList(modalButton("Cancel"), actionButton("stats_go", "Open", class = "btn-primary"))))
  }, once = TRUE)
  observeEvent(input$stats_go, {
    if (unlocked() || !secret_equal(parseQueryString(session$clientData$url_search)$stats, STATS_KEY)) return()
    if (tries >= STATS_MAX_TRIES) { output$stats_pw_msg <- renderUI(note("Too many attempts; reload the page to try again.")); return() }
    tries <<- tries + 1L
    if (!secret_equal(input$stats_pw, STATS_PASSWORD)) {
      Sys.sleep(1); output$stats_pw_msg <- renderUI(note("Wrong password.")); return() }
    removeModal(); unlocked(TRUE)
    nav_insert("nav", stats_page_ui(session$clientData$url_hostname), target = "about", position = "after", select = TRUE, session = session)
  })

  nocount <- reactiveVal(NA)
  observeEvent(input$usage_client, nocount(isTRUE(input$usage_client$nocount)), once = TRUE)
  output$stats_nocount_ui <- renderUI({ req(unlocked())
    on <- isTRUE(nocount())
    span(class = "stats-nocount", actionLink("stats_nocount", if (on) "This browser is not counted (count it again)" else "Don't count my visits from this browser")) })
  observeEvent(input$stats_nocount, { on <- !isTRUE(nocount()); session$sendCustomMessage("usage_nocount", on); nocount(on) })

  u <- reactive({ req(unlocked()); input$stats_refresh
    p <- input$stats_period; days <- if (is.null(p) || p == "all") NULL else as.integer(p)
    x <- usage_read(days); if (is.null(x)) return(NULL); x$days <- days; x })
  tabs <- reactive({ x <- u(); req(x); usage_tables(x) })
  output$stats_tiles <- renderUI({ x <- u(); req(x); k <- usage_kpis(x)
    pct <- function(v) if (is.na(v)) "-" else sprintf("%.0f%%", 100 * v)
    mins <- function(v) if (is.na(v)) "-" else if (v < 1) sprintf("%.0f s", 60 * v) else sprintf("%.1f min", v)
    div(class = "stat-row",
        stats_tile(format(k$visits, big.mark = ","), "Visits", paste("on", k$days, if (k$days == 1) "day" else "days")),
        stats_tile(mins(k$median_min), "Median visit length", "first to last activity"),
        stats_tile(pct(k$engaged), "Visits that used the data", "looked up a gene, compared genes, built a gene set or downloaded"),
        stats_tile(format(k$genes, big.mark = ","), "Different genes looked up"),
        stats_tile(format(k$downloads, big.mark = ","), "Downloads")) })
  output$p_stats_daily <- renderPlotly({ x <- u(); req(x); to_plotly(plot_usage_daily(usage_daily(x, x$days)), height = 320, wrap = 60) })
  dt <- function(df, n = 10) datatable(df, rownames = FALSE, style = "default", options = list(dom = if (nrow(df) > n) "tp" else "t", pageLength = n, ordering = FALSE), class = "compact")
  output$t_stats_ref <- renderDT(dt(tabs()$referrers), server = FALSE)
  output$t_stats_tags <- renderDT(dt(tabs()$tags, 5), server = FALSE)
  output$t_stats_dev <- renderDT(dt(tabs()$devices), server = FALSE)
  output$t_stats_pages <- renderDT(dt(tabs()$pages), server = FALSE)
  output$t_stats_tabs <- renderDT(dt(tabs()$tabs, 12), server = FALSE)
  output$t_stats_genes <- renderDT(dt(tabs()$genes, 15), server = FALSE)
  output$t_stats_dl <- renderDT(dt(tabs()$downloads, 15), server = FALSE)
  output$t_stats_sets <- renderDT(dt(tabs()$sets, 10), server = FALSE)
  output$t_stats_compare <- renderDT(dt(tabs()$compare, 10), server = FALSE)
}
