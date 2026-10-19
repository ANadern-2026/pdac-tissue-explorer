# PDAC Tissue Explorer - Shiny app
# Explore where any gene is expressed in human pancreatic ductal adenocarcinoma (PDAC) tissue, using public
# proteomic, bulk transcriptomic and single-cell data. Data files are built by the scripts in ../data-raw/.
suppressPackageStartupMessages({
  library(shiny); library(bslib); library(ggplot2); library(plotly); library(DT); library(Matrix)
})
# used via pkg:: or indirectly (includeMarkdown); listed here so deployment tools detect them
if (FALSE) { library(markdown); library(patchwork); library(scales) }
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f, local = TRUE)
D <- load_app_data()

EXAMPLE_GENES <- intersect(c("MSLN", "KRT17", "FAP", "GATA6", "CEACAM5", "COL11A1", "CD274", "KRT19"), D$genes$symbol)
EXAMPLE_LISTS <- list(
  "Fibroblast markers" = "COL1A1, DCN, LUM, FAP, ACTA2, PDGFRB, CXCL12",
  "Classical vs basal-like tumour cells" = "GATA6, TFF1, LGALS4, CLDN18, KRT17, KRT5, S100A2, VGLL1",
  "Immune checkpoints" = "CD274, PDCD1LG2, PDCD1, CTLA4, LAG3, HAVCR2, TIGIT")
DEFAULT_GENE <- if (isTRUE(D$has_sc)) { g <- intersect(EXAMPLE_GENES, D$genes$symbol[!is.na(D$genes$sym_sc)]); if (length(g)) g[1] else D$genes$symbol[!is.na(D$genes$sym_sc)][1] } else EXAMPLE_GENES[1]

theme_app <- bs_theme(version = 5, bg = "#fcfcfb", fg = "#0b0b0b", primary = "#2a78d6", secondary = "#52514e",
                      base_font = font_collection("Inter", "system-ui", "-apple-system", "Segoe UI", "Roboto", "Helvetica Neue", "Arial", "sans-serif"),
                      "navbar-bg" = "#ffffff", "card-border-color" = "#e9e8e4")

to_plotly <- function(p, height = NULL) {
  ttl <- p$labels$title; sub <- p$labels$subtitle; faceted <- !inherits(p$facet, "FacetNull")
  p$labels$title <- NULL; p$labels$subtitle <- NULL; p$labels$caption <- NULL
  out <- ggplotly(p, tooltip = "text", height = height)
  out <- plotly_build(out)
  for (i in seq_along(out$x$data)) {
    if (identical(out$x$data[[i]]$type, "box")) { out$x$data[[i]]$marker$opacity <- 0; out$x$data[[i]]$hoverinfo <- "none" }   # hidden outliers
    if (!is.null(out$x$data[[i]]$marker$colorbar)) out$x$data[[i]]$marker$colorbar <- modifyList(out$x$data[[i]]$marker$colorbar,
      list(len = 0.85, lenmode = "fraction", y = 0.5, yanchor = "middle", thickness = 14, title = list(text = "Relative<br>expression")))
  }
  top <- (if (!is.null(sub)) 75 else 45) + (if (faceted) 30 else 0)
  out <- layout(out, showlegend = FALSE,
                title = list(text = paste0("<b>", ttl, "</b>", if (!is.null(sub)) paste0("<br><span style='font-size:12px;color:#52514e'>", sub, "</span>") else ""),
                             x = 0, xanchor = "left", font = list(size = 15)),
                margin = list(t = top), hoverlabel = list(bgcolor = "white", font = list(color = "#0b0b0b")))
  config(out, displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d", "hoverCompareCartesian", "toggleSpikelines"))
}

fig_downloads <- function(id) div(class = "dl-row", downloadButton(paste0(id, "_pdf"), "Figure (PDF)", class = "btn-sm btn-outline-secondary"),
                                  downloadButton(paste0(id, "_png"), "Figure (PNG)", class = "btn-sm btn-outline-secondary"),
                                  downloadButton(paste0(id, "_csv"), "Data (CSV)", class = "btn-sm btn-outline-secondary"))
note <- function(...) p(class = "note", ...)

# ---- UI ------------------------------------------------------------------------------------------------------
ui <- page_navbar(
  title = span(APP_TITLE, span(class = "badge-version", paste0("v", APP_VERSION))), theme = theme_app, id = "nav", fillable = FALSE,
  header = tagList(tags$head(includeCSS("www/style.css"), tags$title(APP_TITLE)),
                   if (length(MOVED_TO)) div(class = "moved-banner", strong("This app has moved: "),
                                             a(href = MOVED_TO[1], target = "_top", sub("/$", "", sub("^https?://", "", MOVED_TO[1]))),
                                             if (length(MOVED_TO) > 1) span(" ", MOVED_TO[2]))),
  nav_panel("Explore a gene", value = "gene",
    layout_sidebar(
      sidebar = sidebar(width = 300,
        selectizeInput("gene", "Gene", choices = NULL, options = list(placeholder = "Symbol, alias or name", maxOptions = 50)),
        div(class = "examples", span("Try: "), lapply(EXAMPLE_GENES[1:min(5, length(EXAMPLE_GENES))], function(g) actionLink(paste0("ex_", g), g))),
        hr(),
        selectInput("context", "Tissue context for the cell-type view", choices = STATE_LEVELS, selected = "Primary tumour"),
        checkboxGroupInput("compartments", "Cell types in the tissue-state view",
                           choices = CELLTYPE_ORDER, selected = c("Fibroblasts", "Ductal: tumour subclusters", "Ductal: other", "Endothelial cells", "Myeloid cells", "T and NK cells")),
        note("Each dot in the single-cell plots is one tissue sample. A sample contributes to a cell type only if it has at least ", MIN_CELLS, " cells of that type.")),
      uiOutput("gene_header"),
      navset_card_underline(id = "gene_tabs",
        nav_panel("Summary", uiOutput("summary_ui")),
        nav_panel("Tumour vs normal",
          layout_columns(col_widths = c(6, 6), plotlyOutput("p_prot", height = "460px"), plotlyOutput("p_rna", height = "460px")),
          fig_downloads("dl_cptac"), hr(), plotlyOutput("p_bulk", height = "440px"), fig_downloads("dl_bulk"),
          note("CPTAC: Cao et al., Cell 2021 (135 PDAC tumours; 64 protein and 19 RNA tumour/adjacent-normal pairs). TCGA/GTEx via recount3. ",
               "Adjacent normal tissue from pancreatic-cancer surgery is fibrotic and inflamed; it is not healthy pancreas.")),
        nav_panel("Cell types",
          plotlyOutput("p_ct", height = "520px"),
          layout_columns(col_widths = c(6, 6), plotlyOutput("p_share", height = "420px"), plotlyOutput("p_frac", height = "420px")),
          fig_downloads("dl_ct"),
          note("Average expression = log1p of the mean normalised expression (CP10k) of all cells of that type in a sample. ",
               "'Ductal: tumour subclusters' are ductal cells assigned to the atlas's malignant subclusters; 'Ductal: other' are remaining ductal cells, mostly non-malignant duct.")),
        nav_panel("Tissue states", plotlyOutput("p_states", height = "640px"), fig_downloads("dl_states"),
          note("Healthy donor pancreas, adjacent normal tissue and treatment-naive primary tumours and liver metastases; treated primary tumours (chemotherapy or chemoradiation) are shown separately. ",
               "Tissue states come from different studies, so study and state are partly confounded.")),
        nav_panel("Subtypes",
          plotlyOutput("p_fibro", height = "460px"), fig_downloads("dl_fibro"), hr(),
          plotlyOutput("p_ductal", height = "440px"), fig_downloads("dl_ductal"),
          note("Fibroblast and tumour-cell subclusters as published with the atlas (Loveless et al., Clin Cancer Res 2025).")),
        nav_panel("Atlas map",
          layout_columns(col_widths = c(6, 6), plotlyOutput("p_umap", height = "560px"), plotlyOutput("p_umap_ct", height = "560px")),
          fig_downloads("dl_umap"),
          note("A stratified subsample of atlas cells, positioned by the published atlas UMAP. Grey = not detected.")),
        nav_panel("Data",
          note("Per-sample single-cell values for the selected gene (all tissue states). Use the buttons to copy or download."),
          DTOutput("t_sc"), hr(), note("CPTAC matched tumour/normal values."), DTOutput("t_cptac"))
      ))),
  nav_panel("Compare genes", value = "multi",
    layout_sidebar(
      sidebar = sidebar(width = 320,
        textAreaInput("genelist", "Genes (comma, space or line separated; up to 40)", value = EXAMPLE_LISTS[[1]], rows = 6),
        div(class = "examples", span("Examples: "), lapply(seq_along(EXAMPLE_LISTS), function(i) actionLink(paste0("exl_", i), names(EXAMPLE_LISTS)[i]))),
        actionButton("go_multi", "Compare", class = "btn-primary w-100"), uiOutput("multi_msg")),
      card(card_header("Cell types"), plotlyOutput("p_mdot", height = "auto"), fig_downloads("dl_mdot")),
      card(card_header("Tumour vs normal"), plotlyOutput("p_mforest", height = "auto"), fig_downloads("dl_mforest")))),
  nav_panel("About", value = "about", includeMarkdown("www/about.md")),
  nav_spacer(),
  nav_item(tags$a(href = REPO_URL, target = "_blank", "Source code"))
)

# ---- server --------------------------------------------------------------------------------------------------
server <- function(input, output, session) {
  if (length(MOVED_TO)) showModal(modalDialog(
    title = paste(APP_TITLE, "has moved"), easyClose = TRUE,
    p("The app now runs at its permanent address:"),
    p(a(href = MOVED_TO[1], target = "_top", style = "font-size: 1.15rem; font-weight: 600;", sub("/$", "", sub("^https?://", "", MOVED_TO[1])))),
    if (length(MOVED_TO) > 1) p(MOVED_TO[2]),
    p(class = "note", "Same content and data; please update your bookmarks."),
    footer = tagList(modalButton("Continue here"), a(class = "btn btn-primary", href = MOVED_TO[1], target = "_top", "Go to the new address"))))
  gene_df <- data.frame(value = D$genes$symbol, label = D$genes$symbol, name = ifelse(is.na(D$genes$name), "", D$genes$name),
                        aliases = D$genes$aliases, stringsAsFactors = FALSE)
  updateSelectizeInput(session, "gene", choices = gene_df, selected = DEFAULT_GENE, server = TRUE,
                       options = list(valueField = "value", labelField = "label", searchField = c("label", "aliases", "name"),
                                      score = I("function(search) { var base = this.getScoreFunction(search); var q = String(typeof search === 'string' ? search : (search.query || '')).toUpperCase().trim();
                                                   return function(item) { var l = item.label.toUpperCase(); var al = (item.aliases || '').toUpperCase().split(/,\\s*/);
                                                     if (l === q) return 1000; if (al.indexOf(q) >= 0) return 500; if (l.indexOf(q) === 0) return 100 + base(item); return base(item); }; }"),
                                      onFocus = I("function() { this._prev = this.getValue(); this.clear(true); }"),
                                      onBlur = I("function() { if (!this.getValue() && this._prev) { this.setValue(this._prev, true); } }"),
                                      render = I("{ option: function(item, escape) { return '<div><b>' + escape(item.label) + '</b> <span class=\"opt-name\">' + escape(item.name) + '</span>' + (item.aliases ? '<div class=\"opt-alias\">' + escape(item.aliases) + '</div>' : '') + '</div>'; } }")))
  for (g in EXAMPLE_GENES) local({ gg <- g; observeEvent(input[[paste0("ex_", gg)]], updateSelectizeInput(session, "gene", choices = gene_df, selected = gg, server = TRUE)) })
  for (i in seq_along(EXAMPLE_LISTS)) local({ ii <- i; observeEvent(input[[paste0("exl_", ii)]], updateTextAreaInput(session, "genelist", value = EXAMPLE_LISTS[[ii]])) })

  gene <- reactive({ req(input$gene, input$gene %in% D$genes$symbol); input$gene })
  pr <- reactive(cptac_pairs(D, gene(), "protein")); rn <- reactive(cptac_pairs(D, gene(), "rna"))
  bk <- reactive(bulk_values(D, gene())); sv <- reactive(sc_values(D, gene()))
  cs <- reactive({ v <- sv(); if (is.null(v)) NULL else celltype_summary(v, input$context) })
  cs_naive <- reactive({ v <- sv(); if (is.null(v)) NULL else celltype_summary(v, "Primary tumour") })

  figs <- list(
    prot = reactive(plot_cptac_paired(pr(), gene(), "protein")), rna = reactive(plot_cptac_paired(rn(), gene(), "RNA")),
    bulk = reactive(plot_bulk(bk(), gene())), ct = reactive(plot_celltypes(cs(), gene(), input$context)),
    share = reactive(plot_share(cs(), gene())), frac = reactive(plot_frac(cs(), gene())),
    states = reactive(plot_states(sv(), gene(), input$compartments)), fibro = reactive(plot_fibro(sv(), gene())),
    ductal = reactive(plot_ductal(sv(), gene())), umap = reactive(plot_umap_static(D$umap, gene(), gene_row(D, gene())$sym_sc)))

  output$gene_header <- renderUI({
    gr <- gene_row(D, gene())
    chip <- function(lbl, ok) span(class = if (ok) "chip chip-on" else "chip chip-off", lbl)
    div(class = "gene-header",
        h2(gr$symbol, span(class = "gene-name", if (!is.na(gr$name)) gr$name else "")),
        if (nzchar(gr$aliases)) p(class = "aliases", "Also known as: ", gr$aliases),
        div(chip("Protein (CPTAC)", !is.na(gr$sym_cptac_protein)), chip("RNA (CPTAC)", !is.na(gr$sym_cptac_rna)),
            chip("RNA (TCGA/GTEx)", !is.na(gr$sym_tcga_gtex)), chip("Single-cell atlas", !is.na(gr$sym_sc) && isTRUE(D$has_sc))))
  })

  output$summary_ui <- renderUI({
    g <- gene(); sp <- paired_stats(pr()); sr <- paired_stats(rn())
    items <- list(txt_cptac(pr(), sp, g, "protein"), txt_cptac(rn(), sr, g, "rna"), txt_bulk(bk(), g), txt_sc(cs_naive(), g), txt_state(sv(), g))
    items <- Filter(Negate(is.null), items)
    tagList(tags$ul(class = "summary-list", lapply(items, tags$li)),
            note("These sentences describe the data; they do not establish mechanism or causality. p-values are unadjusted and refer to the single comparison stated. ",
                 "See the other tabs for the full distributions and the About page for methods and limitations."))
  })

  output$p_prot <- renderPlotly(to_plotly(figs$prot())); output$p_rna <- renderPlotly(to_plotly(figs$rna()))
  output$p_bulk <- renderPlotly(to_plotly(figs$bulk())); output$p_ct <- renderPlotly(to_plotly(figs$ct()))
  output$p_share <- renderPlotly(to_plotly(figs$share())); output$p_frac <- renderPlotly(to_plotly(figs$frac()))
  output$p_states <- renderPlotly(to_plotly(figs$states())); output$p_fibro <- renderPlotly(to_plotly(figs$fibro()))
  output$p_ductal <- renderPlotly(to_plotly(figs$ductal()))

  output$p_umap <- renderPlotly({
    u <- D$umap; validate(need(!is.null(u), "Atlas map not available.")); col <- gene_row(D, gene())$sym_sc
    validate(need(!is.na(col), paste(gene(), "is not in the single-cell data.")))
    val <- as.numeric(u$expr[, col]); o <- order(val)
    plot_ly(x = u$coords$UMAP_1[o], y = u$coords$UMAP_2[o], type = "scattergl", mode = "markers",
            marker = list(size = 3, color = val[o], colorscale = list(c(0, INK[["zero"]]), c(0.001, SEQ_RAMP[1]), c(1, SEQ_RAMP[6])), showscale = TRUE,
                          colorbar = list(title = "log1p\nCP10k", thickness = 10)),
            text = paste0(u$meta$compartment[o], "<br>", gene(), ": ", fmt_num(val[o])), hoverinfo = "text") |>
      layout(title = list(text = paste0("<b>", gene(), " expression</b>"), x = 0), xaxis = list(visible = FALSE), yaxis = list(visible = FALSE, scaleanchor = "x")) |>
      config(displaylogo = FALSE)
  })
  output$p_umap_ct <- renderPlotly({
    u <- D$umap; validate(need(!is.null(u), ""))
    lab <- as.character(u$meta$compartment); lab[grepl("^Ductal", lab)] <- "Ductal cells"
    cen <- stats::aggregate(cbind(UMAP_1, UMAP_2) ~ compartment, data.frame(u$coords, compartment = lab), stats::median)
    plot_ly() |>
      add_trace(x = u$coords$UMAP_1, y = u$coords$UMAP_2, type = "scattergl", mode = "markers", marker = list(size = 3, color = "#d9d8d4"),
                text = paste0(u$meta$compartment, "<br>", u$meta$disease_state), hoverinfo = "text", showlegend = FALSE) |>
      add_annotations(x = cen$UMAP_1, y = cen$UMAP_2, text = as.character(cen$compartment), showarrow = FALSE, font = list(size = 12, color = "#0b0b0b"),
                      bgcolor = "rgba(255,255,255,0.75)") |>
      layout(title = list(text = "<b>Cell types</b>", x = 0), xaxis = list(visible = FALSE), yaxis = list(visible = FALSE, scaleanchor = "x")) |>
      config(displaylogo = FALSE)
  })

  sc_table <- reactive({ v <- sv(); req(v)
    a <- sc_aggregate(v[!is.na(v$state), ], c("sample", "study", "state", "compartment"))
    a$study <- paste0(a$study, " (", STUDY_CITATION[a$study], ")"); a[order(a$compartment, a$state, a$sample), ] })
  cptac_table <- reactive({ rbind(if (!is.null(pr())) data.frame(assay = "protein", pr()$pairs), if (!is.null(rn())) data.frame(assay = "RNA", rn()$pairs)) })
  output$t_sc <- renderDT(datatable(transform(sc_table(), mean_cp10k = round(mean_cp10k, 3), avg_expr = round(avg_expr, 3), frac = round(frac, 4)),
                                    rownames = FALSE, extensions = "Buttons", options = list(pageLength = 10, dom = "Bfrtip", buttons = c("copy", "csv"))))
  output$t_cptac <- renderDT(datatable(cptac_table(), rownames = FALSE, extensions = "Buttons", options = list(pageLength = 10, dom = "Bfrtip", buttons = c("copy", "csv"))) |>
                               formatRound(c("normal", "tumour", "diff"), 3))

  # downloads: figure(s) + underlying data
  fname <- function(stub, ext) function() paste0(gene(), "_", stub, "_", format(Sys.Date()), ".", ext)
  save_figs <- function(file, plots, w, h, dev) {
    plots <- lapply(plots, for_export, width_in = w)
    if (length(plots) == 1) { ggsave(file, plots[[1]], width = w, height = h, device = dev, bg = "white", dpi = 200); return(invisible()) }
    if (dev == "pdf") { grDevices::pdf(file, width = w, height = h); for (p in plots) print(p); grDevices::dev.off() }
    else ggsave(file, patchwork::wrap_plots(plots, ncol = length(plots)), width = w * length(plots), height = h, device = dev, bg = "white", dpi = 200)
  }
  dl <- list(dl_cptac = list(stub = "CPTAC_tumour_vs_normal", plots = function() list(figs$prot(), figs$rna()), data = cptac_table, w = 6, h = 5),
             dl_bulk = list(stub = "TCGA_GTEx", plots = function() list(figs$bulk()), data = bk, w = 7, h = 5),
             dl_ct = list(stub = "cell_types", plots = function() list(figs$ct(), figs$share(), figs$frac()), data = function() cs()$summary, w = 8, h = 6),
             dl_states = list(stub = "tissue_states", plots = function() list(figs$states()), data = sc_table, w = 11, h = 8),
             dl_fibro = list(stub = "fibroblast_subtypes", plots = function() list(figs$fibro()), data = sc_table, w = 9, h = 5),
             dl_ductal = list(stub = "tumour_cell_states", plots = function() list(figs$ductal()), data = sc_table, w = 9, h = 5),
             dl_umap = list(stub = "atlas_map", plots = function() list(figs$umap()), data = function() data.frame(note = "Atlas map values are a subsample; see the Data tab for per-sample values."), w = 7, h = 6))
  for (k in names(dl)) local({ kk <- k; s <- dl[[kk]]
    output[[paste0(kk, "_pdf")]] <- downloadHandler(fname(s$stub, "pdf"), function(file) save_figs(file, s$plots(), s$w, s$h, "pdf"))
    output[[paste0(kk, "_png")]] <- downloadHandler(fname(s$stub, "png"), function(file) save_figs(file, s$plots(), s$w, s$h, "png"))
    output[[paste0(kk, "_csv")]] <- downloadHandler(fname(s$stub, "csv"), function(file) utils::write.csv(s$data(), file, row.names = FALSE))
  })

  # ---- compare genes ----
  multi <- eventReactive(input$go_multi, { pl <- parse_gene_list(input$genelist, D); pl$genes <- utils::head(pl$genes, 40); pl }, ignoreNULL = FALSE)
  output$multi_msg <- renderUI({ m <- multi(); if (length(m$missing)) note("Not found: ", paste(m$missing, collapse = ", ")) })
  mdot <- reactive(multi_dot_table(D, multi()$genes)); mfor <- reactive(multi_forest_table(D, multi()$genes))
  mdot_fig <- reactive(plot_multi_dot(mdot())); mfor_fig <- reactive(plot_multi_forest(mfor()))
  hgt <- reactive(160 + 26 * max(3, length(multi()$genes)))
  output$p_mdot <- renderPlotly(to_plotly(mdot_fig(), height = hgt() + 160))
  output$p_mforest <- renderPlotly(to_plotly(mfor_fig(), height = hgt()))
  mfn <- function(stub, ext) function() paste0("gene_set_", stub, "_", format(Sys.Date()), ".", ext)
  output$dl_mdot_pdf <- downloadHandler(mfn("cell_types", "pdf"), function(f) ggsave(f, for_export(mdot_fig(), 10), width = 10, height = 2 + 0.3 * length(multi()$genes)))
  output$dl_mdot_png <- downloadHandler(mfn("cell_types", "png"), function(f) ggsave(f, for_export(mdot_fig(), 10), width = 10, height = 2 + 0.3 * length(multi()$genes), dpi = 200, bg = "white"))
  output$dl_mdot_csv <- downloadHandler(mfn("cell_types", "csv"), function(f) utils::write.csv(mdot(), f, row.names = FALSE))
  output$dl_mforest_pdf <- downloadHandler(mfn("tumour_vs_normal", "pdf"), function(f) ggsave(f, for_export(mfor_fig(), 11), width = 11, height = 2 + 0.3 * length(multi()$genes)))
  output$dl_mforest_png <- downloadHandler(mfn("tumour_vs_normal", "png"), function(f) ggsave(f, for_export(mfor_fig(), 11), width = 11, height = 2 + 0.3 * length(multi()$genes), dpi = 200, bg = "white"))
  output$dl_mforest_csv <- downloadHandler(mfn("tumour_vs_normal", "csv"), function(f) utils::write.csv(mfor(), f, row.names = FALSE))
}

shinyApp(ui, server)
