# PDAC Tissue Explorer - Shiny app
# Explore where any gene is expressed in human pancreatic ductal adenocarcinoma (PDAC) tissue, using public
# proteomic, bulk transcriptomic and single-cell data. Data files are built by the scripts in ../data-raw/.
suppressPackageStartupMessages({
  library(shiny); library(bslib); library(ggplot2); library(plotly); library(DT)
})
# used via pkg:: or indirectly (includeMarkdown); listed here so deployment tools detect them
if (FALSE) { library(markdown); library(patchwork); library(scales); library(survival) }
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f, local = TRUE)
D <- load_app_data()            # gene index + small tables in memory; per-gene values are read from the on-disk store
onStop(function() close_app_data(D))

EXAMPLE_GENES <- intersect(c("MSLN", "KRT17", "FAP", "GATA6", "CEACAM5", "COL11A1", "CD274", "KRT19"), D$genes$symbol)
EXAMPLE_LISTS <- list(
  "Fibroblast markers" = "COL1A1, DCN, LUM, FAP, ACTA2, PDGFRB, CXCL12",
  "Classical vs basal-like tumour cells" = "GATA6, TFF1, LGALS4, CLDN18, KRT17, KRT5, S100A2, VGLL1",
  "Immune checkpoints" = "CD274, PDCD1LG2, PDCD1, CTLA4, LAG3, HAVCR2, TIGIT")
# start-up gene: a random, well-covered gene per session (measured in CPTAC protein and RNA, TCGA/GTEx and the single-cell atlas,
# and expressed widely enough in at least one cell type for co-expression), so first views are rarely empty
START_POOL <- with(D$genes, symbol[!is.na(sym_cptac_protein) & !is.na(sym_cptac_rna) & !is.na(sym_tcga_gtex) & !is.na(sym_sc)])
if (isTRUE(D$has_coexp)) START_POOL <- START_POOL[gene_row(D, START_POOL)$sym_sc %in% store_ids(D$coexp, "partners")]
if (!length(START_POOL)) START_POOL <- EXAMPLE_GENES[1]
random_start_gene <- function() { set.seed(NULL); sample(START_POOL, 1L) }   # re-seeded: bootstrap CIs set fixed seeds

SURV_CHOICES <- c("TCGA-PAAD RNA: overall survival" = "TCGA-PAAD RNA|os", "TCGA-PAAD RNA: progression-free interval" = "TCGA-PAAD RNA|pfi",
                  "CPTAC protein: overall survival" = "CPTAC protein|os", "CPTAC RNA: overall survival" = "CPTAC RNA|os")
MOFFITT_GENES <- if (isTRUE(D$has_clinical)) unique(unlist(D$clin$moffitt)) else character(0)

theme_app <- bs_theme(version = 5, bg = "#fcfcfb", fg = "#0b0b0b", primary = "#2a78d6", secondary = "#52514e",
                      base_font = font_collection("Inter", "system-ui", "-apple-system", "Segoe UI", "Roboto", "Helvetica Neue", "Arial", "sans-serif"),
                      "navbar-bg" = "#ffffff", "card-border-color" = "#e9e8e4")

to_plotly <- function(p, height = NULL, wrap = NULL) {      # wrap: max characters per subtitle line (plots in narrow columns)
  ttl <- p$labels$title; sub <- p$labels$subtitle
  if (!is.null(sub) && !is.null(wrap)) sub <- paste(unlist(lapply(strsplit(sub, "\n", fixed = TRUE)[[1]], strwrap, width = wrap)), collapse = "\n"); faceted <- !inherits(p$facet, "FacetNull")
  p$labels$title <- NULL; p$labels$subtitle <- NULL; p$labels$caption <- NULL
  out <- ggplotly(p, tooltip = "text", height = height)
  out <- plotly_build(out)
  for (i in seq_along(out$x$data)) {
    if (identical(out$x$data[[i]]$type, "box")) { out$x$data[[i]]$marker$opacity <- 0; out$x$data[[i]]$hoverinfo <- "none" }   # hidden outliers
    if (!is.null(out$x$data[[i]]$marker$colorbar)) out$x$data[[i]]$marker$colorbar <- modifyList(out$x$data[[i]]$marker$colorbar,
      list(len = 0.85, lenmode = "fraction", y = 0.5, yanchor = "middle", thickness = 14, title = list(text = "Relative<br>expression")))
  }
  # multi-line subtitles ("\n") and two-line facet labels get more room; single-line layouts are unchanged
  sub_lines <- if (is.null(sub)) 0L else length(strsplit(sub, "\n", fixed = TRUE)[[1]])
  strip_lines <- if (!faceted) 0L else { v <- names(p$facet$params$facets); if (!length(v) || !v[1] %in% names(p$data)) 1L else
    max(1L, vapply(strsplit(as.character(unique(p$data[[v[1]]])), "\n", fixed = TRUE), length, integer(1))) }
  ttl_lines <- 1L
  if (!is.null(ttl) && !is.null(wrap)) { tw <- strwrap(ttl, width = round(0.8 * wrap)); ttl_lines <- length(tw); ttl <- paste(tw, collapse = "<br>") }   # bold 15 px: fewer characters per line
  top <- 52 + 22 * (ttl_lines - 1L) + 20 * sub_lines + 28 * strip_lines       # title block is anchored at the top (see layout below)
  if (!is.null(sub)) sub <- gsub("\n", "<br>", sub, fixed = TRUE)
  out <- layout(out, showlegend = FALSE,
                title = list(text = paste0("<b>", ttl, "</b>", if (!is.null(sub)) paste0("<br><span style='font-size:12px;color:#52514e'>", sub, "</span>") else ""),
                             x = 0, xanchor = "left", y = 1, yref = "container", yanchor = "top", pad = list(t = 14), font = list(size = 15)),   # anchored at the top: long wrapped titles never reach the plot
                margin = list(t = top), hoverlabel = list(bgcolor = "white", font = list(color = "#0b0b0b")))
  config(out, displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d", "hoverCompareCartesian", "toggleSpikelines"))
}

fig_downloads <- function(id) div(class = "dl-row", downloadButton(paste0(id, "_pdf"), "Figure (PDF)", class = "btn-sm btn-outline-secondary"),
                                  downloadButton(paste0(id, "_png"), "Figure (PNG)", class = "btn-sm btn-outline-secondary"),
                                  downloadButton(paste0(id, "_csv"), "Data (CSV)", class = "btn-sm btn-outline-secondary"))
note <- function(...) p(class = "note", ...)
strip <- function(...) div(class = "settings-strip", div(class = "strip-label", "View settings"), div(class = "strip-controls", ...))
# Data-tab tables: DataTables' own style (style = "default") - the Ubuntu r-cran-dt package used on the server lacks the
# Bootstrap 5 styles that DT would otherwise pick for a bslib page, which left the buttons unstyled. CSV file name: Buttons >= 1.2 reads `filename`, the server's Buttons 1.0.3 reads `title`.
dt_buttons <- function(symbol, stub) { f <- paste0(symbol, "_", stub, "_", format(Sys.Date())); list("copy", list(extend = "csvHtml5", filename = f, title = f)) }

# ---- UI ------------------------------------------------------------------------------------------------------
ui <- page_navbar(
  title = span(APP_TITLE, span(class = "badge-version", paste0("v", APP_VERSION))), theme = theme_app, id = "nav", fillable = FALSE,
  header = tagList(tags$head(includeCSS("www/style.css"), tags$title(APP_TITLE)),
                   if (length(MOVED_TO)) div(class = "moved-banner", strong("This app has moved: "),
                                             a(href = MOVED_TO[1], target = "_top", sub("/$", "", sub("^https?://", "", MOVED_TO[1]))),
                                             if (length(MOVED_TO) > 1) span(" ", MOVED_TO[2]))),
  nav_panel("Explore a gene", value = "gene",
    div(class = "gene-band",
      div(class = "gene-search",
        selectizeInput("gene", "Gene", choices = NULL, width = "100%", options = list(placeholder = "Type a gene symbol, alias or name", maxOptions = 50)),
        div(class = "examples", span("Try: "), lapply(EXAMPLE_GENES[1:min(5, length(EXAMPLE_GENES))], function(g) actionLink(paste0("ex_", g), g)))),
      div(class = "gene-band-header", uiOutput("gene_header"))),
    navset_card_underline(id = "gene_tabs",
        nav_panel("Summary", uiOutput("summary_ui")),
        nav_panel("Tumour vs normal",
          layout_columns(col_widths = breakpoints(sm = 12, lg = c(12, 8), xl = c(8, 4)),
            div(layout_columns(col_widths = breakpoints(sm = 12, md = c(6, 6)), plotlyOutput("p_prot", height = "460px"), plotlyOutput("p_rna", height = "460px")),
                fig_downloads("dl_cptac")),
            div(plotlyOutput("p_bulk", height = "460px"), fig_downloads("dl_bulk"))),
          note("CPTAC: Cao et al., Cell 2021 (135 PDAC tumours; 64 protein and 19 RNA tumour/adjacent-normal pairs). TCGA/GTEx via recount3. ",
               "Adjacent normal tissue from pancreatic-cancer surgery is fibrotic and inflamed; it is not healthy pancreas.")),
        if (isTRUE(D$has_phospho)) nav_panel("Phosphosites",
          strip(radioButtons("ph_mode", "Values", inline = TRUE, choices = c("Site abundance" = "raw", "Relative to the protein level" = "rel"))),
          layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
            plotlyOutput("p_phsites", height = "auto"),
            conditionalPanel("output.ph_has", uiOutput("ph_site_ui"),
              layout_columns(col_widths = breakpoints(sm = 12, md = c(6, 6)), plotlyOutput("p_phsite", height = "440px"), plotlyOutput("p_phsite_rel", height = "440px")))),
          fig_downloads("dl_phospho"), conditionalPanel("output.ph_has", hr(), DTOutput("t_phospho", fill = FALSE)),
          note("CPTAC phosphoproteome (Cao et al., Cell 2021): phosphosite intensities in tumours and matched adjacent normal tissue of the same patients. ",
               "A site is compared only if quantified in both tissues of at least ", PHOSPHO_MIN_PAIRS, " patients (paired Wilcoxon test; false discovery rate across the sites of this gene). ",
               "'Relative to the protein level' subtracts the gene's protein level in the same sample, separating a change in phosphorylation from a change in the amount of protein. ",
               "Positions refer to the RefSeq protein given in the table; a change in phosphosite abundance does not by itself show a change in activity.")),
        if (isTRUE(D$has_clinical)) nav_panel("Clinical context",
          strip(radioButtons("surv_set", "Survival in", inline = TRUE, choices = SURV_CHOICES)),
          layout_columns(col_widths = breakpoints(sm = 12, lg = c(7, 5)), plotlyOutput("p_km", height = "500px"), plotlyOutput("p_hr", height = "500px")),
          fig_downloads("dl_surv"),
          note("Kaplan-Meier curves compare tumours above the median with those at or below it (log-rank test). Hazard ratios come from Cox models of the continuous value ",
               "(per standard deviation within the dataset), unadjusted and adjusted for age and stage. No cut-point was optimised. ",
               "TCGA-PAAD: overall survival and progression-free interval from the TCGA Pan-Cancer Clinical Data Resource (Liu et al., Cell 2018); CPTAC: overall survival from the study's follow-up."),
          hr(), layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
            div(plotlyOutput("p_stage", height = "500px"), fig_downloads("dl_stage")),
            div(plotlyOutput("p_moffitt", height = "500px"), fig_downloads("dl_moffitt"))),
          note("These associations are exploratory and observational. Survival after surgery depends on stage, resection margin, adjuvant therapy and tumour-cell content, ",
               "and with many genes tested, some will show p < 0.05 by chance. Bulk tumour values mix tumour and stromal cells: a gene expressed by stroma can appear ",
               "'classical' or 'basal-like' because the score partly tracks tumour-cell content. Use these views to generate hypotheses, not as evidence of prognostic value.")),
        nav_panel("Cell types",
          strip(radioButtons("context", "Tissue context", choices = STATE_LEVELS, selected = "Primary tumour", inline = TRUE)),
          layout_columns(col_widths = breakpoints(sm = 12, xl = c(6, 6)),
            plotlyOutput("p_ct", height = "580px"), div(plotlyOutput("p_share", height = "290px"), plotlyOutput("p_frac", height = "290px"))),
          fig_downloads("dl_ct"),
          note("Each dot is one tissue sample; a sample contributes to a cell type only if it has at least ", MIN_CELLS, " cells of that type. Average expression = log1p of the mean normalised expression (CP10k) of all cells of that type in a sample. ",
               "'Ductal: tumour subclusters' are ductal cells assigned to the atlas's malignant subclusters; 'Ductal: other' are remaining ductal cells, mostly non-malignant duct.")),
        nav_panel("Tissue states",
          strip(checkboxGroupInput("compartments", "Cell types shown", choices = CELLTYPE_ORDER, inline = TRUE,
                                   selected = c("Fibroblasts", "Ductal: tumour subclusters", "Ductal: other", "Endothelial cells", "Myeloid cells", "T and NK cells"))),
          plotlyOutput("p_states", height = "640px"), fig_downloads("dl_states"),
          note("One dot per tissue sample with at least ", MIN_CELLS, " cells of that type. Healthy donor pancreas, adjacent normal tissue and treatment-naive primary tumours and liver metastases; treated primary tumours (chemotherapy or chemoradiation) are shown separately. ",
               "Tissue states come from different studies, so study and state are partly confounded.")),
        nav_panel("Subtypes",
          layout_columns(col_widths = breakpoints(sm = 12, xl = c(7, 5)),
            div(plotlyOutput("p_fibro", height = "480px"), fig_downloads("dl_fibro")),
            div(plotlyOutput("p_ductal", height = "480px"), fig_downloads("dl_ductal"))),
          note("Fibroblast and tumour-cell subclusters as published with the atlas (Loveless et al., Clin Cancer Res 2025).")),
        if (isTRUE(D$has_coexp)) nav_panel("Co-expression",
          strip(selectInput("cx_cell", "Cell type", choices = NULL, width = "300px"), uiOutput("cx_partner_ui")),
          layout_columns(col_widths = breakpoints(sm = 12, xl = c(5, 7)),
            plotlyOutput("p_coexp", height = "auto"),
            conditionalPanel("output.cx_has", plotlyOutput("p_coexp_pair", height = "460px"), fig_downloads("dl_coexp"), DTOutput("t_coexp", fill = FALSE))),
          note("Genes whose expression rises and falls with the selected gene across patients, within one cell type: for every treatment-naive primary tumour with at least ",
               MIN_CELLS, " cells of that type, the cell type's average expression; each gene is centred within its study to remove study effects, then genes are compared by ",
               "Spearman's rank correlation (rho). Only genes detected in at least ", paste0(round(100 * COEXP_MIN_FRAC), "%"), " of the cell type's cells in at least half of the samples are analysed. ",
               "FDR = false discovery rate across all genes analysed in that cell type. Click a gene in the table to open it. ",
               "Co-expression can reflect shared regulation, a shared cell state or a changing mix of subtypes within the cell type; it does not show that the proteins interact. ",
               "Transcripts typical of another cell type (for example tumour keratins in immune cells) can come from ambient RNA; their co-expression then mostly tracks that contamination.")),
        nav_panel("Atlas map",
          layout_columns(col_widths = c(6, 6), plotlyOutput("p_umap", height = "560px"), plotlyOutput("p_umap_ct", height = "560px")),
          fig_downloads("dl_umap"),
          note("A stratified subsample of atlas cells, positioned by the published atlas UMAP. Grey = not detected.")),
        nav_panel("Data",
          note("Per-sample single-cell values for the selected gene (all tissue states). Use the buttons to copy or download."),
          DTOutput("t_sc", fill = FALSE), hr(), note("CPTAC matched tumour/normal values."), DTOutput("t_cptac", fill = FALSE))
      )),
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
  updateSelectizeInput(session, "gene", choices = gene_df, selected = random_start_gene(), server = TRUE,
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
  fname <- function(stub, ext) function() paste0(gene(), "_", stub, "_", format(Sys.Date()), ".", ext)
  pr <- reactive(cptac_pairs(D, gene(), "protein")); rn <- reactive(cptac_pairs(D, gene(), "rna"))
  bk <- reactive(bulk_values(D, gene())); sv <- reactive(sc_values(D, gene()))
  uv <- reactive(umap_values(D, gene_row(D, gene())$sym_sc))
  cs <- reactive({ v <- sv(); if (is.null(v)) NULL else celltype_summary(v, input$context) })
  cs_naive <- reactive({ v <- sv(); if (is.null(v)) NULL else celltype_summary(v, "Primary tumour") })
  ph <- reactive(phospho_sites(D, gene()))
  ph_raw <- reactive(phospho_table(ph(), relative = FALSE))
  ph_rel <- reactive({ x <- ph(); if (is.null(x) || !x$has_protein) NULL else phospho_table(x, relative = TRUE) })
  ph_rel_p <- reactive({ x <- ph(); if (is.null(x) || !x$has_protein) NULL else phospho_table(x, relative = TRUE, ci = FALSE) })   # summary: p-values only
  ph_relmode <- reactive(identical(input$ph_mode, "rel"))
  ph_tab <- reactive(if (ph_relmode()) ph_rel() else ph_raw())          # NULL: no sites (or, relative, no protein)
  cv <- reactive(clinical_values(D, gene()))
  sv_tab <- reactive(surv_table(cv())); cl_assoc <- reactive(clin_assoc(cv()))
  surv_sel <- reactive({ k <- strsplit(if (is.null(input$surv_set)) SURV_CHOICES[[1]] else input$surv_set, "|", fixed = TRUE)[[1]]; list(dataset = k[1], endpoint = k[2]) })
  surv_x <- reactive({ x <- cv(); if (is.null(x)) NULL else surv_data(x, surv_sel()$dataset, surv_sel()$endpoint) })
  km <- reactive({ x <- surv_x(); if (is.null(x) || !nrow(x)) NULL else km_split(x) })
  in_moffitt <- reactive({ g <- gene(); g %in% MOFFITT_GENES })
  cx <- reactive(coexp_partners(D, gene()))
  cx_default <- reactive({ x <- cx(); if (is.null(x)) return(NULL); have <- unique(as.character(x$compartment))
    s <- cs_naive(); if (!is.null(s)) { sm <- s$summary[!is.na(s$summary$share), ]; o <- as.character(sm$compartment[order(-sm$share)]); o <- o[o %in% have]; if (length(o)) return(o[1]) }
    have[1] })
  cx_cell <- reactive({ x <- cx(); req(x); if (isTRUE(input$cx_cell %in% x$compartment)) input$cx_cell else cx_default() })
  cx_tab <- reactive({ x <- cx(); if (is.null(x)) NULL else x[x$compartment == cx_cell(), ] })
  cx_info <- reactive(D$coexp_cells[D$coexp_cells$compartment == cx_cell(), ])
  cx_partner <- reactive({ t <- cx_tab(); req(t, nrow(t)); if (isTRUE(input$cx_partner %in% t$partner)) input$cx_partner else t$partner[1] })
  cx_pair <- reactive({ t <- cx_tab(); req(t); pt <- cx_partner(); ps <- t$partner_sc[match(pt, t$partner)]; req(!is.na(ps))
    coexp_pair_values(D, gene(), D$genes$symbol[match(ps, D$genes$sym_sc)], cx_cell()) })
  ph_pairs_n <- if (isTRUE(D$has_phospho)) length(intersect(D$phospho$info$colnames$phospho_tumor, D$phospho$info$colnames$phospho_normal)) else 0L

  figs <- list(
    prot = reactive(plot_cptac_paired(pr(), gene(), "protein")), rna = reactive(plot_cptac_paired(rn(), gene(), "RNA")),
    bulk = reactive(plot_bulk(bk(), gene())), ct = reactive(plot_celltypes(cs(), gene(), input$context)),
    share = reactive(plot_share(cs(), gene())), frac = reactive(plot_frac(cs(), gene())),
    states = reactive(plot_states(sv(), gene(), input$compartments)), fibro = reactive(plot_fibro(sv(), gene())),
    ductal = reactive(plot_ductal(sv(), gene())), umap = reactive(plot_umap_static(D$umap, gene(), uv())),
    phsites = reactive(plot_phospho_sites(ph_tab(), gene(), relative = ph_relmode(), has_sites = !is.null(ph()))),
    phsite = reactive(plot_cptac_paired(ph_site_pairs(FALSE), gene(), "phosphosite")),
    km = reactive(plot_km(km(), surv_x(), gene(), surv_sel()$dataset, surv_sel()$endpoint)), hr = reactive(plot_hr(sv_tab(), gene())),
    stage = reactive(plot_stage(cv(), cl_assoc(), gene())), moffitt = reactive(plot_moffitt(cv(), cl_assoc(), gene(), in_moffitt())),
    coexp = reactive(if (is.null(cx())) empty_plot(paste0(gene(), ": not analysed in any cell type (not in the single-cell data or expressed in too few cells)")) else plot_coexp(cx_tab(), gene(), cx_cell(), cx_info())),
    coexp_pair = reactive(if (is.null(cx())) empty_plot("") else plot_coexp_pair(cx_pair(), gene(), cx_partner(), cx_cell())),
    phsite_rel = reactive(if (isTRUE(ph()$has_protein)) plot_cptac_paired(ph_site_pairs(TRUE), gene(), "phosphosite") else empty_plot("Protein not quantified")))

  output$gene_header <- renderUI({
    gr <- gene_row(D, gene())
    chip <- function(lbl, ok) span(class = if (ok) "chip chip-on" else "chip chip-off", lbl)
    div(class = "gene-header",
        h2(gr$symbol, span(class = "gene-name", if (!is.na(gr$name)) gr$name else "")),
        if (nzchar(gr$aliases)) p(class = "aliases", "Also known as: ", gr$aliases),
        div(chip("Protein (CPTAC)", !is.na(gr$sym_cptac_protein)), chip("RNA (CPTAC)", !is.na(gr$sym_cptac_rna)),
            if (isTRUE(D$has_phospho)) chip("Phosphosites (CPTAC)", gr$symbol %in% D$phospho_map$symbol),
            chip("RNA (TCGA/GTEx)", !is.na(gr$sym_tcga_gtex)), chip("Single-cell atlas", !is.na(gr$sym_sc) && isTRUE(D$has_sc))),
        div(class = "links-out", span(class = "links-label", "Elsewhere:"),
            lapply(gene_links(gr), function(l) tags$a(href = l$url, target = "_blank", rel = "noopener noreferrer", title = l$hint, class = "link-out", l$label))))
  })

  output$summary_ui <- renderUI({
    g <- gene(); sp <- paired_stats(pr()); sr <- paired_stats(rn())
    items <- list(txt_cptac(pr(), sp, g, "protein"), txt_cptac(rn(), sr, g, "rna"),
                  if (isTRUE(D$has_phospho)) txt_phospho(ph_raw(), ph_rel_p(), g, ph_pairs_n), txt_bulk(bk(), g),
                  if (isTRUE(D$has_clinical)) txt_survival(sv_tab(), g), if (isTRUE(D$has_clinical)) txt_clin_assoc(cl_assoc(), g, in_moffitt()),
                  txt_sc(cs_naive(), g), txt_state(sv(), g),
                  if (isTRUE(D$has_coexp) && !is.null(cx())) txt_coexp(cx(), g, cx_default(), D$coexp_cells[D$coexp_cells$compartment == cx_default(), ]))
    items <- Filter(Negate(is.null), items)
    tagList(tags$ul(class = "summary-list", lapply(items, tags$li)),
            note("These sentences describe the data; they do not establish mechanism or causality. p-values are unadjusted and refer to the single comparison stated. ",
                 "See the other tabs for the full distributions and the About page for methods and limitations."))
  })

  output$p_prot <- renderPlotly(to_plotly(figs$prot(), wrap = 70)); output$p_rna <- renderPlotly(to_plotly(figs$rna(), wrap = 70))
  output$p_bulk <- renderPlotly(to_plotly(figs$bulk(), wrap = 70)); output$p_ct <- renderPlotly(to_plotly(figs$ct(), wrap = 90))
  output$p_share <- renderPlotly(to_plotly(figs$share())); output$p_frac <- renderPlotly(to_plotly(figs$frac()))
  output$p_states <- renderPlotly(to_plotly(figs$states())); output$p_fibro <- renderPlotly(to_plotly(figs$fibro(), wrap = 100))
  output$p_ductal <- renderPlotly(to_plotly(figs$ductal(), wrap = 75))

  # ---- clinical context ----
  output$p_km <- renderPlotly(to_plotly(figs$km(), wrap = 95)); output$p_hr <- renderPlotly(to_plotly(figs$hr(), wrap = 78))
  output$p_stage <- renderPlotly(to_plotly(figs$stage(), wrap = 85)); output$p_moffitt <- renderPlotly(to_plotly(figs$moffitt(), wrap = 85))
  clin_csv <- reactive({ x <- cv(); if (is.null(x)) return(data.frame(note = paste("No tumour data for", gene())))
    x$value <- round(x$value, 3); names(x)[names(x) == "value"] <- paste0(gene(), "_value")
    names(x)[names(x) == "moffitt"] <- "moffitt_basal_minus_classical"; x })

  # ---- co-expression ----
  observe({ x <- cx(); if (is.null(x)) { updateSelectInput(session, "cx_cell", choices = character(0)); return() }
    have <- intersect(CELLTYPE_ORDER, as.character(x$compartment))
    updateSelectInput(session, "cx_cell", choices = have, selected = if (isTRUE(isolate(input$cx_cell) %in% have)) isolate(input$cx_cell) else cx_default()) })
  output$cx_partner_ui <- renderUI({ t <- cx_tab(); if (is.null(t) || !nrow(t)) return(NULL)
    ch <- stats::setNames(t$partner, paste0(t$partner, "  (rho = ", fmt_num(t$r, 3), ")"))
    selectInput("cx_partner", "Patient samples for", choices = ch, selected = if (isTRUE(isolate(input$cx_partner) %in% t$partner)) isolate(input$cx_partner) else t$partner[1], width = "300px") })
  output$cx_has <- reactive(!is.null(cx())); outputOptions(output, "cx_has", suspendWhenHidden = FALSE)
  cx_n <- reactive({ t <- cx_tab(); if (is.null(t)) 0L else min(sum(t$r > 0), COEXP_SHOW) + min(sum(t$r < 0), 10L) })
  output$p_coexp <- renderPlotly(to_plotly(figs$coexp(), wrap = 80, height = 190 + 20 * max(4, cx_n())))
  output$p_coexp_pair <- renderPlotly({ validate(need(!is.null(cx()), "")); to_plotly(figs$coexp_pair(), wrap = 100) })
  cx_table <- reactive({ x <- cx(); req(x); x <- x[x$compartment == cx_cell(), ]
    data.frame(gene = x$partner, name = D$genes$name[match(x$partner, D$genes$symbol)], rho = x$r, fdr = x$q,
               direction = ifelse(x$r > 0, "positive", "negative"), cell_type = as.character(x$compartment), stringsAsFactors = FALSE) })
  output$t_coexp <- renderDT({ t <- cx_table()
    t$open <- sprintf("<a href=\"#\" class=\"goto-gene\" data-gene=\"%s\">open</a>", htmltools::htmlEscape(t$gene))
    datatable(t[, c("gene", "open", "name", "rho", "fdr", "direction", "cell_type")], rownames = FALSE, escape = -2, style = "default", class = "compact stripe hover",
              extensions = "Buttons", colnames = c("gene", "", "name", "rho", "FDR", "direction", "cell type"),
              callback = JS("table.on('click', 'a.goto-gene', function(e) { e.preventDefault(); Shiny.setInputValue('goto_gene', $(this).data('gene'), {priority: 'event'}); });"),
              options = list(pageLength = 10, dom = "Bfrtip", buttons = dt_buttons(gene(), paste0("coexpression_", gsub("[^A-Za-z]+", "_", cx_cell()))),
                             columnDefs = list(list(orderable = FALSE, targets = 1)))) |>
      formatRound("rho", 3, mark = "") |> formatSignif("fdr", 2) }, server = FALSE)
  observeEvent(input$goto_gene, { g <- input$goto_gene; if (isTRUE(g %in% D$genes$symbol)) updateSelectizeInput(session, "gene", choices = gene_df, selected = g, server = TRUE) })
  cx_csv <- reactive({ x <- cx(); if (is.null(x)) return(data.frame(note = paste(gene(), "was not analysed in any cell type")))
    data.frame(gene = gene(), cell_type = as.character(x$compartment), partner = x$partner, rho = x$r, fdr = signif(x$q, 2)) })

  # ---- phosphosites ----
  ph_n_shown <- reactive({ t <- ph_tab(); if (is.null(t)) 0L else min(sum(!is.na(t$p)), PHOSPHO_MAX_PLOT) })
  output$p_phsites <- renderPlotly(to_plotly(figs$phsites(), wrap = 80, height = 190 + 22 * max(4, ph_n_shown())))
  output$ph_site_ui <- renderUI({
    t <- ph_tab(); if (is.null(t)) t <- ph_raw(); if (is.null(t) || !any(!is.na(t$p))) return(NULL)
    tt <- t[!is.na(t$p), ]; ch <- stats::setNames(tt$label, paste0(tt$label, "  (p = ", fmt_p(tt$p), ", ", tt$n_pairs, " pairs)"))
    sel <- if (isTRUE(isolate(input$ph_site) %in% tt$label)) isolate(input$ph_site) else tt$label[which.min(tt$p)]
    div(class = "ph-site-select", selectInput("ph_site", "Patient-level values for site", choices = ch, selected = sel, width = "360px"))
  })
  ph_site_pairs <- function(rel) {       # the selected site, both value types side by side
    x <- ph(); req(x, input$ph_site); k <- match(input$ph_site, x$sites$label); req(!is.na(k))
    list(pairs = phospho_pairs(x, k, rel), assay = if (rel) "phospho_rel" else "phospho", site = input$ph_site)
  }
  output$ph_has <- reactive({ t <- ph_raw(); !is.null(t) && any(!is.na(t$p)) }); outputOptions(output, "ph_has", suspendWhenHidden = FALSE)
  ph_any <- reactive({ t <- ph_raw(); !is.null(t) && any(!is.na(t$p)) })
  output$p_phsite <- renderPlotly({ validate(need(ph_any(), "")); to_plotly(figs$phsite(), wrap = 42) })
  output$p_phsite_rel <- renderPlotly({ validate(need(ph_any(), "")); to_plotly(figs$phsite_rel(), wrap = 42) })
  ph_table_view <- reactive({ t <- ph_tab(); req(t)
    t[order(is.na(t$p), match(t$label, ph()$sites$label)), c("label", "site_id", "peptide", "n_tumours", "n_pairs", "median_normal", "median_tumour", "median_diff", "ci_low", "ci_high", "p", "q", "n_up", "n_down", "values")] })
  output$t_phospho <- renderDT(datatable(ph_table_view(), rownames = FALSE, style = "default", class = "compact stripe hover", extensions = "Buttons",
                                         colnames = c("site", "site ID", "peptide(s)", "tumours quantified", "pairs", "median normal", "median tumour", "median change", "CI low", "CI high", "p", "FDR", "pairs up", "pairs down", "values"),
                                         options = list(pageLength = 10, dom = "Bfrtip", scrollX = TRUE, buttons = dt_buttons(gene(), "CPTAC_phosphosites"))) |>
                                 formatRound(c("median_normal", "median_tumour", "median_diff", "ci_low", "ci_high"), 3, mark = "") |>
                                 formatSignif(c("p", "q"), 2), server = FALSE)
  ph_long <- reactive({ x <- ph()               # per-patient values of every site (both value types) for the CSV download
    if (is.null(x)) return(data.frame(note = paste("No", gene(), "phosphosites were quantified in the CPTAC data.")))
    do.call(rbind, lapply(seq_len(nrow(x$sites)), function(k) {
      a <- phospho_pairs(x, k, FALSE); if (!nrow(a)) return(NULL)
      b <- if (x$has_protein) phospho_pairs(x, k, TRUE) else a[0, ]
      r <- data.frame(site = x$sites$label[k], site_id = x$sites$site_id[k], a, normal_rel_protein = b$normal[match(a$case_id, b$case_id)],
                      tumour_rel_protein = b$tumour[match(a$case_id, b$case_id)], diff_rel_protein = b$diff[match(a$case_id, b$case_id)])
      num <- vapply(r, is.double, logical(1)); r[num] <- lapply(r[num], round, 3); r       # inputs have 3 decimals
    })) })
  ph_h <- reactive(2.2 + 0.2 * max(4, ph_n_shown()))
  output$dl_phospho_pdf <- downloadHandler(fname("CPTAC_phosphosites", "pdf"), function(file) {
    grDevices::pdf(file, width = 9, height = ph_h()); print(for_export(figs$phsites(), 9))
    if (ph_any()) { print(for_export(figs$phsite(), 9)); if (isTRUE(ph()$has_protein)) print(for_export(figs$phsite_rel(), 9)) }; grDevices::dev.off() })
  output$dl_phospho_png <- downloadHandler(fname("CPTAC_phosphosites", "png"), function(file) ggsave(file, for_export(figs$phsites(), 9), width = 9, height = ph_h(), dpi = 200, bg = "white"))
  output$dl_phospho_csv <- downloadHandler(fname("CPTAC_phosphosites_per_patient", "csv"), function(file) utils::write.csv(ph_long(), file, row.names = FALSE))

  output$p_umap <- renderPlotly({
    u <- D$umap; validate(need(!is.null(u), "Atlas map not available.")); col <- gene_row(D, gene())$sym_sc
    validate(need(!is.na(col), paste(gene(), "is not in the single-cell data.")))
    val <- uv(); validate(need(!is.null(val), paste(gene(), "is not on the atlas map."))); o <- order(val)
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
  # whole table sent to the browser (server = FALSE; <= ~2,500 rows) so Copy/CSV export all rows, not just the visible page;
  # numbers are formatted for display (strings): the server's Buttons 1.0.3 cannot export raw numeric cells
  output$t_sc <- renderDT(datatable(sc_table(), rownames = FALSE, style = "default", class = "compact stripe hover", extensions = "Buttons",
                                    options = list(pageLength = 10, dom = "Bfrtip", buttons = dt_buttons(gene(), "single_cell_per_sample"))) |>
                            formatRound("n_cells", 0, mark = "") |> formatRound(c("mean_cp10k", "avg_expr"), 3, mark = "") |> formatRound("frac", 4, mark = ""), server = FALSE)
  output$t_cptac <- renderDT(datatable(cptac_table(), rownames = FALSE, style = "default", class = "compact stripe hover", extensions = "Buttons", options = list(pageLength = 10, dom = "Bfrtip", buttons = dt_buttons(gene(), "CPTAC_tumour_vs_normal"))) |>
                               formatRound(c("normal", "tumour", "diff"), 3, mark = ""), server = FALSE)

  # downloads: figure(s) + underlying data
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
             dl_surv = list(stub = "survival", plots = function() list(figs$km(), figs$hr()), data = clin_csv, w = 8, h = 6),
             dl_stage = list(stub = "stage", plots = function() list(figs$stage()), data = clin_csv, w = 11, h = 5),
             dl_moffitt = list(stub = "classical_basal", plots = function() list(figs$moffitt()), data = clin_csv, w = 11, h = 5),
             dl_coexp = list(stub = "coexpression", plots = function() list(figs$coexp(), figs$coexp_pair()), data = cx_csv, w = 8, h = 9),
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
