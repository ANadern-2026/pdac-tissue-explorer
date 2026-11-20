# ---- cell populations: two-level labels within atlas compartments (pop_<id>.sqlite, data-raw/07_cell_populations.R) ----
# Level 1 (lineages) are the atlas authors' published labels (myeloid) or named here (fibroblasts); level 2 is named here by
# literature signatures. Every population carries its definition, signature genes, label source and references, which the
# page shows next to every figure. The notes on how the groups were named (why a population is a name, only a score, set
# apart or not used; data-raw/reference/decisions.csv, attached by 05_build_store.R) are shown with their references on the
# cards they concern (collapsed to the headline) and in the section "How the groups were named". Page text is written for
# readers who may want the background but do not need it (user 2026-11-20): plain statements, no labels such as
# "Decision:", no workflow language.
# Plots work per sample (>= MIN_CELLS cells of a population), as in the rest of the app.

POP_MIN_PARENT <- 50L      # composition: a sample counts only with at least this many cells of the parent population
POP_CAT  <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")   # fixed order (dataviz reference)
POP_GREY <- "#a8a7a1"      # unresolved / unlabelled cells
POP_ATLAS <- c("Myeloid cells" = "myeloid", "Mast cells" = "myeloid", "Fibroblasts" = "fibroblasts")   # atlas-map cell type -> population store
POP_METHOD_TXT <- c(published = "Published label", signature = "Named from a published signature", score_only = "Shown as a score",
                    boundary = "Other cell type", unresolved = "Not named", unlabelled = "No published label")
POP_KIND_TXT <- c(structure = "Grouping", rule = "General", signature = "Signature", boundary = "Other cell types",
                  score_only = "Shown as a score", not_labelled = "Not used")
# what the atlas publication's own clusters of a cell type were (shown with the comparison table where groups are named here)
POP_ATLAS_NOTE <- c(Fibroblasts = paste(
  "The atlas publication (Loveless et al., Clin Cancer Res 2025) divided fibroblasts into eight clusters by two pairs of marker genes:",
  "myFb 1 to 5, myofibroblast-like (ACTA2, COL11A1), and iFb 1 to 3, inflammatory (C7, DPT).",
  "The groups on this page were named independently, from the signatures on the cards.",
  "The table shows where the cells of each group fall among the publication's clusters."))
# plain-language versions of the reasons and notes stored with each group of cells (07_cell_populations.R)
pop_plain_reason <- function(x) {
  x <- sub("^small cluster$", "too few cells (fewer than 50)", x)
  x <- sub("^few genes detected per cell for its study \\(possible low-quality cells\\)$", "fewer genes detected per cell than usual for its study (possibly damaged cells)", x)
  x <- sub("^few genes detected per cell \\(possible low-quality cells\\)$", "fewer genes detected per cell than in the rest of these cells (possibly damaged cells)", x)   # stores built before 2026-11-13
  x <- sub("^no signature reaches the minimum score$", "no signature fitted well enough", x)
  x <- sub("^two signatures score almost equally$", "two signatures fitted almost equally well", x)
  x <- sub("^(.*) scores higher than its signature \\(possible technical effect\\)$", "the \\1 was stronger than the signature (possibly from tissue handling)", x)
  x <- sub("^the ([A-Z])([a-z])", "the \\L\\1\\E\\2", x, perl = TRUE)
  sub(" \\(score only\\)", "", x)
}
pop_plain_note <- function(x) {
  parts <- unlist(strsplit(x, "; ", fixed = TRUE)); if (!length(parts)) return(character(0))
  parts <- sub("^resembles (.*) by the level-1 markers$", "looks more like \\1", parts)
  parts <- sub("^strongest programme: (.*) \\(score [-0-9.]+\\)$", "a strong \\1", parts)
  parts <- gsub(" \\(score only\\)", "", parts)
  vapply(parts, function(q) sub("(looks more like |a strong )([A-Z])([a-z])", "\\1\\L\\2\\E\\3", q, perl = TRUE), "", USE.NAMES = FALSE)
}
pop_and <- function(x) if (length(x) < 2) x else paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])
# the notes of one or more groups as one sentence part: "the A and B programmes" / "look more like X"
pop_note_text <- function(parts, verb_show = "show", verb_like = "look more like") {
  parts <- unique(parts); progs <- sub("^a strong ", "", parts[grepl("^a strong ", parts)]); like <- sub("^looks more like ", "", parts[grepl("^looks more like ", parts)])
  re <- "^(.*) programme(, (.+)| \\((.+)\\))?$"                # "X programme", "X programme (Y)" or, in older stores, "X programme, Y"
  if (length(progs) > 1 && all(grepl(re, progs))) {
    alias <- paste0(sub(re, "\\3", progs), sub(re, "\\4", progs))
    progs <- paste0(pop_and(paste0(sub(re, "\\1", progs), ifelse(nzchar(alias), paste0(" (", alias, ")"), ""))), " programmes")
  } else progs <- pop_and(progs)
  c(if (length(progs)) paste(verb_show, "the", progs), if (length(like)) paste(verb_like, pop_and(like)))
}
pop_group_name <- function(id) paste0("Group ", sub("^.*_c", "", id))
POP_DEC_EMPTY <- data.frame(decision_id = character(0), compartment = character(0), applies_to = character(0), kind = character(0),
                            question = character(0), decision = character(0), reason = character(0), refs = character(0), stringsAsFactors = FALSE)

pop_load <- function(dir) {
  out <- list()
  for (f in sort(list.files(dir, pattern = "^pop_.*\\.sqlite$", full.names = TRUE))) {
    st <- store_open(f); i <- st$info
    tax <- i$taxonomy; tax$level <- as.integer(tax$level)
    calls <- i$cluster_calls
    extra_ids <- setdiff(unique(c(calls$call, i$groups$population_id, i$map_meta$l1, i$map_meta$l2)), tax$population_id)
    root <- tax$population_id[tax$level == 0][1]
    extra <- do.call(rbind, lapply(extra_ids, function(id) {
      unr <- grepl("_UNRESOLVED$", id); parent <- if (unr) sub("_UNRESOLVED$", "", id) else root
      plab <- tax$label[match(parent, tax$population_id)]
      data.frame(population_id = id, compartment = i$compartment, level = if (unr && parent != root) 2L else 1L, parent_id = parent,
                 label = if (unr) paste0(plab, ": not named") else paste0(i$compartment, ": no published label"),
                 definition = if (unr) paste0("Cells among the ", sub("^([A-Z])([a-z])", "\\L\\1\\E\\2", plab, perl = TRUE), " that no published signature fitted clearly.")
                              else "Cells of this cell type that the atlas publication did not assign to a lineage.",
                 signature_genes = "", signature_ref = "", definition_refs = "", pdac_evidence = "",
                 label_method = if (unr) "unresolved" else "unlabelled", status = "", also_known_as = "", stringsAsFactors = FALSE)
    }))
    tax <- rbind(tax[, names(tax)], if (!is.null(extra)) extra[, names(tax)])
    tax$colour <- POP_GREY
    for (par in unique(tax$parent_id[tax$level > 0])) {          # fixed colour per population: taxonomy order within its parent
      k <- which(tax$parent_id == par & tax$label_method %in% c("published", "signature"))
      tax$colour[k] <- POP_CAT[(seq_along(k) - 1L) %% length(POP_CAT) + 1L]
    }
    g <- i$groups; g$state <- tissue_state(g$disease_state, g$treatment, g$met_site)
    m <- i$map_meta; m$state <- tissue_state(m$disease_state, m$treatment, m$met_site)
    dec <- i$decisions; if (!is.data.frame(dec) || !nrow(dec)) dec <- POP_DEC_EMPTY     # stores built before the decision record
    l1m <- if (is.null(i$level1_method)) "published" else i$level1_method
    ac <- if (is.data.frame(i$atlas_compare) && nrow(i$atlas_compare)) i$atlas_compare else NULL
    out[[i$id]] <- list(id = i$id, name = i$compartment, st = st, tax = tax, refs = i$references, decisions = dec, groups = g, map = m, scores = i$map_scores,
                        calls = calls, cluster_scores = i$cluster_scores, markers = i$markers, check = i$label_check, params = i$params,
                        built = i$built_07, level1_method = l1m, atlas_compare = ac,
                        words = if (l1m == "signature") c(l1 = "groups", l1_one = "group", l2 = "sub-groups") else c(l1 = "lineages", l1_one = "lineage", l2 = "sub-lineages"))
  }
  out
}

pop_lc <- function(x) ifelse(grepl("^[A-Z][a-z]", x), paste0(tolower(substr(x, 1, 1)), substring(x, 2)), x)   # 'Inflammatory fibroblasts (iCAF family)' -> 'inflammatory ... (iCAF family)'
pop_row <- function(cp, id) cp$tax[match(id, cp$tax$population_id), , drop = FALSE]
pop_label <- function(cp, id) cp$tax$label[match(id, cp$tax$population_id)]
pop_root <- function(cp) cp$tax$population_id[cp$tax$level == 0][1]
# populations shown at a level: lineages of the compartment, or the sub-lineages (labels + unresolved) of one lineage
pop_kids <- function(cp, parent, with_scores = FALSE) {
  t <- cp$tax[cp$tax$parent_id %in% parent & cp$tax$population_id != parent, ]
  have <- unique(c(cp$groups$population_id, cp$map$l1, cp$map$l2))
  keep <- if (with_scores) t$label_method == "score_only" else t$label_method != "score_only" & t$population_id %in% have
  t$population_id[keep]
}
pop_drillable <- function(cp, id) length(pop_kids(cp, id)) > 0 && any(!is.na(cp$map$umap2_x[cp$map$l1 == id]))

# one gene: per-sample values of every population (both levels); NULL if the gene is not in the single-cell data
pop_values <- function(d, cp, symbol) {
  s <- gene_row(d, symbol)$sym_sc; if (length(s) != 1L || is.na(s)) return(NULL)
  rec <- store_get(cp$st, "expr", s); if (is.null(rec)) return(NULL)
  v <- data.frame(cp$groups, mean_cp10k = rec$mean / 100, frac = rec$frac / 1e4)
  v$avg_expr <- log1p(v$mean_cp10k); v$label <- pop_label(cp, v$population_id); v
}
pop_map_values <- function(d, cp, symbol) {
  s <- gene_row(d, symbol)$sym_sc; if (length(s) != 1L || is.na(s)) return(NULL)
  store_get(cp$st, "map", s)
}

# ---- figures ----
pop_caption <- function(cp) src_caption(if (identical(cp$level1_method, "signature"))
  "Atlas: Loveless et al., Clin Cancer Res 2025; groups and sub-groups named here from published signatures (sources on the Cell populations page)" else
  "Atlas: Loveless et al., Clin Cancer Res 2025; lineages as published, sub-lineages named here from published signatures (sources on the Cell populations page)")
pop_boundary <- function(cp) cp$tax$population_id[cp$tax$label_method == "boundary"]   # other cell types: shown, not counted in the make-up
pop_xlab <- function(cp, ids, n) factor(paste0(pop_label(cp, ids), "\n(n = ", n, ")"))

plot_pop_expr <- function(cp, v, kids, symbol, state) {
  if (is.null(v)) return(empty_plot(paste0(symbol, " is not in the single-cell data")))
  if (!length(state)) return(empty_plot("Tick at least one tissue state above the plot"))
  x <- v[v$population_id %in% kids & v$state %in% state & v$n_cells >= MIN_CELLS, ]
  if (!nrow(x)) return(empty_plot(paste0("No samples with at least ", MIN_CELLS, " cells of these populations in this tissue state")))
  n <- table(factor(x$population_id, levels = kids)); kids <- kids[n > 0]
  lv <- paste0(vapply(pop_label(cp, kids), function(l) paste(strwrap(l, 12), collapse = "\n"), ""), "\n(n = ", as.integer(n[kids]), ")")
  x$pop <- factor(lv[match(x$population_id, kids)], levels = lv)
  cols <- stats::setNames(pop_row(cp, kids)$colour, lv)
  x$text <- paste0(x$label, " | sample ", x$sample, " (", x$n_cells, " cells)<br>average expression: ", fmt_num(x$avg_expr), "<br>cells expressing: ", fmt_pct(x$frac))
  ggplot(x, aes(pop, avg_expr)) +
    geom_jitter(aes(colour = pop, text = text), width = 0.2, height = 0, size = 1.3, alpha = 0.6) +
    geom_boxplot(fill = NA, outlier.shape = NA, width = 0.5, colour = INK[["primary"]], linewidth = 0.3) +
    scale_colour_manual(values = cols, guide = "none") +
    labs(x = NULL, y = "Expression per sample (log1p CP10k)", title = paste0(symbol, " in ", pop_lc(pop_label(cp, pop_row(cp, kids[1])$parent_id)), ": ", paste(tolower(state), collapse = ", ")),
         subtitle = paste0("One dot per sample with >= ", MIN_CELLS, " cells of the population; n = samples. Box = median and quartiles."),
         caption = pop_caption(cp)) + theme_px() + theme(axis.text.x = element_text(size = 9))
}

# finest labels of a compartment (sub-groups where they exist, else the group; boundary groups left out): the
# 'Subtypes' tab of the gene page shows these for fibroblasts instead of the atlas's own labels
pop_fine <- function(cp) { l1 <- setdiff(pop_kids(cp, pop_root(cp)), pop_boundary(cp))
  unlist(lapply(l1, function(k) if (pop_drillable(cp, k)) pop_kids(cp, k) else k)) }
plot_pop_fine <- function(cp, v, symbol, state = "Primary tumour") {
  p <- plot_pop_expr(cp, v, pop_fine(cp), symbol, state)
  if (!inherits(p, "ggplot") || is.null(p$labels$subtitle)) return(p)
  p + labs(title = paste0(symbol, " in ", sub("s$", "", pop_lc(cp$name)), " groups (treatment-naive primary tumours)"),
           subtitle = paste0("Finest groups of the Cell populations page; one dot per sample with >= ", MIN_CELLS, " cells; n = samples; grey = not named"))
}
pop_fine_table <- function(cp, v, symbol) {
  if (is.null(v)) return(data.frame(note = paste(symbol, "is not in the single-cell data")))
  x <- v[v$population_id %in% pop_fine(cp), ]
  data.frame(population = x$label, label_source = unname(POP_METHOD_TXT[pop_row(cp, x$population_id)$label_method]), sample = x$sample, study = x$study,
             state = as.character(x$state), n_cells = x$n_cells, mean_cp10k = round(x$mean_cp10k, 3), avg_expr = round(x$avg_expr, 3), frac = round(x$frac, 4),
             stringsAsFactors = FALSE, row.names = NULL)
}

# populations on x, the five tissue states side by side within each population (one panel: reads well on screen and in plotly)
pop_state_layout <- function(cp, x, kids) {
  kids <- kids[kids %in% x$population_id]; ns <- length(STATE_LEVELS); w <- 0.8 / ns
  x$pi <- match(x$population_id, kids); x$si <- as.integer(x$state)
  x$xpos <- x$pi + (x$si - (ns + 1) / 2) * w; x$grp <- paste(x$pi, x$si)
  list(x = x, w = w, breaks = seq_along(kids), labels = vapply(pop_label(cp, kids), function(l) paste(strwrap(l, 18), collapse = "\n"), ""))
}
STATE_KEY <- "Colour = tissue state: green healthy donor, blue adjacent normal, orange primary tumour, violet liver metastasis, pink treated primary tumour"
plot_pop_bystate <- function(L, y, ylab, title, subtitle, cp, yfmt = waiver()) {
  x <- L$x; set.seed(1); x$xj <- x$xpos + stats::runif(nrow(x), -L$w * 0.3, L$w * 0.3)
  ggplot(x) +
    geom_boxplot(aes(x = xpos, y = .data[[y]], group = grp), fill = NA, outlier.shape = NA, width = L$w * 0.85, colour = INK[["primary"]], linewidth = 0.3) +
    geom_point(aes(x = xj, y = .data[[y]], colour = state, text = text), size = 1.1, alpha = 0.65) +
    scale_colour_manual(values = STATE_COLS, guide = "none") + scale_x_continuous(breaks = L$breaks, labels = L$labels) + scale_y_continuous(labels = yfmt) +
    labs(x = NULL, y = paste(strwrap(ylab, 38), collapse = "\n"), title = title, subtitle = subtitle, caption = pop_caption(cp)) + theme_px() +
    theme(axis.text.x = element_text(size = 9), panel.grid.major.x = element_blank())
}
plot_pop_states <- function(cp, v, kids, symbol) {
  if (is.null(v)) return(empty_plot(paste0(symbol, " is not in the single-cell data")))
  x <- v[v$population_id %in% kids & !is.na(v$state) & v$n_cells >= MIN_CELLS, ]
  if (!nrow(x)) return(empty_plot("No samples with enough cells"))
  x$text <- paste0(x$label, " | ", x$state, "<br>sample ", x$sample, " (", x$n_cells, " cells)<br>average expression: ", fmt_num(x$avg_expr), "<br>cells expressing: ", fmt_pct(x$frac))
  plot_pop_bystate(pop_state_layout(cp, x, kids), "avg_expr", "Expression per sample (log1p CP10k)", paste0(symbol, " across tissue states, by population"),
                   paste0("One dot per sample with >= ", MIN_CELLS, " cells of the population. ", STATE_KEY, ".\nTissue states come from different studies (partly confounded)."), cp)
}

# share of each population among the parent's cells per sample (independent of the gene)
pop_composition <- function(cp, parent, kids) {
  lvl <- if (parent == pop_root(cp)) 1L else 2L
  kids <- setdiff(kids, pop_boundary(cp))
  g <- cp$groups[cp$groups$level == lvl & !is.na(cp$groups$state) & !cp$groups$population_id %in% pop_boundary(cp), ]
  if (lvl == 2L) g <- g[g$population_id %in% pop_kids(cp, parent), ]
  tot <- stats::aggregate(n_cells ~ sample, g, sum); tot <- tot[tot$n_cells >= POP_MIN_PARENT, ]
  if (!nrow(tot)) return(NULL)
  smp <- unique(g[g$sample %in% tot$sample, c("sample", "state", "study")])
  x <- merge(expand.grid(sample = smp$sample, population_id = kids, stringsAsFactors = FALSE), g[, c("sample", "population_id", "n_cells")], all.x = TRUE)
  x$n_cells[is.na(x$n_cells)] <- 0L
  x <- merge(merge(x, smp, by = "sample"), stats::setNames(tot, c("sample", "n_parent")), by = "sample")
  x$share <- x$n_cells / x$n_parent; x$label <- pop_label(cp, x$population_id); x
}
plot_pop_composition <- function(cp, parent, kids) {
  x <- pop_composition(cp, parent, kids)
  if (is.null(x)) return(empty_plot(paste0("No samples with at least ", POP_MIN_PARENT, " cells")))
  x$text <- paste0(x$label, " | ", x$state, "<br>sample ", x$sample, ": ", x$n_cells, " of ", x$n_parent, " cells (", fmt_pct(x$share), ")")
  plot_pop_bystate(pop_state_layout(cp, x, kids), "share", "Share of the cells in each sample", paste0("Make-up of ", pop_lc(pop_label(cp, parent)), " by tissue state"),
                   paste0("One dot per sample with >= ", POP_MIN_PARENT, " ", pop_lc(pop_label(cp, parent)), if (length(intersect(kids, pop_boundary(cp)))) " (cells set apart as other cell types not counted)",
                          ". ", STATE_KEY, ".\nShares also depend on how each study dissociated and sorted cells."),
                   cp, scales::percent)
}

# static map for downloads (the page shows an interactive version)
plot_pop_map_static <- function(df, title, colours = NULL, value = NULL) {
  if (!nrow(df)) return(empty_plot("No cells"))
  p <- ggplot(df, aes(x, y))
  if (is.null(value)) {
    cen <- stats::aggregate(cbind(x, y) ~ label, df, stats::median)
    p <- p + geom_point(aes(colour = label), size = 0.3) + scale_colour_manual(values = colours, name = NULL) +
      geom_label(data = cen, aes(label = label), size = 3, label.size = 0, fill = "white", alpha = 0.8) +
      guides(colour = guide_legend(override.aes = list(size = 3)))
  } else {
    df$value <- value; df <- df[order(df$value), ]
    p <- p + geom_point(aes(colour = value), size = 0.3) + scale_colour_gradientn(colours = c(INK[["zero"]], SEQ_RAMP[-1]), name = "value")
  }
  p + coord_equal() + labs(title = title, x = NULL, y = NULL, caption = src_caption("Subsample of cells; map computed in this resource (Harmony-integrated UMAP)")) +
    theme_px() + theme(axis.text = element_blank(), panel.grid = element_blank())
}

# ---- definition cards ----
pop_ref_ids <- function(r) {
  ids <- c(unlist(strsplit(r$definition_refs, "[;,]")), unlist(regmatches(r$signature_ref, gregexpr("[a-z]+[0-9]{4}", r$signature_ref))))
  unique(trimws(ids[nzchar(trimws(ids))]))
}
pop_split <- function(x) { v <- trimws(unlist(strsplit(x, ";", fixed = TRUE))); v[nzchar(v)] }
# notes shown on a population's card: the rows of decisions.csv that name it, collapsed to their headline
pop_decisions_for <- function(cp, id) { d <- cp$decisions; d[vapply(d$applies_to, function(a) id %in% pop_split(a), TRUE), , drop = FALSE] }
pop_note_body <- function(d) tagList(if (nzchar(d$decision)) p(d$decision), if (nzchar(d$reason)) p(d$reason))
pop_why <- function(d) if (nrow(d)) lapply(seq_len(nrow(d)), function(k)
  tags$details(class = "pop-why", tags$summary(d$question[k]), div(class = "pop-why-body", pop_note_body(d[k, ]))))
pop_ref_list <- function(cp, ids) {
  r <- cp$refs[match(ids, cp$refs$ref_id), ]; r <- r[!is.na(r$ref_id), ]
  if (!nrow(r)) return(NULL)
  tags$ul(class = "pop-refs", lapply(seq_len(nrow(r)), function(k) tags$li(r$citation[k], ". ",
    tags$a(href = paste0("https://doi.org/", r$doi[k]), target = "_blank", rel = "noopener noreferrer", paste0("doi:", r$doi[k])))))
}
# where groups are named here: how their cells had been classified in the atlas publication (table with an explanation)
pop_atlas_compare_ui <- function(cp, ids) {
  a <- cp$atlas_compare; if (is.null(a) || !identical(cp$level1_method, "signature") || is.na(POP_ATLAS_NOTE[cp$name])) return(NULL)
  a <- a[a$population_id %in% ids, ]; a <- a[a$level == pop_row(cp, a$population_id)$level, ]; if (!nrow(a)) return(NULL)
  cols <- sort(unique(a$atlas_label)); cols <- c(setdiff(cols, "(no atlas label)"), intersect(cols, "(no atlas label)"))
  rows <- ids[ids %in% a$population_id]
  tab <- vapply(rows, function(i) { r <- a[a$population_id == i, ]; sh <- r$n_cells / sum(r$n_cells); vapply(cols, function(cc) { v <- sum(sh[r$atlas_label == cc]); if (v > 0) fmt_pct(v) else "" }, "") }, character(length(cols)))
  tab <- matrix(tab, nrow = length(cols), dimnames = list(cols, rows))
  tags$details(class = "how pop-compare-box", tags$summary(paste0("How these groups compare with the ", tolower(cp$name), " clusters of the atlas publication")),
    div(class = "how-body", p(POP_ATLAS_NOTE[[cp$name]]),
      tags$table(class = "pop-compare", tags$thead(tags$tr(tags$th("Group on this page"), tags$th("Cells"), lapply(cols, tags$th))),
        tags$tbody(lapply(rows, function(i) tags$tr(tags$td(pop_label(cp, i)), tags$td(format(sum(a$n_cells[a$population_id == i]), big.mark = ",")),
                                                     lapply(cols, function(cc) tags$td(tab[cc, i]))))))))
}
pop_card <- function(cp, id, v = NULL, drill = FALSE) {
  r <- pop_row(cp, id); meth <- r$label_method
  g <- cp$groups[cp$groups$population_id == id, ]
  n_cells <- sum(g$n_cells); n_smp <- sum(g$n_cells >= MIN_CELLS)
  calls <- cp$calls[cp$calls$call == id, ]
  src <- switch(meth, published = "Label published with the atlas (Loveless et al., Clin Cancer Res 2025).",
                boundary = paste0("These cells express the genes of another cell type more strongly than those of ", pop_lc(cp$name),
                                  ". They are shown here but not counted in the make-up."),
                signature = paste0("Named here: ", if (nrow(calls) == 1) "one group" else paste(nrow(calls), "groups"),
                                   " of similar cells fitted this signature clearly better than any other."),
                score_only = "Shown as a score on the map, to see where these genes are active; not used to name cells.",
                unresolved = "", "")
  chk <- if (meth == "published" && id %in% names(cp$check$agreement)) paste0(fmt_pct(cp$check$agreement[[id]]),
           " of these cells express the marker genes of this lineage more strongly than those of the other ", sum(cp$tax$level == 1 & cp$tax$label_method == "published") - 1, " lineages.")
  notes <- unique(pop_plain_note(paste(calls$note[nzchar(calls$note)], collapse = "; ")))
  dec <- pop_decisions_for(cp, id); ref_ids <- unique(c(pop_ref_ids(r), pop_split(dec$refs)))
  
  reasons <- if (meth == "unresolved" && nrow(calls)) tags$details(class = "pop-refbox",
    tags$summary(paste0("Why ", if (nrow(calls) == 1) "this group was" else paste(nrow(calls), "groups were"), " not named")),
    tags$ul(class = "pop-clusters", lapply(seq_len(nrow(calls)), function(k)
      tags$li(paste0(pop_group_name(calls$cluster[k]), " (", format(calls$n_cells[k], big.mark = ","), " cells): ", pop_plain_reason(calls$reason[k]),
                     if (nzchar(calls$note[k])) paste0("; its cells ", paste(pop_note_text(pop_plain_note(calls$note[k]), "show", "looked more like"), collapse = " and ")), ". Closest: ",
                     pop_lc(pop_label(cp, calls$best[k])), ", then ", pop_lc(pop_label(cp, calls$runner_up[k])), ".")))))
  div(class = paste("pop-card", paste0("pop-", meth)),
      div(class = "pop-card-head", span(class = "pop-swatch", style = paste0("background:", r$colour)), span(class = "pop-name", r$label),
          span(class = paste("pop-badge", paste0("badge-", meth)), POP_METHOD_TXT[[meth]])),
      if (meth != "score_only") p(class = "pop-count", format(n_cells, big.mark = ","), " cells; ", n_smp, " samples with >= ", MIN_CELLS, " cells"),
      p(class = "pop-def", r$definition),
      if (nzchar(r$signature_genes)) p(class = "pop-sig", strong(if (meth == "published") "Marker genes: " else "Signature genes: "),
                                       gsub(";", ", ", r$signature_genes), span(class = "pop-src", paste0(" (", r$signature_ref, ")"))),
      if (nzchar(r$pdac_evidence)) p(class = "pop-pdac", strong("In PDAC: "), r$pdac_evidence),
      pop_why(dec),
      if (nzchar(r$also_known_as)) p(class = "pop-aka", strong("Also called: "), r$also_known_as),
      if (nzchar(src)) p(class = "pop-src-line", src), if (!is.null(chk)) p(class = "pop-src-line", chk),
      if (length(notes) && meth == "signature") p(class = "pop-note", paste0("Some of these cells also ", paste(pop_note_text(notes), collapse = ", and some "), ".")),
      reasons,
      if (length(ref_ids)) tags$details(class = "pop-refbox", tags$summary("References"), pop_ref_list(cp, ref_ids)),
      if (drill) tags$a(href = "#", class = "pop-open", onclick = sprintf("Shiny.setInputValue('pg_open', '%s', {priority: 'event'}); return false;", id),
                        HTML(paste0("Open ", htmltools::htmlEscape(pop_lc(r$label)), " &rsaquo;"))))
}

# all notes of a compartment (rows of decisions.csv: question = headline, decision = what was done, reason = background):
# its own notes by kind, then the notes shared by every cell type
pop_decision_item <- function(cp, d) {
  pops <- pop_split(d$applies_to); pops <- pops[pops %in% cp$tax$population_id]
  tags$details(class = "pop-decision",
    tags$summary(span(class = "dec-q", d$question), span(class = paste("pop-badge", paste0("dec-", d$kind)), POP_KIND_TXT[[d$kind]])),
    div(class = "dec-body", pop_note_body(d),
        if (length(pops)) p(class = "dec-applies", "Concerns: ", paste(pop_label(cp, pops), collapse = "; ")),
        if (length(pop_split(d$refs))) pop_ref_list(cp, pop_split(d$refs))))
}
pop_decision_ui <- function(cp) {
  d <- cp$decisions
  if (!nrow(d)) return(p(class = "pop-src-line", "Background notes are not part of this data version; the cards give the definitions and references."))
  own <- d[d$compartment != "all", , drop = FALSE]; own <- own[order(match(own$kind, names(POP_KIND_TXT)), own$decision_id), , drop = FALSE]
  all <- d[d$compartment == "all", , drop = FALSE]
  tagList(if (nrow(own)) tagList(h4(class = "dec-head", cp$name), lapply(seq_len(nrow(own)), function(k) pop_decision_item(cp, own[k, ]))),
          if (nrow(all)) tagList(h4(class = "dec-head", "For every cell type"), lapply(seq_len(nrow(all)), function(k) pop_decision_item(cp, all[k, ]))))
}
pop_decision_table <- function(cp) {
  d <- cp$decisions; if (!nrow(d)) return(NULL)
  data.frame(note = d$decision_id, scope = ifelse(d$compartment == "all", "all cell types", d$compartment), topic = unname(POP_KIND_TXT[d$kind]),
             headline = d$question, summary = d$decision, background = d$reason,
             applies_to = vapply(d$applies_to, function(a) paste(pop_label(cp, pop_split(a)), collapse = "; "), ""),
             references = vapply(d$refs, function(x) { rr <- cp$refs[match(pop_split(x), cp$refs$ref_id), ]; paste(paste0(rr$citation, " doi:", rr$doi), collapse = " | ") }, ""),
             stringsAsFactors = FALSE, row.names = NULL)
}

# ---- interactive maps (plotly; click a population to open it) ----
pop_map_plotly <- function(df, title, source, clickable = character(0)) {
  p <- plot_ly(source = source)
  for (lab in unique(df$label)) {
    s <- df[df$label == lab, ]
    p <- add_trace(p, x = round(s$x, 2), y = round(s$y, 2), type = "scattergl", mode = "markers", name = lab, customdata = s$id,
                   marker = list(size = 3, color = s$colour[1], opacity = if (s$id[1] %in% clickable || !length(clickable)) 0.8 else 0.35),
                   text = s$hover, hoverinfo = "text")
  }
  lab_df <- stats::aggregate(cbind(x, y) ~ label, df, stats::median)
  p |> add_annotations(x = lab_df$x, y = lab_df$y, text = lab_df$label, showarrow = FALSE, font = list(size = 11, color = "#0b0b0b"), bgcolor = "rgba(255,255,255,0.75)") |>
    layout(title = list(text = paste0("<b>", title, "</b>"), x = 0), xaxis = list(visible = FALSE), yaxis = list(visible = FALSE, scaleanchor = "x"),
           legend = list(orientation = "h", y = -0.02, font = list(size = 11)), margin = list(t = 50)) |>
    config(displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d", "toggleSpikelines")) |> event_register("plotly_click") |> pop_slim()
}
pop_slim <- function(p) { p <- plotly_build(p); p$x$attrs <- NULL; p$x$visdat <- NULL; p$x$cur_data <- NULL; p }   # sends the built traces only (about half the size)
pop_value_plotly <- function(df, value, title, cbar) {
  o <- order(value)
  plot_ly(x = round(df$x[o], 2), y = round(df$y[o], 2), type = "scattergl", mode = "markers",
          marker = list(size = 3, color = round(value[o], 2), colorscale = list(c(0, INK[["zero"]]), c(0.001, SEQ_RAMP[1]), c(1, SEQ_RAMP[6])), showscale = TRUE,
                        colorbar = list(title = cbar, thickness = 10)),
          text = paste0(df$label[o], "<br>", cbar, ": ", fmt_num(value[o])), hoverinfo = "text") |>
    layout(title = list(text = paste0("<b>", title, "</b>"), x = 0), xaxis = list(visible = FALSE), yaxis = list(visible = FALSE, scaleanchor = "x"), margin = list(t = 50)) |>
    config(displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d", "toggleSpikelines")) |> pop_slim()
}

# ---- page ----
pop_ui <- function() nav_panel("Cell populations", value = "pops", icon = nav_icon("pops"),
  div(class = "page-intro", h3("Cell populations"),
      p("See a gene across the finer cell populations of the single-cell atlas. Open a cell type on the map, then one of its groups; ",
        "each population has a card with its definition, the genes behind its name and the studies it comes from."),
      div(class = "pop-caveat", "The names come from the literature. Apart from the myeloid lineages, which the atlas authors published, ",
          "groups of similar cells are matched here to the gene signatures of published studies, and each card cites the study and figure behind its name. ",
          "Cells that fit no signature clearly are shown as not named.")),
  div(class = "gene-band",
      div(class = "gene-search", selectizeInput("pg_gene", "Gene", choices = NULL, width = "100%", options = list(placeholder = "Type a gene symbol, alias or name", maxOptions = 50))),
      div(class = "gene-band-header", uiOutput("pg_crumbs"))),
  layout_columns(col_widths = breakpoints(sm = 12, lg = c(6, 6)),
    div(plotlyOutput("pg_map", height = "560px"), note("Click a coloured population (or use the cards below) to open it.")),
    div(plotlyOutput("pg_map_val", height = "560px"), uiOutput("pg_colour_ui"))),
  conditionalPanel("output.pg_in_comp",
    div(class = "sec-head", h3(textOutput("pg_gene_head", inline = TRUE)), p("One dot per tissue sample with at least ", MIN_CELLS, " cells of the population.")),
    layout_columns(col_widths = breakpoints(sm = 12, xl = c(5, 7)),
      div(div(class = "pop-plot-control", checkboxGroupInput("pg_state", "Tissue states in this plot", choices = STATE_LEVELS, selected = "Primary tumour", inline = TRUE, width = "100%")),
          plotlyOutput("pg_expr", height = "500px")),
      div(div(class = "pop-plot-control pop-plot-control-note", "All five tissue states, side by side within each population"),
          plotlyOutput("pg_states", height = "500px"))),
    note("Average expression = log1p of the mean normalised expression (CP10k) of the population's cells in that sample. Grey: cells without a name."),
    div(class = "sec-head", h3("Make-up of the populations"), p("How the cells divide among these populations in each sample, by tissue state. This does not depend on the selected gene.")),
    plotlyOutput("pg_comp", height = "480px"),
    fig_downloads("dl_pop"),
    div(class = "sec-head", h3("About these populations")),
    uiOutput("pg_cards"),
    uiOutput("pg_atlas_cmp"),
    div(class = "sec-head", h3("How the groups were named"),
        p("Background on the choices behind the names, with the studies they rest on.")),
    uiOutput("pg_decisions"),
    div(class = "sec-head", h3("Genes that mark each population")),
    how_to(p("For each population, genes higher in it than in the rest of its parent population, compared within the same samples: median log ratio of the average ",
             "expression across samples with at least ", MIN_CELLS, " cells of both, a one-sided sign test across samples and false discovery rate (FDR) across genes; ",
             "only genes detected in at least 10% of the population's cells; up to 50 genes per population. Click a gene to show it on this page.")),
    DTOutput("pg_markers", fill = FALSE),
    tags$details(class = "how how-table", tags$summary("Per-sample values of the selected gene at this level (table and downloads)"),
                 div(class = "how-body", DTOutput("pg_table", fill = FALSE)))),
  uiOutput("pg_atlas_cards"))

pop_server <- function(input, output, session, D, gene, gene_df, save_figs) {
  cps <- D$pops
  st <- reactiveValues(level = "atlas", cp = names(cps)[1], lin = NULL)
  # own gene box, filled from the "Explore a gene" page: on start and whenever the user comes to this page with another gene there
  # a gene typed into this page's box returns the page to the atlas map (user 2026-11-15); genes set by the app (start, coming from
  # the gene page, clicking a gene in the marker table) keep the current view: `set_by_app` holds the gene the app is about to set
  set_by_app <- reactiveVal(NULL)
  observeEvent(input$pg_gene, ignoreInit = TRUE, { g <- input$pg_gene
    if (!isTRUE(g %in% D$genes$symbol)) return()
    if (identical(g, set_by_app())) { set_by_app(NULL); return() }
    if (st$level != "atlas") go("atlas") })
  observeEvent(gene(), once = TRUE, { set_by_app(gene()); updateSelectizeInput(session, "pg_gene", choices = gene_df, selected = gene(), server = TRUE,
    options = list(valueField = "value", labelField = "label", searchField = c("label", "aliases", "name"),
                   render = I("{ option: function(item, escape) { return '<div><b>' + escape(item.label) + '</b> <span class=\"opt-name\">' + escape(item.name) + '</span></div>'; } }"))) })
  last_main <- reactiveVal(NULL)
  observeEvent(input$nav, { if (!identical(input$nav, "pops")) return(); g <- tryCatch(gene(), error = function(e) NULL)
    if (!is.null(g) && !identical(g, last_main())) { last_main(g); if (!identical(input$pg_gene, g)) { set_by_app(g); updateSelectizeInput(session, "pg_gene", choices = gene_df, selected = g, server = TRUE) } } })
  pg <- reactive({ req(input$pg_gene, input$pg_gene %in% D$genes$symbol); input$pg_gene })
  cp <- reactive(cps[[st$cp]])
  parent <- reactive(if (st$level == "lineage") st$lin else pop_root(cp()))
  kids <- reactive(pop_kids(cp(), parent()))
  v <- reactive(pop_values(D, cp(), pg()))
  mv <- reactive(pop_map_values(D, cp(), pg()))
  output$pg_in_comp <- reactive(st$level != "atlas"); outputOptions(output, "pg_in_comp", suspendWhenHidden = FALSE)

  go <- function(level, lin = NULL) { st$level <- level; st$lin <- lin }
  observeEvent(input$pg_open, { id <- input$pg_open
    if (id %in% names(cps)) { st$cp <- id; go("compartment"); return() }
    if (st$level != "atlas" && pop_drillable(cp(), id)) go("lineage", id) else
      showNotification(paste0(pop_label(cp(), id), ": no ", cp()$words[["l2"]], " are defined"), type = "message", duration = 4) })
  observeEvent(event_data("plotly_click", source = "pg_click"), { e <- event_data("plotly_click", source = "pg_click"); id <- e$customdata
    if (is.null(id) || !length(id)) return()
    if (st$level == "atlas") { k <- POP_ATLAS[id]; if (!is.na(k) && k %in% names(cps)) { st$cp <- k; go("compartment") } else
      showNotification(paste0(id, ": cell populations for this cell type are planned for a later version"), type = "message", duration = 4); return() }
    if (st$level == "compartment" && pop_drillable(cp(), id)) go("lineage", id) })
  observeEvent(input$pg_to_atlas, go("atlas")); observeEvent(input$pg_to_comp, go("compartment"))

  output$pg_crumbs <- renderUI({
    sep <- span(class = "crumb-sep", HTML(" &rsaquo; "))
    div(class = "pop-crumbs",
        if (st$level == "atlas") span(class = "crumb-on", "Atlas") else actionLink("pg_to_atlas", "Atlas"),
        if (st$level != "atlas") tagList(sep, if (st$level == "compartment") span(class = "crumb-on", cp()$name) else actionLink("pg_to_comp", cp()$name)),
        if (st$level == "lineage") tagList(sep, span(class = "crumb-on", pop_label(cp(), st$lin))))
  })

  map_df <- reactive({
    if (st$level == "atlas") {
      u <- D$umap; req(u); lab <- as.character(u$meta$compartment); lab[grepl("^Ductal", lab)] <- "Ductal cells"
      avail <- names(POP_ATLAS)[POP_ATLAS %in% names(cps)]
      col <- ifelse(lab %in% avail, POP_CAT[1], "#d9d8d4")
      return(data.frame(x = u$coords$UMAP_1, y = u$coords$UMAP_2, id = lab, label = lab, colour = col,
                        hover = paste0(lab, if (length(avail)) ifelse(lab %in% avail, "<br>click to open its populations", "")), stringsAsFactors = FALSE))
    }
    m <- cp()$map
    if (st$level == "compartment") { w <- rep(TRUE, nrow(m)); x <- m$umap1_x; y <- m$umap1_y; id <- m$l1 } else {
      w <- m$l1 == st$lin & !is.na(m$umap2_x); x <- m$umap2_x; y <- m$umap2_y; id <- m$l2 }
    r <- pop_row(cp(), id[w])
    data.frame(x = x[w], y = y[w], id = id[w], label = r$label, colour = r$colour, row = which(w),
               hover = paste0(r$label, "<br>", m$state[w], " | sample ", m$sample[w]), stringsAsFactors = FALSE)
  })
  output$pg_map <- renderPlotly({ df <- map_df()
    ttl <- switch(st$level, atlas = "Atlas: click a highlighted cell type", compartment = paste0(cp()$name, ": ", cp()$words[["l1"]], " (click one to open)"),
                  lineage = paste0(pop_label(cp(), st$lin), ": ", cp()$words[["l2"]]))
    clickable <- if (st$level == "atlas") names(POP_ATLAS)[POP_ATLAS %in% names(cps)] else if (st$level == "compartment") Filter(function(i) pop_drillable(cp(), i), unique(df$id)) else character(0)
    pop_map_plotly(df, ttl, "pg_click", clickable) })

  score_choices <- reactive({ if (st$level == "atlas") return(NULL)
    ids <- if (st$level == "compartment") c(pop_kids(cp(), parent()), pop_kids(cp(), parent(), TRUE), unlist(lapply(pop_kids(cp(), parent()), function(k) if (!pop_drillable(cp(), k)) pop_kids(cp(), k, TRUE))))
           else c(pop_kids(cp(), parent()), pop_kids(cp(), parent(), TRUE))
    ids <- intersect(ids, colnames(cp()$scores)); stats::setNames(ids, paste0("Signature score: ", pop_label(cp(), ids))) })
  output$pg_colour_ui <- renderUI({ ch <- score_choices(); if (is.null(ch)) return(NULL)
    selectInput("pg_colour", NULL, choices = c(stats::setNames("gene", paste0("Expression of ", pg())), ch), width = "100%",
                selected = if (isTRUE(isolate(input$pg_colour) %in% c("gene", ch))) isolate(input$pg_colour) else "gene") })
  val_now <- reactive({ df <- map_df()
    if (st$level == "atlas") { val <- umap_values(D, gene_row(D, pg())$sym_sc); return(list(val = val, title = paste0(pg(), " expression"), cbar = "log1p\nCP10k")) }
    col <- if (is.null(input$pg_colour)) "gene" else input$pg_colour
    if (col != "gene" && col %in% colnames(cp()$scores)) return(list(val = cp()$scores[df$row, col], title = paste0("Signature score: ", pop_label(cp(), col)), cbar = "score"))
    val <- mv(); list(val = if (is.null(val)) NULL else val[df$row], title = paste0(pg(), " expression"), cbar = "log1p\nCP10k") })
  output$pg_map_val <- renderPlotly({ df <- map_df(); vn <- val_now()
    validate(need(!is.null(vn$val), paste(pg(), "is not in the single-cell data.")))
    ok <- !is.na(vn$val); pop_value_plotly(df[ok, ], vn$val[ok], vn$title, vn$cbar) })

  figs <- list(expr = reactive(plot_pop_expr(cp(), v(), kids(), pg(), input$pg_state)),
               comp = reactive(plot_pop_composition(cp(), parent(), kids())),
               states = reactive(plot_pop_states(cp(), v(), kids(), pg())),
               map = reactive({ df <- map_df(); plot_pop_map_static(df, paste0(if (st$level == "lineage") pop_label(cp(), st$lin) else cp()$name, ": populations"),
                                                               stats::setNames(df$colour, df$label)[!duplicated(df$label)]) }),
               map_val = reactive({ df <- map_df(); vn <- val_now(); if (is.null(vn$val)) return(empty_plot("Not available"))
                 ok <- !is.na(vn$val); plot_pop_map_static(df[ok, ], vn$title, value = vn$val[ok]) }))
  output$pg_gene_head <- renderText(paste0(pg(), " per sample"))
  output$pg_expr <- renderPlotly(to_plotly(figs$expr(), wrap = 70))
  output$pg_comp <- renderPlotly(to_plotly(figs$comp(), wrap = 170))
  output$pg_states <- renderPlotly(to_plotly(figs$states(), wrap = 100))

  output$pg_cards <- renderUI({ req(st$level != "atlas"); k <- kids(); so <- if (st$level == "lineage") pop_kids(cp(), parent(), TRUE) else
      c(pop_kids(cp(), parent(), TRUE), unlist(lapply(k, function(i) if (!pop_drillable(cp(), i)) pop_kids(cp(), i, TRUE))))
    div(class = "pop-cards", lapply(k, function(i) pop_card(cp(), i, drill = st$level == "compartment" && pop_drillable(cp(), i))),
        lapply(so, function(i) pop_card(cp(), i))) })
  output$pg_atlas_cards <- renderUI({ req(st$level == "atlas")
    div(class = "pop-cards", lapply(names(cps), function(k) { c1 <- cps[[k]]; r <- pop_row(c1, pop_root(c1))
      div(class = "pop-card", div(class = "pop-card-head", span(class = "pop-name", c1$name), span(class = "pop-badge badge-published", "Available")),
          p(class = "pop-def", r$definition), p(class = "pop-count", format(sum(c1$groups$n_cells[c1$groups$level == 1]), big.mark = ","), " cells; ",
                                                length(setdiff(pop_kids(c1, pop_root(c1)), pop_boundary(c1))), " ", c1$words[["l1"]]),
          tags$a(href = "#", class = "pop-open", onclick = sprintf("Shiny.setInputValue('pg_open', '%s', {priority: 'event'}); return false;", k), HTML(paste0("Open ", htmltools::htmlEscape(pop_lc(c1$name)), " &rsaquo;")))) }),
      div(class = "pop-card pop-planned", div(class = "pop-card-head", span(class = "pop-name", "Other cell types"), span(class = "pop-badge badge-unresolved", "Planned")),
          p(class = "pop-def", paste0(paste(setdiff(c("Fibroblasts", "Myeloid cells", "T and NK cells", "Tumour cells", "Endothelial cells and pericytes", "B and plasma cells"),
                                                     c(vapply(cps, `[[`, "", "name"), if ("myeloid" %in% names(cps)) "Myeloid cells")), collapse = ", "),
                                      " will follow, one cell type at a time.")))) })

  output$pg_atlas_cmp <- renderUI({ req(st$level != "atlas"); pop_atlas_compare_ui(cp(), kids()) })
  output$pg_decisions <- renderUI({ req(st$level != "atlas"); div(class = "pop-decisions", pop_decision_ui(cp())) })

  marker_tab <- reactive({ m <- cp()$markers; m <- m[m$population_id %in% kids(), ]; if (!nrow(m)) return(NULL)
    data.frame(population = pop_label(cp(), m$population_id), gene = m$gene, median_log_ratio = round(m$median_log_ratio, 3),
               samples_higher = paste0(m$samples_up, " / ", m$samples), median_frac = round(m$median_frac, 3), fdr = signif(m$q, 2), stringsAsFactors = FALSE) })
  output$pg_markers <- renderDT({ t <- marker_tab(); validate(need(!is.null(t), "No marker genes at this level (too few samples with enough cells)."))
    t$gene <- sprintf("<a href=\"#\" class=\"pg-gene\" data-gene=\"%s\">%s</a>", htmltools::htmlEscape(t$gene), htmltools::htmlEscape(t$gene))
    datatable(t, rownames = FALSE, escape = -2, style = "default", class = "compact stripe hover", extensions = "Buttons",
              colnames = c("population", "gene", "median log ratio", "samples higher", "median fraction expressing", "FDR"),
              callback = JS("table.on('click', 'a.pg-gene', function(e) { e.preventDefault(); Shiny.setInputValue('pg_goto', $(this).data('gene'), {priority: 'event'}); });"),
              options = list(pageLength = 10, dom = "Bfrtip", buttons = dt_buttons(cp()$name, "population_markers"))) }, server = FALSE)
  observeEvent(input$pg_goto, { g <- input$pg_goto; s <- D$genes$symbol[match(g, D$genes$sym_sc)]; if (is.na(s)) s <- g
    if (isTRUE(s %in% D$genes$symbol) && !identical(input$pg_gene, s)) { set_by_app(s); updateSelectizeInput(session, "pg_gene", choices = gene_df, selected = s, server = TRUE) } })

  pop_table <- reactive({ x <- v(); if (is.null(x)) return(data.frame(note = paste(pg(), "is not in the single-cell data")))
    x <- x[x$population_id %in% kids(), ]
    data.frame(population = x$label, label_source = POP_METHOD_TXT[pop_row(cp(), x$population_id)$label_method], sample = x$sample, study = x$study,
               study_ref = unname(STUDY_CITATION[x$study]), state = as.character(x$state), n_cells = x$n_cells, mean_cp10k = round(x$mean_cp10k, 3),
               avg_expr = round(x$avg_expr, 3), frac = round(x$frac, 4), stringsAsFactors = FALSE) })
  output$pg_table <- renderDT(datatable(pop_table(), rownames = FALSE, style = "default", class = "compact stripe hover", extensions = "Buttons",
                                        options = list(pageLength = 10, dom = "Bfrtip", buttons = dt_buttons(pg(), paste0(cp()$id, "_populations_per_sample")))), server = FALSE)

  fn <- function(ext) function() paste0(pg(), "_", cp()$id, "_", if (st$level == "lineage") st$lin else "lineages", "_", format(Sys.Date()), ".", ext)
  dl_plots <- function() list(figs$map(), figs$map_val(), figs$expr(), figs$comp(), figs$states())
  output$dl_pop_pdf <- downloadHandler(fn("pdf"), function(file) save_figs(file, dl_plots(), 10, 7, "pdf"))
  output$dl_pop_png <- downloadHandler(fn("png"), function(file) save_figs(file, dl_plots(), 10, 7, "png", 1))
  output$dl_pop_csv <- downloadHandler(fn("csv"), function(file) {
    t <- pop_table(); r <- cp()$tax[cp()$tax$population_id %in% kids(), ]
    defs <- data.frame(population = r$label, label_source = unname(POP_METHOD_TXT[r$label_method]), definition = r$definition, signature_genes = r$signature_genes,
                       signature_source = r$signature_ref,
                       notes = vapply(r$population_id, function(i) paste(pop_decisions_for(cp(), i)$decision_id, collapse = "; "), ""),
                       references = vapply(seq_len(nrow(r)), function(k) {
                         ids <- unique(c(pop_ref_ids(r[k, ]), pop_split(pop_decisions_for(cp(), r$population_id[k])$refs)))
                         rr <- cp()$refs[match(ids, cp()$refs$ref_id), ]; paste(paste0(rr$citation, " doi:", rr$doi), collapse = " | ") }, ""), row.names = NULL)
    dec <- pop_decision_table(cp())
    con <- file(file, "w"); on.exit(close(con))
    writeLines(paste0("# ", pg(), " in ", cp()$name, if (st$level == "lineage") paste0(" > ", pop_label(cp(), st$lin)),
                      "; per-sample values, then population definitions", if (!is.null(dec)) ", then notes on how the groups were named"), con)
    utils::write.csv(t, con, row.names = FALSE); writeLines("", con); utils::write.csv(defs, con, row.names = FALSE)
    if (!is.null(dec)) { writeLines("", con); utils::write.csv(dec, con, row.names = FALSE) } })
}
