# Cell-population layer: a small synthetic pop_<id>.rds (the shape 07_cell_populations.R writes) goes through
# build_pop() (05) and back through the app's accessors (R/08_populations.R); values must come back unchanged.
root <- normalizePath(file.path(test_path(), "..", ".."))
for (f in sort(list.files(file.path(root, "app", "R"), full.names = TRUE))) source(f)
suppressPackageStartupMessages({ library(Matrix); library(shiny) })
old_opt <- options(pdacx.build_store.autorun = FALSE, pdacx.root = root)
source(file.path(root, "data-raw", "05_build_store.R"))
options(old_opt)

make_pop <- function(dir) {
  set.seed(7)
  tax <- data.frame(population_id = c("TST", "TST_A", "TST_B", "A_1", "A_2", "A_SCORE"), compartment = "Test", level = c(0, 1, 1, 2, 2, 2),
                    parent_id = c("", "TST", "TST", "TST_A", "TST_A", "TST_A"), label = c("Test cells", "Lineage A", "Lineage B", "A one", "A two", "A programme"),
                    definition = paste("Definition of", c("TST", "TST_A", "TST_B", "A_1", "A_2", "A_SCORE")), signature_genes = c("", "G1;G2;G3", "G4;G5;G6", "G1;G7;G8", "G2;G9;G10", "G3;G4;G5"),
                    signature_ref = c("", "ref2020 (Fig. 1)", "ref2020 (Fig. 2)", "ref2021 (Table 1)", "ref2021 (Table 1)", "ref2021 (Table 2)"),
                    definition_refs = c("ref2020", "ref2020", "ref2020", "ref2021;ref2020", "ref2021", "ref2021"), pdac_evidence = "", label_method = c("published", "published", "published", "signature", "signature", "score_only"),
                    status = "verified", also_known_as = "", stringsAsFactors = FALSE)
  refs <- data.frame(ref_id = c("ref2020", "ref2021"), citation = c("Author A. Title one. J 2020", "Author B. Title two. J 2021"), doi = c("10.1/one", "10.1/two"),
                     checked_against = "full text", stringsAsFactors = FALSE)
  samples <- paste0("s", 1:6); state <- c("Primary tumor", "Primary tumor", "Primary tumor", "Donor", "Adjacent normal", "Primary tumor")
  g1 <- expand.grid(sample = samples, population_id = c("TST_A", "TST_B"), stringsAsFactors = FALSE); g1$level <- 1
  g2 <- expand.grid(sample = samples, population_id = c("A_1", "A_2", "TST_A_UNRESOLVED"), stringsAsFactors = FALSE); g2$level <- 2
  groups <- rbind(g1, g2); groups$n_cells <- sample(10:90, nrow(groups), TRUE)
  for (s in samples) groups$n_cells[groups$level == 1 & groups$sample == s & groups$population_id == "TST_A"] <- sum(groups$n_cells[groups$level == 2 & groups$sample == s])
  groups$group_id <- seq_len(nrow(groups)); groups$study <- "GSE1"; groups$disease_state <- state[match(groups$sample, samples)]; groups$met_site <- NA
  groups$treatment <- ifelse(groups$disease_state == "Donor", "Donor (untreated)", "Treatment-naive")
  genes <- paste0("G", 1:12); ng <- nrow(groups)
  expr <- list(genes = genes, mean_cp10k_x100 = matrix(sample(0:50000, length(genes) * ng, TRUE), length(genes), dimnames = list(genes, NULL)),
               frac_x10000 = matrix(sample(0:10000, length(genes) * ng, TRUE), length(genes), dimnames = list(genes, NULL)), note = "synthetic")
  n <- 300; l1 <- sample(c("TST_A", "TST_B"), n, TRUE); l2 <- ifelse(l1 == "TST_A", sample(c("A_1", "A_2", "TST_A_UNRESOLVED"), n, TRUE), l1)
  meta <- data.frame(sample = sample(samples, n, TRUE), study = "GSE1", disease_state = "Primary tumor", met_site = NA, treatment = "Treatment-naive", l1 = l1, l2 = l2,
                     cluster = ifelse(l1 == "TST_A", paste0("TST_A_c", match(l2, c("A_1", "A_2", "TST_A_UNRESOLVED"))), NA), umap1_x = rnorm(n), umap1_y = rnorm(n),
                     umap2_x = ifelse(l1 == "TST_A", rnorm(n), NA), umap2_y = ifelse(l1 == "TST_A", rnorm(n), NA), stringsAsFactors = FALSE)
  X <- rsparsematrix(n, length(genes), 0.2, rand.x = function(k) round(runif(k, 0, 5), 2)); colnames(X) <- genes; X <- as(X, "CsparseMatrix")
  scores <- matrix(round(rnorm(n * 5), 3), n, 5, dimnames = list(NULL, c("TST_A", "TST_B", "A_1", "A_2", "A_SCORE")))
  calls <- data.frame(lineage = "TST_A", cluster = paste0("TST_A_c", 1:3), n_cells = c(60, 50, 40), n_samples = 6, median_genes = 1000,
                      call = c("A_1", "A_2", "TST_A_UNRESOLVED"), reason = c("", "", "two signatures score almost equally"), note = c("", "strongest programme: A programme (score only) (score 0.90)", ""),
                      best = c("A_1", "A_2", "A_1"), best_score = c(0.8, 0.6, 0.3), runner_up = c("A_2", "A_1", "A_2"), margin = c(0.5, 0.3, 0.01), stringsAsFactors = FALSE)
  x <- list(compartment = "Test", taxonomy = tax, references = refs, groups = groups[, c("level", "sample", "population_id", "n_cells", "group_id", "study", "disease_state", "met_site", "treatment")],
            expr = expr, label_check = list(counts = NULL, mean_scores = NULL, agreement = c(TST_A = 0.9, TST_B = 0.8)), clusters = list(calls = calls, scores = list()),
            markers = data.frame(population_id = "A_1", gene = "G7", median_log_ratio = 1.2, samples_up = 5, samples = 6, median_frac = 0.4, p = 0.01, q = 0.03),
            map = list(meta = meta, scores = scores, expr = X, note = "synthetic"), params = list(), built = "now")
  saveRDS(x, file.path(dir, "pop_test.rds")); x
}

# a decision record as data-raw/reference/decisions.csv; ref2022 is only in the current references.csv, not in the 07 output
make_decisions <- function(dir, rows = NULL) {
  d <- data.frame(decision_id = c("ALL_01", "TST_01", "TST_02", "OTHER_01"), compartment = c("all", "Test", "Test", "Other"),
                  applies_to = c("", "A_SCORE", "", "X_1"), kind = c("rule", "score_only", "not_labelled", "signature"),
                  question = c("How are clusters named?", "Why is A programme only a score?", "Why no C cells?", "Other compartment"),
                  decision = c("Highest module score.", "Shown as a score.", "Not used.", "x"), reason = c("Rule text.", "Handling artefact.", "Single study.", "x"),
                  refs = c("", "ref2021", "ref2022", "ref2020"), stringsAsFactors = FALSE)
  if (!is.null(rows)) d <- rows(d)
  utils::write.csv(d, file.path(dir, "decisions.csv"), row.names = FALSE)
  utils::write.csv(data.frame(ref_id = c("ref2020", "ref2021", "ref2022"), citation = c("Author A. Title one. J 2020", "Author B. Title two. J 2021", "Author C. Title three. J 2022"),
                              doi = c("10.1/one", "10.1/two", "10.1/three"), checked_against = "full text"), file.path(dir, "references.csv"), row.names = FALSE)
  dir
}

test_that("population store round-trips every gene exactly and the app reads it", {
  dir <- tempfile("pop"); dir.create(dir); src <- make_pop(dir)
  build_pop(dir, function(...) NULL, file.path(dir, "pop_test.rds"), ref_dir = file.path(dir, "no_reference_dir"))
  cps <- pop_load(dir); expect_named(cps, "test"); cp <- cps$test
  d <- list(genes = data.frame(symbol = c(src$expr$genes, "NOPE"), sym_sc = c(src$expr$genes, NA), stringsAsFactors = FALSE))
  for (g in src$expr$genes) {
    v <- pop_values(d, cp, g)
    expect_equal(v$mean_cp10k, unname(src$expr$mean_cp10k_x100[g, ]) / 100)
    expect_equal(v$frac, unname(src$expr$frac_x10000[g, ]) / 1e4)
    expect_equal(pop_map_values(d, cp, g), unname(as.numeric(src$map$expr[, g])))
  }
  expect_null(pop_values(d, cp, "NOPE"))
  # unresolved cells get a pseudo-population under their lineage, grey, with the cluster reasons on its card
  expect_true("TST_A_UNRESOLVED" %in% cp$tax$population_id)
  expect_equal(pop_row(cp, "TST_A_UNRESOLVED")$parent_id, "TST_A"); expect_equal(pop_row(cp, "TST_A_UNRESOLVED")$colour, POP_GREY)
  expect_equal(pop_kids(cp, "TST"), c("TST_A", "TST_B"))
  expect_equal(pop_kids(cp, "TST_A"), c("A_1", "A_2", "TST_A_UNRESOLVED"))
  expect_equal(pop_kids(cp, "TST_A", with_scores = TRUE), "A_SCORE")
  expect_true(pop_drillable(cp, "TST_A")); expect_false(pop_drillable(cp, "TST_B"))
  # composition: shares of a sample's parent cells sum to 1
  co <- pop_composition(cp, "TST_A", pop_kids(cp, "TST_A")); s <- tapply(co$share, co$sample, sum)
  expect_true(all(abs(s - 1) < 1e-12))
  # every population card shows its references with DOI links
  html <- as.character(pop_card(cp, "A_1")); expect_match(html, "doi.org/10.1/two", fixed = TRUE); expect_match(html, "doi.org/10.1/one", fixed = TRUE)
  expect_match(as.character(pop_card(cp, "TST_A_UNRESOLVED")), "Group 3 (40 cells): two signatures fitted almost equally well. Closest: A one, then A two.", fixed = TRUE)
  expect_match(as.character(pop_card(cp, "A_2")), "Some of these cells also show the A programme.", fixed = TRUE)
  v <- pop_values(d, cp, "G1")
  for (p in suppressWarnings(list(plot_pop_expr(cp, v, pop_kids(cp, "TST_A"), "G1", "Primary tumour"), plot_pop_states(cp, v, pop_kids(cp, "TST_A"), "G1"),
                 plot_pop_composition(cp, "TST_A", pop_kids(cp, "TST_A"))))) expect_s3_class(p, "ggplot")
  for (k in cps) store_close(k$st)
})

test_that("build_pop refuses populations without a definition", {
  dir <- tempfile("pop"); dir.create(dir); x <- make_pop(dir)
  x$map$meta$l2[1] <- "NOT_DEFINED"; saveRDS(x, file.path(dir, "pop_test.rds"))
  expect_error(build_pop(dir, function(...) NULL, file.path(dir, "pop_test.rds")), "without a definition")
})

test_that("the decision record is attached at build time and shown on cards, in its section and in the CSV", {
  dir <- tempfile("pop"); dir.create(dir); make_pop(dir); rd <- make_decisions(file.path(dir, "ref") |> (\(p) { dir.create(p); p })())
  build_pop(dir, function(...) NULL, file.path(dir, "pop_test.rds"), ref_dir = rd)
  cp <- pop_load(dir)$test
  expect_equal(cp$decisions$decision_id, c("ALL_01", "TST_01", "TST_02"))           # own compartment + "all", not other compartments
  expect_true("ref2022" %in% cp$refs$ref_id)                                         # cited reference taken from the current references.csv
  card <- as.character(pop_card(cp, "A_SCORE"))
  expect_match(card, "Why is A programme only a score?", fixed = TRUE); expect_match(card, "Handling artefact.", fixed = TRUE)
  expect_match(card, "doi.org/10.1/two", fixed = TRUE)
  expect_false(grepl("Why no C cells?", as.character(pop_card(cp, "A_1")), fixed = TRUE))   # decisions only on the cards they name
  expect_match(card, "<summary>Why is A programme only a score?</summary>\n    <div class=\"pop-why-body\">", fixed = TRUE)   # collapsed to the headline
  ui <- as.character(pop_decision_ui(cp))
  expect_false(grepl("Decision:", ui, fixed = TRUE)); expect_false(grepl("Why:", ui, fixed = TRUE))
  for (q in c("How are clusters named?", "Why is A programme only a score?", "Why no C cells?", "Not used", "For every cell type", "doi.org/10.1/three"))
    expect_match(ui, q, fixed = TRUE)
  expect_false(grepl("Other compartment", ui, fixed = TRUE))
  expect_lt(regexpr("Why no C cells?", ui, fixed = TRUE), regexpr("How are clusters named?", ui, fixed = TRUE))   # own decisions before shared rules
  tab <- pop_decision_table(cp); expect_equal(nrow(tab), 3); expect_match(tab$references[tab$note == "TST_02"], "doi:10.1/three", fixed = TRUE)
  expect_equal(tab$applies_to[tab$note == "TST_01"], "A programme")
  store_close(cp$st)
})

test_that("a store without a decision record still works", {
  dir <- tempfile("pop"); dir.create(dir); make_pop(dir)
  build_pop(dir, function(...) NULL, file.path(dir, "pop_test.rds"), ref_dir = file.path(dir, "missing"))
  cp <- pop_load(dir)$test; expect_equal(nrow(cp$decisions), 0)
  expect_match(as.character(pop_decision_ui(cp)), "not part of this data version", fixed = TRUE)
  expect_s3_class(pop_card(cp, "A_1"), "shiny.tag"); expect_null(pop_decision_table(cp)); store_close(cp$st)
})

test_that("build_pop refuses a decision record with unknown populations, missing references or no references", {
  dir <- tempfile("pop"); dir.create(dir); make_pop(dir); f <- file.path(dir, "pop_test.rds")
  bad <- function(edit, msg) { rd <- file.path(dir, paste0("ref", sample.int(1e6, 1))); dir.create(rd); make_decisions(rd, edit)
    expect_error(build_pop(dir, function(...) NULL, f, ref_dir = rd), msg) }
  bad(function(d) { d$applies_to[2] <- "NOPE"; d }, "missing from the Test taxonomy")
  bad(function(d) { d$refs[2] <- "ref1999"; d }, "missing from references.csv")
  bad(function(d) { d$refs[3] <- ""; d }, "without references")
  bad(function(d) { d$kind[2] <- "maybe"; d }, "unknown kind")
  bad(function(d) rbind(d, d[2, ]), "duplicated decision_id")
})

test_that("the shipped decision record names existing populations and references", {
  ref <- file.path(root, "data-raw", "reference"); skip_if_not(file.exists(file.path(ref, "decisions.csv")))
  rd <- function(f) utils::read.csv(file.path(ref, f), stringsAsFactors = FALSE, colClasses = "character", encoding = "UTF-8")
  d <- rd("decisions.csv"); pops <- rd("populations.csv"); refs <- rd("references.csv")
  expect_false(anyDuplicated(d$decision_id) > 0)
  expect_true(all(d$kind %in% POP_DECISION_KINDS))
  expect_true(all(d$compartment %in% c("all", pops$compartment)))
  sp <- function(x) { v <- trimws(unlist(strsplit(x, ";", fixed = TRUE))); v[nzchar(v)] }
  expect_true(all(sp(d$applies_to) %in% pops$population_id)); expect_true(all(sp(d$refs) %in% refs$ref_id))
  cited <- unique(c(sp(pops$definition_refs), unlist(regmatches(pops$signature_ref, gregexpr("[a-z]+[0-9]{4}", pops$signature_ref)))))
  expect_true(all(cited %in% refs$ref_id))
  expect_true(all(nzchar(sp(d$refs[!d$kind %in% c("rule", "structure")]))))
})

test_that("level 1 assigned by signature: unresolved at level 1, boundary groups out of the make-up, atlas labels as comparison", {
  dir <- tempfile("pop"); dir.create(dir); x <- make_pop(dir)
  x$taxonomy$label_method[x$taxonomy$population_id == "TST_A"] <- "signature"
  x$taxonomy$label_method[x$taxonomy$population_id == "TST_B"] <- "boundary"
  x$level1_method <- "signature"
  # some level-1 cells unresolved
  g <- x$groups; u <- g[g$level == 1 & g$population_id == "TST_B", ][1:2, ]; u$population_id <- "TST_UNRESOLVED"; u$group_id <- nrow(g) + 1:2
  x$groups <- rbind(g, u); x$expr$mean_cp10k_x100 <- cbind(x$expr$mean_cp10k_x100, x$expr$mean_cp10k_x100[, 1:2]); x$expr$frac_x10000 <- cbind(x$expr$frac_x10000, x$expr$frac_x10000[, 1:2])
  x$clusters$calls <- rbind(x$clusters$calls, data.frame(lineage = "TST", cluster = "TST_c9", n_cells = 40, n_samples = 2, median_genes = 500, call = "TST_UNRESOLVED",
    reason = "no signature reaches the minimum score", note = "", best = "TST_A", best_score = 0.05, runner_up = "TST_B", margin = 0.01))
  x$atlas_compare <- data.frame(level = c(1L, 1L, 2L), population_id = c("TST_A", "TST_A", "A_1"), atlas_label = c("myFb 1", "iFb 2", "myFb 1"), n_cells = c(30L, 10L, 5L))
  saveRDS(x, file.path(dir, "pop_test.rds")); build_pop(dir, function(...) NULL, file.path(dir, "pop_test.rds"), ref_dir = file.path(dir, "none"))
  cp <- pop_load(dir)$test
  expect_equal(cp$level1_method, "signature"); expect_equal(cp$words[["l1"]], "groups")
  r <- pop_row(cp, "TST_UNRESOLVED"); expect_equal(r$level, 1L); expect_equal(r$parent_id, "TST")
  expect_true("TST_UNRESOLVED" %in% pop_kids(cp, "TST"))
  co <- pop_composition(cp, "TST", pop_kids(cp, "TST")); expect_false("TST_B" %in% co$population_id)   # boundary not counted
  s <- tapply(co$share, co$sample, sum); expect_true(all(abs(s - 1) < 1e-12))
  expect_false(grepl("myFb", as.character(pop_card(cp, "TST_A")), fixed = TRUE))              # no atlas labels on the cards
  POP_ATLAS_NOTE[["Test"]] <<- "Explanation of the atlas publication's clusters."
  cmp <- as.character(pop_atlas_compare_ui(cp, pop_kids(cp, "TST")))                              # ... but in an explained comparison table
  expect_match(cmp, "Explanation of the atlas publication", fixed = TRUE); expect_match(cmp, "<td>75%</td>", fixed = TRUE); expect_match(cmp, "myFb 1", fixed = TRUE)
  expect_match(as.character(pop_card(cp, "TST_B")), "Other cell type", fixed = TRUE)
  un <- as.character(pop_card(cp, "TST_UNRESOLVED"))                                               # plain language for groups that were not named
  expect_match(un, "Group 9 (40 cells): no signature fitted well enough. Closest: lineage A, then lineage B.", fixed = TRUE); expect_false(grepl("TST_c9", un, fixed = TRUE))
  expect_match(as.character(pop_caption(cp)), "groups and sub-groups named here from published signatures", fixed = TRUE)
  store_close(cp$st)
})

test_that("the wording of populations comes from the current populations.csv; changed signatures keep the 07 version", {
  dir <- tempfile("pop"); dir.create(dir); x <- make_pop(dir); rd <- file.path(dir, "ref"); dir.create(rd)
  cur <- x$taxonomy; cur$definition[cur$population_id == "A_1"] <- "New wording of A one."; cur$label[cur$population_id == "A_1"] <- "A one, renamed"
  cur$signature_genes[cur$population_id == "A_2"] <- "G99"; cur$definition[cur$population_id == "A_2"] <- "Must not appear."
  utils::write.csv(cur, file.path(rd, "populations.csv"), row.names = FALSE)
  msg <- character(0); build_pop(dir, function(...) msg <<- c(msg, paste0(...)), file.path(dir, "pop_test.rds"), ref_dir = rd)
  cp <- pop_load(dir)$test
  expect_equal(pop_label(cp, "A_1"), "A one, renamed"); expect_equal(pop_row(cp, "A_1")$definition, "New wording of A one.")
  expect_equal(pop_row(cp, "A_2")$definition, x$taxonomy$definition[x$taxonomy$population_id == "A_2"])
  expect_true(any(grepl("signature or method changed for A_2", msg, fixed = TRUE)))
  store_close(cp$st)
})

test_that("notes of groups read as one sentence", {
  expect_equal(pop_and(c("a", "b", "c")), "a, b and c"); expect_equal(pop_and("a"), "a")
  expect_equal(pop_note_text(c("a strong C programme, csX", "a strong M programme")), "show the C (csX) and M programmes")
  expect_equal(pop_note_text(c("a strong X programme", "a strong Y programme", "looks more like Z cells")),
               c("show the X and Y programmes", "look more like Z cells"))
})
