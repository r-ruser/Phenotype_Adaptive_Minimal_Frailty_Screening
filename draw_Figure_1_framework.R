if(.Platform$OS.type=="windows") invisible(try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE))
# Reproducible R-only schematic for the five-cohort FI item-reduction study.
# Run: Rscript scripts/draw_Figure_1_framework.R [revision_root]
options(stringsAsFactors = FALSE, warn = 1)
suppressWarnings(try(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"), silent = TRUE))
required <- c("grid", "svglite", "ragg", "systemfonts")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))

library(grid)
options(stringsAsFactors=FALSE)
ms5_args <- commandArgs(FALSE)
ms5_script <- sub("^--file=", "", ms5_args[grepl("^--file=",ms5_args)][1])
code_dir <- Sys.getenv("MS5_CODE_DIR", unset=if(!is.na(ms5_script)) dirname(normalizePath(ms5_script,winslash="/",mustWork=TRUE)) else getwd())
source(file.path(code_dir,"00_config.R"),encoding="UTF-8")

root <- publication_dir
output_dir <- file.path(root,"figures"); dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
source_dir <- file.path(root,"figure_source"); dir.create(source_dir,recursive=TRUE,showWarnings=FALSE)
# Constants and node specification record the audited fixed-wave analytic framework.
width_mm <- 183
height_mm <- 145
family <- "Arial"
palette <- c(ink = "#203A47", body = "#3D515D", border = "#78929D",
             arrow = "#5B727C", blue = "#276B82", blue_light = "#EDF4F7",
             gold = "#F4DE9A", gold_light = "#FFF7DC", purple = "#EEE7F3",
             purple_ink = "#655575", green_light = "#EDF4EF", neutral = "#F3F5F5")
stem <- file.path(output_dir, "Figure_1")

nodes <- data.frame(
  id = c("cohorts", "inputs", "outer_split", "training", "universal", "phenotype", "lock", "scores",
         "heldout", "apply", "outcomes", "comparisons", "summary"),
  x = c(.1425, .1425, .1425, .49, .394, .586, .49, .49, .8475, .8475, .8475, .8475, .50),
  y = c(.707, .504, .293, .727, .556, .556, .412, .257, .727, .575, .386, .208, .081),
  width = c(.207, .207, .207, .366, .164, .164, .366, .366, .218, .218, .218, .218, .95),
  height = c(.178, .158, .174, .114, .176, .176, .070, .159, .114, .102, .234, .096, .080),
  label = c(
    "Five ageing cohorts | HRS, ELSA, SHARE, CHARLS, MHAS | N = 34,954",
    "Baseline inputs | Complete FI: 39-50 deficits, FI >= 0.25 | Shared bank: 16 questions",
    "Leave-one-cohort-out | 4 training cohorts + 1 held-out cohort | Repeat across 5 outer folds",
    "Four-cohort training only | 3 inner folds; six-item FI-scale mapping",
    "Universal ranking (primary) | XGBoost for complete FI | Cross-fitted TreeSHAP | Rank 16 questions",
    "Phenotype ranking (secondary) | Training-only latent classes | Within-class TreeSHAP",
    "Lock rankings, cutoffs and six-item mapping",
    "Unweighted top-k scores | Mean of observed selected deficits | k = 1-16 | m >= ceiling(0.8k)",
    "One held-out cohort | Untouched during learning",
    "Apply locked scores | Assign trained classes for phenotype scores",
    "Concurrent FI approximation | Follow-up FI frailty | Incident ADL limitation",
    "Matched comparisons | Complete baseline FI | Same outcome-specific participants",
    "Exploratory length–performance trade-off | AUROC, correlation, errors, prospective outcomes, random benchmarks"
  ))
edges <- data.frame(
  from = c("cohorts", "inputs", "outer_split", "outer_split", "training", "training", "universal",
           "phenotype", "lock", "scores", "heldout", "apply", "outcomes", "comparisons"),
  to = c("inputs", "outer_split", "training", "heldout", "universal", "phenotype", "lock",
         "lock", "scores", "apply", "apply", "comparisons", "comparisons", "summary"),
  role = c("baseline inputs", "outer validation split", "four training cohorts", "one held-out cohort",
           "primary ranking", "secondary ranking", "lock universal choices", "lock phenotype choices",
           "construct scores", "transfer locked rules", "evaluate held-out data", "three outcome types",
           "common risk sets", "summarize held-out curves"))
edges$to[edges$from == "apply"] <- "outcomes"

nodes <- rbind(nodes, data.frame(id=c("training_uncertainty","heterogeneity_sources"),x=c(.26,.74),y=c(.115,.115),width=.46,height=.165,label=c("Primary pipeline: 500 coupled bootstrap repetitions; ranks, mapping and cutoffs relearned; training-only and joint evaluation", "FI and candidate distributions; disease/function profiles; relative error improvement; FI-stratified residuals; MHAS reference decomposition")))
edges <- rbind(edges,data.frame(from="summary",to=c("training_uncertainty","heterogeneity_sources"),role=c("development uncertainty","descriptive sources of heterogeneity")))
nodes$label[nodes$id == "cohorts"] <- "Five ageing cohorts; N=34,954; age >=65 and at least one of six specified chronic diseases"
nodes$y[nodes$id == "summary"] <- .075
nodes$height[nodes$id == "summary"] <- .110
main_nodes <- !nodes$id %in% c("training_uncertainty","heterogeneity_sources")
nodes$y[main_nodes] <- .22 + (.78/.9)*nodes$y[main_nodes]
nodes$height[main_nodes] <- (.78/.9)*nodes$height[main_nodes]
qa_text <- list()
txt <- function(label, x, y, size = 7.2, face = "plain", colour = palette["body"],
                just = "centre", max_width = NULL, lineheight = 1.18) {
  prop <- gpar(fontfamily = family, fontsize = size, fontface = face, col = colour, lineheight = lineheight)
  grid.text(label, x = x, y = y, just = just, gp = prop)
  if (!is.null(max_width)) {
    lines <- strsplit(label, "\n", fixed = TRUE)[[1]]
    widths <- vapply(lines, function(z) convertWidth(grobWidth(textGrob(z, gp = prop)), "npc", valueOnly = TRUE), numeric(1))
    qa_text[[length(qa_text) + 1L]] <<- data.frame(label = label, size_pt = size,
      width_npc = max(widths), available_width_npc = max_width,
      fits = max(widths) <= max_width)
  }
}
rect <- function(x, y, w, h, fill = "white", border = palette["border"],
                 dashed = FALSE, radius_mm = 1.4, line_width = .75) {
  grid.roundrect(x = x, y = y, width = w, height = h, r = unit(radius_mm, "mm"),
    gp = gpar(fill = fill, col = border, lwd = line_width, lty = if (dashed) "22" else "solid"))
}
connector <- function(x, y, colour = palette["arrow"], dashed = FALSE, arrow = TRUE, width = .9) {
  grid.lines(x = x, y = y,
    arrow = if (arrow) grid::arrow(length = unit(1.45, "mm"), type = "closed") else NULL,
    gp = gpar(col = colour, fill = colour, lwd = width, lty = if (dashed) "33" else "solid",
              lineend = "round", linejoin = "round"))
}
panel <- function(x0, x1, title, letter, fill) {
  grid.rect(x = (x0 + x1)/2, y = .51, width = x1 - x0, height = .730,
    gp = gpar(fill = "white", col = palette["border"], lwd = .8, lty = "22"))
  rect((x0 + x1)/2, .839, x1 - x0 - .018, .057, fill = fill,
    border = if (fill == palette["gold"]) "#CBB574" else fill, radius_mm = 1.05)
  txt(letter, x0 + .023, .84, 8.5, "bold",
    colour = if (fill == palette["gold"]) palette["ink"] else "white")
  txt(title, (x0 + x1)/2 + .007, .84, 7.9, "bold",
    colour = if (fill == palette["gold"]) palette["ink"] else "white",
    max_width = x1 - x0 - .064)
}
person <- function(x, y, colour, scale = 1) {
  grid.circle(x, y + .013 * scale, r = unit(.85 * scale, "mm"), gp = gpar(fill = colour, col = NA))
  grid.roundrect(x, y - .008 * scale, width = .012 * scale, height = .026 * scale,
    r = unit(.6, "mm"), gp = gpar(fill = colour, col = NA))
}

draw_figure <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  # Remove the two in-figure title lines and trim their 12 mm of headroom.
  # Keep all scientific modules at their original physical sizes.
  pushViewport(viewport(x = 0, y = .22, width = 1, height = .78/.9,
                        just = c("left", "bottom")))
  panel(.025, .260, "Study inputs", "a", palette["blue"])
  panel(.290, .690, "Training & scoring", "b", palette["gold"])
  panel(.720, .975, "Held-out evaluation", "c", palette["blue"])

  # Connectors are drawn before cards, keeping all arrowheads outside card interiors.
  connector(c(.1425, .1425), c(.618, .583))
  connector(c(.1425, .1425), c(.425, .380))
  connector(c(.246, .277, .277, .307), c(.320, .320, .727, .727))
  # Separate data path for the held-out cohort, with no return to training.
  connector(c(.277, .277, .8475, .8475),
            c(.727, .796, .796, .784),
            colour = palette["blue"], dashed = TRUE, width = .75)
  connector(c(.394, .394), c(.670, .644))
  connector(c(.586, .586), c(.670, .644), colour = palette["purple_ink"])
  connector(c(.394, .394), c(.468, .447))
  connector(c(.586, .586), c(.468, .447), colour = palette["purple_ink"])
  connector(c(.49, .49), c(.377, .3365))
  connector(c(.673, .705, .705, .7385), c(.257, .257, .575, .575))
  connector(c(.8475, .8475), c(.670, .626))
  connector(c(.8475, .8475), c(.524, .503))
  connector(c(.8475, .8475), c(.269, .256))
  connector(c(.8475, .8475), c(.160, .131))

  rect(.1425, .707, .207, .178, palette["blue_light"])
  txt("Five ageing cohorts", .1425, .766, 8.3, "bold", palette["ink"], max_width = .185)
  
  txt("HRS · ELSA · SHARE\nCHARLS · MHAS", .1425, .724, 7.35, max_width = .182)
  txt("N = 34,954", .1425, .679, 7.75, "bold", palette["ink"], max_width = .182)

  txt("Age ≥65; ≥1 of six\nspecified chronic diseases", .1425, .649, 6.65, max_width = .185)

  rect(.1425, .504, .207, .158, palette["neutral"])
  txt("Baseline inputs", .1425, .556, 8.0, "bold", palette["ink"], max_width = .185)
  txt("Complete FI: 39–50 deficits\nFrailty: FI ≥ 0.25", .1425, .516, 7.05, max_width = .183)
  txt("Shared bank: 16 questions", .1425, .458, 7.05, "bold", max_width = .183)

  rect(.1425, .293, .207, .174, palette["blue_light"])
  txt("Leave-one-cohort-out", .1425, .346, 7.9, "bold", palette["ink"], max_width = .185)
  for (i in 1:5) {
    rect(.077 + (i - 1)*.033, .298, .023, .028,
      fill = if (i == 5) "#F1DEA8" else "#D2E4EC", border = if (i == 5) "#B49B5C" else "#759CAC", radius_mm = .4)
  }
  txt("4 training + 1 held-out", .1425, .257, 7.1, "bold", max_width = .186)
  txt("Repeat across 5 outer folds", .1425, .220, 6.85, max_width = .186)

  rect(.49, .727, .366, .114, palette["gold_light"], border = "#BDAC78")
  txt("Four-cohort training only", .49, .748, 8.7, "bold", palette["ink"], max_width = .338)
  txt("3 inner folds; six-item FI-scale mapping", .49, .704, 7.2, max_width = .338)

  rect(.394, .556, .164, .176, palette["gold_light"], border = "#BDAC78")
  txt("Universal ranking\n(primary)", .394, .610, 7.65, "bold", palette["ink"], max_width = .146)
  txt("XGBoost for FI\nCross-fitted TreeSHAP\nRank 16 questions", .394, .530, 6.9, max_width = .146)

  rect(.586, .556, .164, .176, palette["purple"], border = "#9A89A5", dashed = TRUE)
  txt("Phenotype ranking\n(secondary)", .586, .610, 7.65, "bold", palette["purple_ink"], max_width = .146)
  txt("Training-only classes\nWithin-class TreeSHAP\nRank 16 questions", .586, .530, 6.9, max_width = .146)

  rect(.49, .412, .366, .070, palette["gold"], border = "#BDAC78", radius_mm = 1)
  txt("Lock ranks, cutoffs & six-item mapping", .49, .413, 7.65, "bold", palette["ink"], max_width = .340)

  rect(.49, .257, .366, .159, palette["blue_light"])
  txt("Unweighted top-k scores", .49, .309, 8.2, "bold", palette["ink"], max_width = .338)
  txt("Score = mean of observed selected deficits", .49, .269, 7.2, max_width = .338)
  txt("k = 1–16;  m ≥ ceiling(0.8k)", .49, .228, 7.3, "bold", max_width = .338)
  txt("m = number of observed selected questions", .49, .196, 6.6, max_width = .338)

  rect(.8475, .727, .218, .114, palette["blue_light"])
  txt("One held-out cohort", .8475, .749, 8.0, "bold", palette["ink"], max_width = .196)
  txt("Untouched during learning", .8475, .706, 6.95, max_width = .196)

  rect(.8475, .575, .218, .102, palette["neutral"])
  txt("Apply locked scores", .8475, .600, 7.9, "bold", palette["ink"], max_width = .196)
  txt("Trained class assignment\nfor phenotype scores", .8475, .555, 6.7, max_width = .196)

  rect(.8475, .386, .218, .234, palette["green_light"], border = "#82998A")
  txt("Evaluation outcomes", .8475, .477, 7.7, "bold", palette["ink"], max_width = .196)
  txt("Concurrent FI approximation", .8475, .434, 6.95, "bold", max_width = .196)
  txt("AUROC + correlation + errors", .8475, .408, 6.6, max_width = .196)
  grid.lines(x = c(.753, .942), y = c(.389, .389), gp = gpar(col = "#D1DFD4", lwd = .55))
  txt("Follow-up FI frailty", .8475, .366, 7.15, "bold", max_width = .196)
  grid.lines(x = c(.753, .942), y = c(.346, .346), gp = gpar(col = "#D1DFD4", lwd = .55))
  txt("Incident ADL limitation", .8475, .323, 7.15, "bold", max_width = .196)
  txt("Baseline limitation-free", .8475, .294, 6.6, max_width = .196)

  rect(.8475, .208, .218, .096, palette["neutral"])
  txt("Matched comparisons", .8475, .239, 7.15, "bold", palette["ink"], max_width = .196)
  txt("Complete baseline FI\nSame outcome-specific risk sets", .8475, .196, 6.6, max_width = .196)

  rect(.50, .075, .95, .110, palette["neutral"], border = "#97AAB2", radius_mm = 1)
  txt("Exploratory length–performance trade-off", .50, .110, 8.0, "bold", palette["ink"], max_width = .915)
  txt("Discrimination · FI correlation · MAE/RMSE · future outcomes", .50, .076, 7.1, max_width = .915)
  txt("Simple correlation · 1,000 random six-item sets · sensitivity analyses", .50, .043, 6.9, max_width = .915)
  popViewport()
  connector(c(.50,.50,.26,.26),c(.237,.220,.220,.198),width=.75)
  connector(c(.50,.50,.74,.74),c(.237,.220,.220,.198),width=.75)
  rect(.26,.115,.46,.165,palette["gold_light"],border="#BDAC78")
  txt("Primary training uncertainty",.26,.174,8.0,"bold",palette["ink"],max_width=.43)
  txt("500 coupled resamples; 7,500 XGBoost fits\nRelearn ranks, mappings & cutoffs\nOriginal + resampled held-out evaluation",.26,.107,7.0,max_width=.43,lineheight=1.3)
  rect(.74,.115,.46,.165,palette["green_light"],border="#82998A")
  txt("Heterogeneity & error sources",.74,.174,8.0,"bold",palette["ink"],max_width=.43)
  txt("Cohort FI, diagnoses, function & item profiles\nMapped MAE/RMSE gains; FI-level residuals\nMatched-score MHAS reference decomposition",.74,.107,7.0,max_width=.43,lineheight=1.3)
}

write.csv(nodes, file.path(source_dir, "Figure_1_framework_nodes.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(edges, file.path(source_dir, "Figure_1_framework_edges.csv"), row.names = FALSE, fileEncoding = "UTF-8")
w <- width_mm/25.4
h <- height_mm/25.4
svglite::svglite(paste0(stem, ".svg"), width = w, height = h, bg = "white")
draw_figure(); dev.off()

ragg::agg_png(paste0(stem, ".png"), width = width_mm, height = height_mm, units = "mm", res = 600, background = "white")
draw_figure(); dev.off()
ragg::agg_tiff(paste0(stem, ".tiff"), width = width_mm, height = height_mm, units = "mm", res = 600, compression = "lzw", background = "white")
draw_figure(); dev.off()
ragg::agg_png(file.path(output_dir, "Figure_1_preview.png"), width = width_mm, height = height_mm, units = "mm", res = 220, background = "white")
qa_text <- list(); draw_figure(); dev.off()
qa <- do.call(rbind, qa_text)
write.csv(qa, file.path(output_dir, "Figure_1_text_fit_QA.csv"), row.names = FALSE, fileEncoding = "UTF-8")
if (any(!qa$fits)) {
  print(qa[!qa$fits, ], row.names = FALSE)
  stop("One or more labels exceed their planned text boxes.")
}
legend <- paste0(
 "Figure 1. Five-cohort item-reduction framework. The study included 34,954 adults aged ≥65 years with at least one of hypertension, diabetes, heart disease, cancer, stroke, or arthritis. Each outer fold used four cohorts for learning and one held-out cohort for evaluation. A shared bank contained 16 questions, and complete reference FIs contained 39–50 deficits (frailty threshold ≥0.25). Three inner folds supplied cross-fitted TreeSHAP rankings for the primary universal scores; training-only latent classes supplied secondary phenotype rankings. Rankings, brief-score cutoffs, and six-item linear mappings were estimated in training and locked before held-out evaluation. Equal-weight top-k scores (k=1–16) required at least ceiling(0.8k) observed questions. Evaluation assessed concurrent FI discrimination and numerical agreement, follow-up FI frailty, and incident ADL limitation in baseline limitation-free participants, using common outcome-specific risk sets. Complete FI, simple correlation selection, and 1,000 random six-item combinations supplied benchmarks. Held-out length curves identified six items as an exploratory practical candidate. A further 500 coupled participant bootstraps repeated the primary universal training pipeline, including 7,500 XGBoost fits, ranking, mapping and cutoff estimation; duplicate participants stayed together within inner folds. Evaluation used original and jointly resampled held-out participants. The question bank and exploratory six-item length remained fixed during this resampling. Heterogeneity analyses compared FI distributions, diagnoses, functional limitations and candidate responses, quantified relative MAE/RMSE improvements and FI-stratified residual patterns, and decomposed MHAS bias at matched raw scores. Alt text: Three columns show cohort inputs, training-only score development and held-out evaluation. The lower modules show exploratory length trade-offs, primary training uncertainty and descriptive sources of heterogeneity."
)
writeLines(legend, file.path(output_dir, "Figure_1_legend.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"), useBytes = TRUE)
manifest <- data.frame(file = file.path("figures",list.files(output_dir, pattern="^Figure_1([.]|_preview|_legend|_text_fit)")), purpose = "Figure 1 revision output")
write.csv(manifest, file.path(source_dir, "Figure_1_output_manifest.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("PASS: R-only PNG/TIFF 600 dpi, editable SVG and preview exported.\n")
cat("PASS: All ", nrow(qa), " measured text lines fit their planned widths.\n", sep = "")
cat("Output directory: ", output_dir, "\n", sep = "")

