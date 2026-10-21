if(.Platform$OS.type=="windows") invisible(try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE))
# Reproduce the six-panel length comparison from saved aggregate results.
# Backend: R; no model refitting. Run with --revision-dir <submission directory>.
options(stringsAsFactors = FALSE, warn = 1)
invisible(try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"), silent = TRUE))
required <- c("ggplot2", "patchwork", "ragg", "svglite")
stopifnot(all(vapply(required, requireNamespace, logical(1), quietly = TRUE)))
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
options(stringsAsFactors=FALSE)
ms5_args <- commandArgs(FALSE)
ms5_script <- sub("^--file=", "", ms5_args[grepl("^--file=",ms5_args)][1])
code_dir <- Sys.getenv("MS5_CODE_DIR", unset=if(!is.na(ms5_script)) dirname(normalizePath(ms5_script,winslash="/",mustWork=TRUE)) else getwd())
source(file.path(code_dir,"00_config.R"),encoding="UTF-8")

revision_dir <- publication_dir
results_dir <- score_source
figure_dir <- file.path(publication_dir,"figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
cohorts <- c("HRS", "ELSA", "SHARE", "CHARLS", "MHAS")
colors <- c(HRS = "#4E79A7", ELSA = "#59A14F", SHARE = "#E15759", CHARLS = "#B07AA1", MHAS = "#F28E2B")
ink <- "#202B33"
p <- read.csv(file.path(results_dir, "FI_item_number_performance.csv"), check.names = FALSE, fileEncoding = "UTF-8")
p <- p[p$strategy == "Universal" & p$k %in% 1:16, ]
stopifnot(nrow(p) == 240L, all(p$rule == "80pct"))
meta <- read.csv(file.path(figure_dir, "Figure_2_pooled_source.csv"), check.names = FALSE, fileEncoding = "UTF-8")
stopifnot(nrow(meta) == 64L, all(meta$n_cohorts == 5L))
macro <- read.csv(file.path(results_dir, "length_curve", "FI_1_to_16_continuous_error_cohort_mean.csv"), check.names = FALSE, fileEncoding = "UTF-8")
macro <- macro[macro$strategy == "Universal", ]
stopifnot(nrow(macro) == 16L)

# Scientific contract: measure the changing balance of discrimination, rank
# correlation, raw-scale agreement and future outcome discrimination at every
# length. The first four metrics share common concurrent participants; future
# outcomes use their own common risk sets. Six items remains an exploratory length.
spec <- data.frame(panel = letters[1:6], outcome = c("concurrent", "concurrent", "concurrent", "concurrent", "followup_frail", "incident_adl"),
                   metric = c("AUROC", "Spearman", "MAE", "RMSE", "AUROC", "AUROC"),
                   title = c("Concurrent FI frailty", "Continuous FI correlation", "Raw-scale absolute error", "Raw-scale root mean squared error", "Follow-up FI frailty", "Incident ADL limitation"),
                   ylab = c("AUROC", "Spearman correlation", "MAE", "RMSE", "AUROC", "AUROC"))
cohort_rows <- list(); summary_rows <- list()
for (j in seq_len(nrow(spec))) {
  s <- spec[j, ]; z <- p[p$outcome == s$outcome, ]
  z <- z[order(match(z$cohort, cohorts), z$k), ]
  stopifnot(nrow(z) == 80L, all(is.finite(z[[s$metric]])))
  z$panel <- s$panel; z$metric <- s$metric; z$value <- z[[s$metric]]
  cohort_rows[[j]] <- z[, c("panel", "cohort", "k", "strategy", "outcome", "metric", "value", "n", "events", "rule")]
  if (s$metric %in% c("MAE", "RMSE")) {
    mm <- data.frame(k = macro$k, pooled = macro[[paste0(s$metric, "_cohort_mean")]], low = NA_real_, high = NA_real_,
                     aggregation = "Unweighted arithmetic mean across the five held-out cohorts")
    direct <- vapply(mm$k, function(k) mean(z[[s$metric]][z$k == k]), numeric(1))
    stopifnot(max(abs(direct - mm$pooled)) < 1e-12)
  } else {
    mm <- meta[meta$outcome == s$outcome & meta$metric == s$metric, c("k", "pooled", "low", "high")]
    mm$aggregation <- if (s$metric == "AUROC") "Saved pointwise REML random-effects summary on the logit AUROC scale" else "Saved pointwise REML random-effects summary on the Fisher correlation scale"
  }
  stopifnot(nrow(mm) == 16L, setequal(mm$k, 1:16))
  mm$panel <- s$panel; mm$outcome <- s$outcome; mm$metric <- s$metric
  summary_rows[[j]] <- mm[, c("panel", "outcome", "metric", "k", "pooled", "low", "high", "aggregation")]
}
cohort_data <- do.call(rbind, cohort_rows)
summary_data <- do.call(rbind, summary_rows)
write.csv(cohort_data, file.path(figure_dir, "Figure_2_cohort_source.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(summary_data, file.path(figure_dir, "Figure_2_summary_source.csv"), row.names = FALSE, fileEncoding = "UTF-8", na = "")

th <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = .3, colour = ink), axis.ticks = element_line(linewidth = .3, colour = ink),
        axis.text = element_text(size = 7.2, colour = ink), axis.title = element_text(size = 8, colour = ink),
        plot.title = element_text(size = 8.6, face = "bold", colour = ink),
        panel.grid.major.y = element_line(colour = "#E8EBED", linewidth = .2),
        plot.margin = margin(3, 5, 3, 3), legend.position = "bottom", legend.text = element_text(size = 7.3),
        legend.key.width = unit(11, "mm"), legend.margin = margin(1, 0, 1, 0))
panels <- list()
for (j in seq_len(nrow(spec))) {
  s <- spec[j, ]; z <- cohort_data[cohort_data$panel == s$panel, ]; mm <- summary_data[summary_data$panel == s$panel, ]
  z$cohort <- factor(z$cohort, levels = cohorts)
  g <- ggplot(z, aes(k, value, group = cohort, colour = cohort))
  band <- s$metric %in% c("AUROC", "Spearman")
  if (band) g <- g + geom_ribbon(data = mm, aes(x = k, ymin = low, ymax = high), inherit.aes = FALSE, fill = ink, alpha = .11)
  g <- g + geom_line(linewidth = .47, alpha = .90) + geom_point(size = .85, alpha = .95) +
    geom_line(data = mm, aes(x = k, y = pooled, group = 1), inherit.aes = FALSE, colour = ink, linewidth = .9) +
    scale_colour_manual(values = colors, limits = cohorts, drop = FALSE) +
    scale_x_continuous(breaks = c(1, 4, 6, 8, 12, 16), limits = c(1, 16), expand = expansion(mult = .02)) +
    labs(title = paste0(s$panel, "  ", s$title), x = "Number of items", y = s$ylab, colour = NULL) + th
  if (s$metric == "MAE") g <- g + scale_y_continuous(limits = c(0, .4), breaks = seq(0, .4, .1), expand = expansion(mult = .015))
  if (s$metric == "RMSE") g <- g + scale_y_continuous(limits = c(0, .45), breaks = seq(0, .4, .1), expand = expansion(mult = .015))
  if (s$metric %in% c("AUROC", "Spearman")) g <- g + scale_y_continuous(limits = c(.5, 1), breaks = seq(.5, 1, .1), expand = expansion(mult = .015))
  points <- ggplot_build(g)$data[[if (band) 3L else 2L]]
  key <- paste(as.character(z$cohort), z$k)
  idx <- match(paste(names(colors)[match(points$colour, colors)], points$x), key)
  stopifnot(!anyNA(idx), nrow(points) == 80L, max(abs(points$y - z$value[idx])) < 1e-14)
  panels[[j]] <- g
}
fig <- wrap_plots(panels, ncol = 2, byrow = TRUE, guides = "collect") & theme(legend.position = "bottom")
stem <- file.path(figure_dir, "Figure_2")
w <- 183; h <- 210
ggsave(paste0(stem, ".png"), fig, width = w, height = h, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(paste0(stem, ".tiff"), fig, width = w, height = h, units = "mm", dpi = 600, compression = "lzw", device = ragg::agg_tiff, bg = "white")
ggsave(paste0(stem, ".svg"), fig, width = w, height = h, units = "mm", device = svglite::svglite, bg = "white")
ggsave(paste0(stem, "_preview.png"), fig, width = w, height = h, units = "mm", dpi = 180, device = ragg::agg_png, bg = "white")
legend <- paste0("Figure 2. Length-dependent discrimination, rank correlation, numerical error and prospective performance of universal equal-weight scores. ",
                 "Panels show (a) concurrent complete-FI frailty AUROC, (b) continuous complete-FI Spearman correlation, (c) raw-scale MAE, (d) raw-scale RMSE, ",
                 "(e) follow-up FI frailty AUROC and (f) incident ADL limitation AUROC. Colored lines show the five held-out cohorts. ",
                 "Dark lines in panels a, b, e and f show the original pointwise random-effects summaries with 95% intervals (gray bands); dark lines in panels c and d show unweighted cohort means. ",
                 "No uncertainty band is assigned to the descriptive error means. Concurrent metrics use common participants calculable under the 80%-observed rule across both strategies and all lengths ",
                 "(HRS 8,075; ELSA 4,196; SHARE 12,552; CHARLS 3,946; MHAS 5,433; total 34,202). ",
                 "The common outcome-specific totals are 21,861 for follow-up FI and 18,503 for incident ADL limitation. ",
                 "MAE and RMSE compare each raw equal-weight score with the cohort-specific 39-50-deficit FI. Length comparisons are descriptive and six items represents an exploratory candidate length.")
writeLines(legend, file.path(figure_dir, "Figure_2_legend.txt"), useBytes = TRUE)
png_dimensions <- function(path) {con <- file(path, "rb"); on.exit(close(con)); seek(con, 16); readBin(con, "integer", n = 2, size = 4, endian = "big")}
dims <- png_dimensions(paste0(stem, ".png"))
stopifnot(all(abs(dims - c(w, h)/25.4*600) <= 1L))
svg <- readLines(paste0(stem, ".svg"), encoding = "UTF-8", warn = FALSE)
stopifnot(any(grepl("<text", svg, fixed = TRUE)), any(grepl("Arial", svg, fixed = TRUE)))
qa <- c("PASS: all 480 cohort point values match saved all-length results.",
        "PASS: the 64 original random-effects summaries are preserved; 32 descriptive error means match the saved five-cohort means.",
        "Figure dimensions: 183 x 210 mm; PNG/TIFF 600 dpi; SVG has editable Arial text.",
        paste("PNG pixels:", paste(dims, collapse = " x ")), "No models were refitted; only saved aggregate results were rendered.")
writeLines(qa, file.path(figure_dir, "Figure_2_QA.txt"), useBytes = TRUE)
cat(paste(qa, collapse = "\n"), "\n")



# Official caption is shared with the submitted DOCX.
canonical_path <- file.path(revision_dir,"figure_source/figure_legends_master.csv")
if(file.exists(canonical_path)) {
  canonical <- read.csv(canonical_path,fileEncoding="UTF-8",stringsAsFactors=FALSE)
  writeLines(canonical$legend[canonical$figure_id=="2"],file.path(revision_dir,"figures/Figure_2_legend.txt"),useBytes=TRUE)
}
