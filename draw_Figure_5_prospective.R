if(.Platform$OS.type=="windows") invisible(try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE))
# Reproduce prospective comparisons from saved paired results.
# Backend: R; run with --revision-dir <submission directory>.
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
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
cohorts <- c("HRS", "ELSA", "SHARE", "CHARLS", "MHAS")
re <- read.csv(file.path(runtime_dir,"extensions/approximation", "FI6_paired_prospective_retention.csv"), fileEncoding = "UTF-8", check.names = FALSE)
p <- read.csv(file.path(results_dir, "FI_item_number_performance.csv"), fileEncoding = "UTF-8", check.names = FALSE)
stopifnot(nrow(re) == 10L, !anyDuplicated(paste(re$cohort, re$outcome)))
long <- do.call(rbind, lapply(c("Full_FI", "Universal"), function(strategy) {
  stem <- if (strategy == "Full_FI") "AUROC_full" else "AUROC_FI6"
  data.frame(cohort = re$cohort, outcome = re$outcome, strategy = strategy, k = if (strategy == "Full_FI") 0L else 6L,
             n = re$n, events = re$events, AUROC = re[[stem]], AUROC_low = re[[paste0(stem, "_low")]], AUROC_high = re[[paste0(stem, "_high")]],
             interval_definition = "500 paired participant-bootstrap resamples; conditional on fixed scoring rules")
}))
# Each pair was computed on the same outcome-specific common participants.
orig <- p[p$outcome %in% c("followup_frail", "incident_adl") & ((p$k == 6L & p$strategy == "Universal") | p$strategy == "Full_FI"), ]
ix <- match(paste(long$cohort, long$outcome, long$strategy), paste(orig$cohort, orig$outcome, orig$strategy))
stopifnot(nrow(long) == 20L, nrow(orig) == 20L, !anyNA(ix), all(long$n == orig$n[ix]), all(long$events == orig$events[ix]),
          max(abs(long$AUROC - orig$AUROC[ix])) < 1e-12, all(long$AUROC_low <= long$AUROC), all(long$AUROC <= long$AUROC_high))
write.csv(long, file.path(figure_dir, "Figure_5_source_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")
colors <- c(Full_FI = "#263D34", Universal = "#176B8A")
labels <- c(Full_FI = "Complete baseline FI", Universal = "Universal six-item score")
th <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = .3, colour = "#202B33"), axis.ticks = element_line(linewidth = .3),
        axis.text = element_text(size = 7.7, colour = "#202B33"), axis.title = element_text(size = 8),
        plot.title = element_text(size = 8.7, face = "bold", colour = "#202B33"),
        panel.grid.major.y = element_line(colour = "#EBEEF0", linewidth = .2),
        plot.margin = margin(3, 5, 3, 3), legend.text = element_text(size = 7.8),
        legend.key.width = unit(11, "mm"), legend.position = "bottom", legend.margin = margin(1, 0, 1, 0))
panels <- list()
outcomes <- c("followup_frail", "incident_adl")
titles <- c("a  Follow-up FI frailty", "b  Incident ADL limitation")
for (j in seq_along(outcomes)) {
  z <- long[long$outcome == outcomes[j], ]
  z$cohort <- factor(z$cohort, levels = rev(cohorts)); z$strategy <- factor(z$strategy, levels = names(colors))
  z$y <- as.numeric(z$cohort) + ifelse(z$strategy == "Full_FI", .13, -.13)
  n_by_cohort <- re[re$outcome == outcomes[j], ]; counts <- n_by_cohort$n[match(rev(cohorts), n_by_cohort$cohort)]
  ylabels <- paste0(rev(cohorts), "   (n = ", format(counts, big.mark = ",", trim = TRUE), ")")
  g <- ggplot(z, aes(AUROC, y, colour = strategy, shape = strategy)) +
    geom_segment(aes(x = AUROC_low, xend = AUROC_high, yend = y), linewidth = .65) +
    geom_point(size = 2.0) + scale_colour_manual(values = colors, limits = names(colors), labels = labels, drop = FALSE) +
    scale_shape_manual(values = c(Full_FI = 16, Universal = 17), limits = names(colors), labels = labels, drop = FALSE) +
    scale_y_continuous(breaks = 1:5, labels = ylabels, limits = c(.5, 5.5), expand = c(0, 0)) +
    scale_x_continuous(limits = c(.5, 1), breaks = seq(.5, 1, .1), expand = expansion(mult = .015)) +
    labs(title = titles[j], x = "AUROC", y = NULL, colour = NULL, shape = NULL) + th
  point_data <- ggplot_build(g)$data[[2L]]
  stopifnot(nrow(point_data) == 10L, max(abs(point_data$x - z$AUROC)) < 1e-14)
  panels[[j]] <- g
}
fig <- wrap_plots(panels, ncol = 1, guides = "collect") & theme(legend.position = "bottom")
stem <- file.path(figure_dir, "Figure_5"); w <- 183; h <- 160
ggsave(paste0(stem, ".png"), fig, width = w, height = h, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(paste0(stem, ".tiff"), fig, width = w, height = h, units = "mm", dpi = 600, compression = "lzw", device = ragg::agg_tiff, bg = "white")
ggsave(paste0(stem, ".svg"), fig, width = w, height = h, units = "mm", device = svglite::svglite, bg = "white")
ggsave(paste0(stem, "_preview.png"), fig, width = w, height = h, units = "mm", dpi = 180, device = ragg::agg_png, bg = "white")
legend <- paste0("Figure 5. Prospective discrimination of complete baseline FI and the universal six-item score on identical outcome-specific common participants. ",
                 "Panels show (a) follow-up FI frailty and (b) incident ADL limitation. Points indicate AUROC; horizontal intervals are percentile 95% intervals from 500 paired participant-bootstrap resamples, ",
                 "conditional on the fixed scoring rules. Both scores use the same participants within each cohort-outcome comparison; cohort-specific sample sizes appear beside each row. ",
                 "Follow-up FI frailty is defined by complete follow-up FI >=0.25. The totals are 21,861 for follow-up FI and 18,503 for incident ADL limitation. ",
                 "The displayed intervals describe each score's AUROC; paired AUROC differences and excess-over-chance AUROC retention are reported in Supplementary Table 7.")
writeLines(legend, file.path(figure_dir, "Figure_5_legend.txt"), useBytes = TRUE)
png_dimensions <- function(path) {con <- file(path, "rb"); on.exit(close(con)); seek(con, 16); readBin(con, "integer", n = 2, size = 4, endian = "big")}
dims <- png_dimensions(paste0(stem, ".png")); stopifnot(all(abs(dims - c(w, h)/25.4*600) <= 1L))
svg <- readLines(paste0(stem, ".svg"), encoding = "UTF-8", warn = FALSE)
stopifnot(any(grepl("<text", svg, fixed = TRUE)), any(grepl("Arial", svg, fixed = TRUE)))
qa <- c("PASS: 20 AUROC points match saved main-analysis estimates; all cohort-outcome n and event counts match.",
        "PASS: displayed 95% intervals exactly reuse the saved 500-resample paired participant-bootstrap outputs.",
        "Figure dimensions: 183 x 160 mm; PNG/TIFF 600 dpi; editable Arial SVG.", paste("PNG pixels:", paste(dims, collapse = " x ")),
        "Phenotype-specific results remain in the supplementary evidence; no analysis was refitted.")
writeLines(qa, file.path(figure_dir, "Figure_5_QA.txt"), useBytes = TRUE)
cat(paste(qa, collapse = "\n"), "\n")

# Official caption is shared with the submitted DOCX.
canonical_path <- file.path(revision_dir,"figure_source/figure_legends_master.csv")
if(file.exists(canonical_path)) {
  canonical <- read.csv(canonical_path,fileEncoding="UTF-8",stringsAsFactors=FALSE)
  writeLines(canonical$legend[canonical$figure_id=="5"],file.path(revision_dir,"figures/Figure_5_legend.txt"),useBytes=TRUE)
}
