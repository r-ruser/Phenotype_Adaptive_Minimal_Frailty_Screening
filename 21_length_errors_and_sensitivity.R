if(.Platform$OS.type=="windows") invisible(try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE))
options(stringsAsFactors = FALSE, warn = 1)
try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"), silent = TRUE)
library(ggplot2)
library(patchwork)

# R-only derivation from previously executed held-out analyses; no model refitting.
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)[1])
package_dir <- Sys.getenv("MS5_PACKAGE_DIR", unset="")
if (!nzchar(package_dir)) package_dir <- dirname(normalizePath(script_file,winslash="/",mustWork=TRUE))
source(file.path(package_dir,"00_config.R"),encoding="UTF-8")
source_dir <- score_source
base_dir <- source_dir
output_root <- score_source
if (!nzchar(output_root)) stop("Set MS5_AGGREGATE_OUTPUT_DIR to a writable output directory")
out_dir <- file.path(output_root,"length_curve")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(out_dir), file.exists(file.path(source_dir, "FI_item_number_performance.csv")))
log_path <- file.path(out_dir, "length_curve_run.log")
sink(log_path, split = TRUE)
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")

required <- c("ggplot2", "patchwork", "ragg", "svglite")
stopifnot(all(vapply(required, requireNamespace, logical(1), quietly = TRUE)))
cohorts <- c("HRS", "ELSA", "SHARE", "CHARLS", "MHAS")
strategies <- c("Universal", "Phenotype")
outcomes <- c("concurrent", "followup_frail", "incident_adl")
keys <- c("cohort", "k", "strategy", "outcome")
metric_columns <- c("AUROC", "Spearman", "MAE", "RMSE")
perf_path <- file.path(source_dir, "FI_item_number_performance.csv")
cc_path <- file.path(source_dir, "FI_item_number_complete_case.csv")
perf <- read.csv(perf_path, check.names = FALSE)
cc <- read.csv(cc_path, check.names = FALSE)
perf <- perf[perf$strategy %in% strategies & perf$k %in% 1:16 & perf$outcome %in% outcomes, ]
cc <- cc[cc$strategy %in% strategies & cc$k %in% 1:16 & cc$outcome %in% outcomes, ]
sort_rows <- function(x) x[order(match(x$cohort, cohorts), match(x$strategy, strategies), match(x$outcome, outcomes), x$k), ]
perf <- sort_rows(perf); cc <- sort_rows(cc)
stopifnot(nrow(perf) == 480L, nrow(cc) == 480L, all(perf$rule == "80pct"), all(cc$rule == "complete_items"))
row_key <- function(x) do.call(paste, c(x[keys], sep = "|"))
stopifnot(!anyDuplicated(row_key(perf)), !anyDuplicated(row_key(cc)), setequal(row_key(perf), row_key(cc)))

inventory <- data.frame(
  source = c("80%-observed complete length curves", "Complete-selected-item sensitivity", "Paired selection-strategy comparisons", "Item stability summary", "Item stability rankings", "Exploratory length-rule results", "Exploratory length-rule summary", "Reference-composition sensitivity"),
  path = c(perf_path, cc_path, file.path(source_dir, "FI_universal_vs_phenotype.csv"), file.path(base_dir, "FI_item_stability.csv"), file.path(base_dir, "FI_item_stability_matrix.csv"), file.path(source_dir, "FI_optimal_item_number.csv"), file.path(source_dir, "FI_optimal_item_number_summary.csv"), file.path(source_dir, "JGMS_reference_overlap_sensitivity.csv"))
)
stopifnot(all(file.exists(inventory$path)))
inventory$rows <- vapply(inventory$path, function(p) nrow(read.csv(p)), integer(1))
inventory$columns <- vapply(inventory$path, function(p) paste(names(read.csv(p)), collapse = "; "), character(1))
write.csv(inventory, file.path(out_dir, "source_file_inventory.csv"), row.names = FALSE, na = "")

raw <- perf[perf$outcome == "concurrent", c("cohort", "k", "strategy", "n", "AUROC", "AUROC_low", "AUROC_high", "Spearman", "MAE", "RMSE")]
raw$score_scale <- "Raw equal-weight mean of observed selected items"
raw$reference_scale <- "Complete baseline FI; cohort-specific 39-50 deficits"
raw$risk_set <- "Common participants calculable under the 80%-observed rule for both strategies at every length 1-16"
write.csv(raw, file.path(out_dir, "FI_1_to_16_continuous_error_summary.csv"), row.names = FALSE, na = "")

macro_rows <- list()
for (strategy in strategies) for (k in 1:16) {
  z <- raw[raw$strategy == strategy & raw$k == k, ]
  stopifnot(nrow(z) == 5L)
  out <- data.frame(k = k, strategy = strategy, n_cohorts = 5L, n_total = sum(z$n), n_min = min(z$n), n_max = max(z$n))
  for (metric in metric_columns) {
    out[[paste0(metric, "_cohort_mean")]] <- mean(z[[metric]])
    out[[paste0(metric, "_cohort_min")]] <- min(z[[metric]])
    out[[paste0(metric, "_cohort_max")]] <- max(z[[metric]])
  }
  out$aggregation <- "Unweighted arithmetic mean across five held-out cohorts"
  macro_rows[[length(macro_rows) + 1L]] <- out
}
macro <- do.call(rbind, macro_rows)
write.csv(macro, file.path(out_dir, "FI_1_to_16_continuous_error_cohort_mean.csv"), row.names = FALSE)
key_lengths <- c(4L, 5L, 6L, 8L, 16L)
write.csv(macro[macro$k %in% key_lengths, ], file.path(out_dir, "FI_key_lengths_cohort_mean.csv"), row.names = FALSE)
write.csv(raw[raw$k %in% key_lengths, ], file.path(out_dir, "FI_key_lengths_cohort_results.csv"), row.names = FALSE)
table13 <- macro[macro$strategy == "Universal", c("k", "n_cohorts", "n_total", "AUROC_cohort_mean", "Spearman_cohort_mean", "MAE_cohort_mean", "RMSE_cohort_mean")]
names(table13)[names(table13) == "n_total"] <- "concurrent_n_total"
names(table13)[names(table13) == "AUROC_cohort_mean"] <- "concurrent_AUROC_cohort_mean"
for (outcome in c("followup_frail", "incident_adl")) {
  table13[[paste0(outcome, "_n_total")]] <- vapply(table13$k, function(k) sum(perf$n[perf$strategy == "Universal" & perf$outcome == outcome & perf$k == k]), numeric(1))
  table13[[paste0(outcome, "_AUROC_cohort_mean")]] <- vapply(table13$k, function(k) mean(perf$AUROC[perf$strategy == "Universal" & perf$outcome == outcome & perf$k == k]), numeric(1))
}
table13$aggregation <- "Unweighted arithmetic mean across five held-out cohorts; sample totals sum the cohort-specific common risk sets"
write.csv(table13, file.path(out_dir, "Supplementary_Table_13_universal_length_curve_summary.csv"), row.names = FALSE)

# A matched results-row comparison is distinct from a paired participant analysis.
# The complete-item rule retains a nested subset of each main analysis risk set.
paired <- merge(perf[, c(keys, "n", "events", "calculable_n", metric_columns)],
                cc[, c(keys, "n", "events", "calculable_n", metric_columns)],
                by = keys, suffixes = c("_80pct", "_complete_items"), sort = FALSE)
paired <- sort_rows(paired)
stopifnot(nrow(paired) == 480L, all(paired$n_complete_items <= paired$n_80pct))
paired$n_change <- paired$n_complete_items - paired$n_80pct
paired$n_retained_fraction <- paired$n_complete_items / paired$n_80pct
paired$risk_set_identical <- paired$n_change == 0L
for (metric in metric_columns) paired[[paste0("delta_", metric)]] <- paired[[paste0(metric, "_complete_items")]] - paired[[paste0(metric, "_80pct")]]
paired$delta_definition <- "Complete-selected-item result minus 80%-observed result"
paired$comparison_scope <- "Complete-item subset nested within common all-length main risk set; sample composition varies with length and strategy"
write.csv(paired, file.path(out_dir, "FI_complete_item_vs_80pct_by_length.csv"), row.names = FALSE, na = "")

range_stat <- function(x, fun) if (all(!is.finite(x))) NA_real_ else fun(x[is.finite(x)])
summary_rows <- list()
for (cohort in cohorts) for (strategy in strategies) for (outcome in outcomes) {
  z <- paired[paired$cohort == cohort & paired$strategy == strategy & paired$outcome == outcome, ]
  stopifnot(nrow(z) == 16L, length(unique(z$n_80pct)) == 1L)
  out <- data.frame(cohort = cohort, strategy = strategy, outcome = outcome, n_lengths = 16L,
                    n_80pct = unique(z$n_80pct), n_complete_min = min(z$n_complete_items), n_complete_max = max(z$n_complete_items),
                    n_retained_fraction_min = min(z$n_retained_fraction), n_retained_fraction_max = max(z$n_retained_fraction),
                    lengths_with_same_n = paste(z$k[z$risk_set_identical], collapse = ";"),
                    lengths_with_fewer_participants = paste(z$k[!z$risk_set_identical], collapse = ";"))
  for (metric in metric_columns) {
    d <- z[[paste0("delta_", metric)]]
    out[[paste0("delta_", metric, "_min")]] <- range_stat(d, min)
    out[[paste0("delta_", metric, "_max")]] <- range_stat(d, max)
    out[[paste0("abs_delta_", metric, "_max")]] <- range_stat(abs(d), max)
  }
  out$interpretation <- "Descriptive robustness across answering rules; stricter-rule subsets also change sample composition"
  summary_rows[[length(summary_rows) + 1L]] <- out
}
sensitivity <- do.call(rbind, summary_rows)
stopifnot(nrow(sensitivity) == 30L)
write.csv(sensitivity, file.path(out_dir, "FI_complete_item_sensitivity_summary_30_rows.csv"), row.names = FALSE, na = "")

# Figure contract: numerical agreement varies with length and strategy, with
# transparent cohort-specific raw-scale errors. Two quantitative panels show MAE/RMSE.
colors <- c(HRS = "#4E79A7", ELSA = "#59A14F", SHARE = "#E15759", CHARLS = "#B07AA1", MHAS = "#F28E2B")
fig_data <- do.call(rbind, lapply(c("MAE", "RMSE"), function(metric) {
  data.frame(cohort = raw$cohort, k = raw$k, strategy = raw$strategy, n = raw$n, metric = metric, value = raw[[metric]])
}))
fig_data$cohort <- factor(fig_data$cohort, levels = cohorts)
fig_data$strategy <- factor(fig_data$strategy, levels = strategies)
stopifnot(nrow(fig_data) == 320L, all(is.finite(fig_data$value)))
stopifnot(max(fig_data$value[fig_data$metric == "MAE"]) < .40,
          max(fig_data$value[fig_data$metric == "RMSE"]) < .50)
write.csv(fig_data, file.path(out_dir, "Figure_S4_source_data.csv"), row.names = FALSE)

base_theme <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = .35, colour = "#28343D"),
        axis.ticks = element_line(linewidth = .35, colour = "#28343D"),
        axis.text = element_text(size = 7, colour = "#28343D"),
        axis.title = element_text(size = 8, colour = "#28343D"),
        legend.position = "bottom", legend.title = element_text(size = 7), legend.text = element_text(size = 7),
        legend.key.width = grid::unit(9, "mm"), legend.key.height = grid::unit(4, "mm"),
        legend.box = "vertical", legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(4, 6, 3, 5), plot.background = element_rect(fill = "white", colour = NA),
        panel.grid = element_blank())
one_panel <- function(metric, ylab) {
  z <- fig_data[fig_data$metric == metric, ]
  ggplot(z, aes(k, value, colour = cohort, linetype = strategy, group = interaction(cohort, strategy))) +
    geom_line(linewidth = .55, alpha = .94) +
    scale_colour_manual(values = colors, name = "Held-out cohort") +
    scale_linetype_manual(values = c(Universal = "solid", Phenotype = "22"), name = "Selection") +
    scale_x_continuous(breaks = c(1, 4, 6, 8, 12, 16), limits = c(1, 16), expand = expansion(mult = c(.02, .03))) +
    scale_y_continuous(limits = c(0, if (metric == "MAE") .40 else .50), breaks = if (metric == "MAE") seq(0, .40, .10) else seq(0, .50, .10), expand = expansion(mult = c(0, .02))) +
    labs(x = "Number of selected items", y = ylab) +
    guides(colour = guide_legend(order = 1, nrow = 1, override.aes = list(linetype = "solid", linewidth = .75)),
           linetype = guide_legend(order = 2, nrow = 1, override.aes = list(colour = "#263238", linewidth = .75))) +
    base_theme
}
fig <- (one_panel("MAE", "Mean absolute error") | one_panel("RMSE", "Root mean squared error")) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom", legend.box = "vertical", plot.tag = element_text(size = 10, face = "bold", family = "Arial"))

width_mm <- 183; height_mm <- 105; dpi <- 600
prefix <- file.path(out_dir, "Figure_S4_continuous_error_by_length")
svglite::svglite(paste0(prefix, ".svg"), width = width_mm / 25.4, height = height_mm / 25.4)
print(fig); dev.off()
ragg::agg_png(paste0(prefix, ".png"), width = width_mm, height = height_mm, units = "mm", res = dpi, background = "white")
print(fig); dev.off()
ragg::agg_tiff(paste0(prefix, ".tiff"), width = width_mm, height = height_mm, units = "mm", res = dpi, compression = "lzw", background = "white")
print(fig); dev.off()
ragg::agg_png(paste0(prefix, "_preview.png"), width = width_mm, height = height_mm, units = "mm", res = 180, background = "white")
print(fig); dev.off()

paired_strategy_path <- file.path(source_dir, "FI_universal_vs_phenotype.csv")
paired_strategy <- read.csv(paired_strategy_path)
rules <- read.csv(file.path(source_dir, "FI_optimal_item_number_summary.csv"))
stability <- read.csv(file.path(base_dir, "FI_item_stability.csv"))
overlap <- read.csv(file.path(source_dir, "JGMS_reference_overlap_sensitivity.csv"))
analysis_inventory <- data.frame(
  analysis = c("Complete length curves", "Complete-selected-item sensitivity", "Universal versus phenotype-selected scores", "Paired universal versus phenotype selection", "Forced-distinct question combinations", "Selection stability", "Exploratory length criteria", "Reference composition", "Simple-correlation and random-combination benchmark"),
  scope = c("1-16 items; 5 cohorts; 2 strategies; concurrent FI, follow-up FI, incident ADL",
            "1-16 items; 5 cohorts; 2 strategies; 3 outcomes",
            "1-16 items; 5 cohorts; 3 outcomes",
            "5, 6, and 8 items; 5 cohorts; 3 outcomes",
            "5 items only; 5 cohorts; 3 outcomes",
            "16 candidate-item rankings across 5 universal training folds and 29 training phenotypes",
            "1-16 items; concurrent AUROC fractions 90%, 95%, 97.5%; follow-up FI alone or follow-up FI plus ADL; 4/5 or 5/5 cohorts",
            "Universal 6-item score versus complete FI and FI omitting all 16 candidate items; paired complete-reference subset",
            "6 items only; 5 outer folds; 1,000 random 6-item combinations; training-only simple correlation"),
  method = c("Previously locked training-selected equal-weight scores evaluated on common all-length held-out risk sets",
             "Re-evaluate the locked scores among participants with every selected item observed; compare with 80%-observed rule",
             "Direct comparison of previously executed curves; common all-length risk sets",
             "500 paired participant-bootstrap samples for AUROC differences",
             "Previously executed rule-enforced different phenotype question sets; 500 paired participant-bootstrap samples",
             "Rank summaries and top-3/top-5/top-8 inclusion frequencies from training-only selections",
             "Post-evaluation descriptive criteria combining concurrent discrimination, Spearman >=0.80, and prospective excess-AUROC retention >=0.85",
             "Evaluate altered reference definition; report AUROC, Spearman, and reference-positive fractions",
             "Fixed-seed random combinations and simple-correlation ranking evaluated at the same six-item length on aligned benchmark risk sets"),
  result_rows = c(nrow(perf), nrow(cc), nrow(perf), sum(paired_strategy$comparison == "Phenotype_minus_Universal"), sum(paired_strategy$comparison == "ForcedDistinct_minus_Universal"), nrow(stability), nrow(rules), nrow(overlap), nrow(read.csv(file.path(source_dir, "benchmarks", "six_item_benchmark_performance.csv")))),
  source_files = c("analysis/results/FI_item_number_performance.csv",
                   "analysis/results/FI_item_number_complete_case.csv",
                   "analysis/results/FI_item_number_performance.csv",
                   "analysis/results/FI_universal_vs_phenotype.csv",
                   "analysis/results/FI_universal_vs_phenotype.csv",
                   "FI_item_stability.csv; FI_item_stability_matrix.csv",
                   "analysis/results/FI_optimal_item_number.csv; FI_optimal_item_number_summary.csv",
                   "analysis/results/JGMS_reference_overlap_sensitivity.csv",
                   "submission_revision_20261001/benchmarks/six_item_benchmark_performance.csv; SHAP_vs_simple_correlation_paired_bootstrap.csv; random_combination_distribution_summary.csv")
)
write.csv(analysis_inventory, file.path(out_dir, "FI_length_selection_sensitivity_inventory.csv"), row.names = FALSE)

fmt <- function(x, digits = 3L) formatC(x, digits = digits, format = "f")
u <- macro[macro$strategy == "Universal" & macro$k %in% key_lengths, ]
u6 <- raw[raw$strategy == "Universal" & raw$k == 6L, ]
concurrent_pairs <- paired[paired$outcome == "concurrent", ]
complete_ranges <- vapply(outcomes, function(o) {
  z <- paired[paired$outcome == o, ]
  paste0(o, ": AUROC differences ", fmt(min(z$delta_AUROC), 4), " to ", fmt(max(z$delta_AUROC), 4))
}, character(1))
key_sentence <- paste(vapply(seq_len(nrow(u)), function(i) paste0(u$k[i], " items: MAE ", fmt(u$MAE_cohort_mean[i]), ", RMSE ", fmt(u$RMSE_cohort_mean[i])), character(1)), collapse = "; ")
text <- c(
  "# Length curves, continuous error, and sensitivity analyses",
  "",
  "## Manuscript Methods insertion",
  "We evaluated the locked universal and phenotype-selected equal-weight scores at every length from 1 to 16 items. Concurrent evaluation included AUROC, Spearman correlation, mean absolute error (MAE), and root mean squared error (RMSE) against the complete baseline FI. MAE and RMSE use the raw mean of observed selected items. Each cohort's primary curves used a common participant set with calculable scores under the 80%-observed rule for both selection strategies at all 16 lengths. Complete-selected-item sensitivity analyses used the subset with every selected item observed at each length. This descriptive comparison captures changes in both answering rules and participant composition.",
  "",
  "## Manuscript Results insertion",
  paste0("Across all lengths, the concurrent error analyses covered 160 cohort-by-length-by-strategy results. The unweighted five-cohort means for the universal score were ", key_sentence, ". At six items, cohort-specific raw MAE ranged from ", fmt(min(u6$MAE)), " to ", fmt(max(u6$MAE)), " and RMSE from ", fmt(min(u6$RMSE)), " to ", fmt(max(u6$RMSE)), ". Error changed non-monotonically with length; numerical agreement and discrimination therefore describe complementary aspects of compression (Supplementary Figure S4)."),
  paste0("Complete-selected-item sensitivity analyses covered all 16 lengths, both selection strategies, and all three outcomes (480 result comparisons). ", paste(complete_ranges, collapse = "; "), ". Across concurrent evaluations, differences in MAE ranged from ", fmt(min(concurrent_pairs$delta_MAE)), " to ", fmt(max(concurrent_pairs$delta_MAE)), ", and differences in RMSE from ", fmt(min(concurrent_pairs$delta_RMSE)), " to ", fmt(max(concurrent_pairs$delta_RMSE)), ". These differences describe stricter-rule subsets nested within the main common risk sets."),
  "The analysis sequence comprised full 1-16-item length curves, complete-selected-item sensitivity, universal versus phenotype-selected curves, paired strategy comparisons at 5, 6, and 8 items, forced-distinct five-item combinations, training-ranking stability across 29 phenotypes, alternative descriptive length criteria, and reference-composition sensitivity. The simple-correlation and random-combination benchmarks specifically assessed the six-item candidate.",
  "",
  "## Interpretation of candidate length",
  "Length criteria were applied descriptively to the five held-out performance curves. At the 90% and 95% concurrent AUROC criteria, the combination of Spearman correlation >=0.80 and follow-up FI excess-AUROC retention >=0.85 selected five items for four of five cohorts and six items for all five cohorts. Adding the ADL retention condition selected 11 items for all five cohorts. The 97.5% criterion selected six items for four of five cohorts under the same prospective definitions. These results identify outcome-specific candidate lengths and support six items as a promising practical candidate for further prospective validation.",
  "",
  "## Supplementary Figure S4 legend",
  "Supplementary Figure S4. Raw numerical error across score lengths. (a) Mean absolute error and (b) root mean squared error of 1-16-item scores against each cohort's complete baseline FI. Colours identify the held-out cohorts; solid lines indicate universal selection and dashed lines indicate phenotype-specific selection. Scores equal the mean of the observed selected items, with at least 80% of selected items observed. Each cohort uses the common participant set calculable for both strategies at every length: HRS n=8,075; ELSA n=4,196; SHARE n=12,552; CHARLS n=3,946; MHAS n=5,433. Lines show point estimates from previously locked scores. Source data comprise 320 metric-by-cohort-by-length-by-strategy values. The raw scores have their own numerical scale; the six-item training-mapped agreement analysis appears separately in Figure 3.",
  "",
  "## Full-data notes",
  "The five-cohort means are unweighted descriptive averages, with cohort-level estimates supplied in the accompanying CSV. Sensitivity deltas use complete-selected-item minus 80%-observed results. The available complete-item sensitivity tables provide point estimates; the paired comparison here matches result rows and has a participant-composition component. All output is derived from existing executed results, with the original locked models and scores preserved."
)
writeLines(text, file.path(out_dir, "length_curve_methods_results.md"), useBytes = TRUE)
writeLines("Supplementary Figure S4. Raw numerical error across score lengths. (a) Mean absolute error and (b) root mean squared error of 1-16-item scores against each cohort's complete baseline FI. Colours identify the held-out cohorts; solid lines indicate universal selection and dashed lines indicate phenotype-specific selection. Scores equal the mean of the observed selected items, with at least 80% of selected items observed. Each cohort uses the common participant set calculable for both strategies at every length: HRS n=8,075; ELSA n=4,196; SHARE n=12,552; CHARLS n=3,946; MHAS n=5,433. Lines show point estimates from previously locked scores. Source data comprise 320 metric-by-cohort-by-length-by-strategy values. The raw scores have their own numerical scale; the six-item training-mapped agreement analysis appears separately in Figure 3.", file.path(out_dir, "Figure_S4_legend.txt"))

qa <- c(
  paste0("PASS: main 1-16-item rows = ", nrow(perf)),
  paste0("PASS: complete-selected-item rows = ", nrow(cc)),
  paste0("PASS: concurrent raw-error rows = ", nrow(raw)),
  paste0("PASS: row-matched answering-rule comparisons = ", nrow(paired)),
  paste0("PASS: cohort x strategy x outcome summaries = ", nrow(sensitivity)),
  paste0("PASS: complete-item n <= main n in all ", nrow(paired), " comparisons"),
  paste0("PASS: figure source values = ", nrow(fig_data)),
  "PASS: five-cohort palette matches the original length-curve figure",
  "PASS: plotting and all visual exports generated in R",
  "PASS: export size 183 x 105 mm; PNG and TIFF 600 dpi; editable vector SVG and SVG",
  "PASS: source files preserved; no raw model refitting",
  "Visual inspection status: recorded separately in Figure_S4_visual_QA.txt and length_curve_final_QA.txt"
)
writeLines(qa, file.path(out_dir, "length_curve_numerical_QA.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_sessionInfo.txt"))
files <- list.files(out_dir, full.names = TRUE)
manifest <- data.frame(file = basename(files), bytes = file.info(files)$size, purpose = ifelse(grepl("Figure_S4", files), "Continuous-error figure and source data", "Length and sensitivity derivation, results, and audit"))
write.csv(manifest, file.path(out_dir, "output_manifest.csv"), row.names = FALSE)
cat("Key universal length values:\n"); print(u[, c("k", "AUROC_cohort_mean", "Spearman_cohort_mean", "MAE_cohort_mean", "RMSE_cohort_mean")], row.names = FALSE)
cat("Complete-item AUROC delta range by outcome:\n"); cat(paste(complete_ranges, collapse = "\n"), "\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
sink()
