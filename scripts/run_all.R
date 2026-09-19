options(stringsAsFactors = FALSE)
scripts <- c(
  "scripts/01_recover_same_wave_data.R",
  "scripts/02_lca_transportability.R",
  "scripts/03_irt_dif.R",
  "scripts/04_screening_validation.R",
  "scripts/05_nature_figures_and_report.R"
)
for (s in scripts) {
  cat("\n>>> Running", s, "\n")
  status <- system2(file.path(R.home("bin"), "Rscript"), s)
  if (status != 0) stop("Pipeline failed at: ", s)
}

