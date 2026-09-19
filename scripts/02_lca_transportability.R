###############################################################################
# Stage 2: Clinical-social phenotype discovery and LOCO transportability
###############################################################################

options(stringsAsFactors = FALSE, warn = 1)
try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"), silent = TRUE)
suppressPackageStartupMessages({
  library(poLCA)
  library(dplyr)
  library(tidyr)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (basename(project_dir) == "scripts") project_dir <- dirname(project_dir)
workspace_dir <- normalizePath(dirname(project_dir), winslash = "/", mustWork = TRUE)
output_dir <- file.path(workspace_dir, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260918)
d <- readRDS(file.path(workspace_dir, "eligible_with_bmi.rds"))
context_vars <- c("hypertension", "diabetes", "heart_disease", "cancer", "stroke",
                  "arthritis", "not_married", "low_education")

cat("=== LCA phenotype discovery ===\n")
cat("Eligible:", nrow(d), "\n")

fit_lca <- function(data, k, nrep = 3, seed = 1) {
  x <- data[, context_vars, drop = FALSE]
  x <- x[complete.cases(x), , drop = FALSE]
  for (v in context_vars) x[[v]] <- factor(x[[v]], levels = c(0, 1))
  f <- as.formula(paste0("cbind(", paste(context_vars, collapse = ","), ") ~ 1"))
  set.seed(seed)
  poLCA(f, data = x, nclass = k, nrep = nrep, maxiter = 1500,
        tol = 1e-8, verbose = FALSE, calc.se = FALSE)
}

model_stats <- function(m, n) {
  post <- m$posterior
  k <- ncol(post)
  entropy <- 1 + sum(post * log(post + 1e-12)) / (nrow(post) * log(k))
  assigned <- max.col(post, ties.method = "first")
  avg_post <- sapply(seq_len(k), function(i) mean(post[assigned == i, i], na.rm = TRUE))
  npar <- if (!is.null(m$npar)) m$npar else length(unlist(coef(m)))
  data.frame(AIC = m$aic, BIC = m$bic,
             aBIC = -2 * m$llik + log((n + 2) / 24) * npar,
             logLik = m$llik, entropy = entropy,
             min_class_pct = 100 * min(tabulate(assigned, nbins = k)) / nrow(post),
             mean_assigned_posterior = mean(avg_post), npar = npar)
}

profile_matrix <- function(m) {
  t(vapply(context_vars, function(v) m$probs[[v]][, 2], numeric(nrow(m$probs[[context_vars[1]]])))) |>
    t()
}

predict_lca <- function(data, m) {
  prof <- profile_matrix(m)
  k <- nrow(prof)
  logp <- matrix(log(pmax(m$P, 1e-12)), nrow = nrow(data), ncol = k, byrow = TRUE)
  for (j in seq_along(context_vars)) {
    x <- data[[context_vars[j]]]
    obs <- which(!is.na(x) & x %in% c(0, 1))
    if (!length(obs)) next
    for (cl in seq_len(k)) {
      pr <- pmin(pmax(prof[cl, j], 1e-8), 1 - 1e-8)
      logp[obs, cl] <- logp[obs, cl] + x[obs] * log(pr) + (1 - x[obs]) * log(1 - pr)
    }
  }
  row_max <- apply(logp, 1, max)
  p <- exp(logp - row_max)
  p <- p / rowSums(p)
  list(class = max.col(p, ties.method = "first"), posterior = p,
       max_posterior = apply(p, 1, max))
}

comparison <- data.frame()
models <- list()
complete_n <- sum(complete.cases(d[, context_vars]))
for (k in 3:7) {
  cat("Fitting", k, "classes...\n")
  t0 <- proc.time()[3]
  m <- fit_lca(d, k, nrep = 3, seed = 20260918 + k)
  st <- model_stats(m, complete_n)
  comparison <- bind_rows(comparison, cbind(n_class = k, st, elapsed_sec = proc.time()[3] - t0))
  models[[as.character(k)]] <- m
  cat(sprintf("  BIC %.1f | entropy %.3f | min class %.1f%%\n", st$BIC, st$entropy, st$min_class_pct))
}

valid <- comparison$n_class <= 6 & comparison$entropy >= 0.60 & comparison$min_class_pct >= 5
if (!any(valid)) valid <- comparison$min_class_pct >= 3
candidate <- comparison[valid, , drop = FALSE]
selected_k <- candidate$n_class[which.min(candidate$BIC)]
cat("Selected class count:", selected_k, "\n")
cat("Refitting selected model with 10 random starts...\n")
primary_model <- fit_lca(d, selected_k, nrep = 10, seed = 20261001)
models[[as.character(selected_k)]] <- primary_model
primary_profile <- profile_matrix(primary_model)
colnames(primary_profile) <- context_vars

name_classes <- function(profile) {
  disease_cols <- match(c("hypertension","diabetes","heart_disease","cancer","stroke","arthritis"), colnames(profile))
  social_cols <- match(c("not_married","low_education"), colnames(profile))
  out <- character(nrow(profile))
  for (i in seq_len(nrow(profile))) {
    p <- profile[i, ]
    disease_mean <- mean(p[disease_cols])
    if (disease_mean >= 0.48) label <- "Multisystem burden"
    else if (p["diabetes"] >= 0.55) label <- "Metabolic"
    else if (max(p[c("heart_disease","stroke")]) >= 0.45) label <- "Cardiovascular"
    else if (p["arthritis"] >= 0.70) label <- "Musculoskeletal"
    else if (mean(p[social_cols]) >= 0.65) label <- "Social vulnerability"
    else if (p["cancer"] >= 0.35) label <- "Cancer-enriched"
    else if (p["hypertension"] >= 0.75) label <- "Hypertension-dominant"
    else label <- "Lower-burden mixed"
    out[i] <- paste0(label, " (P", i, ")")
  }
  out
}

class_names <- name_classes(primary_profile)
pred <- predict_lca(d, primary_model)
d$phenotype <- pred$class
d$phenotype_name <- factor(class_names[d$phenotype], levels = class_names)
d$max_posterior <- pred$max_posterior
for (k in seq_len(selected_k)) d[[paste0("phenotype_p", k)]] <- pred$posterior[, k]

prob_table <- data.frame(phenotype = class_names, primary_profile, check.names = FALSE)
class_summary <- d |>
  group_by(phenotype, phenotype_name) |>
  summarise(N = n(), prevalence_pct = 100 * n() / nrow(d),
            mean_max_posterior = mean(max_posterior),
            pct_posterior_ge_070 = 100 * mean(max_posterior >= 0.70), .groups = "drop")
cohort_prevalence <- d |>
  count(cohort, phenotype, phenotype_name, name = "N") |>
  group_by(cohort) |>
  mutate(prevalence_pct = 100 * N / sum(N)) |>
  ungroup()

permutations <- function(v) {
  if (length(v) == 1) return(matrix(v, nrow = 1))
  do.call(rbind, lapply(seq_along(v), function(i) {
    cbind(v[i], permutations(v[-i]))
  }))
}

align_classes <- function(source_profile, target_profile) {
  perms <- permutations(seq_len(nrow(target_profile)))
  loss <- apply(perms, 1, function(p) sum((source_profile - target_profile[p, , drop = FALSE])^2))
  perms[which.min(loss), ]
}

cat("Running leave-one-cohort-out phenotype transportability...\n")
loco_models <- list()
loco_assignments <- list()
loco_transport <- data.frame()
for (heldout in unique(d$cohort)) {
  cat("  Held out:", heldout, "\n")
  train <- d[d$cohort != heldout, , drop = FALSE]
  m <- fit_lca(train, selected_k, nrep = 3, seed = 20262000 + match(heldout, unique(d$cohort)))
  source_profile <- profile_matrix(m)
  colnames(source_profile) <- context_vars
  map_raw_to_primary <- align_classes(source_profile, primary_profile)
  pa <- predict_lca(d, m)
  aligned_class <- map_raw_to_primary[pa$class]
  loco_models[[heldout]] <- list(model = m, raw_to_primary = map_raw_to_primary,
                                 aligned_profile = source_profile[order(map_raw_to_primary), , drop = FALSE])
  loco_assignments[[heldout]] <- data.frame(uid = d$uid, cohort = d$cohort,
                                            phenotype = aligned_class,
                                            max_posterior = pa$max_posterior)
  aligned_profile <- source_profile[match(seq_len(selected_k), map_raw_to_primary), , drop = FALSE]
  class_cor <- sapply(seq_len(selected_k), function(i) cor(aligned_profile[i, ], primary_profile[i, ]))
  test_idx <- d$cohort == heldout
  loco_transport <- bind_rows(loco_transport, data.frame(
    held_out_cohort = heldout,
    N_train = sum(!test_idx), N_test = sum(test_idx),
    profile_correlation = mean(class_cor, na.rm = TRUE),
    profile_rmse = sqrt(mean((aligned_profile - primary_profile)^2)),
    test_mean_max_posterior = mean(pa$max_posterior[test_idx]),
    test_pct_posterior_ge_070 = 100 * mean(pa$max_posterior[test_idx] >= 0.70)
  ))
}

lca_results <- list(
  context_vars = context_vars, comparison = comparison, selected_k = selected_k,
  primary_model = primary_model, primary_profile = primary_profile,
  class_names = class_names, class_summary = class_summary,
  cohort_prevalence = cohort_prevalence,
  loco_models = loco_models, loco_assignments = loco_assignments,
  loco_transport = loco_transport
)

saveRDS(lca_results, file.path(workspace_dir, "lca_results.rds"))
saveRDS(lca_results, file.path(workspace_dir, "lca_validation_results.rds"))
saveRDS(d, file.path(workspace_dir, "eligible_with_phenotypes.rds"))
write.csv(comparison, file.path(output_dir, "LCA_MODEL_COMPARISON.csv"), row.names = FALSE)
write.csv(prob_table, file.path(output_dir, "LCA_CONDITIONAL_PROBABILITIES.csv"), row.names = FALSE)
write.csv(class_summary, file.path(output_dir, "LCA_CLASS_SUMMARY.csv"), row.names = FALSE)
write.csv(cohort_prevalence, file.path(output_dir, "LCA_COHORT_PREVALENCE.csv"), row.names = FALSE)
write.csv(loco_transport, file.path(output_dir, "LCA_LOCO_TRANSPORTABILITY.csv"), row.names = FALSE)

cat("\nSelected:", selected_k, "classes\n")
print(class_summary)
print(loco_transport)
cat("=== LCA complete ===\n")
