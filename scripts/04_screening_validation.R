###############################################################################
# Stage 4: Phenotype-adaptive 4/5-item screening and locked LOCO validation
###############################################################################

options(stringsAsFactors = FALSE, warn = 1)
try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"), silent = TRUE)
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(pROC)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (basename(project_dir) == "scripts") project_dir <- dirname(project_dir)
workspace_dir <- normalizePath(dirname(project_dir), winslash = "/", mustWork = TRUE)
output_dir <- file.path(workspace_dir, "output")

set.seed(20260918)
d <- readRDS(file.path(workspace_dir, "eligible_with_phenotypes.rds"))
lca <- readRDS(file.path(workspace_dir, "lca_results.rds"))
item_vars <- c("scr_fatigue", "scr_srh", "scr_walk", "scr_adl", "scr_iadl")
item_labels <- c(scr_fatigue="Fatigue/effort", scr_srh="Self-rated health",
                 scr_walk="Walking", scr_adl="ADL", scr_iadl="IADL")

d$incident_adl <- ifelse(is.na(d$adl) | is.na(d$fu_adl) | d$adl >= 1,
                         NA, as.numeric(d$fu_adl >= 1))
d$frail_concurrent <- ifelse(is.na(d$fi_screen_excluded), NA, as.numeric(d$fi_screen_excluded >= 0.25))
d$frail_followup <- ifelse(is.na(d$fu_fi_core13), NA, as.numeric(d$fu_fi_core13 >= 0.25))
d$fi_progression <- ifelse(is.na(d$fi_core13) | is.na(d$fu_fi_core13), NA,
                           as.numeric(d$fu_fi_core13 - d$fi_core13 >= 0.05))

cat("=== Screening validation ===\n")
cat("Eligible:", nrow(d), "| five-item complete:", sum(complete.cases(d[, item_vars])), "\n")

fit_logit <- function(data, outcome, items) {
  x <- data[, c(outcome, items), drop = FALSE]
  names(x)[1] <- "y"
  x <- x[complete.cases(x), , drop = FALSE]
  if (nrow(x) < 100 || length(unique(x$y)) < 2) return(NULL)
  glm(reformulate(items, response = "y"), data = x, family = binomial)
}

safe_auc <- function(y, p) {
  ok <- !is.na(y) & !is.na(p)
  if (sum(ok) < 50 || length(unique(y[ok])) < 2) return(NA_real_)
  as.numeric(auc(roc(y[ok], p[ok], quiet = TRUE, direction = "<")))
}

choose_items <- function(data, outcome, size = 4, seed = 1) {
  combos <- combn(item_vars, size, simplify = FALSE)
  x <- data[complete.cases(data[, c(outcome, item_vars)]), , drop = FALSE]
  if (nrow(x) < 300 || length(unique(x[[outcome]])) < 2) return(combos[[1]])
  set.seed(seed)
  val_idx <- sample(seq_len(nrow(x)), size = max(100, floor(0.25 * nrow(x))))
  tr <- x[-val_idx, , drop = FALSE]
  va <- x[val_idx, , drop = FALSE]
  scores <- sapply(combos, function(items) {
    m <- fit_logit(tr, outcome, items)
    if (is.null(m)) return(NA_real_)
    safe_auc(va[[outcome]], predict(m, newdata = va, type = "response"))
  })
  if (all(is.na(scores))) combos[[1]] else combos[[which.max(scores)]]
}

screen_threshold <- function(y, p, sensitivity_target = 0.85) {
  ok <- !is.na(y) & !is.na(p)
  y <- y[ok]; p <- p[ok]
  if (sum(y == 1) < 20) return(0.5)
  unname(quantile(p[y == 1], probs = 1 - sensitivity_target, na.rm = TRUE, type = 1))
}

fit_universal <- function(train, outcome, size, seed) {
  items <- if (size == 5) item_vars else choose_items(train, outcome, size, seed)
  m <- fit_logit(train, outcome, items)
  p <- if (is.null(m)) rep(NA_real_, nrow(train)) else predict(m, newdata = train, type = "response")
  list(model = m, items = items, threshold = screen_threshold(train[[outcome]], p))
}

fit_adaptive <- function(train, outcome, size, phenotype, seed) {
  fits <- list()
  for (cl in sort(unique(phenotype[!is.na(phenotype)]))) {
    idx <- phenotype == cl & !is.na(phenotype)
    sub <- train[idx, , drop = FALSE]
    items <- if (size == 5) item_vars else choose_items(sub, outcome, size, seed + cl)
    m <- fit_logit(sub, outcome, items)
    if (is.null(m)) next
    p <- predict(m, newdata = sub, type = "response")
    fits[[as.character(cl)]] <- list(model = m, items = items,
                                     threshold = screen_threshold(sub[[outcome]], p), N = nrow(sub))
  }
  fits
}

predict_universal <- function(fit, newdata) {
  p <- predict(fit$model, newdata = newdata, type = "response")
  data.frame(prob = p, threshold = fit$threshold,
             screen_positive = as.numeric(p >= fit$threshold))
}

predict_adaptive <- function(fits, fallback, newdata, phenotype) {
  p <- rep(NA_real_, nrow(newdata)); thr <- rep(NA_real_, nrow(newdata))
  used <- rep("", nrow(newdata))
  for (i in seq_len(nrow(newdata))) {
    key <- as.character(phenotype[i])
    f <- fits[[key]]
    if (is.null(f)) {
      f <- fallback
      used[i] <- paste(f$items, collapse = "+")
    } else used[i] <- paste(f$items, collapse = "+")
    p[i] <- predict(f$model, newdata = newdata[i, , drop = FALSE], type = "response")
    thr[i] <- f$threshold
  }
  data.frame(prob = p, threshold = thr, screen_positive = as.numeric(p >= thr), items = used)
}

calibration_stats <- function(y, p) {
  ok <- !is.na(y) & !is.na(p)
  y <- y[ok]; p <- pmin(pmax(p[ok], 1e-6), 1 - 1e-6)
  if (length(y) < 100 || length(unique(y)) < 2) return(c(intercept=NA, slope=NA))
  lp <- qlogis(p)
  intercept <- tryCatch(coef(glm(y ~ 1 + offset(lp), family = binomial))[1], error = function(e) NA_real_)
  slope <- tryCatch(coef(glm(y ~ lp, family = binomial))[2], error = function(e) NA_real_)
  c(intercept = intercept, slope = slope)
}

calc_metrics <- function(x, outcome, outcome_label) {
  y <- x[[outcome]]
  ok <- !is.na(y) & !is.na(x$prob) & !is.na(x$screen_positive)
  y <- y[ok]; p <- x$prob[ok]; z <- x$screen_positive[ok]
  if (length(y) < 100 || length(unique(y)) < 2) return(NULL)
  roc_obj <- roc(y, p, quiet = TRUE, direction = "<")
  ci <- as.numeric(ci.auc(roc_obj, method = "delong"))
  tp <- sum(z == 1 & y == 1); tn <- sum(z == 0 & y == 0)
  fp <- sum(z == 1 & y == 0); fn <- sum(z == 0 & y == 1)
  cal <- calibration_stats(y, p)
  data.frame(
    outcome = outcome_label, N = length(y), events = sum(y == 1), prevalence = mean(y),
    AUROC = as.numeric(auc(roc_obj)), AUROC_low = ci[1], AUROC_high = ci[3],
    Brier = mean((p - y)^2), calibration_intercept = cal[1], calibration_slope = cal[2],
    sensitivity = tp / (tp + fn), specificity = tn / (tn + fp),
    PPV = tp / (tp + fp), NPV = tn / (tn + fn),
    false_negatives_per_1000 = 1000 * fn / length(y),
    screen_positive_pct = 100 * mean(z == 1), stringsAsFactors = FALSE)
}

item_sets <- data.frame()
predictions <- data.frame()
models_by_holdout <- list()

for (heldout in unique(d$cohort)) {
  cat("Held out:", heldout, "\n")
  train_idx <- d$cohort != heldout
  test_idx <- d$cohort == heldout
  complete_bank <- complete.cases(d[, item_vars])
  train_keep <- train_idx & complete_bank
  test_keep <- test_idx & complete_bank
  train <- d[train_keep, , drop = FALSE]
  test <- d[test_keep, , drop = FALSE]
  loco_all <- lca$loco_assignments[[heldout]]
  loc_class <- loco_all$phenotype[match(d$uid, loco_all$uid)]
  train_class <- loc_class[train_keep]
  test_class <- loc_class[test_keep]

  u4 <- fit_universal(train, "frail_concurrent", 4, 100 + match(heldout, unique(d$cohort)))
  u5 <- fit_universal(train, "frail_concurrent", 5, 200 + match(heldout, unique(d$cohort)))
  p4 <- fit_adaptive(train, "frail_concurrent", 4, train_class, 300 + match(heldout, unique(d$cohort)))
  p5 <- fit_adaptive(train, "frail_concurrent", 5, train_class, 400 + match(heldout, unique(d$cohort)))

  pu4 <- predict_universal(u4, test); pu4$model <- "Universal 4-item"
  pu5 <- predict_universal(u5, test); pu5$model <- "Universal 5-item"
  pp4 <- predict_adaptive(p4, u4, test, test_class); pp4$model <- "Adaptive 4-item"
  pp5 <- predict_adaptive(p5, u5, test, test_class); pp5$model <- "Adaptive 5-item"
  pred_list <- list(pu4, pp4, pu5, pp5)
  for (pp in pred_list) {
    predictions <- bind_rows(predictions, data.frame(
      uid = test$uid, cohort = test$cohort, held_out_cohort = heldout,
      phenotype = test_class, model = pp$model, prob = pp$prob,
      threshold = pp$threshold, screen_positive = pp$screen_positive,
      frail_concurrent = test$frail_concurrent,
      frail_followup = test$frail_followup,
      fi_progression = test$fi_progression,
      incident_adl = test$incident_adl,
      stringsAsFactors = FALSE))
  }

  item_sets <- bind_rows(item_sets,
    data.frame(held_out_cohort=heldout, model="Universal 4-item", phenotype="All", items=paste(u4$items, collapse=" + ")),
    data.frame(held_out_cohort=heldout, model="Universal 5-item", phenotype="All", items=paste(u5$items, collapse=" + ")))
  for (cl in names(p4)) item_sets <- bind_rows(item_sets, data.frame(
    held_out_cohort=heldout, model="Adaptive 4-item", phenotype=cl,
    items=paste(p4[[cl]]$items, collapse=" + ")))
  for (cl in names(p5)) item_sets <- bind_rows(item_sets, data.frame(
    held_out_cohort=heldout, model="Adaptive 5-item", phenotype=cl,
    items=paste(p5[[cl]]$items, collapse=" + ")))
  models_by_holdout[[heldout]] <- list(u4=u4, u5=u5, p4=p4, p5=p5)
}

outcomes <- c(frail_concurrent="Concurrent FI >=0.25",
              frail_followup="Follow-up FI >=0.25",
              fi_progression="FI progression >=0.05",
              incident_adl="Incident ADL limitation")
performance <- data.frame()
for (heldout in unique(predictions$held_out_cohort)) {
  for (model in unique(predictions$model)) {
    x <- predictions[predictions$held_out_cohort == heldout & predictions$model == model, , drop = FALSE]
    for (o in names(outcomes)) {
      met <- calc_metrics(x, o, outcomes[[o]])
      if (!is.null(met)) performance <- bind_rows(performance, cbind(held_out_cohort=heldout, model=model, met))
    }
  }
}
for (model in unique(predictions$model)) {
  x <- predictions[predictions$model == model, , drop = FALSE]
  for (o in names(outcomes)) {
    met <- calc_metrics(x, o, outcomes[[o]])
    if (!is.null(met)) performance <- bind_rows(performance, cbind(held_out_cohort="Pooled LOCO", model=model, met))
  }
}

# Decision-curve source data for concurrent and prospective frailty.
dca <- data.frame()
for (model in unique(predictions$model)) {
  x <- predictions[predictions$model == model, , drop = FALSE]
  for (o in c("frail_concurrent", "frail_followup")) {
    ok <- !is.na(x[[o]]) & !is.na(x$prob)
    y <- x[[o]][ok]; p <- x$prob[ok]; n <- length(y)
    for (pt in seq(0.05, 0.75, by = 0.025)) {
      tp <- sum(p >= pt & y == 1); fp <- sum(p >= pt & y == 0)
      nb <- tp / n - fp / n * pt / (1 - pt)
      dca <- bind_rows(dca, data.frame(model=model, outcome=unname(outcomes[o]), threshold=pt,
                                       net_benefit=nb, N=n))
    }
  }
}

# Final full-data models for deployable scores and intervention-relevance linkage.
full <- d[complete.cases(d[, item_vars]) & !is.na(d$frail_concurrent), , drop = FALSE]
full_class <- full$phenotype
final_u4 <- fit_universal(full, "frail_concurrent", 4, 999)
final_u5 <- fit_universal(full, "frail_concurrent", 5, 1000)
final_p4 <- fit_adaptive(full, "frail_concurrent", 4, full_class, 1001)
final_p5 <- fit_adaptive(full, "frail_concurrent", 5, full_class, 1002)
full_u4 <- predict_universal(final_u4, full)
full_p4 <- predict_adaptive(final_p4, final_u4, full, full_class)
full_scores <- data.frame(uid=full$uid, cohort=full$cohort, phenotype=full$phenotype,
                          universal4_prob=full_u4$prob, adaptive4_prob=full_p4$prob,
                          universal4_positive=full_u4$screen_positive,
                          adaptive4_positive=full_p4$screen_positive)

cate_path <- file.path(dirname(workspace_dir), "MS4", "output", "LOCO_CATE_PREDICTIONS_v2.csv")
intervention_metrics <- data.frame()
benefit_gradient <- data.frame()
if (file.exists(cate_path)) {
  cate <- read.csv(cate_path, stringsAsFactors = FALSE, check.names = FALSE)
  cate <- cate[cate$comparison == "A_vs_C", , drop = FALSE]
  linked <- merge(full_scores, cate, by.x = "uid", by.y = "person_id")
  if (nrow(linked) >= 100) {
    benefit_cut <- quantile(linked$predicted_benefit, 0.80, na.rm = TRUE)
    linked$high_benefit <- as.numeric(linked$predicted_benefit >= benefit_cut)
    intervention_metrics <- data.frame(
      N_linked = nrow(linked),
      AUROC_universal4 = safe_auc(linked$high_benefit, linked$universal4_prob),
      AUROC_adaptive4 = safe_auc(linked$high_benefit, linked$adaptive4_prob),
      Spearman_universal4 = cor(linked$predicted_benefit, linked$universal4_prob, method="spearman"),
      Spearman_adaptive4 = cor(linked$predicted_benefit, linked$adaptive4_prob, method="spearman"),
      high_benefit_capture_universal4 = mean(linked$universal4_positive[linked$high_benefit == 1]),
      high_benefit_capture_adaptive4 = mean(linked$adaptive4_positive[linked$high_benefit == 1])
    )
    for (score in c("universal4_prob", "adaptive4_prob")) {
      linked$stratum <- cut(linked[[score]], unique(quantile(linked[[score]], probs=seq(0,1,0.25), na.rm=TRUE)),
                            include.lowest=TRUE, labels=FALSE)
      benefit_gradient <- bind_rows(benefit_gradient, linked |>
        group_by(stratum) |>
        summarise(N=n(), mean_predicted_benefit=mean(predicted_benefit),
                  mean_cate=mean(cate), high_benefit_pct=100*mean(high_benefit), .groups="drop") |>
        mutate(model=ifelse(score=="universal4_prob","Universal 4-item","Adaptive 4-item")))
    }
    write.csv(linked, file.path(output_dir, "INTERVENTION_RELEVANCE_LINKED.csv"), row.names=FALSE)
  }
}

analysis_results <- list(
  item_vars=item_vars, item_labels=item_labels,
  models_by_holdout=models_by_holdout, final_models=list(u4=final_u4,u5=final_u5,p4=final_p4,p5=final_p5),
  item_sets=item_sets, performance=performance, predictions=predictions, dca=dca,
  full_scores=full_scores, intervention_metrics=intervention_metrics,
  benefit_gradient=benefit_gradient
)
saveRDS(analysis_results, file.path(workspace_dir, "analysis_results.rds"))
saveRDS(predictions, file.path(workspace_dir, "evaluation_data.rds"))
write.csv(item_sets, file.path(output_dir, "FINAL_ITEM_SETS.csv"), row.names=FALSE)
write.csv(performance, file.path(output_dir, "LOCO_SCREENING_PERFORMANCE.csv"), row.names=FALSE)
write.csv(predictions, file.path(output_dir, "LOCO_PREDICTIONS.csv"), row.names=FALSE)
write.csv(dca, file.path(output_dir, "DECISION_CURVE_SOURCE_DATA.csv"), row.names=FALSE)
write.csv(full_scores, file.path(output_dir, "FINAL_SCREENING_SCORES.csv"), row.names=FALSE)
write.csv(intervention_metrics, file.path(output_dir, "INTERVENTION_RELEVANCE_METRICS.csv"), row.names=FALSE)
write.csv(benefit_gradient, file.path(output_dir, "INTERVENTION_BENEFIT_GRADIENT.csv"), row.names=FALSE)

cat("\nPooled LOCO performance:\n")
print(performance[performance$held_out_cohort == "Pooled LOCO", ])
if (nrow(intervention_metrics)) print(intervention_metrics)
cat("=== Screening validation complete ===\n")
