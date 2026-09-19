###############################################################################
# Stage 3: IRT item information and DIF magnitude
###############################################################################

options(stringsAsFactors = FALSE, warn = 1)
try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"), silent = TRUE)
suppressPackageStartupMessages({
  library(ltm)
  library(dplyr)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (basename(project_dir) == "scripts") project_dir <- dirname(project_dir)
workspace_dir <- normalizePath(dirname(project_dir), winslash = "/", mustWork = TRUE)
output_dir <- file.path(workspace_dir, "output")

d <- readRDS(file.path(workspace_dir, "eligible_with_phenotypes.rds"))
item_vars <- c("scr_fatigue", "scr_srh", "scr_walk", "scr_adl", "scr_iadl")
item_labels <- c(scr_fatigue="Fatigue/effort", scr_srh="Fair/poor self-rated health",
                 scr_walk="Walking difficulty", scr_adl="ADL limitation",
                 scr_iadl="IADL limitation")

cat("=== IRT and DIF ===\n")
irt_cc <- d[complete.cases(d[, item_vars]), item_vars, drop = FALSE]
cat("Complete five-item IRT sample:", nrow(irt_cc), "\n")
if (nrow(irt_cc) < 500) stop("Insufficient complete five-item sample for IRT")

set.seed(20260918)
# Randomized starts avoid the boundary solution seen with the default start for
# the high-prevalence self-rated-health item.
irt_model <- ltm(irt_cc ~ z1, IRT.param = TRUE, start.val = "random",
                 control = list(iter.qN = 200))
params <- coef(irt_model)
param_table <- data.frame(
  item = rownames(params), label = unname(item_labels[rownames(params)]),
  difficulty = params[, "Dffclt"], discrimination = params[, "Dscrmn"],
  stringsAsFactors = FALSE
)

theta_grid <- seq(-3, 3, by = 0.05)
item_information <- bind_rows(lapply(seq_len(nrow(param_table)), function(i) {
  a <- param_table$discrimination[i]
  b <- param_table$difficulty[i]
  p <- 1 / (1 + exp(-a * (theta_grid - b)))
  data.frame(item = param_table$item[i], label = param_table$label[i],
             theta = theta_grid, information = a^2 * p * (1 - p))
}))
param_table$information_theta0 <- sapply(seq_len(nrow(param_table)), function(i) {
  a <- param_table$discrimination[i]; b <- param_table$difficulty[i]
  p <- 1 / (1 + exp(-a * (0 - b)))
  a^2 * p * (1 - p)
})

pseudo_r2 <- function(model) {
  ll <- as.numeric(logLik(model))
  null <- glm(model$y ~ 1, family = binomial)
  1 - ll / as.numeric(logLik(null))
}

run_dif <- function(data, group_var) {
  ans <- list()
  for (item in item_vars) {
    others <- setdiff(item_vars, item)
    x <- data[, c(item, others, group_var), drop = FALSE]
    x$rest_score <- rowMeans(x[, others, drop = FALSE], na.rm = TRUE)
    x$rest_n <- rowSums(!is.na(x[, others, drop = FALSE]))
    x <- x[!is.na(x[[item]]) & x$rest_n >= 3 & !is.na(x[[group_var]]), , drop = FALSE]
    x$group <- factor(x[[group_var]])
    x$rest_z <- as.numeric(scale(x$rest_score))
    if (nrow(x) < 500 || nlevels(x$group) < 2) next
    m0 <- glm(x[[item]] ~ rest_z, data = x, family = binomial)
    m1 <- glm(x[[item]] ~ rest_z + group, data = x, family = binomial)
    m2 <- glm(x[[item]] ~ rest_z * group, data = x, family = binomial)
    ll0 <- as.numeric(logLik(m0)); ll1 <- as.numeric(logLik(m1)); ll2 <- as.numeric(logLik(m2))
    r20 <- 1 - ll0 / as.numeric(logLik(glm(x[[item]] ~ 1, data = x, family = binomial)))
    r21 <- 1 - ll1 / as.numeric(logLik(glm(x[[item]] ~ 1, data = x, family = binomial)))
    r22 <- 1 - ll2 / as.numeric(logLik(glm(x[[item]] ~ 1, data = x, family = binomial)))
    p_uniform <- anova(m0, m1, test = "Chisq")$`Pr(>Chi)`[2]
    p_nonuniform <- anova(m1, m2, test = "Chisq")$`Pr(>Chi)`[2]
    delta_total <- r22 - r20
    ans[[item]] <- data.frame(
      grouping = group_var, item = item, label = unname(item_labels[item]), N = nrow(x),
      delta_r2_uniform = r21 - r20, delta_r2_nonuniform = r22 - r21,
      delta_r2_total = delta_total, p_uniform = p_uniform, p_nonuniform = p_nonuniform,
      magnitude = cut(delta_total, breaks = c(-Inf, 0.02, 0.05, Inf),
                      labels = c("Negligible", "Moderate", "Large"), right = FALSE),
      stringsAsFactors = FALSE
    )
  }
  out <- bind_rows(ans)
  out$q_uniform <- p.adjust(out$p_uniform, method = "BH")
  out$q_nonuniform <- p.adjust(out$p_nonuniform, method = "BH")
  out
}

dif_phenotype <- run_dif(d, "phenotype")
dif_cohort <- run_dif(d, "cohort")
dif_results <- bind_rows(dif_phenotype, dif_cohort)

result <- list(
  item_vars = item_vars, item_labels = item_labels, irt_model = irt_model,
  parameters = param_table, item_information = item_information,
  dif_results = dif_results, N_irt = nrow(irt_cc)
)
saveRDS(result, file.path(workspace_dir, "irt_dif_results.rds"))
write.csv(param_table, file.path(output_dir, "IRT_ITEM_PARAMETERS.csv"), row.names = FALSE)
write.csv(item_information, file.path(output_dir, "IRT_ITEM_INFORMATION_CURVES.csv"), row.names = FALSE)
write.csv(dif_results, file.path(output_dir, "DIF_MAGNITUDE.csv"), row.names = FALSE)

print(param_table)
print(dif_results)
cat("=== IRT/DIF complete ===\n")
