###############################################################################
# Stage 1: Recover harmonized same-wave baseline and follow-up data
###############################################################################

options(stringsAsFactors = FALSE, warn = 1)
try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"), silent = TRUE)

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
})

project_dir <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
if (basename(project_dir) == "scripts") project_dir <- dirname(project_dir)
workspace_dir <- normalizePath(dirname(project_dir), winslash = "/", mustWork = TRUE)
data_root <- normalizePath(dirname(workspace_dir), winslash = "/", mustWork = TRUE)
output_dir <- file.path(workspace_dir, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Same-wave data recovery ===\n")
cat("Project:", project_dir, "\n")
cat("Workspace:", workspace_dir, "\n")

clean_num <- function(x, valid = NULL) {
  x <- suppressWarnings(as.numeric(x))
  x[is.finite(x) & x < 0] <- NA_real_
  if (!is.null(valid)) x[!is.na(x) & !(x %in% valid)] <- NA_real_
  x
}

pull_num <- function(d, name, valid = NULL) {
  if (is.null(name) || is.na(name) || !name %in% names(d)) return(rep(NA_real_, nrow(d)))
  clean_num(d[[name]], valid)
}

pull_id <- function(d, name) {
  x <- d[[name]]
  if (is.numeric(x)) format(x, scientific = FALSE, trim = TRUE) else as.character(x)
}

recode_married <- function(x) {
  x <- clean_num(x)
  ifelse(is.na(x), NA_real_, ifelse(x %in% c(1, 2, 3), 1, ifelse(x %in% c(4, 5, 6, 7, 8), 0, NA)))
}

recode_fatigue <- function(x, type) {
  x <- clean_num(x)
  if (type == "binary") return(ifelse(is.na(x), NA_real_, as.numeric(x >= 1)))
  if (type == "effort_binary") return(ifelse(is.na(x), NA_real_, as.numeric(x == 1)))
  if (type == "effort_ordinal") return(ifelse(is.na(x), NA_real_, as.numeric(x >= 3)))
  stop("Unknown fatigue type: ", type)
}

recode_srh <- function(x) clean_num(x, 1:5)
recode_any <- function(x) {
  x <- clean_num(x)
  ifelse(is.na(x), NA_real_, as.numeric(x >= 1))
}

normalize_height <- function(x) {
  x <- clean_num(x)
  cm_idx <- which(!is.na(x) & x > 3 & x <= 250)
  x[cm_idx] <- x[cm_idx] / 100
  x[which(!is.na(x) & (x < 1.0 | x > 2.3))] <- NA_real_
  x
}

normalize_weight <- function(x) {
  x <- clean_num(x)
  x[which(!is.na(x) & (x < 25 | x > 250))] <- NA_real_
  x
}

make_standard <- function(d, cohort, country, id_col, gender_col, education_col,
                          baseline_wave, followup_wave, baseline_map, followup_map,
                          fatigue_type) {
  get_map <- function(map, key) if (key %in% names(map)) unname(map[[key]]) else NA_character_
  out <- data.frame(
    id = pull_id(d, id_col),
    cohort = cohort,
    country = country,
    baseline_wave = baseline_wave,
    followup_wave = followup_wave,
    gender = pull_num(d, gender_col, c(1, 2)),
    education = pull_num(d, education_col, c(1, 2, 3)),
    stringsAsFactors = FALSE
  )

  extract_wave <- function(map, prefix = "") {
    vals <- list(
      age = pull_num(d, get_map(map, "age")),
      iwstat = pull_num(d, get_map(map, "iwstat")),
      hypertension = pull_num(d, get_map(map, "hypertension"), c(0, 1)),
      diabetes = pull_num(d, get_map(map, "diabetes"), c(0, 1)),
      heart_disease = pull_num(d, get_map(map, "heart_disease"), c(0, 1)),
      cancer = pull_num(d, get_map(map, "cancer"), c(0, 1)),
      lung_disease = pull_num(d, get_map(map, "lung_disease"), c(0, 1)),
      stroke = pull_num(d, get_map(map, "stroke"), c(0, 1)),
      arthritis = pull_num(d, get_map(map, "arthritis"), c(0, 1)),
      married = recode_married(pull_num(d, get_map(map, "married"))),
      srh = recode_srh(pull_num(d, get_map(map, "srh"))),
      walk_difficulty = recode_any(pull_num(d, get_map(map, "walk_difficulty"))),
      fatigue = recode_fatigue(pull_num(d, get_map(map, "fatigue")), fatigue_type),
      fall = recode_any(pull_num(d, get_map(map, "fall"))),
      adl = clean_num(pull_num(d, get_map(map, "adl"))),
      iadl = clean_num(pull_num(d, get_map(map, "iadl"))),
      sight = recode_srh(pull_num(d, get_map(map, "sight"))),
      grip = clean_num(pull_num(d, get_map(map, "grip"))),
      height_m = normalize_height(pull_num(d, get_map(map, "height"))),
      weight_kg = normalize_weight(pull_num(d, get_map(map, "weight")))
    )
    names(vals) <- paste0(prefix, names(vals))
    as.data.frame(vals, stringsAsFactors = FALSE)
  }

  bind_cols(out, extract_wave(baseline_map), extract_wave(followup_map, "fu_"))
}

read_selected <- function(path, cols) {
  cols <- unique(cols[!is.na(cols) & nzchar(cols)])
  read_dta(path, col_select = any_of(cols))
}

map_cols <- function(...) unlist(list(...), use.names = TRUE)

# HRS: clinical context and mobility/SRH come from RAND; performance items from H_HRS.
hrs_main_path <- file.path(data_root, "HRS_USA", "H_HRS_d.dta")
hrs_rand_path <- file.path(data_root, "HRS_USA", "randhrs1992_2022v1.dta")
hrs_b <- map_cols(age="r14agey_b", iwstat="r14iwstat", hypertension="r14hibpe",
                  diabetes="r14diabe", heart_disease="r14hearte", cancer="r14cancre",
                  lung_disease="r14lunge", stroke="r14stroke", arthritis="r14arthre",
                  married="r14mstat", srh="r14shlt", walk_difficulty="r14walkra",
                  fatigue="r14effort", fall="r14fall", adl="r14adl5a",
                  iadl="r14iadl5a", sight="r14sight", grip="r14lgrip1",
                  height="r14mheight", weight="r14mweight")
hrs_f <- map_cols(age="r15agey_b", iwstat="r15iwstat", hypertension="r15hibpe",
                  diabetes="r15diabe", heart_disease="r15hearte", cancer="r15cancre",
                  lung_disease="r15lunge", stroke="r15stroke", arthritis="r15arthre",
                  married="r15mstat", srh="r15shlt", walk_difficulty="r15walkra",
                  fatigue="r15effort", fall="r15fall", adl="r15adl5a",
                  iadl="r15iadl5a", sight="r15sight")
hrs_rand_cols <- c("hhidpn", "ragender", setdiff(c(hrs_b, hrs_f), c("r14fall","r14sight","r14lgrip1","r14mheight","r14mweight","r15fall","r15sight")))
hrs_main_cols <- c("hhidpn", "raeducl", intersect(c(hrs_b, hrs_f), c("r14fall","r14sight","r14lgrip1","r14mheight","r14mweight","r15fall","r15sight")))
cat("Reading HRS selected columns...\n")
hrs_rand <- read_selected(hrs_rand_path, hrs_rand_cols)
hrs_main <- read_selected(hrs_main_path, hrs_main_cols)
hrs <- merge(hrs_rand, hrs_main, by = "hhidpn", all = TRUE, sort = FALSE)
hrs_std <- make_standard(hrs, "HRS", "USA", "hhidpn", "ragender", "raeducl",
                         14, 15, hrs_b, hrs_f, "effort_binary")

# ELSA waves 9 -> 10.
elsa_path <- file.path(data_root, "ELSA_UK", "UKDA-5050-stata", "stata", "stata13_se", "gh_elsa_h.dta")
elsa_b <- map_cols(age="r9agey", iwstat="r9iwstat", hypertension="r9hibpe",
                   diabetes="r9diabe", heart_disease="r9hearte", cancer="r9cancre",
                   lung_disease="r9lunge", stroke="r9stroke", arthritis="r9arthre",
                   married="r9mstat", srh="r9shlt", walk_difficulty="r9walkra",
                   fatigue="r9effort", fall="r9fall", adl="r9adlfive",
                   iadl="r9iadlfour", sight="r9sight", weight="r9mweight")
elsa_f <- sub("r9", "r10", elsa_b, fixed = TRUE)
cat("Reading ELSA selected columns...\n")
elsa <- read_selected(elsa_path, c("idauniq", "ragender", "raeducl", elsa_b, elsa_f))
elsa_std <- make_standard(elsa, "ELSA", "England", "idauniq", "ragender", "raeducl",
                          9, 10, elsa_b, elsa_f, "effort_binary")

# SHARE waves 2 -> 4; wave 2 is the latest SHARE wave with the common effort item.
share_path <- file.path(data_root, "SHARE_Europe", "GH_SHARE_g.dta")
share_b <- map_cols(age="r2agey", iwstat="r2iwstat", hypertension="r2hibpe",
                    diabetes="r2diabe", heart_disease="r2hearte", cancer="r2cancre",
                    lung_disease="r2lunge", stroke="r2stroke", arthritis="r2arthre",
                    married="r2mstat", srh="r2shlt", walk_difficulty="r2walkra",
                    fatigue="r2effort", fall="r2fall_s", adl="r2adlfive",
                    iadl="r2iadlfour", grip="r2lgrip1", height="r2height", weight="r2weight")
share_f <- map_cols(age="r4agey", iwstat="r4iwstat", hypertension="r4hibpe",
                    diabetes="r4diabe", heart_disease="r4hearte", cancer="r4cancre",
                    lung_disease="r4lunge", stroke="r4stroke", arthritis="r4arthre",
                    married="r4mstat", srh="r4shlt", walk_difficulty="r4walkra",
                    fall="r4fall_s", adl="r4adlfive", iadl="r4iadlfour",
                    grip="r4lgrip1", height="r4height", weight="r4weight")
cat("Reading SHARE selected columns...\n")
share <- read_selected(share_path, c("mergeid", "ragender", "raeducl", share_b, share_f))
share_std <- make_standard(share, "SHARE", "Europe", "mergeid", "ragender", "raeducl",
                           2, 4, share_b, share_f, "effort_binary")

# CHARLS waves 3 -> 4.
charls_path <- file.path(data_root, "CHARLS_China", "H_CHARLS_D_Data.dta")
charls_b <- map_cols(age="r3agey", iwstat="r3iwstat", hypertension="r3hibpe",
                     diabetes="r3diabe", heart_disease="r3hearte", cancer="r3cancre",
                     lung_disease="r3lunge", stroke="r3stroke", arthritis="r3arthre",
                     married="r3mstat", srh="r3shlt", walk_difficulty="r3walk100a",
                     fatigue="r3effortl", adl="r3adlfive", iadl="r3iadla",
                     grip="r3lgrip1", height="r3mheight", weight="r3mweight")
charls_f <- map_cols(age="r4agey", iwstat="r4iwstat", hypertension="r4hibpe",
                     diabetes="r4diabe", heart_disease="r4hearte", cancer="r4cancre",
                     lung_disease="r4lunge", stroke="r4stroke", arthritis="r4arthre",
                     married="r4mstat", srh="r4shlta", walk_difficulty="r4walk100a",
                     fatigue="r4effortl", adl="r4adlfive", iadl="r4iadla")
cat("Reading CHARLS selected columns...\n")
charls <- read_selected(charls_path, c("ID", "ragender", "raeducl", charls_b, charls_f))
charls_std <- make_standard(charls, "CHARLS", "China", "ID", "ragender", "raeducl",
                            3, 4, charls_b, charls_f, "effort_ordinal")

# MHAS waves 5 -> 6.
mhas_path <- file.path(data_root, "MHAS_Mexico", "H_MHAS_d.dta")
mhas_b <- map_cols(age="r5agey", iwstat="r5iwstat", hypertension="r5hibpe",
                   diabetes="r5diabe", heart_disease="r5hearte", cancer="r5cancre",
                   stroke="r5stroke", arthritis="r5arthre", married="r5mstat",
                   srh="r5shlt", walk_difficulty="r5walkra", fatigue="r5fatigue",
                   fall="r5fall", adl="r5adlfive", iadl="r5iadlfour", sight="r5sight",
                   height="r5height", weight="r5weight")
mhas_f <- sub("r5", "r6", mhas_b, fixed = TRUE)
cat("Reading MHAS selected columns...\n")
mhas <- read_selected(mhas_path, c("unhhidnp", "ragender", "raeducl", mhas_b, mhas_f))
mhas_std <- make_standard(mhas, "MHAS", "Mexico", "unhhidnp", "ragender", "raeducl",
                          5, 6, mhas_b, mhas_f, "binary")

combined <- bind_rows(hrs_std, elsa_std, share_std, charls_std, mhas_std)
combined$uid <- paste(combined$cohort, combined$id, sep = ":")
combined$sex <- factor(ifelse(combined$gender == 1, "Male", ifelse(combined$gender == 2, "Female", NA)),
                       levels = c("Male", "Female"))
combined$low_education <- ifelse(is.na(combined$education), NA_real_, as.numeric(combined$education == 1))
combined$not_married <- ifelse(is.na(combined$married), NA_real_, 1 - combined$married)
combined$bmi_continuous <- combined$weight_kg / combined$height_m^2
combined$bmi_continuous[which(!is.na(combined$bmi_continuous) & (combined$bmi_continuous < 10 | combined$bmi_continuous > 80))] <- NA_real_
combined$fu_bmi_continuous <- combined$fu_weight_kg / combined$fu_height_m^2
combined$fu_bmi_continuous[which(!is.na(combined$fu_bmi_continuous) & (combined$fu_bmi_continuous < 10 | combined$fu_bmi_continuous > 80))] <- NA_real_

# Follow-up education is treated as time-invariant; marital status remains wave-specific.
combined$fu_low_education <- combined$low_education
combined$fu_not_married <- ifelse(is.na(combined$fu_married), NA_real_, 1 - combined$fu_married)

# Grip weakness: cohort- and sex-specific lowest quintile among baseline respondents.
combined$low_grip <- NA_real_
for (coh in unique(combined$cohort)) {
  for (sx in levels(combined$sex)) {
    idx <- combined$cohort == coh & combined$sex == sx & !is.na(combined$grip)
    if (sum(idx) >= 100) {
      cut <- unname(quantile(combined$grip[idx], 0.20, na.rm = TRUE, type = 7))
      combined$low_grip[idx] <- as.numeric(combined$grip[idx] <= cut)
    }
  }
}

make_deficits <- function(d, prefix = "") {
  g <- function(v) d[[paste0(prefix, v)]]
  defs <- data.frame(
    hypertension = g("hypertension"), diabetes = g("diabetes"),
    heart = g("heart_disease"), cancer = g("cancer"), stroke = g("stroke"),
    arthritis = g("arthritis"), not_married = g("not_married"),
    low_education = g("low_education"),
    srh = ifelse(is.na(g("srh")), NA, (g("srh") - 1) / 4),
    walk = g("walk_difficulty"), fatigue = g("fatigue"),
    adl = ifelse(is.na(g("adl")), NA, as.numeric(g("adl") >= 1)),
    iadl = ifelse(is.na(g("iadl")), NA, as.numeric(g("iadl") >= 1)),
    stringsAsFactors = FALSE
  )
  defs
}

calc_fi <- function(defs, min_fraction = 0.80) {
  n_avail <- rowSums(!is.na(defs))
  ans <- rowSums(defs, na.rm = TRUE) / n_avail
  ans[n_avail < ceiling(ncol(defs) * min_fraction)] <- NA_real_
  ans
}

base_defs <- make_deficits(combined)
fu_defs <- make_deficits(combined, "fu_")
combined$fi_core13 <- calc_fi(base_defs)
combined$fu_fi_core13 <- calc_fi(fu_defs)

bmi_deficit <- function(x) ifelse(is.na(x), NA_real_,
  ifelse(x < 18.5, 0.75, ifelse(x < 25, 0,
  ifelse(x < 30, 0.25, ifelse(x < 35, 0.50, ifelse(x < 40, 0.75, 1))))))
base_extra <- data.frame(
  fall = combined$fall, low_grip = combined$low_grip,
  bmi = bmi_deficit(combined$bmi_continuous),
  sight = ifelse(is.na(combined$sight), NA_real_, (combined$sight - 1) / 4)
)
fu_extra <- data.frame(
  fall = combined$fu_fall, low_grip = NA_real_,
  bmi = bmi_deficit(combined$fu_bmi_continuous),
  sight = ifelse(is.na(combined$fu_sight), NA_real_, (combined$fu_sight - 1) / 4)
)
# Primary concurrent criterion excludes all five candidate screening questions.
# At least 9 of 12 non-screen deficits must be observed.
screen_excluded_base <- cbind(base_defs[, 1:8, drop = FALSE], base_extra)
screen_excluded_fu <- cbind(fu_defs[, 1:8, drop = FALSE], fu_extra)
combined$fi_screen_excluded <- rowSums(screen_excluded_base, na.rm = TRUE) / rowSums(!is.na(screen_excluded_base))
combined$fi_screen_excluded[rowSums(!is.na(screen_excluded_base)) < 9] <- NA_real_
combined$fu_fi_screen_excluded <- rowSums(screen_excluded_fu, na.rm = TRUE) / rowSums(!is.na(screen_excluded_fu))
combined$fu_fi_screen_excluded[rowSums(!is.na(screen_excluded_fu)) < 9] <- NA_real_
combined$fi_full <- calc_fi(cbind(base_defs, base_extra))
combined$fu_fi_full <- calc_fi(cbind(fu_defs, fu_extra))

combined$age_65plus <- !is.na(combined$age) & combined$age >= 65
chronic_vars <- c("hypertension", "diabetes", "heart_disease", "cancer", "stroke", "arthritis")
combined$has_chronic <- rowSums(combined[, chronic_vars, drop = FALSE] == 1, na.rm = TRUE) >= 1
combined$has_social <- !is.na(combined$married) | !is.na(combined$education)
combined$eligible <- combined$age_65plus & combined$has_chronic & combined$has_social
eligible <- combined[combined$eligible %in% TRUE, , drop = FALSE]
eligible$age_group <- cut(eligible$age, breaks = c(65, 75, 85, Inf), right = FALSE,
                          labels = c("65-74", "75-84", "85+"))
eligible$frail_concurrent <- ifelse(is.na(eligible$fi_screen_excluded), NA, eligible$fi_screen_excluded >= 0.25)
eligible$frail_followup <- ifelse(is.na(eligible$fu_fi_core13), NA, eligible$fu_fi_core13 >= 0.25)
eligible$fi_progression <- ifelse(is.na(eligible$fi_core13) | is.na(eligible$fu_fi_core13), NA,
                                  eligible$fu_fi_core13 - eligible$fi_core13 >= 0.05)
eligible$incident_adl <- ifelse(is.na(eligible$adl) | is.na(eligible$fu_adl) | eligible$adl >= 1,
                                NA, as.numeric(eligible$fu_adl >= 1))

# Five-item harmonized screening bank.
eligible$scr_fatigue <- eligible$fatigue
eligible$scr_srh <- ifelse(is.na(eligible$srh), NA, as.numeric(eligible$srh >= 3))
eligible$scr_walk <- eligible$walk_difficulty
eligible$scr_adl <- ifelse(is.na(eligible$adl), NA, as.numeric(eligible$adl >= 1))
eligible$scr_iadl <- ifelse(is.na(eligible$iadl), NA, as.numeric(eligible$iadl >= 1))

screen_vars <- c("scr_fatigue", "scr_srh", "scr_walk", "scr_adl", "scr_iadl")
context_vars <- c("hypertension", "diabetes", "heart_disease", "cancer", "stroke",
                  "arthritis", "not_married", "low_education")

availability <- bind_rows(lapply(split(eligible, eligible$cohort), function(x) {
  data.frame(cohort = x$cohort[1], variable = c(context_vars, screen_vars, "fi_core13", "fu_fi_core13"),
             N = nrow(x), available_n = sapply(x[, c(context_vars, screen_vars, "fi_core13", "fu_fi_core13"), drop = FALSE], function(z) sum(!is.na(z))),
             available_pct = 100 * sapply(x[, c(context_vars, screen_vars, "fi_core13", "fu_fi_core13"), drop = FALSE], function(z) mean(!is.na(z))),
             stringsAsFactors = FALSE)
}))

sample_flow <- bind_rows(
  data.frame(stage = "Same-wave records", cohort = names(table(combined$cohort)), N = as.integer(table(combined$cohort))),
  data.frame(stage = "Eligible age 65+ with chronic disease", cohort = names(table(eligible$cohort)), N = as.integer(table(eligible$cohort))),
  data.frame(stage = "Complete five-item bank", cohort = names(table(eligible$cohort[complete.cases(eligible[, screen_vars])])), N = as.integer(table(eligible$cohort[complete.cases(eligible[, screen_vars])]))),
  data.frame(stage = "Follow-up FI available", cohort = names(table(eligible$cohort[!is.na(eligible$fu_fi_core13)])), N = as.integer(table(eligible$cohort[!is.na(eligible$fu_fi_core13)])))
)

# Explicitly overwrite the prior intermediate data and QC outputs.
saveRDS(combined, file.path(workspace_dir, "harmonized_data.rds"))
saveRDS(combined, file.path(workspace_dir, "combined_with_bmi.rds"))
saveRDS(eligible, file.path(workspace_dir, "eligible_with_bmi.rds"))
saveRDS(eligible, file.path(workspace_dir, "eligible_with_phenotypes.rds"))
write.csv(availability, file.path(output_dir, "RECOVERED_VARIABLE_AVAILABILITY.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(sample_flow, file.path(output_dir, "RECOVERED_SAMPLE_FLOW.csv"), row.names = FALSE, fileEncoding = "UTF-8")

qc <- list(
  generated_at = Sys.time(),
  baseline_followup_waves = data.frame(cohort=c("HRS","ELSA","SHARE","CHARLS","MHAS"), baseline=c(14,9,2,3,5), followup=c(15,10,4,4,6)),
  N_combined = nrow(combined), N_eligible = nrow(eligible),
  duplicate_uid = sum(duplicated(combined$uid)),
  five_item_complete = sum(complete.cases(eligible[, screen_vars])),
  followup_fi_available = sum(!is.na(eligible$fu_fi_core13))
)
saveRDS(qc, file.path(workspace_dir, "qc_results.rds"))

cat("Combined:", nrow(combined), "\n")
cat("Eligible:", nrow(eligible), "\n")
print(table(eligible$cohort))
cat("Complete five-item bank:", qc$five_item_complete, "\n")
cat("Follow-up FI available:", qc$followup_fi_available, "\n")
cat("=== Recovery complete ===\n")
