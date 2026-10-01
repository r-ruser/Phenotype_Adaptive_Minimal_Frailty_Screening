options(stringsAsFactors=FALSE,warn=1)
try(Sys.setlocale("LC_ALL","Chinese (Simplified)_China.utf8"),silent=TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)[1])
package_dir <- Sys.getenv("MS5_PACKAGE_DIR", unset="")
if (!nzchar(package_dir)) package_dir <- dirname(dirname(normalizePath(script_file, winslash="/", mustWork=TRUE)))
source(file.path(package_dir,"code/00_config.R"), encoding="UTF-8")
suppressPackageStartupMessages({library(haven);library(dplyr)})
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


combined$age_65plus <- !is.na(combined$age) & combined$age >= 65
chronic_vars <- c("hypertension","diabetes","heart_disease","cancer","stroke","arthritis")
combined$has_chronic <- rowSums(combined[,chronic_vars,drop=FALSE] == 1,na.rm=TRUE)>=1
combined$has_social <- !is.na(combined$married) | !is.na(combined$education)
combined$eligible <- combined$age_65plus & combined$has_chronic & combined$has_social
eligible <- combined[combined$eligible %in% TRUE,,drop=FALSE]
eligible$incident_adl <- ifelse(is.na(eligible$adl)|is.na(eligible$fu_adl)|eligible$adl>=1, NA_real_,as.numeric(eligible$fu_adl>=1))
saveRDS(eligible,eligible_input)
flow <- data.frame(cohort=names(table(eligible$cohort)),eligible_n=as.integer(table(eligible$cohort)))
write.csv(flow,file.path(out,"baseline_eligible_counts.csv"),row.names=FALSE)
print(flow)
