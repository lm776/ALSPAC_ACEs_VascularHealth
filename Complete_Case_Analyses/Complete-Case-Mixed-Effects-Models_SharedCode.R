# Complete-Case Longitudinal Mixed-Effects Models
# Study: Adverse Childhood Experiences and Changes in Vascular Health from 
#    Childhood to Mid-Adulthood: Cross-Sectional and Longitudinal Evidence from
#    the ALSPAC Study
#
# Description:
#   Fits linear mixed-effects models (random intercept per child) for
#   pulse wave velocity (PWV) and carotid intima-media thickness (cIMT)
#   outcomes using the Classic ACEs exposure score. Models M1-M5 cover
#   unadjusted through fully adjusted specifications (whole sample only).
#   One Word document per model specification is written to disk.
#
# Data:
#   ALSPAC (Avon Longitudinal Study of Parents and Children).
#   Data access: https://www.bristol.ac.uk/alspac/researchers/access/
#   Please cite: Boyd A et al. (2013) Cohort Profile: the 'Children of the
#   90s'—the index offspring of the Avon Longitudinal Study of Parents and
#   Children. Int J Epidemiol, 42(1):111-127.
#   https://doi.org/10.1093/ije/dys064
#
# Dependency:
#   This script expects the object `cov_complete_both_classic` to exist in
#   the R environment. That object is created in
#   Complete-Case-Regression-Analyses-SharedCode.R, which must be run first.
#   `cov_complete_both_classic` is the subset of participants with complete
#   data on Classic ACEs, both vascular outcomes (17y and 24y), and all
#   covariates.
#
# Variable notes:
#   - cidB4619: ALSPAC pregnancy/mother identifier (publicly documented in
#               the ALSPAC data dictionary). Used here to define the random
#               intercept grouping structure (one row per child per time point
#               after pivoting to long format).
#   - age17_c / age24_c: clinic age centred at the grand mean (created in
#                         the regression script).
#   - Age_PHV_c: age at peak height velocity, centred at the sample mean
#                (created in the regression script).
#   - BP_systolic in model formulas: this is the *time-varying* systolic BP
#                  variable created during pivoting to long format (takes
#                  BP_systolic_17 at age 17 and BP_systolic_24 at age 24).
#
# Author:  Laura Macro
# Date:    June 2026


# Required packages ----

library(dplyr)
library(tidyr)
library(tibble)
library(officer)
library(flextable)
library(lme4)
library(lmerTest)   # extends lme4 to provide Satterthwaite p-values
library(broom.mixed)


# Section 1: Construct child IDs ----
# cidB4619 is the ALSPAC pregnancy ID (one per mother). Sibling pairs within
# a family share the same cidB4619; we append a letter suffix (A, B, ...)
# to create a unique child-level ID for the random intercept.

make_child_id_wide <- function(dat,
                               mother_id = "cidB4619",
                               order_vars = c("qlet", "Child_sex")) {
  dat %>%
    dplyr::arrange(dplyr::across(dplyr::all_of(c(mother_id, order_vars)))) %>%
    dplyr::group_by(.data[[mother_id]]) %>%
    dplyr::mutate(
      sib_index  = dplyr::row_number(),
      sib_suffix = LETTERS[sib_index],
      child_id   = paste0(.data[[mother_id]], sib_suffix)
    ) %>%
    dplyr::ungroup()
}

cov_complete_both_classic <- make_child_id_wide(
  cov_complete_both_classic,
  mother_id  = "cidB4619",
  order_vars = c("qlet", "Child_sex")
)


# Section 2: Pivot to long format ----
# Creates time-varying columns: age_clinic, age_clinic_centred, BP_systolic,
# BMI — each taking the 17y or 24y value depending on the time point.

make_classic_long_cc <- function(dat) {
  dat %>%
    tidyr::pivot_longer(
      cols      = intersect(c("PWV_17", "PWV_24", "cIMT_17", "cIMT_24"),
                            names(dat)),
      names_to  = c("Measure", "Age_lbl"),
      names_sep = "_",
      values_to = "Value"
    ) %>%
    dplyr::mutate(
      age_clinic = dplyr::case_when(
        Age_lbl == "17" ~ Age_17_clinic_years,
        Age_lbl == "24" ~ Age_24_clinic_years,
        TRUE ~ NA_real_
      ),
      age_clinic_centred = dplyr::case_when(
        Age_lbl == "17" ~ age17_c,
        Age_lbl == "24" ~ age24_c,
        TRUE ~ NA_real_
      ),
      BP_systolic = dplyr::case_when(
        Age_lbl == "17" ~ BP_systolic_17,
        Age_lbl == "24" ~ BP_systolic_24,
        TRUE ~ NA_real_
      ),
      BMI = dplyr::case_when(
        Age_lbl == "17" ~ BMI_17,
        Age_lbl == "24" ~ BMI_24,
        TRUE ~ NA_real_
      )
    )
}


# Section 3: Subset long data by outcome ----

get_outcome_subset_cc <- function(long_dat,
                                  outcome = c("PWV", "cIMT")) {
  outcome <- match.arg(outcome)
  long_dat %>%
    dplyr::filter(Measure == outcome) %>%
    dplyr::mutate(Y = Value)
}


# Section 4: Select formula list by outcome ----

pick_forms <- function(outcome) {
  if (outcome == "PWV") forms_whole_PWV else forms_whole_cIMT
}


# Section 5: Fit a single mixed-effects model ----
# Uses lmerTest::lmer for Satterthwaite denominator df and p-values.
# Random intercept is added automatically from the fixed formula.

fit_cc_one <- function(fixed_formula, df, id = "cidB4619") {

  if (!id %in% names(df)) stop("ID column missing: ", id)
  df[[id]] <- as.factor(df[[id]])

  full_fml <- update(
    fixed_formula,
    as.formula(paste0(". ~ . + (1 | ", id, ")"))
  )

  tryCatch(
    lmerTest::lmer(
      full_fml,
      data      = df,
      na.action = na.omit,
      control   = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl   = list(maxfun = 2e5)
      )
    ),
    error = function(e) e
  )
}


# Section 6: Build a results table for one model ----

cc_results_table <- function(mod, id = "cidB4619") {

  if (!inherits(mod, "merMod")) {
    return(tibble::tibble(
      term      = "MODEL FAILED",
      estimate  = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      df        = NA_character_,
      p.value   = NA_real_,
      CI        = conditionMessage(mod)
    ))
  }

  vc        <- as.data.frame(lme4::VarCorr(mod))
  resid_var <- stats::sigma(mod)^2

  ri_var <- vc %>%
    dplyr::filter(grp == id, var1 == "(Intercept)", is.na(var2)) %>%
    dplyr::pull(vcov)

  ri_var    <- if (length(ri_var) == 0) NA_real_ else ri_var[1]
  total_var <- sum(c(ri_var, resid_var), na.rm = TRUE)

  fit_part <- tibble::tibble(
    term      = c("Total variance (sum components)",
                  "Var(random intercept)",
                  "Var(residual)",
                  "AIC"),
    estimate  = c(total_var, ri_var, resid_var, AIC(mod)),
    std.error = NA_real_,
    statistic = NA_real_,
    df        = NA_character_,
    p.value   = NA_real_,
    CI        = NA_character_
  )

  fe_raw <- broom.mixed::tidy(mod, effects = "fixed", conf.int = TRUE)

  df_out   <- if ("df" %in% names(fe_raw))
    as.character(round(as.numeric(fe_raw$df), 1)) else NA_character_

  stat_out <- if ("statistic" %in% names(fe_raw))
    round(as.numeric(fe_raw$statistic), 2) else NA_real_

  p_out    <- if ("p.value" %in% names(fe_raw))
    round(as.numeric(fe_raw$p.value), 3) else NA_real_

  fe <- fe_raw %>%
    dplyr::transmute(
      term      = term,
      estimate  = round(as.numeric(estimate), 3),
      std.error = round(as.numeric(std.error), 3),
      statistic = stat_out,
      df        = df_out,
      p.value   = p_out,
      CI        = paste0("[", round(conf.low, 3), ", ",
                         round(conf.high, 3), "]")
    )

  dplyr::bind_rows(fit_part, fe)
}


# Section 7: Write one Word document per model specification ----
# Whole-sample models only. Outputs to `word_results_CC/` by default.

write_one_doc_per_model_cc <- function(dat_cc,
                                       model_id = c("M1", "M2", "M3",
                                                    "M4", "M5"),
                                       doc_dir  = "word_results_CC",
                                       id       = "cidB4619") {

  model_id <- match.arg(model_id)
  dir.create(doc_dir, showWarnings = FALSE, recursive = TRUE)

  long <- make_classic_long_cc(dat_cc)

  doc <- officer::read_docx() %>%
    officer::body_add_par(
      paste0("Trajectory mixed-effects models (Complete-case, whole sample): ",
             model_id),
      style = "heading 1"
    )

  for (outcome in c("PWV", "cIMT")) {

    doc   <- doc %>% officer::body_add_par(outcome, style = "heading 1")
    forms <- pick_forms(outcome)
    fixed <- forms[[model_id]]
    df    <- get_outcome_subset_cc(long, outcome = outcome)

    message("Fitting (CC whole): ", outcome, " ", model_id)

    mod <- fit_cc_one(fixed, df, id = id)
    tab <- cc_results_table(mod, id = id)

    fit_rows <- which(is.na(tab$statistic))
    ft       <- flextable::flextable(tab) %>% flextable::autofit()
    if (length(fit_rows) > 0)
      ft <- ft %>%
        flextable::bold(i   = fit_rows, bold   = TRUE) %>%
        flextable::italic(i = fit_rows, italic = TRUE)

    doc <- doc %>%
      officer::body_add_par(model_id, style = "heading 2")
    doc <- flextable::body_add_flextable(doc, value = ft)
    doc <- doc %>%
      officer::body_add_par("", style = "Normal")

    rm(mod, tab, ft); gc()
  }

  out_doc <- file.path(doc_dir,
                       paste0("CC_mixed_effects_whole_", model_id, ".docx"))
  print(doc, target = out_doc)
  rm(doc); gc()
  out_doc
}


# Section 8: Define model formulas ----

# --- PWV (includes BP_systolic as time-varying covariate) ---

# M1: age-by-ACEs trajectory, adjusted for BP
Model1_PWV_Classic <- Y ~ age_clinic_centred * Classic_ACEs + BP_systolic

# M2 whole sample: adds sex interaction + full covariate set
Model2_PWV_Classic_whole <- Y ~ age_clinic_centred * Classic_ACEs * Child_sex +
  BP_systolic + Townsend_sum + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
  Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI

# M2 sex-stratified: no sex interaction term (used when running within sex subsets)
Model2_PWV_Classic_sex <- Y ~ age_clinic_centred * Classic_ACEs +
  BP_systolic + Townsend_sum + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
  Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI

# M3: adds lifestyle mediators (PA, diet, sleep, alcohol)
Model3_PWV_Classic <- Y ~ age_clinic_centred * Classic_ACEs +
  BP_systolic + Townsend_sum + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
  Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI +
  Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
  Non_milk_sugar_13 + Weekday_sleep_duration_15y +
  Weekend_sleep_duration_15y + Child_alc_15

# M4: adds biological mediators (cortisol, CRP, IL-6)
Model4_PWV_Classic <- Y ~ age_clinic_centred * Classic_ACEs +
  BP_systolic + Townsend_sum + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
  Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI +
  Cortisol_15 + CRP_15 + IL6_24

# M5: fully adjusted (lifestyle + biological mediators)
Model5_PWV_Classic <- Y ~ age_clinic_centred * Classic_ACEs +
  BP_systolic + Townsend_sum + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
  Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI +
  Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
  Non_milk_sugar_13 + Weekday_sleep_duration_15y +
  Weekend_sleep_duration_15y + Child_alc_15 +
  Cortisol_15 + CRP_15 + IL6_24


# --- cIMT (BP systolic not included) ---

# M1: age-by-ACEs trajectory only
Model1_cIMT_Classic <- Y ~ age_clinic_centred * Classic_ACEs

# M2 whole sample: adds sex interaction + full covariate set
Model2_cIMT_Classic_whole <- Y ~ age_clinic_centred * Classic_ACEs * Child_sex +
  Townsend_sum + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
  Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI

# M2 sex-stratified
Model2_cIMT_Classic_sex <- Y ~ age_clinic_centred * Classic_ACEs +
  Townsend_sum + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
  Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI

# M3: adds lifestyle mediators
Model3_cIMT_Classic <- Y ~ age_clinic_centred * Classic_ACEs +
  Townsend_sum + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
  Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI +
  Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
  Non_milk_sugar_13 + Weekday_sleep_duration_15y +
  Weekend_sleep_duration_15y + Child_alc_15

# M4: adds biological mediators
Model4_cIMT_Classic <- Y ~ age_clinic_centred * Classic_ACEs +
  Townsend_sum + Parent_edu + Mat_PND_gest +
  Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke +
  Mat_preg_alc + Family_CVD + Age_PHV_c + BMI +
  Cortisol_15 + CRP_15 + IL6_24

# M5: fully adjusted
Model5_cIMT_Classic <- Y ~ age_clinic_centred * Classic_ACEs +
  Townsend_sum + Parent_edu + Mat_PND_gest +
  Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke +
  Mat_preg_alc + Family_CVD + Age_PHV_c + BMI +
  Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
  Non_milk_sugar_13 + Weekday_sleep_duration_15y +
  Weekend_sleep_duration_15y +
  Cortisol_15 + CRP_15 + IL6_24


# Section 9: Assemble formula lists ----

forms_whole_PWV <- list(
  M1 = Model1_PWV_Classic,
  M2 = Model2_PWV_Classic_whole,
  M3 = Model3_PWV_Classic,
  M4 = Model4_PWV_Classic,
  M5 = Model5_PWV_Classic
)

# Sex-stratified list: M2 uses the sex-stratified formula (no sex interaction)
# Note: forms_sex_PWV and forms_sex_cIMT are defined here for completeness
# and can be used for sex-stratified runs if required, but are not called in
# the main run loop below (which is whole-sample only).
forms_sex_PWV      <- forms_whole_PWV
forms_sex_PWV$M2   <- Model2_PWV_Classic_sex

forms_whole_cIMT <- list(
  M1 = Model1_cIMT_Classic,
  M2 = Model2_cIMT_Classic_whole,
  M3 = Model3_cIMT_Classic,
  M4 = Model4_cIMT_Classic,
  M5 = Model5_cIMT_Classic
)

forms_sex_cIMT      <- forms_whole_cIMT
forms_sex_cIMT$M2   <- Model2_cIMT_Classic_sex


# Section 10: Run models and write output ----
# Produces 5 Word documents (one per model specification) in word_results_CC/

files_written_cc <- lapply(c("M1", "M2", "M3", "M4", "M5"), function(m) {
  write_one_doc_per_model_cc(
    dat_cc   = cov_complete_both_classic,
    model_id = m,
    doc_dir  = "word_results_CC",
    id       = "cidB4619"   # ALSPAC pregnancy ID; publicly documented variable
  )
})

files_written_cc
