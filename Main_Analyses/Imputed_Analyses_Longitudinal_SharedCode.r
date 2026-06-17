# File header ----
# ACEs and Vascular Health — Imputed Mixed-Effects Trajectory Analyses (Classic ACEs)
# Study:   Adverse Childhood Experiences and Changes in Vascular Health from
#          Childhood to Mid-Adulthood: Cross-Sectional and Longitudinal
#          Evidence from the ALSPAC Study
# Cohort:  Avon Longitudinal Study of Parents and Children (ALSPAC)
#          Boyd A, et al. (2012) Cohort Profile: the 'children of the 90s'.
#          International Journal of Epidemiology, 42(1):111-127.
#          https://doi.org/10.1093/ije/dys064
#
# NOTE:    ALSPAC data are available to approved researchers via application.
#          This code is shared for transparency and reproducibility only.
#          The underlying data cannot be shared and are not included.
#          Data access: https://www.bristol.ac.uk/alspac/researchers/access/
#
# ADAPTING FOR EXTENDED ACEs:
#   The Extended ACEs section of this script (Section 3) uses the same pipeline
#   functions defined in Section 1. The only differences are:
#
#   1. Exposure variable names:
#        Classic continuous  → Classic_ACEs
#        Classic categorical → Classic_ACEs_cat
#        Extended continuous → Extended_ACEs
#      (No binary equivalent is used in the Extended analyses.)
#
#   2. Formula object names:
#        Classic  → forms_whole_PWV,     forms_sex_PWV,
#                   forms_whole_cIMT,    forms_sex_cIMT
#        Cat      → forms_whole_PWV_cat, forms_sex_PWV_cat,
#                   forms_whole_cIMT_cat,forms_sex_cIMT_cat
#        Extended → forms_whole_PWV_ext, forms_sex_PWV_ext,
#                   forms_whole_cIMT_ext,forms_sex_cIMT_ext
#
#   3. Input mids object:
#        Classic/Cat → Study1_Imp_Classic_Final_Transformed_V2
#        Extended    → Study1_Imp_Ext_Trns_Reduced
#
#   4. pick_forms helper:
#        Classic  → pick_forms()
#        Cat      → pick_forms_cat()
#        Extended → pick_forms_ext()
#
# Author:  Laura Macro
# Date:    2026


# Section 1: Required packages ----

library(mice)
library(dplyr)
library(lme4)
library(tidyr)


# Section 2: Pipeline functions ----

# --- Helper: create unique child ID from mother ID + sibling order ---
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


# --- Helper: reshape wide imputed data to long format ---
# Creates time-varying Y (outcome), age_clinic_centered, BP_systolic, and BMI
# from the wide-format imputed dataset.
make_classic_long_imp <- function(dat) {
  dat %>%
    pivot_longer(
      cols = intersect(c("PWV_17", "PWV_24", "cIMT_17", "cIMT_24"), names(dat)),
      names_to  = c("Measure", "Age_lbl"),
      names_sep = "_",
      values_to = "Value"
    ) %>%
    mutate(
      age_clinic = case_when(
        Age_lbl == "17" ~ Age_17_clinic_years,
        Age_lbl == "24" ~ Age_24_clinic_years,
        TRUE ~ NA_real_
      ),

      # age17_c / age24_c are centered age variables created during imputation
      age_clinic_centred = case_when(
        Age_lbl == "17" ~ age17_c,
        Age_lbl == "24" ~ age24_c,
        TRUE ~ NA_real_
      ),

      BP_systolic = case_when(
        Age_lbl == "17" ~ BP_systolic_17,
        Age_lbl == "24" ~ BP_systolic_24,
        TRUE ~ NA_real_
      ),

      BMI = case_when(
        Age_lbl == "17" ~ BMI_17,
        Age_lbl == "24" ~ BMI_24,
        TRUE ~ NA_real_
      )
    )
}


# --- Helper: filter long data to outcome and sex subgroup ---
get_outcome_subset_imp <- function(long_dat,
                                    outcome = c("PWV", "cIMT"),
                                    sex = c("All", "Female", "Male"),
                                    restrict_PA_Diet = TRUE) {

  outcome <- match.arg(outcome)
  sex     <- match.arg(sex)

  df <- long_dat %>%
    dplyr::filter(Measure == outcome) %>%
    dplyr::mutate(Y = Value)

  if (sex != "All") df <- df %>% dplyr::filter(Child_sex == sex)

  if (restrict_PA_Diet) {
    df <- df %>%
      dplyr::filter(
        PA_include_15 == TRUE,
        Diet_include_13 == TRUE
      )
  }

  df
}


# --- Core fitting function: runs lmer across all imputations ---
# Returns pooled results, a mira object, or a list of fits depending on
# the 'return' argument.
fit_and_pool <- function(mids_obj,
                         outcome = c("PWV", "cIMT"),
                         sex = c("All", "Female", "Male"),
                         fixed_formula,
                         mother_id = "cidB4619",
                         child_id  = "child_id",
                         random_slope_age = FALSE,
                         return = c("pooled", "mira", "fits"),
                         optimizer = "bobyqa",
                         maxfun = 2e5,
                         verbose_imp = TRUE,
                         restrict_PA_Diet = TRUE) {

  outcome <- match.arg(outcome)
  sex     <- match.arg(sex)
  return  <- match.arg(return)

  m <- mids_obj$m

  fits <- lapply(seq_len(m), function(i) {

    if (verbose_imp) message("  Imputation ", i, "/", m)

    dat_i <- mice::complete(mids_obj, i) %>%
      make_child_id_wide(mother_id = mother_id,
                         order_vars = c("qlet", "Child_sex"))

    long_i <- make_classic_long_imp(dat_i)

    df_i <- get_outcome_subset_imp(
      long_i,
      outcome = outcome,
      sex = sex,
      restrict_PA_Diet = restrict_PA_Diet
    ) %>%
      dplyr::mutate(
        !!child_id := as.factor(.data[[child_id]])
      )

    re_term <- if (random_slope_age) {
      paste0("(1 + age_clinic_centred | ", child_id, ")")
    } else {
      paste0("(1 | ", child_id, ")")
    }

    full_formula <- update(fixed_formula, paste(". ~ . +", re_term))

    lme4::lmer(
      full_formula,
      data = df_i,
      control = lme4::lmerControl(
        optimizer = optimizer,
        optCtrl = list(maxfun = maxfun)
      )
    )
  })

  if (return == "fits") return(fits)

  mira <- mice::as.mira(fits)
  if (return == "mira") return(mira)

  mice::pool(mira)
}


# --- Table builders: variance components, AIC, pooled coefficients ---

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

ensure_mira <- function(x) {
  if (inherits(x, "mira")) return(x)
  if (is.list(x) && length(x) > 0 && inherits(x[[1]], "merMod")) return(mice::as.mira(x))
  stop("Expected a 'mira' object or a list of lmer fits (merMod).")
}

pool_coef_tbl_lmer <- function(mira_obj) {
  s <- summary(pool(mira_obj), conf.int = TRUE) |> as.data.frame()

  CI <- dplyr::case_when(
    all(c("2.5 %", "97.5 %") %in% names(s)) ~
      paste0("[", round(as_num(s$`2.5 %`), 3), ", ", round(as_num(s$`97.5 %`), 3), "]"),
    all(c("conf.low", "conf.high") %in% names(s)) ~
      paste0("[", round(as_num(s$conf.low), 3), ", ", round(as_num(s$conf.high), 3), "]"),
    TRUE ~ NA_character_
  )

  data.frame(
    term      = s$term,
    estimate  = round(as_num(s$estimate), 3),
    std.error = round(as_num(s$std.error), 3),
    statistic = round(as_num(s$statistic), 2),
    df        = as.character(round(as_num(s$df), 1)),
    p.value   = round(as_num(s$p.value), 3),
    CI        = CI,
    stringsAsFactors = FALSE
  )
}

var_components_tbl <- function(mira_obj, groups = "child_id") {

  by_imp <- dplyr::bind_rows(lapply(seq_along(mira_obj$analyses), function(i) {
    mod <- mira_obj$analyses[[i]]
    vc  <- as.data.frame(VarCorr(mod))
    resid_var <- sigma(mod)^2

    ri <- vc %>%
      dplyr::filter(var1 == "(Intercept)", is.na(var2), grp %in% groups) %>%
      dplyr::transmute(term = paste0("Var(random intercept: ", grp, ")"),
                       value = vcov)

    rs <- vc %>%
      dplyr::filter(!is.na(var1), var1 != "(Intercept)", is.na(var2), grp %in% groups) %>%
      dplyr::transmute(term = paste0("Var(random slope: ", grp, " - ", var1, ")"),
                       value = vcov)

    base <- data.frame(term = "Var(residual)", value = resid_var,
                       stringsAsFactors = FALSE)

    total_var <- sum(c(ri$value, rs$value, resid_var), na.rm = TRUE)

    dplyr::bind_rows(
      ri, rs, base,
      data.frame(term = "Total variance (sum components)", value = total_var,
                 stringsAsFactors = FALSE)
    ) %>%
      dplyr::mutate(imp = i)
  }))

  by_imp %>%
    dplyr::group_by(term) %>%
    dplyr::summarise(
      estimate  = round(mean(value, na.rm = TRUE), 3),
      std.error = round(sd(value, na.rm = TRUE), 3),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      statistic = NA_real_,
      df        = NA_character_,
      p.value   = NA_real_,
      CI        = NA_character_
    ) %>%
    dplyr::select(term, estimate, std.error, statistic, df, p.value, CI)
}

pool_AIC_tbl <- function(mira_obj) {
  aics <- vapply(mira_obj$analyses, AIC, numeric(1))
  data.frame(
    term      = "AIC (mean across imputations)",
    estimate  = round(mean(aics, na.rm = TRUE), 2),
    std.error = round(sd(aics, na.rm = TRUE), 2),
    statistic = NA_real_,
    df        = NA_character_,
    p.value   = NA_real_,
    CI        = NA_character_,
    stringsAsFactors = FALSE
  )
}

build_mixed_table <- function(mira_or_fits, groups = "child_id") {
  mira <- ensure_mira(mira_or_fits)

  dplyr::bind_rows(
    var_components_tbl(mira, groups = groups),
    pool_AIC_tbl(mira),
    pool_coef_tbl_lmer(mira)
  )
}


# --- Helper: print results for one outcome/sex/model to console ---
print_mixed_results <- function(mira, label, groups = "child_id") {
  cat("\n====", label, "====\n")
  tab <- build_mixed_table(mira, groups = groups)
  print(tab)
}


# --- Helper: select formula list by outcome and group ---
sex_label <- function(sex) if (sex == "All") "Whole sample" else paste0("Sex-stratified: ", sex)

pick_forms <- function(outcome, group) {
  if (outcome == "PWV") {
    if (group == "Whole") forms_whole_PWV else forms_sex_PWV
  } else {
    if (group == "Whole") forms_whole_cIMT else forms_sex_cIMT
  }
}

pick_forms_cat <- function(outcome, group) {
  if (outcome == "PWV") {
    if (group == "Whole") forms_whole_PWV_cat else forms_sex_PWV_cat
  } else {
    if (group == "Whole") forms_whole_cIMT_cat else forms_sex_cIMT_cat
  }
}

pick_forms_ext <- function(outcome, group) {
  if (outcome == "PWV") {
    if (group == "Whole") forms_whole_PWV_ext else forms_sex_PWV_ext
  } else {
    if (group == "Whole") forms_whole_cIMT_ext else forms_sex_cIMT_ext
  }
}


# Section 3: Classic ACEs (continuous) ----

# --- Fixed effects formulas ---

# PWV (includes BP_systolic)
Model1_PWV_Classic <- Y ~ age_clinic_centred * Classic_ACEs + BP_systolic

Model2__PWV_Classic_whole <- Y ~ age_clinic_centred * Classic_ACEs * Child_sex + BP_systolic +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI

Model2_PWV_Classic_sex <- Y ~ age_clinic_centred * Classic_ACEs + BP_systolic +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI

Model3_PWV_Classic <- Y ~ age_clinic_centred * Classic_ACEs + BP_systolic +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI + Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
  non_milk_sugar_13 + Weekday_sleep_duration_15y + Weekend_sleep_duration_15y +
  Child_alc_15 + Child_smoke_15

Model4_PWV_Classic <- Y ~ age_clinic_centred * Classic_ACEs + BP_systolic +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke +
  Mat_preg_alc + Family_CVD + Age_PHV_c + BMI + Cortisol_15 + CRP_15 + IL6_24

Model5_PWV_Classic <- Y ~ age_clinic_centred * Classic_ACEs + BP_systolic +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI + Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
  non_milk_sugar_13 + Weekday_sleep_duration_15y + Weekend_sleep_duration_15y +
  Child_alc_15 + Child_smoke_15 + Cortisol_15 + CRP_15 + IL6_24

# cIMT (excludes BP_systolic)
Model1_cIMT_Classic <- Y ~ age_clinic_centred * Classic_ACEs

Model2_cIMT_Classic_whole <- Y ~ age_clinic_centred * Classic_ACEs * Child_sex +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI

Model2_cIMT_Classic_sex <- Y ~ age_clinic_centred * Classic_ACEs +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI

Model3_cIMT_Classic <- Y ~ age_clinic_centred * Classic_ACEs +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI + Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
  non_milk_sugar_13 + Weekday_sleep_duration_15y + Weekend_sleep_duration_15y +
  Child_alc_15 + Child_smoke_15

Model4_cIMT_Classic <- Y ~ age_clinic_centred * Classic_ACEs +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu + Mat_PND_gest +
  Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI + Cortisol_15 + CRP_15 + IL6_24

Model5_cIMT_Classic <- Y ~ age_clinic_centred * Classic_ACEs +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu + Mat_PND_gest +
  Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI + Daily_MVPA_15 + Daily_Light_PA_15 +
  Diet_pattern_13_calib + non_milk_sugar_13 + Weekday_sleep_duration_15y +
  Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15 + Cortisol_15 + CRP_15 + IL6_24

# Formula lists
forms_whole_PWV <- list(
  M1 = Model1_PWV_Classic,
  M2 = Model2__PWV_Classic_whole,
  M3 = Model3_PWV_Classic,
  M4 = Model4_PWV_Classic,
  M5 = Model5_PWV_Classic
)

forms_sex_PWV <- forms_whole_PWV
forms_sex_PWV$M2 <- Model2_PWV_Classic_sex

forms_whole_cIMT <- list(
  M1 = Model1_cIMT_Classic,
  M2 = Model2_cIMT_Classic_whole,
  M3 = Model3_cIMT_Classic,
  M4 = Model4_cIMT_Classic,
  M5 = Model5_cIMT_Classic
)

forms_sex_cIMT <- forms_whole_cIMT
forms_sex_cIMT$M2 <- Model2_cIMT_Classic_sex


# --- Run models: fit across all imputations, print results to console ---

for (model_id in c("M1", "M2", "M3", "M4", "M5")) {
  for (outcome in c("PWV", "cIMT")) {
    for (sex in c("All", "Female", "Male")) {

      group <- if (sex == "All") "Whole" else "Sex"
      fml   <- pick_forms(outcome, group)[[model_id]]

      message("Fitting: Classic ", outcome, " ", model_id, " ", sex)

      mira <- fit_and_pool(
        mids_obj      = Study1_Imp_Classic_Final_Transformed_V2,
        outcome       = outcome,
        sex           = sex,
        fixed_formula = fml,
        child_id      = "child_id",
        random_slope_age = FALSE,
        return        = "mira",
        optimizer     = "bobyqa",
        maxfun        = 2e5,
        verbose_imp   = TRUE
      )

      print_mixed_results(
        mira,
        label  = paste("Classic ACEs (continuous) |", outcome, "|", model_id, "|", sex_label(sex)),
        groups = "child_id"
      )

      rm(mira)
      gc()
    }
  }
}


# Section 4: Classic ACEs (categorical) ----

# --- Fixed effects formulas ---

# PWV (includes BP_systolic)
Model1_PWV_Classic_cat <- Y ~ age_clinic_centred * Classic_ACEs_cat + BP_systolic

Model2__PWV_Classic_cat_whole <- Y ~ age_clinic_centred * Classic_ACEs_cat * Child_sex +
  BP_systolic + Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI

Model2_PWV_Classic_cat_sex <- Y ~ age_clinic_centred * Classic_ACEs_cat + BP_systolic +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI

Model3_PWV_Classic_cat <- Y ~ age_clinic_centred * Classic_ACEs_cat + BP_systolic +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI + Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
  non_milk_sugar_13 + Weekday_sleep_duration_15y + Weekend_sleep_duration_15y +
  Child_alc_15 + Child_smoke_15

Model4_PWV_Classic_cat <- Y ~ age_clinic_centred * Classic_ACEs_cat + BP_systolic +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke +
  Mat_preg_alc + Family_CVD + Age_PHV_c + BMI + Cortisol_15 + CRP_15 + IL6_24

Model5_PWV_Classic_cat <- Y ~ age_clinic_centred * Classic_ACEs_cat + BP_systolic +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI + Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
  non_milk_sugar_13 + Weekday_sleep_duration_15y + Weekend_sleep_duration_15y +
  Child_alc_15 + Child_smoke_15 + Cortisol_15 + CRP_15 + IL6_24

# cIMT (excludes BP_systolic)
Model1_cIMT_Classic_cat <- Y ~ age_clinic_centred * Classic_ACEs_cat

Model2_cIMT_Classic_cat_whole <- Y ~ age_clinic_centred * Classic_ACEs_cat * Child_sex +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI

Model2_cIMT_Classic_cat_sex <- Y ~ age_clinic_centred * Classic_ACEs_cat +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI

Model3_cIMT_Classic_cat <- Y ~ age_clinic_centred * Classic_ACEs_cat +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI + Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
  non_milk_sugar_13 + Weekday_sleep_duration_15y + Weekend_sleep_duration_15y +
  Child_alc_15 + Child_smoke_15

Model4_cIMT_Classic_cat <- Y ~ age_clinic_centred * Classic_ACEs_cat +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu + Mat_PND_gest +
  Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI + Cortisol_15 + CRP_15 + IL6_24

Model5_cIMT_Classic_cat <- Y ~ age_clinic_centred * Classic_ACEs_cat +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu + Mat_PND_gest +
  Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI + Daily_MVPA_15 + Daily_Light_PA_15 +
  Diet_pattern_13_calib + non_milk_sugar_13 + Weekday_sleep_duration_15y +
  Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15 + Cortisol_15 + CRP_15 + IL6_24

# Formula lists
forms_whole_PWV_cat <- list(
  M1 = Model1_PWV_Classic_cat,
  M2 = Model2__PWV_Classic_cat_whole,
  M3 = Model3_PWV_Classic_cat,
  M4 = Model4_PWV_Classic_cat,
  M5 = Model5_PWV_Classic_cat
)

forms_sex_PWV_cat <- forms_whole_PWV_cat
forms_sex_PWV_cat$M2 <- Model2_PWV_Classic_cat_sex

forms_whole_cIMT_cat <- list(
  M1 = Model1_cIMT_Classic_cat,
  M2 = Model2_cIMT_Classic_cat_whole,
  M3 = Model3_cIMT_Classic_cat,
  M4 = Model4_cIMT_Classic_cat,
  M5 = Model5_cIMT_Classic_cat
)

forms_sex_cIMT_cat <- forms_whole_cIMT_cat
forms_sex_cIMT_cat$M2 <- Model2_cIMT_Classic_cat_sex


# --- Run models: fit across all imputations, print results to console ---

for (model_id in c("M1", "M2", "M3", "M4", "M5")) {
  for (outcome in c("PWV", "cIMT")) {
    for (sex in c("All", "Female", "Male")) {

      group <- if (sex == "All") "Whole" else "Sex"
      fml   <- pick_forms_cat(outcome, group)[[model_id]]

      message("Fitting: Classic categorical ", outcome, " ", model_id, " ", sex)

      mira <- fit_and_pool(
        mids_obj      = Study1_Imp_Classic_Final_Transformed_V2,
        outcome       = outcome,
        sex           = sex,
        fixed_formula = fml,
        child_id      = "cidB4619",
        random_slope_age = FALSE,
        return        = "mira",
        optimizer     = "bobyqa",
        maxfun        = 2e5,
        verbose_imp   = TRUE
      )

      print_mixed_results(
        mira,
        label  = paste("Classic ACEs (categorical) |", outcome, "|", model_id, "|", sex_label(sex)),
        groups = "cidB4619"
      )

      rm(mira)
      gc()
    }
  }
}


# Section 5: Extended ACEs (continuous) ----

# --- Fixed effects formulas ---

# PWV (includes BP_systolic)
Model1_PWV_Extended <- Y ~ age_clinic_centred * Extended_ACEs + BP_systolic

Model2__PWV_Extended_whole <- Y ~ age_clinic_centred * Extended_ACEs * Child_sex +
  BP_systolic + Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI

Model2_PWV_Extended_sex <- Y ~ age_clinic_centred * Extended_ACEs + BP_systolic +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI

Model3_PWV_Extended <- Y ~ age_clinic_centred * Extended_ACEs + BP_systolic +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI + Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
  non_milk_sugar_13 + Weekday_sleep_duration_15y + Weekend_sleep_duration_15y +
  Child_alc_15 + Child_smoke_15

Model4_PWV_Extended <- Y ~ age_clinic_centred * Extended_ACEs + BP_systolic +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke +
  Mat_preg_alc + Family_CVD + Age_PHV_c + BMI + Cortisol_15 + CRP_15 + IL6_24

Model5_PWV_Extended <- Y ~ age_clinic_centred * Extended_ACEs + BP_systolic +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI + Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
  non_milk_sugar_13 + Weekday_sleep_duration_15y + Weekend_sleep_duration_15y +
  Child_alc_15 + Child_smoke_15 + Cortisol_15 + CRP_15 + IL6_24

# cIMT (excludes BP_systolic)
Model1_cIMT_Extended <- Y ~ age_clinic_centred * Extended_ACEs

Model2_cIMT_Extended_whole <- Y ~ age_clinic_centred * Extended_ACEs * Child_sex +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI

Model2_cIMT_Extended_sex <- Y ~ age_clinic_centred * Extended_ACEs +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI

Model3_cIMT_Extended <- Y ~ age_clinic_centred * Extended_ACEs +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
  Mat_PND_gest + Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI + Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
  non_milk_sugar_13 + Weekday_sleep_duration_15y + Weekend_sleep_duration_15y +
  Child_alc_15 + Child_smoke_15

Model4_cIMT_Extended <- Y ~ age_clinic_centred * Extended_ACEs +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu + Mat_PND_gest +
  Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI + Cortisol_15 + CRP_15 + IL6_24

Model5_cIMT_Extended <- Y ~ age_clinic_centred * Extended_ACEs +
  Child_ethnicity + Townsend_sum + Marital_status + Parent_edu + Mat_PND_gest +
  Mat_age_delivery + Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
  Family_CVD + Age_PHV_c + BMI + Daily_MVPA_15 + Daily_Light_PA_15 +
  Diet_pattern_13_calib + non_milk_sugar_13 + Weekday_sleep_duration_15y +
  Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15 + Cortisol_15 + CRP_15 + IL6_24

# Formula lists
forms_whole_PWV_ext <- list(
  M1 = Model1_PWV_Extended,
  M2 = Model2__PWV_Extended_whole,
  M3 = Model3_PWV_Extended,
  M4 = Model4_PWV_Extended,
  M5 = Model5_PWV_Extended
)

forms_sex_PWV_ext <- forms_whole_PWV_ext
forms_sex_PWV_ext$M2 <- Model2_PWV_Extended_sex

forms_whole_cIMT_ext <- list(
  M1 = Model1_cIMT_Extended,
  M2 = Model2_cIMT_Extended_whole,
  M3 = Model3_cIMT_Extended,
  M4 = Model4_cIMT_Extended,
  M5 = Model5_cIMT_Extended
)

forms_sex_cIMT_ext <- forms_whole_cIMT_ext
forms_sex_cIMT_ext$M2 <- Model2_cIMT_Extended_sex


# --- Run models: fit across all imputations, print results to console ---

for (model_id in c("M1", "M2", "M3", "M4", "M5")) {
  for (outcome in c("PWV", "cIMT")) {
    for (sex in c("All", "Female", "Male")) {

      group <- if (sex == "All") "Whole" else "Sex"
      fml   <- pick_forms_ext(outcome, group)[[model_id]]

      message("Fitting: Extended ", outcome, " ", model_id, " ", sex)

      mira <- fit_and_pool(
        mids_obj      = Study1_Imp_Ext_Trns_Reduced,
        outcome       = outcome,
        sex           = sex,
        fixed_formula = fml,
        child_id      = "cidB4619",
        random_slope_age = FALSE,
        return        = "mira",
        optimizer     = "bobyqa",
        maxfun        = 2e5,
        verbose_imp   = TRUE
      )

      print_mixed_results(
        mira,
        label  = paste("Extended ACEs (continuous) |", outcome, "|", model_id, "|", sex_label(sex)),
        groups = "cidB4619"
      )

      rm(mira)
      gc()
    }
  }
}


# Section 6: Clean up ----

rm(forms_whole_PWV,     forms_whole_cIMT,     forms_sex_PWV,     forms_sex_cIMT,
   forms_whole_PWV_cat, forms_whole_cIMT_cat, forms_sex_PWV_cat, forms_sex_cIMT_cat,
   forms_whole_PWV_ext, forms_whole_cIMT_ext, forms_sex_PWV_ext, forms_sex_cIMT_ext)
