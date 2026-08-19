# File header ----
# ACEs and Vascular Health — Imputed Regression Analyses (Classic ACEs)
# Study:   Adverse Childhood Experiences and Changes in Vascular Health from
#          Adolescence to Early Adulthood: Cross-Sectional and Longitudinal
#          Evidence from the ALSPAC Study
# Cohort:  Avon Longitudinal Study of Parents and Children (ALSPAC)
#          Boyd A, et al. (2012) Cohort Profile: the 'children of the 90s'.
#          International Journal of Epidemiology, 42(1):111-127.
#          https://doi.org/10.1093/eje/dys064
#
# NOTE:    ALSPAC data are available to approved researchers via application.
#          This code is shared for transparency and reproducibility only.
#          The underlying data cannot be shared and are not included.
#          Data access: https://www.bristol.ac.uk/alspac/researchers/access/
#
# ADAPTING FOR EXTENDED ACEs:
#   This script runs all analyses using the Classic ACE score as the exposure.
#   This script can be amended to run the same analyses using the Extended ACE 
#   score. The two scripts are structurally identical except for the following:
#
#   1. Exposure variable names:
#        Classic  → Classic_ACEs  (continuous),  Classic_ACEs_cat  (categorical)
#        Extended → Extended_ACEs (continuous),  Ext_ACEs_cat      (categorical)
#      Note: Classic_ACEs_binary is retained in the variable subset (line ~85)
#      for completeness; no binary equivalent is used in the Extended script.
#
#   2. Input mids object names:
#        Classic  → Study1_Imp_Classic_Final_Transformed_V2
#                   Study1_Imp_Classic_Trns_Reduced
#        Extended → Study1_Imp_Ext_Final_Transformed_V2
#                   Study1_Imp_Ext_Trns_Reduced
#
#   3. Long-format intermediate object names:
#        Classic  → Study1_Classic_Long, Study1_Classic_Long_V2,
#                   Study1_Classic_Long_Subset
#        Extended → Study1_Ext_Long, Study1_Ext_Long_V2,
#                   Study1_Ext_Long_Subset
#
#   4. Model object and list names:
#        All "_Classic_" / "_classic" suffixes → "_Ext_" / "_ext"
#        e.g. PWV_Model1_Classic_17_Whole → PWV_Model1_Ext_17_Whole
#             Model1_fits_classic         → Model1_fits_ext
#
#   5. Saved .rds filenames:
#        e.g. "Model1_fits_classic.rds" → "Model1_fits_ext.rds"
#
#   6. Model assumption check section (Section 5):
#        The mids object used to extract imputation 1 differs:
#        Classic  → Study1_Imp_Classic_Trns_Reduced
#        Extended → Study1_Imp_Ext_Trns_Reduced
#
# Author:  Laura Macro
# Date:    2026


# Section 1: Required packages ----

library(mice)
library(miceadds)
library(car)
library(performance)
library(dplyr)
library(tidyr)
library(ggplot2)


# Section 2: Data preparation ----

# --- Step 1: Convert mids to long ---
Study1_Classic_Long <- complete(
  Study1_Imp_Classic_Final_Transformed_V2,
  action  = "long",
  include = TRUE
)

# --- Step 2: Force ordinal variables to unordered factors ---
# (prevents polynomial contrasts in regression models)
vars_to_unorder <- c(
  "Parent_edu",
  "Mat_preg_smoke",
  "Mat_preg_alc",
  "Child_ethnicity",
  "Child_alc_15"
)

vars_to_unorder <- intersect(vars_to_unorder, names(Study1_Classic_Long))

Study1_Classic_Long[vars_to_unorder] <- lapply(
  Study1_Classic_Long[vars_to_unorder],
  function(x) factor(x, ordered = FALSE)
)

# --- Step 3: Set reference categories ---
if ("Parent_edu" %in% names(Study1_Classic_Long)) {
  Study1_Classic_Long$Parent_edu <- relevel(
    Study1_Classic_Long$Parent_edu,
    ref = "No qualifications"
  )
}

if ("Child_ethnicity" %in% names(Study1_Classic_Long)) {
  Study1_Classic_Long$Child_ethnicity <- relevel(
    Study1_Classic_Long$Child_ethnicity,
    ref = "White"
  )
}

if ("Child_alc_15" %in% names(Study1_Classic_Long)) {
  Study1_Classic_Long$Child_alc_15 <- factor(
    Study1_Classic_Long$Child_alc_15,
    levels = c("None", "1 or 2", "3 to 5", "6 to 9",
               "10 to 19", "20 to 39", "40+"),
    ordered = FALSE
  )
}

# --- Step 4: Optional level checks ---
levels(Study1_Classic_Long$Parent_edu)
levels(Study1_Classic_Long$Child_ethnicity)
levels(Study1_Classic_Long$Child_alc_15)

# --- Step 5: Convert long back to mids ---
Study1_Imp_Classic_Final_Transformed_V2 <- as.mids(Study1_Classic_Long)

saveRDS(Study1_Imp_Classic_Final_Transformed_V2,
        "Study1_Imp_Classic_Final_Transformed_V2.rds")

# --- Step 6: Subset to analysis variables only (speeds up model fitting) ---
# cidB4619 = ALSPAC pregnancy identifier (publicly documented)
all_vars_analyses_Classic <- c(
  "cidB4619",
  "Classic_ACEs", "Classic_ACEs_cat", "Classic_ACEs_binary",
  "Child_sex", "BP_systolic_17", "BP_systolic_24",
  "Age_17_clinic_years", "Age_24_clinic_years",
  "Child_ethnicity", "Townsend_sum", "Marital_status", "Parent_edu",
  "Mat_PND_gest", "Mat_age_delivery", "Birth_weight_kg",
  "Mat_preg_smoke", "Mat_preg_alc", "Family_CVD", "Age_PHV_c",
  "BMI_17", "BMI_24", "Cortisol_15", "CRP_15", "IL6_24",
  "Daily_MVPA_15", "Daily_Light_PA_15", "Diet_pattern_13_calib",
  "non_milk_sugar_13", "Weekday_sleep_duration_15y",
  "Weekend_sleep_duration_15y", "Child_alc_15", "Child_smoke_15",
  "PWV_17", "PWV_24", "cIMT_17", "cIMT_24",
  "Arterial_Dist_17", "PA_include_15", "Diet_include_13",
  ".imp", ".id"
)

Study1_Classic_Long_V2 <- complete(
  Study1_Imp_Classic_Final_Transformed_V2,
  action  = "long",
  include = TRUE
)

Study1_Classic_Long_Subset <- Study1_Classic_Long_V2 %>%
  select(all_of(all_vars_analyses_Classic))

Study1_Imp_Classic_Trns_Reduced <- as.mids(Study1_Classic_Long_Subset)

saveRDS(Study1_Imp_Classic_Trns_Reduced, "Study1_Imp_Classic_Trns_Reduced.rds")

rm(Study1_Classic_Long, Study1_Classic_Long_V2, Study1_Classic_Long_Subset,
   Study1_Imp_Classic_Final_Transformed_V2)


# Section 3: Continuous ACEs — regression models ----

# Model structure (all models):
#   Whole group: Classic_ACEs (or * Child_sex for interaction models)
#   Split by sex: Classic_ACEs only (no interaction term needed)
#   Subset: Models 3 and 5 (which include physical activity and diet as
#           lifestyle covariates) restrict to PA_include_15 == TRUE &
#           Diet_include_13 == TRUE via the `subset =` argument. Models 1,
#           2, and 4 use the full available (unfiltered) sample.
#
# Models 1–5 are progressively adjusted:
#   Model 1: Minimal adjustment (age, BP)
#   Model 2: + sociodemographic and perinatal covariates (with ACE × sex interaction)
#   Model 3: Model 2 + health risk behaviour covariates
#   Model 4: Model 2 + stress/inflammatory biomarkers
#   Model 5: Model 2 + all lifestyle and biomarker covariates (full adjustment)


# Section 3.1: Model 1 (minimal adjustment) ----

# --- PWV ---

PWV_Model1_Classic_17_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years
  )
)

PWV_Model1_Classic_24_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_24 ~ Classic_ACEs + BP_systolic_24 + Age_24_clinic_years
  )
)

PWV_Model1_Classic_17_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years,
    subset = (Child_sex == "Female")
  )
)

PWV_Model1_Classic_17_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years,
    subset = (Child_sex == "Male")
  )
)

PWV_Model1_Classic_24_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_24 ~ Classic_ACEs + BP_systolic_24 + Age_24_clinic_years,
    subset = (Child_sex == "Female")
  )
)

PWV_Model1_Classic_24_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_24 ~ Classic_ACEs + BP_systolic_24 + Age_24_clinic_years,
    subset = (Child_sex == "Male")
  )
)

# --- cIMT ---

cIMT_Model1_Classic_17_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_17 ~ Classic_ACEs + Age_17_clinic_years
  )
)

cIMT_Model1_Classic_24_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_24 ~ Classic_ACEs + Age_24_clinic_years
  )
)

cIMT_Model1_Classic_17_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_17 ~ Classic_ACEs + Age_17_clinic_years,
    subset = (Child_sex == "Female")
  )
)

cIMT_Model1_Classic_17_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_17 ~ Classic_ACEs + Age_17_clinic_years,
    subset = (Child_sex == "Male")
  )
)

cIMT_Model1_Classic_24_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_24 ~ Classic_ACEs + Age_24_clinic_years,
    subset = (Child_sex == "Female")
  )
)

cIMT_Model1_Classic_24_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_24 ~ Classic_ACEs + Age_24_clinic_years,
    subset = (Child_sex == "Male")
  )
)

# --- Arterial distensibility ---

Dist_Model1_Classic_17_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    Arterial_Dist_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years
  )
)

Dist_Model1_Classic_17_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    Arterial_Dist_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years,
    subset = (Child_sex == "Female")
  )
)

Dist_Model1_Classic_17_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    Arterial_Dist_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years,
    subset = (Child_sex == "Male")
  )
)

Model1_fits_classic <- list(
  PWV_17_Whole  = PWV_Model1_Classic_17_Whole,
  PWV_17_F      = PWV_Model1_Classic_17_F,
  PWV_17_M      = PWV_Model1_Classic_17_M,
  PWV_24_Whole  = PWV_Model1_Classic_24_Whole,
  PWV_24_F      = PWV_Model1_Classic_24_F,
  PWV_24_M      = PWV_Model1_Classic_24_M,
  cIMT_17_Whole = cIMT_Model1_Classic_17_Whole,
  cIMT_17_F     = cIMT_Model1_Classic_17_F,
  cIMT_17_M     = cIMT_Model1_Classic_17_M,
  cIMT_24_Whole = cIMT_Model1_Classic_24_Whole,
  cIMT_24_F     = cIMT_Model1_Classic_24_F,
  cIMT_24_M     = cIMT_Model1_Classic_24_M,
  Dist_17_Whole = Dist_Model1_Classic_17_Whole,
  Dist_17_F     = Dist_Model1_Classic_17_F,
  Dist_17_M     = Dist_Model1_Classic_17_M
)

saveRDS(Model1_fits_classic, "Model1_fits_classic.rds", compress = "gzip")

rm(PWV_Model1_Classic_17_Whole,  PWV_Model1_Classic_17_F,  PWV_Model1_Classic_17_M,
   PWV_Model1_Classic_24_Whole,  PWV_Model1_Classic_24_F,  PWV_Model1_Classic_24_M,
   cIMT_Model1_Classic_17_Whole, cIMT_Model1_Classic_17_F, cIMT_Model1_Classic_17_M,
   cIMT_Model1_Classic_24_Whole, cIMT_Model1_Classic_24_F, cIMT_Model1_Classic_24_M,
   Dist_Model1_Classic_17_Whole, Dist_Model1_Classic_17_F, Dist_Model1_Classic_17_M)


# Section 3.2: Pooled results helper functions ----
# These functions are defined once here and reused for all models.

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

pool_overall_F_tbl <- function(mira_obj) {
  fs   <- lapply(mira_obj$analyses, function(mod) summary(mod)$fstatistic)
  Fm   <- vapply(fs, function(x) as_num(x["value"]), numeric(1))
  df1m <- vapply(fs, function(x) as_num(x["numdf"]),  numeric(1))
  df1  <- df1m[1]
  m    <- length(Fm)
  Wi   <- df1 * Fm
  Wbar <- mean(Wi, na.rm = TRUE)
  B    <- stats::var(Wi, na.rm = TRUE)
  r    <- (1 + 1/m) * (B / Wbar)
  F_pooled <- Wbar / (df1 * (1 + r))
  df2_pool <- (m - 1) * (1 + 1/r)^2
  p        <- stats::pf(F_pooled, df1 = df1, df2 = df2_pool, lower.tail = FALSE)
  data.frame(
    term      = "Overall model (all slopes = 0)",
    estimate  = NA_real_,
    std.error = NA_real_,
    statistic = round(F_pooled, 2),
    df        = paste0(round(df1, 0), ", ", round(df2_pool, 1)),
    p.value   = round(p, 3),
    CI        = NA_character_,
    stringsAsFactors = FALSE
  )
}

pool_r2_tbl <- function(mira_obj) {
  r2      <- mice::pool.r.squared(mira_obj, adjusted = TRUE) |> as.data.frame()
  if (!"term" %in% names(r2)) r2 <- cbind(term = rownames(r2), r2)
  est_col <- intersect(c("estimate", "est", "r.squared", "R2"), names(r2))[1]
  if (is.na(est_col))
    stop("Could not find an R² estimate column in pool.r.squared() output.")
  data.frame(
    term      = r2$term,
    estimate  = round(as_num(r2[[est_col]]), 3),
    std.error = if ("std.error" %in% names(r2)) round(as_num(r2[["std.error"]]), 3) else NA_real_,
    statistic = NA_real_,
    df        = NA_character_,
    p.value   = NA_real_,
    CI        = NA_character_,
    stringsAsFactors = FALSE
  )
}

pool_coef_tbl <- function(mira_obj) {
  s  <- summary(pool(mira_obj), conf.int = TRUE) |> as.data.frame()
  CI <- if (all(c("2.5 %", "97.5 %") %in% names(s)))
    paste0("[", round(as_num(s$`2.5 %`), 3), ", ", round(as_num(s$`97.5 %`), 3), "]")
  else if (all(c("conf.low", "conf.high") %in% names(s)))
    paste0("[", round(as_num(s$conf.low),  3), ", ", round(as_num(s$conf.high), 3), "]")
  else NA_character_
  data.frame(
    term      = s$term,
    estimate  = round(as_num(s$estimate),  3),
    std.error = round(as_num(s$std.error), 3),
    statistic = round(as_num(s$statistic), 2),
    df        = as.character(round(as_num(s$df), 1)),
    p.value   = round(as_num(s$p.value), 3),
    CI        = CI,
    stringsAsFactors = FALSE
  )
}

# Helper: print pooled results for a named list of mira objects
print_pooled_results <- function(fits_list, label) {
  cat("\n====", label, "====\n")
  for (nm in names(fits_list)) {
    cat("\n---", nm, "---\n")
    tbl <- rbind(
      pool_overall_F_tbl(fits_list[[nm]]),
      pool_r2_tbl(fits_list[[nm]]),
      pool_coef_tbl(fits_list[[nm]])
    )
    print(tbl)
  }
}

# Print Model 1 results
Model1_fits_classic <- readRDS("Model1_fits_classic.rds")
print_pooled_results(Model1_fits_classic, "Model 1 — Classic ACEs (continuous)")
rm(Model1_fits_classic)


# Section 3.3: Model 2 (sociodemographic + ACE × sex interaction) ----

# --- PWV ---

PWV_Model2_Classic_17_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_17 ~ Classic_ACEs * Child_sex + BP_systolic_17 + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17
  )
)

PWV_Model2_Classic_24_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_24 ~ Classic_ACEs * Child_sex + BP_systolic_24 + Age_24_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_24
  )
)

PWV_Model2_Classic_17_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17,
    subset = (Child_sex == "Female")
  )
)

PWV_Model2_Classic_17_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17,
    subset = (Child_sex == "Male")
  )
)

PWV_Model2_Classic_24_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_24 ~ Classic_ACEs + BP_systolic_24 + Age_24_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_24,
    subset = (Child_sex == "Female")
  )
)

PWV_Model2_Classic_24_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_24 ~ Classic_ACEs + BP_systolic_24 + Age_24_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_24,
    subset = (Child_sex == "Male")
  )
)

# --- cIMT ---

cIMT_Model2_Classic_17_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_17 ~ Classic_ACEs * Child_sex + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17
  )
)

cIMT_Model2_Classic_24_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_24 ~ Classic_ACEs * Child_sex + Age_24_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_24
  )
)

cIMT_Model2_Classic_17_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_17 ~ Classic_ACEs + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17,
    subset = (Child_sex == "Female")
  )
)

cIMT_Model2_Classic_17_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_17 ~ Classic_ACEs + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17,
    subset = (Child_sex == "Male")
  )
)

cIMT_Model2_Classic_24_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_24 ~ Classic_ACEs + Age_24_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_24,
    subset = (Child_sex == "Female")
  )
)

cIMT_Model2_Classic_24_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_24 ~ Classic_ACEs + Age_24_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_24,
    subset = (Child_sex == "Male")
  )
)

# --- Arterial distensibility ---

Dist_Model2_Classic_17_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    Arterial_Dist_17 ~ Classic_ACEs * Child_sex + BP_systolic_17 + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17
  )
)

Dist_Model2_Classic_17_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    Arterial_Dist_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17,
    subset = (Child_sex == "Female")
  )
)

Dist_Model2_Classic_17_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    Arterial_Dist_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17,
    subset = (Child_sex == "Male")
  )
)

Model2_fits_classic <- list(
  PWV_17_Whole  = PWV_Model2_Classic_17_Whole,
  PWV_17_F      = PWV_Model2_Classic_17_F,
  PWV_17_M      = PWV_Model2_Classic_17_M,
  PWV_24_Whole  = PWV_Model2_Classic_24_Whole,
  PWV_24_F      = PWV_Model2_Classic_24_F,
  PWV_24_M      = PWV_Model2_Classic_24_M,
  cIMT_17_Whole = cIMT_Model2_Classic_17_Whole,
  cIMT_17_F     = cIMT_Model2_Classic_17_F,
  cIMT_17_M     = cIMT_Model2_Classic_17_M,
  cIMT_24_Whole = cIMT_Model2_Classic_24_Whole,
  cIMT_24_F     = cIMT_Model2_Classic_24_F,
  cIMT_24_M     = cIMT_Model2_Classic_24_M,
  Dist_17_Whole = Dist_Model2_Classic_17_Whole,
  Dist_17_F     = Dist_Model2_Classic_17_F,
  Dist_17_M     = Dist_Model2_Classic_17_M
)

saveRDS(Model2_fits_classic, "Model2_fits_classic.rds", compress = "gzip")

rm(PWV_Model2_Classic_17_Whole,  PWV_Model2_Classic_17_F,  PWV_Model2_Classic_17_M,
   PWV_Model2_Classic_24_Whole,  PWV_Model2_Classic_24_F,  PWV_Model2_Classic_24_M,
   cIMT_Model2_Classic_17_Whole, cIMT_Model2_Classic_17_F, cIMT_Model2_Classic_17_M,
   cIMT_Model2_Classic_24_Whole, cIMT_Model2_Classic_24_F, cIMT_Model2_Classic_24_M,
   Dist_Model2_Classic_17_Whole, Dist_Model2_Classic_17_F, Dist_Model2_Classic_17_M)

Model2_fits_classic <- readRDS("Model2_fits_classic.rds")
print_pooled_results(Model2_fits_classic, "Model 2 — Classic ACEs (continuous)")
rm(Model2_fits_classic)


# Section 3.4: Model 3 (+ lifestyle / behavioural covariates) ----

# --- PWV ---

PWV_Model3_Classic_17_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15,
    subset = (PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

PWV_Model3_Classic_24_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_24 ~ Classic_ACEs + BP_systolic_24 + Age_24_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_24 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15,
    subset = (PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

PWV_Model3_Classic_17_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15,
    subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

PWV_Model3_Classic_17_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15,
    subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

PWV_Model3_Classic_24_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_24 ~ Classic_ACEs + BP_systolic_24 + Age_24_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_24 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15,
    subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

PWV_Model3_Classic_24_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_24 ~ Classic_ACEs + BP_systolic_24 + Age_24_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_24 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15,
    subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

# --- cIMT ---

cIMT_Model3_Classic_17_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_17 ~ Classic_ACEs + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15,
    subset = (PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

cIMT_Model3_Classic_24_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_24 ~ Classic_ACEs + Age_24_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_24 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15,
    subset = (PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

cIMT_Model3_Classic_17_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_17 ~ Classic_ACEs + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15,
    subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

cIMT_Model3_Classic_17_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_17 ~ Classic_ACEs + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15,
    subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

cIMT_Model3_Classic_24_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_24 ~ Classic_ACEs + Age_24_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_24 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15,
    subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

cIMT_Model3_Classic_24_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_24 ~ Classic_ACEs + Age_24_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_24 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15,
    subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

# --- Arterial distensibility ---

Dist_Model3_Classic_17_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    Arterial_Dist_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15,
    subset = (PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

Dist_Model3_Classic_17_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    Arterial_Dist_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15,
    subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

Dist_Model3_Classic_17_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    Arterial_Dist_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15,
    subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

Model3_fits_classic <- list(
  PWV_17_Whole  = PWV_Model3_Classic_17_Whole,
  PWV_17_F      = PWV_Model3_Classic_17_F,
  PWV_17_M      = PWV_Model3_Classic_17_M,
  PWV_24_Whole  = PWV_Model3_Classic_24_Whole,
  PWV_24_F      = PWV_Model3_Classic_24_F,
  PWV_24_M      = PWV_Model3_Classic_24_M,
  cIMT_17_Whole = cIMT_Model3_Classic_17_Whole,
  cIMT_17_F     = cIMT_Model3_Classic_17_F,
  cIMT_17_M     = cIMT_Model3_Classic_17_M,
  cIMT_24_Whole = cIMT_Model3_Classic_24_Whole,
  cIMT_24_F     = cIMT_Model3_Classic_24_F,
  cIMT_24_M     = cIMT_Model3_Classic_24_M,
  Dist_17_Whole = Dist_Model3_Classic_17_Whole,
  Dist_17_F     = Dist_Model3_Classic_17_F,
  Dist_17_M     = Dist_Model3_Classic_17_M
)

saveRDS(Model3_fits_classic, "Model3_fits_classic.rds", compress = "gzip")

rm(PWV_Model3_Classic_17_Whole,  PWV_Model3_Classic_17_F,  PWV_Model3_Classic_17_M,
   PWV_Model3_Classic_24_Whole,  PWV_Model3_Classic_24_F,  PWV_Model3_Classic_24_M,
   cIMT_Model3_Classic_17_Whole, cIMT_Model3_Classic_17_F, cIMT_Model3_Classic_17_M,
   cIMT_Model3_Classic_24_Whole, cIMT_Model3_Classic_24_F, cIMT_Model3_Classic_24_M,
   Dist_Model3_Classic_17_Whole, Dist_Model3_Classic_17_F, Dist_Model3_Classic_17_M)

Model3_fits_classic <- readRDS("Model3_fits_classic.rds")
print_pooled_results(Model3_fits_classic, "Model 3 — Classic ACEs (continuous)")
rm(Model3_fits_classic)


# Section 3.5: Model 4 (+ neuroendocrine / inflammatory biomarkers) ----
# Note: IL6_24 is included in 24y models only (not available at 17y).

# --- PWV ---

PWV_Model4_Classic_17_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
      Cortisol_15 + CRP_15
  )
)

PWV_Model4_Classic_24_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_24 ~ Classic_ACEs + BP_systolic_24 + Age_24_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_24 +
      Cortisol_15 + CRP_15 + IL6_24
  )
)

PWV_Model4_Classic_17_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
      Cortisol_15 + CRP_15,
    subset = (Child_sex == "Female")
  )
)

PWV_Model4_Classic_17_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
      Cortisol_15 + CRP_15,
    subset = (Child_sex == "Male")
  )
)

PWV_Model4_Classic_24_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_24 ~ Classic_ACEs + BP_systolic_24 + Age_24_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_24 +
      Cortisol_15 + CRP_15 + IL6_24,
    subset = (Child_sex == "Female")
  )
)

PWV_Model4_Classic_24_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_24 ~ Classic_ACEs + BP_systolic_24 + Age_24_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_24 +
      Cortisol_15 + CRP_15 + IL6_24,
    subset = (Child_sex == "Male")
  )
)

# --- cIMT ---

cIMT_Model4_Classic_17_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_17 ~ Classic_ACEs + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
      Cortisol_15 + CRP_15
  )
)

cIMT_Model4_Classic_24_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_24 ~ Classic_ACEs + Age_24_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_24 +
      Cortisol_15 + CRP_15 + IL6_24
  )
)

cIMT_Model4_Classic_17_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_17 ~ Classic_ACEs + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
      Cortisol_15 + CRP_15,
    subset = (Child_sex == "Female")
  )
)

cIMT_Model4_Classic_17_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_17 ~ Classic_ACEs + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
      Cortisol_15 + CRP_15,
    subset = (Child_sex == "Male")
  )
)

cIMT_Model4_Classic_24_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_24 ~ Classic_ACEs + Age_24_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_24 +
      Cortisol_15 + CRP_15 + IL6_24,
    subset = (Child_sex == "Female")
  )
)

cIMT_Model4_Classic_24_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_24 ~ Classic_ACEs + Age_24_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_24 +
      Cortisol_15 + CRP_15 + IL6_24,
    subset = (Child_sex == "Male")
  )
)

# --- Arterial distensibility ---

Dist_Model4_Classic_17_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    Arterial_Dist_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
      Cortisol_15 + CRP_15
  )
)

Dist_Model4_Classic_17_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    Arterial_Dist_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
      Cortisol_15 + CRP_15,
    subset = (Child_sex == "Female")
  )
)

Dist_Model4_Classic_17_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    Arterial_Dist_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
      Cortisol_15 + CRP_15,
    subset = (Child_sex == "Male")
  )
)

Model4_fits_classic <- list(
  PWV_17_Whole  = PWV_Model4_Classic_17_Whole,
  PWV_17_F      = PWV_Model4_Classic_17_F,
  PWV_17_M      = PWV_Model4_Classic_17_M,
  PWV_24_Whole  = PWV_Model4_Classic_24_Whole,
  PWV_24_F      = PWV_Model4_Classic_24_F,
  PWV_24_M      = PWV_Model4_Classic_24_M,
  cIMT_17_Whole = cIMT_Model4_Classic_17_Whole,
  cIMT_17_F     = cIMT_Model4_Classic_17_F,
  cIMT_17_M     = cIMT_Model4_Classic_17_M,
  cIMT_24_Whole = cIMT_Model4_Classic_24_Whole,
  cIMT_24_F     = cIMT_Model4_Classic_24_F,
  cIMT_24_M     = cIMT_Model4_Classic_24_M,
  Dist_17_Whole = Dist_Model4_Classic_17_Whole,
  Dist_17_F     = Dist_Model4_Classic_17_F,
  Dist_17_M     = Dist_Model4_Classic_17_M
)

rm(PWV_Model4_Classic_17_Whole,  PWV_Model4_Classic_17_F,  PWV_Model4_Classic_17_M,
   PWV_Model4_Classic_24_Whole,  PWV_Model4_Classic_24_F,  PWV_Model4_Classic_24_M,
   cIMT_Model4_Classic_17_Whole, cIMT_Model4_Classic_17_F, cIMT_Model4_Classic_17_M,
   cIMT_Model4_Classic_24_Whole, cIMT_Model4_Classic_24_F, cIMT_Model4_Classic_24_M,
   Dist_Model4_Classic_17_Whole, Dist_Model4_Classic_17_F, Dist_Model4_Classic_17_M)

saveRDS(Model4_fits_classic, "Model4_fits_classic.rds", compress = "gzip")

Model4_fits_classic <- readRDS("Model4_fits_classic.rds")
print_pooled_results(Model4_fits_classic, "Model 4 — Classic ACEs (continuous)")
rm(Model4_fits_classic)


# Section 3.6: Model 5 (full adjustment: Models 3 + 4 combined) ----

# --- PWV ---

PWV_Model5_Classic_17_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15 +
      Cortisol_15 + CRP_15, 
    subset = (PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

PWV_Model5_Classic_24_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_24 ~ Classic_ACEs + BP_systolic_24 + Age_24_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_24 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15 +
      Cortisol_15 + CRP_15 + IL6_24,
    subset = (PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

PWV_Model5_Classic_17_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15 +
      Cortisol_15 + CRP_15,
    subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

PWV_Model5_Classic_17_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15 +
      Cortisol_15 + CRP_15,
    subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

PWV_Model5_Classic_24_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_24 ~ Classic_ACEs + BP_systolic_24 + Age_24_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_24 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15 +
      Cortisol_15 + CRP_15 + IL6_24,
    subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

PWV_Model5_Classic_24_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_24 ~ Classic_ACEs + BP_systolic_24 + Age_24_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_24 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15 +
      Cortisol_15 + CRP_15 + IL6_24,
    subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

# --- cIMT ---

cIMT_Model5_Classic_17_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_17 ~ Classic_ACEs + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15 +
      Cortisol_15 + CRP_15,
    subset = (PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

cIMT_Model5_Classic_24_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_24 ~ Classic_ACEs + Age_24_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_24 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15 +
      Cortisol_15 + CRP_15 + IL6_24,
    subset = (PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

cIMT_Model5_Classic_17_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_17 ~ Classic_ACEs + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15 +
      Cortisol_15 + CRP_15,
    subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

cIMT_Model5_Classic_17_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_17 ~ Classic_ACEs + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15 +
      Cortisol_15 + CRP_15,
    subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

cIMT_Model5_Classic_24_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_24 ~ Classic_ACEs + Age_24_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_24 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15 +
      Cortisol_15 + CRP_15 + IL6_24,
    subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

cIMT_Model5_Classic_24_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_24 ~ Classic_ACEs + Age_24_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_24 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15 +
      Cortisol_15 + CRP_15 + IL6_24,
    subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

# --- Arterial distensibility ---

Dist_Model5_Classic_17_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    Arterial_Dist_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15 +
      Cortisol_15 + CRP_15,
    subset = (PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

Dist_Model5_Classic_17_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    Arterial_Dist_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15 +
      Cortisol_15 + CRP_15,
    subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

Dist_Model5_Classic_17_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    Arterial_Dist_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years +
      Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
      Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
      Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
      Daily_MVPA_15 + Daily_Light_PA_15 + Diet_pattern_13_calib +
      non_milk_sugar_13 + Weekday_sleep_duration_15y +
      Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15 +
      Cortisol_15 + CRP_15,
    subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE)
  )
)

Model5_fits_Classic <- list(
  PWV_17_Whole  = PWV_Model5_Classic_17_Whole,
  PWV_17_F      = PWV_Model5_Classic_17_F,
  PWV_17_M      = PWV_Model5_Classic_17_M,
  PWV_24_Whole  = PWV_Model5_Classic_24_Whole,
  PWV_24_F      = PWV_Model5_Classic_24_F,
  PWV_24_M      = PWV_Model5_Classic_24_M,
  cIMT_17_Whole = cIMT_Model5_Classic_17_Whole,
  cIMT_17_F     = cIMT_Model5_Classic_17_F,
  cIMT_17_M     = cIMT_Model5_Classic_17_M,
  cIMT_24_Whole = cIMT_Model5_Classic_24_Whole,
  cIMT_24_F     = cIMT_Model5_Classic_24_F,
  cIMT_24_M     = cIMT_Model5_Classic_24_M,
  Dist_17_Whole = Dist_Model5_Classic_17_Whole,
  Dist_17_F     = Dist_Model5_Classic_17_F,
  Dist_17_M     = Dist_Model5_Classic_17_M
)

rm(PWV_Model5_Classic_17_Whole,  PWV_Model5_Classic_17_F,  PWV_Model5_Classic_17_M,
   PWV_Model5_Classic_24_Whole,  PWV_Model5_Classic_24_F,  PWV_Model5_Classic_24_M,
   cIMT_Model5_Classic_17_Whole, cIMT_Model5_Classic_17_F, cIMT_Model5_Classic_17_M,
   cIMT_Model5_Classic_24_Whole, cIMT_Model5_Classic_24_F, cIMT_Model5_Classic_24_M,
   Dist_Model5_Classic_17_Whole, Dist_Model5_Classic_17_F, Dist_Model5_Classic_17_M)

saveRDS(Model5_fits_Classic, "Model5_fits_Classic.rds", compress = "gzip")

Model5_fits_classic <- readRDS("Model5_fits_Classic.rds")
print_pooled_results(Model5_fits_classic, "Model 5 — Classic ACEs (continuous)")
rm(Model5_fits_classic)


# Section 4: Categorical ACEs — regression models (Models 1–5) ----
# Identical structure to Section 3 but using Classic_ACEs_cat (0 / 1–3 / 4+)
# as the exposure. The ACE × sex interaction in Model 2 (Whole group) uses
# Classic_ACEs_cat * Child_sex.
#
# ADAPTING FOR EXTENDED ACEs:
#   Replace Classic_ACEs_cat with Ext_ACEs_cat throughout this section.
#   All object names follow the same _CatClassic → _CatExt convention.


# Section 4.1: Categorical Model 1 ----

PWV_Model1_CatClassic_17_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_17 ~ Classic_ACEs_cat + BP_systolic_17 + Age_17_clinic_years
  )
)

PWV_Model1_CatClassic_24_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_24 ~ Classic_ACEs_cat + BP_systolic_24 + Age_24_clinic_years
  )
)

PWV_Model1_CatClassic_17_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_17 ~ Classic_ACEs_cat + BP_systolic_17 + Age_17_clinic_years,
    subset = (Child_sex == "Female")
  )
)

PWV_Model1_CatClassic_17_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_17 ~ Classic_ACEs_cat + BP_systolic_17 + Age_17_clinic_years,
    subset = (Child_sex == "Male")
  )
)

PWV_Model1_CatClassic_24_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_24 ~ Classic_ACEs_cat + BP_systolic_24 + Age_24_clinic_years,
    subset = (Child_sex == "Female")
  )
)

PWV_Model1_CatClassic_24_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    PWV_24 ~ Classic_ACEs_cat + BP_systolic_24 + Age_24_clinic_years,
    subset = (Child_sex == "Male")
  )
)

cIMT_Model1_CatClassic_17_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_17 ~ Classic_ACEs_cat + Age_17_clinic_years
  )
)

cIMT_Model1_CatClassic_24_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_24 ~ Classic_ACEs_cat + Age_24_clinic_years
  )
)

cIMT_Model1_CatClassic_17_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_17 ~ Classic_ACEs_cat + Age_17_clinic_years,
    subset = (Child_sex == "Female")
  )
)

cIMT_Model1_CatClassic_17_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_17 ~ Classic_ACEs_cat + Age_17_clinic_years,
    subset = (Child_sex == "Male")
  )
)

cIMT_Model1_CatClassic_24_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_24 ~ Classic_ACEs_cat + Age_24_clinic_years,
    subset = (Child_sex == "Female")
  )
)

cIMT_Model1_CatClassic_24_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    cIMT_24 ~ Classic_ACEs_cat + Age_24_clinic_years,
    subset = (Child_sex == "Male")
  )
)

Dist_Model1_CatClassic_17_Whole <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    Arterial_Dist_17 ~ Classic_ACEs_cat + BP_systolic_17 + Age_17_clinic_years
  )
)

Dist_Model1_CatClassic_17_F <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    Arterial_Dist_17 ~ Classic_ACEs_cat + BP_systolic_17 + Age_17_clinic_years,
    subset = (Child_sex == "Female")
  )
)

Dist_Model1_CatClassic_17_M <- with(
  Study1_Imp_Classic_Trns_Reduced,
  lm(
    Arterial_Dist_17 ~ Classic_ACEs_cat + BP_systolic_17 + Age_17_clinic_years,
    subset = (Child_sex == "Male")
  )
)

Model1_fits_CatClassic <- list(
  PWV_17_Whole  = PWV_Model1_CatClassic_17_Whole,
  PWV_17_F      = PWV_Model1_CatClassic_17_F,
  PWV_17_M      = PWV_Model1_CatClassic_17_M,
  PWV_24_Whole  = PWV_Model1_CatClassic_24_Whole,
  PWV_24_F      = PWV_Model1_CatClassic_24_F,
  PWV_24_M      = PWV_Model1_CatClassic_24_M,
  cIMT_17_Whole = cIMT_Model1_CatClassic_17_Whole,
  cIMT_17_F     = cIMT_Model1_CatClassic_17_F,
  cIMT_17_M     = cIMT_Model1_CatClassic_17_M,
  cIMT_24_Whole = cIMT_Model1_CatClassic_24_Whole,
  cIMT_24_F     = cIMT_Model1_CatClassic_24_F,
  cIMT_24_M     = cIMT_Model1_CatClassic_24_M,
  Dist_17_Whole = Dist_Model1_CatClassic_17_Whole,
  Dist_17_F     = Dist_Model1_CatClassic_17_F,
  Dist_17_M     = Dist_Model1_CatClassic_17_M
)

saveRDS(Model1_fits_CatClassic, "Model1_fits_CatClassic.rds", compress = "gzip")

rm(PWV_Model1_CatClassic_17_Whole,  PWV_Model1_CatClassic_17_F,  PWV_Model1_CatClassic_17_M,
   PWV_Model1_CatClassic_24_Whole,  PWV_Model1_CatClassic_24_F,  PWV_Model1_CatClassic_24_M,
   cIMT_Model1_CatClassic_17_Whole, cIMT_Model1_CatClassic_17_F, cIMT_Model1_CatClassic_17_M,
   cIMT_Model1_CatClassic_24_Whole, cIMT_Model1_CatClassic_24_F, cIMT_Model1_CatClassic_24_M,
   Dist_Model1_CatClassic_17_Whole, Dist_Model1_CatClassic_17_F, Dist_Model1_CatClassic_17_M)

Model1_fits_CatClassic <- readRDS("Model1_fits_CatClassic.rds")
print_pooled_results(Model1_fits_CatClassic, "Model 1 — Classic ACEs (categorical)")
rm(Model1_fits_CatClassic)


# Section 4.2: Categorical Models 2–5 ----
# Models 2–5 follow the same covariate structure as Sections 3.3–3.6 above,
# with Classic_ACEs replaced by Classic_ACEs_cat throughout, and
# Classic_ACEs_cat * Child_sex in the Whole-group Model 2 interaction.
# Code is omitted here for brevity but follows identically from the continuous
# ACEs models above. Replace exposure term and update object/file names.
# As in Section 3, only Models 3 and 5 restrict to
# PA_include_15 == TRUE & Diet_include_13 == TRUE via `subset =`; Models 2
# and 4 use the full available (unfiltered) sample.
#
# Example for Model 3 (filtered), PWV, whole group:
#
#   PWV_Model3_CatClassic_17_Whole <- with(
#     Study1_Imp_Classic_Trns_Reduced,
#     lm(
#       PWV_17 ~ Classic_ACEs_cat + BP_systolic_17 + Age_17_clinic_years +
#         Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
#         Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
#         Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17 +
#         Weekend_sleep_duration_15y + Child_alc_15 + Child_smoke_15,
#       subset = (PA_include_15 == TRUE & Diet_include_13 == TRUE)
#     )
#   )
#
# Example for Model 2 (unfiltered, with ACE × sex interaction), PWV, whole group:
#
#   PWV_Model2_CatClassic_17_Whole <- with(
#     Study1_Imp_Classic_Trns_Reduced,
#     lm(
#       PWV_17 ~ Classic_ACEs_cat * Child_sex + BP_systolic_17 + Age_17_clinic_years +
#         Child_ethnicity + Townsend_sum + Marital_status + Parent_edu +
#         Mat_PND_gest + Mat_age_delivery + Birth_weight_kg +
#         Mat_preg_smoke + Mat_preg_alc + Family_CVD + Age_PHV_c + BMI_17
#     )
#   )
#
# And so on for Models 4 and 5 following the same pattern as their continuous
# ACEs counterparts in Sections 3.5–3.6.


# Section 5: Model assumption checks (single imputation) ----
# Checks are run on imputation 1 only as a visual diagnostic.
# These plots assess linearity, homoscedasticity, and normality of residuals
# for the continuous ACEs models.

imp1 <- complete(Study1_Imp_Classic_Trns_Reduced, action = 1)

# --- Model 1: PWV ---

PWV_m1_17 <- lm(PWV_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years,
                 data = imp1)
ggplot(imp1, aes(Classic_ACEs, PWV_17)) +
  geom_point() + geom_smooth(method = "lm")
check_model(PWV_m1_17)

PWV_m1_24 <- lm(PWV_24 ~ Classic_ACEs + BP_systolic_24 + Age_24_clinic_years,
                 data = imp1)
ggplot(imp1, aes(Classic_ACEs, PWV_24)) +
  geom_point() + geom_smooth(method = "lm")
check_model(PWV_m1_24)

PWV_m1_17_F <- lm(PWV_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years,
                   data = imp1, subset = (Child_sex == "Female"))
ggplot(subset(imp1, Child_sex == "Female"), aes(Classic_ACEs, PWV_17)) +
  geom_point() + geom_smooth(method = "lm")
check_model(PWV_m1_17_F)

PWV_m1_17_M <- lm(PWV_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years,
                   data = imp1, subset = (Child_sex == "Male"))
ggplot(subset(imp1, Child_sex == "Male"), aes(Classic_ACEs, PWV_17)) +
  geom_point() + geom_smooth(method = "lm")
check_model(PWV_m1_17_M)

PWV_m1_24_F <- lm(PWV_24 ~ Classic_ACEs + BP_systolic_24 + Age_24_clinic_years,
                   data = imp1, subset = (Child_sex == "Female"))
ggplot(subset(imp1, Child_sex == "Female"), aes(Classic_ACEs, PWV_24)) +
  geom_point() + geom_smooth(method = "lm")
check_model(PWV_m1_24_F)

PWV_m1_24_M <- lm(PWV_24 ~ Classic_ACEs + BP_systolic_24 + Age_24_clinic_years,
                   data = imp1, subset = (Child_sex == "Male"))
ggplot(subset(imp1, Child_sex == "Male"), aes(Classic_ACEs, PWV_24)) +
  geom_point() + geom_smooth(method = "lm")
check_model(PWV_m1_24_M)

# --- Model 1: cIMT ---

cIMT_m1_17 <- lm(cIMT_17 ~ Classic_ACEs + Age_17_clinic_years, data = imp1)
ggplot(imp1, aes(Classic_ACEs, cIMT_17)) +
  geom_point() + geom_smooth(method = "lm")
check_model(cIMT_m1_17)

cIMT_m1_24 <- lm(cIMT_24 ~ Classic_ACEs + Age_24_clinic_years, data = imp1)
ggplot(imp1, aes(Classic_ACEs, cIMT_24)) +
  geom_point() + geom_smooth(method = "lm")
check_model(cIMT_m1_24)

cIMT_m1_17_F <- lm(cIMT_17 ~ Classic_ACEs + Age_17_clinic_years,
                    data = imp1, subset = (Child_sex == "Female"))
ggplot(subset(imp1, Child_sex == "Female"), aes(Classic_ACEs, cIMT_17)) +
  geom_point() + geom_smooth(method = "lm")
check_model(cIMT_m1_17_F)

cIMT_m1_17_M <- lm(cIMT_17 ~ Classic_ACEs + Age_17_clinic_years,
                    data = imp1, subset = (Child_sex == "Male"))
ggplot(subset(imp1, Child_sex == "Male"), aes(Classic_ACEs, cIMT_17)) +
  geom_point() + geom_smooth(method = "lm")
check_model(cIMT_m1_17_M)

cIMT_m1_24_F <- lm(cIMT_24 ~ Classic_ACEs + Age_24_clinic_years,
                    data = imp1, subset = (Child_sex == "Female"))
ggplot(subset(imp1, Child_sex == "Female"), aes(Classic_ACEs, cIMT_24)) +
  geom_point() + geom_smooth(method = "lm")
check_model(cIMT_m1_24_F)

cIMT_m1_24_M <- lm(cIMT_24 ~ Classic_ACEs + Age_24_clinic_years,
                    data = imp1, subset = (Child_sex == "Male"))
ggplot(subset(imp1, Child_sex == "Male"), aes(Classic_ACEs, cIMT_24)) +
  geom_point() + geom_smooth(method = "lm")
check_model(cIMT_m1_24_M)

# --- Model 1: Arterial distensibility ---

Dist_m1_17 <- lm(Arterial_Dist_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years,
                  data = imp1)
ggplot(imp1, aes(Classic_ACEs, Arterial_Dist_17)) +
  geom_point() + geom_smooth(method = "lm")
check_model(Dist_m1_17)

Dist_m1_17_F <- lm(Arterial_Dist_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years,
                    data = imp1, subset = (Child_sex == "Female"))
ggplot(subset(imp1, Child_sex == "Female"), aes(Classic_ACEs, Arterial_Dist_17)) +
  geom_point() + geom_smooth(method = "lm")
check_model(Dist_m1_17_F)

Dist_m1_17_M <- lm(Arterial_Dist_17 ~ Classic_ACEs + BP_systolic_17 + Age_17_clinic_years,
                    data = imp1, subset = (Child_sex == "Male"))
ggplot(subset(imp1, Child_sex == "Male"), aes(Classic_ACEs, Arterial_Dist_17)) +
  geom_point() + geom_smooth(method = "lm")
check_model(Dist_m1_17_M)
