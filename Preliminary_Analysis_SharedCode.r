# File header ----
# ACEs and Cardiovascular Disease Risk Factors — Preliminary Analysis
# Study:   Adverse Childhood Experiences and Changes in Vascular Health from
#          Adolescence to Early Adulthood: Cross-Sectional and Longitudinal
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
# Author:  Laura Macro
# Date:    2026


# Section 1: Setup ----

# Load all required packages
library(haven)
library(labelled)
library(mice)
library(ggmice)
library(ggplot2)
library(dplyr)
library(tidyr)
library(psych)
library(forcats)
library(stringr)
library(sandwich)
library(lmtest)
library(car)

# Set your working directory and load the ALSPAC dataset here.
# The raw data object is referred to as `raw_data` throughout this script.
# Users must have approved ALSPAC access and load their own extract.

setwd("Your_directory_path_here")  # <-- Update this to your working directory
raw_data <- haven::read_spss("yout_ALSPAC_data_extract.sav")  # <-- Update this to your ALSPAC data file


# Section 2: Variable dictionary ----
# This section documents the mapping between ALSPAC variable codes and the
# analysis variable names used throughout this script.

# --- Exposures ---
# Classic_ACEs   : Classic ACE score            (ALSPAC: clon122)
# Extended_ACEs  : Extended ACE score            (ALSPAC: clon120)

# --- Outcomes @ 17y ---
# PWV_17         : Pulse wave velocity           (ALSPAC: FJAR083d)
# Arterial_Dist_17: Arterial distensibility      (ALSPAC: FJAR077d)
# cIMT_17        : Carotid intima-media thickness(ALSPAC: FJAR079d)

# --- Outcomes @ 24y ---
# PWV_24         : Pulse wave velocity           (ALSPAC: FKCV4200)
# cIMT_24        : Carotid intima-media thickness(ALSPAC: FKCV1131)

# --- Covariates ---
# Family_CVD              : Family history of CVD                         (ALSPAC: FJAR049)
# Parent_edu              : Parental education, NQF level                 (ALSPAC: xa520)
# Age_13_clinic           : Age at 13y clinic visit                       (ALSPAC: fh0011a)
# Age_15_clinic           : Age at 15y clinic visit                       (ALSPAC: fg0011a)
# Age_17_clinic           : Age at 17y clinic visit                       (ALSPAC: FJ003a)
# Age_24_clinic           : Age at 24y clinic visit                       (ALSPAC: FKAR0010)
# BP_systolic_17          : Systolic blood pressure @ 17y                 (ALSPAC: FJAR020a)
# BP_diastolic_17         : Diastolic blood pressure @ 17y                (ALSPAC: FJAR020b)
# BP_systolic_24          : Systolic blood pressure @ 24y                 (ALSPAC: FKBP1030)
# BP_diastolic_24         : Diastolic blood pressure @ 24y                (ALSPAC: FKBP1031)
# Child_sex               : Child sex at birth                            (ALSPAC: kz021)
# Child_ethnicity         : Child ethnic background                       (ALSPAC: c804)
# Marital_status          : Marital status of mother at birth             (ALSPAC: a525)
# Mat_PND_gest            : Maternal postnatal depression, 18wk gestation (ALSPAC: b370)
# Mat_age_delivery        : Maternal age at delivery                      (ALSPAC: mz028b)
# Birth_weight_grams      : Child birth weight (grams)                    (ALSPAC: kz030)
# Mat_preg_smoke          : Maternal smoking during pregnancy             (ALSPAC: b670)
# Mat_preg_alc            : Maternal alcohol use during pregnancy         (ALSPAC: b721)
# Townsend_Gest           : Townsend deprivation score, gestation         (ALSPAC: bTownsendq5)
# Townsend_7 to _15       : Townsend deprivation score at each follow-up  (ALSPAC: f7–tf3Townsendq5)
# Daily_MVPA_15           : Mean daily MVPA, 15y                          (ALSPAC: fh5110)
# Daily_Light_PA_15       : Mean daily light PA, 15y                      (ALSPAC: fh5108)
# Valid_Days_Wk_15        : Valid accelerometry weekdays (>=600 min), 15y (ALSPAC: fh5011)
# Valid_Days_We_15        : Valid accelerometry weekends (>600 min), 15y  (ALSPAC: fh5012)
# Cortisol_15             : Plasma cortisol, 15y                          (ALSPAC: Cortisol_TF3)
# CRP_15                  : C-reactive protein, 15y                       (ALSPAC: crp_TF3)
# IL6_24                  : Interleukin-6, 24y                            (ALSPAC: IL6_F24)
# Sleep_weekday_15        : Sleep duration weekdays, 15y (derived)
# Sleep_weekend_15        : Sleep duration weekends, 15y (derived)
# Age_PHV                 : Age at peak height velocity                   (ALSPAC: clon063)
# Diet_pattern_13         : RRR dietary pattern z-score, 13y              (ALSPAC: fg1600)
# Non_milk_sugar_13       : Non-milk extrinsic sugars, 13y                (ALSPAC: fg1555)
# Days_diet_report_13     : Days of dietary data reported, 13y            (ALSPAC: fg1400)
# Diet_report_plausible_13: Diet reporting validity, 13y                  (ALSPAC: fg1590)
# Child_alc_15            : Child alcohol use in past 30 days, 15y        (ALSPAC: fh8542)
# Child_smoke_15          : Child smoking status (binary), 15y            (ALSPAC: fh8456)
# Child_alc_24            : Past-year alcohol consumption, 24y            (ALSPAC: FKAL1021)
# Child_smoke_24          : Child smoking status (binary), 24y            (ALSPAC: FKSM1090)
# BMI_17                  : Body mass index, 17y                          (ALSPAC: FJMR022a)
# BMI_24                  : Body mass index, 24y                          (ALSPAC: FKMS1040)
# Trig_17 / Trig_24       : Triglycerides, 17y / 24y                     (ALSPAC: TRIG_TF4 / Trig_F24)
# HDL_17 / HDL_24         : HDL cholesterol, 17y / 24y                   (ALSPAC: HDL_TF4 / HDL_F24)
# Glucose_17 / Glucose_24 : Plasma glucose, 17y / 24y                    (ALSPAC: glucose_TF4 / Glucose_F24)
# Insulin_17 / Insulin_24 : Plasma insulin, 17y / 24y                    (ALSPAC: insulin_TF4 / Insulin_F24)
# Child_ICD10_15          : ICD-10 depression diagnosis, 15y              (ALSPAC: fh6892)
# Child_ICD10_17          : ICD-10 depression diagnosis, 17y              (ALSPAC: FJCI1001)
# Child_GAD_15            : Generalised anxiety disorder, 15y             (ALSPAC: FJCI602)
# Child_GAD_24            : Generalised anxiety disorder, 24y             (ALSPAC: FKDQ1030)


# Section 3: Data tidying ----

# Convert SPSS user-missing values to NA
data_na <- user_na_to_na(raw_data)


# --- Derive sleep duration variables ---

# 15y: simple addition of hours and minutes
data_na$Weekday_sleep_duration_15y <- round(data_na$fh5440 + data_na$fh5441 / 60, 2)
data_na$Weekend_sleep_duration_15y <- round(data_na$fh5480 + data_na$fh5481 / 60, 2)

# 13y: wake and bed times are separate; convert to sleep duration
# Step 1 – combine hours and minutes into decimal time
data_na$Weekday_Wakeup_13y  <- round(data_na$ku340a + data_na$ku340b / 60, 2)
data_na$Weekday_Bedtime_13y <- round(data_na$ku341a + data_na$ku341b / 60, 2)
data_na$Weekend_Wakeup_13y  <- round(data_na$ku342a + data_na$ku342b / 60, 2)
data_na$Weekend_Bedtime_13y <- round(data_na$ku343a + data_na$ku343b / 60, 2)

# Step 2 – convert bedtime to 24-h format (assume PM) and wakeup to 24-h (assume AM)
bedtime_24_wkday_13y <- ifelse(data_na$Weekday_Bedtime_13y < 12,
                               data_na$Weekday_Bedtime_13y + 12,
                               data_na$Weekday_Bedtime_13y)
bedtime_24_wkend_13y <- ifelse(data_na$Weekend_Bedtime_13y < 12,
                               data_na$Weekend_Bedtime_13y + 12,
                               data_na$Weekend_Bedtime_13y)
wakeup_24_wkday_13y  <- ifelse(data_na$Weekday_Wakeup_13y == 12, 0, data_na$Weekday_Wakeup_13y)
wakeup_24_wkend_13y  <- ifelse(data_na$Weekend_Wakeup_13y == 12, 0, data_na$Weekend_Wakeup_13y)

# Step 3 – calculate sleep duration (hours)
data_na$Weekday_sleep_duration_13y <- (24 - bedtime_24_wkday_13y) + wakeup_24_wkday_13y
data_na$Weekend_sleep_duration_13y <- (24 - bedtime_24_wkend_13y) + wakeup_24_wkend_13y


# --- Transform smoking / alcohol variables to categorical ---

data_na <- data_na %>%
  mutate(
    Child_smoke_15 = case_when(
      fh8456 > 0    ~ "Yes",
      fh8456 == -8  ~ "No",
      fh8456 == 0   ~ "No",
      fh8456 == 7710 ~ "Yes",
      fh8456 == 7750 ~ "Yes",
      fh8456 == 7770 ~ "Yes",
      TRUE ~ NA_character_
    ),
    Child_smoke_24 = case_when(
      FKSM1090 > 0  ~ "Yes",
      FKSM1090 == -6 ~ "Yes",
      FKSM1090 == -5 ~ "No",
      FKSM1090 == -4 ~ "No",
      FKSM1090 == -3 ~ "No",
      TRUE ~ NA_character_
    ),
    Child_alc_24 = case_when(
      FKAL1021 == -4 ~ "None",
      FKAL1021 == -3 ~ "None",
      FKAL1021 == 0  ~ "1 or 2",
      FKAL1021 == 1  ~ "3 or 4",
      FKAL1021 == 2  ~ "5 or 6",
      FKAL1021 == 3  ~ "7 to 9",
      FKAL1021 == 4  ~ "10+",
      TRUE ~ NA_character_
    ),
    Child_alc_15 = case_when(
      fh8542 == 1 ~ "None",
      fh8542 == 2 ~ "1 or 2",
      fh8542 == 3 ~ "3 to 5",
      fh8542 == 4 ~ "6 to 9",
      fh8542 == 5 ~ "10 to 19",
      fh8542 == 6 ~ "20 to 39",
      fh8542 == 7 ~ "40 to 99",
      fh8542 == 8 ~ "100+",
      TRUE ~ NA_character_
    ),
    across(c("Child_smoke_15", "Child_smoke_24", "Child_alc_24", "Child_alc_15"), as.factor)
  )


# --- Rename variables to analysis names ---

data_renamed <- data_na %>%
  rename(
    # Exposures
    Classic_ACEs          = clon122,
    Extended_ACEs         = clon120,
    # Outcomes
    PWV_17                = FJAR083d,
    Arterial_Dist_17      = FJAR077d,
    cIMT_17               = FJAR079d,
    PWV_24                = FKCV4200,
    cIMT_24               = FKCV1131,
    # CVD risk factors
    BMI_17                = FJMR022a,
    Trig_17               = TRIG_TF4,
    HDL_17                = HDL_TF4,
    Glucose_17            = glucose_TF4,
    Insulin_17            = insulin_TF4,
    BMI_24                = FKMS1040,
    Trig_24               = Trig_F24,
    HDL_24                = HDL_F24,
    Glucose_24            = Glucose_F24,
    Insulin_24            = Insulin_F24,
    # Covariates
    Family_CVD            = FJAR049,
    Parent_edu            = xa520,
    Age_13_clinic         = fh0011a,
    Age_15_clinic         = fg0011a,
    Age_17_clinic         = FJ003a,
    Age_24_clinic         = FKAR0010,
    BP_systolic_17        = FJAR020a,
    BP_diastolic_17       = FJAR020b,
    BP_systolic_24        = FKBP1030,
    BP_diastolic_24       = FKBP1031,
    Child_sex             = kz021,
    Child_ethnicity       = c804,
    Marital_status        = a525,
    Mat_PND_gest          = b370,
    Mat_age_delivery      = mz028b,
    Birth_weight_grams    = kz030,
    Mat_preg_smoke        = b670,
    Mat_preg_alc          = b721,
    Townsend_Gest         = bTownsendq5,
    Townsend_7            = f7Townsendq5,
    Townsend_9            = f9Townsendq5,
    Townsend_10           = f10Townsendq5,
    Townsend_11           = f11Townsendq5,
    Townsend_12           = tf1Townsendq5,
    Townsend_13           = tf2Townsendq5,
    Townsend_15           = tf3Townsendq5,
    Daily_MVPA_15         = fh5110,
    Daily_Light_PA_15     = fh5108,
    Valid_Days_Wk_15      = fh5011,
    Valid_Days_We_15      = fh5012,
    Cortisol_15           = Cortisol_TF3,
    CRP_15                = crp_TF3,
    IL6_24                = IL6_F24,
    Age_PHV               = clon063,
    Diet_pattern_13       = fg1600,
    Non_milk_sugar_13     = fg1555,
    Days_diet_report_13   = fg1400,
    Diet_report_plausible_13 = fg1590,
    # Mental health
    Child_ICD10_15        = fh6892,
    Child_ICD10_17        = FJCI1001,
    Child_GAD_15          = FJCI602,
    Child_GAD_24          = FKDQ1030
  )


# --- Convert SPSS labelled variables to R factor / numeric types ---

data_converted <- data_renamed %>%
  mutate(

    # Exposures (continuous)
    Classic_ACEs          = as.numeric(zap_labels(Classic_ACEs)),
    Extended_ACEs         = as.numeric(zap_labels(Extended_ACEs)),

    # Vascular outcomes (continuous)
    PWV_17                = as.numeric(zap_labels(PWV_17)),
    Arterial_Dist_17      = as.numeric(zap_labels(Arterial_Dist_17)),
    cIMT_17               = as.numeric(zap_labels(cIMT_17)),
    PWV_24                = as.numeric(zap_labels(PWV_24)),
    cIMT_24               = as.numeric(zap_labels(cIMT_24)),

    # CVD risk factors (continuous)
    BMI_17                = as.numeric(zap_labels(BMI_17)),
    Insulin_17            = as.numeric(zap_labels(Insulin_17)),
    Glucose_17            = as.numeric(zap_labels(Glucose_17)),
    Trig_17               = as.numeric(zap_labels(Trig_17)),
    HDL_17                = as.numeric(zap_labels(HDL_17)),
    BMI_24                = as.numeric(zap_labels(BMI_24)),
    Insulin_24            = as.numeric(zap_labels(Insulin_24)),
    Glucose_24            = as.numeric(zap_labels(Glucose_24)),
    Trig_24               = as.numeric(zap_labels(Trig_24)),
    HDL_24                = as.numeric(zap_labels(HDL_24)),

    # Binary nominal covariates (reference level first)
    Family_CVD            = fct_relevel(haven::as_factor(Family_CVD), "No", "Yes"),
    Child_ethnicity       = fct_relevel(haven::as_factor(Child_ethnicity), "White", "Non-white"),
    Child_sex             = fct_relevel(haven::as_factor(Child_sex), "Male", "Female", "Not known"),

    # Mental health: ICD-10 depression
    Child_ICD10_15_num    = as.numeric(Child_ICD10_15),
    Child_ICD10_15        = case_when(
      Child_ICD10_15_num == 2 ~ "No Depression",
      Child_ICD10_15_num == 1 ~ "Yes Depression",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("No Depression", "Yes Depression")),

    Child_ICD10_17_num    = as.numeric(Child_ICD10_17),
    Child_ICD10_17        = case_when(
      Child_ICD10_17_num == 0 ~ "No Depression",
      Child_ICD10_17_num == 1 ~ "Yes Depression",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("No Depression", "Yes Depression")),

    # Mental health: GAD
    Child_GAD_15_num      = as.numeric(Child_GAD_15),
    Child_GAD_15          = case_when(
      Child_GAD_15_num == 0 ~ "No",
      Child_GAD_15_num == 1 ~ "Yes",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("No", "Yes")),

    Child_GAD_24_num      = as.numeric(Child_GAD_24),
    Child_GAD_24          = case_when(
      Child_GAD_24_num == 0 ~ "No",
      Child_GAD_24_num == 1 ~ "Yes",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("No", "Yes")),

    # Ordinal categorical covariates
    Parent_edu            = factor(
      fct_relevel(haven::as_factor(Parent_edu),
                  "Level 1", "Level 2", "Level 3", "Level 4", "Level 5"),
      ordered = TRUE),

    Mat_preg_smoke        = factor(
      fct_relevel(haven::as_factor(Mat_preg_smoke),
                  "0", "1-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30+"),
      ordered = TRUE),

    Mat_preg_alc          = factor(
      fct_relevel(haven::as_factor(Mat_preg_alc),
                  "never", "<1 glass PWK", "1+ glasses PWK",
                  "1-2 glasses PDAY", "3-9 glasses PDAY", "10+ glasses PDAY", "DK"),
      ordered = TRUE),

    Townsend_Gest         = factor(haven::as_factor(Townsend_Gest), ordered = TRUE),
    Townsend_7            = factor(haven::as_factor(Townsend_7),    ordered = TRUE),
    Townsend_9            = factor(haven::as_factor(Townsend_9),    ordered = TRUE),
    Townsend_10           = factor(haven::as_factor(Townsend_10),   ordered = TRUE),
    Townsend_11           = factor(haven::as_factor(Townsend_11),   ordered = TRUE),
    Townsend_12           = factor(haven::as_factor(Townsend_12),   ordered = TRUE),
    Townsend_13           = factor(haven::as_factor(Townsend_13),   ordered = TRUE),
    Townsend_15           = factor(haven::as_factor(Townsend_15),   ordered = TRUE),

    Diet_report_plausible_13 = haven::as_factor(Diet_report_plausible_13),

    # Nominal categorical (unordered)
    Marital_status        = fct_relevel(
      haven::as_factor(Marital_status),
      "Never married", "Widowed", "Divorced", "Separated",
      "1st marriage", "Marriage 2 or 3"),

    # Continuous covariates
    Age_13_clinic         = as.numeric(zap_labels(Age_13_clinic)),
    Age_15_clinic         = as.numeric(zap_labels(Age_15_clinic)),
    Age_17_clinic         = as.numeric(zap_labels(Age_17_clinic)),
    Age_24_clinic         = as.numeric(zap_labels(Age_24_clinic)),
    BP_systolic_17        = as.numeric(zap_labels(BP_systolic_17)),
    BP_diastolic_17       = as.numeric(zap_labels(BP_diastolic_17)),
    BP_systolic_24        = as.numeric(zap_labels(BP_systolic_24)),
    BP_diastolic_24       = as.numeric(zap_labels(BP_diastolic_24)),
    Mat_PND_gest          = as.numeric(zap_labels(Mat_PND_gest)),
    Mat_age_delivery      = as.numeric(zap_labels(Mat_age_delivery)),
    Birth_weight_grams    = as.numeric(zap_labels(Birth_weight_grams)),
    Daily_MVPA_15         = as.numeric(zap_labels(Daily_MVPA_15)),
    Daily_Light_PA_15     = as.numeric(zap_labels(Daily_Light_PA_15)),
    Valid_Days_Wk_15      = as.numeric(zap_labels(Valid_Days_Wk_15)),
    Valid_Days_We_15      = as.numeric(zap_labels(Valid_Days_We_15)),
    Cortisol_15           = as.numeric(zap_labels(Cortisol_15)),
    CRP_15                = as.numeric(zap_labels(CRP_15)),
    IL6_24                = as.numeric(zap_labels(IL6_24)),
    Age_PHV               = as.numeric(zap_labels(Age_PHV)),
    Diet_pattern_13       = as.numeric(zap_labels(Diet_pattern_13)),
    Non_milk_sugar_13     = as.numeric(zap_labels(Non_milk_sugar_13)),
    Days_diet_report_13   = as.numeric(zap_labels(Days_diet_report_13))
  )


# --- Dietary intake calibration for reporting bias ---
# Fit calibration model: adjust dietary pattern score by reporting validity
# (under-reporters and over-reporters are corrected to the valid-reporter mean)

Diet_calib_fit <- lm(
  Diet_pattern_13 ~ stats::relevel(as.factor(Diet_report_plausible_13), ref = "Valid reporter"),
  data = data_converted
)

beta_under <- coef(Diet_calib_fit)["stats::relevel(as.factor(Diet_report_plausible_13), ref = \"Valid reporter\")Under reporter"]
beta_over  <- coef(Diet_calib_fit)["stats::relevel(as.factor(Diet_report_plausible_13), ref = \"Valid reporter\")Over reporter"]


# --- Helper: sex-stratified age-adjusted standardised residuals ---
# Used to construct the cardiometabolic risk (CMR) score

std_resid <- function(formula, data) {
  mf <- model.frame(formula, data = data, na.action = na.pass)
  ok <- stats::complete.cases(mf)
  if (sum(ok) < 2) return(rep(NA_real_, nrow(data)))
  fit <- stats::lm(formula, data = data[ok, , drop = FALSE])
  out <- rep(NA_real_, nrow(data))
  r   <- stats::residuals(fit)
  s   <- stats::sd(r, na.rm = TRUE)
  out[ok] <- if (is.na(s) || s == 0) 0 else (r - mean(r, na.rm = TRUE)) / s
  out
}


# --- Derive cardiometabolic risk scores and additional variables ---

data_transformed <- data_converted %>%
  mutate(Child_sex_num = as.numeric(Child_sex) - 1) %>%
  group_by(Child_sex) %>%
  group_modify(~ {
    .x %>% mutate(
      # Age-adjusted standardised residuals @ 17y
      Z_SBP_17      = std_resid(BP_systolic_17 ~ Age_17_clinic, .x),
      Z_BMI_17      = std_resid(BMI_17         ~ Age_17_clinic, .x),
      Z_insulin_17  = std_resid(Insulin_17     ~ Age_17_clinic, .x),
      Z_trig_17     = std_resid(Trig_17        ~ Age_17_clinic, .x),
      Z_glucose_17  = std_resid(Glucose_17     ~ Age_17_clinic, .x),
      Z_HDL_17      = std_resid(HDL_17         ~ Age_17_clinic, .x),
      Z_HDL_inv_17  = -1 * Z_HDL_17,
      # CMR @ 17y: NA unless all components present
      CMR_17 = ifelse(
        stats::complete.cases(Z_SBP_17, Z_BMI_17, Z_insulin_17,
                              Z_trig_17, Z_glucose_17, Z_HDL_inv_17),
        Z_SBP_17 + Z_BMI_17 + Z_insulin_17 + Z_trig_17 + Z_glucose_17 + Z_HDL_inv_17,
        NA_real_
      ),
      # Age-adjusted standardised residuals @ 24y
      Z_SBP_24      = std_resid(BP_systolic_24 ~ Age_24_clinic, .x),
      Z_BMI_24      = std_resid(BMI_24         ~ Age_24_clinic, .x),
      Z_insulin_24  = std_resid(Insulin_24     ~ Age_24_clinic, .x),
      Z_trig_24     = std_resid(Trig_24        ~ Age_24_clinic, .x),
      Z_glucose_24  = std_resid(Glucose_24     ~ Age_24_clinic, .x),
      Z_HDL_24      = std_resid(HDL_24         ~ Age_24_clinic, .x),
      Z_HDL_inv_24  = -1 * Z_HDL_24,
      # CMR @ 24y
      CMR_24 = ifelse(
        stats::complete.cases(Z_SBP_24, Z_BMI_24, Z_insulin_24,
                              Z_trig_24, Z_glucose_24, Z_HDL_inv_24),
        Z_SBP_24 + Z_BMI_24 + Z_insulin_24 + Z_trig_24 + Z_glucose_24 + Z_HDL_inv_24,
        NA_real_
      )
    )
  }) %>%
  ungroup() %>%
  mutate(
    # Calibrated dietary pattern score
    Diet_pattern_13_calib = case_when(
      Diet_report_plausible_13 == "Under reporter" ~ Diet_pattern_13 - beta_under,
      Diet_report_plausible_13 == "Over reporter"  ~ Diet_pattern_13 - beta_over,
      TRUE ~ Diet_pattern_13
    ),

    # ACE exposure categories
    Ext_ACEs_cat = factor(
      cut(Extended_ACEs, breaks = c(-Inf, 0, 3, Inf),
          labels = c("0", "1-3", "4+"), include.lowest = TRUE),
      levels = c("0", "1-3", "4+")),

    Ext_ACEs_binary = factor(
      ifelse(Extended_ACEs > 0, "Yes", "No"), levels = c("No", "Yes")),

    Classic_ACEs_Cat = factor(
      cut(Classic_ACEs, breaks = c(-Inf, 0, 3, Inf),
          labels = c("0", "1-3", "4+"), include.lowest = TRUE),
      levels = c("0", "1-3", "4+")),

    Classic_ACEs_binary = factor(
      ifelse(Classic_ACEs > 0, "Yes", "No"), levels = c("No", "Yes")),

    # Accelerometry inclusion flags
    PA_include_15   = Valid_Days_Wk_15 >= 3 & Valid_Days_We_15 >= 1,
    Diet_include_13 = Days_diet_report_13 >= 3,

    # Unit conversions
    Birth_weight_kg     = Birth_weight_grams / 1000,
    Age_17_clinic_years = Age_17_clinic / 12,
    Age_24_clinic_years = Age_24_clinic / 12,

    # Summed Townsend deprivation index across all time points
    Townsend_sum = as.numeric(Townsend_Gest) + as.numeric(Townsend_7) +
      as.numeric(Townsend_9)  + as.numeric(Townsend_10) +
      as.numeric(Townsend_11) + as.numeric(Townsend_12) +
      as.numeric(Townsend_13) + as.numeric(Townsend_15)
  )

# Restrict to participants with known sex and ALSPAC core membership
# kz011b is the ALSPAC "core" sample indicator (singletons alive at 1 year)
data_transformed <- data_transformed %>%
  filter(!is.na(Child_sex)) %>%
  mutate(kz011b = haven::as_factor(kz011b)) %>%
  filter(haven::as_factor(kz011b) == "Yes")

# Centred age variables for interaction / sensitivity analyses
# Centred at 17 years (clinic visit reference age)
data_transformed$age17_c <- data_transformed$Age_17_clinic_years - 17
data_transformed$age24_c <- data_transformed$Age_24_clinic_years - 17

# Save processed dataset
saveRDS(data_transformed, file = "data_transformed.rds")

rm(data_converted, Diet_calib_fit)


# Section 4: Missing data patterns ----

# Visual missingness patterns for key exposure–outcome pairs
plot_pattern(data_transformed, c("Classic_ACEs", "PWV_17", "Arterial_Dist_17", "cIMT_17"), rotate = TRUE)
plot_pattern(data_transformed, c("Extended_ACEs", "PWV_17", "Arterial_Dist_17", "cIMT_17"), rotate = TRUE)
plot_pattern(data_transformed, c("Classic_ACEs",  "PWV_24", "cIMT_24"), rotate = TRUE)
plot_pattern(data_transformed, c("Extended_ACEs", "PWV_24", "cIMT_24"), rotate = TRUE)


# Section 5: Missing data summary table ----

# Variable sets
vars_exposures       <- c("Classic_ACEs", "Extended_ACEs")
vascular_outcomes_17 <- c("PWV_17", "Arterial_Dist_17", "cIMT_17")
vascular_outcomes_24 <- c("PWV_24", "cIMT_24")
vars_outcomes_CVD    <- c("BMI_17", "Insulin_17", "Glucose_17", "Trig_17", "HDL_17",
                          "BMI_24", "Insulin_24", "Glucose_24", "Trig_24", "HDL_24")
vars_covariates      <- c("Child_sex", "Age_17_clinic_years", "Age_24_clinic_years",
                          "BP_systolic_17", "BP_systolic_24", "BP_diastolic_17", "BP_diastolic_24",
                          "Child_ethnicity", "Marital_status", "Mat_PND_gest", "Mat_age_delivery",
                          "Birth_weight_kg", "Mat_preg_smoke", "Mat_preg_alc", "Townsend_sum",
                          "Age_PHV", "Family_CVD", "Parent_edu", "BMI_17", "BMI_24",
                          "Daily_MVPA_15", "Daily_Light_PA_15", "Diet_pattern_13_calib",
                          "Non_milk_sugar_13", "Weekday_sleep_duration_15y",
                          "Weekend_sleep_duration_15y", "Child_alc_15", "Child_smoke_15",
                          "Child_alc_24", "Child_smoke_24", "Cortisol_15", "CRP_15",
                          "IL6_24", "age17_c", "age24_c")

# Function: % and n missing per variable
get_missing_summary <- function(data, vars) {
  data.frame(
    Variable        = vars,
    Percent_Missing = colMeans(is.na(data[, vars])) * 100,
    N_Missing       = colSums(is.na(data[, vars])),
    row.names       = NULL
  )
}

# Summarise by variable block
combined_summary_missingdata <- rbind(
  cbind(Block = "Classic/Extended ACEs", get_missing_summary(data_transformed, vars_exposures)),
  cbind(Block = "Outcomes @ 17",         get_missing_summary(data_transformed, vascular_outcomes_17)),
  cbind(Block = "Outcomes @ 24",         get_missing_summary(data_transformed, vascular_outcomes_24)),
  cbind(Block = "CVD Outcomes",          get_missing_summary(data_transformed, vars_outcomes_CVD)),
  cbind(Block = "Covariates",            get_missing_summary(data_transformed, vars_covariates))
)

print(combined_summary_missingdata)
rm(combined_summary_missingdata)


# Section 6: Descriptive statistics ----

# Variable groupings
exposures    <- c("Classic_ACEs", "Extended_ACEs", "Classic_ACEs_Cat", "Classic_ACEs_binary")
cov_socio    <- c("Family_CVD", "Parent_edu", "Child_ethnicity", "Marital_status",
                  "Mat_PND_gest", "Mat_age_delivery", "Birth_weight_kg",
                  "Mat_preg_smoke", "Mat_preg_alc", "Townsend_sum")
cov_physio   <- c("Child_sex", "Age_17_clinic_years", "Age_24_clinic_years",
                  "age17_c", "age24_c", "BP_systolic_17", "BP_diastolic_17",
                  "BMI_17", "BMI_24", "BP_systolic_24", "BP_diastolic_24",
                  "Cortisol_15", "CRP_15", "IL6_24", "Age_PHV")
cov_behav    <- c("Daily_MVPA_15", "Daily_Light_PA_15", "Diet_pattern_13",
                  "Non_milk_sugar_13", "Weekday_sleep_duration_15y",
                  "Weekend_sleep_duration_15y", "Child_alc_15", "Child_smoke_15",
                  "Child_alc_24", "Child_smoke_24")

# Descriptive statistics for all variable blocks
psych::describe(data_transformed %>% dplyr::select(all_of(exposures)))
psych::describe(data_transformed %>% dplyr::select(all_of(c(vascular_outcomes_17, vascular_outcomes_24))))
psych::describe(data_transformed %>% dplyr::select(where(is.numeric) & all_of(cov_socio)))
psych::describe(data_transformed %>% dplyr::select(where(is.numeric) & all_of(cov_physio)))
psych::describe(data_transformed %>% dplyr::select(where(is.numeric) & all_of(cov_behav)))

# Histograms / bar charts grouped by variable type
# NOTE: histograms_or_bars() is a custom helper function that produces
# histograms for numeric variables and bar charts for factors.
# Define this function in your environment before running these lines.
# Example signature: histograms_or_bars(data, vars, title)

var_groups <- list(
  exposures    = exposures,
  outcomes_17  = vascular_outcomes_17,
  outcomes_24  = vascular_outcomes_24,
  CVD_outcomes = vars_outcomes_CVD,
  cov_socio    = cov_socio,
  cov_physio   = cov_physio,
  cov_behav    = cov_behav
)

for (grp_name in names(var_groups)) {
  png(file = paste0(grp_name, "_histograms_and_boxplots_%03d.png"),
      width = 3200, height = 1600, res = 200)
  par(mfrow = c(3, 3))
  histograms_or_bars(data_transformed, var_groups[[grp_name]], grp_name)
  dev.off()
}

# Note: visual inspection identified variables needing transformation:
#   Log transform:  Daily_MVPA_15, CRP_15, Non_milk_sugar_13, Valid_Days_Wk_15, Valid_Days_We_15
#   Sqrt transform: Mat_PND_gest, Classic_ACEs, Extended_ACEs (zero-inflated)


# Section 7: Complete-case and missing-case descriptive statistics ----

# Create missingness indicators for each exposure–outcome analysis
data_transformed <- data_transformed %>%
  mutate(
    y17_complete_classic    = factor(complete.cases(select(., all_of(c("Classic_ACEs",  vascular_outcomes_17)))),
                                     levels = c(FALSE, TRUE), labels = c("incomplete", "complete")),
    y17_complete_extended   = factor(complete.cases(select(., all_of(c("Extended_ACEs", vascular_outcomes_17)))),
                                     levels = c(FALSE, TRUE), labels = c("incomplete", "complete")),
    y24_complete_classic    = factor(complete.cases(select(., all_of(c("Classic_ACEs",  vascular_outcomes_24)))),
                                     levels = c(FALSE, TRUE), labels = c("incomplete", "complete")),
    y24_complete_extended   = factor(complete.cases(select(., all_of(c("Extended_ACEs", vascular_outcomes_24)))),
                                     levels = c(FALSE, TRUE), labels = c("incomplete", "complete")),

    y17_cov_complete_classic  = factor(complete.cases(select(., all_of(c("Classic_ACEs",  vascular_outcomes_17, vars_covariates)))),
                                       levels = c(FALSE, TRUE), labels = c("incomplete", "complete")),
    y17_cov_complete_extended = factor(complete.cases(select(., all_of(c("Extended_ACEs", vascular_outcomes_17, vars_covariates)))),
                                       levels = c(FALSE, TRUE), labels = c("incomplete", "complete")),
    y24_cov_complete_classic  = factor(complete.cases(select(., all_of(c("Classic_ACEs",  vascular_outcomes_24, vars_covariates)))),
                                       levels = c(FALSE, TRUE), labels = c("incomplete", "complete")),
    y24_cov_complete_extended = factor(complete.cases(select(., all_of(c("Extended_ACEs", vascular_outcomes_24, vars_covariates)))),
                                       levels = c(FALSE, TRUE), labels = c("incomplete", "complete")),

    y17_missing_classic    = factor(!complete.cases(select(., all_of(c("Classic_ACEs",  vascular_outcomes_17)))),
                                    levels = c(FALSE, TRUE), labels = c("not missing", "missing")),
    y17_missing_extended   = factor(!complete.cases(select(., all_of(c("Extended_ACEs", vascular_outcomes_17)))),
                                    levels = c(FALSE, TRUE), labels = c("not missing", "missing")),
    y24_missing_classic    = factor(!complete.cases(select(., all_of(c("Classic_ACEs",  vascular_outcomes_24)))),
                                    levels = c(FALSE, TRUE), labels = c("not missing", "missing")),
    y24_missing_extended   = factor(!complete.cases(select(., all_of(c("Extended_ACEs", vascular_outcomes_24)))),
                                    levels = c(FALSE, TRUE), labels = c("not missing", "missing")),

    y17_cov_missing_classic   = factor(!complete.cases(select(., all_of(c("Classic_ACEs",  vascular_outcomes_17, vars_covariates)))),
                                       levels = c(FALSE, TRUE), labels = c("not missing", "missing")),
    y17_cov_missing_extended  = factor(!complete.cases(select(., all_of(c("Extended_ACEs", vascular_outcomes_17, vars_covariates)))),
                                       levels = c(FALSE, TRUE), labels = c("not missing", "missing")),
    y24_cov_missing_classic   = factor(!complete.cases(select(., all_of(c("Classic_ACEs",  vascular_outcomes_24, vars_covariates)))),
                                       levels = c(FALSE, TRUE), labels = c("not missing", "missing")),
    y24_cov_missing_extended  = factor(!complete.cases(select(., all_of(c("Extended_ACEs", vascular_outcomes_24, vars_covariates)))),
                                       levels = c(FALSE, TRUE), labels = c("not missing", "missing"))
  )

saveRDS(data_transformed, file = "data_transformed.rds")

# Subset datasets for complete and missing cases
complete_classic_17     <- data_transformed %>% filter(y17_complete_classic    == "complete")
complete_ext_17         <- data_transformed %>% filter(y17_complete_extended   == "complete")
cov_complete_classic_17 <- data_transformed %>% filter(y17_cov_complete_classic  == "complete")
cov_complete_ext_17     <- data_transformed %>% filter(y17_cov_complete_extended == "complete")

complete_classic_24     <- data_transformed %>% filter(y24_complete_classic    == "complete")
complete_ext_24         <- data_transformed %>% filter(y24_complete_extended   == "complete")
cov_complete_classic_24 <- data_transformed %>% filter(y24_cov_complete_classic  == "complete")
cov_complete_ext_24     <- data_transformed %>% filter(y24_cov_complete_extended == "complete")

missing_classic_17 <- data_transformed %>% filter(y17_missing_classic  == "missing")
missing_classic_24 <- data_transformed %>% filter(y24_missing_classic  == "missing")
missing_ext_17     <- data_transformed %>% filter(y17_missing_extended == "missing")
missing_ext_24     <- data_transformed %>% filter(y24_missing_extended == "missing")


# --- Descriptive summaries for complete-case and missing-case subsets ---

all_vars <- c(vars_exposures, vars_covariates, vascular_outcomes_17, vascular_outcomes_24)
num_vars <- all_vars[sapply(data_transformed[all_vars], is.numeric)]

datasets <- list(
  complete_classic_17 = complete_classic_17, complete_ext_17 = complete_ext_17,
  complete_classic_24 = complete_classic_24, complete_ext_24 = complete_ext_24,
  cov_complete_classic_17 = cov_complete_classic_17, cov_complete_ext_17 = cov_complete_ext_17,
  cov_complete_classic_24 = cov_complete_classic_24, cov_complete_ext_24 = cov_complete_ext_24,
  missing_classic_17 = missing_classic_17, missing_ext_17 = missing_ext_17,
  missing_classic_24 = missing_classic_24, missing_ext_24 = missing_ext_24
)

for (ds_name in names(datasets)) {
  cat("\n--- Numeric descriptives:", ds_name, "---\n")
  print(psych::describe(datasets[[ds_name]][num_vars]))
}

rm(datasets)


# --- Histograms / bar charts for complete-case and missing-case subsets ---

subset_plot_list <- list(
  cov_complete_classic_17 = cov_complete_classic_17,
  cov_complete_ext_17     = cov_complete_ext_17,
  cov_complete_classic_24 = cov_complete_classic_24,
  cov_complete_ext_24     = cov_complete_ext_24,
  complete_classic_17     = complete_classic_17,
  complete_ext_17         = complete_ext_17,
  complete_classic_24     = complete_classic_24,
  complete_ext_24         = complete_ext_24
)

plot_var_groups <- list(
  exposures   = vars_exposures,
  outcomes_17 = vascular_outcomes_17,
  outcomes_24 = vascular_outcomes_24,
  cov_socio   = cov_socio,
  cov_physio  = cov_physio,
  cov_behav   = cov_behav
)

for (ds_name in names(subset_plot_list)) {
  ds <- subset_plot_list[[ds_name]]
  for (grp_name in names(plot_var_groups)) {
    png(file = paste0(grp_name, "_Plots_", ds_name, "_%03d.png"),
        width = 3200, height = 1600, res = 200)
    par(mfrow = c(3, 3))
    histograms_or_bars(ds, plot_var_groups[[grp_name]], grp_name)
    dev.off()
  }
}


# Section 8: Log / sqrt transformations ----
# Variables identified as non-normally distributed via histogram inspection:
#   Log transform:  CRP_15, Non_milk_sugar_13
#   Sqrt transform: Mat_PND_gest, Classic_ACEs, Extended_ACEs, Daily_MVPA_15

data_transformed <- data_transformed %>%
  mutate(
    across(all_of(c("CRP_15", "Non_milk_sugar_13")),
           ~ ifelse(is.na(.x), NA, log(.x)), .names = "ln_{.col}"),
    across(all_of(c("Mat_PND_gest", "Classic_ACEs", "Extended_ACEs", "Daily_MVPA_15")),
           ~ ifelse(is.na(.x), NA, sqrt(.x)), .names = "sqrt_{.col}")
  )

transformed_variables <- c("ln_CRP_15", "ln_Non_milk_sugar_13",
                            "sqrt_Mat_PND_gest", "sqrt_Classic_ACEs",
                            "sqrt_Extended_ACEs", "sqrt_Daily_MVPA_15")

# Descriptive statistics for transformed variables
descriptives_transformed <- psych::describe(
  data_transformed %>% dplyr::select(all_of(transformed_variables))
)
print(descriptives_transformed)

# Visual check of transformed distributions
png("transformed_vars_plots%03d.png", width = 1600, height = 800, res = 150)
print(histograms_or_bars(data_transformed, transformed_variables, "Transformed Variables"))
dev.off()


# Section 9: Helper functions for complete-case association analyses ----
# Used in Sections 10 and 11 below.

# Recode marital status to binary (married vs not married)
recode_marital_status_bin <- function(df,
                                      var = "Marital_status",
                                      new_var = "Marital_status_bin",
                                      married_levels     = c("1st marriage", "Marriage 2 or 3"),
                                      not_married_levels = c("Never married", "Divorced", "Widowed_Separated"),
                                      married_patterns     = c("marriage", "married"),
                                      not_married_patterns = c("never", "divorc", "widow", "separ")) {
  if (!var %in% names(df)) return(df)
  x_chr   <- as.character(df[[var]])
  x_clean <- stringr::str_squish(x_chr)
  out     <- rep(NA_character_, length(x_clean))
  out[x_clean %in% married_levels]     <- "Married"
  out[x_clean %in% not_married_levels] <- "Not married"
  still <- is.na(out) & !is.na(x_clean)
  if (any(still)) {
    lower      <- stringr::str_to_lower(x_clean[still])
    is_married <- Reduce(`|`, lapply(married_patterns,     function(p) stringr::str_detect(lower, p)))
    is_not     <- Reduce(`|`, lapply(not_married_patterns, function(p) stringr::str_detect(lower, p)))
    out[still][is_married & !is_not] <- "Married"
    out[still][is_not]               <- "Not married"
  }
  df[[new_var]] <- factor(out, levels = c("Not married", "Married"))
  df
}

# Apply to complete-case datasets
cov_complete_classic_17 <- recode_marital_status_bin(cov_complete_classic_17)
cov_complete_classic_24 <- recode_marital_status_bin(cov_complete_classic_24)
cov_complete_ext_17     <- recode_marital_status_bin(cov_complete_ext_17)
cov_complete_ext_24     <- recode_marital_status_bin(cov_complete_ext_24)

# Collapse rare levels before regression (prevents fitting issues with sparse cells)
collapse_rare_levels <- function(x, min_level_n = 5, other_level = "Other") {
  if (is.character(x) || is.logical(x)) x <- factor(x)
  if (!is.factor(x)) return(list(x = x, collapsed = FALSE))
  if (nlevels(droplevels(x)) < 2) return(list(x = x, collapsed = FALSE))
  before <- levels(droplevels(x))
  x2     <- forcats::fct_lump_min(x, min = min_level_n, other_level = other_level)
  after  <- levels(droplevels(x2))
  if (other_level %in% levels(x2)) x2 <- forcats::fct_relevel(x2, other_level, after = Inf)
  list(x = x2, collapsed = !identical(before, after))
}

# Robust variance–covariance: HC3 with HC2 fallback
robust_vcov <- function(fit) {
  V  <- tryCatch(sandwich::vcovHC(fit, type = "HC3"),
                 warning = function(w) sandwich::vcovHC(fit, type = "HC2"),
                 error   = function(e) sandwich::vcovHC(fit, type = "HC2"))
  vc <- "HC3"
  if (any(!is.finite(diag(V)))) { V <- sandwich::vcovHC(fit, type = "HC2"); vc <- "HC2" }
  list(V = V, vcov_type = vc)
}

# p-value formatting helper
fmt_p <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)))
}


# Section 10: Covariate–exposure associations ----
# For each covariate, fits: ACE_count ~ covariate
# Continuous predictors: beta (mean difference in ACE count per unit), HC3 robust SE
# Binary factors: beta + 95% CI; multi-level factors: global robust F-test p-value only

run_cc_covariate_models_global <- function(df, exposure, covars,
                                           conf = 0.95, min_n = 20, min_level_n = 5,
                                           other_level = "Other") {
  if (!exposure %in% names(df)) stop("Exposure not found: ", exposure)

  lapply(covars, function(cv) {
    if (!cv %in% names(df))
      return(data.frame(Covariate = cv, N = NA_integer_, Type = "",  Beta = "",
                        `95% CI` = "", `p-value` = "", Notes = "variable missing",
                        check.names = FALSE))

    dat <- df[!is.na(df[[exposure]]) & !is.na(df[[cv]]), ]
    dat <- data.frame(y = dat[[exposure]], x_raw = dat[[cv]])
    n <- nrow(dat)
    if (n < min_n)
      return(data.frame(Covariate = cv, N = n, Type = "", Beta = "", `95% CI` = "",
                        `p-value` = "", Notes = paste0("insufficient n (<", min_n, ")"),
                        check.names = FALSE))

    x         <- dat$x_raw
    lumped    <- collapse_rare_levels(x, min_level_n = min_level_n, other_level = other_level)
    x         <- lumped$x
    lump_note <- if (is.factor(x) && lumped$collapsed)
      paste0("Rare levels collapsed (<", min_level_n, " -> ", other_level, ")") else ""

    if (is.factor(x)) {
      x <- droplevels(x)
      if (length(table(x)) < 2 || any(table(x) < min_level_n))
        return(data.frame(Covariate = cv, N = n, Type = "Categorical (skipped)", Beta = "",
                          `95% CI` = "", `p-value` = "",
                          Notes = paste0("sparse factor levels (min n < ", min_level_n, ")",
                                         if (lump_note != "") paste0("; ", lump_note) else ""),
                          check.names = FALSE))
    }

    dat2 <- data.frame(y = dat$y, x = x)
    fit  <- lm(y ~ x, data = dat2)
    rv   <- robust_vcov(fit); V <- rv$V; vc <- rv$vcov_type
    notes <- paste(c(lump_note,
                     if (vc == "HC2") "SEs use HC2 (HC3 unstable)" else "")[
                       c(lump_note, if (vc == "HC2") "SEs use HC2 (HC3 unstable)" else "") != ""],
                   collapse = "; ")

    if (!is.factor(dat2$x)) {
      ct   <- lmtest::coeftest(fit, vcov. = V)
      est  <- unname(ct["x", 1]); se <- unname(ct["x", 2]); p <- unname(ct["x", 4])
      crit <- qt(1 - (1 - conf) / 2, df = fit$df.residual)
      return(data.frame(Covariate = cv, N = n, Type = "Continuous",
                        Beta = formatC(est, format = "f", digits = 3),
                        `95% CI` = paste0(formatC(est - crit * se, format = "f", digits = 3), ", ",
                                          formatC(est + crit * se, format = "f", digits = 3)),
                        `p-value` = fmt_p(p), Notes = notes, check.names = FALSE))
    }

    k        <- nlevels(dat2$x)
    lh       <- car::linearHypothesis(fit, matchCoefs(fit, "x"), vcov = V, test = "F")
    p_global <- lh$`Pr(>F)`[2]

    if (k == 2) {
      ct     <- lmtest::coeftest(fit, vcov. = V)
      x_term <- rownames(ct)[grepl("^x", rownames(ct))][1]
      est    <- unname(ct[x_term, 1]); se <- unname(ct[x_term, 2])
      crit   <- qt(1 - (1 - conf) / 2, df = fit$df.residual)
      return(data.frame(Covariate = cv, N = n, Type = "Categorical (binary)",
                        Beta = formatC(est, format = "f", digits = 3),
                        `95% CI` = paste0(formatC(est - crit * se, format = "f", digits = 3), ", ",
                                          formatC(est + crit * se, format = "f", digits = 3)),
                        `p-value` = fmt_p(p_global), Notes = notes, check.names = FALSE))
    }

    data.frame(Covariate = cv, N = n, Type = paste0("Categorical (", k, " levels)"),
               Beta = "", `95% CI` = "", `p-value` = fmt_p(p_global), Notes = notes,
               check.names = FALSE)
  }) |> do.call(what = rbind)
}

# Covariates for association analyses
covars <- c(
  "Age_17_clinic_years", "Age_24_clinic_years", "age17_c", "age24_c",
  "BP_systolic_17", "BP_systolic_24", "BP_diastolic_17", "BP_diastolic_24",
  "Child_ethnicity", "Marital_status_bin", "Mat_PND_gest", "Mat_age_delivery",
  "Birth_weight_kg", "Mat_preg_smoke", "Mat_preg_alc", "Townsend_sum",
  "Age_PHV", "Family_CVD", "Parent_edu", "BMI_17", "BMI_24",
  "Daily_MVPA_15", "Daily_Light_PA_15", "Diet_pattern_13_calib", "Non_milk_sugar_13",
  "Weekday_sleep_duration_15y", "Weekend_sleep_duration_15y",
  "Child_alc_15", "Child_smoke_15", "Child_alc_24", "Child_smoke_24",
  "Cortisol_15", "CRP_15", "IL6_24"
)

# Run for Classic ACEs
cat("\n--- Classic ACEs ~ covariates (17y) ---\n")
print(run_cc_covariate_models_global(cov_complete_classic_17, "Classic_ACEs", covars))
cat("\n--- Classic ACEs ~ covariates (24y) ---\n")
print(run_cc_covariate_models_global(cov_complete_classic_24, "Classic_ACEs", covars))

# Run for Extended ACEs
cat("\n--- Extended ACEs ~ covariates (17y) ---\n")
print(run_cc_covariate_models_global(cov_complete_ext_17, "Extended_ACEs", covars, min_n = 15))
cat("\n--- Extended ACEs ~ covariates (24y) ---\n")
print(run_cc_covariate_models_global(cov_complete_ext_24, "Extended_ACEs", covars, min_n = 15))


# Section 11: Covariate–outcome associations ----
# For each covariate, fits: outcome ~ covariate
# Same robust regression approach as Section 10.

run_cc_covariate_models_global_y <- function(df, outcome, covars,
                                              conf = 0.95, min_n = 20, min_level_n = 5,
                                              other_level = "Other") {
  if (!outcome %in% names(df)) stop("Outcome not found: ", outcome)

  lapply(covars, function(cv) {
    if (!cv %in% names(df))
      return(data.frame(Outcome = outcome, Covariate = cv, N = NA_integer_, Type = "",
                        Beta = "", `95% CI` = "", `p-value` = "", Notes = "variable missing",
                        check.names = FALSE))

    dat <- df[!is.na(df[[outcome]]) & !is.na(df[[cv]]), ]
    dat <- data.frame(y = dat[[outcome]], x_raw = dat[[cv]])
    n <- nrow(dat)
    if (n < min_n)
      return(data.frame(Outcome = outcome, Covariate = cv, N = n, Type = "", Beta = "",
                        `95% CI` = "", `p-value` = "",
                        Notes = paste0("insufficient n (<", min_n, ")"),
                        check.names = FALSE))

    x         <- dat$x_raw
    lumped    <- collapse_rare_levels(x, min_level_n = min_level_n, other_level = other_level)
    x         <- lumped$x
    lump_note <- if (is.factor(x) && lumped$collapsed)
      paste0("Rare levels collapsed (<", min_level_n, " -> ", other_level, ")") else ""

    if (is.factor(x)) {
      x <- droplevels(x)
      if (length(table(x)) < 2 || any(table(x) < min_level_n))
        return(data.frame(Outcome = outcome, Covariate = cv, N = n,
                          Type = "Categorical (skipped)", Beta = "", `95% CI` = "", `p-value` = "",
                          Notes = paste0("sparse factor levels (min n < ", min_level_n, ")",
                                         if (lump_note != "") paste0("; ", lump_note) else ""),
                          check.names = FALSE))
    }

    dat2 <- data.frame(y = dat$y, x = x)
    fit  <- lm(y ~ x, data = dat2)
    rv   <- robust_vcov(fit); V <- rv$V; vc <- rv$vcov_type
    notes <- paste(c(lump_note,
                     if (vc == "HC2") "SEs use HC2 (HC3 unstable)" else "")[
                       c(lump_note, if (vc == "HC2") "SEs use HC2 (HC3 unstable)" else "") != ""],
                   collapse = "; ")

    if (!is.factor(dat2$x)) {
      ct   <- lmtest::coeftest(fit, vcov. = V)
      est  <- unname(ct["x", 1]); se <- unname(ct["x", 2]); p <- unname(ct["x", 4])
      crit <- qt(1 - (1 - conf) / 2, df = fit$df.residual)
      return(data.frame(Outcome = outcome, Covariate = cv, N = n, Type = "Continuous",
                        Beta = formatC(est, format = "f", digits = 3),
                        `95% CI` = paste0(formatC(est - crit * se, format = "f", digits = 3), ", ",
                                          formatC(est + crit * se, format = "f", digits = 3)),
                        `p-value` = fmt_p(p), Notes = notes, check.names = FALSE))
    }

    k        <- nlevels(dat2$x)
    lh       <- car::linearHypothesis(fit, matchCoefs(fit, "x"), vcov = V, test = "F")
    p_global <- lh$`Pr(>F)`[2]

    if (k == 2) {
      ct     <- lmtest::coeftest(fit, vcov. = V)
      x_term <- rownames(ct)[grepl("^x", rownames(ct))][1]
      est    <- unname(ct[x_term, 1]); se <- unname(ct[x_term, 2])
      crit   <- qt(1 - (1 - conf) / 2, df = fit$df.residual)
      return(data.frame(Outcome = outcome, Covariate = cv, N = n,
                        Type = "Categorical (binary)",
                        Beta = formatC(est, format = "f", digits = 3),
                        `95% CI` = paste0(formatC(est - crit * se, format = "f", digits = 3), ", ",
                                          formatC(est + crit * se, format = "f", digits = 3)),
                        `p-value` = fmt_p(p_global), Notes = notes, check.names = FALSE))
    }

    data.frame(Outcome = outcome, Covariate = cv, N = n,
               Type = paste0("Categorical (", k, " levels)"),
               Beta = "", `95% CI` = "", `p-value` = fmt_p(p_global), Notes = notes,
               check.names = FALSE)
  }) |> do.call(what = rbind)
}

# Run for each vascular outcome
outcome_map <- data.frame(
  label = c("PWV",  "cIMT",  "AIx"),
  y17   = c("PWV_17", "cIMT_17", "Arterial_Dist_17"),
  y24   = c("PWV_24", "cIMT_24", NA_character_),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(outcome_map))) {
  cat("\n---", outcome_map$label[i], "~ covariates (17y, Classic) ---\n")
  print(run_cc_covariate_models_global_y(cov_complete_classic_17,
                                         outcome_map$y17[i], covars))
  if (!is.na(outcome_map$y24[i])) {
    cat("\n---", outcome_map$label[i], "~ covariates (24y, Classic) ---\n")
    print(run_cc_covariate_models_global_y(cov_complete_classic_24,
                                           outcome_map$y24[i], covars))
  }
}


# Section 12: Missingness association analyses ----
# Tests whether observed covariates and outcomes predict missingness in the
# primary analytic variables (to inform the missing-data mechanism).

# Recode marital status and lump sparse levels in the full dataset
data_transformed_bin <- recode_marital_status_bin(data_transformed)

data_transformed_bin$Child_alc_24 <- data_transformed_bin$Child_alc_24 %>%
  as.factor() %>% fct_drop() %>% forcats::fct_lump_min(min = 5)
data_transformed_bin$Mat_preg_alc <- data_transformed_bin$Mat_preg_alc %>%
  as.factor() %>% fct_drop() %>% forcats::fct_lump_min(min = 5)


# Logistic regression: categorical covariate vs missingness indicator (robust SEs)
test_logit_cat_vs_missing <- function(cat_df, missing_vec) {
  y <- as.numeric(missing_vec == "missing")
  do.call(rbind, lapply(names(cat_df), function(var) {
    x <- as.factor(cat_df[[var]])
    if (nlevels(x) < 2) return(NULL)
    contrasts(x) <- contr.treatment(nlevels(x))
    dat <- data.frame(y = y, x = x)
    fit <- tryCatch(glm(y ~ x, data = dat, family = binomial), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    coefs        <- summary(fit)$coefficients
    robust_coefs <- lmtest::coeftest(fit, vcov = sandwich::vcovHC(fit, type = "HC0"))
    do.call(rbind, lapply(2:nrow(coefs), function(i) {
      data.frame(
        Covariate  = var,
        Level      = levels(x)[i],
        OR         = exp(coefs[i, 1]),
        CI_Lower   = exp(coefs[i, 1] - 1.96 * coefs[i, 2]),
        CI_Upper   = exp(coefs[i, 1] + 1.96 * coefs[i, 2]),
        P_Value    = coefs[i, 4],
        Robust_SE  = robust_coefs[i, 2]
      )
    }))
  }))
}

# Logistic regression: numeric covariate vs missingness indicator (robust SEs)
test_logit_num_vs_missing <- function(num_df, missing_vec) {
  y <- as.numeric(missing_vec == "missing")
  do.call(rbind, lapply(names(num_df), function(var) {
    x <- as.numeric(num_df[[var]])
    if (all(is.na(x))) return(NULL)
    dat <- data.frame(y = y, x = x)
    fit <- tryCatch(glm(y ~ x, data = dat, family = binomial), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    est          <- summary(fit)$coefficients
    robust_coefs <- lmtest::coeftest(fit, vcov = sandwich::vcovHC(fit, type = "HC0"))
    data.frame(
      Covariate = var,
      OR        = exp(est[2, 1]),
      CI_Lower  = exp(est[2, 1] - 1.96 * est[2, 2]),
      CI_Upper  = exp(est[2, 1] + 1.96 * est[2, 2]),
      P_Value   = est[2, 4],
      Robust_SE = robust_coefs[2, 2]
    )
  }))
}


# Define variable selections for missingness analyses
# (covariates split into numeric and categorical for the two logistic functions above)
covariates_num_vars <- c(
  "Child_sex", "Age_17_clinic_years", "Age_24_clinic_years",
  "BP_systolic_17", "BP_systolic_24", "BP_diastolic_17", "BP_diastolic_24",
  "Mat_PND_gest", "Mat_age_delivery", "Birth_weight_kg",
  "Daily_MVPA_15", "Daily_Light_PA_15", "Cortisol_15", "CRP_15", "IL6_24",
  "Weekday_sleep_duration_15y", "Weekend_sleep_duration_15y",
  "Age_PHV", "Diet_pattern_13", "Non_milk_sugar_13", "Townsend_sum"
)

# Helper to run a full missingness analysis for one missingness indicator
run_missingness_analysis <- function(data, missing_var, label) {
  missing_vec <- data[[missing_var]]

  exp_num <- data[, c("Classic_ACEs", "Extended_ACEs")]
  exp_num[] <- lapply(exp_num, function(x) if (is.factor(x)) as.numeric(as.character(x)) else x)

  outcomes_num <- data[, c("PWV_17", "Arterial_Dist_17", "cIMT_17", "PWV_24", "cIMT_24")]
  outcomes_num[] <- lapply(outcomes_num, function(x) if (is.factor(x)) as.numeric(as.character(x)) else x)

  cov_num <- data[, covariates_num_vars]
  cov_num[] <- lapply(cov_num, function(x) if (is.factor(x)) as.numeric(as.character(x)) else x)

  cov_cat_cols <- vars_covariates[sapply(data[vars_covariates],
                                         function(x) is.factor(x) || is.character(x))]
  cov_cat <- data[, cov_cat_cols]
  cov_cat[] <- lapply(cov_cat, as.factor)

  list(
    label    = label,
    cat_cov  = test_logit_cat_vs_missing(cov_cat,      missing_vec),
    num_cov  = test_logit_num_vs_missing(cov_num,      missing_vec),
    outcome  = test_logit_num_vs_missing(outcomes_num, missing_vec),
    exposure = test_logit_num_vs_missing(exp_num,      missing_vec)
  )
}

# Run for all four missingness indicators
missingness_analyses <- list(
  classic_17   = run_missingness_analysis(data_transformed_bin, "y17_missing_classic",  "Classic 17y"),
  classic_24   = run_missingness_analysis(data_transformed_bin, "y24_missing_classic",  "Classic 24y"),
  extended_17  = run_missingness_analysis(data_transformed_bin, "y17_missing_extended", "Extended 17y"),
  extended_24  = run_missingness_analysis(data_transformed_bin, "y24_missing_extended", "Extended 24y")
)

# Print results to console
for (nm in names(missingness_analyses)) {
  res <- missingness_analyses[[nm]]
  cat("\n===", res$label, "===\n")
  cat("-- Categorical covariates --\n"); print(res$cat_cov)
  cat("-- Numeric covariates --\n");     print(res$num_cov)
  cat("-- Outcomes --\n");               print(res$outcome)
  cat("-- Exposures --\n");              print(res$exposure)
}

rm(missingness_analyses, data_transformed_bin)
