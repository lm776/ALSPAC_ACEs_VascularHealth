# File header ----
# ACEs and Vascular Health — Multiple Imputation Setup
# Study:  Adverse Childhood Experiences and Changes in Vascular Health from 
#         Childhood to Mid-Adulthood: Cross-Sectional and Longitudinal Evidence 
#         from the ALSPAC Study
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


# Section 1: Required packages ----

library(haven)
library(dplyr)
library(forcats)
library(ggmice)
library(sandwich)
library(lmtest)
library(labelled)
library(mice)
library(miceadds)


# Section 2: Variable dictionary ----
# Documents the mapping between analysis variable names and ALSPAC variable
# codes used in the MI setup.

# --- Exposures ---
# Classic_ACEs  = clon122
# Extended_ACEs = clon120

# --- Outcomes @ 17 ---
# PWV_17             : Pulse wave velocity @ 17
# Arterial_Dist_17   : Arterial distensibility @ 17
# cIMT_17            : Carotid intima-media thickness @ 17
# BMI_17             : Body mass index @ 17
# Trig_17            : Triglycerides @ 17
# HDL_17             : HDL Cholesterol @ 17
# Glucose_17         : Blood glucose @ 17
# Insulin_17         : Insulin @ 17

# --- Outcomes @ 24 ---
# PWV_24             : Pulse wave velocity @ 24
# cIMT_24            : Carotid intima-media thickness @ 24
# BMI_24             : Body mass index @ 24
# Trig_24            : Triglycerides @ 24
# HDL_24             : HDL Cholesterol @ 24
# Glucose_24         : Blood glucose @ 24
# Insulin_24         : Insulin @ 24

# --- Covariates ---
# Family_CVD                   : Family history of CVD
# Parent_edu                   : Parental education, NQF level
# Age_13_clinic_years          : Age at visit for 13y clinic
# Age_15_clinic_years          : Age at visit for 15y clinic
# Age_17_clinic_years          : Age at visit for 17y clinic
# Age_24_clinic_years          : Age at visit for 24y clinic
# age17_c                      : Centred age at 17 (Age_17_clinic_years - 17)
# age24_c                      : Centred age at 24 (Age_24_clinic_years - 17)
# BP_systolic_17               : Systolic blood pressure @ 17y
# BP_systolic_24               : Systolic blood pressure @ 24y
# BP_diastolic_17              : Diastolic blood pressure @ 17y
# BP_diastolic_24              : Diastolic blood pressure @ 24y
# Child_ethnicity              : Child ethnic background
# Marital_status               : Marital status of mother at birth
# Mat_PND_gest                 : Maternal post-natal depression, 18 weeks gestation
# Mat_age_delivery             : Maternal age at child delivery
# Birth_weight_kg              : Child birth weight (kg)
# Mat_preg_smoke               : Maternal smoking during pregnancy
# Mat_preg_alc                 : Maternal alcohol consumption during pregnancy
# Townsend_sum                 : Townsend score — sum of 7 separate timepoints
# Daily_MVPA_15                : Mean daily moderate-to-vigorous physical activity, 15y
# Daily_Light_PA_15            : Mean daily light physical activity, 15y
# Valid_Days_Wk_15             : Number of valid weekdays (>=600 mins), 15y
# Valid_Days_We_15             : Number of valid weekend days (>=600 mins), 15y
# Cortisol_15                  : Plasma cortisol, 15y
# CRP_15                       : C-reactive protein, 15y
# IL6_24                       : Interleukin-6, 24y
# Weekday_sleep_duration_15y   : Sleep duration weekdays, 15y
# Weekend_sleep_duration_15y   : Sleep duration weekends, 15y
# Age_PHV                      : Age at peak height velocity
# Diet_pattern_13              : RRR dietary pattern z-score, 13y
# Non_milk_sugar_13            : Non-milk extrinsic sugars intake, 13y
# Days_diet_report_13          : Number of days dietary data reported, 13y
# Diet_report_plausible_13     : Valid/under/over reporter of dietary data, 13y
# Child_alc_15                 : Child alcohol consumption/30 days, 15y
# Child_smoke_15               : Child smoking status (binary), 15y
# Child_alc_24                 : Past year alcohol consumption, 24y
# Child_smoke_24               : Child smoking status (binary), 24y
# Child_ICD10_15               : Child ICD-10 depression diagnosis, 15y
# Child_ICD10_17               : Child ICD-10 depression diagnosis, 17y
# Child_GAD_15                 : Child GAD-7 anxiety score, 15y
# Child_GAD_24                 : Child GAD-7 anxiety score, 24y


# Section 3: Variable sets ----

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
                          "Non_milk_sugar_13", "Weekday_sleep_duration_15y", "Weekend_sleep_duration_15y",
                          "Child_alc_15", "Child_smoke_15", "Child_alc_24", "Child_smoke_24",
                          "Cortisol_15", "CRP_15", "IL6_24", "age17_c", "age24_c")


# Section 4: Data preparation ----
# Assumes data_transformed has been loaded from the preliminary analysis script.
# cidB4619 = ALSPAC pregnancy identifier (publicly documented)
# kz011b   = ALSPAC core sample membership flag

# Collapse sparse factor levels identified from contingency table inspection
data_transformed <- data_transformed %>%
  filter(!(Marital_status %in% c("Missing", "Consent withdrawn")))

data_transformed$Marital_status <- droplevels(data_transformed$Marital_status)

data_transformed <- data_transformed %>%
  mutate(Marital_status = fct_collapse(Marital_status,
                                       Widowed_Separated = c("Widowed", "Separated")))

data_transformed <- data_transformed %>%
  mutate(Mat_preg_smoke = fct_collapse(Mat_preg_smoke, "25+" = c("25-29", "30+")))

data_transformed <- data_transformed %>%
  mutate(Mat_preg_alc = fct_collapse(Mat_preg_alc,
                                     "3+ glasses PDAY" = c("3-9 glasses PDAY", "10+ glasses PDAY")))

data_transformed$Child_alc_15 <- droplevels(data_transformed$Child_alc_15)

data_transformed <- data_transformed %>%
  mutate(Child_alc_15 = fct_collapse(Child_alc_15, "40+" = c("40 to 99", "100+")))

# Convert age-13 and age-15 clinic variables from months to years
data_transformed <- data_transformed %>%
  mutate(Age_13_clinic_years = Age_13_clinic / 12,
         Age_15_clinic_years = Age_15_clinic / 12)

# Convert all variables to factor or numeric (from SPSS labelled format)
data_transformed <- data_transformed %>%
  mutate(

    # Classic ACE individual items
    clon100 = haven::as_factor(clon100),
    clon101 = haven::as_factor(clon101),
    clon102 = haven::as_factor(clon102),
    clon103 = haven::as_factor(clon103),
    clon104 = haven::as_factor(clon104),
    clon105 = haven::as_factor(clon105),
    clon106 = haven::as_factor(clon106),
    clon107 = haven::as_factor(clon107),
    clon108 = haven::as_factor(clon108),
    clon109 = haven::as_factor(clon109),

    # Extended ACE individual items
    clon111 = haven::as_factor(clon111),
    clon112 = haven::as_factor(clon112),
    clon113 = haven::as_factor(clon113),
    clon114 = haven::as_factor(clon114),
    clon115 = haven::as_factor(clon115),
    clon116 = haven::as_factor(clon116),
    clon117 = haven::as_factor(clon117),
    clon118 = haven::as_factor(clon118),
    clon119 = haven::as_factor(clon119),

    # ACEs general
    e390    = as.numeric(zap_labels(e390)),
    b650    = haven::as_factor(b650),
    b669    = factor(
      fct_relevel(haven::as_factor(b669),
                  "0", "1-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30+"),
      ordered = TRUE),
    clon166 = haven::as_factor(clon166),
    clon167 = as.numeric(zap_labels(clon167)),
    clon168 = as.numeric(zap_labels(clon168)),
    clon169 = as.numeric(zap_labels(clon169)),
    clon170 = as.numeric(zap_labels(clon170)),
    qlet    = haven::as_factor(qlet),

    # CVD outcomes / auxiliary variables
    bestgest    = as.numeric(zap_labels(bestgest)),
    CRP_f9      = as.numeric(zap_labels(CRP_f9)),
    CRP_TF4     = as.numeric(zap_labels(CRP_TF4)),
    CRP_F24     = as.numeric(zap_labels(CRP_F24)),
    IL6_F9      = as.numeric(zap_labels(IL6_F9)),
    f7ms026a    = as.numeric(zap_labels(f7ms026a)),
    f9ms026a    = as.numeric(zap_labels(f9ms026a)),
    fdms026a    = as.numeric(zap_labels(fdms026a)),
    fems026a    = as.numeric(zap_labels(fems026a)),
    ff2039      = as.numeric(zap_labels(ff2039)),
    fg3139      = as.numeric(zap_labels(fg3139)),
    fh3019      = as.numeric(zap_labels(fh3019)),
    f9ms018     = as.numeric(zap_labels(f9ms018)),
    fdms018     = as.numeric(zap_labels(fdms018)),
    fems018     = as.numeric(zap_labels(fems018)),
    ff2020      = as.numeric(zap_labels(ff2020)),
    fg3120      = as.numeric(zap_labels(fg3120)),
    fh4020      = as.numeric(zap_labels(fh4020)),
    FKMS1052    = as.numeric(zap_labels(FKMS1052)),
    f9dx135     = as.numeric(zap_labels(f9dx135)),
    fedx135     = as.numeric(zap_labels(fedx135)),
    fg3254      = as.numeric(zap_labels(fg3254)),
    insulin_F9  = as.numeric(zap_labels(insulin_F9)),
    insulin_TF3 = as.numeric(zap_labels(insulin_TF3)),
    HDL_f9      = as.numeric(zap_labels(HDL_f9)),
    hdl_TF3     = as.numeric(zap_labels(hdl_TF3)),
    LDL_f9      = as.numeric(zap_labels(LDL_f9)),
    glucose_TF3 = as.numeric(zap_labels(glucose_TF3)),
    trig_f9     = as.numeric(zap_labels(trig_f9)),
    trig_TF3    = as.numeric(zap_labels(trig_TF3)),
    FJAR040     = fct_relevel(haven::as_factor(FJAR040), "No", "Yes"),
    FJAR043     = fct_relevel(haven::as_factor(FJAR043), "No", "Yes"),
    FJAR047     = fct_relevel(haven::as_factor(FJAR047), "No", "Yes"),
    FJAR048     = fct_relevel(haven::as_factor(FJAR048), "No", "Yes"),
    YPB1213     = fct_relevel(haven::as_factor(YPB1213), "No", "Not Sure", "Yes, by Self", "Yes, by a Doctor"),
    YPB1214     = fct_relevel(haven::as_factor(YPB1214), "No", "Not Sure", "Yes, by Self", "Yes, by a Doctor"),
    YPB1230     = fct_relevel(haven::as_factor(YPB1230), "No", "Not Sure", "Yes, by Self", "Yes, by a Doctor"),
    fdar117     = as.numeric(zap_labels(fdar117)),
    fdar118     = as.numeric(zap_labels(fdar118)),
    FJAR020c    = as.numeric(zap_labels(FJAR020c)),
    fdar111     = as.numeric(zap_labels(fdar111)),
    fdar113     = as.numeric(zap_labels(fdar113)),
    fdar114     = as.numeric(zap_labels(fdar114)),
    fdar115     = as.numeric(zap_labels(fdar115)),
    fdar116     = as.numeric(zap_labels(fdar116)),
    fdar119     = as.numeric(zap_labels(fdar119)),
    FKBP1032    = as.numeric(zap_labels(FKBP1032)),

    # MVPA auxiliary variables
    feag210 = as.numeric(zap_labels(feag210)),
    fg1310  = as.numeric(zap_labels(fg1310)),
    fg1302  = as.numeric(zap_labels(fg1302)),
    fg1206  = as.numeric(zap_labels(fg1206)),
    fg1207  = as.numeric(zap_labels(fg1207)),
    fg1308  = as.numeric(zap_labels(fg1308)),
    feag101 = as.numeric(zap_labels(feag101)),
    feag102 = as.numeric(zap_labels(feag102)),
    feag208 = as.numeric(zap_labels(feag208)),

    # Dietary auxiliary variables
    f7dd400 = as.numeric(zap_labels(f7dd400)),
    fddd400 = as.numeric(zap_labels(fddd400)),
    fg1560  = as.numeric(zap_labels(fg1560)),
    f7dd100 = as.numeric(zap_labels(f7dd100)),
    fddd100 = as.numeric(zap_labels(fddd100)),
    f7dd507 = as.numeric(zap_labels(f7dd507)),
    f7dd346 = as.numeric(zap_labels(f7dd346)),
    f7dd600 = as.numeric(zap_labels(f7dd600)),
    fddd506 = as.numeric(zap_labels(fddd506)),
    fddd345 = as.numeric(zap_labels(fddd345)),
    fddd600 = haven::as_factor(fddd600),
    fg1585  = as.numeric(zap_labels(fg1585)),

    # Puberty auxiliary variables
    pub430 = factor(haven::as_factor(pub430), ordered = TRUE),
    pub435 = factor(haven::as_factor(pub435), ordered = TRUE),
    pub450 = factor(haven::as_factor(pub450), ordered = TRUE),
    pub455 = factor(haven::as_factor(pub455), ordered = TRUE),
    pub530 = factor(haven::as_factor(pub530), ordered = TRUE),
    pub535 = factor(haven::as_factor(pub535), ordered = TRUE),
    pub550 = factor(haven::as_factor(pub550), ordered = TRUE),
    pub555 = factor(haven::as_factor(pub555), ordered = TRUE),
    pub630 = factor(haven::as_factor(pub630), ordered = TRUE),
    pub635 = factor(haven::as_factor(pub635), ordered = TRUE),
    pub650 = factor(haven::as_factor(pub650), ordered = TRUE),
    clon070 = as.numeric(zap_labels(clon070)),

    # Sleep auxiliary variables
    Weekday_sleep_duration_13y = as.numeric(zap_labels(Weekday_sleep_duration_13y)),
    Weekend_sleep_duration_13y = as.numeric(zap_labels(Weekend_sleep_duration_13y)),
    FJCI250  = as.numeric(zap_labels(FJCI250)),
    FKDQ4000 = haven::as_factor(FKDQ4000),
    FKDQ4100 = haven::as_factor(FKDQ4100),

    # Mental health auxiliary variables
    kv8603 = as.numeric(zap_labels(kv8603)),
    kv8613 = as.numeric(zap_labels(kv8613)),
    kv8617 = haven::as_factor(kv8617),
    kv8618 = haven::as_factor(kv8618),
    tb8603 = as.numeric(zap_labels(tb8603)),
    tb8619 = haven::as_factor(tb8619),
    tb8614 = as.numeric(zap_labels(tb8614)),
    tb8618 = haven::as_factor(tb8618),

    # Age auxiliary variables
    f7003c  = as.numeric(zap_labels(f7003c)),
    f9003c  = as.numeric(zap_labels(f9003c)),
    fd003c  = as.numeric(zap_labels(fd003c)),
    fe003c  = as.numeric(zap_labels(fe003c)),
    ff0011a = as.numeric(zap_labels(ff0011a)),
    YPB9992 = as.numeric(zap_labels(YPB9992)),
    kv9991a = as.numeric(zap_labels(kv9991a)),
    tb9991a = as.numeric(zap_labels(tb9991a)),
    ku991a  = as.numeric(zap_labels(ku991a))
  )

saveRDS(data_transformed, file = "data_transformed.rds")


# Section 5: Auxiliary variable selection ----
# Candidate auxiliary variables were evaluated for inclusion in the predictor
# matrix based on: (1) % missing lower than the variable to be imputed,
# (2) association with variables to be imputed, and (3) association with
# missingness indicators. The functions and variable vectors below document
# this selection process.

# Candidate auxiliary variable blocks
vars_aux_ACEs <- c("Townsend_Gest", "Townsend_7", "Townsend_9", "Townsend_10",
                   "Townsend_11", "Townsend_12", "Townsend_13", "Townsend_15",
                   "clon100", "clon101", "clon102", "clon103", "clon104",
                   "clon105", "clon106", "clon107", "clon108", "clon109",
                   "clon111", "clon112", "clon113", "clon114", "clon115",
                   "clon116", "clon117", "clon118", "clon119",
                   "e390", "b650", "b669", "clon166", "clon167", "clon168",
                   "clon169", "clon170", "qlet")
vars_aux_CVD <- c("bestgest", "CRP_f9", "CRP_TF4", "CRP_F24", "IL6_F9",
                  "f7ms026a", "f9ms026a", "fdms026a", "fems026a", "ff2039",
                  "fg3139", "fh3019", "f9ms018", "fdms018", "fems018", "ff2020",
                  "fg3120", "fh4020", "FKMS1052", "f9dx135", "fedx135", "fg3254",
                  "insulin_F9", "insulin_TF3", "HDL_f9", "hdl_TF3", "LDL_f9",
                  "glucose_TF3", "trig_f9", "trig_TF3", "FJAR040", "FJAR043",
                  "FJAR047", "FJAR048", "YPB1213", "YPB1214", "YPB1230",
                  "fdar117", "fdar118", "FJAR020c", "fdar111", "fdar113",
                  "fdar114", "fdar115", "fdar116", "fdar119", "FKBP1032")
vars_aux_MVPA    <- c("feag210", "fg1310", "fg1302", "fg1206", "fg1207",
                      "fg1308", "feag101", "feag102", "feag208")
vars_aux_dietary <- c("f7dd400", "fddd400", "fg1560", "f7dd100", "fddd100",
                      "f7dd507", "f7dd346", "f7dd600", "fddd506", "fddd345",
                      "fddd600", "fg1585")
vars_aux_puberty <- c("Child_sex", "pub435", "pub450", "pub455", "pub530",
                      "pub535", "pub550", "pub555", "pub630", "pub635",
                      "pub650", "clon070")
vars_aux_sleep   <- c("Weekday_sleep_duration_13y", "Weekend_sleep_duration_13y",
                      "FJCI250", "FKDQ4000", "FKDQ4100")
vars_aux_mental  <- c("kv8603", "kv8618", "kv8613", "kv8617", "tb8603",
                      "tb8619", "tb8614", "tb8618")
vars_aux_age     <- c("f7003c", "f9003c", "fd003c", "fe003c", "ff0011a",
                      "YPB9992", "kv9991a", "tb9991a", "ku991a")

# Variables to be imputed — by type
ACEs_general_num <- c("Classic_ACEs", "Extended_ACEs", "Mat_PND_gest")
ACEs_general_cat <- c("Mat_preg_smoke", "Mat_preg_alc", "Townsend_sum", "Marital_status")
CVD_outcomes_num <- c(vascular_outcomes_17, vascular_outcomes_24,
                      "BP_systolic_17", "BP_systolic_24", "BP_diastolic_17", "BP_diastolic_24",
                      "Birth_weight_kg", "Cortisol_15", "CRP_15", "IL6_24")
age_covs_num     <- c("Age_13_clinic_years", "Age_15_clinic_years",
                      "Age_17_clinic_years", "Age_24_clinic_years", "age17_c", "age24_c")
diet_covs_num    <- c("Diet_pattern_13", "Non_milk_sugar_13", "Days_diet_report_13")
mental_covs_cat  <- c("Child_ICD10_15", "Child_ICD10_17", "Child_GAD_15", "Child_GAD_24")
MVPA_covs_num    <- c("Daily_MVPA_15", "Daily_Light_PA_15", "Valid_Days_Wk_15", "Valid_Days_We_15")
sleep_covs_num   <- c("Weekday_sleep_duration_15y", "Weekend_sleep_duration_15y")

# Final set of candidate auxiliaries used in association tests
test_aux_all <- c(
  "Townsend_Gest", "Townsend_7", "Townsend_9", "Townsend_10", "Townsend_11",
  "Townsend_12", "Townsend_13", "Townsend_15",
  "clon100", "clon101", "clon102", "clon103", "clon104", "clon105", "clon106",
  "clon107", "clon108", "clon109", "clon111", "clon112", "clon113", "clon114",
  "clon115", "clon116", "clon117", "clon118", "clon119",
  "e390", "b650", "clon167", "clon168", "clon169", "clon170",
  "bestgest", "qlet", "CRP_f9", "IL6_F9", "f9ms018", "insulin_TF3",
  "HDL_f9", "glucose_TF3", "trig_f9", "FJAR040",
  "fdar117", "fdar114", "fdar116", "fdar118",
  "feag210", "feag101", "feag102", "feag208",
  "fddd100", "f7dd600", "fddd345", "fddd600",
  "pubic_hair_stage",
  "Weekday_sleep_duration_13y", "Weekend_sleep_duration_13y",
  "tb8603", "tb8614", "fd003c", "ff0011a"
)

# Combine puberty variables (male: pub455, female: pub435) into one column
data_transformed <- data_transformed %>%
  mutate(
    pubic_hair_stage = case_when(
      !is.na(pub455) & !is.na(pub435) ~
        ifelse(as.numeric(as.character(pub455)) >= as.numeric(as.character(pub435)),
               pub455, pub435),
      !is.na(pub455) ~ pub455,
      !is.na(pub435) ~ pub435,
      TRUE           ~ NA_character_
    )
  ) %>%
  mutate(pubic_hair_stage = factor(pubic_hair_stage, ordered = TRUE))

saveRDS(data_transformed, file = "data_transformed.rds")

# --- Association test functions ---

# Numeric vs numeric: Spearman correlation matrix
test_num_vs_num_fast <- function(num_df, outcome_df) {
  mat <- suppressWarnings(
    cor(num_df, outcome_df, method = "spearman", use = "pairwise.complete.obs")
  )
  expand.grid(Var1 = rownames(mat), Var2 = colnames(mat),
              stringsAsFactors = FALSE) |>
    dplyr::mutate(Correlation = as.vector(mat))
}

# Categorical vs numeric: t-test (binary) or ANOVA/Kruskal-Wallis (multi-level)
test_cat_vs_num <- function(cat_df, outcome_df) {
  results <- list()
  for (cat_var in names(cat_df)) {
    for (out_var in names(outcome_df)) {
      cat       <- cat_df[[cat_var]]
      out       <- outcome_df[[out_var]]
      valid_idx <- complete.cases(cat, out)
      cat       <- cat[valid_idx]
      out       <- out[valid_idx]
      n_groups  <- length(unique(cat))
      if (n_groups < 2) next

      p_val_ttest <- p_val_wilcox <- p_val_aov <- p_val_kruskal <- NA
      try({
        if (n_groups == 2) {
          p_val_ttest  <- t.test(out ~ cat)$p.value
          p_val_wilcox <- wilcox.test(out ~ cat, exact = FALSE)$p.value
        } else {
          aov_sum <- summary(aov(out ~ cat))[[1]]
          if (nrow(aov_sum) > 0) p_val_aov <- aov_sum[1, "Pr(>F)"]
          p_val_kruskal <- kruskal.test(out ~ cat)$p.value
        }
      }, silent = TRUE)

      test_type <- ifelse(n_groups == 2, "t-test / Wilcoxon", "ANOVA / Kruskal-Wallis")
      results[[length(results) + 1]] <- data.frame(
        Var1 = cat_var, Var2 = out_var, Test = test_type,
        TTest_P = p_val_ttest, Wilcox_P = p_val_wilcox,
        ANOVA_P = p_val_aov, Kruskal_P = p_val_kruskal
      )
    }
  }
  bind_rows(results)
}

# Categorical vs categorical: chi-squared (and Fisher for 2x2)
test_cat_vs_cat <- function(cat_df1, cat_df2) {
  results <- list()
  for (cat_var1 in names(cat_df1)) {
    for (cat_var2 in names(cat_df2)) {
      x         <- cat_df1[[cat_var1]]
      y         <- cat_df2[[cat_var2]]
      valid_idx <- complete.cases(x, y)
      x         <- droplevels(x[valid_idx])
      y         <- droplevels(y[valid_idx])
      tbl       <- table(x, y)
      keep_rows <- rowSums(tbl) > 0
      keep_cols <- colSums(tbl) > 0
      tbl       <- tbl[keep_rows, keep_cols, drop = FALSE]
      x         <- factor(x, levels = rownames(tbl))
      y         <- factor(y, levels = colnames(tbl))

      p_val_chisq <- p_val_fisher <- NA
      test_type   <- "Chi-squared"
      if (nrow(tbl) > 1 && ncol(tbl) > 1) {
        if (all(dim(tbl) == 2)) {
          fisher_res <- tryCatch(fisher.test(tbl), error = function(e) NULL)
          if (!is.null(fisher_res)) p_val_fisher <- fisher_res$p.value
          test_type <- "Chi-squared / Fisher"
        }
        chisq_res <- tryCatch(chisq.test(tbl), error = function(e) NULL)
        if (!is.null(chisq_res)) p_val_chisq <- chisq_res$p.value
      }
      results[[length(results) + 1]] <- data.frame(
        Var1 = cat_var1, Var2 = cat_var2, Test = test_type,
        Chisq_P = p_val_chisq, Fisher_P = p_val_fisher
      )
    }
  }
  dplyr::bind_rows(results)
}

# Faster screening versions for high-dimensional combinations
test_cat_vs_cat_fast <- function(cat_df1, cat_df2, max_levels = 10, do_fisher = TRUE) {
  results <- vector("list", length(names(cat_df1)) * length(names(cat_df2)))
  k <- 0
  for (cat_var1 in names(cat_df1)) {
    x0 <- if (!is.factor(cat_df1[[cat_var1]])) factor(cat_df1[[cat_var1]]) else cat_df1[[cat_var1]]
    if (nlevels(droplevels(x0)) < 2) next
    for (cat_var2 in names(cat_df2)) {
      y0 <- if (!is.factor(cat_df2[[cat_var2]])) factor(cat_df2[[cat_var2]]) else cat_df2[[cat_var2]]
      if (cat_var1 == cat_var2) next
      valid_idx <- stats::complete.cases(x0, y0)
      x <- droplevels(x0[valid_idx]); y <- droplevels(y0[valid_idx])
      if (length(x) < 10 || nlevels(x) < 2 || nlevels(y) < 2) next
      if (nlevels(x) > max_levels || nlevels(y) > max_levels) next
      tbl <- table(x, y)
      if (nrow(tbl) < 2 || ncol(tbl) < 2) next

      p_val_chisq <- p_val_fisher <- NA_real_
      test_type   <- "Chi-squared"
      if (do_fisher && all(dim(tbl) == 2)) {
        fisher_res <- tryCatch(fisher.test(tbl), error = function(e) NULL)
        if (!is.null(fisher_res)) p_val_fisher <- fisher_res$p.value
        test_type <- "Chi-squared / Fisher"
      }
      chisq_res <- tryCatch(chisq.test(tbl), error = function(e) NULL)
      if (!is.null(chisq_res)) p_val_chisq <- chisq_res$p.value

      k <- k + 1
      results[[k]] <- data.frame(Var1 = cat_var1, Var2 = cat_var2, Test = test_type,
                                 Chisq_P = p_val_chisq, Fisher_P = p_val_fisher,
                                 stringsAsFactors = FALSE)
    }
  }
  if (k == 0) return(dplyr::tibble(Var1 = character(), Var2 = character(),
                                   Test = character(), Chisq_P = numeric(), Fisher_P = numeric()))
  dplyr::bind_rows(results[seq_len(k)])
}

test_cat_vs_num_fast <- function(cat_df, outcome_df, min_n = 20, max_levels = 10) {
  results <- vector("list", length(names(cat_df)) * length(names(outcome_df)))
  k <- 0
  for (cat_var in names(cat_df)) {
    cat0 <- if (!is.factor(cat_df[[cat_var]])) factor(cat_df[[cat_var]]) else cat_df[[cat_var]]
    if (sum(!is.na(cat0)) < min_n) next
    if (nlevels(droplevels(cat0)) < 2 || nlevels(droplevels(cat0)) > max_levels) next
    for (out_var in names(outcome_df)) {
      out0 <- outcome_df[[out_var]]
      if (!is.numeric(out0)) next
      idx  <- stats::complete.cases(cat0, out0)
      cat  <- droplevels(cat0[idx]); out <- out0[idx]
      if (length(out) < min_n || nlevels(cat) < 2) next

      grand_mean  <- mean(out)
      group_means <- tapply(out, cat, mean)
      group_ns    <- table(cat)
      ss_between  <- sum(group_ns * (group_means - grand_mean)^2)
      ss_total    <- sum((out - grand_mean)^2)
      eta2        <- if (ss_total > 0) ss_between / ss_total else NA_real_

      k <- k + 1
      results[[k]] <- data.frame(Var1 = cat_var, Var2 = out_var,
                                 Test = "Eta2_screen", Stat = eta2,
                                 stringsAsFactors = FALSE)
    }
  }
  if (k == 0) return(dplyr::tibble(Var1 = character(), Var2 = character(),
                                   Test = character(), Stat = numeric()))
  dplyr::bind_rows(results[seq_len(k)])
}

# Run association tests between candidate auxiliaries and variables to be imputed
aux_ACEs_results_all <- bind_rows(
  test_num_vs_num_fast(
    data_transformed %>% select(all_of(ACEs_general_num)),
    data_transformed %>% select(all_of(test_aux_all)) %>% select(where(is.numeric))
  ),
  test_cat_vs_num(
    data_transformed %>% select(all_of(test_aux_all)) %>% select(where(is.factor)),
    data_transformed %>% select(all_of(ACEs_general_num))
  ),
  test_cat_vs_num(
    data_transformed %>% select(all_of(ACEs_general_cat)),
    data_transformed %>% select(all_of(test_aux_all)) %>% select(where(is.numeric))
  ),
  test_cat_vs_cat(
    data_transformed %>% select(all_of(ACEs_general_cat)),
    data_transformed %>% select(all_of(test_aux_all)) %>% select(where(is.factor))
  )
)

aux_age_results_all <- test_num_vs_num_fast(
  data_transformed %>% select(all_of(age_covs_num)),
  data_transformed %>% select(all_of(test_aux_all)) %>% select(where(is.numeric))
)

aux_CVD_results_all <- bind_rows(
  test_num_vs_num_fast(
    data_transformed %>% select(all_of(CVD_outcomes_num)),
    data_transformed %>% select(all_of(test_aux_all)) %>% select(where(is.numeric))
  ),
  test_cat_vs_num_fast(
    data_transformed %>% select("Family_CVD"),
    data_transformed %>% select(all_of(test_aux_all)) %>% select(where(is.numeric))
  ),
  test_cat_vs_num_fast(
    data_transformed %>% select(all_of(test_aux_all)) %>% select(where(is.factor)),
    data_transformed %>% select(all_of(CVD_outcomes_num))
  ),
  test_cat_vs_cat_fast(
    data_transformed %>% select("Family_CVD"),
    data_transformed %>% select(all_of(test_aux_all)) %>% select(where(is.factor))
  )
)

aux_diet_results_all <- bind_rows(
  test_num_vs_num_fast(
    data_transformed %>% select(all_of(diet_covs_num)),
    data_transformed %>% select(all_of(test_aux_all)) %>% select(where(is.numeric))
  ),
  test_cat_vs_num(
    data_transformed %>% select("Diet_report_plausible_13"),
    data_transformed %>% select(all_of(test_aux_all)) %>% select(where(is.numeric))
  ),
  test_cat_vs_num(
    data_transformed %>% select(all_of(test_aux_all)) %>% select(where(is.factor)),
    data_transformed %>% select(all_of(diet_covs_num))
  ),
  test_cat_vs_cat(
    data_transformed %>% select("Diet_report_plausible_13"),
    data_transformed %>% select(all_of(test_aux_all)) %>% select(where(is.factor))
  )
)

aux_mental_results_all <- bind_rows(
  test_cat_vs_num(
    data_transformed %>% select(all_of(mental_covs_cat)),
    data_transformed %>% select(all_of(test_aux_all)) %>% select(where(is.numeric))
  ),
  test_cat_vs_cat(
    data_transformed %>% select(all_of(mental_covs_cat)),
    data_transformed %>% select(all_of(test_aux_all)) %>% select(where(is.factor))
  )
)

aux_MVPA_results_all <- test_num_vs_num_fast(
  data_transformed %>% select(all_of(MVPA_covs_num)),
  data_transformed %>% select(all_of(test_aux_all)) %>% select(where(is.numeric))
)

aux_puberty_results_all <- bind_rows(
  test_num_vs_num_fast(
    data_transformed %>% select("Age_PHV"),
    data_transformed %>% select(all_of(test_aux_all)) %>% select(where(is.numeric))
  ),
  test_cat_vs_num(
    data_transformed %>% select(all_of(test_aux_all)) %>% select(where(is.factor)),
    data_transformed %>% select("Age_PHV")
  )
)

aux_sleep_results_all <- test_num_vs_num_fast(
  data_transformed %>% select(all_of(sleep_covs_num)),
  data_transformed %>% select(all_of(test_aux_all)) %>% select(where(is.numeric))
)

# --- Auxiliary–missingness association tests ---
# Logistic regression testing whether observed auxiliaries predict missingness
# in the primary analytic variables (informs the missing-data mechanism).

# Categorical predictor vs missingness indicator (robust SE)
test_logit_cat_vs_missing <- function(cat_df, missing_vec) {
  results <- list()
  y <- as.numeric(missing_vec == "missing")
  for (var in names(cat_df)) {
    x   <- as.factor(cat_df[[var]])
    if (nlevels(x) < 2) next
    dat <- data.frame(y = y, x = x)
    fit <- tryCatch(glm(y ~ x, data = dat, family = binomial), error = function(e) NULL)
    if (!is.null(fit)) {
      coefs        <- summary(fit)$coefficients
      robust_coefs <- coeftest(fit, vcov = vcovHC(fit, type = "HC0"))
      for (i in 2:nrow(coefs)) {
        results[[length(results) + 1]] <- data.frame(
          Variable  = var,
          Level     = gsub("^x", "", rownames(coefs)[i]),
          OR        = exp(coefs[i, 1]),
          CI_Lower  = exp(coefs[i, 1] - 1.96 * coefs[i, 2]),
          CI_Upper  = exp(coefs[i, 1] + 1.96 * coefs[i, 2]),
          P_Value   = coefs[i, 4],
          Robust_SE = robust_coefs[i, 2]
        )
      }
    }
  }
  dplyr::bind_rows(results)
}

# Numeric predictor vs missingness indicator (robust SE)
test_logit_num_vs_missing <- function(num_df, missing_vec) {
  results <- list()
  y <- as.numeric(missing_vec == "missing")
  for (var in names(num_df)) {
    x   <- as.numeric(num_df[[var]])
    if (all(is.na(x))) next
    dat <- data.frame(y = y, x = x)
    fit <- tryCatch(glm(y ~ x, data = dat, family = binomial), error = function(e) NULL)
    if (!is.null(fit)) {
      est          <- summary(fit)$coefficients
      robust_coefs <- coeftest(fit, vcov = vcovHC(fit, type = "HC0"))
      results[[length(results) + 1]] <- data.frame(
        Variable  = var,
        OR        = exp(est[2, 1]),
        CI_Lower  = exp(est[2, 1] - 1.96 * est[2, 2]),
        CI_Upper  = exp(est[2, 1] + 1.96 * est[2, 2]),
        P_Value   = est[2, 4],
        Robust_SE = robust_coefs[2, 2]
      )
    }
  }
  dplyr::bind_rows(results)
}

possible_aux_numeric <- data_transformed %>%
  select(all_of(test_aux_all)) %>% select(where(is.numeric))
possible_aux_factor  <- data_transformed %>%
  select(all_of(test_aux_all)) %>% select(where(is.factor))

aux_missingness_17 <- bind_rows(
  test_logit_num_vs_missing(possible_aux_numeric, data_transformed$y17_missing_classic),
  test_logit_cat_vs_missing(possible_aux_factor,  data_transformed$y17_missing_classic)
) %>% mutate(year = "17")

aux_missingness_24 <- bind_rows(
  test_logit_num_vs_missing(possible_aux_numeric, data_transformed$y24_missing_classic),
  test_logit_cat_vs_missing(possible_aux_factor,  data_transformed$y24_missing_classic)
) %>% mutate(year = "24")

aux_missingness_17_ext <- bind_rows(
  test_logit_num_vs_missing(possible_aux_numeric, data_transformed$y17_missing_extended),
  test_logit_cat_vs_missing(possible_aux_factor,  data_transformed$y17_missing_extended)
) %>% mutate(year = "17ext")

aux_missingness_24_ext <- bind_rows(
  test_logit_num_vs_missing(possible_aux_numeric, data_transformed$y24_missing_extended),
  test_logit_cat_vs_missing(possible_aux_factor,  data_transformed$y24_missing_extended)
) %>% mutate(year = "24ext")

rm(possible_aux_numeric, possible_aux_factor,
   aux_ACEs_results_all, aux_CVD_results_all, aux_age_results_all, aux_diet_results_all,
   aux_mental_results_all, aux_puberty_results_all, aux_MVPA_results_all, aux_sleep_results_all,
   aux_missingness_17, aux_missingness_17_ext, aux_missingness_24, aux_missingness_24_ext)


# Section 6: MICE data preparation ----

# Drop invalid levels and standardise factor coding
data_transformed <- data_transformed %>%
  mutate(
    Family_CVD      = if_else(Family_CVD      %in% c("No", "Yes"),                        as.character(Family_CVD),      NA_character_),
    Child_ethnicity = if_else(Child_ethnicity %in% c("White", "Non-white"),               as.character(Child_ethnicity), NA_character_),
    Child_ICD10_15  = if_else(Child_ICD10_15  %in% c("No Depression", "Yes Depression"),  as.character(Child_ICD10_15),  NA_character_),
    Child_ICD10_17  = if_else(Child_ICD10_17  %in% c("No Depression", "Yes Depression"),  as.character(Child_ICD10_17),  NA_character_),
    Child_GAD_15    = if_else(Child_GAD_15    %in% c("No", "Yes"),                        as.character(Child_GAD_15),    NA_character_),
    Child_GAD_24    = if_else(Child_GAD_24    %in% c("No", "Yes"),                        as.character(Child_GAD_24),    NA_character_),
    clon100 = if_else(clon100 %in% c("No", "Yes"), as.character(clon100), NA_character_),
    clon101 = if_else(clon101 %in% c("No", "Yes"), as.character(clon101), NA_character_),
    clon102 = if_else(clon102 %in% c("No", "Yes"), as.character(clon102), NA_character_),
    clon103 = if_else(clon103 %in% c("No", "Yes"), as.character(clon103), NA_character_),
    clon104 = if_else(clon104 %in% c("No", "Yes"), as.character(clon104), NA_character_),
    clon105 = if_else(clon105 %in% c("No", "Yes"), as.character(clon105), NA_character_),
    clon106 = if_else(clon106 %in% c("No", "Yes"), as.character(clon106), NA_character_),
    clon107 = if_else(clon107 %in% c("No", "Yes"), as.character(clon107), NA_character_),
    clon108 = if_else(clon108 %in% c("No", "Yes"), as.character(clon108), NA_character_),
    clon109 = if_else(clon109 %in% c("No", "Yes"), as.character(clon109), NA_character_),
    clon111 = if_else(clon111 %in% c("No", "Yes"), as.character(clon111), NA_character_),
    clon112 = if_else(clon112 %in% c("No", "Yes"), as.character(clon112), NA_character_),
    clon113 = if_else(clon113 %in% c("No", "Yes"), as.character(clon113), NA_character_),
    clon114 = if_else(clon114 %in% c("No", "Yes"), as.character(clon114), NA_character_),
    clon115 = if_else(clon115 %in% c("No", "Yes"), as.character(clon115), NA_character_),
    clon116 = if_else(clon116 %in% c("No", "Yes"), as.character(clon116), NA_character_),
    clon117 = if_else(clon117 %in% c("No", "Yes"), as.character(clon117), NA_character_),
    clon118 = if_else(clon118 %in% c("No", "Yes"), as.character(clon118), NA_character_),
    clon119 = if_else(clon119 %in% c("No", "Yes"), as.character(clon119), NA_character_),
    b650    = if_else(b650    %in% c("N", "Y"),     as.character(b650),    NA_character_),
    FJAR040 = if_else(FJAR040 %in% c("No", "Yes"), as.character(FJAR040), NA_character_),
    FJAR048 = if_else(FJAR048 %in% c("No", "Yes"), as.character(FJAR048), NA_character_),
    fddd600 = if_else(fddd600 %in% c("Under reporter", "Valid reporter", "Over reporter"),
                      as.character(fddd600), NA_character_),
    pub435  = if_else(pub435  %in% c("Stage 1", "Stage 2", "Stage 3", "Stage 4", "Stage 5"),
                      as.character(pub435),  NA_character_),
    pub455  = if_else(pub455  %in% c("Stage 1", "Stage 2", "Stage 3", "Stage 4", "Stage 5"),
                      as.character(pub455),  NA_character_),
    Child_sex = if_else(Child_sex %in% c("Male", "Female"), as.character(Child_sex), NA_character_),
    Diet_report_plausible_13 = if_else(
      Diet_report_plausible_13 %in% c("Under reporter", "Valid reporter", "Over reporter"),
      as.character(Diet_report_plausible_13), NA_character_)
  ) %>%
  mutate(
    across(
      c(Family_CVD, Child_ethnicity, Child_ICD10_15, Child_ICD10_17,
        Child_GAD_15, Child_GAD_24,
        clon100, clon101, clon102, clon103, clon105, clon106, clon107,
        clon108, clon109, clon112, clon113, clon114, clon115, clon117,
        clon118, clon119, b650, FJAR040, FJAR048, pub435, pub455, Child_sex),
      ~ factor(.x) %>% droplevels()
    ),
    fddd600 = factor(fddd600,
                     levels  = c("Under reporter", "Valid reporter", "Over reporter"),
                     ordered = TRUE),
    Diet_report_plausible_13 = factor(Diet_report_plausible_13,
                                      levels  = c("Under reporter", "Valid reporter", "Over reporter"),
                                      ordered = TRUE)
  )

# Sanity checks
is.ordered(data_transformed$fddd600)
is.ordered(data_transformed$Diet_report_plausible_13)
table(data_transformed$fddd600, useNA = "ifany")
table(data_transformed$Diet_report_plausible_13, useNA = "ifany")

# Convert puberty variable to ordered factor
data_transformed <- data_transformed %>%
  mutate(across(pubic_hair_stage,
                ~ factor(haven::as_factor(.x), ordered = TRUE),
                .names = "{col}"))

# Variable vectors for imputation subsets
all_vars_main <- c(
  "cidB4619",
  "Classic_ACEs", "Extended_ACEs", "Child_sex",
  "PWV_17", "Arterial_Dist_17", "cIMT_17",
  "BMI_17", "Trig_17", "HDL_17", "Glucose_17", "Insulin_17",
  "PWV_24", "cIMT_24", "BMI_24", "Trig_24", "HDL_24", "Glucose_24", "Insulin_24",
  "Daily_MVPA_15", "Daily_Light_PA_15", "Valid_Days_Wk_15", "Valid_Days_We_15",
  "Diet_pattern_13", "Non_milk_sugar_13", "Days_diet_report_13", "Diet_report_plausible_13",
  "Weekday_sleep_duration_15y", "Weekend_sleep_duration_15y",
  "Child_alc_15", "Child_smoke_15", "Child_alc_24", "Child_smoke_24",
  "Child_ICD10_15", "Child_ICD10_17", "Child_GAD_15", "Child_GAD_24",
  "Age_13_clinic_years", "Age_15_clinic_years", "Age_17_clinic_years", "Age_24_clinic_years",
  "age17_c", "age24_c",
  "BP_systolic_17", "BP_systolic_24", "BP_diastolic_17", "BP_diastolic_24",
  "Cortisol_15", "CRP_15", "IL6_24", "Age_PHV",
  "Family_CVD", "Parent_edu", "Child_ethnicity", "Marital_status",
  "Mat_PND_gest", "Mat_age_delivery", "Birth_weight_kg",
  "Mat_preg_smoke", "Mat_preg_alc",
  "Townsend_Gest", "Townsend_7", "Townsend_9", "Townsend_10",
  "Townsend_11", "Townsend_12", "Townsend_13", "Townsend_15"
)

final_aux_all <- c(
  "clon100", "clon101", "clon102", "clon103", "clon104", "clon105", "clon106",
  "clon107", "clon108", "clon109", "clon111", "clon112", "clon113", "clon114",
  "clon115", "clon116", "clon117", "clon118", "clon119",
  "e390", "b650", "clon167", "clon168", "clon169", "clon170",
  "bestgest", "qlet", "CRP_f9", "IL6_F9", "f9ms018", "insulin_TF3",
  "HDL_f9", "glucose_TF3", "trig_f9", "FJAR040",
  "fdar117", "fdar118", "fdar114", "fdar116",
  "feag210", "feag101", "feag102", "feag208",
  "fddd100", "f7dd600", "fddd345", "fddd600",
  "pubic_hair_stage",
  "Weekday_sleep_duration_13y", "Weekend_sleep_duration_13y",
  "tb8603", "tb8614", "fd003c", "ff0011a"
)

all_vars_main_classic <- setdiff(all_vars_main, "Extended_ACEs")
all_vars_main_ext     <- setdiff(all_vars_main, "Classic_ACEs")

# Filter to core ALSPAC sample (kz011b = "Yes") with Child_sex observed
# kz011b: ALSPAC core sample membership flag (publicly documented)
data_sub <- data_transformed %>%
  filter(!is.na(Child_sex)) %>%
  mutate(kz011b = haven::as_factor(kz011b)) %>%
  dplyr::filter(haven::as_factor(kz011b) == "Yes")

# Make Townsend components numeric for passive imputation of Townsend sum
townsend_cols <- c("Townsend_Gest", "Townsend_7", "Townsend_9", "Townsend_10",
                   "Townsend_11", "Townsend_12", "Townsend_13", "Townsend_15")
data_sub <- data_sub %>%
  mutate(across(all_of(townsend_cols), ~ as.numeric(as.character(.x))))

# Create exposure-specific subsets
Study1_Sub_Classic <- data_sub %>% select(all_of(c(all_vars_main_classic, final_aux_all)))
Study1_Sub_Ext     <- data_sub %>% select(all_of(c(all_vars_main_ext,     final_aux_all)))

# Numeric sex variable (Female = 0, Male = 1) for interaction terms
Study1_Sub_Classic$Child_sex_num <- dplyr::case_when(
  Study1_Sub_Classic$Child_sex == "Female" ~ 0,
  Study1_Sub_Classic$Child_sex == "Male"   ~ 1,
  TRUE ~ NA_real_
)
Study1_Sub_Ext$Child_sex_num <- dplyr::case_when(
  Study1_Sub_Ext$Child_sex == "Female" ~ 0,
  Study1_Sub_Ext$Child_sex == "Male"   ~ 1,
  TRUE ~ NA_real_
)

# Interaction terms (passively imputed)
Study1_Sub_Classic$ACE_sex       <- Study1_Sub_Classic$Classic_ACEs  * Study1_Sub_Classic$Child_sex_num
Study1_Sub_Classic$age17_ACE     <- Study1_Sub_Classic$age17_c       * Study1_Sub_Classic$Classic_ACEs
Study1_Sub_Classic$age24_ACE     <- Study1_Sub_Classic$age24_c       * Study1_Sub_Classic$Classic_ACEs
Study1_Sub_Classic$age17_ACE_sex <- Study1_Sub_Classic$age17_c * Study1_Sub_Classic$Classic_ACEs * Study1_Sub_Classic$Child_sex_num
Study1_Sub_Classic$age24_ACE_sex <- Study1_Sub_Classic$age24_c * Study1_Sub_Classic$Classic_ACEs * Study1_Sub_Classic$Child_sex_num

Study1_Sub_Ext$ACE_sex       <- Study1_Sub_Ext$Extended_ACEs * Study1_Sub_Ext$Child_sex_num
Study1_Sub_Ext$age17_ACE     <- Study1_Sub_Ext$age17_c       * Study1_Sub_Ext$Extended_ACEs
Study1_Sub_Ext$age24_ACE     <- Study1_Sub_Ext$age24_c       * Study1_Sub_Ext$Extended_ACEs
Study1_Sub_Ext$age17_ACE_sex <- Study1_Sub_Ext$age17_c * Study1_Sub_Ext$Extended_ACEs * Study1_Sub_Ext$Child_sex_num
Study1_Sub_Ext$age24_ACE_sex <- Study1_Sub_Ext$age24_c * Study1_Sub_Ext$Extended_ACEs * Study1_Sub_Ext$Child_sex_num

# Townsend_sum initialised as NA — passively imputed as sum of components
Study1_Sub_Classic$Townsend_sum <- NA
Study1_Sub_Ext$Townsend_sum     <- NA

# Convert individual ACE item variables to binary numeric (0/1) for imputation
clon_vars <- c("clon100", "clon101", "clon102", "clon103", "clon104",
               "clon105", "clon106", "clon107", "clon108", "clon109",
               "clon111", "clon112", "clon113", "clon114", "clon115",
               "clon116", "clon117", "clon118", "clon119")

recode_clon <- function(df, vars) {
  df %>% mutate(across(all_of(vars), ~ case_when(
    . %in% c("yes", "Yes", "Y", 1) ~ 1,
    . %in% c("no",  "No",  "N", 0) ~ 0,
    TRUE ~ NA_real_
  )))
}
Study1_Sub_Classic <- recode_clon(Study1_Sub_Classic, clon_vars)
Study1_Sub_Ext     <- recode_clon(Study1_Sub_Ext,     clon_vars)

# Remove any remaining SPSS labels
Study1_Sub_Classic <- labelled::unlabelled(Study1_Sub_Classic)
Study1_Sub_Ext     <- labelled::unlabelled(Study1_Sub_Ext)

# Verify no labelled columns remain
sum(sapply(Study1_Sub_Classic, haven::is.labelled))
sum(sapply(Study1_Sub_Ext,     haven::is.labelled))

# Log-transform skewed variables for better imputation convergence
Study1_Sub_Classic <- Study1_Sub_Classic %>%
  mutate(log_Insulin_17        = log(Insulin_17        + 1),
         log_Insulin_24        = log(Insulin_24        + 1),
         log_Trig_17           = log(Trig_17           + 1),
         log_Trig_24           = log(Trig_24           + 1),
         log_Glucose_17        = log(Glucose_17        + 1),
         log_Glucose_24        = log(Glucose_24        + 1),
         log_Non_milk_sugar_13 = log(Non_milk_sugar_13 + 1)) %>%
  select(-Insulin_17, -Insulin_24, -Trig_17, -Trig_24,
         -Glucose_17, -Glucose_24, -Non_milk_sugar_13)

Study1_Sub_Ext <- Study1_Sub_Ext %>%
  mutate(log_Insulin_17        = log(Insulin_17        + 1),
         log_Insulin_24        = log(Insulin_24        + 1),
         log_Trig_17           = log(Trig_17           + 1),
         log_Trig_24           = log(Trig_24           + 1),
         log_Glucose_17        = log(Glucose_17        + 1),
         log_Glucose_24        = log(Glucose_24        + 1),
         log_Non_milk_sugar_13 = log(Non_milk_sugar_13 + 1)) %>%
  select(-Insulin_17, -Insulin_24, -Trig_17, -Trig_24,
         -Glucose_17, -Glucose_24, -Non_milk_sugar_13)

saveRDS(Study1_Sub_Classic, file = "Study1_Sub_Classic.rds")
saveRDS(Study1_Sub_Ext,     file = "Study1_Sub_Ext.rds")


# Section 7: MICE predictor matrix and method setup ----
# The predictor matrix was initially derived from quickpred() (mincor = 0.3, including all
# main variables) and then manually amended based on the auxiliary association
# results from Section 5. The final matrix is read from a saved CSV.

# --- Classic ACEs ---

Study1_Imp_Classic_Setup <- mice(Study1_Sub_Classic, maxit = 0)
method_classic_final     <- Study1_Imp_Classic_Setup$method

# Ordinal variables: proportional-odds logistic regression
method_classic_final["Child_alc_15"] <- "polr"
method_classic_final["Child_alc_24"] <- "polr"
method_classic_final["fddd600"]      <- "polr"

# Passive imputation: Classic ACEs score = sum of individual items
method_classic_final["Classic_ACEs"] <-
  "~ I(clon100 + clon101 + clon102 + clon103 + clon104 +
       clon105 + clon106 + clon107 + clon108 + clon109)"

# Passive imputation: auxiliary ACE count
method_classic_final["clon170"] <- "~I(clon167 + clon168 + clon169)"

# Passive imputation: interaction terms
method_classic_final["ACE_sex"]       <- "~ I(Classic_ACEs * Child_sex_num)"
method_classic_final["age17_ACE"]     <- "~ I(age17_c * Classic_ACEs)"
method_classic_final["age24_ACE"]     <- "~ I(age24_c * Classic_ACEs)"
method_classic_final["age17_ACE_sex"] <- "~ I(age17_c * Classic_ACEs * Child_sex_num)"
method_classic_final["age24_ACE_sex"] <- "~ I(age24_c * Classic_ACEs * Child_sex_num)"

# Passive imputation: centred age variables
method_classic_final["age17_c"] <- "~ I(Age_17_clinic_years - 17)"
method_classic_final["age24_c"] <- "~ I(Age_24_clinic_years - 17)"

# Child_sex_num: deterministic from Child_sex — do not impute
method_classic_final["Child_sex_num"] <- ""

# Townsend components: pmm; Townsend_sum: passive sum
for (col in townsend_cols) method_classic_final[col] <- "pmm"
method_classic_final["Townsend_sum"] <- paste0(
  "~I(", paste(townsend_cols, collapse = " + "), ")"
)

# Read final predictor matrix (amended from quickpred output)
pred_classic_final <- as.matrix(read.csv("quickpred_classic_final_v7.csv", row.names = 1))

# Alignment checks
stopifnot(all(colnames(pred_classic_final) == rownames(pred_classic_final)))
stopifnot(ncol(pred_classic_final) == ncol(Study1_Sub_Classic))
cat("Dims match:",     ncol(Study1_Sub_Classic) == ncol(pred_classic_final), "\n")
cat("Colnames match:", all(colnames(Study1_Sub_Classic) == colnames(pred_classic_final)), "\n")
cat("Rownames match:", all(colnames(Study1_Sub_Classic) == rownames(pred_classic_final)), "\n")

# --- Extended ACEs ---

Study1_Imp_Ext_Setup <- mice(Study1_Sub_Ext, maxit = 0)
pred_ext_final       <- as.matrix(read.csv("quickpred_ext_final_v7.csv", row.names = 1))
method_ext_final     <- Study1_Imp_Ext_Setup$method

method_ext_final["Child_alc_15"] <- "polr"
method_ext_final["Child_alc_24"] <- "polr"
method_ext_final["fddd600"]      <- "polr"

method_ext_final["Extended_ACEs"] <-
  "~I(clon100 + clon101 + clon102 + clon103 + clon104 + clon105 + clon106 +
      clon107 + clon108 + clon109 + clon111 + clon112 + clon113 + clon114 +
      clon115 + clon116 + clon117 + clon118 + clon119)"

method_ext_final["clon170"]       <- "~I(clon167 + clon168 + clon169)"
method_ext_final["ACE_sex"]       <- "~ I(Extended_ACEs * Child_sex_num)"
method_ext_final["age17_ACE"]     <- "~ I(age17_c * Extended_ACEs)"
method_ext_final["age24_ACE"]     <- "~ I(age24_c * Extended_ACEs)"
method_ext_final["age17_ACE_sex"] <- "~ I(age17_c * Extended_ACEs * Child_sex_num)"
method_ext_final["age24_ACE_sex"] <- "~ I(age24_c * Extended_ACEs * Child_sex_num)"
method_ext_final["age17_c"]       <- "~ I(Age_17_clinic_years - 17)"
method_ext_final["age24_c"]       <- "~ I(Age_24_clinic_years - 17)"
method_ext_final["Child_sex_num"] <- ""
for (col in townsend_cols) method_ext_final[col] <- "pmm"
method_ext_final["Townsend_sum"]  <- paste0("~I(", paste(townsend_cols, collapse = " + "), ")")


# Section 8: Final imputation runs ----
# Two independent runs per exposure model (different seeds) to assess convergence
# and allow pooling across a combined m = 100 datasets.

# --- Classic ACEs: Run 1 (seed 1912) ---
Study1_Imp_Classic_Final_Run1 <- mice(
  data            = Study1_Sub_Classic,
  m               = 50,
  maxit           = 50,
  predictorMatrix = pred_classic_final,
  method          = method_classic_final,
  seed            = 1912
)
saveRDS(Study1_Imp_Classic_Final_Run1, file = "Study1_Imp_Classic_Final_Run1.rds")

# --- Classic ACEs: Run 2 (seed 2112) ---
Study1_Imp_Classic_Final_Run2 <- mice(
  data            = Study1_Sub_Classic,
  m               = 50,
  maxit           = 50,
  predictorMatrix = pred_classic_final,
  method          = method_classic_final,
  seed            = 2112
)
saveRDS(Study1_Imp_Classic_Final_Run2, file = "Study1_Imp_Classic_Final_Run2.rds")

# --- Extended ACEs: Run 1 (seed 2312) ---
Study1_Imp_Ext_Final_Run1 <- mice(
  data            = Study1_Sub_Ext,
  m               = 50,
  maxit           = 50,
  predictorMatrix = pred_ext_final,
  method          = method_ext_final,
  seed            = 2312
)
saveRDS(Study1_Imp_Ext_Final_Run1, file = "Study1_Imp_Ext_Final_Run1.rds")

# --- Extended ACEs: Run 2 (seed 3112) ---
Study1_Imp_Ext_Final_Run2 <- mice(
  data            = Study1_Sub_Ext,
  m               = 50,
  maxit           = 50,
  predictorMatrix = pred_ext_final,
  method          = method_ext_final,
  seed            = 3112
)
saveRDS(Study1_Imp_Ext_Final_Run2, file = "Study1_Imp_Ext_Final_Run2.rds")


# Section 9: Post-imputation transformations ----
# Applied identically across all completed datasets before pooling.
# Includes back-transformation of log variables, diet calibration,
# CMR score construction, and derived categorical/flag variables.

# Helper: standardised residuals (for CMR z-score construction)
std_resid <- function(formula, data) {
  as.numeric(scale(residuals(lm(formula, data = data))))
}

transform_completed_datasets <- function(completed_list,
                                         exposure_type = c("Classic", "Extended"),
                                         beta_under, beta_over) {
  exposure_type <- match.arg(exposure_type)
  lapply(completed_list, function(df) {

    # Birth weight: grams → kg
    df <- df %>% mutate(Birth_weight_kg = Birth_weight_kg / 1000)

    # Back-transform log variables
    df <- df %>%
      mutate(
        Insulin_17        = exp(log_Insulin_17),
        Insulin_24        = exp(log_Insulin_24),
        Trig_17           = exp(log_Trig_17),
        Trig_24           = exp(log_Trig_24),
        Glucose_17        = exp(log_Glucose_17),
        Glucose_24        = exp(log_Glucose_24),
        Non_milk_sugar_13 = exp(log_Non_milk_sugar_13)
      )

    # Diet calibration (adjust for under/over reporting)
    df <- df %>%
      mutate(
        Diet_pattern_13_calib = dplyr::case_when(
          Diet_report_plausible_13 == "Under reporter" ~ Diet_pattern_13 - beta_under,
          Diet_report_plausible_13 == "Over reporter"  ~ Diet_pattern_13 - beta_over,
          TRUE                                          ~ Diet_pattern_13
        )
      )

    # CMR score at 17 (sex-stratified z-scores, age- and ethnicity-adjusted)
    df <- df %>%
      group_by(Child_sex) %>%
      group_modify(~ {
        dat <- .x
        dat %>% mutate(
          Z_SBP_17    = std_resid(BP_systolic_17 ~ Age_17_clinic_years + Child_ethnicity, dat),
          Z_BMI_17    = std_resid(BMI_17         ~ Age_17_clinic_years + Child_ethnicity, dat),
          Z_INS_17    = std_resid(Insulin_17     ~ Age_17_clinic_years + Child_ethnicity, dat),
          Z_TRIG_17   = std_resid(Trig_17        ~ Age_17_clinic_years + Child_ethnicity, dat),
          Z_GLU_17    = std_resid(Glucose_17     ~ Age_17_clinic_years + Child_ethnicity, dat),
          Z_HDL_17    = std_resid(HDL_17         ~ Age_17_clinic_years + Child_ethnicity, dat),
          Z_HDLinv_17 = -1 * Z_HDL_17,
          CMR_17      = Z_SBP_17 + Z_BMI_17 + Z_INS_17 + Z_TRIG_17 + Z_GLU_17 + Z_HDLinv_17
        )
      }) %>% ungroup()

    # CMR score at 24
    df <- df %>%
      group_by(Child_sex) %>%
      group_modify(~ {
        dat <- .x
        dat %>% mutate(
          Z_SBP_24    = std_resid(BP_systolic_24 ~ Age_24_clinic_years + Child_ethnicity, dat),
          Z_BMI_24    = std_resid(BMI_24         ~ Age_24_clinic_years + Child_ethnicity, dat),
          Z_INS_24    = std_resid(Insulin_24     ~ Age_24_clinic_years + Child_ethnicity, dat),
          Z_TRIG_24   = std_resid(Trig_24        ~ Age_24_clinic_years + Child_ethnicity, dat),
          Z_GLU_24    = std_resid(Glucose_24     ~ Age_24_clinic_years + Child_ethnicity, dat),
          Z_HDL_24    = std_resid(HDL_24         ~ Age_24_clinic_years + Child_ethnicity, dat),
          Z_HDLinv_24 = -1 * Z_HDL_24,
          CMR_24      = Z_SBP_24 + Z_BMI_24 + Z_INS_24 + Z_TRIG_24 + Z_GLU_24 + Z_HDLinv_24
        )
      }) %>% ungroup()

    # ACE categorical variables
    if (exposure_type == "Classic") {
      df <- df %>% mutate(
        Classic_ACEs_cat = cut(Classic_ACEs,
                               breaks = c(-Inf, 0, 3, Inf),
                               labels = c("0", "1-3", "4+"),
                               include.lowest = TRUE),
        Classic_ACEs_binary = factor(ifelse(Classic_ACEs > 0, "Yes", "No"),
                                     levels = c("No", "Yes"))
      )
    } else {
      df <- df %>% mutate(
        Ext_ACEs_cat = cut(Extended_ACEs,
                           breaks = c(-Inf, 0, 3, 5, Inf),
                           labels = c("0", "1-3", "4-5", "6+"),
                           include.lowest = TRUE),
        Ext_ACEs_binary = factor(ifelse(Extended_ACEs > 0, "Yes", "No"),
                                 levels = c("No", "Yes"))
      )
    }

    # PA and diet inclusion flags
    df <- df %>%
      mutate(
        PA_include_15   = Valid_Days_Wk_15 >= 3 & Valid_Days_We_15 >= 1,
        Diet_include_13 = Days_diet_report_13 >= 3
      )

    # Age at PHV centred at sample mean
    df <- df %>% mutate(Age_PHV_c = Age_PHV - mean(Age_PHV, na.rm = TRUE))

    # Hypertension flags
    df <- df %>%
      mutate(
        HYP_17 = BP_systolic_17 >= 140 | BP_diastolic_17 >= 90,
        HYP_24 = BP_systolic_24 >= 140 | BP_diastolic_24 >= 90
      )

    return(df)
  })
}

# Helper: add numeric sex to completed datasets
add_sex_num <- function(df_list) {
  lapply(df_list, function(df) { df$Child_sex_num <- as.numeric(df$Child_sex) - 1; df })
}

# --- Classic ACEs ---
completed_dfs_run1_Classic <- add_sex_num(
  lapply(1:Study1_Imp_Classic_Final_Run1$m, function(i) complete(Study1_Imp_Classic_Final_Run1, i))
)
completed_dfs_run2_Classic <- add_sex_num(
  lapply(1:Study1_Imp_Classic_Final_Run2$m, function(i) complete(Study1_Imp_Classic_Final_Run2, i))
)

# Pool diet calibration coefficients across imputed datasets
temp_mids_classic       <- miceadds::datalist2mids(c(completed_dfs_run1_Classic, completed_dfs_run2_Classic))
Diet_calib_fits_classic <- with(temp_mids_classic, {
  Diet_rep <- factor(Diet_report_plausible_13,
                     levels = c("Valid reporter", "Under reporter", "Over reporter"),
                     ordered = TRUE)
  lm(Diet_pattern_13 ~ Diet_rep)
})
beta_vec_classic   <- setNames(pool(Diet_calib_fits_classic)$pooled$estimate,
                                pool(Diet_calib_fits_classic)$pooled$term)
beta_under_classic <- beta_vec_classic["Diet_repUnder reporter"]
beta_over_classic  <- beta_vec_classic["Diet_repOver reporter"]

completed_dfs_list_Classic <- transform_completed_datasets(
  c(completed_dfs_run1_Classic, completed_dfs_run2_Classic),
  exposure_type = "Classic",
  beta_under    = beta_under_classic,
  beta_over     = beta_over_classic
)

Study1_Imp_Classic_Final_Transformed <- miceadds::datalist2mids(completed_dfs_list_Classic)
saveRDS(Study1_Imp_Classic_Final_Transformed, file = "Study1_Imp_Classic_Final_Transformed.rds")

rm(completed_dfs_run1_Classic, completed_dfs_run2_Classic, temp_mids_classic,
   Study1_Imp_Classic_Final_Run1, Study1_Imp_Classic_Final_Run2,
   Diet_calib_fits_classic, beta_vec_classic, beta_under_classic, beta_over_classic)

# --- Extended ACEs ---
completed_dfs_run1_Ext <- add_sex_num(
  lapply(1:Study1_Imp_Ext_Final_Run1$m, function(i) complete(Study1_Imp_Ext_Final_Run1, i))
)
completed_dfs_run2_Ext <- add_sex_num(
  lapply(1:Study1_Imp_Ext_Final_Run2$m, function(i) complete(Study1_Imp_Ext_Final_Run2, i))
)

temp_mids_ext       <- miceadds::datalist2mids(c(completed_dfs_run1_Ext, completed_dfs_run2_Ext))
Diet_calib_fits_ext <- with(temp_mids_ext, {
  Diet_rep <- factor(Diet_report_plausible_13,
                     levels = c("Valid reporter", "Under reporter", "Over reporter"),
                     ordered = TRUE)
  lm(Diet_pattern_13 ~ Diet_rep)
})
beta_vec_ext   <- setNames(pool(Diet_calib_fits_ext)$pooled$estimate,
                            pool(Diet_calib_fits_ext)$pooled$term)
beta_under_ext <- beta_vec_ext["Diet_repUnder reporter"]
beta_over_ext  <- beta_vec_ext["Diet_repOver reporter"]

completed_dfs_list_Ext <- transform_completed_datasets(
  c(completed_dfs_run1_Ext, completed_dfs_run2_Ext),
  exposure_type = "Extended",
  beta_under    = beta_under_ext,
  beta_over     = beta_over_ext
)

Study1_Imp_Ext_Final_Transformed <- miceadds::datalist2mids(completed_dfs_list_Ext)
saveRDS(Study1_Imp_Ext_Final_Transformed, file = "Study1_Imp_Ext_Final_Transformed.rds")

rm(completed_dfs_run1_Ext, completed_dfs_run2_Ext, temp_mids_ext,
   Study1_Imp_Ext_Final_Run1, Study1_Imp_Ext_Final_Run2,
   Diet_calib_fits_ext, beta_vec_ext, beta_under_ext, beta_over_ext)


# Section 10: Imputation diagnostics ----
# Trace plots, box-and-whisker plots, and density plots to assess convergence
# and plausibility of imputed values.

Study1_Imp_Classic_Final_Run1        <- readRDS("Study1_Imp_Classic_Final_Run1.rds")
Study1_Imp_Classic_Final_Run2        <- readRDS("Study1_Imp_Classic_Final_Run2.rds")
Study1_Imp_Ext_Final_Run1            <- readRDS("Study1_Imp_Ext_Final_Run1.rds")
Study1_Imp_Ext_Final_Run2            <- readRDS("Study1_Imp_Ext_Final_Run2.rds")
Study1_Imp_Classic_Final_Transformed <- readRDS("Study1_Imp_Classic_Final_Transformed.rds")
Study1_Imp_Ext_Final_Transformed     <- readRDS("Study1_Imp_Ext_Final_Transformed.rds")

# Subset to first 10 imputations for plotting
subset_mids_10 <- function(mids_obj) {
  mids_obj$m   <- 10
  mids_obj$imp <- lapply(mids_obj$imp, function(x) x[, 1:min(10, ncol(x)), drop = FALSE])
  mids_obj
}

# --- Trace plots (convergence) ---
pdf("Study1_Imp_Classic_TracePlots_Run1.pdf", width = 10, height = 8)
print(plot(Study1_Imp_Classic_Final_Run1, layout = c(2, 2)))
dev.off()

pdf("Study1_Imp_Classic_TracePlots_Run2.pdf", width = 10, height = 8)
print(plot(Study1_Imp_Classic_Final_Run2, layout = c(2, 2)))
dev.off()

pdf("Study1_Imp_Ext_TracePlots_Run1.pdf", width = 10, height = 8)
print(plot(Study1_Imp_Ext_Final_Run1, layout = c(2, 2)))
dev.off()

pdf("Study1_Imp_Ext_TracePlots_Run2.pdf", width = 10, height = 8)
print(plot(Study1_Imp_Ext_Final_Run2, layout = c(2, 2)))
dev.off()

# --- Box-and-whisker plots (distribution across imputations) ---
Classic_Run1_10 <- subset_mids_10(Study1_Imp_Classic_Final_Run1)
Classic_Run2_10 <- subset_mids_10(Study1_Imp_Classic_Final_Run2)
Ext_Run1_10     <- subset_mids_10(Study1_Imp_Ext_Final_Run1)
Ext_Run2_10     <- subset_mids_10(Study1_Imp_Ext_Final_Run2)

pdf("Study1_Imp_Classic_BWPlots_Run1.pdf", width = 10, height = 8)
print(bwplot(Classic_Run1_10, layout = c(3, 2)))
dev.off()

pdf("Study1_Imp_Classic_BWPlots_Run2.pdf", width = 10, height = 8)
print(bwplot(Classic_Run2_10, layout = c(3, 2)))
dev.off()

pdf("Study1_Imp_Ext_BWPlots_Run1.pdf", width = 10, height = 8)
print(bwplot(Ext_Run1_10, layout = c(3, 2)))
dev.off()

pdf("Study1_Imp_Ext_BWPlots_Run2.pdf", width = 10, height = 8)
print(bwplot(Ext_Run2_10, layout = c(3, 2)))
dev.off()

# --- Density plots (key outcomes and exposures) ---
Key_vars_classic <- c("Classic_ACEs", "PWV_17", "cIMT_17", "Arterial_Dist_17",
                      "BMI_17", "HDL_17", "Trig_17", "Insulin_17",
                      "PWV_24", "cIMT_24", "BMI_24", "Trig_24",
                      "HDL_24", "Glucose_24", "Insulin_24")
key_vars_ext     <- c("Extended_ACEs", "PWV_17", "cIMT_17", "Arterial_Dist_17",
                      "BMI_17", "HDL_17", "Trig_17", "Insulin_17",
                      "PWV_24", "cIMT_24", "BMI_24", "Trig_24",
                      "HDL_24", "Glucose_24", "Insulin_24")

Classic_Trans_10 <- subset_mids_10(Study1_Imp_Classic_Final_Transformed)
Ext_Trans_10     <- subset_mids_10(Study1_Imp_Ext_Final_Transformed)

pdf("Study1_Imp_Classic_DensityPlots_Outcomes.pdf", width = 10, height = 8)
print(densityplot(Classic_Trans_10,
                  as.formula(paste("~", paste(Key_vars_classic, collapse = " + "))),
                  layout = c(3, 2)))
dev.off()

pdf("Study1_Imp_Ext_DensityPlots_Outcomes.pdf", width = 10, height = 8)
print(densityplot(Ext_Trans_10,
                  as.formula(paste("~", paste(key_vars_ext, collapse = " + "))),
                  layout = c(3, 2)))
dev.off()

# --- Density plots (key covariates) ---
key_covariates_num <- c("BP_systolic_17", "BP_systolic_24", "Cortisol_15",
                        "CRP_15", "IL6_24", "Age_PHV", "Townsend_sum",
                        "Daily_MVPA_15", "Daily_Light_PA_15",
                        "Weekday_sleep_duration_15y", "Weekend_sleep_duration_15y",
                        "Diet_pattern_13")

pdf("Study1_Imp_Classic_DensityPlots_Covariates.pdf", width = 10, height = 8)
print(densityplot(Classic_Trans_10,
                  as.formula(paste("~", paste(key_covariates_num, collapse = " + "))),
                  layout = c(3, 2)))
dev.off()

pdf("Study1_Imp_Ext_DensityPlots_Covariates.pdf", width = 10, height = 8)
print(densityplot(Ext_Trans_10,
                  as.formula(paste("~", paste(key_covariates_num, collapse = " + "))),
                  layout = c(3, 2)))
dev.off()

rm(Classic_Run1_10, Classic_Run2_10, Ext_Run1_10, Ext_Run2_10,
   Classic_Trans_10, Ext_Trans_10)


# Section 11: Monte Carlo error check ----
# Verifies that m = 100 imputations is sufficient by checking MCE relative to SE.

mc_error_check <- function(mids_obj, outcome, exposure, extra_covar, label) {
  model     <- with(data = mids_obj,
                    exp = lm(as.formula(paste(outcome, "~", exposure, "+", extra_covar))))
  m         <- mids_obj$m
  estimates <- sapply(model$analyses, coef)
  mce       <- apply(estimates, 1, sd) / sqrt(m)
  pooled    <- summary(pool(model))
  results   <- data.frame(term              = pooled$term,
                          estimate          = pooled$estimate,
                          std.error         = pooled$std.error,
                          monte_carlo_error = mce)
  cat("\nMCE check —", label, "\n")
  print(results)
  invisible(results)
}

mc_error_check(Study1_Imp_Classic_Final_Transformed,
               outcome = "PWV_17", exposure = "Classic_ACEs",
               extra_covar = "BP_systolic_17", label = "Classic ACEs")

mc_error_check(Study1_Imp_Ext_Final_Transformed,
               outcome = "PWV_17", exposure = "Extended_ACEs",
               extra_covar = "BP_systolic_17", label = "Extended ACEs")
