# File header ----
# Complete-Case OLS Regression Analyses
# Study: Adverse Childhood Experiences and Changes in Vascular Health from
#        Adolescence to Early Adulthood: Cross-Sectional and Longitudinal
#        Evidence from the ALSPAC Study
#
# Description:
#   Data preparation (complete-case subsetting) followed by ordinary
#   least squares (OLS) regression for Classic ACEs (continuous and
#   categorical) across three vascular outcomes:
#     - Pulse Wave Velocity (PWV) at 17y and 24y
#     - Carotid Intima-Media Thickness (cIMT) at 17y and 24y
#     - Arterial Distensibility at 17y
#   Models M1-M5 run for each outcome x stratum (whole sample, females,
#   males). Model assumption checks (avPlots, check_model, plot) are
#   included throughout. Results are exported as Word documents.
#
# Data:
#   ALSPAC (Avon Longitudinal Study of Parents and Children).
#   Data access: https://www.bristol.ac.uk/alspac/researchers/access/
#   Please cite: Boyd A et al. (2013) Cohort Profile: the 'Children of the
#   90s'—the index offspring of the Avon Longitudinal Study of Parents and
#   Children. Int J Epidemiol, 42(1):111-127.
#   https://doi.org/10.1093/ije/dys064
#
# Input:
#   `data_transformed`: the cleaned, renamed, and transformed ALSPAC dataset
#   produced by the preliminary analysis script. Must be loaded into the
#   R environment before running this script.
#
# Variable notes:
#   - Age_PHV_c: age at peak height velocity, grand-mean centred. Created
#                below via mutate().
#   - Child_sex reference level: Female (set explicitly below).
#   - Parent_edu reference level: "No qualifications" (set explicitly below).
#   - PA_include_15 / Diet_include_13: quality-control flags applied only
#     as a `subset =` restriction on the Model 3 and Model 5 lm() calls
#     (which include PA/diet as lifestyle mediators). Models 1, 2, and 4
#     use the full available sample.
#
# Author:  Laura Macro
# Date:    June 2026


# Required packages ----

library(car)          # avPlots()
library(performance)  # check_model()
library(ggplot2)      # assumption scatter plots
library(see)          # required by performance for check_model() plots
library(broom)        # tidy(), confint()
library(dplyr)
library(tidyr)
library(tibble)
library(flextable)
library(purrr)
library(officer)
library(parameters)


# DATA PREPARATION ----

# Save a copy of the transformed dataset for reproducibility records.
# Replace the path/filename as appropriate for your project directory.
saveRDS(data_transformed, file = "data_transformed.rds")

# PA/diet quality-control flags ----
# PA_include_15 and Diet_include_13 are quality-control flags in ALSPAC
# indicating whether the accelerometry and diet recall data meet minimum
# wear/completion thresholds. The working dataset below is NOT filtered
# on these flags up front: Models 1, 2, and 4 use the full available
# sample. Only Model 3 and Model 5 (which include PA/diet as lifestyle
# mediators) restrict to PA_include_15 == TRUE & Diet_include_13 == TRUE,
# applied via the `subset =` argument on those specific lm() calls below.

data_transformed_CC <- data_transformed

saveRDS(data_transformed_CC, file = "data_transformed_CC.rds")

# Centre age at peak height velocity ----
data_transformed_CC <- data_transformed_CC %>%
  mutate(Age_PHV_c = Age_PHV - mean(Age_PHV, na.rm = TRUE))

# Unorder ordered factors ----
# Some factor variables were imported as ordered; converting to unordered
# for standard dummy-coded treatment contrasts in lm().
vars_to_unorder <- c(
  "Parent_edu",
  "Mat_preg_smoke",
  "Mat_preg_alc",
  "Child_ethnicity",
  "Child_alc_15"
)

vars_to_unorder <- intersect(vars_to_unorder, names(data_transformed_CC))

data_transformed_CC[vars_to_unorder] <- lapply(
  data_transformed_CC[vars_to_unorder],
  function(x) factor(x, ordered = FALSE)
)

# Set reference categories ----

# Parent education: reference = "No qualifications"
if ("Parent_edu" %in% names(data_transformed_CC)) {
  data_transformed_CC$Parent_edu <- relevel(
    data_transformed_CC$Parent_edu,
    ref = "No qualifications"
  )
}

# Child ethnicity: reference = "White"
if ("Child_ethnicity" %in% names(data_transformed_CC)) {
  data_transformed_CC$Child_ethnicity <- relevel(
    data_transformed_CC$Child_ethnicity,
    ref = "White"
  )
}

# Child alcohol at 15y: set explicit level order (unordered factor)
if ("Child_alc_15" %in% names(data_transformed_CC)) {
  data_transformed_CC$Child_alc_15 <- factor(
    data_transformed_CC$Child_alc_15,
    levels = c("None", "1 or 2", "3 to 5", "6 to 9",
               "10 to 19", "20 to 39", "40+"),
    ordered = FALSE
  )
}

# Child sex: restrict to Male/Female, set Female as reference
if ("Child_sex" %in% names(data_transformed_CC)) {
  data_transformed_CC$Child_sex <- ifelse(
    data_transformed_CC$Child_sex %in% c("Male", "Female"),
    as.character(data_transformed_CC$Child_sex),
    NA_character_
  )
  data_transformed_CC$Child_sex <- factor(
    data_transformed_CC$Child_sex,
    levels  = c("Female", "Male"),
    ordered = FALSE
  )
}

# Define variable lists ----

vars_exposures_classic <- c("Classic_ACEs", "Classic_ACEs_Cat",
                             "Classic_ACEs_binary")
vars_exposures_ext     <- c("Extended_ACEs", "Ext_ACEs_cat",
                             "Ext_ACEs_binary")

vars_outcome_17 <- c("PWV_17", "cIMT_17", "Arterial_Dist_17")
vars_outcome_24 <- c("PWV_24", "cIMT_24")

vars_covariates <- c(
  "Age_17_clinic_years", "Age_24_clinic_years",
  "BP_systolic_17", "BP_systolic_24",
  "Child_ethnicity", "Marital_status",
  "Mat_PND_gest", "Mat_age_delivery", "Birth_weight_kg",
  "Mat_preg_smoke", "Mat_preg_alc",
  "Townsend_sum", "Age_PHV_c", "Family_CVD", "Parent_edu",
  "BMI_17", "BMI_24",
  "Daily_MVPA_15", "Daily_Light_PA_15",
  "Diet_pattern_13_calib", "Non_milk_sugar_13",
  "Weekday_sleep_duration_15y", "Weekend_sleep_duration_15y",
  "Child_alc_15", "Child_smoke_15",
  "Child_alc_24", "Child_smoke_24",
  "Cortisol_15", "CRP_15", "IL6_24"
)

# Create complete-case indicators ----
# Flags for whether each participant has complete data on exposures +
# outcomes (y*_complete_*) and exposures + outcomes + all covariates
# (y*_cov_complete_*).

data_transformed_CC <- data_transformed_CC %>%
  mutate(
    y17_complete_classic = factor(
      complete.cases(select(., all_of(c(vars_exposures_classic,
                                        vars_outcome_17)))),
      levels = c(FALSE, TRUE), labels = c("incomplete", "complete")
    ),
    y17_complete_extended = factor(
      complete.cases(select(., all_of(c(vars_exposures_ext,
                                        vars_outcome_17)))),
      levels = c(FALSE, TRUE), labels = c("incomplete", "complete")
    ),
    y24_complete_classic = factor(
      complete.cases(select(., all_of(c(vars_exposures_classic,
                                        vars_outcome_24)))),
      levels = c(FALSE, TRUE), labels = c("incomplete", "complete")
    ),
    y24_complete_extended = factor(
      complete.cases(select(., all_of(c(vars_exposures_ext,
                                        vars_outcome_24)))),
      levels = c(FALSE, TRUE), labels = c("incomplete", "complete")
    ),
    y17_cov_complete_classic = factor(
      complete.cases(select(., all_of(c(vars_exposures_classic,
                                        vars_outcome_17, vars_covariates)))),
      levels = c(FALSE, TRUE), labels = c("incomplete", "complete")
    ),
    y17_cov_complete_extended = factor(
      complete.cases(select(., all_of(c(vars_exposures_ext,
                                        vars_outcome_17, vars_covariates)))),
      levels = c(FALSE, TRUE), labels = c("incomplete", "complete")
    ),
    y24_cov_complete_classic = factor(
      complete.cases(select(., all_of(c(vars_exposures_classic,
                                        vars_outcome_24, vars_covariates)))),
      levels = c(FALSE, TRUE), labels = c("incomplete", "complete")
    ),
    y24_cov_complete_extended = factor(
      complete.cases(select(., all_of(c(vars_exposures_ext,
                                        vars_outcome_24, vars_covariates)))),
      levels = c(FALSE, TRUE), labels = c("incomplete", "complete")
    )
  )

# Create complete-case analysis subsets ----

cov_complete_classic_17 <- data_transformed_CC %>%
  filter(!!sym("y17_cov_complete_classic") == "complete")

cov_complete_ext_17 <- data_transformed_CC %>%
  filter(!!sym("y17_cov_complete_extended") == "complete")

cov_complete_classic_24 <- data_transformed_CC %>%
  filter(!!sym("y24_cov_complete_classic") == "complete")

cov_complete_ext_24 <- data_transformed_CC %>%
  filter(!!sym("y24_cov_complete_extended") == "complete")

# Participants complete at both 17y and 24y (used in mixed-effects script)
cov_complete_both_classic <- data_transformed_CC %>%
  filter(
    !!sym("y17_cov_complete_classic") == "complete",
    !!sym("y24_cov_complete_classic") == "complete"
  )


# CLASSIC ACEs — CONTINUOUS EXPOSURE ----

# Model 1: unadjusted (age + BP only) ----

# PWV ----

PWV_Model1_Classic_17_Whole <- lm(PWV_17 ~ Classic_ACEs + BP_systolic_17 +
                                    Age_17_clinic_years,
                                  cov_complete_classic_17)
summary(PWV_Model1_Classic_17_Whole)

PWV_Model1_Classic_24_Whole <- lm(PWV_24 ~ Classic_ACEs + BP_systolic_24 +
                                    Age_24_clinic_years,
                                  cov_complete_classic_24)
summary(PWV_Model1_Classic_24_Whole)

PWV_Model1_Classic_17_F <- lm(PWV_17 ~ Classic_ACEs + BP_systolic_17 +
                                 Age_17_clinic_years,
                               cov_complete_classic_17,
                               subset = (Child_sex == "Female"))
summary(PWV_Model1_Classic_17_F)

PWV_Model1_Classic_17_M <- lm(PWV_17 ~ Classic_ACEs + BP_systolic_17 +
                                 Age_17_clinic_years,
                               cov_complete_classic_17,
                               subset = (Child_sex == "Male"))
summary(PWV_Model1_Classic_17_M)

PWV_Model1_Classic_24_F <- lm(PWV_24 ~ Classic_ACEs + BP_systolic_24 +
                                 Age_24_clinic_years,
                               cov_complete_classic_24,
                               subset = (Child_sex == "Female"))
summary(PWV_Model1_Classic_24_F)

PWV_Model1_Classic_24_M <- lm(PWV_24 ~ Classic_ACEs + BP_systolic_24 +
                                 Age_24_clinic_years,
                               cov_complete_classic_24,
                               subset = (Child_sex == "Male"))
summary(PWV_Model1_Classic_24_M)

# Assumption checks (linearity, homoscedasticity, normality)
ggplot(cov_complete_classic_17,
       aes(x = Classic_ACEs, y = PWV_17)) +
  geom_point() + geom_smooth(method = "lm")
check_model(PWV_Model1_Classic_17_Whole)

ggplot(cov_complete_classic_24,
       aes(x = Classic_ACEs, y = PWV_24)) +
  geom_point() + geom_smooth(method = "lm")
check_model(PWV_Model1_Classic_24_Whole)

ggplot(subset(cov_complete_classic_17, Child_sex == "Female"),
       aes(x = Classic_ACEs, y = PWV_17)) +
  geom_point() + geom_smooth(method = "lm")
check_model(PWV_Model1_Classic_17_F)

ggplot(subset(cov_complete_classic_17, Child_sex == "Male"),
       aes(x = Classic_ACEs, y = PWV_17)) +
  geom_point() + geom_smooth(method = "lm")
check_model(PWV_Model1_Classic_17_M)

ggplot(subset(cov_complete_classic_24, Child_sex == "Female"),
       aes(x = Classic_ACEs, y = PWV_24)) +
  geom_point() + geom_smooth(method = "lm")
check_model(PWV_Model1_Classic_24_F)

ggplot(subset(cov_complete_classic_24, Child_sex == "Male"),
       aes(x = Classic_ACEs, y = PWV_24)) +
  geom_point() + geom_smooth(method = "lm")
check_model(PWV_Model1_Classic_24_M)

# cIMT ----

cIMT_Model1_Classic_17_Whole <- lm(cIMT_17 ~ Classic_ACEs +
                                      Age_17_clinic_years,
                                    cov_complete_classic_17)
summary(cIMT_Model1_Classic_17_Whole)

cIMT_Model1_Classic_24_Whole <- lm(cIMT_24 ~ Classic_ACEs +
                                      Age_24_clinic_years,
                                    cov_complete_classic_24)
summary(cIMT_Model1_Classic_24_Whole)

cIMT_Model1_Classic_17_F <- lm(cIMT_17 ~ Classic_ACEs + Age_17_clinic_years,
                                cov_complete_classic_17,
                                subset = (Child_sex == "Female"))
summary(cIMT_Model1_Classic_17_F)

cIMT_Model1_Classic_17_M <- lm(cIMT_17 ~ Classic_ACEs + Age_17_clinic_years,
                                cov_complete_classic_17,
                                subset = (Child_sex == "Male"))
summary(cIMT_Model1_Classic_17_M)

cIMT_Model1_Classic_24_F <- lm(cIMT_24 ~ Classic_ACEs + Age_24_clinic_years,
                                cov_complete_classic_24,
                                subset = (Child_sex == "Female"))
summary(cIMT_Model1_Classic_24_F)

cIMT_Model1_Classic_24_M <- lm(cIMT_24 ~ Classic_ACEs + Age_24_clinic_years,
                                cov_complete_classic_24,
                                subset = (Child_sex == "Male"))
summary(cIMT_Model1_Classic_24_M)

# Assumption checks
ggplot(cov_complete_classic_17,
       aes(x = Classic_ACEs, y = cIMT_17)) +
  geom_point() + geom_smooth(method = "lm")
check_model(cIMT_Model1_Classic_17_Whole)

ggplot(cov_complete_classic_24,
       aes(x = Classic_ACEs, y = cIMT_24)) +
  geom_point() + geom_smooth(method = "lm")
check_model(cIMT_Model1_Classic_24_Whole)

ggplot(subset(cov_complete_classic_17, Child_sex == "Female"),
       aes(x = Classic_ACEs, y = cIMT_17)) +
  geom_point() + geom_smooth(method = "lm")
check_model(cIMT_Model1_Classic_17_F)

ggplot(subset(cov_complete_classic_17, Child_sex == "Male"),
       aes(x = Classic_ACEs, y = cIMT_17)) +
  geom_point() + geom_smooth(method = "lm")
check_model(cIMT_Model1_Classic_17_M)

ggplot(subset(cov_complete_classic_24, Child_sex == "Female"),
       aes(x = Classic_ACEs, y = cIMT_24)) +
  geom_point() + geom_smooth(method = "lm")
check_model(cIMT_Model1_Classic_24_F)

ggplot(subset(cov_complete_classic_24, Child_sex == "Male"),
       aes(x = Classic_ACEs, y = cIMT_24)) +
  geom_point() + geom_smooth(method = "lm")
check_model(cIMT_Model1_Classic_24_M)

# Arterial Distensibility ----

Dist_Model1_Classic_17_Whole <- lm(Arterial_Dist_17 ~ Classic_ACEs +
                                     BP_systolic_17 + Age_17_clinic_years,
                                   cov_complete_classic_17)
summary(Dist_Model1_Classic_17_Whole)

Dist_Model1_Classic_17_F <- lm(Arterial_Dist_17 ~ Classic_ACEs +
                                  BP_systolic_17 + Age_17_clinic_years,
                                cov_complete_classic_17,
                                subset = (Child_sex == "Female"))
summary(Dist_Model1_Classic_17_F)

Dist_Model1_Classic_17_M <- lm(Arterial_Dist_17 ~ Classic_ACEs +
                                  BP_systolic_17 + Age_17_clinic_years,
                                cov_complete_classic_17,
                                subset = (Child_sex == "Male"))
summary(Dist_Model1_Classic_17_M)

# Assumption checks
ggplot(cov_complete_classic_17,
       aes(x = Classic_ACEs, y = Arterial_Dist_17)) +
  geom_point() + geom_smooth(method = "lm")
check_model(Dist_Model1_Classic_17_Whole)

ggplot(subset(cov_complete_classic_17, Child_sex == "Female"),
       aes(x = Classic_ACEs, y = Arterial_Dist_17)) +
  geom_point() + geom_smooth(method = "lm")
check_model(Dist_Model1_Classic_17_F)

ggplot(subset(cov_complete_classic_17, Child_sex == "Male"),
       aes(x = Classic_ACEs, y = Arterial_Dist_17)) +
  geom_point() + geom_smooth(method = "lm")
check_model(Dist_Model1_Classic_17_M)


# Model 2: sociodemographic covariates ----
# Note: ethnicity excluded due to insufficient level variation in the
# complete-case subset.

# PWV ----

PWV_Model2_Classic_17_Whole <- lm(PWV_17 ~ Classic_ACEs * Child_sex +
                                    BP_systolic_17 + Age_17_clinic_years +
                                    Townsend_sum + Marital_status + Parent_edu +
                                    Mat_PND_gest + Mat_age_delivery +
                                    Birth_weight_kg + Mat_preg_smoke +
                                    Mat_preg_alc + Family_CVD + Age_PHV_c +
                                    BMI_17,
                                  cov_complete_classic_17)
summary(PWV_Model2_Classic_17_Whole)

PWV_Model2_Classic_24_Whole <- lm(PWV_24 ~ Classic_ACEs * Child_sex +
                                    BP_systolic_24 + Age_24_clinic_years +
                                    Townsend_sum + Marital_status + Parent_edu +
                                    Mat_PND_gest + Mat_age_delivery +
                                    Birth_weight_kg + Mat_preg_smoke +
                                    Mat_preg_alc + Family_CVD + Age_PHV_c +
                                    BMI_24,
                                  cov_complete_classic_24)
summary(PWV_Model2_Classic_24_Whole)

PWV_Model2_Classic_17_F <- lm(PWV_17 ~ Classic_ACEs + BP_systolic_17 +
                                 Age_17_clinic_years + Townsend_sum +
                                 Marital_status + Parent_edu +
                                 Mat_PND_gest + Mat_age_delivery +
                                 Birth_weight_kg + Mat_preg_smoke +
                                 Mat_preg_alc + Family_CVD + Age_PHV_c +
                                 BMI_17,
                               cov_complete_classic_17,
                               subset = (Child_sex == "Female"))
summary(PWV_Model2_Classic_17_F)

PWV_Model2_Classic_17_M <- lm(PWV_17 ~ Classic_ACEs + BP_systolic_17 +
                                 Age_17_clinic_years + Townsend_sum +
                                 Marital_status + Parent_edu +
                                 Mat_PND_gest + Mat_age_delivery +
                                 Birth_weight_kg + Mat_preg_smoke +
                                 Mat_preg_alc + Family_CVD + Age_PHV_c +
                                 BMI_17,
                               cov_complete_classic_17,
                               subset = (Child_sex == "Male"))
summary(PWV_Model2_Classic_17_M)

PWV_Model2_Classic_24_F <- lm(PWV_24 ~ Classic_ACEs + BP_systolic_24 +
                                 Age_24_clinic_years + Townsend_sum +
                                 Marital_status + Parent_edu +
                                 Mat_PND_gest + Mat_age_delivery +
                                 Birth_weight_kg + Mat_preg_smoke +
                                 Mat_preg_alc + Family_CVD + Age_PHV_c +
                                 BMI_24,
                               cov_complete_classic_24,
                               subset = (Child_sex == "Female"))
summary(PWV_Model2_Classic_24_F)

PWV_Model2_Classic_24_M <- lm(PWV_24 ~ Classic_ACEs + BP_systolic_24 +
                                 Age_24_clinic_years + Townsend_sum +
                                 Marital_status + Parent_edu +
                                 Mat_PND_gest + Mat_age_delivery +
                                 Birth_weight_kg + Mat_preg_smoke +
                                 Mat_preg_alc + Family_CVD + Age_PHV_c +
                                 BMI_24,
                               cov_complete_classic_24,
                               subset = (Child_sex == "Male"))
summary(PWV_Model2_Classic_24_M)

# Assumption checks
avPlots(PWV_Model2_Classic_17_Whole); plot(PWV_Model2_Classic_17_Whole)
avPlots(PWV_Model2_Classic_24_Whole); plot(PWV_Model2_Classic_24_Whole)
avPlots(PWV_Model2_Classic_17_F);     plot(PWV_Model2_Classic_17_F)
avPlots(PWV_Model2_Classic_17_M);     plot(PWV_Model2_Classic_17_M)
avPlots(PWV_Model2_Classic_24_F);     plot(PWV_Model2_Classic_24_F)
avPlots(PWV_Model2_Classic_24_M);     plot(PWV_Model2_Classic_24_M)

# cIMT ----

cIMT_Model2_Classic_17_Whole <- lm(cIMT_17 ~ Classic_ACEs * Child_sex +
                                     Age_17_clinic_years + Townsend_sum +
                                     Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_17,
                                   cov_complete_classic_17)
summary(cIMT_Model2_Classic_17_Whole)

cIMT_Model2_Classic_24_Whole <- lm(cIMT_24 ~ Classic_ACEs * Child_sex +
                                     Age_24_clinic_years + Townsend_sum +
                                     Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_24,
                                   cov_complete_classic_24)
summary(cIMT_Model2_Classic_24_Whole)

cIMT_Model2_Classic_17_F <- lm(cIMT_17 ~ Classic_ACEs + Age_17_clinic_years +
                                  Townsend_sum + Marital_status + Parent_edu +
                                  Mat_PND_gest + Mat_age_delivery +
                                  Birth_weight_kg + Mat_preg_smoke +
                                  Mat_preg_alc + Family_CVD + Age_PHV_c +
                                  BMI_17,
                                cov_complete_classic_17,
                                subset = (Child_sex == "Female"))
summary(cIMT_Model2_Classic_17_F)

cIMT_Model2_Classic_17_M <- lm(cIMT_17 ~ Classic_ACEs + Age_17_clinic_years +
                                  Townsend_sum + Marital_status + Parent_edu +
                                  Mat_PND_gest + Mat_age_delivery +
                                  Birth_weight_kg + Mat_preg_smoke +
                                  Mat_preg_alc + Family_CVD + Age_PHV_c +
                                  BMI_17,
                                cov_complete_classic_17,
                                subset = (Child_sex == "Male"))
summary(cIMT_Model2_Classic_17_M)

cIMT_Model2_Classic_24_F <- lm(cIMT_24 ~ Classic_ACEs + Age_24_clinic_years +
                                  Townsend_sum + Marital_status + Parent_edu +
                                  Mat_PND_gest + Mat_age_delivery +
                                  Birth_weight_kg + Mat_preg_smoke +
                                  Mat_preg_alc + Family_CVD + Age_PHV_c +
                                  BMI_24,
                                cov_complete_classic_24,
                                subset = (Child_sex == "Female"))
summary(cIMT_Model2_Classic_24_F)

cIMT_Model2_Classic_24_M <- lm(cIMT_24 ~ Classic_ACEs + Age_24_clinic_years +
                                  Townsend_sum + Marital_status + Parent_edu +
                                  Mat_PND_gest + Mat_age_delivery +
                                  Birth_weight_kg + Mat_preg_smoke +
                                  Mat_preg_alc + Family_CVD + Age_PHV_c +
                                  BMI_24,
                                cov_complete_classic_24,
                                subset = (Child_sex == "Male"))
summary(cIMT_Model2_Classic_24_M)

# Assumption checks
avPlots(cIMT_Model2_Classic_17_Whole); plot(cIMT_Model2_Classic_17_Whole)
avPlots(cIMT_Model2_Classic_24_Whole); plot(cIMT_Model2_Classic_24_Whole)
avPlots(cIMT_Model2_Classic_17_F);     plot(cIMT_Model2_Classic_17_F)
avPlots(cIMT_Model2_Classic_17_M);     plot(cIMT_Model2_Classic_17_M)
avPlots(cIMT_Model2_Classic_24_F);     plot(cIMT_Model2_Classic_24_F)
avPlots(cIMT_Model2_Classic_24_M);     plot(cIMT_Model2_Classic_24_M)

# Arterial Distensibility ----

Dist_Model2_Classic_17_Whole <- lm(Arterial_Dist_17 ~ Classic_ACEs * Child_sex +
                                     BP_systolic_17 + Age_17_clinic_years +
                                     Townsend_sum + Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_17,
                                   cov_complete_classic_17)
summary(Dist_Model2_Classic_17_Whole)

Dist_Model2_Classic_17_F <- lm(Arterial_Dist_17 ~ Classic_ACEs +
                                  BP_systolic_17 + Age_17_clinic_years +
                                  Townsend_sum + Marital_status + Parent_edu +
                                  Mat_PND_gest + Mat_age_delivery +
                                  Birth_weight_kg + Mat_preg_smoke +
                                  Mat_preg_alc + Family_CVD + Age_PHV_c +
                                  BMI_17,
                                cov_complete_classic_17,
                                subset = (Child_sex == "Female"))
summary(Dist_Model2_Classic_17_F)

Dist_Model2_Classic_17_M <- lm(Arterial_Dist_17 ~ Classic_ACEs +
                                  BP_systolic_17 + Age_17_clinic_years +
                                  Townsend_sum + Marital_status + Parent_edu +
                                  Mat_PND_gest + Mat_age_delivery +
                                  Birth_weight_kg + Mat_preg_smoke +
                                  Mat_preg_alc + Family_CVD + Age_PHV_c +
                                  BMI_17,
                                cov_complete_classic_17,
                                subset = (Child_sex == "Male"))
summary(Dist_Model2_Classic_17_M)

# Assumption checks
avPlots(Dist_Model2_Classic_17_Whole); plot(Dist_Model2_Classic_17_Whole)
avPlots(Dist_Model2_Classic_17_F);     plot(Dist_Model2_Classic_17_F)
avPlots(Dist_Model2_Classic_17_M);     plot(Dist_Model2_Classic_17_M)


# Model 3: + lifestyle mediators (PA, diet, sleep, alcohol) ----

# PWV ----
# Note: child smoking removed due to insufficient level variation in the
# complete-case subset.

PWV_Model3_Classic_17_Whole <- lm(PWV_17 ~ Classic_ACEs + BP_systolic_17 +
                                    Age_17_clinic_years + Townsend_sum +
                                    Marital_status + Parent_edu +
                                    Mat_PND_gest + Mat_age_delivery +
                                    Birth_weight_kg + Mat_preg_smoke +
                                    Mat_preg_alc + Family_CVD + Age_PHV_c +
                                    BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                    Diet_pattern_13_calib + Non_milk_sugar_13 +
                                    Weekday_sleep_duration_15y +
                                    Weekend_sleep_duration_15y + Child_alc_15,
                                  cov_complete_classic_17, subset = PA_include_15 == TRUE & Diet_include_13 == TRUE)
summary(PWV_Model3_Classic_17_Whole)

PWV_Model3_Classic_24_Whole <- lm(PWV_24 ~ Classic_ACEs + BP_systolic_24 +
                                    Age_24_clinic_years + Townsend_sum +
                                    Marital_status + Parent_edu +
                                    Mat_PND_gest + Mat_age_delivery +
                                    Birth_weight_kg + Mat_preg_smoke +
                                    Mat_preg_alc + Family_CVD + Age_PHV_c +
                                    BMI_24 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                    Diet_pattern_13_calib + Non_milk_sugar_13 +
                                    Weekday_sleep_duration_15y +
                                    Weekend_sleep_duration_15y + Child_alc_15,
                                  cov_complete_classic_24, subset = PA_include_15 == TRUE & Diet_include_13 == TRUE)
summary(PWV_Model3_Classic_24_Whole)

PWV_Model3_Classic_17_F <- lm(PWV_17 ~ Classic_ACEs + BP_systolic_17 +
                                 Age_17_clinic_years + Townsend_sum +
                                 Marital_status + Parent_edu +
                                 Mat_PND_gest + Mat_age_delivery +
                                 Birth_weight_kg + Mat_preg_smoke +
                                 Mat_preg_alc + Family_CVD + Age_PHV_c +
                                 BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                 Diet_pattern_13_calib + Non_milk_sugar_13 +
                                 Weekday_sleep_duration_15y +
                                 Weekend_sleep_duration_15y + Child_alc_15,
                               cov_complete_classic_17,
                               subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(PWV_Model3_Classic_17_F)

PWV_Model3_Classic_17_M <- lm(PWV_17 ~ Classic_ACEs + BP_systolic_17 +
                                 Age_17_clinic_years + Townsend_sum +
                                 Marital_status + Parent_edu +
                                 Mat_PND_gest + Mat_age_delivery +
                                 Birth_weight_kg + Mat_preg_smoke +
                                 Mat_preg_alc + Family_CVD + Age_PHV_c +
                                 BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                 Diet_pattern_13_calib + Non_milk_sugar_13 +
                                 Weekday_sleep_duration_15y +
                                 Weekend_sleep_duration_15y + Child_alc_15,
                               cov_complete_classic_17,
                               subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(PWV_Model3_Classic_17_M)

PWV_Model3_Classic_24_F <- lm(PWV_24 ~ Classic_ACEs + BP_systolic_24 +
                                 Age_24_clinic_years + Townsend_sum +
                                 Marital_status + Parent_edu +
                                 Mat_PND_gest + Mat_age_delivery +
                                 Birth_weight_kg + Mat_preg_smoke +
                                 Mat_preg_alc + Family_CVD + Age_PHV_c +
                                 BMI_24 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                 Diet_pattern_13_calib + Non_milk_sugar_13 +
                                 Weekday_sleep_duration_15y +
                                 Weekend_sleep_duration_15y + Child_alc_15,
                               cov_complete_classic_24,
                               subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(PWV_Model3_Classic_24_F)

PWV_Model3_Classic_24_M <- lm(PWV_24 ~ Classic_ACEs + BP_systolic_24 +
                                 Age_24_clinic_years + Townsend_sum +
                                 Marital_status + Parent_edu +
                                 Mat_PND_gest + Mat_age_delivery +
                                 Birth_weight_kg + Mat_preg_smoke +
                                 Mat_preg_alc + Family_CVD + Age_PHV_c +
                                 BMI_24 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                 Diet_pattern_13_calib + Non_milk_sugar_13 +
                                 Weekday_sleep_duration_15y +
                                 Weekend_sleep_duration_15y + Child_alc_15,
                               cov_complete_classic_24,
                               subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(PWV_Model3_Classic_24_M)

# Assumption checks
avPlots(PWV_Model3_Classic_17_F); check_model(PWV_Model3_Classic_17_F)
avPlots(PWV_Model3_Classic_17_M); check_model(PWV_Model3_Classic_17_M)
avPlots(PWV_Model3_Classic_24_F); check_model(PWV_Model3_Classic_24_F)
avPlots(PWV_Model3_Classic_24_M); check_model(PWV_Model3_Classic_24_M)

# cIMT ----

cIMT_Model3_Classic_17_Whole <- lm(cIMT_17 ~ Classic_ACEs +
                                     Age_17_clinic_years + Townsend_sum +
                                     Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                     Diet_pattern_13_calib + Non_milk_sugar_13 +
                                     Weekday_sleep_duration_15y +
                                     Weekend_sleep_duration_15y + Child_alc_15,
                                   data = cov_complete_classic_17, subset = PA_include_15 == TRUE & Diet_include_13 == TRUE)
summary(cIMT_Model3_Classic_17_Whole)

cIMT_Model3_Classic_24_Whole <- lm(cIMT_24 ~ Classic_ACEs +
                                     Age_24_clinic_years + Townsend_sum +
                                     Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_24 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                     Diet_pattern_13_calib + Non_milk_sugar_13 +
                                     Weekday_sleep_duration_15y +
                                     Weekend_sleep_duration_15y + Child_alc_15,
                                   data = cov_complete_classic_24, subset = PA_include_15 == TRUE & Diet_include_13 == TRUE)
summary(cIMT_Model3_Classic_24_Whole)

cIMT_Model3_Classic_17_F <- lm(cIMT_17 ~ Classic_ACEs + Age_17_clinic_years +
                                  Townsend_sum + Marital_status + Parent_edu +
                                  Mat_PND_gest + Mat_age_delivery +
                                  Birth_weight_kg + Mat_preg_smoke +
                                  Mat_preg_alc + Family_CVD + Age_PHV_c +
                                  BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                  Diet_pattern_13_calib + Non_milk_sugar_13 +
                                  Weekday_sleep_duration_15y +
                                  Weekend_sleep_duration_15y + Child_alc_15,
                                data = cov_complete_classic_17,
                                subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(cIMT_Model3_Classic_17_F)

# Note: Marital_status removed for male subset due to insufficient variation
cIMT_Model3_Classic_17_M <- lm(cIMT_17 ~ Classic_ACEs + Age_17_clinic_years +
                                  Townsend_sum + Parent_edu +
                                  Mat_PND_gest + Mat_age_delivery +
                                  Birth_weight_kg + Mat_preg_smoke +
                                  Mat_preg_alc + Family_CVD + Age_PHV_c +
                                  BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                  Diet_pattern_13_calib + Non_milk_sugar_13 +
                                  Weekday_sleep_duration_15y +
                                  Weekend_sleep_duration_15y + Child_alc_15,
                                data = cov_complete_classic_17,
                                subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(cIMT_Model3_Classic_17_M)

cIMT_Model3_Classic_24_F <- lm(cIMT_24 ~ Classic_ACEs + Age_24_clinic_years +
                                  Townsend_sum + Marital_status + Parent_edu +
                                  Mat_PND_gest + Mat_age_delivery +
                                  Birth_weight_kg + Mat_preg_smoke +
                                  Mat_preg_alc + Family_CVD + Age_PHV_c +
                                  BMI_24 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                  Diet_pattern_13_calib + Non_milk_sugar_13 +
                                  Weekday_sleep_duration_15y +
                                  Weekend_sleep_duration_15y + Child_alc_15,
                                data = cov_complete_classic_24,
                                subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(cIMT_Model3_Classic_24_F)

cIMT_Model3_Classic_24_M <- lm(cIMT_24 ~ Classic_ACEs + Age_24_clinic_years +
                                  Townsend_sum + Marital_status + Parent_edu +
                                  Mat_PND_gest + Mat_age_delivery +
                                  Birth_weight_kg + Mat_preg_smoke +
                                  Mat_preg_alc + Family_CVD + Age_PHV_c +
                                  BMI_24 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                  Diet_pattern_13_calib + Non_milk_sugar_13 +
                                  Weekday_sleep_duration_15y +
                                  Weekend_sleep_duration_15y + Child_alc_15,
                                data = cov_complete_classic_24,
                                subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(cIMT_Model3_Classic_24_M)

# Assumption checks
avPlots(cIMT_Model3_Classic_17_F); plot(cIMT_Model3_Classic_17_F)
avPlots(cIMT_Model3_Classic_17_M); plot(cIMT_Model3_Classic_17_M)
avPlots(cIMT_Model3_Classic_24_F); plot(cIMT_Model3_Classic_24_F)
avPlots(cIMT_Model3_Classic_24_M); plot(cIMT_Model3_Classic_24_M)

# Arterial Distensibility ----

Dist_Model3_Classic_17_Whole <- lm(Arterial_Dist_17 ~ Classic_ACEs +
                                     BP_systolic_17 + Age_17_clinic_years +
                                     Townsend_sum + Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                     Diet_pattern_13_calib + Non_milk_sugar_13 +
                                     Weekday_sleep_duration_15y +
                                     Weekend_sleep_duration_15y + Child_alc_15,
                                   cov_complete_classic_17, subset = PA_include_15 == TRUE & Diet_include_13 == TRUE)
summary(Dist_Model3_Classic_17_Whole)

Dist_Model3_Classic_17_F <- lm(Arterial_Dist_17 ~ Classic_ACEs +
                                  BP_systolic_17 + Age_17_clinic_years +
                                  Townsend_sum + Marital_status + Parent_edu +
                                  Mat_PND_gest + Mat_age_delivery +
                                  Birth_weight_kg + Mat_preg_smoke +
                                  Mat_preg_alc + Family_CVD + Age_PHV_c +
                                  BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                  Diet_pattern_13_calib + Non_milk_sugar_13 +
                                  Weekday_sleep_duration_15y +
                                  Weekend_sleep_duration_15y + Child_alc_15,
                                cov_complete_classic_17,
                                subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(Dist_Model3_Classic_17_F)

Dist_Model3_Classic_17_M <- lm(Arterial_Dist_17 ~ Classic_ACEs +
                                  BP_systolic_17 + Age_17_clinic_years +
                                  Townsend_sum + Marital_status + Parent_edu +
                                  Mat_PND_gest + Mat_age_delivery +
                                  Birth_weight_kg + Mat_preg_smoke +
                                  Mat_preg_alc + Family_CVD + Age_PHV_c +
                                  BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                  Diet_pattern_13_calib + Non_milk_sugar_13 +
                                  Weekday_sleep_duration_15y +
                                  Weekend_sleep_duration_15y + Child_alc_15,
                                cov_complete_classic_17,
                                subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(Dist_Model3_Classic_17_M)

# Assumption checks
avPlots(Dist_Model3_Classic_17_F); plot(Dist_Model3_Classic_17_F)
avPlots(Dist_Model3_Classic_17_M); plot(Dist_Model3_Classic_17_M)


# Model 4: + biological mediators (cortisol, CRP, IL-6) ----

# PWV ----

PWV_Model4_Classic_17_Whole <- lm(PWV_17 ~ Classic_ACEs + BP_systolic_17 +
                                    Age_17_clinic_years + Townsend_sum +
                                    Marital_status + Parent_edu +
                                    Mat_PND_gest + Mat_age_delivery +
                                    Birth_weight_kg + Mat_preg_smoke +
                                    Mat_preg_alc + Family_CVD + Age_PHV_c +
                                    BMI_17 + Cortisol_15 + CRP_15 + IL6_24,
                                  cov_complete_classic_17)
summary(PWV_Model4_Classic_17_Whole)

PWV_Model4_Classic_24_Whole <- lm(PWV_24 ~ Classic_ACEs + BP_systolic_24 +
                                    Age_24_clinic_years + Townsend_sum +
                                    Marital_status + Parent_edu +
                                    Mat_PND_gest + Mat_age_delivery +
                                    Birth_weight_kg + Mat_preg_smoke +
                                    Mat_preg_alc + Family_CVD + Age_PHV_c +
                                    BMI_24 + Cortisol_15 + CRP_15 + IL6_24,
                                  cov_complete_classic_24)
summary(PWV_Model4_Classic_24_Whole)

PWV_Model4_Classic_17_F <- lm(PWV_17 ~ Classic_ACEs + BP_systolic_17 +
                                 Age_17_clinic_years + Townsend_sum +
                                 Marital_status + Parent_edu +
                                 Mat_PND_gest + Mat_age_delivery +
                                 Birth_weight_kg + Mat_preg_smoke +
                                 Mat_preg_alc + Family_CVD + Age_PHV_c +
                                 BMI_17 + Cortisol_15 + CRP_15 + IL6_24,
                               cov_complete_classic_17,
                               subset = (Child_sex == "Female"))
summary(PWV_Model4_Classic_17_F)

PWV_Model4_Classic_17_M <- lm(PWV_17 ~ Classic_ACEs + BP_systolic_17 +
                                 Age_17_clinic_years + Townsend_sum +
                                 Marital_status + Parent_edu +
                                 Mat_PND_gest + Mat_age_delivery +
                                 Birth_weight_kg + Mat_preg_smoke +
                                 Mat_preg_alc + Family_CVD + Age_PHV_c +
                                 BMI_17 + Cortisol_15 + CRP_15 + IL6_24,
                               cov_complete_classic_17,
                               subset = (Child_sex == "Male"))
summary(PWV_Model4_Classic_17_M)

PWV_Model4_Classic_24_F <- lm(PWV_24 ~ Classic_ACEs + BP_systolic_24 +
                                 Age_24_clinic_years + Townsend_sum +
                                 Marital_status + Parent_edu +
                                 Mat_PND_gest + Mat_age_delivery +
                                 Birth_weight_kg + Mat_preg_smoke +
                                 Mat_preg_alc + Family_CVD + Age_PHV_c +
                                 BMI_24 + Cortisol_15 + CRP_15 + IL6_24,
                               cov_complete_classic_24,
                               subset = (Child_sex == "Female"))
summary(PWV_Model4_Classic_24_F)

PWV_Model4_Classic_24_M <- lm(PWV_24 ~ Classic_ACEs + BP_systolic_24 +
                                 Age_24_clinic_years + Townsend_sum +
                                 Marital_status + Parent_edu +
                                 Mat_PND_gest + Mat_age_delivery +
                                 Birth_weight_kg + Mat_preg_smoke +
                                 Mat_preg_alc + Family_CVD + Age_PHV_c +
                                 BMI_24 + Cortisol_15 + CRP_15 + IL6_24,
                               cov_complete_classic_24,
                               subset = (Child_sex == "Male"))
summary(PWV_Model4_Classic_24_M)

# Assumption checks
avPlots(PWV_Model4_Classic_17_F); plot(PWV_Model4_Classic_17_F)
avPlots(PWV_Model4_Classic_17_M); plot(PWV_Model4_Classic_17_M)
avPlots(PWV_Model4_Classic_24_F); plot(PWV_Model4_Classic_24_F)
avPlots(PWV_Model4_Classic_24_M); plot(PWV_Model4_Classic_24_M)

# cIMT ----

cIMT_Model4_Classic_17_Whole <- lm(cIMT_17 ~ Classic_ACEs +
                                     Age_17_clinic_years + Townsend_sum +
                                     Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_17 + Cortisol_15 + CRP_15 + IL6_24,
                                   cov_complete_classic_17)
summary(cIMT_Model4_Classic_17_Whole)

cIMT_Model4_Classic_24_Whole <- lm(cIMT_24 ~ Classic_ACEs +
                                     Age_24_clinic_years + Townsend_sum +
                                     Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_24 + Cortisol_15 + CRP_15 + IL6_24,
                                   cov_complete_classic_24)
summary(cIMT_Model4_Classic_24_Whole)

cIMT_Model4_Classic_17_F <- lm(cIMT_17 ~ Classic_ACEs + Age_17_clinic_years +
                                  Townsend_sum + Marital_status + Parent_edu +
                                  Mat_PND_gest + Mat_age_delivery +
                                  Birth_weight_kg + Mat_preg_smoke +
                                  Mat_preg_alc + Family_CVD + Age_PHV_c +
                                  BMI_17 + Cortisol_15 + CRP_15 + IL6_24,
                                cov_complete_classic_17,
                                subset = (Child_sex == "Female"))
summary(cIMT_Model4_Classic_17_F)

cIMT_Model4_Classic_17_M <- lm(cIMT_17 ~ Classic_ACEs + Age_17_clinic_years +
                                  Townsend_sum + Marital_status + Parent_edu +
                                  Mat_PND_gest + Mat_age_delivery +
                                  Birth_weight_kg + Mat_preg_smoke +
                                  Mat_preg_alc + Family_CVD + Age_PHV_c +
                                  BMI_17 + Cortisol_15 + CRP_15 + IL6_24,
                                cov_complete_classic_17,
                                subset = (Child_sex == "Male"))
summary(cIMT_Model4_Classic_17_M)

cIMT_Model4_Classic_24_F <- lm(cIMT_24 ~ Classic_ACEs + Age_24_clinic_years +
                                  Townsend_sum + Marital_status + Parent_edu +
                                  Mat_PND_gest + Mat_age_delivery +
                                  Birth_weight_kg + Mat_preg_smoke +
                                  Mat_preg_alc + Family_CVD + Age_PHV_c +
                                  BMI_24 + Cortisol_15 + CRP_15 + IL6_24,
                                cov_complete_classic_24,
                                subset = (Child_sex == "Female"))
summary(cIMT_Model4_Classic_24_F)

cIMT_Model4_Classic_24_M <- lm(cIMT_24 ~ Classic_ACEs + Age_24_clinic_years +
                                  Townsend_sum + Marital_status + Parent_edu +
                                  Mat_PND_gest + Mat_age_delivery +
                                  Birth_weight_kg + Mat_preg_smoke +
                                  Mat_preg_alc + Family_CVD + Age_PHV_c +
                                  BMI_24 + Cortisol_15 + CRP_15 + IL6_24,
                                cov_complete_classic_24,
                                subset = (Child_sex == "Male"))
summary(cIMT_Model4_Classic_24_M)

# Assumption checks
avPlots(cIMT_Model4_Classic_17_F); plot(cIMT_Model4_Classic_17_F)
avPlots(cIMT_Model4_Classic_17_M); plot(cIMT_Model4_Classic_17_M)
avPlots(cIMT_Model4_Classic_24_F); plot(cIMT_Model4_Classic_24_F)
avPlots(cIMT_Model4_Classic_24_M); plot(cIMT_Model4_Classic_24_M)

# Arterial Distensibility ----

Dist_Model4_Classic_17_Whole <- lm(Arterial_Dist_17 ~ Classic_ACEs +
                                     BP_systolic_17 + Age_17_clinic_years +
                                     Townsend_sum + Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_17 + Cortisol_15 + CRP_15 + IL6_24,
                                   cov_complete_classic_17)
summary(Dist_Model4_Classic_17_Whole)

Dist_Model4_Classic_17_F <- lm(Arterial_Dist_17 ~ Classic_ACEs +
                                  BP_systolic_17 + Age_17_clinic_years +
                                  Townsend_sum + Marital_status + Parent_edu +
                                  Mat_PND_gest + Mat_age_delivery +
                                  Birth_weight_kg + Mat_preg_smoke +
                                  Mat_preg_alc + Family_CVD + Age_PHV_c +
                                  BMI_17 + Cortisol_15 + CRP_15 + IL6_24,
                                cov_complete_classic_17,
                                subset = (Child_sex == "Female"))
summary(Dist_Model4_Classic_17_F)

Dist_Model4_Classic_17_M <- lm(Arterial_Dist_17 ~ Classic_ACEs +
                                  BP_systolic_17 + Age_17_clinic_years +
                                  Townsend_sum + Marital_status + Parent_edu +
                                  Mat_PND_gest + Mat_age_delivery +
                                  Birth_weight_kg + Mat_preg_smoke +
                                  Mat_preg_alc + Family_CVD + Age_PHV_c +
                                  BMI_17 + Cortisol_15 + CRP_15 + IL6_24,
                                cov_complete_classic_17,
                                subset = (Child_sex == "Male"))
summary(Dist_Model4_Classic_17_M)

# Assumption checks
avPlots(Dist_Model4_Classic_17_F); plot(Dist_Model4_Classic_17_F)
avPlots(Dist_Model4_Classic_17_M); plot(Dist_Model4_Classic_17_M)


# Model 5: fully adjusted (lifestyle + biological mediators) ----
# Note: issues with convergence/levels when adding smoking/alcohol in some
# subgroups; retained where possible.

# PWV ----

PWV_Model5_Classic_17_Whole <- lm(PWV_17 ~ Classic_ACEs + BP_systolic_17 +
                                    Age_17_clinic_years + Townsend_sum +
                                    Marital_status + Parent_edu +
                                    Mat_PND_gest + Mat_age_delivery +
                                    Birth_weight_kg + Mat_preg_smoke +
                                    Mat_preg_alc + Family_CVD + Age_PHV_c +
                                    BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                    Diet_pattern_13_calib + Non_milk_sugar_13 +
                                    Weekday_sleep_duration_15y +
                                    Weekend_sleep_duration_15y + Child_alc_15 +
                                    Cortisol_15 + CRP_15 + IL6_24,
                                  cov_complete_classic_17, subset = PA_include_15 == TRUE & Diet_include_13 == TRUE)
summary(PWV_Model5_Classic_17_Whole)

PWV_Model5_Classic_24_Whole <- lm(PWV_24 ~ Classic_ACEs + BP_systolic_24 +
                                    Age_24_clinic_years + Townsend_sum +
                                    Marital_status + Parent_edu +
                                    Mat_PND_gest + Mat_age_delivery +
                                    Birth_weight_kg + Mat_preg_smoke +
                                    Mat_preg_alc + Family_CVD + Age_PHV_c +
                                    BMI_24 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                    Diet_pattern_13_calib + Non_milk_sugar_13 +
                                    Weekday_sleep_duration_15y +
                                    Weekend_sleep_duration_15y + Child_alc_15 +
                                    Cortisol_15 + CRP_15 + IL6_24,
                                  cov_complete_classic_24, subset = PA_include_15 == TRUE & Diet_include_13 == TRUE)
summary(PWV_Model5_Classic_24_Whole)

PWV_Model5_Classic_17_F <- lm(PWV_17 ~ Classic_ACEs + BP_systolic_17 +
                                 Age_17_clinic_years + Townsend_sum +
                                 Marital_status + Parent_edu +
                                 Mat_PND_gest + Mat_age_delivery +
                                 Birth_weight_kg + Mat_preg_smoke +
                                 Mat_preg_alc + Family_CVD + Age_PHV_c +
                                 BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                 Diet_pattern_13_calib + Non_milk_sugar_13 +
                                 Weekday_sleep_duration_15y +
                                 Weekend_sleep_duration_15y + Child_alc_15 +
                                 Cortisol_15 + CRP_15 + IL6_24,
                               cov_complete_classic_17,
                               subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(PWV_Model5_Classic_17_F)

PWV_Model5_Classic_17_M <- lm(PWV_17 ~ Classic_ACEs + BP_systolic_17 +
                                 Age_17_clinic_years + Townsend_sum +
                                 Marital_status + Parent_edu +
                                 Mat_PND_gest + Mat_age_delivery +
                                 Birth_weight_kg + Mat_preg_smoke +
                                 Mat_preg_alc + Family_CVD + Age_PHV_c +
                                 BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                 Diet_pattern_13_calib + Non_milk_sugar_13 +
                                 Weekday_sleep_duration_15y +
                                 Weekend_sleep_duration_15y + Child_alc_15 +
                                 Cortisol_15 + CRP_15 + IL6_24,
                               cov_complete_classic_17,
                               subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(PWV_Model5_Classic_17_M)

PWV_Model5_Classic_24_F <- lm(PWV_24 ~ Classic_ACEs + BP_systolic_24 +
                                 Age_24_clinic_years + Townsend_sum +
                                 Marital_status + Parent_edu +
                                 Mat_PND_gest + Mat_age_delivery +
                                 Birth_weight_kg + Mat_preg_smoke +
                                 Mat_preg_alc + Family_CVD + Age_PHV_c +
                                 BMI_24 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                 Diet_pattern_13_calib + Non_milk_sugar_13 +
                                 Weekday_sleep_duration_15y +
                                 Weekend_sleep_duration_15y + Child_alc_15 +
                                 Cortisol_15 + CRP_15 + IL6_24,
                               cov_complete_classic_24,
                               subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(PWV_Model5_Classic_24_F)

PWV_Model5_Classic_24_M <- lm(PWV_24 ~ Classic_ACEs + BP_systolic_24 +
                                 Age_24_clinic_years + Townsend_sum +
                                 Marital_status + Parent_edu +
                                 Mat_PND_gest + Mat_age_delivery +
                                 Birth_weight_kg + Mat_preg_smoke +
                                 Mat_preg_alc + Family_CVD + Age_PHV_c +
                                 BMI_24 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                 Diet_pattern_13_calib + Non_milk_sugar_13 +
                                 Weekday_sleep_duration_15y +
                                 Weekend_sleep_duration_15y + Child_alc_15 +
                                 Cortisol_15 + CRP_15 + IL6_24,
                               cov_complete_classic_24,
                               subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(PWV_Model5_Classic_24_M)

# Assumption checks
avPlots(PWV_Model5_Classic_17_F); plot(PWV_Model5_Classic_17_F)
avPlots(PWV_Model5_Classic_17_M); plot(PWV_Model5_Classic_17_M)
avPlots(PWV_Model5_Classic_24_F); plot(PWV_Model5_Classic_24_F)
avPlots(PWV_Model5_Classic_24_M); plot(PWV_Model5_Classic_24_M)

# cIMT ----

cIMT_Model5_Classic_17_Whole <- lm(cIMT_17 ~ Classic_ACEs +
                                     Age_17_clinic_years + Townsend_sum +
                                     Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                     Diet_pattern_13_calib + Non_milk_sugar_13 +
                                     Weekday_sleep_duration_15y +
                                     Weekend_sleep_duration_15y + Child_alc_15 +
                                     Cortisol_15 + CRP_15 + IL6_24,
                                   cov_complete_classic_17, subset = PA_include_15 == TRUE & Diet_include_13 == TRUE)
summary(cIMT_Model5_Classic_17_Whole)

cIMT_Model5_Classic_24_Whole <- lm(cIMT_24 ~ Classic_ACEs +
                                     Age_24_clinic_years + Townsend_sum +
                                     Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_24 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                     Diet_pattern_13_calib + Non_milk_sugar_13 +
                                     Weekday_sleep_duration_15y +
                                     Weekend_sleep_duration_15y + Child_alc_15 +
                                     Cortisol_15 + CRP_15 + IL6_24,
                                   cov_complete_classic_24, subset = PA_include_15 == TRUE & Diet_include_13 == TRUE)
summary(cIMT_Model5_Classic_24_Whole)

cIMT_Model5_Classic_17_F <- lm(cIMT_17 ~ Classic_ACEs + Age_17_clinic_years +
                                  Townsend_sum + Marital_status + Parent_edu +
                                  Mat_PND_gest + Mat_age_delivery +
                                  Birth_weight_kg + Mat_preg_smoke +
                                  Mat_preg_alc + Family_CVD + Age_PHV_c +
                                  BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                  Diet_pattern_13_calib + Non_milk_sugar_13 +
                                  Weekday_sleep_duration_15y +
                                  Weekend_sleep_duration_15y + Child_alc_15 +
                                  Cortisol_15 + CRP_15 + IL6_24,
                                cov_complete_classic_17,
                                subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(cIMT_Model5_Classic_17_F)

cIMT_Model5_Classic_17_M <- lm(cIMT_17 ~ Classic_ACEs + Age_17_clinic_years +
                                  Townsend_sum + Marital_status + Parent_edu +
                                  Mat_PND_gest + Mat_age_delivery +
                                  Birth_weight_kg + Mat_preg_smoke +
                                  Mat_preg_alc + Family_CVD + Age_PHV_c +
                                  BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                  Diet_pattern_13_calib + Non_milk_sugar_13 +
                                  Weekday_sleep_duration_15y +
                                  Weekend_sleep_duration_15y + Child_alc_15 +
                                  Cortisol_15 + CRP_15 + IL6_24,
                                cov_complete_classic_17,
                                subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(cIMT_Model5_Classic_17_M)

cIMT_Model5_Classic_24_F <- lm(cIMT_24 ~ Classic_ACEs + Age_24_clinic_years +
                                  Townsend_sum + Marital_status + Parent_edu +
                                  Mat_PND_gest + Mat_age_delivery +
                                  Birth_weight_kg + Mat_preg_smoke +
                                  Mat_preg_alc + Family_CVD + Age_PHV_c +
                                  BMI_24 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                  Diet_pattern_13_calib + Non_milk_sugar_13 +
                                  Weekday_sleep_duration_15y +
                                  Weekend_sleep_duration_15y + Child_alc_15 +
                                  Cortisol_15 + CRP_15 + IL6_24,
                                cov_complete_classic_24,
                                subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(cIMT_Model5_Classic_24_F)

cIMT_Model5_Classic_24_M <- lm(cIMT_24 ~ Classic_ACEs + Age_24_clinic_years +
                                  Townsend_sum + Marital_status + Parent_edu +
                                  Mat_PND_gest + Mat_age_delivery +
                                  Birth_weight_kg + Mat_preg_smoke +
                                  Mat_preg_alc + Family_CVD + Age_PHV_c +
                                  BMI_24 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                  Diet_pattern_13_calib + Non_milk_sugar_13 +
                                  Weekday_sleep_duration_15y +
                                  Weekend_sleep_duration_15y + Child_alc_15 +
                                  Cortisol_15 + CRP_15 + IL6_24,
                                cov_complete_classic_24,
                                subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(cIMT_Model5_Classic_24_M)

# Assumption checks
avPlots(cIMT_Model5_Classic_17_F); plot(cIMT_Model5_Classic_17_F)
avPlots(cIMT_Model5_Classic_17_M); plot(cIMT_Model5_Classic_17_M)
avPlots(cIMT_Model5_Classic_24_F); plot(cIMT_Model5_Classic_24_F)
avPlots(cIMT_Model5_Classic_24_M); plot(cIMT_Model5_Classic_24_M)

# Arterial Distensibility ----

Dist_Model5_Classic_17_Whole <- lm(Arterial_Dist_17 ~ Classic_ACEs +
                                     BP_systolic_17 + Age_17_clinic_years +
                                     Townsend_sum + Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                     Diet_pattern_13_calib + Non_milk_sugar_13 +
                                     Weekday_sleep_duration_15y +
                                     Weekend_sleep_duration_15y + Child_alc_15 +
                                     Cortisol_15 + CRP_15 + IL6_24,
                                   cov_complete_classic_17, subset = PA_include_15 == TRUE & Diet_include_13 == TRUE)
summary(Dist_Model5_Classic_17_Whole)

Dist_Model5_Classic_17_F <- lm(Arterial_Dist_17 ~ Classic_ACEs +
                                  BP_systolic_17 + Age_17_clinic_years +
                                  Townsend_sum + Marital_status + Parent_edu +
                                  Mat_PND_gest + Mat_age_delivery +
                                  Birth_weight_kg + Mat_preg_smoke +
                                  Mat_preg_alc + Family_CVD + Age_PHV_c +
                                  BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                  Diet_pattern_13_calib + Non_milk_sugar_13 +
                                  Weekday_sleep_duration_15y +
                                  Weekend_sleep_duration_15y + Child_alc_15 +
                                  Cortisol_15 + CRP_15 + IL6_24,
                                cov_complete_classic_17,
                                subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(Dist_Model5_Classic_17_F)

Dist_Model5_Classic_17_M <- lm(Arterial_Dist_17 ~ Classic_ACEs +
                                  BP_systolic_17 + Age_17_clinic_years +
                                  Townsend_sum + Marital_status + Parent_edu +
                                  Mat_PND_gest + Mat_age_delivery +
                                  Birth_weight_kg + Mat_preg_smoke +
                                  Mat_preg_alc + Family_CVD + Age_PHV_c +
                                  BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                  Diet_pattern_13_calib + Non_milk_sugar_13 +
                                  Weekday_sleep_duration_15y +
                                  Weekend_sleep_duration_15y + Child_alc_15 +
                                  Cortisol_15 + CRP_15 + IL6_24,
                                cov_complete_classic_17,
                                subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(Dist_Model5_Classic_17_M)

# Assumption checks
avPlots(Dist_Model5_Classic_17_F); plot(Dist_Model5_Classic_17_F)
avPlots(Dist_Model5_Classic_17_M); plot(Dist_Model5_Classic_17_M)


# CLASSIC ACEs — CATEGORICAL EXPOSURE ----

# Model 1: unadjusted (age + BP only) ----

# PWV ----

PWV_Model1_CatClassic_17_Whole <- lm(PWV_17 ~ Classic_ACEs_Cat +
                                        BP_systolic_17 + Age_17_clinic_years,
                                      cov_complete_classic_17)
summary(PWV_Model1_CatClassic_17_Whole)

PWV_Model1_CatClassic_24_Whole <- lm(PWV_24 ~ Classic_ACEs_Cat +
                                        BP_systolic_24 + Age_24_clinic_years,
                                      cov_complete_classic_24)
summary(PWV_Model1_CatClassic_24_Whole)

PWV_Model1_CatClassic_17_F <- lm(PWV_17 ~ Classic_ACEs_Cat + BP_systolic_17 +
                                    Age_17_clinic_years,
                                  cov_complete_classic_17,
                                  subset = (Child_sex == "Female"))
summary(PWV_Model1_CatClassic_17_F)

PWV_Model1_CatClassic_17_M <- lm(PWV_17 ~ Classic_ACEs_Cat + BP_systolic_17 +
                                    Age_17_clinic_years,
                                  cov_complete_classic_17,
                                  subset = (Child_sex == "Male"))
summary(PWV_Model1_CatClassic_17_M)

PWV_Model1_CatClassic_24_F <- lm(PWV_24 ~ Classic_ACEs_Cat + BP_systolic_24 +
                                    Age_24_clinic_years,
                                  cov_complete_classic_24,
                                  subset = (Child_sex == "Female"))
summary(PWV_Model1_CatClassic_24_F)

PWV_Model1_CatClassic_24_M <- lm(PWV_24 ~ Classic_ACEs_Cat + BP_systolic_24 +
                                    Age_24_clinic_years,
                                  cov_complete_classic_24,
                                  subset = (Child_sex == "Male"))
summary(PWV_Model1_CatClassic_24_M)

# Assumption checks
ggplot(cov_complete_classic_17,
       aes(x = Classic_ACEs_Cat, y = PWV_17)) +
  geom_point() + geom_smooth(method = "lm")
check_model(PWV_Model1_CatClassic_17_Whole)

ggplot(cov_complete_classic_24,
       aes(x = Classic_ACEs_Cat, y = PWV_24)) +
  geom_point() + geom_smooth(method = "lm")
check_model(PWV_Model1_CatClassic_24_Whole)

ggplot(subset(cov_complete_classic_17, Child_sex == "Female"),
       aes(x = Classic_ACEs_Cat, y = PWV_17)) +
  geom_point() + geom_smooth(method = "lm")
check_model(PWV_Model1_CatClassic_17_F)

ggplot(subset(cov_complete_classic_17, Child_sex == "Male"),
       aes(x = Classic_ACEs_Cat, y = PWV_17)) +
  geom_point() + geom_smooth(method = "lm")
check_model(PWV_Model1_CatClassic_17_M)

ggplot(subset(cov_complete_classic_24, Child_sex == "Female"),
       aes(x = Classic_ACEs_Cat, y = PWV_24)) +
  geom_point() + geom_smooth(method = "lm")
check_model(PWV_Model1_CatClassic_24_F)

ggplot(subset(cov_complete_classic_24, Child_sex == "Male"),
       aes(x = Classic_ACEs_Cat, y = PWV_24)) +
  geom_point() + geom_smooth(method = "lm")
check_model(PWV_Model1_CatClassic_24_M)

# cIMT ----

cIMT_Model1_CatClassic_17_Whole <- lm(cIMT_17 ~ Classic_ACEs_Cat +
                                         Age_17_clinic_years,
                                       cov_complete_classic_17)
summary(cIMT_Model1_CatClassic_17_Whole)

cIMT_Model1_CatClassic_24_Whole <- lm(cIMT_24 ~ Classic_ACEs_Cat +
                                         Age_24_clinic_years,
                                       cov_complete_classic_24)
summary(cIMT_Model1_CatClassic_24_Whole)

cIMT_Model1_CatClassic_17_F <- lm(cIMT_17 ~ Classic_ACEs_Cat +
                                     Age_17_clinic_years,
                                   cov_complete_classic_17,
                                   subset = (Child_sex == "Female"))
summary(cIMT_Model1_CatClassic_17_F)

cIMT_Model1_CatClassic_17_M <- lm(cIMT_17 ~ Classic_ACEs_Cat +
                                     Age_17_clinic_years,
                                   cov_complete_classic_17,
                                   subset = (Child_sex == "Male"))
summary(cIMT_Model1_CatClassic_17_M)

cIMT_Model1_CatClassic_24_F <- lm(cIMT_24 ~ Classic_ACEs_Cat +
                                     Age_24_clinic_years,
                                   cov_complete_classic_24,
                                   subset = (Child_sex == "Female"))
summary(cIMT_Model1_CatClassic_24_F)

cIMT_Model1_CatClassic_24_M <- lm(cIMT_24 ~ Classic_ACEs_Cat +
                                     Age_24_clinic_years,
                                   cov_complete_classic_24,
                                   subset = (Child_sex == "Male"))
summary(cIMT_Model1_CatClassic_24_M)

# Assumption checks
ggplot(cov_complete_classic_17,
       aes(x = Classic_ACEs_Cat, y = cIMT_17)) +
  geom_point() + geom_smooth(method = "lm")
check_model(cIMT_Model1_CatClassic_17_Whole)

ggplot(cov_complete_classic_24,
       aes(x = Classic_ACEs_Cat, y = cIMT_24)) +
  geom_point() + geom_smooth(method = "lm")
check_model(cIMT_Model1_CatClassic_24_Whole)

ggplot(subset(cov_complete_classic_17, Child_sex == "Female"),
       aes(x = Classic_ACEs_Cat, y = cIMT_17)) +
  geom_point() + geom_smooth(method = "lm")
check_model(cIMT_Model1_CatClassic_17_F)

ggplot(subset(cov_complete_classic_17, Child_sex == "Male"),
       aes(x = Classic_ACEs_Cat, y = cIMT_17)) +
  geom_point() + geom_smooth(method = "lm")
check_model(cIMT_Model1_CatClassic_17_M)

ggplot(subset(cov_complete_classic_24, Child_sex == "Female"),
       aes(x = Classic_ACEs_Cat, y = cIMT_24)) +
  geom_point() + geom_smooth(method = "lm")
check_model(cIMT_Model1_CatClassic_24_F)

ggplot(subset(cov_complete_classic_24, Child_sex == "Male"),
       aes(x = Classic_ACEs_Cat, y = cIMT_24)) +
  geom_point() + geom_smooth(method = "lm")
check_model(cIMT_Model1_CatClassic_24_M)

# Arterial Distensibility ----

Dist_Model1_CatClassic_17_Whole <- lm(Arterial_Dist_17 ~ Classic_ACEs_Cat +
                                         BP_systolic_17 + Age_17_clinic_years,
                                       cov_complete_classic_17)
summary(Dist_Model1_CatClassic_17_Whole)

Dist_Model1_CatClassic_17_F <- lm(Arterial_Dist_17 ~ Classic_ACEs_Cat +
                                     BP_systolic_17 + Age_17_clinic_years,
                                   cov_complete_classic_17,
                                   subset = (Child_sex == "Female"))
summary(Dist_Model1_CatClassic_17_F)


Dist_Model1_CatClassic_17_M <- lm(Arterial_Dist_17 ~ Classic_ACEs_Cat +
                                     BP_systolic_17 + Age_17_clinic_years,
                                   cov_complete_classic_17,
                                   subset = (Child_sex == "Male"))
summary(Dist_Model1_CatClassic_17_M)

# Assumption checks
ggplot(cov_complete_classic_17,
       aes(x = Classic_ACEs_Cat, y = Arterial_Dist_17)) +
  geom_point() + geom_smooth(method = "lm")
check_model(Dist_Model1_CatClassic_17_Whole)

ggplot(subset(cov_complete_classic_17, Child_sex == "Female"),
       aes(x = Classic_ACEs_Cat, y = Arterial_Dist_17)) +
  geom_point() + geom_smooth(method = "lm")
check_model(Dist_Model1_CatClassic_17_F)

ggplot(subset(cov_complete_classic_17, Child_sex == "Male"),
       aes(x = Classic_ACEs_Cat, y = Arterial_Dist_17)) +
  geom_point() + geom_smooth(method = "lm")
# Note: check_model applied to male model (corrected from female by copy-paste)
check_model(Dist_Model1_CatClassic_17_M)


# Model 2: sociodemographic covariates ----
# Note: ethnicity excluded due to insufficient level variation.

# PWV ----

PWV_Model2_CatClassic_17_Whole <- lm(PWV_17 ~ Classic_ACEs_Cat * Child_sex +
                                        BP_systolic_17 + Age_17_clinic_years +
                                        Townsend_sum + Marital_status +
                                        Parent_edu + Mat_PND_gest +
                                        Mat_age_delivery + Birth_weight_kg +
                                        Mat_preg_smoke + Mat_preg_alc +
                                        Family_CVD + Age_PHV_c + BMI_17,
                                      cov_complete_classic_17)
summary(PWV_Model2_CatClassic_17_Whole)

PWV_Model2_CatClassic_24_Whole <- lm(PWV_24 ~ Classic_ACEs_Cat * Child_sex +
                                       BP_systolic_24 + Age_24_clinic_years +
                                        Townsend_sum + Marital_status +
                                        Parent_edu + Mat_PND_gest +
                                        Mat_age_delivery + Birth_weight_kg +
                                        Mat_preg_smoke + Mat_preg_alc +
                                        Family_CVD + Age_PHV_c + BMI_24,
                                      cov_complete_classic_24)
summary(PWV_Model2_CatClassic_24_Whole)

PWV_Model2_CatClassic_17_F <- lm(PWV_17 ~ Classic_ACEs_Cat + BP_systolic_17 +
                                    Age_17_clinic_years + Townsend_sum +
                                    Marital_status + Parent_edu +
                                    Mat_PND_gest + Mat_age_delivery +
                                    Birth_weight_kg + Mat_preg_smoke +
                                    Mat_preg_alc + Family_CVD + Age_PHV_c +
                                    BMI_17,
                                  cov_complete_classic_17,
                                  subset = (Child_sex == "Female"))
summary(PWV_Model2_CatClassic_17_F)

PWV_Model2_CatClassic_17_M <- lm(PWV_17 ~ Classic_ACEs_Cat + BP_systolic_17 +
                                    Age_17_clinic_years + Townsend_sum +
                                    Marital_status + Parent_edu +
                                    Mat_PND_gest + Mat_age_delivery +
                                    Birth_weight_kg + Mat_preg_smoke +
                                    Mat_preg_alc + Family_CVD + Age_PHV_c +
                                    BMI_17,
                                  cov_complete_classic_17,
                                  subset = (Child_sex == "Male"))
summary(PWV_Model2_CatClassic_17_M)

PWV_Model2_CatClassic_24_F <- lm(PWV_24 ~ Classic_ACEs_Cat + BP_systolic_24 +
                                    Age_24_clinic_years + Townsend_sum +
                                    Marital_status + Parent_edu +
                                    Mat_PND_gest + Mat_age_delivery +
                                    Birth_weight_kg + Mat_preg_smoke +
                                    Mat_preg_alc + Family_CVD + Age_PHV_c +
                                    BMI_24,
                                  cov_complete_classic_24,
                                  subset = (Child_sex == "Female"))
summary(PWV_Model2_CatClassic_24_F)

PWV_Model2_CatClassic_24_M <- lm(PWV_24 ~ Classic_ACEs_Cat + BP_systolic_24 +
                                    Age_24_clinic_years + Townsend_sum +
                                    Marital_status + Parent_edu +
                                    Mat_PND_gest + Mat_age_delivery +
                                    Birth_weight_kg + Mat_preg_smoke +
                                    Mat_preg_alc + Family_CVD + Age_PHV_c +
                                    BMI_24,
                                  cov_complete_classic_24,
                                  subset = (Child_sex == "Male"))
summary(PWV_Model2_CatClassic_24_M)

# Assumption checks
avPlots(PWV_Model2_CatClassic_17_Whole); plot(PWV_Model2_CatClassic_17_Whole)
avPlots(PWV_Model2_CatClassic_24_Whole); plot(PWV_Model2_CatClassic_24_Whole)
avPlots(PWV_Model2_CatClassic_17_F);     plot(PWV_Model2_CatClassic_17_F)
avPlots(PWV_Model2_CatClassic_17_M);     plot(PWV_Model2_CatClassic_17_M)
avPlots(PWV_Model2_CatClassic_24_F);     plot(PWV_Model2_CatClassic_24_F)
avPlots(PWV_Model2_CatClassic_24_M);     plot(PWV_Model2_CatClassic_24_M)

# cIMT ----

cIMT_Model2_CatClassic_17_Whole <- lm(cIMT_17 ~ Classic_ACEs_Cat * Child_sex +
                                         Age_17_clinic_years + Townsend_sum +
                                         Marital_status + Parent_edu +
                                         Mat_PND_gest + Mat_age_delivery +
                                         Birth_weight_kg + Mat_preg_smoke +
                                         Mat_preg_alc + Family_CVD + Age_PHV_c +
                                         BMI_17,
                                       cov_complete_classic_17)
summary(cIMT_Model2_CatClassic_17_Whole)

cIMT_Model2_CatClassic_24_Whole <- lm(cIMT_24 ~ Classic_ACEs_Cat * Child_sex +
                                         Age_24_clinic_years + Townsend_sum +
                                         Marital_status + Parent_edu +
                                         Mat_PND_gest + Mat_age_delivery +
                                         Birth_weight_kg + Mat_preg_smoke +
                                         Mat_preg_alc + Family_CVD + Age_PHV_c +
                                         BMI_24,
                                       cov_complete_classic_24)
summary(cIMT_Model2_CatClassic_24_Whole)

cIMT_Model2_CatClassic_17_F <- lm(cIMT_17 ~ Classic_ACEs_Cat +
                                     Age_17_clinic_years + Townsend_sum +
                                     Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_17,
                                   cov_complete_classic_17,
                                   subset = (Child_sex == "Female"))
summary(cIMT_Model2_CatClassic_17_F)

cIMT_Model2_CatClassic_17_M <- lm(cIMT_17 ~ Classic_ACEs_Cat +
                                     Age_17_clinic_years + Townsend_sum +
                                     Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_17,
                                   cov_complete_classic_17,
                                   subset = (Child_sex == "Male"))
summary(cIMT_Model2_CatClassic_17_M)

cIMT_Model2_CatClassic_24_F <- lm(cIMT_24 ~ Classic_ACEs_Cat +
                                     Age_24_clinic_years + Townsend_sum +
                                     Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_24,
                                   cov_complete_classic_24,
                                   subset = (Child_sex == "Female"))
summary(cIMT_Model2_CatClassic_24_F)

cIMT_Model2_CatClassic_24_M <- lm(cIMT_24 ~ Classic_ACEs_Cat +
                                     Age_24_clinic_years + Townsend_sum +
                                     Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_24,
                                   cov_complete_classic_24,
                                   subset = (Child_sex == "Male"))
summary(cIMT_Model2_CatClassic_24_M)

# Assumption checks
avPlots(cIMT_Model2_CatClassic_17_Whole); plot(cIMT_Model2_CatClassic_17_Whole)
avPlots(cIMT_Model2_CatClassic_24_Whole); plot(cIMT_Model2_CatClassic_24_Whole)
avPlots(cIMT_Model2_CatClassic_17_F);     plot(cIMT_Model2_CatClassic_17_F)
avPlots(cIMT_Model2_CatClassic_17_M);     plot(cIMT_Model2_CatClassic_17_M)
avPlots(cIMT_Model2_CatClassic_24_F);     plot(cIMT_Model2_CatClassic_24_F)
avPlots(cIMT_Model2_CatClassic_24_M);     plot(cIMT_Model2_CatClassic_24_M)

# Arterial Distensibility ----

Dist_Model2_CatClassic_17_Whole <- lm(Arterial_Dist_17 ~
                                         Classic_ACEs_Cat * Child_sex +
                                         BP_systolic_17 + Age_17_clinic_years +
                                         Townsend_sum + Marital_status +
                                         Parent_edu + Mat_PND_gest +
                                         Mat_age_delivery + Birth_weight_kg +
                                         Mat_preg_smoke + Mat_preg_alc +
                                         Family_CVD + Age_PHV_c + BMI_17,
                                       cov_complete_classic_17)
summary(Dist_Model2_CatClassic_17_Whole)

Dist_Model2_CatClassic_17_F <- lm(Arterial_Dist_17 ~ Classic_ACEs_Cat +
                                     BP_systolic_17 + Age_17_clinic_years +
                                     Townsend_sum + Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_17,
                                   cov_complete_classic_17,
                                   subset = (Child_sex == "Female"))
summary(Dist_Model2_CatClassic_17_F)

Dist_Model2_CatClassic_17_M <- lm(Arterial_Dist_17 ~ Classic_ACEs_Cat +
                                     BP_systolic_17 + Age_17_clinic_years +
                                     Townsend_sum + Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_17,
                                   cov_complete_classic_17,
                                   subset = (Child_sex == "Male"))
summary(Dist_Model2_CatClassic_17_M)

# Assumption checks
avPlots(Dist_Model2_CatClassic_17_Whole); plot(Dist_Model2_CatClassic_17_Whole)
avPlots(Dist_Model2_CatClassic_17_F);     plot(Dist_Model2_CatClassic_17_F)
avPlots(Dist_Model2_CatClassic_17_M);     plot(Dist_Model2_CatClassic_17_M)


# Model 3: + lifestyle mediators ----

# PWV ----
# Note: BP_systolic and age excluded from M3 CatClassic whole-sample due to
# issues with smoking/alcohol convergence in this specification.

PWV_Model3_CatClassic_17_Whole <- lm(PWV_17 ~ Classic_ACEs_Cat +
                                        Townsend_sum + Marital_status +
                                        Parent_edu + Mat_PND_gest +
                                        Mat_age_delivery + Birth_weight_kg +
                                        Mat_preg_smoke + Mat_preg_alc +
                                        Family_CVD + Age_PHV_c + BMI_17 +
                                        Daily_MVPA_15 + Daily_Light_PA_15 +
                                        Diet_pattern_13_calib + Non_milk_sugar_13 +
                                        Weekday_sleep_duration_15y +
                                        Weekend_sleep_duration_15y + Child_alc_15,
                                      cov_complete_classic_17, subset = PA_include_15 == TRUE & Diet_include_13 == TRUE)
summary(PWV_Model3_CatClassic_17_Whole)

PWV_Model3_CatClassic_24_Whole <- lm(PWV_24 ~ Classic_ACEs_Cat +
                                        Townsend_sum + Marital_status +
                                        Parent_edu + Mat_PND_gest +
                                        Mat_age_delivery + Birth_weight_kg +
                                        Mat_preg_smoke + Mat_preg_alc +
                                        Family_CVD + Age_PHV_c + BMI_24 +
                                        Daily_MVPA_15 + Daily_Light_PA_15 +
                                        Diet_pattern_13_calib,
                                      cov_complete_classic_24, subset = PA_include_15 == TRUE & Diet_include_13 == TRUE)
summary(PWV_Model3_CatClassic_24_Whole)

PWV_Model3_CatClassic_17_F <- lm(PWV_17 ~ Classic_ACEs_Cat +
                                    Townsend_sum + Marital_status + Parent_edu +
                                    Mat_PND_gest + Mat_age_delivery +
                                    Birth_weight_kg + Mat_preg_smoke +
                                    Mat_preg_alc + Family_CVD + Age_PHV_c +
                                    BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                    Diet_pattern_13_calib + Non_milk_sugar_13 +
                                    Weekday_sleep_duration_15y +
                                    Weekend_sleep_duration_15y + Child_alc_15,
                                  cov_complete_classic_17,
                                  subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(PWV_Model3_CatClassic_17_F)

PWV_Model3_CatClassic_17_M <- lm(PWV_17 ~ Classic_ACEs_Cat +
                                    Townsend_sum + Marital_status + Parent_edu +
                                    Mat_PND_gest + Mat_age_delivery +
                                    Birth_weight_kg + Mat_preg_smoke +
                                    Mat_preg_alc + Family_CVD + Age_PHV_c +
                                    BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                    Diet_pattern_13_calib + Non_milk_sugar_13 +
                                    Weekday_sleep_duration_15y +
                                    Weekend_sleep_duration_15y + Child_alc_15,
                                  cov_complete_classic_17,
                                  subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(PWV_Model3_CatClassic_17_M)

PWV_Model3_CatClassic_24_F <- lm(PWV_24 ~ Classic_ACEs_Cat +
                                    Townsend_sum + Marital_status + Parent_edu +
                                    Mat_PND_gest + Mat_age_delivery +
                                    Birth_weight_kg + Mat_preg_smoke +
                                    Mat_preg_alc + Family_CVD + Age_PHV_c +
                                    BMI_24 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                    Diet_pattern_13_calib + Non_milk_sugar_13 +
                                    Weekday_sleep_duration_15y +
                                    Weekend_sleep_duration_15y + Child_alc_15,
                                  cov_complete_classic_24,
                                  subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(PWV_Model3_CatClassic_24_F)

PWV_Model3_CatClassic_24_M <- lm(PWV_24 ~ Classic_ACEs_Cat +
                                    Townsend_sum + Marital_status + Parent_edu +
                                    Mat_PND_gest + Mat_age_delivery +
                                    Birth_weight_kg + Mat_preg_smoke +
                                    Mat_preg_alc + Family_CVD + Age_PHV_c +
                                    BMI_24 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                    Diet_pattern_13_calib + Non_milk_sugar_13 +
                                    Weekday_sleep_duration_15y +
                                    Weekend_sleep_duration_15y + Child_alc_15,
                                  cov_complete_classic_24,
                                  subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(PWV_Model3_CatClassic_24_M)

# Assumption checks
avPlots(PWV_Model3_CatClassic_17_F); check_model(PWV_Model3_CatClassic_17_F)
avPlots(PWV_Model3_CatClassic_17_M); check_model(PWV_Model3_CatClassic_17_M)
avPlots(PWV_Model3_CatClassic_24_F); check_model(PWV_Model3_CatClassic_24_F)
avPlots(PWV_Model3_CatClassic_24_M); check_model(PWV_Model3_CatClassic_24_M)

# cIMT ----

cIMT_Model3_CatClassic_17_Whole <- lm(cIMT_17 ~ Classic_ACEs_Cat +
                                         Age_17_clinic_years + Townsend_sum +
                                         Marital_status + Parent_edu +
                                         Mat_PND_gest + Mat_age_delivery +
                                         Birth_weight_kg + Mat_preg_smoke +
                                         Mat_preg_alc + Family_CVD + Age_PHV_c +
                                         BMI_17 + Daily_MVPA_15 +
                                         Daily_Light_PA_15 +
                                         Diet_pattern_13_calib + Non_milk_sugar_13 +
                                         Weekday_sleep_duration_15y +
                                         Weekend_sleep_duration_15y + Child_alc_15,
                                       data = cov_complete_classic_17, subset = PA_include_15 == TRUE & Diet_include_13 == TRUE)
summary(cIMT_Model3_CatClassic_17_Whole)

cIMT_Model3_CatClassic_24_Whole <- lm(cIMT_24 ~ Classic_ACEs_Cat +
                                         Age_24_clinic_years + Townsend_sum +
                                         Marital_status + Parent_edu +
                                         Mat_PND_gest + Mat_age_delivery +
                                         Birth_weight_kg + Mat_preg_smoke +
                                         Mat_preg_alc + Family_CVD + Age_PHV_c +
                                         BMI_24 + Daily_MVPA_15 +
                                         Daily_Light_PA_15 +
                                         Diet_pattern_13_calib + Non_milk_sugar_13 +
                                         Weekday_sleep_duration_15y +
                                         Weekend_sleep_duration_15y + Child_alc_15,
                                       data = cov_complete_classic_24, subset = PA_include_15 == TRUE & Diet_include_13 == TRUE)
summary(cIMT_Model3_CatClassic_24_Whole)

cIMT_Model3_CatClassic_17_F <- lm(cIMT_17 ~ Classic_ACEs_Cat +
                                     Age_17_clinic_years + Townsend_sum +
                                     Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                     Diet_pattern_13_calib + Non_milk_sugar_13 +
                                     Weekday_sleep_duration_15y +
                                     Weekend_sleep_duration_15y + Child_alc_15,
                                   data = cov_complete_classic_17,
                                   subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(cIMT_Model3_CatClassic_17_F)

cIMT_Model3_CatClassic_17_M <- lm(cIMT_17 ~ Classic_ACEs_Cat +
                                     Age_17_clinic_years + Townsend_sum +
                                     Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                     Diet_pattern_13_calib + Non_milk_sugar_13 +
                                     Weekday_sleep_duration_15y +
                                     Weekend_sleep_duration_15y + Child_alc_15,
                                   data = cov_complete_classic_17,
                                   subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(cIMT_Model3_CatClassic_17_M)

cIMT_Model3_CatClassic_24_F <- lm(cIMT_24 ~ Classic_ACEs_Cat +
                                     Age_24_clinic_years + Townsend_sum +
                                     Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_24 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                     Diet_pattern_13_calib + Non_milk_sugar_13 +
                                     Weekday_sleep_duration_15y +
                                     Weekend_sleep_duration_15y + Child_alc_15,
                                   data = cov_complete_classic_24,
                                   subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(cIMT_Model3_CatClassic_24_F)

cIMT_Model3_CatClassic_24_M <- lm(cIMT_24 ~ Classic_ACEs_Cat +
                                     Age_24_clinic_years + Townsend_sum +
                                     Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_24 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                     Diet_pattern_13_calib + Non_milk_sugar_13 +
                                     Weekday_sleep_duration_15y +
                                     Weekend_sleep_duration_15y + Child_alc_15,
                                   data = cov_complete_classic_24,
                                   subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(cIMT_Model3_CatClassic_24_M)

# Assumption checks
avPlots(cIMT_Model3_CatClassic_17_F); plot(cIMT_Model3_CatClassic_17_F)
avPlots(cIMT_Model3_CatClassic_17_M); plot(cIMT_Model3_CatClassic_17_M)
avPlots(cIMT_Model3_CatClassic_24_F); plot(cIMT_Model3_CatClassic_24_F)
avPlots(cIMT_Model3_CatClassic_24_M); plot(cIMT_Model3_CatClassic_24_M)

# Arterial Distensibility ----

Dist_Model3_CatClassic_17_Whole <- lm(Arterial_Dist_17 ~ Classic_ACEs_Cat +
                                         BP_systolic_17 + Age_17_clinic_years +
                                         Townsend_sum + Marital_status +
                                         Parent_edu + Mat_PND_gest +
                                         Mat_age_delivery + Birth_weight_kg +
                                         Mat_preg_smoke + Mat_preg_alc +
                                         Family_CVD + Age_PHV_c + BMI_17 +
                                         Daily_MVPA_15 + Daily_Light_PA_15 +
                                         Diet_pattern_13_calib + Non_milk_sugar_13 +
                                         Weekday_sleep_duration_15y +
                                         Weekend_sleep_duration_15y + Child_alc_15,
                                       cov_complete_classic_17, subset = PA_include_15 == TRUE & Diet_include_13 == TRUE)
summary(Dist_Model3_CatClassic_17_Whole)

Dist_Model3_CatClassic_17_F <- lm(Arterial_Dist_17 ~ Classic_ACEs_Cat +
                                     BP_systolic_17 + Age_17_clinic_years +
                                     Townsend_sum + Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                     Diet_pattern_13_calib + Non_milk_sugar_13 +
                                     Weekday_sleep_duration_15y +
                                     Weekend_sleep_duration_15y + Child_alc_15,
                                   cov_complete_classic_17,
                                   subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(Dist_Model3_CatClassic_17_F)

Dist_Model3_CatClassic_17_M <- lm(Arterial_Dist_17 ~ Classic_ACEs_Cat +
                                     BP_systolic_17 + Age_17_clinic_years +
                                     Townsend_sum + Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                     Diet_pattern_13_calib + Non_milk_sugar_13 +
                                     Weekday_sleep_duration_15y +
                                     Weekend_sleep_duration_15y + Child_alc_15,
                                   cov_complete_classic_17,
                                   subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(Dist_Model3_CatClassic_17_M)

# Assumption checks
avPlots(Dist_Model3_CatClassic_17_F); plot(Dist_Model3_CatClassic_17_F)
avPlots(Dist_Model3_CatClassic_17_M); plot(Dist_Model3_CatClassic_17_M)


# Model 4: + biological mediators ----

# PWV ----

PWV_Model4_CatClassic_17_Whole <- lm(PWV_17 ~ Classic_ACEs_Cat +
                                        BP_systolic_17 + Age_17_clinic_years +
                                        Townsend_sum + Marital_status +
                                        Parent_edu + Mat_PND_gest +
                                        Mat_age_delivery + Birth_weight_kg +
                                        Mat_preg_smoke + Mat_preg_alc +
                                        Family_CVD + Age_PHV_c + BMI_17 +
                                        Cortisol_15 + CRP_15 + IL6_24,
                                      cov_complete_classic_17)
summary(PWV_Model4_CatClassic_17_Whole)

PWV_Model4_CatClassic_24_Whole <- lm(PWV_24 ~ Classic_ACEs_Cat +
                                        BP_systolic_24 + Age_24_clinic_years +
                                        Townsend_sum + Marital_status +
                                        Parent_edu + Mat_PND_gest +
                                        Mat_age_delivery + Birth_weight_kg +
                                        Mat_preg_smoke + Mat_preg_alc +
                                        Family_CVD + Age_PHV_c + BMI_24 +
                                        Cortisol_15 + CRP_15 + IL6_24,
                                      cov_complete_classic_24)
summary(PWV_Model4_CatClassic_24_Whole)

PWV_Model4_CatClassic_17_F <- lm(PWV_17 ~ Classic_ACEs_Cat + BP_systolic_17 +
                                    Age_17_clinic_years + Townsend_sum +
                                    Marital_status + Parent_edu +
                                    Mat_PND_gest + Mat_age_delivery +
                                    Birth_weight_kg + Mat_preg_smoke +
                                    Mat_preg_alc + Family_CVD + Age_PHV_c +
                                    BMI_17 + Cortisol_15 + CRP_15 + IL6_24,
                                  cov_complete_classic_17,
                                  subset = (Child_sex == "Female"))
summary(PWV_Model4_CatClassic_17_F)

PWV_Model4_CatClassic_17_M <- lm(PWV_17 ~ Classic_ACEs_Cat + BP_systolic_17 +
                                    Age_17_clinic_years + Townsend_sum +
                                    Marital_status + Parent_edu +
                                    Mat_PND_gest + Mat_age_delivery +
                                    Birth_weight_kg + Mat_preg_smoke +
                                    Mat_preg_alc + Family_CVD + Age_PHV_c +
                                    BMI_17 + Cortisol_15 + CRP_15 + IL6_24,
                                  cov_complete_classic_17,
                                  subset = (Child_sex == "Male"))
summary(PWV_Model4_CatClassic_17_M)

PWV_Model4_CatClassic_24_F <- lm(PWV_24 ~ Classic_ACEs_Cat + BP_systolic_24 +
                                    Age_24_clinic_years + Townsend_sum +
                                    Marital_status + Parent_edu +
                                    Mat_PND_gest + Mat_age_delivery +
                                    Birth_weight_kg + Mat_preg_smoke +
                                    Mat_preg_alc + Family_CVD + Age_PHV_c +
                                    BMI_24 + Cortisol_15 + CRP_15 + IL6_24,
                                  cov_complete_classic_24,
                                  subset = (Child_sex == "Female"))
summary(PWV_Model4_CatClassic_24_F)

PWV_Model4_CatClassic_24_M <- lm(PWV_24 ~ Classic_ACEs_Cat + BP_systolic_24 +
                                    Age_24_clinic_years + Townsend_sum +
                                    Marital_status + Parent_edu +
                                    Mat_PND_gest + Mat_age_delivery +
                                    Birth_weight_kg + Mat_preg_smoke +
                                    Mat_preg_alc + Family_CVD + Age_PHV_c +
                                    BMI_24 + Cortisol_15 + CRP_15 + IL6_24,
                                  cov_complete_classic_24,
                                  subset = (Child_sex == "Male"))
summary(PWV_Model4_CatClassic_24_M)

# Assumption checks
avPlots(PWV_Model4_CatClassic_17_F); plot(PWV_Model4_CatClassic_17_F)
avPlots(PWV_Model4_CatClassic_17_M); plot(PWV_Model4_CatClassic_17_M)
avPlots(PWV_Model4_CatClassic_24_F); plot(PWV_Model4_CatClassic_24_F)
avPlots(PWV_Model4_CatClassic_24_M); plot(PWV_Model4_CatClassic_24_M)

# cIMT ----

cIMT_Model4_CatClassic_17_Whole <- lm(cIMT_17 ~ Classic_ACEs_Cat +
                                         Age_17_clinic_years + Townsend_sum +
                                         Marital_status + Parent_edu +
                                         Mat_PND_gest + Mat_age_delivery +
                                         Birth_weight_kg + Mat_preg_smoke +
                                         Mat_preg_alc + Family_CVD + Age_PHV_c +
                                         BMI_17 + Cortisol_15 + CRP_15 + IL6_24,
                                       cov_complete_classic_17)
summary(cIMT_Model4_CatClassic_17_Whole)

cIMT_Model4_CatClassic_24_Whole <- lm(cIMT_24 ~ Classic_ACEs_Cat +
                                         Age_24_clinic_years + Townsend_sum +
                                         Marital_status + Parent_edu +
                                         Mat_PND_gest + Mat_age_delivery +
                                         Birth_weight_kg + Mat_preg_smoke +
                                         Mat_preg_alc + Family_CVD + Age_PHV_c +
                                         BMI_24 + Cortisol_15 + CRP_15 + IL6_24,
                                       cov_complete_classic_24)
summary(cIMT_Model4_CatClassic_24_Whole)

cIMT_Model4_CatClassic_17_F <- lm(cIMT_17 ~ Classic_ACEs_Cat +
                                     Age_17_clinic_years + Townsend_sum +
                                     Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_17 + Cortisol_15 + CRP_15 + IL6_24,
                                   cov_complete_classic_17,
                                   subset = (Child_sex == "Female"))
summary(cIMT_Model4_CatClassic_17_F)

cIMT_Model4_CatClassic_17_M <- lm(cIMT_17 ~ Classic_ACEs_Cat +
                                     Age_17_clinic_years + Townsend_sum +
                                     Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_17 + Cortisol_15 + CRP_15 + IL6_24,
                                   cov_complete_classic_17,
                                   subset = (Child_sex == "Male"))
summary(cIMT_Model4_CatClassic_17_M)

cIMT_Model4_CatClassic_24_F <- lm(cIMT_24 ~ Classic_ACEs_Cat +
                                     Age_24_clinic_years + Townsend_sum +
                                     Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_24 + Cortisol_15 + CRP_15 + IL6_24,
                                   cov_complete_classic_24,
                                   subset = (Child_sex == "Female"))
summary(cIMT_Model4_CatClassic_24_F)

cIMT_Model4_CatClassic_24_M <- lm(cIMT_24 ~ Classic_ACEs_Cat +
                                     Age_24_clinic_years + Townsend_sum +
                                     Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_24 + Cortisol_15 + CRP_15 + IL6_24,
                                   cov_complete_classic_24,
                                   subset = (Child_sex == "Male"))
summary(cIMT_Model4_CatClassic_24_M)

# Assumption checks
avPlots(cIMT_Model4_CatClassic_17_F); plot(cIMT_Model4_CatClassic_17_F)
avPlots(cIMT_Model4_CatClassic_17_M); plot(cIMT_Model4_CatClassic_17_M)
avPlots(cIMT_Model4_CatClassic_24_F); plot(cIMT_Model4_CatClassic_24_F)
avPlots(cIMT_Model4_CatClassic_24_M); plot(cIMT_Model4_CatClassic_24_M)

# Arterial Distensibility ----

Dist_Model4_CatClassic_17_Whole <- lm(Arterial_Dist_17 ~ Classic_ACEs_Cat +
                                         BP_systolic_17 + Age_17_clinic_years +
                                         Townsend_sum + Marital_status +
                                         Parent_edu + Mat_PND_gest +
                                         Mat_age_delivery + Birth_weight_kg +
                                         Mat_preg_smoke + Mat_preg_alc +
                                         Family_CVD + Age_PHV_c + BMI_17 +
                                         Cortisol_15 + CRP_15 + IL6_24,
                                       cov_complete_classic_17)
summary(Dist_Model4_CatClassic_17_Whole)

Dist_Model4_CatClassic_17_F <- lm(Arterial_Dist_17 ~ Classic_ACEs_Cat +
                                     BP_systolic_17 + Age_17_clinic_years +
                                     Townsend_sum + Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_17 + Cortisol_15 + CRP_15 + IL6_24,
                                   cov_complete_classic_17,
                                   subset = (Child_sex == "Female"))
summary(Dist_Model4_CatClassic_17_F)

Dist_Model4_CatClassic_17_M <- lm(Arterial_Dist_17 ~ Classic_ACEs_Cat +
                                     BP_systolic_17 + Age_17_clinic_years +
                                     Townsend_sum + Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_17 + Cortisol_15 + CRP_15 + IL6_24,
                                   cov_complete_classic_17,
                                   subset = (Child_sex == "Male"))
summary(Dist_Model4_CatClassic_17_M)

# Assumption checks
avPlots(Dist_Model4_CatClassic_17_F); plot(Dist_Model4_CatClassic_17_F)
avPlots(Dist_Model4_CatClassic_17_M); plot(Dist_Model4_CatClassic_17_M)


# Model 5: fully adjusted ----
# Note: issues with smoking/alcohol convergence in some subgroups.

# PWV ----

PWV_Model5_CatClassic_17_Whole <- lm(PWV_17 ~ Classic_ACEs_Cat +
                                        BP_systolic_17 + Age_17_clinic_years +
                                        Townsend_sum + Marital_status +
                                        Parent_edu + Mat_PND_gest +
                                        Mat_age_delivery + Birth_weight_kg +
                                        Mat_preg_smoke + Mat_preg_alc +
                                        Family_CVD + Age_PHV_c + BMI_17 +
                                        Daily_MVPA_15 + Daily_Light_PA_15 +
                                        Diet_pattern_13_calib + Non_milk_sugar_13 +
                                        Weekday_sleep_duration_15y +
                                        Weekend_sleep_duration_15y + Child_alc_15 +
                                        Cortisol_15 + CRP_15 + IL6_24,
                                      cov_complete_classic_17, subset = PA_include_15 == TRUE & Diet_include_13 == TRUE)
summary(PWV_Model5_CatClassic_17_Whole)

PWV_Model5_CatClassic_24_Whole <- lm(PWV_24 ~ Classic_ACEs_Cat +
                                        BP_systolic_24 + Age_24_clinic_years +
                                        Townsend_sum + Marital_status +
                                        Parent_edu + Mat_PND_gest +
                                        Mat_age_delivery + Birth_weight_kg +
                                        Mat_preg_smoke + Mat_preg_alc +
                                        Family_CVD + Age_PHV_c + BMI_24 +
                                        Daily_MVPA_15 + Daily_Light_PA_15 +
                                        Diet_pattern_13_calib + Non_milk_sugar_13 +
                                        Weekday_sleep_duration_15y +
                                        Weekend_sleep_duration_15y + Child_alc_15 +
                                        Cortisol_15 + CRP_15 + IL6_24,
                                      cov_complete_classic_24, subset = PA_include_15 == TRUE & Diet_include_13 == TRUE)
summary(PWV_Model5_CatClassic_24_Whole)

PWV_Model5_CatClassic_17_F <- lm(PWV_17 ~ Classic_ACEs_Cat + BP_systolic_17 +
                                    Age_17_clinic_years + Townsend_sum +
                                    Marital_status + Parent_edu +
                                    Mat_PND_gest + Mat_age_delivery +
                                    Birth_weight_kg + Mat_preg_smoke +
                                    Mat_preg_alc + Family_CVD + Age_PHV_c +
                                    BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                    Diet_pattern_13_calib +
                                    Cortisol_15 + CRP_15 + IL6_24,
                                  cov_complete_classic_17,
                                  subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(PWV_Model5_CatClassic_17_F)

PWV_Model5_CatClassic_17_M <- lm(PWV_17 ~ Classic_ACEs_Cat + BP_systolic_17 +
                                    Age_17_clinic_years + Townsend_sum +
                                    Marital_status + Parent_edu +
                                    Mat_PND_gest + Mat_age_delivery +
                                    Birth_weight_kg + Mat_preg_smoke +
                                    Mat_preg_alc + Family_CVD + Age_PHV_c +
                                    BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                    Diet_pattern_13_calib + Non_milk_sugar_13 +
                                    Weekday_sleep_duration_15y +
                                    Weekend_sleep_duration_15y + Child_alc_15 +
                                    Cortisol_15 + CRP_15 + IL6_24,
                                  cov_complete_classic_17,
                                  subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(PWV_Model5_CatClassic_17_M)

PWV_Model5_CatClassic_24_F <- lm(PWV_24 ~ Classic_ACEs_Cat + BP_systolic_24 +
                                    Age_24_clinic_years + Townsend_sum +
                                    Marital_status + Parent_edu +
                                    Mat_PND_gest + Mat_age_delivery +
                                    Birth_weight_kg + Mat_preg_smoke +
                                    Mat_preg_alc + Family_CVD + Age_PHV_c +
                                    BMI_24 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                    Diet_pattern_13_calib + Non_milk_sugar_13 +
                                    Weekday_sleep_duration_15y +
                                    Weekend_sleep_duration_15y + Child_alc_15 +
                                    Cortisol_15 + CRP_15 + IL6_24,
                                  cov_complete_classic_24,
                                  subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(PWV_Model5_CatClassic_24_F)

PWV_Model5_CatClassic_24_M <- lm(PWV_24 ~ Classic_ACEs_Cat + BP_systolic_24 +
                                    Age_24_clinic_years + Townsend_sum +
                                    Marital_status + Parent_edu +
                                    Mat_PND_gest + Mat_age_delivery +
                                    Birth_weight_kg + Mat_preg_smoke +
                                    Mat_preg_alc + Family_CVD + Age_PHV_c +
                                    BMI_24 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                    Diet_pattern_13_calib + Non_milk_sugar_13 +
                                    Weekday_sleep_duration_15y +
                                    Weekend_sleep_duration_15y + Child_alc_15 +
                                    Cortisol_15 + CRP_15 + IL6_24,
                                  cov_complete_classic_24,
                                  subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(PWV_Model5_CatClassic_24_M)

# Assumption checks
avPlots(PWV_Model5_CatClassic_17_F); plot(PWV_Model5_CatClassic_17_F)
avPlots(PWV_Model5_CatClassic_17_M); plot(PWV_Model5_CatClassic_17_M)
avPlots(PWV_Model5_CatClassic_24_F); plot(PWV_Model5_CatClassic_24_F)
avPlots(PWV_Model5_CatClassic_24_M); plot(PWV_Model5_CatClassic_24_M)

# cIMT ----

cIMT_Model5_CatClassic_17_Whole <- lm(cIMT_17 ~ Classic_ACEs_Cat +
                                         Age_17_clinic_years + Townsend_sum +
                                         Marital_status + Parent_edu +
                                         Mat_PND_gest + Mat_age_delivery +
                                         Birth_weight_kg + Mat_preg_smoke +
                                         Mat_preg_alc + Family_CVD + Age_PHV_c +
                                         BMI_17 + Daily_MVPA_15 +
                                         Daily_Light_PA_15 +
                                         Diet_pattern_13_calib + Non_milk_sugar_13 +
                                         Weekday_sleep_duration_15y +
                                         Weekend_sleep_duration_15y + Child_alc_15 +
                                         Cortisol_15 + CRP_15 + IL6_24,
                                       cov_complete_classic_17, subset = PA_include_15 == TRUE & Diet_include_13 == TRUE)
summary(cIMT_Model5_CatClassic_17_Whole)

cIMT_Model5_CatClassic_24_Whole <- lm(cIMT_24 ~ Classic_ACEs_Cat +
                                         Age_24_clinic_years + Townsend_sum +
                                         Marital_status + Parent_edu +
                                         Mat_PND_gest + Mat_age_delivery +
                                         Birth_weight_kg + Mat_preg_smoke +
                                         Mat_preg_alc + Family_CVD + Age_PHV_c +
                                         BMI_24 + Daily_MVPA_15 +
                                         Daily_Light_PA_15 +
                                         Diet_pattern_13_calib + Non_milk_sugar_13 +
                                         Weekday_sleep_duration_15y +
                                         Weekend_sleep_duration_15y + Child_alc_15 +
                                         Cortisol_15 + CRP_15 + IL6_24,
                                       cov_complete_classic_24, subset = PA_include_15 == TRUE & Diet_include_13 == TRUE)
summary(cIMT_Model5_CatClassic_24_Whole)

cIMT_Model5_CatClassic_17_F <- lm(cIMT_17 ~ Classic_ACEs_Cat +
                                     Age_17_clinic_years + Townsend_sum +
                                     Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                     Diet_pattern_13_calib + Non_milk_sugar_13 +
                                     Weekday_sleep_duration_15y +
                                     Weekend_sleep_duration_15y + Child_alc_15 +
                                     Cortisol_15 + CRP_15 + IL6_24,
                                   cov_complete_classic_17,
                                   subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(cIMT_Model5_CatClassic_17_F)

cIMT_Model5_CatClassic_17_M <- lm(cIMT_17 ~ Classic_ACEs_Cat +
                                     Age_17_clinic_years + Townsend_sum +
                                     Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                     Diet_pattern_13_calib + Non_milk_sugar_13 +
                                     Weekday_sleep_duration_15y +
                                     Weekend_sleep_duration_15y + Child_alc_15 +
                                     Cortisol_15 + CRP_15 + IL6_24,
                                   cov_complete_classic_17,
                                   subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(cIMT_Model5_CatClassic_17_M)

cIMT_Model5_CatClassic_24_F <- lm(cIMT_24 ~ Classic_ACEs_Cat +
                                     Age_24_clinic_years + Townsend_sum +
                                     Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_24 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                     Diet_pattern_13_calib + Non_milk_sugar_13 +
                                     Weekday_sleep_duration_15y +
                                     Weekend_sleep_duration_15y + Child_alc_15 +
                                     Cortisol_15 + CRP_15 + IL6_24,
                                   cov_complete_classic_24,
                                   subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(cIMT_Model5_CatClassic_24_F)

cIMT_Model5_CatClassic_24_M <- lm(cIMT_24 ~ Classic_ACEs_Cat +
                                     Age_24_clinic_years + Townsend_sum +
                                     Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_24 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                     Diet_pattern_13_calib +
                                     Cortisol_15 + CRP_15 + IL6_24,
                                   cov_complete_classic_24,
                                   subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(cIMT_Model5_CatClassic_24_M)

# Assumption checks
avPlots(cIMT_Model5_CatClassic_17_F); plot(cIMT_Model5_CatClassic_17_F)
avPlots(cIMT_Model5_CatClassic_17_M); plot(cIMT_Model5_CatClassic_17_M)
avPlots(cIMT_Model5_CatClassic_24_F); plot(cIMT_Model5_CatClassic_24_F)
avPlots(cIMT_Model5_CatClassic_24_M); plot(cIMT_Model5_CatClassic_24_M)

# Arterial Distensibility ----

Dist_Model5_CatClassic_17_Whole <- lm(Arterial_Dist_17 ~ Classic_ACEs_Cat +
                                         BP_systolic_17 + Age_17_clinic_years +
                                         Townsend_sum + Marital_status +
                                         Parent_edu + Mat_PND_gest +
                                         Mat_age_delivery + Birth_weight_kg +
                                         Mat_preg_smoke + Mat_preg_alc +
                                         Family_CVD + Age_PHV_c + BMI_17 +
                                         Daily_MVPA_15 + Daily_Light_PA_15 +
                                         Diet_pattern_13_calib + Non_milk_sugar_13 +
                                         Weekday_sleep_duration_15y +
                                         Weekend_sleep_duration_15y + Child_alc_15 +
                                         Cortisol_15 + CRP_15 + IL6_24,
                                       cov_complete_classic_17, subset = PA_include_15 == TRUE & Diet_include_13 == TRUE)
summary(Dist_Model5_CatClassic_17_Whole)

Dist_Model5_CatClassic_17_F <- lm(Arterial_Dist_17 ~ Classic_ACEs_Cat +
                                     BP_systolic_17 + Age_17_clinic_years +
                                     Townsend_sum + Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                     Diet_pattern_13_calib + Non_milk_sugar_13 +
                                     Weekday_sleep_duration_15y +
                                     Weekend_sleep_duration_15y + Child_alc_15 +
                                     Cortisol_15 + CRP_15 + IL6_24,
                                   cov_complete_classic_17,
                                   subset = (Child_sex == "Female" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(Dist_Model5_CatClassic_17_F)

Dist_Model5_CatClassic_17_M <- lm(Arterial_Dist_17 ~ Classic_ACEs_Cat +
                                     BP_systolic_17 + Age_17_clinic_years +
                                     Townsend_sum + Marital_status + Parent_edu +
                                     Mat_PND_gest + Mat_age_delivery +
                                     Birth_weight_kg + Mat_preg_smoke +
                                     Mat_preg_alc + Family_CVD + Age_PHV_c +
                                     BMI_17 + Daily_MVPA_15 + Daily_Light_PA_15 +
                                     Diet_pattern_13_calib + Non_milk_sugar_13 +
                                     Weekday_sleep_duration_15y +
                                     Weekend_sleep_duration_15y + Child_alc_15 +
                                     Cortisol_15 + CRP_15 + IL6_24,
                                   cov_complete_classic_17,
                                   subset = (Child_sex == "Male" & PA_include_15 == TRUE & Diet_include_13 == TRUE))
summary(Dist_Model5_CatClassic_17_M)

# Assumption checks
avPlots(Dist_Model5_CatClassic_17_F); plot(Dist_Model5_CatClassic_17_F)
avPlots(Dist_Model5_CatClassic_17_M); plot(Dist_Model5_CatClassic_17_M)
