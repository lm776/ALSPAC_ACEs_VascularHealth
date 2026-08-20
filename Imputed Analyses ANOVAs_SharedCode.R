# ACEs and Vascular Health — Imputed MI Cross-Sectional Group Comparisons (D1 Tests)
# Study: Adverse Childhood Experiences and Changes in Vascular Health from
#        Adolescence to Early Adulthood: Cross-Sectional and Longitudinal
#        Evidence from the ALSPAC Study
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
# Overview:
#   Outcomes: PWV, Arterial_Dist, cIMT at ages 17 and 24
#   Groups:   Whole sample, Female, Male
#   Exposure: Classic_ACEs_cat (categorical, already present in mids object)
#
#   Model 1:
#     PWV/Arterial_Dist: Y ~ Classic_ACEs_cat + age_clinic + BP_systolic
#     cIMT:              Y ~ Classic_ACEs_cat + age_clinic
#
#   Model 2:
#     PWV/Arterial_Dist: + Child_ethnicity + Townsend_sum + Marital_status +
#                          Parent_edu + Mat_PND_gest + Mat_age_delivery +
#                          Birth_weight_kg + Mat_preg_smoke + Mat_preg_alc +
#                          Family_CVD + Age_PHV_c + BMI_(wave) + BP_systolic + age_clinic
#     cIMT:              same covariates but NO BP_systolic
#
#   D1 tests:
#     - Overall ACE group effect: compare full vs reduced model (drop Classic_ACEs_cat)
#     - (Optional) Interaction in whole-sample Model 2: Classic_ACEs_cat * Child_sex
#
# Author:  Laura Macro
# Date:    2026


# Section 1: Required packages ----

library(emmeans)
library(mice)
library(mitools)
library(tidyr)
library(dplyr)

emm_options(rg.limit = 50000)


# Section 2: Helper functions ----

# Wave-specific variable name constructors
y_var   <- function(outcome, age) paste0(outcome, "_", age)          # e.g. PWV_17
age_var <- function(age)          paste0("Age_", age, "_clinic_years") # e.g. Age_17_clinic_years
bp_var  <- function(age)          paste0("BP_systolic_", age)         # e.g. BP_systolic_17
bmi_var <- function(age)          paste0("BMI_", age)                 # e.g. BMI_17

# Model 2 base covariates (non-wave-specific)
m2_covars_base <- c(
  "Child_ethnicity", "Townsend_sum", "Marital_status", "Parent_edu",
  "Mat_PND_gest", "Mat_age_delivery", "Birth_weight_kg", "Mat_preg_smoke", "Mat_preg_alc",
  "Family_CVD", "Age_PHV_c"
)

# Partial eta^2 from D1 result
d1_to_eta2 <- function(d1) {
  if (is.null(d1) || is.null(d1$result)) return(NULL)
  df <- as.data.frame(d1$result)
  names(df) <- c("F_stat", "df1", "df2", "p_value", "riv")
  df$eta2_p <- (df$F_stat * df$df1) / (df$F_stat * df$df1 + df$df2)
  df
}

# Pooled residual SD across imputations (used for Cohen's d)
get_pooled_sigma <- function(mira) {
  sigmas <- sapply(mira$analyses, function(m) summary(m)$sigma)
  mean(sigmas)
}


# Section 3: emmeans and pooled contrasts ----

# Computes emmeans and pairwise contrasts across imputations using Rubin's rules.
# Returns pooled emmeans and pooled contrasts (including Cohen's d).
compute_mi_pairwise_emmeans <- function(mira, interaction = FALSE,
                                        sigma_pooled,
                                        age) {
  bmi    <- bmi_var(age)
  bp     <- bp_var(age)
  agec   <- age_var(age)

  nuisance_base <- c(
    "Marital_status", "Parent_edu",
    "Mat_PND_gest", "Mat_age_delivery", "Birth_weight_kg",
    "Mat_preg_smoke", "Mat_preg_alc", "Family_CVD", "Age_PHV_c"
  )
  if (outcome %in% c("PWV", "Arterial_Dist")) {
    nuisance_wave <- c(bmi, agec, bp)
  } else {  # cIMT
    nuisance_wave <- c(bmi, agec)
  }  
  nuisance_vars <- c(nuisance_base, nuisance_wave)

  emm_list <- lapply(mira$analyses, function(fit) {
    emm <- suppressMessages(
      if (interaction) {
        emmeans::emmeans(
          fit,
          ~ Classic_ACEs_cat | Child_sex,
          nuisance = nuisance_vars,
          rg.limit = 50000
        )
      } else {
        emmeans::emmeans(
          fit,
          ~ Classic_ACEs_cat,
          nuisance = nuisance_vars,
          rg.limit = 50000
        )
      }
    )

    emm_df <- as.data.frame(emm) |>
      dplyr::rename(
        estimate = emmean,
        CI_low   = lower.CL,
        CI_high  = upper.CL
      )

    contr_df <- as.data.frame(contrast(emm, method = "pairwise", adjust = "none"))

    list(
      emmeans   = emm_df,
      contrasts = contr_df
    )
  })

  # Pool EMMs over imputations (Rubin's rules)
  emm_all <- dplyr::bind_rows(lapply(emm_list, `[[`, "emmeans"), .id = "imp")

  emm_pooled <- emm_all |>
    dplyr::group_by(Classic_ACEs_cat, Child_sex) |>
    dplyr::summarise(
      m   = dplyr::n(),
      est = mean(estimate),
      W   = mean(SE^2),
      B   = var(estimate),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      T_var   = W + (1 + 1/m) * B,
      SE      = sqrt(T_var),
      t       = est / SE,
      df      = dplyr::if_else(B == 0, 9999, (m - 1) * (1 + 1/m * W / B)^2),
      p_value = 2 * pt(-abs(t), df = df),
      CI_low  = est - qt(0.975, df) * SE,
      CI_high = est + qt(0.975, df) * SE
    )

  # Pool contrasts over imputations (Rubin's rules)
  contr_all <- dplyr::bind_rows(lapply(emm_list, `[[`, "contrasts"), .id = "imp")

  contrast_pooled <- contr_all |>
    dplyr::group_by(contrast) |>
    dplyr::summarise(
      m   = dplyr::n(),
      est = mean(estimate),
      W   = mean(SE^2),
      B   = var(estimate),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      T_var    = W + (1 + 1/m) * B,
      SE       = sqrt(T_var),
      t        = est / SE,
      df       = dplyr::if_else(B == 0, 9999, (m - 1) * (1 + 1/m * W / B)^2),
      p_value  = 2 * pt(-abs(t), df = df),
      CI_low   = est - qt(0.975, df) * SE,
      CI_high  = est + qt(0.975, df) * SE,
      Cohens_d = est / sigma_pooled
    )

  list(
    emmeans   = emm_pooled,
    contrasts = contrast_pooled
  )
}


# Section 4: Formula builder ----

build_fml <- function(outcome, age,
                      model = c("M1", "M2"),
                      group = c("Whole", "Female", "Male"),
                      include_sex_interaction_M2 = TRUE) {

  model <- match.arg(model)
  group <- match.arg(group)

  y    <- y_var(outcome, age)
  agec <- age_var(age)
  bp   <- bp_var(age)
  bmi  <- bmi_var(age)

  ace <- "Classic_ACEs_cat"

  # Model 1 RHS
  rhs_m1 <- if (outcome %in% c("PWV", "Arterial_Dist")) {
    paste(c(ace, agec, bp), collapse = " + ")
  } else { # cIMT
    paste(c(ace, agec), collapse = " + ")
  }

  if (model == "M1") {
    return(as.formula(paste0(y, " ~ ", rhs_m1)))
  }

  # Model 2 RHS
  rhs_m2 <- c(m2_covars_base, bmi, agec)

  # PWV/Arterial_Dist include BP; cIMT does not
  if (outcome %in% c("PWV", "Arterial_Dist")) rhs_m2 <- c(rhs_m2, bp)

  # Whole group: optionally include ACE x sex interaction in Model 2
  if (group == "Whole") {
    if (include_sex_interaction_M2) {
      rhs <- paste(c(paste0(ace, " * Child_sex"), rhs_m2), collapse = " + ")
    } else {
      rhs <- paste(c(ace, "Child_sex", rhs_m2), collapse = " + ")
    }
  } else {
    # Sex-stratified: do NOT include Child_sex or interaction
    rhs <- paste(c(ace, rhs_m2), collapse = " + ")
  }

  as.formula(paste0(y, " ~ ", rhs))
}


# Section 5: MI fitting and D1 tests ----

fit_mi <- function(mids_obj, fml, subset_expr = NULL) {
  fml_chr <- paste(deparse(fml), collapse = " ")
  with(
    data = mids_obj,
    expr = {
      fml2 <- stats::as.formula(fml_chr)
      environment(fml2) <- environment()
      if (is.null(subset_expr)) {
        lm(fml2)
      } else {
        lm(fml2, subset = eval(subset_expr))
      }
    }
  )
}

d1_test_ace_fast <- function(mira_full, fml_full, mids_obj, subset_expr = NULL) {
  fml_red  <- update(fml_full, . ~ . - Classic_ACEs_cat)
  mira_red <- fit_mi(mids_obj, fml_red, subset_expr)
  mice::D1(mira_full, mira_red)
}

d1_test_interaction_fast <- function(mira_full, fml_full, mids_obj) {
  fml_red  <- update(fml_full, . ~ . - Classic_ACEs_cat:Child_sex)
  mira_red <- fit_mi(mids_obj, fml_red)
  mice::D1(mira_full, mira_red)
}


# Section 6: Run plan ----

run_mi_d1_suite_clean <- function(mids_obj,
                                  outcomes = c("PWV", "Arterial_Dist", "cIMT"),
                                  ages     = c(17, 24),
                                  models   = c("M1", "M2"),
                                  groups   = c("Whole", "Female", "Male"),
                                  include_sex_interaction_M2 = FALSE) {

  plan <- tidyr::expand_grid(
    outcome = outcomes,
    age     = ages,
    model   = models,
    group   = groups
  )

  res <- lapply(seq_len(nrow(plan)), function(i) {

    outcome <- plan$outcome[i]
    age     <- plan$age[i]
    model   <- plan$model[i]
    group   <- plan$group[i]

    subset_expr <- switch(group,
                          Whole  = NULL,
                          Female = quote(Child_sex == "Female"),
                          Male   = quote(Child_sex == "Male"))

    fml <- build_fml(outcome, age, model, group, include_sex_interaction_M2)

    message("Fitting: ", outcome, "_", age, " ", model, " ", group)

    mira   <- fit_mi(mids_obj, fml, subset_expr)
    pooled <- summary(mice::pool(mira), conf.int = TRUE)

    d1_ace <- d1_test_ace_fast(mira, fml, mids_obj, subset_expr)
    d1_eta <- d1_to_eta2(d1_ace)

    d1_int <- NULL
    if (group == "Whole" && model == "M2" && include_sex_interaction_M2) {
      d1_int <- d1_test_interaction_fast(mira, fml, mids_obj)
    }

    sigma_pooled <- get_pooled_sigma(mira)

    emm_results <- NULL
    if (model == "M2" && group == "Whole") {
      emm_results <- compute_mi_pairwise_emmeans(
        mira,
        interaction   = include_sex_interaction_M2,
        sigma_pooled  = sigma_pooled,
        age           = age,
        outcome       = outcome   # <-- pass outcome in
      )
    }

    list(
      spec      = data.frame(outcome = outcome, age = age, model = model, group = group,
                             formula = deparse(fml), stringsAsFactors = FALSE),
      pooled    = pooled,
      D1_ACE    = d1_ace,
      D1_ETA    = d1_eta,
      D1_INT    = d1_int,
      emmeans   = emm_results$emmeans,
      contrasts = emm_results$contrasts
    )
  })

  names(res) <- with(plan, paste(outcome, age, model, group, sep = "_"))
  res
}


# Section 7: Print results to console ----

print_anova_results <- function(results_list, ace_only = FALSE) {
  for (nm in names(results_list)) {
    res <- results_list[[nm]]
    cat("\n\n====", nm, "====\n")

    if (!is.null(res$spec))   { cat("\n-- Specification --\n");    print(res$spec) }
    if (!is.null(res$pooled)) { cat("\n-- Pooled regression --\n"); print(as.data.frame(res$pooled)) }

    if (!is.null(res$D1_ACE)) {
      cat("\n-- D1 test: ACE effect --\n")
      print(as.data.frame(res$D1_ACE$result))
    }
    if (!is.null(res$D1_ETA)) {
      cat("\n-- Partial eta^2 --\n")
      print(as.data.frame(res$D1_ETA))
    }
    if (!is.null(res$D1_INT)) {
      cat("\n-- D1 test: ACE x Sex interaction --\n")
      print(as.data.frame(res$D1_INT$result))
    }
    if (!is.null(res$emmeans)) {
      cat("\n-- Pooled emmeans --\n");   print(as.data.frame(res$emmeans))
      cat("\n-- Pairwise contrasts --\n"); print(as.data.frame(res$contrasts))
    }
  }
}


# Section 8: Run models ----

# Set TRUE to include Classic_ACEs_cat * Child_sex interaction in whole-sample Model 2
include_sex_interaction_M2 <- TRUE

Results_D1_Imputed_17_PWV <- run_mi_d1_suite_clean(
  mids_obj   = Study1_Imp_Classic_Final_Transformed_V2,
  outcomes   = "PWV",
  ages       = 17,
  models     = c("M1", "M2"),
  groups     = c("Whole", "Female", "Male"),
  include_sex_interaction_M2 = include_sex_interaction_M2
)
saveRDS(Results_D1_Imputed_17_PWV, file = "Results_D1_Imputed_17_PWV.rds")
rm(Results_D1_Imputed_17_PWV)

Results_D1_Imputed_17_Dist <- run_mi_d1_suite_clean(
  mids_obj   = Study1_Imp_Classic_Final_Transformed_V2,
  outcomes   = "Arterial_Dist",
  ages       = 17,
  models     = c("M1", "M2"),
  groups     = c("Whole", "Female", "Male"),
  include_sex_interaction_M2 = include_sex_interaction_M2
)
saveRDS(Results_D1_Imputed_17_Dist, file = "Results_D1_Imputed_17_Dist.rds")
rm(Results_D1_Imputed_17_Dist)

Results_D1_Imputed_17_cIMT <- run_mi_d1_suite_clean(
  mids_obj   = Study1_Imp_Classic_Final_Transformed_V2,
  outcomes   = "cIMT",
  ages       = 17,
  models     = c("M1", "M2"),
  groups     = c("Whole", "Female", "Male"),
  include_sex_interaction_M2 = include_sex_interaction_M2
)
saveRDS(Results_D1_Imputed_17_cIMT, file = "Results_D1_Imputed_17_cIMT.rds")
rm(Results_D1_Imputed_17_cIMT)

Results_D1_Imputed_24_PWV <- run_mi_d1_suite_clean(
  mids_obj   = Study1_Imp_Classic_Final_Transformed_V2,
  outcomes   = "PWV",
  ages       = 24,
  models     = c("M1", "M2"),
  groups     = c("Whole", "Female", "Male"),
  include_sex_interaction_M2 = include_sex_interaction_M2
)
saveRDS(Results_D1_Imputed_24_PWV, file = "Results_D1_Imputed_24_PWV.rds")
rm(Results_D1_Imputed_24_PWV)

Results_D1_Imputed_24_cIMT <- run_mi_d1_suite_clean(
  mids_obj   = Study1_Imp_Classic_Final_Transformed_V2,
  outcomes   = "cIMT",
  ages       = 24,
  models     = c("M1", "M2"),
  groups     = c("Whole", "Female", "Male"),
  include_sex_interaction_M2 = include_sex_interaction_M2
)
saveRDS(Results_D1_Imputed_24_cIMT, file = "Results_D1_Imputed_24_cIMT.rds")
rm(Results_D1_Imputed_24_cIMT)


# Section 9: Print results ----

Results_D1_Imputed_17_PWV  <- readRDS("Results_D1_Imputed_17_PWV.rds")
print_anova_results(Results_D1_Imputed_17_PWV)
rm(Results_D1_Imputed_17_PWV)

Results_D1_Imputed_17_Dist <- readRDS("Results_D1_Imputed_17_Dist.rds")
print_anova_results(Results_D1_Imputed_17_Dist)
rm(Results_D1_Imputed_17_Dist)

Results_D1_Imputed_17_cIMT <- readRDS("Results_D1_Imputed_17_cIMT.rds")
print_anova_results(Results_D1_Imputed_17_cIMT)
rm(Results_D1_Imputed_17_cIMT)

Results_D1_Imputed_24_PWV  <- readRDS("Results_D1_Imputed_24_PWV.rds")
print_anova_results(Results_D1_Imputed_24_PWV)
rm(Results_D1_Imputed_24_PWV)

Results_D1_Imputed_24_cIMT <- readRDS("Results_D1_Imputed_24_cIMT.rds")
print_anova_results(Results_D1_Imputed_24_cIMT)
rm(Results_D1_Imputed_24_cIMT)
