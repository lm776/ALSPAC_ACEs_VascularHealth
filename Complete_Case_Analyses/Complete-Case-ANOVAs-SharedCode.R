# ACEs and Vascular Health — Complete-Case Cross-Sectional Group Comparisons (ANOVAs)
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
# Overview:
#   Outcomes: PWV, Arterial_Dist, cIMT at ages 17 and 24
#   Groups:   Whole sample, Female, Male
#   Exposure: Classic_ACEs_cat (categorical: 0, 1-3, 4+)
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
#   Nested F-tests:
#     - Overall ACE group effect: compare full vs reduced model (drop Classic_ACEs_cat)
#     - (Optional) Interaction: Classic_ACEs_cat * Child_sex in whole-sample Model 2
#
#   Input data objects:
#     cov_complete_classic_17  — complete-case dataset at age 17
#     cov_complete_classic_24  — complete-case dataset at age 24
#
# Author:  Laura Macro
# Date:    2026


# Section 1: Required packages ----

library(emmeans)
library(dplyr)
library(tidyr)


# Section 2: Create categorical ACEs variable ----

cov_complete_classic_17 <- cov_complete_classic_17 |>
  mutate(Classic_ACEs_cat = case_when(
    Classic_ACEs == 0                          ~ "0",
    Classic_ACEs >= 1 & Classic_ACEs <= 3      ~ "1-3",
    Classic_ACEs >= 4                          ~ "4+"
  )) |>
  mutate(Classic_ACEs_cat = factor(Classic_ACEs_cat, levels = c("0", "1-3", "4+")))

cov_complete_classic_24 <- cov_complete_classic_24 |>
  mutate(Classic_ACEs_cat = case_when(
    Classic_ACEs == 0                          ~ "0",
    Classic_ACEs >= 1 & Classic_ACEs <= 3      ~ "1-3",
    Classic_ACEs >= 4                          ~ "4+"
  )) |>
  mutate(Classic_ACEs_cat = factor(Classic_ACEs_cat, levels = c("0", "1-3", "4+")))


# Section 3: Helper functions ----

# Wave-specific variable name constructors
y_var   <- function(outcome, age) paste0(outcome, "_", age)            # e.g. PWV_17
age_var <- function(age)          paste0("Age_", age, "_clinic_years") # e.g. Age_17_clinic_years
bp_var  <- function(age)          paste0("BP_systolic_", age)          # e.g. BP_systolic_17
bmi_var <- function(age)          paste0("BMI_", age)                  # e.g. BMI_17

# Model 2 base covariates (non-wave-specific)
m2_covars_base <- c(
  "Child_ethnicity", "Townsend_sum", "Marital_status", "Parent_edu",
  "Mat_PND_gest", "Mat_age_delivery", "Birth_weight_kg", "Mat_preg_smoke", "Mat_preg_alc",
  "Family_CVD", "Age_PHV_c"
)

# Drop covariates with fewer than 2 observed levels in the analysis subset.
# This avoids rank-deficiency errors in lm() for sex-stratified models.
drop_single_level_terms <- function(data, fml, subset_expr = NULL, verbose = TRUE) {
  vars_in_fml <- all.vars(fml)

  dat <- data
  if (!is.null(subset_expr)) {
    keep <- eval(subset_expr, envir = dat)
    keep <- ifelse(is.na(keep), FALSE, keep)
    dat  <- dat[keep, , drop = FALSE]
  }

  dropped <- character(0)

  for (v in vars_in_fml) {
    if (!v %in% names(dat)) next
    x <- dat[[v]]
    if (is.character(x) || is.factor(x)) {
      n_lev <- length(unique(stats::na.omit(as.character(x))))
      if (n_lev < 2) dropped <- c(dropped, v)
    }
  }

  outcome <- all.vars(stats::terms(fml))[1]
  dropped <- setdiff(unique(dropped), outcome)

  fml_new <- fml
  if (length(dropped) > 0) {
    for (v in dropped) fml_new <- update(fml_new, paste(". ~ . -", v))
  }

  if (verbose && length(dropped) > 0) {
    message("Dropped due to <2 observed levels: ", paste(dropped, collapse = ", "))
  }

  list(formula = fml_new, dropped = dropped, data_used = dat)
}

# Fit a single complete-case lm, dropping any single-level terms first
fit_cc <- function(data, fml, subset_expr = NULL, verbose = TRUE) {
  chk  <- drop_single_level_terms(data, fml, subset_expr = subset_expr, verbose = verbose)
  mod  <- lm(chk$formula, data = chk$data_used)
  list(model = mod, formula = chk$formula, dropped = chk$dropped, n = nrow(chk$data_used))
}

# Nested F-test for the overall ACE group effect
cc_test_ace <- function(data, fml_full, subset_expr = NULL, verbose = TRUE) {
  term_labels <- attr(terms(fml_full), "term.labels")
  has_int     <- any(grepl("^Classic_ACEs_cat:Child_sex$", term_labels))

  fml_red <- if (has_int) {
    update(fml_full, . ~ . - Classic_ACEs_cat - Classic_ACEs_cat:Child_sex)
  } else {
    update(fml_full, . ~ . - Classic_ACEs_cat)
  }

  fit_full <- fit_cc(data, fml_full, subset_expr, verbose = verbose)
  fit_red  <- fit_cc(data, fml_red,  subset_expr, verbose = verbose)

  list(
    anova        = anova(fit_red$model, fit_full$model),
    dropped_full = fit_full$dropped,
    dropped_red  = fit_red$dropped,
    formula_full = fit_full$formula,
    formula_red  = fit_red$formula
  )
}

# Nested F-test for the ACE x Sex interaction
cc_test_interaction <- function(data, fml_full, subset_expr = NULL, verbose = TRUE) {
  fml_red  <- update(fml_full, . ~ . - Classic_ACEs_cat:Child_sex)
  fit_full <- fit_cc(data, fml_full, subset_expr, verbose = verbose)
  fit_red  <- fit_cc(data, fml_red,  subset_expr, verbose = verbose)

  list(
    anova        = anova(fit_red$model, fit_full$model),
    dropped_full = fit_full$dropped,
    dropped_red  = fit_red$dropped,
    formula_full = fit_full$formula,
    formula_red  = fit_red$formula
  )
}

# Partial eta^2 from a nested anova() result
anova_to_eta2 <- function(aov_obj) {
  df <- as.data.frame(aov_obj$anova)
  if (nrow(df) < 2L) return(data.frame(note = "No nested comparison row"))

  row    <- df[2, ]
  F_stat <- row$F
  df1    <- row$Df
  df2    <- row$Res.Df
  p_val  <- row$`Pr(>F)`
  eta2_p <- (F_stat * df1) / (F_stat * df1 + df2)

  data.frame(F_stat = F_stat, df1 = df1, df2 = df2, p_value = p_val, eta2_p = eta2_p)
}

# Compute emmeans and pairwise contrasts (including Cohen's d) for a single lm fit
compute_cc_pairwise_emmeans <- function(mod, interaction = FALSE, age) {
  bmi  <- bmi_var(age)
  bp   <- bp_var(age)
  agec <- age_var(age)

  nuisance_base <- c(
    "Marital_status", "Parent_edu",
    "Mat_PND_gest", "Mat_age_delivery", "Birth_weight_kg",
    "Mat_preg_smoke", "Mat_preg_alc", "Family_CVD", "Age_PHV_c"
  )
  nuisance_vars <- c(nuisance_base, bmi, agec, bp)

  emm <- suppressMessages(
    if (interaction) {
      emmeans::emmeans(mod, ~ Classic_ACEs_cat | Child_sex,
                       nuisance = nuisance_vars, rg.limit = 50000)
    } else {
      emmeans::emmeans(mod, ~ Classic_ACEs_cat,
                       nuisance = nuisance_vars, rg.limit = 50000)
    }
  )

  emm_df <- as.data.frame(emm) |>
    dplyr::rename(estimate = emmean, CI_low = lower.CL, CI_high = upper.CL)

  contr_df           <- as.data.frame(emmeans::contrast(emm, method = "pairwise", adjust = "none"))
  contr_df$Cohens_d  <- contr_df$estimate / summary(mod)$sigma

  list(emmeans = emm_df, contrasts = contr_df)
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
  ace  <- "Classic_ACEs_cat"

  # Model 1 RHS
  rhs_m1 <- if (outcome %in% c("PWV", "Arterial_Dist")) {
    paste(c(ace, agec, bp), collapse = " + ")
  } else {
    paste(c(ace, agec), collapse = " + ")
  }

  if (model == "M1") return(as.formula(paste0(y, " ~ ", rhs_m1)))

  # Model 2 RHS
  rhs_m2 <- c(m2_covars_base, bmi, agec)
  if (outcome %in% c("PWV", "Arterial_Dist")) rhs_m2 <- c(rhs_m2, bp)

  rhs <- if (group == "Whole") {
    if (include_sex_interaction_M2) {
      paste(c(paste0(ace, " * Child_sex"), rhs_m2), collapse = " + ")
    } else {
      paste(c(ace, "Child_sex", rhs_m2), collapse = " + ")
    }
  } else {
    paste(c(ace, rhs_m2), collapse = " + ")
  }

  as.formula(paste0(y, " ~ ", rhs))
}


# Section 5: Run plan ----

run_cc_anova_suite_clean <- function(data,
                                     outcomes = c("PWV", "Arterial_Dist", "cIMT"),
                                     ages     = c(17, 24),
                                     models   = c("M1", "M2"),
                                     groups   = c("Whole", "Female", "Male"),
                                     include_sex_interaction_M2 = FALSE,
                                     verbose  = TRUE) {

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

    message("Fitting CC: ", outcome, "_", age, " ", model, " ", group)

    fit_obj   <- fit_cc(data, fml, subset_expr, verbose = verbose)
    mod       <- fit_obj$model
    summ      <- summary(mod)

    anova_ace <- cc_test_ace(data, fml, subset_expr, verbose = verbose)
    eta_ace   <- anova_to_eta2(anova_ace)

    anova_int <- NULL
    eta_int   <- NULL
    if (group == "Whole" && include_sex_interaction_M2) {
      anova_int <- cc_test_interaction(data, fml, subset_expr, verbose = verbose)
      eta_int   <- anova_to_eta2(anova_int)
    }

    emm_results <- NULL
    if (model == "M2" && group == "Whole") {
      emm_results <- compute_cc_pairwise_emmeans(
        mod,
        interaction = include_sex_interaction_M2,
        age         = age
      )
    }

    list(
      spec = data.frame(
        outcome           = outcome,
        age               = age,
        model             = model,
        group             = group,
        formula_requested = paste(deparse(fml), collapse = " "),
        formula_fitted    = paste(deparse(stats::formula(mod)), collapse = " "),
        n                 = stats::nobs(mod),
        stringsAsFactors  = FALSE
      ),
      dropped_terms = fit_obj$dropped,
      model         = mod,
      summary       = summ,
      anova_ACE     = anova_ace,
      eta_ACE       = eta_ace,
      anova_INT     = anova_int,
      eta_INT       = eta_int,
      emmeans       = if (!is.null(emm_results)) emm_results$emmeans   else NULL,
      contrasts     = if (!is.null(emm_results)) emm_results$contrasts else NULL
    )
  })

  names(res) <- with(plan, paste(outcome, age, model, group, sep = "_"))
  res
}


# Section 6: Run models ----

include_sex_interaction_M2 <- TRUE

Results_CC_17_PWV <- run_cc_anova_suite_clean(
  data     = cov_complete_classic_17,
  outcomes = "PWV",
  ages     = 17,
  models   = c("M1", "M2"),
  groups   = c("Whole", "Female", "Male"),
  include_sex_interaction_M2 = include_sex_interaction_M2,
  verbose  = TRUE
)
saveRDS(Results_CC_17_PWV, file = "Results_CC_17_PWV.rds")

Results_CC_17_Arterial_Dist <- run_cc_anova_suite_clean(
  data     = cov_complete_classic_17,
  outcomes = "Arterial_Dist",
  ages     = 17,
  models   = c("M1", "M2"),
  groups   = c("Whole", "Female", "Male"),
  include_sex_interaction_M2 = include_sex_interaction_M2,
  verbose  = TRUE
)
saveRDS(Results_CC_17_Arterial_Dist, file = "Results_CC_17_Arterial_Dist.rds")

Results_CC_17_cIMT <- run_cc_anova_suite_clean(
  data     = cov_complete_classic_17,
  outcomes = "cIMT",
  ages     = 17,
  models   = c("M1", "M2"),
  groups   = c("Whole", "Female", "Male"),
  include_sex_interaction_M2 = include_sex_interaction_M2,
  verbose  = TRUE
)
saveRDS(Results_CC_17_cIMT, file = "Results_CC_17_cIMT.rds")

Results_CC_24_PWV <- run_cc_anova_suite_clean(
  data     = cov_complete_classic_24,
  outcomes = "PWV",
  ages     = 24,
  models   = c("M1", "M2"),
  groups   = c("Whole", "Female", "Male"),
  include_sex_interaction_M2 = include_sex_interaction_M2,
  verbose  = TRUE
)
saveRDS(Results_CC_24_PWV, file = "Results_CC_24_PWV.rds")

Results_CC_24_cIMT <- run_cc_anova_suite_clean(
  data     = cov_complete_classic_24,
  outcomes = "cIMT",
  ages     = 24,
  models   = c("M1", "M2"),
  groups   = c("Whole", "Female", "Male"),
  include_sex_interaction_M2 = include_sex_interaction_M2,
  verbose  = TRUE
)
saveRDS(Results_CC_24_cIMT, file = "Results_CC_24_cIMT.rds")


# Section 7: Print results to console ----

print_cc_results <- function(results_list) {
  for (nm in names(results_list)) {
    res <- results_list[[nm]]
    cat("\n\n====", nm, "====\n")

    if (!is.null(res$spec)) {
      cat("\n-- Specification --\n")
      print(res$spec)
    }
    if (length(res$dropped_terms) > 0) {
      cat("\n-- Dropped terms --\n")
      print(res$dropped_terms)
    }
    if (!is.null(res$summary)) {
      cat("\n-- Model summary --\n")
      print(res$summary)
    }
    if (!is.null(res$anova_ACE)) {
      cat("\n-- Nested F-test: ACE effect --\n")
      print(res$anova_ACE$anova)
      cat("\n-- Partial eta^2: ACE effect --\n")
      print(as.data.frame(res$eta_ACE))
    }
    if (!is.null(res$anova_INT)) {
      cat("\n-- Nested F-test: ACE x Sex interaction --\n")
      print(res$anova_INT$anova)
      cat("\n-- Partial eta^2: ACE x Sex interaction --\n")
      print(as.data.frame(res$eta_INT))
    }
    if (!is.null(res$emmeans)) {
      cat("\n-- Estimated marginal means --\n")
      print(as.data.frame(res$emmeans))
      cat("\n-- Pairwise contrasts --\n")
      print(as.data.frame(res$contrasts))
    }
  }
}

print_cc_results(Results_CC_17_PWV)
print_cc_results(Results_CC_17_Arterial_Dist)
print_cc_results(Results_CC_17_cIMT)
print_cc_results(Results_CC_24_PWV)
print_cc_results(Results_CC_24_cIMT)


# Section 8: Sensitivity table (for supplementary materials)----

# Extracts F, df, p, and partial eta^2 for the ACE main effect and ACE x Sex
# interaction for each whole-sample model, and combines into a wide comparison
# table (M1 vs M2 side by side).

extract_sensitivity_row <- function(res, effect = c("ACE", "ACExSex")) {
  effect <- match.arg(effect)

  spec  <- as.data.frame(res$spec)
  outc  <- spec$outcome[1]
  age   <- spec$age[1]
  model <- spec$model[1]
  group <- spec$group[1]
  n     <- spec$n[1]

  if (effect == "ACE") {
    aov_obj <- res$anova_ACE
    eta_obj <- res$eta_ACE
    eff_lab <- "ACE main effect"
  } else {
    aov_obj <- res$anova_INT
    eta_obj <- res$eta_INT
    eff_lab <- "ACE x Sex interaction"
  }

  na_row <- data.frame(
    Outcome = outc, Age = age, Model = model, Group = group, n = n,
    Effect = eff_lab, df1 = NA_real_, df2 = NA_real_,
    F = NA_real_, p_value = NA_real_, eta2_p = NA_real_,
    stringsAsFactors = FALSE
  )

  if (is.null(aov_obj)) return(na_row)

  aov_tab <- if (!is.null(aov_obj$anova)) aov_obj$anova else aov_obj
  df_aov  <- as.data.frame(aov_tab)
  if (nrow(df_aov) < 2) return(na_row)

  row2    <- df_aov[2, , drop = FALSE]
  eta_val <- NA_real_
  if (!is.null(eta_obj)) {
    eta_df  <- as.data.frame(eta_obj)
    if ("eta2_p" %in% names(eta_df)) eta_val <- eta_df$eta2_p[1]
  }

  data.frame(
    Outcome = outc, Age = age, Model = model, Group = group, n = n,
    Effect  = eff_lab,
    df1     = row2$Df,
    df2     = row2$Res.Df,
    F       = row2$F,
    p_value = row2$`Pr(>F)`,
    eta2_p  = eta_val,
    stringsAsFactors = FALSE
  )
}

filter_whole <- function(results_list) {
  keep <- vapply(results_list, function(res) {
    as.data.frame(res$spec)$group[1] == "Whole"
  }, logical(1))
  results_list[keep]
}

build_sensitivity_for_results <- function(results_list) {
  results_whole <- filter_whole(results_list)
  do.call(rbind, lapply(results_whole, function(res) {
    rbind(
      extract_sensitivity_row(res, effect = "ACE"),
      extract_sensitivity_row(res, effect = "ACExSex")
    )
  }))
}

sens_17_PWV  <- build_sensitivity_for_results(Results_CC_17_PWV)
sens_17_Dist <- build_sensitivity_for_results(Results_CC_17_Arterial_Dist)
sens_17_cIMT <- build_sensitivity_for_results(Results_CC_17_cIMT)
sens_24_PWV  <- build_sensitivity_for_results(Results_CC_24_PWV)
sens_24_cIMT <- build_sensitivity_for_results(Results_CC_24_cIMT)

sensitivity_all <- rbind(
  sens_17_PWV,
  sens_17_Dist,
  sens_17_cIMT,
  sens_24_PWV,
  sens_24_cIMT
)

# Reshape to wide: M1 and M2 side by side
sens_compare <- sensitivity_all |>
  mutate(
    F       = round(F, 2),
    eta2_p  = round(eta2_p, 3),
    p_value = case_when(
      is.na(p_value)  ~ NA_character_,
      p_value < 0.001 ~ "<0.001",
      TRUE            ~ formatC(p_value, digits = 3, format = "f")
    ),
    df1 = as.character(df1),
    df2 = as.character(df2)
  ) |>
  select(Outcome, Age, Model, Effect, df1, df2, F, p_value, eta2_p) |>
  pivot_wider(
    names_from  = Model,
    values_from = c(df1, df2, F, p_value, eta2_p),
    names_glue  = "{Model}_{.value}"
  ) |>
  arrange(Outcome, Age, Effect)

# Print supplementary sensitivity table to console
cat("\n\n==== Supplementary sensitivity ANOVA table ====\n")
print(as.data.frame(sens_compare))
