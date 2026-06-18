# ALSPAC_ACEs_VascularHealth

This repository shares the analysis code used in the study:

> **"Adverse Childhood Experiences and Changes in Vascular Health from Childhood to Mid-Adulthood: Cross-Sectional and Longitudinal Evidence from the ALSPAC Study"**
> (in preparation - authors hidden for blinded review purposes)

The study was pre-registered via OSF: [https://osf.io/6gxn8/overview?view_only=1e2e34729194476cb8c74380d653e696](https://osf.io/6gxn8/overview?view_only=1e2e34729194476cb8c74380d653e696)

This code repository is archived on Zenodo: [https://doi.org/10.5281/zenodo.20744231](https://doi.org/10.5281/zenodo.20744231)

---

## Study Aims

1. Assess whether ACEs are associated with pulse wave velocity (PWV), arterial distensibility, and carotid intima-media thickness (cIMT) measured in adolescence and early adulthood.
2. Determine whether ACE exposure is associated with the level of change in these vascular markers between adolescence and early adulthood.
3. Explore whether relationships between ACEs and vascular health differ by sex, and after adjusting for sociodemographic variables, child health risk behaviours, and physiological biomarkers.

---

## Data Access

The data used in these analyses were obtained from the **Avon Longitudinal Study of Parents and Children (ALSPAC)**, a UK birth cohort study. The underlying data cannot be shared as they are held under a data access agreement with the ALSPAC Executive Committee.

Researchers wishing to use ALSPAC data may apply via the ALSPAC data access process:
[https://www.bristol.ac.uk/alspac/researchers/access/](https://www.bristol.ac.uk/alspac/researchers/access/)

ALSPAC variable names used in this code are publicly documented via the ALSPAC variable search tool:
[https://variables.alspac.bris.ac.uk/](https://variables.alspac.bris.ac.uk/)

Please cite ALSPAC as:
> Boyd A, et al. (2012). Cohort Profile: the 'children of the 90s' — the index offspring of the Avon Longitudinal Study of Parents and Children. *International Journal of Epidemiology*, 42(1):111–127. https://doi.org/10.1093/ije/dys064
>
> Fraser A, et al. (2012). Cohort Profile: the Avon Longitudinal Study of Parents and Children: ALSPAC mothers cohort. *International Journal of Epidemiology*, 42(1):97–110. https://doi.org/10.1093/ije/dys066

---

## Ethical Approval

Ethical approval for the ALSPAC study was obtained from the ALSPAC Ethics and Law Committee and the Local Research Ethics Committees. Ethical approval for this analysis was obtained in accordance with the terms of the ALSPAC data access agreement.

---

## Repository Structure

Scripts are numbered in the order they should be run. Complete-case analyses (scripts 07–09) are sensitivity analyses and can be run independently of the multiple imputation pipeline.

```
ALSPAC_ACEs_VascularHealth/
│
├── 01_Preliminary_Analysis_PublicCode.R
│       Data exploration and preliminary analysis. Descriptive statistics,
│       missing data summaries, and covariate distributions.
│
├── 02_MI-Analysis_Setup_SharedCode.R
│       Multiple imputation setup using the mice package. Includes predictor
│       matrix specification, imputation model selection, and generation of
│       the final imputed datasets for Classic and Extended ACE exposures.
│
├── 03_Imputed-Analyses-Participant-Characteristics-SharedCode.R
│       Derives pooled participant characteristics tables (Table 1) from
│       multiply-imputed data using Rubin's rules. Produces whole-sample
│       and sex-stratified tables for both ACE exposures.
│
├── 04_Imputed-Analyses-ANOVAs-SharedCode.R
│       Cross-sectional group comparisons across ACE categories using D1
│       pooled F-tests (mice). Includes estimated marginal means, pairwise
│       contrasts, and partial eta-squared for PWV, arterial distensibility,
│       and cIMT at ages 17 and 24.
│
├── 05_Imputed-Analyses-Regressions-Classic-SharedCode.R
│       Pooled linear regression analyses (Models 1–5) for Classic ACEs
│       (continuous and categorical) across vascular outcomes at ages 17
│       and 24. Includes adaptation notes at the top of the file explaining
│       how to substitute Classic ACEs variables and objects for Extended ACEs.
│
├── 06_Imputed-Analyses-Mixed-Effects-Classic-SharedCode.R
│       Pooled mixed-effects trajectory models (lme4, fitted across
│       imputations) examining change in vascular outcomes between ages
│       17 and 24, for Classic and Extended ACE exposures.
│
├── 07_Complete-Case-ANOVAs-SharedCode.R
│       Complete-case sensitivity analyses mirroring script 04. Includes
│       nested F-tests, emmeans, pairwise contrasts, and a supplementary
│       sensitivity comparison table (M1 vs M2).
│
├── 08_Complete-Case-Regression-Analyses-SharedCode.R
│       Complete-case sensitivity analyses mirroring script 05.
│       Linear regression models (Models 1–5) for Classic ACEs.
│
├── 09_Complete-Case-Mixed-Effects-Models-SharedCode.R
│       Complete-case sensitivity analyses mirroring script 06.
│       Mixed-effects trajectory models for Classic ACEs.
│
└── README.md
```

---

## Adapting the Code for Extended ACEs

Scripts 05 and 06 cover both Classic and Extended ACE exposures. The Classic ACEs scripts contain an **adaptation note at the top of the file** detailing the specific variable names, mids object names, and helper functions to substitute when running the Extended ACEs analyses. The key substitutions are:

| Classic ACEs | Extended ACEs |
|---|---|
| `Classic_ACEs` | `Extended_ACEs` |
| `Classic_ACEs_cat` | `Ext_ACEs_cat` |
| `Study1_Imp_Classic_Final_Transformed_V2` | `Study1_Imp_Ext_Trns_Reduced` |
| `pick_forms()` | `pick_forms_ext()` |

---

## Software Requirements

All analyses were conducted in **R** (version 4.2.0 or later). The following packages are required:

### Core analysis
| Package | Purpose |
|---|---|
| `mice` | Multiple imputation and pooling |
| `miceadds` | `datalist2mids()` for combining imputation runs |
| `lme4` | Linear mixed-effects models |
| `lmerTest` | Satterthwaite p-values for mixed-effects models |
| `emmeans` | Estimated marginal means and pairwise contrasts |
| `car` | Added-variable plots (`avPlots()`) |
| `performance` | Model assumption checks (`check_model()`) |
| `see` | Required by `performance` for assumption plots |
| `parameters` | Model parameter summaries |
| `broom` | Tidying linear model output |
| `broom.mixed` | Tidying mixed-effects model output |

### Data wrangling and visualisation
| Package | Purpose |
|---|---|
| `dplyr` | Data wrangling |
| `tidyr` | Data reshaping |
| `forcats` | Factor handling |
| `stringr` | String manipulation |
| `ggplot2` | Data visualisation |
| `ggmice` | Visualising imputation diagnostics |
| `psych` | Descriptive statistics (`describe()`) |

### Data import and labelling
| Package | Purpose |
|---|---|
| `haven` | Reading SPSS/Stata data files |
| `labelled` | Handling labelled variable data |
| `lmtest` | Coefficient tests |
| `sandwich` | Robust standard errors |

Install all required packages with:

```r
install.packages(c(
  # Core analysis
  "mice", "miceadds", "lme4", "lmerTest", "emmeans",
  "car", "performance", "see", "parameters", "broom", "broom.mixed",
  # Data wrangling and visualisation
  "dplyr", "tidyr", "forcats", "stringr", "ggplot2", "ggmice", "psych",
  # Data import and labelling
  "haven", "labelled", "lmtest", "sandwich"
))
```

---

## Important Notes for Code Users

- **No data are included in this repository.** Scripts will not run without access to the ALSPAC dataset obtained via the formal application process described above.
- All scripts print results to the R console. No output files are written — this was a deliberate choice to avoid sharing derived data that could indirectly identify participants.
- Scripts assume that the input data object is named `data_transformed` (or the relevant mids object as specified in each script header). Users will need to rename their own data objects accordingly, or update the object names at the top of each script.
- Raw ALSPAC variable names are used throughout and are publicly documented via the ALSPAC variable search tool linked above.

---

## Citation

The code itself can be cited directly via Zenodo:

> ALSPAC ACEs and Vascular Health Analysis (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.20744231

Please also cite ALSPAC as detailed in the Data Access section above.
