# =============================================================================
# Risk Factors for Diabetes and Hypertension (Non-Communicable Diseases)
# Cross-sectional survey, Tete, Mozambique (2023)
# =============================================================================
#
# Goal:     Estimate the prevalence of diabetes, pre-diabetes, hypertension and
#           overweight/obesity, and identify associated factors with Fisher's
#           exact tests and logistic regression (crude and adjusted odds ratios).
#
# Input:    A semicolon-separated CSV with one row per participant
#           (see README for the columns). Data are NOT included in this repository.
#
# Sections:
#   1. Setup (packages, paths, helper functions)
#   2. Data import
#   3. Variable recoding
#   4. Descriptive statistics and prevalence
#   5. Bivariate analysis (Fisher's exact test)
#   6. Crude odds ratios
#   7. Multivariable model: diabetes
#   8. Multivariable model: hypertension
#   9. Correlation (age and height)
#
# Run the script top to bottom, section by section, in RStudio.
# =============================================================================


# =============================================================================
# 1. SETUP ----
# =============================================================================

# --- 1.1 Packages ------------------------------------------------------------
pkgs <- c("tidyverse", "broom", "knitr", "pROC", "ResourceSelection", "car")
missing_pkgs <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(tidyverse)          # data wrangling and ggplot2
  library(broom)              # tidy model output
  library(knitr)              # kable() tables
  library(pROC)               # ROC curves and AUC
  library(ResourceSelection)  # Hosmer-Lemeshow goodness-of-fit test
  library(car)                # vif()
})

# --- 1.2 Paths ---------------------------------------------------------------
# Edit these lines for your machine. Do NOT commit the data file.
DATA_FILE <- "data/Dataset_DCNT_Tete.csv"
OUT_DIR   <- "output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- 1.3 Helper functions ----------------------------------------------------

# Prevalence (%) with exact binomial 95% CI for each category of a variable.
prevalence_table <- function(x, label = deparse(substitute(x))) {
  x <- x[!is.na(x)]
  N <- length(x)
  map_dfr(sort(unique(x)), function(level) {
    n  <- sum(x == level)
    ci <- binom.test(n, N)$conf.int * 100
    tibble(variable = label, category = as.character(level), n = n, N = N,
           prevalence = round(n / N * 100, 2),
           ci_low = round(ci[1], 2), ci_high = round(ci[2], 2))
  })
}

# Cross-tabulation of an exposure by a binary outcome, with row percentages
# and Fisher's exact test. Returns the p-value.
crosstab_fisher <- function(data, exposure, outcome) {
  tab <- table(data[[exposure]], data[[outcome]], dnn = c(exposure, outcome))
  cat("\n=====", exposure, "x", outcome, "=====\n")
  print(tab)
  print(round(prop.table(tab, margin = 1) * 100, 1))   # row percentages
  p <- fisher.test(tab)$p.value
  cat("Fisher exact p-value:", signif(p, 3), "\n")
  tibble(exposure = exposure, outcome = outcome, p_fisher = signif(p, 3))
}

# Odds ratios with Wald 95% CIs from a logistic model (intercept dropped).
or_table <- function(fit) {
  co <- summary(fit)$coefficients
  z  <- qnorm(0.975)
  tibble(term      = rownames(co),
         OR        = exp(co[, "Estimate"]),
         conf.low  = exp(co[, "Estimate"] - z * co[, "Std. Error"]),
         conf.high = exp(co[, "Estimate"] + z * co[, "Std. Error"]),
         p.value   = co[, ncol(co)]) |>
    filter(term != "(Intercept)") |>
    mutate(across(where(is.numeric), ~ round(.x, 3)))
}

# Crude (unadjusted) odds ratios: one logistic model per exposure.
crude_or <- function(data, exposures, outcome) {
  map_dfr(exposures, function(v) {
    fit <- glm(as.formula(paste(outcome, "~", v)), data = data, family = binomial)
    or_table(fit) |> mutate(exposure = v, .before = 1)
  })
}

# Forest plot of odds ratios (log scale).
plot_or <- function(fit, title) {
  or_table(fit) |>
    ggplot(aes(x = OR, y = reorder(term, OR), xmin = conf.low, xmax = conf.high)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "red") +
    geom_pointrange() +
    scale_x_log10() +
    labs(title = title, x = "Odds ratio (log scale, 95% CI)", y = NULL) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
}

# Model performance: n, AIC, AUC (95% CI), Hosmer-Lemeshow test, and VIF.
# Uses the model's own data, so rows with missing values are handled correctly.
evaluate_model <- function(fit, label) {
  roc_obj <- roc(fit$y, fitted(fit), quiet = TRUE)
  hl      <- hoslem.test(fit$y, fitted(fit), g = 10)

  plot(roc_obj, main = paste("ROC curve:", label), col = "blue", lwd = 2)

  # With factor predictors vif() returns a GVIF matrix; use GVIF^(1/(2*Df)).
  v <- vif(fit)
  vif_mean <- if (is.matrix(v)) mean(v[, "GVIF^(1/(2*Df))"]) else mean(v)

  tibble(model = label, n = length(fit$y), AIC = round(AIC(fit), 1),
         AUC = round(as.numeric(auc(roc_obj)), 3),
         AUC_low = round(ci.auc(roc_obj)[1], 3), AUC_high = round(ci.auc(roc_obj)[3], 3),
         HL_p = round(hl$p.value, 3), mean_adj_GVIF = round(vif_mean, 2))
}


# =============================================================================
# 2. DATA IMPORT ----
# =============================================================================

# read.csv2() reads semicolon-separated files (decimal comma)
Banco <- read.csv2(DATA_FILE)
names(Banco)
cat("Participants:", nrow(Banco), "\n")


# =============================================================================
# 3. VARIABLE RECODING ----
# =============================================================================
# Raw questionnaire items used:
#   Pergunta_2  = takes medication for hypertension (1 yes, 2 no, 3 sometimes)
#   Pergunta_3  = family history of diabetes (1 yes, 2 no, 3 do not know)
#   Pergunta_5  = family history of hypertension (1 yes, 2 no, 3 do not know)
#   Pergunta_9  = time spent on physical activity (1: 5-10 min, 2: 10-20, 3: 20-30, 4: 30-40)
#   Pergunta_12 = sugar intake (2 = <= 5 tablespoons, 3 = > 5 tablespoons)
#   Pergunta_14 = salt intake (1 = < 1 tablespoon, 2 = 1 tablespoon, 3 = > 1 tablespoon)
#   Pergunta_16 = fruit intake frequency (0 = never)
#   Pergunta_17 = vegetable intake frequency (1 = 1-2x, 2 = 3-4x, 3-4 = > 4x)
#   Diagnostigo_do_diabetes = 2 normal, 3 pre-diabetic, 4 diabetic
#   Sexo = 1 male, 2 female

Banco <- Banco |>
  mutate(
    # --- Age categories ---
    Idade_cat = cut(Idade, breaks = c(-Inf, 30, 40, Inf),
                    labels = c("< 30", "30-40", "> 40"), include.lowest = TRUE),

    # --- Body mass index (kg/m2; height in meters) ---
    IMC_1 = Peso / Altura^2,

    # --- BMI category (WHO) from the BMI column in the dataset ---
    IMC_cat = cut(IMC, breaks = c(-Inf, 18.5, 25, 30, Inf), right = FALSE,
                  labels = c("Underweight", "Normal weight", "Overweight", "Obese")),

    # --- Sex ---
    Sexo = factor(Sexo, levels = c(1, 2), labels = c("Male", "Female")),

    # --- Diabetes: status (3 levels) and binary outcome (diabetic vs. not) ---
    Diabetes_status = factor(Diagnostigo_do_diabetes, levels = c(2, 3, 4),
                             labels = c("Normal", "Pre-diabetic", "Diabetic")),
    Diabetes_cat = if_else(Diagnostigo_do_diabetes == 4, 1, 0),
    Diabetes_f   = factor(Diabetes_cat, levels = c(0, 1), labels = c("No", "Yes")),

    # --- Hypertension: binary outcome ---
    Hipertensao_cat = if_else(Hipertensao == 1, 1, 0),
    Hipertensao_f   = factor(Hipertensao_cat, levels = c(0, 1), labels = c("No", "Yes")),

    # --- Control of hypertension (non-hypertensive/missing -> "Not_hypertensive") ---
    Controla_hipertensao = coalesce(
      case_when(Pergunta_2 == 1 ~ "Yes",
                Pergunta_2 == 2 ~ "No",
                Pergunta_2 == 3 ~ "Sometimes"),
      "Not_hypertensive"),

    # --- Family history of diabetes / hypertension ---
    Historico_fam_diabetes = factor(Pergunta_3, levels = c(2, 1, 3),
                                    labels = c("No", "Yes", "Do not know")),
    Historico_fam_htn      = factor(Pergunta_5, levels = c(2, 1, 3),
                                    labels = c("No", "Yes", "Do not know")),

    # --- Smoking ---
    Tabagismo = factor(Tabagismo, levels = c(2, 1), labels = c("No", "Yes")),

    # --- Physical activity (code 8 = none) ---
    Ativ_fisica = factor(if_else(Atividade_fisica == 8, "No", "Yes"),
                         levels = c("No", "Yes")),

    # --- Minutes of physical activity per session (missing = none) ---
    tempo_prat_ativ = factor(
      case_when(Pergunta_9 %in% c(1, 2) ~ "<= 20",
                Pergunta_9 %in% c(3, 4) ~ "> 20",
                is.na(Pergunta_9)       ~ "None"),
      levels = c("> 20", "<= 20", "None")),

    # --- Sugar and salt intake ---
    Tablespoons_sugar = factor(Pergunta_12, levels = c(2, 3),
                               labels = c("<= 5 tablespoons", "> 5 tablespoons")),
    Tablespoons_salt  = factor(Pergunta_14, levels = c(1, 2, 3),
                               labels = c("< 1 tablespoon", "1 tablespoon", "> 1 tablespoon")),

    # --- Fruit and vegetable intake ---
    Eat_fuits = factor(if_else(Pergunta_16 == 0, "No", "Yes"), levels = c("Yes", "No")),
    Eat_Vegetables = factor(
      case_when(Pergunta_17 == 1            ~ "1-2 times",
                Pergunta_17 == 2            ~ "3-4 times",
                Pergunta_17 %in% c(3, 4)    ~ "> 4 times"),
      levels = c("1-2 times", "3-4 times", "> 4 times"))
  )

# Reference categories: the first level of each factor (e.g. "< 30", "Underweight"
# for IMC_cat is overridden below to "Normal weight").
Banco$IMC_cat <- relevel(Banco$IMC_cat, ref = "Normal weight")

# Quick checks of the recoded variables
recoded_vars <- c("Idade_cat", "IMC_cat", "Sexo", "Diabetes_status", "Diabetes_cat",
                  "Hipertensao_cat", "Controla_hipertensao", "Historico_fam_diabetes",
                  "Historico_fam_htn", "Tabagismo", "Ativ_fisica", "tempo_prat_ativ",
                  "Tablespoons_sugar", "Tablespoons_salt", "Eat_fuits", "Eat_Vegetables")
lapply(Banco[recoded_vars], table, useNA = "ifany")

# Variables used in the analysis
#   Continuous : Idade, Peso, Altura, Glicemia, IMC_1
#   Categorical: Idade_cat, Sexo, IMC_cat, Ativ_fisica, tempo_prat_ativ, Tabagismo,
#                Historico_fam_htn, Historico_fam_diabetes, Eat_fuits, Eat_Vegetables,
#                Tablespoons_sugar, Tablespoons_salt
#   Outcomes   : Diabetes_cat, Hipertensao_cat

# Exposures used in the bivariate and crude analyses
FISHER_VARS <- c("Idade_cat", "Sexo", "IMC_cat", "Ativ_fisica", "tempo_prat_ativ",
                 "Tabagismo", "Historico_fam_htn", "Historico_fam_diabetes",
                 "Eat_fuits", "Eat_Vegetables", "Tablespoons_salt", "Tablespoons_sugar")

# Physical activity (yes/no) is left out of the crude ORs, as in the original analysis.
CRUDE_VARS_DM  <- c("Idade_cat", "Sexo", "IMC_cat", "tempo_prat_ativ", "Tabagismo",
                    "Hipertensao_f", "Historico_fam_htn", "Historico_fam_diabetes",
                    "Eat_fuits", "Eat_Vegetables", "Tablespoons_sugar", "Tablespoons_salt")
CRUDE_VARS_HTN <- c("Idade_cat", "Sexo", "IMC_cat", "tempo_prat_ativ", "Tabagismo",
                    "Historico_fam_htn", "Diabetes_f", "Eat_fuits", "Eat_Vegetables",
                    "Tablespoons_sugar", "Tablespoons_salt")


# =============================================================================
# 4. DESCRIPTIVE STATISTICS AND PREVALENCE ----
# =============================================================================

# --- 4.1 BMI (continuous) ----------------------------------------------------
summary(Banco$IMC_1)
c(mean = mean(Banco$IMC_1, na.rm = TRUE), sd = sd(Banco$IMC_1, na.rm = TRUE))
shapiro.test(Banco$IMC_1)   # normality test

# --- 4.2 Prevalence with exact binomial 95% CIs ------------------------------
prevalence <- bind_rows(
  prevalence_table(Banco$Diabetes_status, "Diabetes status"),
  prevalence_table(factor(Banco$Hipertensao_f), "Hypertension"),
  prevalence_table(Banco$IMC_cat, "BMI category")
)
prevalence
write_csv(prevalence, file.path(OUT_DIR, "prevalence.csv"))


# =============================================================================
# 5. BIVARIATE ANALYSIS (FISHER'S EXACT TEST) ----
# =============================================================================

# --- 5.1 Diabetes ------------------------------------------------------------
fisher_dm <- map_dfr(FISHER_VARS, ~ crosstab_fisher(Banco, .x, "Diabetes_cat"))

# --- 5.2 Hypertension --------------------------------------------------------
fisher_htn <- map_dfr(FISHER_VARS, ~ crosstab_fisher(Banco, .x, "Hipertensao_cat"))

write_csv(bind_rows(fisher_dm, fisher_htn), file.path(OUT_DIR, "fisher_tests.csv"))


# =============================================================================
# 6. CRUDE ODDS RATIOS ----
# =============================================================================

# --- 6.1 Diabetes as outcome -------------------------------------------------
crude_dm <- crude_or(Banco, CRUDE_VARS_DM, "Diabetes_cat")
print(crude_dm, n = Inf)
write_csv(crude_dm, file.path(OUT_DIR, "crude_OR_diabetes.csv"))

# --- 6.2 Hypertension as outcome ---------------------------------------------
crude_htn <- crude_or(Banco, CRUDE_VARS_HTN, "Hipertensao_cat")
print(crude_htn, n = Inf)
write_csv(crude_htn, file.path(OUT_DIR, "crude_OR_hypertension.csv"))


# =============================================================================
# 7. MULTIVARIABLE MODEL: DIABETES ----
# =============================================================================

# --- 7.1 Association between categorical explanatory variables ---------------
# Proportion of each category by physical-activity time (visual collinearity check).
Banco |>
  select(Idade_cat, Sexo, IMC_cat, tempo_prat_ativ, Tabagismo, Hipertensao_f,
         Historico_fam_htn, Historico_fam_diabetes, Eat_fuits, Eat_Vegetables,
         Tablespoons_salt, Tablespoons_sugar) |>
  mutate(across(everything(), as.character)) |>
  pivot_longer(-tempo_prat_ativ) |>
  ggplot(aes(value, fill = tempo_prat_ativ)) +
  geom_bar(position = "fill") +
  ylab("proportion") +
  facet_wrap(~ name, scales = "free", ncol = 2) +
  coord_flip()

# --- 7.2 Full model ----------------------------------------------------------
fit_dm <- glm(
  Diabetes_cat ~ Idade_cat + Sexo + IMC_cat + tempo_prat_ativ + Tabagismo +
    Hipertensao_f + Historico_fam_htn + Historico_fam_diabetes +
    Eat_fuits + Eat_Vegetables + Tablespoons_salt + Tablespoons_sugar,
  data = Banco, family = binomial
)
summary(fit_dm)
kable(or_table(fit_dm), format = "markdown", align = "c")

# Performance (ROC/AUC, Hosmer-Lemeshow) and multicollinearity (VIF)
vif(fit_dm)
eval_dm <- evaluate_model(fit_dm, "Diabetes - full model")

# Forest plot of adjusted odds ratios
p_or_dm <- plot_or(fit_dm, "Adjusted odds ratios: diabetes")
p_or_dm
ggsave(file.path(OUT_DIR, "forest_OR_diabetes.png"), p_or_dm,
       width = 9, height = 7, dpi = 300, bg = "white")

# --- 7.3 Reduced model -------------------------------------------------------
fit_dm_red <- glm(
  Diabetes_cat ~ Idade_cat + IMC_cat + Tabagismo + Historico_fam_htn +
    Historico_fam_diabetes + Eat_fuits,
  data = Banco, family = binomial
)
summary(fit_dm_red)
kable(or_table(fit_dm_red), format = "markdown", align = "c")
vif(fit_dm_red)
eval_dm_red <- evaluate_model(fit_dm_red, "Diabetes - reduced model")


# =============================================================================
# 8. MULTIVARIABLE MODEL: HYPERTENSION ----
# =============================================================================

# --- 8.1 Full model ----------------------------------------------------------
fit_htn <- glm(
  Hipertensao_cat ~ Idade_cat + Sexo + IMC_cat + tempo_prat_ativ + Tabagismo +
    Diabetes_f + Historico_fam_htn + Eat_fuits + Eat_Vegetables +
    Tablespoons_salt + Tablespoons_sugar,
  data = Banco, family = binomial
)
summary(fit_htn)
kable(or_table(fit_htn), format = "markdown", align = "c")

vif(fit_htn)
eval_htn <- evaluate_model(fit_htn, "Hypertension - full model")

p_or_htn <- plot_or(fit_htn, "Adjusted odds ratios: hypertension")
p_or_htn
ggsave(file.path(OUT_DIR, "forest_OR_hypertension.png"), p_or_htn,
       width = 9, height = 7, dpi = 300, bg = "white")

# --- 8.2 Reduced model -------------------------------------------------------
fit_htn_red <- glm(
  Hipertensao_cat ~ Sexo + Tabagismo + Eat_Vegetables,
  data = Banco, family = binomial
)
summary(fit_htn_red)
kable(or_table(fit_htn_red), format = "markdown", align = "c")
vif(fit_htn_red)
eval_htn_red <- evaluate_model(fit_htn_red, "Hypertension - reduced model")

# --- 8.3 Model comparison ----------------------------------------------------
model_comparison <- bind_rows(eval_dm, eval_dm_red, eval_htn, eval_htn_red)
model_comparison
write_csv(model_comparison, file.path(OUT_DIR, "model_comparison.csv"))


# =============================================================================
# 9. CORRELATION (AGE AND HEIGHT) ----
# =============================================================================

# Spearman rank correlation (non-parametric)
cor_age_height <- cor.test(Banco$Idade, Banco$Altura, method = "spearman")
cor_age_height
cat("Spearman rho:", round(cor_age_height$estimate, 3),
    "| p-value:", signif(cor_age_height$p.value, 3), "\n")
