# Diabetes and Hypertension Risk Factors

R analysis of a cross-sectional survey on non-communicable diseases (NCDs). It estimates the prevalence of
diabetes, pre-diabetes, hypertension and excess weight, and identifies associated lifestyle and family-history factors.

## Contents
- `ncd_diabetes_hypertension.R`: full pipeline (recoding, prevalence, Fisher tests, crude and adjusted logistic models, ROC)

## Methods
1. **Prevalence** with exact binomial 95% CIs
2. **Bivariate:** Fisher's exact test of each factor against diabetes and hypertension
3. **Crude odds ratios** (one logistic model per factor)
4. **Adjusted logistic models** (full and reduced) for diabetes and hypertension
5. **Model checks:** AUC, Hosmer-Lemeshow test, variance inflation factors

## Data (not included)
Participant-level data are confidential and excluded by `.gitignore`. Place the CSV at
`data/Dataset_DCNT_Tete.csv` (semicolon-separated) or edit `DATA_FILE`. Expected columns include:
`Idade, Peso, Altura, IMC, Sexo, Diagnostigo_do_diabetes, Hipertensao, Tabagismo, Atividade_fisica,
Pergunta_2, Pergunta_3, Pergunta_5, Pergunta_9, Pergunta_12, Pergunta_14, Pergunta_16, Pergunta_17`.

## Run
Open the script in RStudio and run it section by section. Packages are installed automatically
(`tidyverse, broom, knitr, pROC, ResourceSelection, car`). Results are written to `output/`.
