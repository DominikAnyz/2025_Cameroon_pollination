###* Tables for the manuscript
###* 
###* FIRST table, plant species summary
###* 
###* I would like to create a table in which I have information about the amounts 
###* of observations from distinct treatmen
# Load packages
pacman::p_load(tidyverse, glmmTMB, DHARMa, gt)

# Ensure dplyr::select is preferred
select <- dplyr::select

# Set seed for reproducibility
set.seed(1234)

# Load and preprocess data
seed.data <- read.delim("data/clean_seeds.txt", na = c("na"))

# Calculate indices
seed.indices <- 
  seed.data %>%
  group_by(elevation, plant_number, species) %>%
  mutate(O.mean.plantnumber = mean(seedset[treatment == "outcrossing"], na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(elevation, species) %>%
  mutate(O.mean.elevation = mean(seedset[treatment == "outcrossing"], na.rm = TRUE)) %>%
  mutate(O.mean = case_when(
    is.nan(O.mean.plantnumber) ~ O.mean.elevation,
    .default = O.mean.plantnumber
  )) %>%
  mutate(index = seedset / O.mean) %>%
  mutate(seedset = if_else(seedset != round(seedset), round(seedset), seedset))

# Clean index values and structure
seed.indices <- seed.indices %>%
  filter(!(index %in% c("NA", "Inf") | is.na(seedset))) %>%
  mutate(index = ifelse(is.finite(index), index, 0))

# Ensure factor variables
seed.indices$elevation <- as.factor(seed.indices$elevation)
seed.indices$species <- as.factor(seed.indices$species)

# Add elevation.species column
seed.indices <- seed.indices %>%
  mutate(elevation.species = paste0(elevation, species))

# Create replicate summary for table
replicate_summary_2 <- seed.indices %>%
  mutate(treatment = recode(treatment,
                            "control" = "C",
                            "autogamy" = "A",
                            "geitonogamy" = "G",
                            "outcrossing" = "O")) %>%
  filter(treatment %in% c("A", "G", "O", "C")) %>%
  group_by(species, elevation, treatment) %>%
  summarise(n_replicates = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = treatment,
    values_from = n_replicates,
    values_fill = 0
  ) %>%
  arrange(species, elevation)








###* TABLE 1
# Recode incorrect elevation values
replicate_summary_2 <- replicate_summary_2 %>%
  mutate(elevation = recode(elevation,
                            `3500` = "3400",
                            `4000` = "3800"))

replicate_summary_2$elevation <- factor(replicate_summary_2$elevation,
                                        levels = c("2300", "2800", "3400", "3800"))


# Format species names with italics and family names
rep_table <- replicate_summary_2 %>%
  mutate(Species = case_when(
    species == "Clematis s"  ~ "<i>Clematis simensis</i> (Ranunculaceae)",
    species == "Crepis h"    ~ "<i>Crepis hypochoeridea</i> (Asteraceae)",
    species == "Geranium a"  ~ "<i>Geranium arabicum</i> (Geraniaceae)",
    species == "Hypericum r" ~ "<i>Hypericum revolutum</i> (Hypericaceae)",
    species == "Lactuca i"   ~ "<i>Lactuca inermis</i> (Asteraceae)",
    species == "Senecio b"   ~ "<i>Senecio burtonii</i> (Asteraceae)",
    species == "Senecio p"   ~ "<i>Senecio purpureus</i> (Asteraceae)"
  )) %>%
  arrange(Species, elevation) %>%
  group_by(Species) %>%
  mutate(Species = ifelse(row_number() == 1, Species, "")) %>%
  ungroup() %>%
  select(Species, Elevation = elevation, Autogamy = A, Geitonogamy = G,
         Outcrossing = O, Control = C)

# Build and display the table
gt(rep_table) %>%
  fmt_markdown(columns = "Species") %>%
  tab_options(
    table.font.size = "small",
    data_row.padding = px(2)
  )

# Title for table
tab_header(
  title = md("**Table 1.** Overview of experimental species and sample sizes per treatment across elevations.")
)

# Build and store the gt table
rep_table_gt <- gt(rep_table) %>%
  fmt_markdown(columns = "Species") %>%
  tab_header(
    title = md("**Table 1.** Overview of experimental species and sample sizes per treatment across elevations.")
  ) %>%
  tab_options(
    table.font.size = "small",
    data_row.padding = px(2),
    column_labels.font.weight = "bold"
  )

gtsave(rep_table_gt, "tables/table1.html")








###* TABLE 2

c_model <- readRDS("glm_outputs/c_model.rds")
c_null  <- readRDS("glm_outputs/c_null.rds")
pl_model <- readRDS("glm_outputs/pl_model.rds")
pl_null  <- readRDS("glm_outputs/pl_null.rds")
ao_model <- readRDS("glm_outputs/ao_model.rds")
ao_null  <- readRDS("glm_outputs/ao_null.rds")
go_model <- readRDS("glm_outputs/go_model.rds")
go_null  <- readRDS("glm_outputs/go_null.rds")

extract_model_summary <- function(index_name_display, model, null_model) {
  # Extract distribution family
  dist_family <- as.character(model$call$family)[1]
  # Check zero-inflation
  ziform <- model$call$ziformula
  zero_inflated <- if (is.null(ziform)) {
    "No"
  } else {
    zi_char <- as.character(ziform)[2]
    if (zi_char == "0") "No" else "Yes"
  }
  # Number of observations
  n_obs <- nobs(model)
  # ΔAIC
  delta_aic <- round(AIC(null_model) - AIC(model), 3)
  # Likelihood ratio test for elevation effect
  anova_res <- tryCatch({
    anova(null_model, model)$`Pr(>Chisq)`[2]
  }, error = function(e) NA)
  elevation_sig <- if (!is.na(anova_res) && anova_res < 0.05) "Yes" else "No"
  # Format p-value
  p_value_fmt <- case_when(
    is.na(anova_res)        ~ "NA",
    anova_res < 0.001       ~ "< 0.001",
    TRUE                    ~ formatC(anova_res, format = "f", digits = 3)
  )
  tibble(
    Index = index_name_display,
    Distribution = dist_family,
    Zero_Inflated = zero_inflated,
    Observations = n_obs,
    Delta_AIC = delta_aic,
    Elevation_Effect = elevation_sig,
    P_value = p_value_fmt
  )
}

table2 <- bind_rows(
  extract_model_summary("Natural seed-set", c_model, c_null),
  extract_model_summary("Pollen limitaion", pl_model, pl_null),
  extract_model_summary("Autogamy", ao_model, ao_null),
  extract_model_summary("Geitonogamy", go_model, go_null)
)

# View table
print(table2)

# Create the gt table
gt_table2 <- table2 %>%
  gt() %>%
  cols_label(
    Index = "Index",
    Distribution = "Distribution",
    Zero_Inflated = "Zero-inflated",
    Observations = "n",
    Delta_AIC = "ΔAIC",
    Elevation_Effect = "Elevation effect",
    P_value = "p-value"
  ) %>%
  cols_align(align = "center", columns = everything()) %>%
  tab_options(
    table.font.size = "small",
    data_row.padding = px(2),
    column_labels.font.weight = "bold"
  ) %>%
  tab_header(
    title = md("**Table 2.** Generalized linear mixed models testing the effect of elevation on reproductive indices.")
  )

gt_table2

gtsave(gt_table2, "tables/table2.html")











###* TABLE 3
v_model <- readRDS("glm_outputs/v_model.rds")
v_null  <- readRDS("glm_outputs/v_null.rds")
m_model <- readRDS("glm_outputs/m_model.rds")
m_null  <- readRDS("glm_outputs/m_null.rds")
f_model <- readRDS("glm_outputs/f_model.rds")
f_null  <- readRDS("glm_outputs/f_null.rds")

extract_model_summary <- function(index_name_display, model, null_model) {
  # Distribution family
  dist_family <- as.character(model$call$family)[1]
  # Zero-inflation check
  ziform <- model$call$ziformula
  zero_inflated <- if (is.null(ziform)) {
    "No"
  } else {
    zi_char <- as.character(ziform)[2]
    if (zi_char == "0") "No" else "Yes"
  }
  # Number of observations
  n_obs <- nobs(model)
  # ΔAIC
  delta_aic <- round(AIC(null_model) - AIC(model), 3)
  # Elevation effect (likelihood ratio test)
  anova_res <- tryCatch({
    anova(null_model, model)$`Pr(>Chisq)`[2]
  }, error = function(e) NA)
  elevation_sig <- if (!is.na(anova_res) && anova_res < 0.05) "Yes" else "No"
  # Format p-value
  p_value_fmt <- case_when(
    is.na(anova_res)        ~ "NA",
    anova_res < 0.001       ~ "< 0.001",
    TRUE                    ~ formatC(anova_res, format = "f", digits = 3)
  )
  tibble(
    Response = index_name_display,
    Distribution = dist_family,
    Zero_Inflated = zero_inflated,
    Observations = n_obs,
    Delta_AIC = delta_aic,
    Elevation_Effect = elevation_sig,
    P_value = p_value_fmt
  )
}

table3 <- bind_rows(
  extract_model_summary("Visitation frequency", v_model, v_null),
  extract_model_summary("Morphospecies richness", m_model, m_null),
  extract_model_summary("Functional group richness", f_model, f_null)
)

gt_table3 <- table3 %>%
  gt() %>%
  cols_label(
    Response = "Response variable",
    Distribution = "Distribution",
    Zero_Inflated = "Zero-inflated",
    Observations = "n",
    Delta_AIC = "ΔAIC",
    Elevation_Effect = "Elevation effect",
    P_value = "p-value"
  ) %>%
  cols_align(align = "center", columns = everything()) %>%
  tab_options(
    table.font.size = "small",
    data_row.padding = px(2),
    column_labels.font.weight = "bold"
  ) %>%
  tab_header(
    title = md("**Table 3.** Generalized linear mixed models testing the effect of elevation on pollinator visitation metrics.")
  )

gtsave(gt_table3, "tables/table3.html")























###* TABLE 4
source("scripts/Bayesian setup June.R")

c.pl.final.table.5 <- c.pl.final.table.4 %>%
  filter(species != "Hypericum r" | elevation != 4000)

extract_bayes_model_info <- function(index_name, model = NULL, null_model = NULL, distribution = NA) {
  fit_success <- !is.null(model)
  
  if (!fit_success) {
    return(tibble(
      Index = index_name,
      Model_Fit = "No",
      "ΔR²" = "–",
      Distribution = distribution,
      Best_Predictors = "–",
      Estimate = "–",
      CI_95 = "–",
      Significant = "–"
    ))
  }
  
  # Bayes R² and delta
  r2_val <- bayes_R2(model)[1, "Estimate"]
  null_r2 <- if (!is.null(null_model)) bayes_R2(null_model)[1, "Estimate"] else NA
  delta_r2 <- if (!is.na(null_r2)) round(r2_val - null_r2, 2) else NA
  
  # If delta R2 is missing or zero or negative, treat as model not better than null
  if (is.na(delta_r2) || delta_r2 <= 0) {
    return(tibble(
      Index = index_name,
      Model_Fit = "Yes",
      "ΔR²" = as.character(delta_r2),
      Distribution = distribution,
      Best_Predictors = "–",
      Estimate = "–",
      CI_95 = "–",
      Significant = "–"
    ))
  }
  
  # Extract fixed effects (excluding intercept)
  coefs <- fixef(model)
  predictors <- rownames(coefs)[rownames(coefs) != "Intercept"]
  
  predictor_labels <- predictors %>%
    str_replace_all(c(
      "memean\\.visited\\.flowers\\.scaled.*" = "visitation",
      "memean\\.morpho\\.scaled.*" = "morphospecies",
      "memean\\.func\\.scaled.*" = "functional group"
    ))
  
  estimates <- round(coefs[predictors, "Estimate"], 2)
  estimates_str <- paste(estimates, collapse = "; ")
  
  ci_95 <- apply(coefs[predictors, , drop = FALSE], 1, function(row) {
    paste0("[", round(row["Q2.5"], 2), ", ", round(row["Q97.5"], 2), "]")
  }) %>% paste(collapse = "; ")
  
  ci_flags <- coefs[predictors, "Q2.5"] < 0 & coefs[predictors, "Q97.5"] > 0
  significant <- ifelse(any(!ci_flags), "Yes", "No")
  
  return(tibble(
    Index = index_name,
    Model_Fit = "Yes",
    "ΔR²" = as.character(delta_r2),
    Distribution = distribution,
    Best_Predictors = paste(predictor_labels, collapse = " + "),
    Estimate = estimates_str,
    CI_95 = ci_95,
    Significant = significant
  ))
}


# Load your models
pl_model <- readRDS("../brms_models/bayesian_pl_noel_w12_p1_vmf.rds")
pl_null  <- readRDS("../brms_models/bayesian_pl_noel_w12_p1_null.rds")

ao_model <- readRDS("../brms_models/bayesian_ao_noel_w12_p1_f.rds")
ao_null  <- readRDS("../brms_models/bayesian_ao_noel_w12_p1_null.rds")

go_model <- readRDS("../brms_models/bayesian_go_noel_w12_p1_v.rds")
go_null  <- readRDS("../brms_models/bayesian_go_noel_w12_p1_null.rds")

# For control model, which failed
c_summary <- extract_bayes_model_info("Natural seed-set", model = NULL, distribution = "–")
pl_summary <- extract_bayes_model_info("Pollen Limitation", model = pl_model, null_model = pl_null, distribution = "zero one inflated beta")
ao_summary <- extract_bayes_model_info("Autogamy", model = ao_model, null_model = ao_null, distribution = "zero inflated beta")
go_summary <- extract_bayes_model_info("Geitonogamy", model = go_model, null_model = go_null, distribution = "–")

# Combine into one table
table4 <- bind_rows(c_summary, pl_summary, ao_summary, go_summary)

gt_table4 <- table4 %>%
  mutate(
    Best_Predictors = gsub(" \\+ ", " +<br>", Best_Predictors),
    Estimate = gsub("; ", ";<br>", Estimate),
    CI_95 = gsub("; ", ";<br>", CI_95)
  ) %>%
  gt() %>%
  fmt_markdown(columns = c("Best_Predictors", "Estimate", "CI_95")) %>%
  cols_label(
    Index = "Index",
    Model_Fit = "Model fit",
    `ΔR²` = "ΔR²",
    Distribution = "Distribution",
    Best_Predictors = "Best predictors",
    Estimate = "Estimate",
    CI_95 = "95% CI",
    Significant = "Significant"
  ) %>%
  cols_align(align = "center", columns = everything()) %>%
  tab_options(
    table.font.size = "small",
    data_row.padding = px(4),
    column_labels.font.weight = "bold"
  ) %>%
  tab_header(
    title = md("**Table 4.** Bayesian models testing the relationship between reproductive indices and pollinator visitation metrics.")
  )

# Save
gtsave(gt_table4, "tables/table4.html")
#gtsave(gt_table4, "tables/table4.pdf")





