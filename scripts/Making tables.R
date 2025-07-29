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

# # Build and display the table
# gt(rep_table) %>%
#   fmt_markdown(columns = "Species") %>%
#   tab_options(
#     table.font.size = "small",
#     data_row.padding = px(2)
#   )
# 
# # Title for table
# tab_header(
#   title = md("**Table 1.** Overview of experimental species and sample sizes per treatment across elevations.")
# )
# 
# # Build and store the gt table
# rep_table_gt <- gt(rep_table) %>%
#   fmt_markdown(columns = "Species") %>%
#   tab_header(
#     title = md("**Table 1.** Overview of experimental species and sample sizes per treatment across elevations.")
#   ) %>%
#   tab_options(
#     table.font.size = "small",
#     data_row.padding = px(2),
#     column_labels.font.weight = "bold"
#   )
# 
# gtsave(rep_table_gt, "tables/table1.html")

readr::write_csv(rep_table, "tables/table1.csv")









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
  #n_obs <- nobs(model)
  # ΔAIC
  delta_aic <- round(AIC(null_model) - AIC(model), 3)
  # Likelihood ratio test for elevation effect
  anova_res <- tryCatch({
    anova(null_model, model)$`Pr(>Chisq)`[2]
  }, error = function(e) NA)
  elevation_sig <- if (!is.na(anova_res) && anova_res < 0.05) "Yes" else "No"
  # Format p-value
  # p_value_fmt <- case_when(
  #   is.na(anova_res)        ~ "NA",
  #   anova_res < 0.001       ~ "< 0.001",
  #   TRUE                    ~ formatC(anova_res, format = "f", digits = 3)
  #)
  tibble(
    Index = index_name_display,
    Distribution = dist_family,
    Zero_Inflated = zero_inflated,
    #Observations = n_obs,
    Delta_AIC = delta_aic,
    Elevation_Effect = elevation_sig,
    #P_value = p_value_fmt
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

readr::write_csv(table2, "tables/table2.csv")











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
  #n_obs <- nobs(model)
  # ΔAIC
  delta_aic <- round(AIC(null_model) - AIC(model), 3)
  # Elevation effect (likelihood ratio test)
  anova_res <- tryCatch({
    anova(null_model, model)$`Pr(>Chisq)`[2]
  }, error = function(e) NA)
  elevation_sig <- if (!is.na(anova_res) && anova_res < 0.05) "Yes" else "No"
  # Format p-value
  # p_value_fmt <- case_when(
  #   is.na(anova_res)        ~ "NA",
  #   anova_res < 0.001       ~ "< 0.001",
  #   TRUE                    ~ formatC(anova_res, format = "f", digits = 3)
  # )
  tibble(
    Response = index_name_display,
    Distribution = dist_family,
    Zero_Inflated = zero_inflated,
    #Observations = n_obs,
    Delta_AIC = delta_aic,
    Elevation_Effect = elevation_sig,
    #P_value = p_value_fmt
  )
}

table3 <- bind_rows(
  extract_model_summary("Visitation frequency", v_model, v_null),
  extract_model_summary("Morphospecies richness", m_model, m_null),
  extract_model_summary("Functional group richness", f_model, f_null)
)


readr::write_csv(table3, "tables/table3.csv")







#TABLE 3.5
# Rename column in table3 to match table2
table3 <- table3 %>%
  rename(Index = Response)

# Combine the two tables
combined_table <- bind_rows(table2, table3)

# View or export the combined table
print(combined_table)

readr::write_csv(combined_table, "tables/table3.5.csv")
























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

readr::write_csv(table4, "tables/table4.csv")

# Save
gtsave(gt_table4, "tables/table4.html")
#gtsave(gt_table4, "tables/table4.pdf")






















library(dplyr)
library(emmeans)

# Load models
c_model <- readRDS("glm_outputs/c_model.rds")
pl_model <- readRDS("glm_outputs/pl_model.rds")
v_model <- readRDS("glm_outputs/v_model.rds")
m_model <- readRDS("glm_outputs/m_model.rds")
f_model <- readRDS("glm_outputs/f_model.rds")

# Extract pairwise post-hoc comparisons
c_emm <- pairs(emmeans(c_model, ~ elevation))
pl_emm <- pairs(emmeans(pl_model, ~ elevation))
v_emm <- pairs(emmeans(v_model, ~ elevation))
m_emm <- pairs(emmeans(m_model, ~ elevation))
f_emm <- pairs(emmeans(f_model, ~ elevation))

# Convert to data frames
c_df <- as.data.frame(c_emm)[, c("contrast", "p.value")]
pl_df <- as.data.frame(pl_emm)[, c("contrast", "p.value")]
v_df <- as.data.frame(v_emm)[, c("contrast", "p.value")]
m_df <- as.data.frame(m_emm)[, c("contrast", "p.value")]
f_df <- as.data.frame(f_emm)[, c("contrast", "p.value")]

# Create full contrast list
all_contrasts <- unique(c(
  c_df$contrast, pl_df$contrast, v_df$contrast, m_df$contrast, f_df$contrast
))

format_pval <- function(p) {
  if (is.na(p)) {
    return("")
  } else if (p < 0.001) {
    return("p < 0.001")
  } else if (p < 0.01) {
    return("p < 0.01")
  } else if (p < 0.05) {
    return(paste0("p = ", format(round(p, 3), nsmall = 3)))
  } else {
    return(paste0("p = ", format(round(p, 2), nsmall = 2)))
  }
}


# Initialize result table
posthoc_table <- data.frame(Contrast = all_contrasts)

# Merge p-values into the table
posthoc_table$`Natural seed-set` <- c_df$p.value[match(posthoc_table$Contrast, c_df$contrast)]
posthoc_table$`Pollen limitation` <- pl_df$p.value[match(posthoc_table$Contrast, pl_df$contrast)]
posthoc_table$`Visitation frequency` <- v_df$p.value[match(posthoc_table$Contrast, v_df$contrast)]
posthoc_table$`Morphospecies richness` <- m_df$p.value[match(posthoc_table$Contrast, m_df$contrast)]
posthoc_table$`Functional group richness` <- f_df$p.value[match(posthoc_table$Contrast, f_df$contrast)]

posthoc_table_formatted <- posthoc_table
posthoc_table_formatted[, -1] <- lapply(posthoc_table_formatted[, -1], function(col) sapply(col, format_pval))


# View the final table
print(posthoc_table_formatted)

readr::write_csv(posthoc_table_formatted, "tables/table5.csv")




















# Convert to data frames including estimate
c_df <- as.data.frame(c_emm)[, c("contrast", "estimate", "p.value")]
pl_df <- as.data.frame(pl_emm)[, c("contrast", "estimate", "p.value")]
v_df <- as.data.frame(v_emm)[, c("contrast", "estimate", "p.value")]
m_df <- as.data.frame(m_emm)[, c("contrast", "estimate", "p.value")]
f_df <- as.data.frame(f_emm)[, c("contrast", "estimate", "p.value")]

posthoc_table <- data.frame(Contrast = all_contrasts)

# Helper function for estimates
format_estimate <- function(e) ifelse(is.na(e), "", format(round(e, 2), nsmall = 2))

# Add values to the table
posthoc_table$`Seedset_estimate` <- c_df$estimate[match(posthoc_table$Contrast, c_df$contrast)]
posthoc_table$`Seedset_pval`     <- c_df$p.value[match(posthoc_table$Contrast, c_df$contrast)]

posthoc_table$`PL_estimate` <- pl_df$estimate[match(posthoc_table$Contrast, pl_df$contrast)]
posthoc_table$`PL_pval`     <- pl_df$p.value[match(posthoc_table$Contrast, pl_df$contrast)]

posthoc_table$`VF_estimate` <- v_df$estimate[match(posthoc_table$Contrast, v_df$contrast)]
posthoc_table$`VF_pval`     <- v_df$p.value[match(posthoc_table$Contrast, v_df$contrast)]

posthoc_table$`MR_estimate` <- m_df$estimate[match(posthoc_table$Contrast, m_df$contrast)]
posthoc_table$`MR_pval`     <- m_df$p.value[match(posthoc_table$Contrast, m_df$contrast)]

posthoc_table$`FR_estimate` <- f_df$estimate[match(posthoc_table$Contrast, f_df$contrast)]
posthoc_table$`FR_pval`     <- f_df$p.value[match(posthoc_table$Contrast, f_df$contrast)]


# Apply formatting
posthoc_table_formatted <- posthoc_table
posthoc_table_formatted[grep("_pval", names(posthoc_table_formatted))] <- lapply(
  posthoc_table_formatted[grep("_pval", names(posthoc_table_formatted))],
  function(col) sapply(col, format_pval)
)

posthoc_table_formatted[grep("_estimate", names(posthoc_table_formatted))] <- lapply(
  posthoc_table_formatted[grep("_estimate", names(posthoc_table_formatted))],
  function(col) sapply(col, format_estimate)
)

readr::write_csv(posthoc_table_formatted, "tables/table3_with_estimates.csv")





















