###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* Making tables - supplementary
###* 
library(loo)
library(dplyr)
library(tibble)



###* 
###* 
###* 
###* 
###* Table S1
c_brm_model_zinb_new <- readRDS("brms_models/c_brm_model_zinb_new.rds")
pl_ord_model_new <- readRDS("brms_models/pl_ord_model_new.rds")
ao_brm_model_zinb_new <- readRDS("brms_models/ao_brm_model_zinb_new.rds")
go_brm_model_zinb_new <- readRDS("brms_models/go_brm_model_zinb_new.rds")

library(emmeans)
library(dplyr)

# Example: contrasts for Natural seed-set model
c_emm <- emmeans(c_brm_model_zinb_new, ~ elevation)
c_contr <- contrast(c_emm, method = "pairwise") %>%
  as.data.frame()

# Format columns
c_contr <- c_contr %>%
  mutate(
    Metric = "Natural seed-set",
    Estimate = round(estimate, 2),
    l95 = round(lower.HPD, 2),
    u95 = round(upper.HPD, 2)
  ) %>%
  select(Metric, contrast, Estimate, l95, u95)

# Pollen limitation
pl_emm <- emmeans(pl_ord_model_new, ~ elevation)
pl_contr <- contrast(pl_emm, method = "pairwise") %>%
  as.data.frame() %>%
  mutate(Metric = "Pollen limitation",
         Estimate = round(estimate, 2),
         l95 = round(lower.HPD, 2),
         u95 = round(upper.HPD, 2)) %>%
  select(Metric, contrast, Estimate, l95, u95)

# ... repeat for AO, GO, M, F, Visitation

# Autogamy index
ao_emm <- emmeans(ao_brm_model_zinb_new, ~ elevation)
ao_contr <- contrast(ao_emm, method = "pairwise") %>%
  as.data.frame() %>%
  mutate(Metric = "Autogamy",
         Estimate = round(estimate, 2),
         l95 = round(lower.HPD, 2),
         u95 = round(upper.HPD, 2)) %>%
  select(Metric, contrast, Estimate, l95, u95)

# Geitonogamy index
summary(go_brm_model_zinb_new)

go_emm <- emmeans(go_brm_model_zinb_new, ~ elevation)
go_contr <- contrast(go_emm, method = "pairwise") %>%
  as.data.frame() %>%
  mutate(Metric = "Geitonogamy",
         Estimate = round(estimate, 2),
         l95 = round(lower.HPD, 2),
         u95 = round(upper.HPD, 2)) %>%
  select(Metric, contrast, Estimate, l95, u95)

go_contr
# Morpho richness
m_emm <- emmeans(m_brm_model_nb, ~ elevation)
m_contr <- contrast(m_emm, method = "pairwise") %>%
  as.data.frame() %>%
  mutate(Metric = "Morphospecies richness",
         Estimate = round(estimate, 2),
         l95 = round(lower.HPD, 2),
         u95 = round(upper.HPD, 2)) %>%
  select(Metric, contrast, Estimate, l95, u95)

# Functional richness
f_emm <- emmeans(f_brm_model_po, ~ elevation)
f_contr <- contrast(f_emm, method = "pairwise") %>%
  as.data.frame() %>%
  mutate(Metric = "Functional group richness",
         Estimate = round(estimate, 2),
         l95 = round(lower.HPD, 2),
         u95 = round(upper.HPD, 2)) %>%
  select(Metric, contrast, Estimate, l95, u95)

# Visitation frequency
v_emm <- emmeans(vis_ord_model, ~ elevation)
v_contr <- contrast(v_emm, method = "pairwise") %>%
  as.data.frame() %>%
  mutate(Metric = "Visitation frequency",
         Estimate = round(estimate, 2),
         l95 = round(lower.HPD, 2),
         u95 = round(upper.HPD, 2)) %>%
  select(Metric, contrast, Estimate, l95, u95)

all_contrasts <- bind_rows(c_contr, pl_contr, ao_contr, m_contr, f_contr, v_contr)

all_contrasts <- all_contrasts %>%
  mutate(
    CI_95 = paste0("[", l95, ", ", u95, "]")
  ) %>%
  select(Metric, contrast, Estimate, CI_95)

print(all_contrasts)

readr::write_csv(all_contrasts, "tables/new_bayes_post.csv")
###*
###*
###*
###*
###*
###*









###* COntrol seedset
###* Prepared loo images for models
c_brm_model_zinb_new_k <- readRDS("brms_models/c_brm_model_zinb_new_k.rds")
c_brm_model_po_new_k <- readRDS("brms_models/c_brm_model_po_new_k.rds")
c_brm_model_zipo_new_k <- readRDS("brms_models/c_brm_model_zipo_new_k.rds")
c_brm_model_nb_new_k <- readRDS("brms_models/c_brm_model_nb_new_k.rds")

# Put the k-fold objects in a named list
mods <- list(
  po   = c_brm_model_po_new_k,
  zipo = c_brm_model_zipo_new_k,
  nb   = c_brm_model_nb_new_k,
  zinb = c_brm_model_zinb_new_k
)

# 1) Run loo_compare on these k-fold objects
lc <- loo_compare(mods)   # uses elpd_kfold internally

# 2) Extract absolute ELPD (elpd_kfold) for each model
elpd_vec <- sapply(mods, function(x) x$estimates["elpd_kfold", "Estimate"])

# 3) Build the summary table
tab_models <- as.data.frame(lc) %>%
  rownames_to_column("model") %>%
  mutate(
    Distribution = case_when(
      model == "po"   ~ "Poisson (po)",
      model == "zipo" ~ "Zero-inflated Poisson (zipo)",
      model == "nb"   ~ "Negative binomial (nb)",
      model == "zinb" ~ "Zero-inflated negative binomial (zinb)",
      TRUE ~ model
    ),
    ELPD   = elpd_vec[model],                              # absolute elpd_kfold
    z_ELPD = if_else(se_diff > 0, elpd_diff / se_diff, NA_real_)  # best model = NA
  ) %>%
  select(Distribution, ELPD, elpd_diff, se_diff, z_ELPD)

tab_models

write.csv(tab_models,
          "tables/tab_models_seedset_distributions.csv",
          row.names = FALSE)





###* Autogamy
ao_brm_model_zinb_new_k <- readRDS("brms_models/ao_brm_model_zinb_new_k.rds")
ao_brm_model_po_new_k <- readRDS("brms_models/ao_brm_model_po_new_k.rds")
ao_brm_model_zipo_new_k <- readRDS("brms_models/ao_brm_model_zipo_new_k.rds")
ao_brm_model_nb_new_k <- readRDS("brms_models/ao_brm_model_nb_new_k.rds")


# Put the k-fold objects in a named list
mods <- list(
  po   = ao_brm_model_po_new_k,
  zipo = ao_brm_model_zipo_new_k,
  nb   = ao_brm_model_nb_new_k,
  zinb = ao_brm_model_zinb_new_k
)

# 1) Run loo_compare on these k-fold objects
lc <- loo_compare(mods)   # uses elpd_kfold internally

# 2) Extract absolute ELPD (elpd_kfold) for each model
elpd_vec <- sapply(mods, function(x) x$estimates["elpd_kfold", "Estimate"])

# 3) Build the summary table
tab_models <- as.data.frame(lc) %>%
  rownames_to_column("model") %>%
  mutate(
    Distribution = case_when(
      model == "po"   ~ "Poisson (po)",
      model == "zipo" ~ "Zero-inflated Poisson (zipo)",
      model == "nb"   ~ "Negative binomial (nb)",
      model == "zinb" ~ "Zero-inflated negative binomial (zinb)",
      TRUE ~ model
    ),
    ELPD   = elpd_vec[model],                              # absolute elpd_kfold
    z_ELPD = if_else(se_diff > 0, elpd_diff / se_diff, NA_real_)  # best model = NA
  ) %>%
  select(Distribution, ELPD, elpd_diff, se_diff, z_ELPD)

tab_models

write.csv(tab_models,
          "tables/tab_models_autogamy_distributions.csv",
          row.names = FALSE)



###* Geitonogamy
go_brm_model_zinb_new_k <- readRDS("brms_models/go_brm_model_zinb_new_k.rds")
go_brm_model_po_new_k <- readRDS("brms_models/go_brm_model_po_new_k.rds")
go_brm_model_zipo_new_k <- readRDS("brms_models/go_brm_model_zipo_new_k.rds")
go_brm_model_nb_new_k <- readRDS("brms_models/go_brm_model_nb_new_k.rds")


# Put the k-fold objects in a named list
mods <- list(
  po   = go_brm_model_po_new_k,
  zipo = go_brm_model_zipo_new_k,
  nb   = go_brm_model_nb_new_k,
  zinb = go_brm_model_zinb_new_k
)

# 1) Run loo_compare on these k-fold objects
lc <- loo_compare(mods)   # uses elpd_kfold internally

# 2) Extract absolute ELPD (elpd_kfold) for each model
elpd_vec <- sapply(mods, function(x) x$estimates["elpd_kfold", "Estimate"])

# 3) Build the summary table
tab_models <- as.data.frame(lc) %>%
  rownames_to_column("model") %>%
  mutate(
    Distribution = case_when(
      model == "po"   ~ "Poisson (po)",
      model == "zipo" ~ "Zero-inflated Poisson (zipo)",
      model == "nb"   ~ "Negative binomial (nb)",
      model == "zinb" ~ "Zero-inflated negative binomial (zinb)",
      TRUE ~ model
    ),
    ELPD   = elpd_vec[model],                              # absolute elpd_kfold
    z_ELPD = if_else(se_diff > 0, elpd_diff / se_diff, NA_real_)  # best model = NA
  ) %>%
  select(Distribution, ELPD, elpd_diff, se_diff, z_ELPD)

tab_models

write.csv(tab_models,
          "tables/tab_models_geitonogamy_distributions.csv",
          row.names = FALSE)
























###* MORPHOSPECIES RICHNESS
m_brm_model_po <- readRDS("brms_models/m_brm_model_po.rds")
m_brm_model_zipo <- readRDS("brms_models/m_brm_model_zipo.rds")
m_brm_model_nb <- readRDS("brms_models/m_brm_model_nb.rds")
m_brm_model_zinb <- readRDS("brms_models/m_brm_model_zinb.rds")

# 1) LOO for the best model (zinb) to get its absolute ELPD
loo_nb <- loo(m_brm_model_nb)
best_elpd <- loo_nb$estimates["elpd_loo", "Estimate"]

loo_zinb <- loo(m_brm_model_zinb)
loo_po <- loo(m_brm_model_po)
loo_zipo <- loo(m_brm_model_zipo)

# 2) Compare all candidate models
lc <- loo_compare(
  loo_zipo,
  loo_po,
  loo_nb,
  loo_zinb
)

lc


# 3) Turn the compare object into a tidy data frame
tab_models <- as.data.frame(lc) %>%
  rownames_to_column("model") %>%
  mutate(
    Distribution = case_when(
      model == "loo_po"   ~ "Poisson (po)",
      model == "loo_zipo" ~ "Zero-inflated Poisson (zipo)",
      model == "loo_nb"   ~ "Negative binomial (nb)",
      model == "loo_zinb" ~ "Zero-inflated negative binomial (zinb)",
      TRUE ~ model
    ),
    # reconstruct absolute ELPD using the best model's elpd
    ELPD   = best_elpd + elpd_diff,
    z_ELPD = if_else(se_diff > 0, elpd_diff / se_diff, NA_real_)  # best model gets NA
  ) %>%
  select(Distribution, ELPD, elpd_diff, se_diff, z_ELPD)

tab_models

write.csv(tab_models,
          "tables/tab_models_morpho_distributions.csv",
          row.names = FALSE)












###* FUNCTIONAL GROUP RICHNESS
f_brm_model_po <- readRDS("brms_models/f_brm_model_po.rds")
f_brm_model_zipo <- readRDS("brms_models/f_brm_model_zipo.rds")
f_brm_model_nb <- readRDS("brms_models/f_brm_model_nb.rds")
f_brm_model_zinb <- readRDS("brms_models/f_brm_model_zinb.rds")

# 1) LOO for the best model (zinb) to get its absolute ELPD
loo_po <- loo(f_brm_model_po)
best_elpd <- loo_po$estimates["elpd_loo", "Estimate"]

loo_zinb <- loo(f_brm_model_zinb)
loo_nb <- loo(f_brm_model_nb)
loo_zipo <- loo(f_brm_model_zipo)

# 2) Compare all candidate models
lc <- loo_compare(
  loo_zipo,
  loo_po,
  loo_nb,
  loo_zinb
)

lc


# 3) Turn the compare object into a tidy data frame
tab_models <- as.data.frame(lc) %>%
  rownames_to_column("model") %>%
  mutate(
    Distribution = case_when(
      model == "loo_po"   ~ "Poisson (po)",
      model == "loo_zipo" ~ "Zero-inflated Poisson (zipo)",
      model == "loo_nb"   ~ "Negative binomial (nb)",
      model == "loo_zinb" ~ "Zero-inflated negative binomial (zinb)",
      TRUE ~ model
    ),
    # reconstruct absolute ELPD using the best model's elpd
    ELPD   = best_elpd + elpd_diff,
    z_ELPD = if_else(se_diff > 0, elpd_diff / se_diff, NA_real_)  # best model gets NA
  ) %>%
  select(Distribution, ELPD, elpd_diff, se_diff, z_ELPD)

tab_models

write.csv(tab_models,
          "tables/tab_models_func_distributions.csv",
          row.names = FALSE)



























kfold_13 <- readRDS("brms_models/kfold_bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt.rds")
kfold_1 <- readRDS("brms_models/kfold_bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt.rds")
kfold_3 <- readRDS("brms_models/kfold_bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt.rds")
kfold_5 <- readRDS("brms_models/kfold_bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt.rds")
kfold_35 <- readRDS("brms_models/kfold_bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt.rds")
kfold_135 <- readRDS("brms_models/kfold_bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt.rds")
kfold_15 <- readRDS("brms_models/kfold_bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt.rds")
kfold_null <- readRDS("brms_models/kfold_bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt.rds")



# 1) Put all K-fold objects (EXCEPT the null model) into a named list ----
# Replace the object names here with your actual k-fold objects
mods_pred <- list(
  m15  = kfold_15,    # visitation (1) + functional-group richness (5)
  m135 = kfold_135,   # 1 + 3 + 5
  m13  = kfold_13,    # 1 + 3
  m35  = kfold_35,    # 3 + 5
  m1   = kfold_1,     # 1 only
  m3   = kfold_3,     # 3 only
  m5   = kfold_5,      # 5 only
  mnull= kfold_null
)

# 2) Compare all candidate predictor sets
lc_pred <- loo_compare(mods_pred)   # works on k-fold objects, uses elpd_kfold
lc_pred

# 3) Get absolute ELPD (elpd_kfold) for each predictor set
elpd_vec_pred <- sapply(mods_pred, function(x) x$estimates["elpd_kfold", "Estimate"])

# 4) Build the summary table
tab_predictors <- as.data.frame(lc_pred) %>%
  rownames_to_column("model") %>%
  mutate(
    Predictor_combination = case_when(
      model == "m15"  ~ "1 + 5  (visitation + functional-group richness)",
      model == "m135" ~ "1 + 3 + 5 (visitation + morphospecies + functional-group)",
      model == "m13"  ~ "1 + 3  (visitation + morphospecies richness)",
      model == "m35"  ~ "3 + 5  (morphospecies + functional-group richness)",
      model == "m1"   ~ "1  (visitation only)",
      model == "m3"   ~ "3  (morphospecies richness only)",
      model == "m5"   ~ "5  (functional-group richness only)",
      model == "mnull"~ "null (intercept only)",
      TRUE ~ model
    ),
    ELPD   = elpd_vec_pred[model],                         # absolute elpd_kfold
    z_ELPD = if_else(se_diff > 0, elpd_diff / se_diff, NA_real_)  # best model = NA
  ) %>%
  select(Predictor_combination, ELPD, elpd_diff, se_diff, z_ELPD)

# 5) Reorder rows: 15, 135, 13, 35, then the rest
order_levels <- c(
  "1 + 5  (visitation + functional-group richness)",
  "1 + 3 + 5 (visitation + morphospecies + functional-group)",
  "1 + 3  (visitation + morphospecies richness)",
  "3 + 5  (morphospecies + functional-group richness)"
)

tab_predictors <- tab_predictors %>%
  mutate(
    Predictor_combination = factor(Predictor_combination,
                                   levels = c(order_levels,
                                              setdiff(Predictor_combination, order_levels))
    )
  ) %>%
  arrange(Predictor_combination)

tab_predictors

write.csv(tab_predictors,
          "tables/tab_models_control_variations.csv",
          row.names = FALSE)




###* POLLEN LIMITATION
kfold_pl_13 <- readRDS("brms_models/kfold_bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2.rds")
kfold_pl_null <- readRDS("brms_models/kfold_bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2.rds")
kfold_pl_135 <- readRDS("brms_models/kfold_bayesian_pl_mean_corrected_scaled_135_zoib_k_opt_2.rds")
kfold_pl_1 <- readRDS("brms_models/kfold_bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2.rds")
kfold_pl_3 <- readRDS("brms_models/kfold_bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2.rds")
kfold_pl_5 <- readRDS("brms_models/kfold_bayesian_pl_mean_corrected_scaled_5_zoib_k_opt_2.rds")
kfold_pl_35 <- readRDS("brms_models/kfold_bayesian_pl_mean_corrected_scaled_35_zoib_k_opt_2.rds")
kfold_pl_15 <- readRDS("brms_models/kfold_bayesian_pl_mean_corrected_scaled_15_zoib_k_opt_2.rds")


# 1) Put all K-fold objects (EXCEPT the null model) into a named list ----
# Replace the object names here with your actual k-fold objects
mods_pred_pl <- list(
  m15  = kfold_pl_15,    # visitation (1) + functional-group richness (5)
  m135 = kfold_pl_135,   # 1 + 3 + 5
  m13  = kfold_pl_13,    # 1 + 3
  m35  = kfold_pl_35,    # 3 + 5
  m1   = kfold_pl_1,     # 1 only
  m3   = kfold_pl_3,     # 3 only
  m5   = kfold_pl_5,      # 5 only
  mnull= kfold_pl_null
)

# 2) Compare all candidate predictor sets
lc_pred_pl <- loo_compare(mods_pred_pl)   # works on k-fold objects, uses elpd_kfold
lc_pred_pl

# 3) Get absolute ELPD (elpd_kfold) for each predictor set
elpd_vec_pred_pl <- sapply(mods_pred_pl, function(x) x$estimates["elpd_kfold", "Estimate"])

# 4) Build the summary table
tab_predictors_pl <- as.data.frame(lc_pred_pl) %>%
  rownames_to_column("model") %>%
  mutate(
    Predictor_combination = case_when(
      model == "m15"  ~ "1 + 5  (visitation + functional-group richness)",
      model == "m135" ~ "1 + 3 + 5 (visitation + morphospecies + functional-group)",
      model == "m13"  ~ "1 + 3  (visitation + morphospecies richness)",
      model == "m35"  ~ "3 + 5  (morphospecies + functional-group richness)",
      model == "m1"   ~ "1  (visitation only)",
      model == "m3"   ~ "3  (morphospecies richness only)",
      model == "m5"   ~ "5  (functional-group richness only)",
      model == "mnull"~ "null (intercept only)",
      TRUE ~ model
    ),
    ELPD   = elpd_vec_pred_pl[model],                         # absolute elpd_kfold
    z_ELPD = if_else(se_diff > 0, elpd_diff / se_diff, NA_real_)  # best model = NA
  ) %>%
  select(Predictor_combination, ELPD, elpd_diff, se_diff, z_ELPD)

tab_predictors_pl

write.csv(tab_predictors_pl,
          "tables/tab_models_pl_variations.csv",
          row.names = FALSE)





###* AUTOGAMY
kfold_ao_13 <- readRDS("brms_models/kfold_bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt.rds")
kfold_ao_null <- readRDS("brms_models/kfold_bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt.rds")
kfold_ao_135 <- readRDS("brms_models/kfold_bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt.rds")
kfold_ao_1 <- readRDS("brms_models/kfold_bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt.rds")
kfold_ao_3 <- readRDS("brms_models/kfold_bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt.rds")
kfold_ao_5 <- readRDS("brms_models/kfold_bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt.rds")
kfold_ao_35 <- readRDS("brms_models/kfold_bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt.rds")
kfold_ao_15 <- readRDS("brms_models/kfold_bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt.rds")

loo_compare(kfold_ao_13,
            kfold_ao_135,
            kfold_ao_1,
            kfold_ao_3,
            kfold_ao_35,
            kfold_ao_5,
            kfold_ao_15,
            kfold_ao_null)

# 1) Put all K-fold objects (EXCEPT the null model) into a named list ----
# Replace the object names here with your actual k-fold objects
mods_pred_ao <- list(
  m15  = kfold_ao_15,    # visitation (1) + functional-group richness (5)
  m135 = kfold_ao_135,   # 1 + 3 + 5
  m13  = kfold_ao_13,    # 1 + 3
  m35  = kfold_ao_35,    # 3 + 5
  m1   = kfold_ao_1,     # 1 only
  m3   = kfold_ao_3,     # 3 only
  m5   = kfold_ao_5,      # 5 only
  mnull= kfold_ao_null
)

loo_compare(mods_pred_ao)

# 2) Compare all candidate predictor sets
lc_pred_ao <- loo_compare(mods_pred_ao)   # works on k-fold objects, uses elpd_kfold
lc_pred_ao

# 3) Get absolute ELPD (elpd_kfold) for each predictor set
elpd_vec_pred_ao <- sapply(mods_pred_ao, function(x) x$estimates["elpd_kfold", "Estimate"])

# 4) Build the summary table
tab_predictors_ao <- as.data.frame(lc_pred_ao) %>%
  rownames_to_column("model") %>%
  mutate(
    Predictor_combination = case_when(
      model == "m15"  ~ "1 + 5  (visitation + functional-group richness)",
      model == "m135" ~ "1 + 3 + 5 (visitation + morphospecies + functional-group)",
      model == "m13"  ~ "1 + 3  (visitation + morphospecies richness)",
      model == "m35"  ~ "3 + 5  (morphospecies + functional-group richness)",
      model == "m1"   ~ "1  (visitation only)",
      model == "m3"   ~ "3  (morphospecies richness only)",
      model == "m5"   ~ "5  (functional-group richness only)",
      model == "mnull"~ "null (intercept only)",
      TRUE ~ model
    ),
    ELPD   = elpd_vec_pred_ao[model],                         # absolute elpd_kfold
    z_ELPD = if_else(se_diff > 0, elpd_diff / se_diff, NA_real_)  # best model = NA
  ) %>%
  select(Predictor_combination, ELPD, elpd_diff, se_diff, z_ELPD)

tab_predictors_ao

write.csv(tab_predictors_ao,
          "tables/tab_models_ao_variations.csv",
          row.names = FALSE)



###* GEITONOGAMY
kfold_go_13 <- readRDS("brms_models/kfold_bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt.rds")
kfold_go_null <- readRDS("brms_models/kfold_bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt.rds")
kfold_go_135 <- readRDS("brms_models/kfold_bayesian_go_mean_corrected_scaled_135_2_zoib_k_opt.rds")
kfold_go_1 <- readRDS("brms_models/kfold_bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt.rds")
kfold_go_3 <- readRDS("brms_models/kfold_bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt.rds")
kfold_go_5 <- readRDS("brms_models/kfold_bayesian_go_mean_corrected_scaled_5_2_zoib_k_opt.rds")
kfold_go_35 <- readRDS("brms_models/kfold_bayesian_go_mean_corrected_scaled_35_2_zoib_k_opt.rds")
kfold_go_15 <- readRDS("brms_models/kfold_bayesian_go_mean_corrected_scaled_15_2_zoib_k_opt.rds")

loo_compare(kfold_go_13,
            kfold_go_135,
            kfold_go_1,
            kfold_go_3,
            kfold_go_35,
            kfold_go_5,
            kfold_go_15
            #kfold_go_null
            )

# 1) Put all K-fold objects (EXCEPT the null model) into a named list ----
# Replace the object names here with your actual k-fold objects
mods_pred_go <- list(
  m15  = kfold_go_15,    # visitation (1) + functional-group richness (5)
  m135 = kfold_go_135,   # 1 + 3 + 5
  m13  = kfold_go_13,    # 1 + 3
  m35  = kfold_go_35,    # 3 + 5
  m1   = kfold_go_1,     # 1 only
  m3   = kfold_go_3,     # 3 only
  m5   = kfold_go_5      # 5 only
  #mnull= kfold_go_null
)

loo_compare(mods_pred_go)

# 2) Compare all candidate predictor sets
lc_pred_go <- loo_compare(mods_pred_go)   # works on k-fold objects, uses elpd_kfold
lc_pred_go

# 3) Get absolute ELPD (elpd_kfold) for each predictor set
elpd_vec_pred_go <- sapply(mods_pred_go, function(x) x$estimates["elpd_kfold", "Estimate"])

# 4) Build the summary table
tab_predictors_go <- as.data.frame(lc_pred_go) %>%
  rownames_to_column("model") %>%
  mutate(
    Predictor_combination = case_when(
      model == "m15"  ~ "1 + 5  (visitation + functional-group richness)",
      model == "m135" ~ "1 + 3 + 5 (visitation + morphospecies + functional-group)",
      model == "m13"  ~ "1 + 3  (visitation + morphospecies richness)",
      model == "m35"  ~ "3 + 5  (morphospecies + functional-group richness)",
      model == "m1"   ~ "1  (visitation only)",
      model == "m3"   ~ "3  (morphospecies richness only)",
      model == "m5"   ~ "5  (functional-group richness only)",
      #model == "mnull"~ "null (intercept only)",
      TRUE ~ model
    ),
    ELPD   = elpd_vec_pred_go[model],                         # absolute elpd_kfold
    z_ELPD = if_else(se_diff > 0, elpd_diff / se_diff, NA_real_)  # best model = NA
  ) %>%
  select(Predictor_combination, ELPD, elpd_diff, se_diff, z_ELPD)

tab_predictors_go

write.csv(tab_predictors_go,
          "tables/tab_models_go_variations.csv",
          row.names = FALSE)











packageVersion("emmeans")
