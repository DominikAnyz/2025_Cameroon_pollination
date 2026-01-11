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
pacman::p_load(dplyr, loo, emmeans)
###* 
###* 
###* 
###* 
###* Table S1
c_brm_model_zinb_new <- readRDS("brms_models/c_brm_model_zinb_new.rds")
pl_ord_model_new <- readRDS("brms_models/pl_ord_model_new.rds")
ao_brm_model_zinb_new <- readRDS("brms_models/ao_brm_model_zinb_new.rds")
go_brm_model_zinb_new <- readRDS("brms_models/go_brm_model_zinb_new.rds")
vis_ord_model <- readRDS("brms_models/vis_ord_model.rds")
m_brm_model_nb <- readRDS("brms_models/m_brm_model_nb.rds")
f_brm_model_po <- readRDS("brms_models/f_brm_model_po.rds")

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
###* 
###* 
###* 
###* TABLE S2
###* 
###* 
###* 
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

loo_compare(c_brm_model_po_new_k, c_brm_model_zipo_new_k, c_brm_model_nb_new_k, c_brm_model_zinb_new_k)


###* Autogamy
ao_brm_model_zinb_new_k <- readRDS("brms_models/ao_brm_model_zinb_new_k.rds")
ao_brm_model_po_new_k <- readRDS("brms_models/ao_brm_model_po_new_k.rds")
ao_brm_model_zipo_new_k <- readRDS("brms_models/ao_brm_model_zipo_new_k.rds")
ao_brm_model_nb_new_k <- readRDS("brms_models/ao_brm_model_nb_new_k.rds")

loo_compare(ao_brm_model_zinb_new_k, ao_brm_model_po_new_k, ao_brm_model_zipo_new_k, ao_brm_model_nb_new_k)


###* Geitonogamy
go_brm_model_zinb_new_k <- readRDS("brms_models/go_brm_model_zinb_new_k.rds")
go_brm_model_po_new_k <- readRDS("brms_models/go_brm_model_po_new_k.rds")
go_brm_model_zipo_new_k <- readRDS("brms_models/go_brm_model_zipo_new_k.rds")
go_brm_model_nb_new_k <- readRDS("brms_models/go_brm_model_nb_new_k.rds")

loo(go_brm_model_po_new, go_brm_model_zipo_new, go_brm_model_nb_new, go_brm_model_zinb_new)


###* MORPHOSPECIES RICHNESS
m_brm_model_po <- readRDS("brms_models/m_brm_model_po.rds")
m_brm_model_zipo <- readRDS("brms_models/m_brm_model_zipo.rds")
m_brm_model_nb <- readRDS("brms_models/m_brm_model_nb.rds")
m_brm_model_zinb <- readRDS("brms_models/m_brm_model_zinb.rds")

loo(m_brm_model_po, m_brm_model_zipo, m_brm_model_nb, m_brm_model_zinb)


###* FUNCTIONAL GROUP RICHNESS
f_brm_model_po <- readRDS("brms_models/f_brm_model_po.rds")
f_brm_model_zipo <- readRDS("brms_models/f_brm_model_zipo.rds")
f_brm_model_nb <- readRDS("brms_models/f_brm_model_nb.rds")
f_brm_model_zinb <- readRDS("brms_models/f_brm_model_zinb.rds")

loo(f_brm_model_po, f_brm_model_zipo, f_brm_model_zinb)
###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* TABLE S3
###* 
###* 
###* 
###* 
###* 
###* 
###* Natural seed set
kfold_bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt <- readRDS("brms_models/kfold_bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt <- readRDS("brms_models/kfold_bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt <- readRDS("brms_models/kfold_bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt <- readRDS("brms_models/kfold_bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt.rds")

###* Comparing existing kfolded object to see, which model is best
loo_compare(kfold_bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt,
            kfold_bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt,
            kfold_bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt,
            kfold_bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt)


###* Pollen limitation
kfold_bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2 <- readRDS("brms_models/kfold_bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2.rds")
kfold_bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2 <- readRDS("brms_models/kfold_bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2.rds")
kfold_bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2 <- readRDS("brms_models/kfold_bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2.rds")
kfold_bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2 <- readRDS("brms_models/kfold_bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2.rds")

###* COomparison of individual kfolded objects
loo_compare(kfold_bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2,
            kfold_bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2,
            kfold_bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2,
            kfold_bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2)

###* Autogamy
kfold_bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt <- readRDS("brms_models/kfold_bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt <- readRDS("brms_models/kfold_bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt <- readRDS("brms_models/kfold_bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt <- readRDS("brms_models/kfold_bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt.rds")

loo_compare(kfold_bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt,
            kfold_bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt,
            kfold_bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt,
            kfold_bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt)

###* Geitonogamy
kfold_bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt <- readRDS("brms_models/kfold_bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt <- readRDS("brms_models/kfold_bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt <- readRDS("brms_models/kfold_bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt <- readRDS("brms_models/kfold_bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt.rds")


loo_compare(kfold_bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt,
            kfold_bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt,
            kfold_bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt,
            kfold_bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt)