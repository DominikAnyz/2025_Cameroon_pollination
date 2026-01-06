###* This script is to test the effect of the predictors "total.morphospecies", 
###* "total.functional.groups" and "visited.flowers.per.minute"

###* Loading the necessary packages
pacman::p_load(tidyverse, brms, ordbetareg, scales, loo)

set.seed(1234)

source("scripts/05.Setup_for_visitation_indices.R")

final.table$elevation <- as.factor(final.table$elevation)

# Recode incorrect elevation values
final.table <- final.table %>%
  mutate(elevation = recode(elevation,
                            `3500` = "3400",
                            `4000` = "3800"))


###* For each response index, I fit a set of four candidate models that reflect
###* the observed distribution of the response (morpho): Poisson, zero-inflated
###* Poisson, negative binomial, and zero-inflated negative binomial.
###*
###* Across candidates, the model structure is identical; only these elements 
###* vary:
###*   (i) model name (po / zipo / nb / zinb),
###*   (ii) whether a zero-inflation component is included,
###*   (iii) the family (poisson / zero_inflated_poisson / negbinomial /
###*        zero_inflated_negbinomial),
###*   (iv) the prior on the zero-inflation intercept (for ZI families).
###*
###* All fitted models are saved as .rds files. The paths are kept in the script
###* so that models can be loaded quickly without refitting. (This is helpful
###* because model fitting and, where needed, k-fold cross-validation are 
###* time-consuming.)
###*
###* Model comparison is done using LOO-CV by default. If LOO diagnostics indicate
###* unreliable importance sampling (e.g., Pareto k > 0.7), I re-fit comparisons
###* using k-fold CV. (Saved objects are also kept for transparency and reuse.)

###* MORPHOSPECIES RICHNESS
m_brm_model_po <- readRDS("brms_models/m_brm_model_po.rds")
m_brm_model_zipo <- readRDS("brms_models/m_brm_model_zipo.rds")
m_brm_model_nb <- readRDS("brms_models/m_brm_model_nb.rds")
m_brm_model_zinb <- readRDS("brms_models/m_brm_model_zinb.rds")

loo(m_brm_model_po, m_brm_model_zipo, m_brm_model_nb, m_brm_model_zinb)

###* clearly, the negative binomial model is the best fit
###* The code for it is:
m_brm_model_nb <- brm(
  formula = bf(total.morpho ~ elevation  + (1|plant.species)),
  data = final.table,
  family = negbinomial(),
  prior = c(
    set_prior("normal(0, 2)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd")
  ),
  chains = 6, cores = 6,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  iter = 5000, warmup = 2000,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)

###* Null model to check if model explains more than random
m_brm_null_nb <- brm(
  formula = bf(total.morpho ~ 1  + (1|plant.species)),
  data = final.table,
  family = negbinomial(),
  prior = c(
    set_prior("student_t(3, 0, 5)", class = "sd")
  ),
  chains = 6, cores = 6,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  iter = 8000, warmup = 3000,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)

saveRDS(m_brm_null_nb, file = "brms_models/m_brm_null_nb.rds")
m_brm_null_nb <- readRDS("brms_models/m_brm_null_nb.rds")

loo(m_brm_null_nb, m_brm_model_nb)
bayes_R2(m_brm_null_nb)
bayes_R2(m_brm_model_nb)

###* Clearly, the model explains more than the null model
###* For clarification, how does loo work?
###* 
###* For every observation, loo asks "how well can the model predict this point,
###* if it had been left out when fitting"?
###* 
###* The results is expected log predictive density (elpd_loo)
###* this is the sum of log predictive densities across all points under loo CV
###* SE is the uncertainty for the estinate
###* -> the higher the eldp_loo, the better predictive accuracy
###* 
###* p_loo / effective number of parameters (model complexity penalty)
###* 
###* elpd_diff, se_diff
###* if Δelpd > ~2 × se_diff, then evidence for the better model is strong
###* 
###* pareto k diagnostics
###* Each observation has a “k” diagnostic telling whether the LOO estimate is reliable
###* k < 0.7 = good; 0.7–1 = problematic; >1 = unreliable.
###* 
###* 
###* So now that we know which model is the best and that it is more reliable
###* than the null, can it be trusted?

summary(m_brm_model_nb)
###* RHAT should be 1 (good if <1.01)
###* bulk_ESS and tail_ESS should be reasonable large (ideally > 1000)

###* Do simulated datasets from posterioirs look like real data?
pp_check(m_brm_model_nb, type = "dens_overlay")
pp_check(m_brm_model_nb, type = "hist")
pp_check(m_brm_model_nb, type = "stat", stat = "mean")
pp_check(m_brm_model_nb, type = "stat", stat = "sd")

###* compute residuals and check patterns unexplained by model
resids <- residuals(m_brm_model_nb, type = "pearson")
plot(resids)      # residual distribution
hist(resids)      # check symmetry/variance

###* model is good fit! 
###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* 

###* FUNCTIONAL GROUP RICHNESS
f_brm_model_po <- readRDS("brms_models/f_brm_model_po.rds")
f_brm_model_zipo <- readRDS("brms_models/f_brm_model_zipo.rds")
f_brm_model_nb <- readRDS("brms_models/f_brm_model_nb.rds")
f_brm_model_zinb <- readRDS("brms_models/f_brm_model_zinb.rds")

loo(f_brm_model_po, f_brm_model_zipo, f_brm_model_zinb)

f_brm_model_nb <- brm(
  formula = bf(total.func ~ elevation  + (1|plant.species)),
  data = final.table,
  family = negbinomial(),
  prior = c(
    set_prior("normal(0, 2)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd")
  ),
  chains = 5, cores = 5,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  iter = 5000, warmup = 2000,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 123
)

saveRDS(f_brm_model_nb, file = "brms_models/f_brm_model_nb.rds")

###* Null model to check if model explains more than random

f_brm_null_po <- brm(
  formula = bf(total.func ~ 1  + (1|plant.species)),
  data = final.table,
  family = poisson(),
  prior = c(
    set_prior("student_t(3, 0, 5)", class = "sd")
  ),
  chains = 5, cores = 5,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  iter = 5000, warmup = 2000,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)

saveRDS(f_brm_null_po, file = "brms_models/f_brm_null_po.rds")
f_brm_null_po <- readRDS("brms_models/f_brm_null_po.rds")

loo(f_brm_null_po, f_brm_model_po)
bayes_R2(f_brm_null_po)
bayes_R2(f_brm_model_po)

###* Validation
summary(f_brm_model_po)

###* Do simulated datasets from posterioirs look like real data?
pp_check(f_brm_model_po, type = "dens_overlay")
pp_check(f_brm_model_po, type = "hist")
pp_check(f_brm_model_po, type = "stat", stat = "mean")
pp_check(f_brm_model_po, type = "stat", stat = "sd")

###* compute residuals and check patterns unexplained by model
resids <- residuals(f_brm_model_po, type = "pearson")
plot(resids)      # residual distribution
hist(resids)      # check symmetry/variance

bayes_R2(f_brm_model_po)


###* VISITATION FREQUENCY
final.table <- final.table %>%
  ungroup() %>%  # Ensure no grouping is active
  mutate(visited.0.1.scaled = rescale(visited.flowers.per.minute, to = c(0, 1), na.rm = TRUE))

vis_ord_model <- ordbetareg(
  formula = visited.0.1.scaled ~ elevation  + (1|plant.species),
  data = final.table,
  chains = 5, iter = 5000, warmup = 2000, cores = 5,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)

saveRDS(vis_ord_model, file = "brms_models/vis_ord_model.rds")
vis_ord_model <- readRDS("brms_models/vis_ord_model.rds")

final.table$intercept_only <- 1

vis_ord_null <- ordbetareg(
  formula = visited.0.1.scaled ~ intercept_only + (1|plant.species),
  data = final.table,
  chains = 5, iter = 5000, warmup = 2000, cores = 5,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)

saveRDS(vis_ord_null, file = "brms_models/vis_ord_null.rds")
vis_ord_null <- readRDS("brms_models/vis_ord_null.rds")

loo(vis_ord_model, vis_ord_null)
bayes_R2(vis_ord_model)
bayes_R2(vis_ord_null)

summary(vis_ord_model)

###* Do simulated datasets from posterioirs look like real data?
pp_check(vis_ord_model, type = "dens_overlay")
pp_check(vis_ord_model, type = "hist")
pp_check(vis_ord_model, type = "stat", stat = "mean")
pp_check(vis_ord_model, type = "stat", stat = "sd")

###* compute residuals and check patterns unexplained by model
resids <- residuals(vis_ord_model, type = "pearson")
plot(resids)      # residual distribution
hist(resids)      # check symmetry/variance

bayes_R2(vis_ord_model)
