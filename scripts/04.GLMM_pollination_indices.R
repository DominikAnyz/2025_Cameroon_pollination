###* This script contains all of the code used to generate the bayesian glmm
###* models for the pollination indices, aquired from the hand-pollination 
###* experiment

pacman::p_load(tidyverse, loo, ordbetareg, brms)

source("scripts/03.Setup_for_pollination_indices.R")

###* For each response index, I fit a set of four candidate models that reflect
###* the observed distribution of the response (seedset): Poisson, zero-inflated
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
###*
###* MODEL FOR CONTROL SEEDSET
###* 
###* ZINB
c_brm_model_nb_new <- brm(
  formula = bf(seedset ~ elevation + (1|species) + (1|plant.id)#,
               #zi ~ elevation
               ),
  data = c.index,
  family = negbinomial(),
  prior = c(
    set_prior("normal(0, 1.5)", class = "b"),
    #set_prior("normal(0, 2)", class = "b", dpar = "zi"),
    set_prior("student_t(3, 0, 5)", class = "sd")
  ),
  chains = 5, cores = 5,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  iter = 2000, warmup = 1000,
  control = list(max_treedepth = 15, adapt_delta = 0.95),
  seed = 1234
)

###* Save models with different distributions
# saveRDS(c_brm_model_zinb_new, "brms_models/c_brm_model_zinb_new.rds")
# saveRDS(c_brm_model_po_new, "brms_models/c_brm_model_po_new.rds")
# saveRDS(c_brm_model_zipo_new, "brms_models/c_brm_model_zipo_new.rds")
# saveRDS(c_brm_model_nb_new, "brms_models/c_brm_model_nb_new.rds")

###* Load models with different distributions
c_brm_model_po_new <- readRDS("brms_models/c_brm_model_po_new.rds")
c_brm_model_zipo_new <- readRDS("brms_models/c_brm_model_zipo_new.rds")
c_brm_model_nb_new <- readRDS("brms_models/c_brm_model_nb_new.rds")
c_brm_model_zinb_new <- readRDS("brms_models/c_brm_model_zinb_new.rds")

###* Check loo's of models
loo(c_brm_model_po_new, c_brm_model_zipo_new, c_brm_model_nb_new, c_brm_model_zinb_new)
loo(c_brm_model_zinb_new)
###* zinb seems to be best, but there are pareto warings
###* 
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
###* Getting rid of pareto warnings by using kfold
###* 
# c_brm_model_zinb_new_k <- kfold(c_brm_model_zinb_new, K = 5, cores = 5)
# saveRDS(c_brm_model_zinb_new_k, "brms_models/c_brm_model_zinb_new_k.rds")
# c_brm_model_po_new_k <- kfold(c_brm_model_po_new, K = 5, cores = 5)
# saveRDS(c_brm_model_po_new_k, "brms_models/c_brm_model_po_new_k.rds")
# c_brm_model_zipo_new_k <- kfold(c_brm_model_zipo_new, K = 5, cores = 5)
# saveRDS(c_brm_model_zipo_new_k, "brms_models/c_brm_model_zipo_new_k.rds")
# c_brm_model_nb_new_k <- kfold(c_brm_model_nb_new, K = 5, cores = 5)
# saveRDS(c_brm_model_nb_new_k, "brms_models/c_brm_model_nb_new_k.rds")

###* Prepared loo images for models
c_brm_model_zinb_new_k <- readRDS("brms_models/c_brm_model_zinb_new_k.rds")
c_brm_model_po_new_k <- readRDS("brms_models/c_brm_model_po_new_k.rds")
c_brm_model_zipo_new_k <- readRDS("brms_models/c_brm_model_zipo_new_k.rds")
c_brm_model_nb_new_k <- readRDS("brms_models/c_brm_model_nb_new_k.rds")

###* Comparing models after kfolding
loo_compare(c_brm_model_po_new_k, c_brm_model_zipo_new_k, c_brm_model_nb_new_k, c_brm_model_zinb_new_k)
###* ZINB model comes out on top

###* Creating null model to my zinb model to see, whether full model shows
###* better support
###* NULL
c_brm_null_zinb_new <- brm(
  formula = bf(seedset ~ 1 + (1|species) + (1|plant.id),
               zi ~ elevation),
  data = c.index,
  family = zero_inflated_negbinomial(),
  prior = c(
    set_prior("student_t(3, 0, 5)", class = "sd"),
    set_prior("normal(0, 2)", class = "b", dpar = "zi")
  ),
  chains = 5, cores = 5,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  iter = 2000, warmup = 1000,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)

###* Save model
saveRDS(c_brm_null_zinb_new, file = "brms_models/c_brm_null_zinb_new.rds")
c_brm_null_zinb_new <- readRDS("brms_models/c_brm_null_zinb_new.rds")
loo(c_brm_null_zinb_new)
###* Also problem with pareto

###* Fix pareto problem in null
c_brm_null_zinb_new_k <- kfold(c_brm_null_zinb_new, K = 5, cores = 5)
#saveRDS(c_brm_null_zinb_new_k, "brms_models/c_brm_null_zinb_new_k.rds")

###* Fit kfold fixed model
c_brm_null_zinb_new_k <- readRDS("brms_models/c_brm_null_zinb_new_k.rds")

c_brm_model_zinb_new_k
c_brm_null_zinb_new_k

loo_compare(c_brm_null_zinb_new_k, c_brm_model_zinb_new_k)

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
###* So now that we know which model is the best and that it is more reliable
###* than the null, can it be trusted?
###* 
c_brm_model_zinb_new <- readRDS("brms_models/c_brm_model_zinb_new.rds")
summary(c_brm_model_zinb_new)

pp_check(c_brm_model_zinb_new, type = "dens_overlay")
pp_check(c_brm_model_zinb_new, type = "hist")
pp_check(c_brm_model_zinb_new, type = "stat", stat = "mean")
pp_check(c_brm_model_zinb_new, type = "stat", stat = "sd")

###* compute residuals and check patterns unexplained by model
resids <- residuals(c_brm_model_zinb_new, type = "pearson")
plot(resids)      
hist(resids)

bayes_R2(c_brm_model_zinb_new)
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
###* MODEL FOR POLLEN LIMITATION

###* Exclude Hypericum from highest eelvation, as it produced no seeds whatsover
c.index <- c.index %>%
  filter(species != "Hypericum r" | elevation != 3800)

###* Use ordbetareg package, since the data ahve a zero-one inflated beta distribution
pl_ord_model_new <- ordbetareg(
  formula = PL.index ~ elevation + (1|species) + (1|plant.id),
  data = c.index,
  chains = 5, cores = 5,
  iter = 2000, warmup = 1000,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  control = list(adapt_delta = 0.99),
  seed = 1234
)

###* Save model and check loo
#saveRDS(pl_ord_model_new, file = "brms_models/pl_ord_model_new.rds")
pl_ord_model_new <- readRDS("brms_models/pl_ord_model_new.rds")
loo(pl_ord_model_new)

###* K-fold, since pareto warnings
pl_ord_model_new_k <- kfold(pl_ord_model_new, K = 5, cores = 5)
#saveRDS(pl_ord_model_new_k, "brms_models/pl_ord_model_new_k.rds")
pl_ord_model_new_k <- readRDS("brms_models/pl_ord_model_new_k.rds")

###* Add "intercept = 1" to dataset, in order to be able to run null model. 
###* 
###* ordbetareg requires an explicit predictor term; adding a constant column
###* provides an intercept-only fixed effect while keeping the same random 
###* effects.
c.index$intercept_only <- 1

###* Run null model
pl_ord_null_new <- ordbetareg(
  formula = PL.index ~ intercept_only + (1|species) + (1|plant.id),
  data = c.index,
  chains = 5, cores = 5,
  iter = 2000, warmup = 1000, 
  save_pars = save_pars(latent = TRUE, all = TRUE), 
  control = list(adapt_delta = 0.99),
  seed = 1234
)

#saveRDS(pl_ord_null_new, file = "brms_models/pl_ord_null_new.rds")
pl_ord_null_new <- readRDS("brms_models/pl_ord_null_new.rds")

###* K-fold null because of parteto warnings
pl_ord_null_new_k <- kfold(pl_ord_null_new, K = 5, cores = 5)
#saveRDS(pl_ord_null_new_k, "brms_models/pl_ord_null_new_k.rds")
pl_ord_null_new_k <- readRDS("brms_models/pl_ord_null_new_k.rds")

###* Compare loo§s of models
loo_compare(pl_ord_null_new_k, pl_ord_model_new_k)
###* yes, model is better than null

###* MODEL OUTPUT
pl_ord_model_new <- readRDS("brms_models/pl_ord_model_new.rds")
summary(pl_ord_model_new)

pp_check(pl_ord_model_new, type = "dens_overlay")
pp_check(pl_ord_model_new, type = "hist")
pp_check(pl_ord_model_new, type = "stat", stat = "mean")
pp_check(pl_ord_model_new, type = "stat", stat = "sd")

###* compute residuals and check patterns unexplained by model
resids <- residuals(pl_ord_model_new, type = "pearson")
plot(resids)      # residual distribution
hist(resids)

bayes_R2(pl_ord_model_new)
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
###* MODEL FOR AO INDEX
###* 
###* Run all model variations
ao_brm_model_zipo_new <- brm(
  formula = bf(index_trans ~ elevation + (1|species) + (1|plant.id),
               zi ~ elevation
               ),
  data = ao.index,
  family = zero_inflated_poisson(),
  prior = c(
    set_prior("normal(0, 1.5)", class = "b"),
    set_prior("normal(0, 2)", class = "b", dpar = "zi"),
    set_prior("student_t(3, 0, 5)", class = "sd")
  ),
  chains = 5, cores = 5,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  iter = 2000, warmup = 1000,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)


###* Save models with different distributions
# saveRDS(ao_brm_model_po_new, "brms_models/ao_brm_model_po_new.rds")
# saveRDS(ao_brm_model_zipo_new, "brms_models/ao_brm_model_zipo_new.rds")
# saveRDS(ao_brm_model_nb_new, "brms_models/ao_brm_model_nb_new.rds")
# saveRDS(ao_brm_model_zinb_new, "brms_models/ao_brm_model_zinb_new.rds")

###* Load models with different distributions
ao_brm_model_po_new <- readRDS("brms_models/ao_brm_model_po_new.rds")
ao_brm_model_zipo_new <- readRDS("brms_models/ao_brm_model_zipo_new.rds")
ao_brm_model_nb_new <- readRDS("brms_models/ao_brm_model_nb_new.rds")
ao_brm_model_zinb_new <- readRDS("brms_models/ao_brm_model_zinb_new.rds")

###* Check loo's of models
loo(ao_brm_model_po_new, ao_brm_model_zipo_new, ao_brm_model_nb_new, ao_brm_model_zinb_new)
loo(ao_brm_model_zinb_new)
###* zinb seems to be best, but there are pareto warings

###* Getting rid of pareto warnings by using kfold
###* 
# ao_brm_model_po_new_k <- kfold(ao_brm_model_po_new, K = 2, cores = 5)
# ao_brm_model_po_new_k
# saveRDS(ao_brm_model_po_new_k, "brms_models/ao_brm_model_po_new_k.rds")
# ao_brm_model_zipo_new_k <- kfold(ao_brm_model_zipo_new, K = 3, cores = 5)
# ao_brm_model_zipo_new_k
# saveRDS(ao_brm_model_zipo_new_k, "brms_models/ao_brm_model_zipo_new_k.rds")
# ao_brm_model_nb_new_k <- kfold(ao_brm_model_nb_new, K = 3, cores = 5)
# ao_brm_model_nb_new_k
# saveRDS(ao_brm_model_nb_new_k, "brms_models/ao_brm_model_nb_new_k.rds")
# ao_brm_model_zinb_new_k <- kfold(ao_brm_model_zinb_new, K = 3, cores = 5)
# ao_brm_model_zinb_new_k
# saveRDS(ao_brm_model_zinb_new_k, "brms_models/ao_brm_model_zinb_new_k.rds")


###* Prepared loo images for models
ao_brm_model_zinb_new_k <- readRDS("brms_models/ao_brm_model_zinb_new_k.rds")
ao_brm_model_po_new_k <- readRDS("brms_models/ao_brm_model_po_new_k.rds")
ao_brm_model_zipo_new_k <- readRDS("brms_models/ao_brm_model_zipo_new_k.rds")
ao_brm_model_nb_new_k <- readRDS("brms_models/ao_brm_model_nb_new_k.rds")

###* Comopare loo for models
loo_compare(ao_brm_model_zinb_new_k, ao_brm_model_po_new_k, ao_brm_model_zipo_new_k, ao_brm_model_nb_new_k)
###* ZINB is the best
 
###* Run  null model for zinb
ao_brm_null_zinb <- brm(
  formula = bf(index_trans ~ 1 + (1|species) + (1|plant.id),
               zi ~ elevation),
  data = ao.index,
  family = zero_inflated_negbinomial(),
  prior = c(
    set_prior("student_t(3, 0, 5)", class = "sd"),
    set_prior("normal(0, 2)", class = "b", dpar = "zi")
  ),
  chains = 5, cores = 5,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  iter = 5000, warmup = 2000,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)

#saveRDS(ao_brm_null_zinb, file = "brms_models/ao_brm_null_zinb.rds")
ao_brm_null_zinb <- readRDS("brms_models/ao_brm_null_zinb.rds")

###* Check loo
loo(ao_brm_null_zinb)

###* K-fold model to get rid of pareto warnings
ao_brm_null_zinb_k <- kfold(ao_brm_null_zinb, K = 2, cores = 5)
ao_brm_null_zinb_k
#saveRDS(ao_brm_null_zinb_k, "brms_models/ao_brm_null_zinb_k.rds")

###* Compare full and null model
loo_compare(ao_brm_null_zinb_k, ao_brm_model_zinb_new_k)
###* model betterthan null
###* 
###* 
###* 
###* 
###* MODEL OUTPUT
ao_brm_model_zinb <- readRDS("brms_models/ao_brm_model_zinb.rds")
summary(ao_brm_model_zinb)

###* PP check
pp_check(ao_brm_model_zinb, type = "dens_overlay")
pp_check(ao_brm_model_zinb, type = "hist")
pp_check(ao_brm_model_zinb, type = "stat", stat = "mean")
pp_check(ao_brm_model_zinb, type = "stat", stat = "sd")

###* compute residuals and check patterns unexplained by model
resids <- residuals(ao_brm_model_zinb, type = "pearson")
plot(resids)      # residual distribution
hist(resids)

bayes_R2(ao_brm_model_zinb)
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
###* MODEL FOR GEITONOGAMY

###* Create versions of geitonogamy model
go_brm_model_zinb_new <- brm(
  formula = bf(index_trans ~ elevation + (1|species) + (1|plant.id),
               zi ~ elevation
  ),
  data = go.index,
  family = zero_inflated_negbinomial(),
  prior = c(
    set_prior("normal(0, 1.5)", class = "b"),
    set_prior("normal(0, 2)", class = "b", dpar = "zi"),
    set_prior("student_t(3, 0, 5)", class = "sd")
  ),
  chains = 5, cores = 5,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  iter = 2000, warmup = 1000,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)

###* Save models with different distributions
# saveRDS(go_brm_model_po_new, "brms_models/go_brm_model_po_new.rds")
# saveRDS(go_brm_model_zipo_new, "brms_models/go_brm_model_zipo_new.rds")
# saveRDS(go_brm_model_nb_new, "brms_models/go_brm_model_nb_new.rds")
# saveRDS(go_brm_model_zinb_new, "brms_models/go_brm_model_zinb_new.rds")

###* Load models with different distributions
go_brm_model_po_new <- readRDS("brms_models/go_brm_model_po_new.rds")
go_brm_model_zipo_new <- readRDS("brms_models/go_brm_model_zipo_new.rds")
go_brm_model_nb_new <- readRDS("brms_models/go_brm_model_nb_new.rds")
go_brm_model_zinb_new <- readRDS("brms_models/go_brm_model_zinb_new.rds")

###* Check loo's of models
loo(go_brm_model_po_new, go_brm_model_zipo_new, go_brm_model_nb_new, go_brm_model_zinb_new)
loo(go_brm_model_zinb_new)
###* zinb seems to be best, but there are pareto warings

###* Getting rid of pareto warnings by using kfold
###* 
# go_brm_model_po_new_k <- kfold(go_brm_model_po_new, K = 2, cores = 5)
# go_brm_model_po_new_k
# saveRDS(go_brm_model_po_new_k, "brms_models/go_brm_model_po_new_k.rds")
# go_brm_model_zipo_new_k <- kfold(go_brm_model_zipo_new, K = 3, cores = 5)
# go_brm_model_zipo_new_k
# saveRDS(go_brm_model_zipo_new_k, "brms_models/go_brm_model_zipo_new_k.rds")
# go_brm_model_nb_new_k <- kfold(go_brm_model_nb_new, K = 3, cores = 5)
# go_brm_model_nb_new_k
# saveRDS(go_brm_model_nb_new_k, "brms_models/go_brm_model_nb_new_k.rds")
# go_brm_model_zinb_new_k <- kfold(go_brm_model_zinb_new, K = 3, cores = 5)
# go_brm_model_zinb_new_k
# saveRDS(go_brm_model_zinb_new_k, "brms_models/go_brm_model_zinb_new_k.rds")


###* Prepared loo images for models
go_brm_model_zinb_new_k <- readRDS("brms_models/go_brm_model_zinb_new_k.rds")
go_brm_model_po_new_k <- readRDS("brms_models/go_brm_model_po_new_k.rds")
go_brm_model_zipo_new_k <- readRDS("brms_models/go_brm_model_zipo_new_k.rds")
go_brm_model_nb_new_k <- readRDS("brms_models/go_brm_model_nb_new_k.rds")

loo_compare(go_brm_model_po_new_k, go_brm_model_zipo_new_k, go_brm_model_nb_new_k, go_brm_model_zinb_new_k)

###* Running null model
go_brm_null_zinb_new <- brm(
  formula = bf(index_trans ~ 1 + (1|species) + (1|plant.id),
               zi ~ elevation),
  data = go.index,
  family = zero_inflated_negbinomial(),
  prior = c(
    set_prior("student_t(3, 0, 5)", class = "sd"),
    set_prior("normal(0, 2)", class = "b", dpar = "zi")
  ),
  chains = 5, cores = 5,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  iter = 2000, warmup = 1000,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)

saveRDS(go_brm_null_zinb_new, "brms_models/go_brm_null_zinb_new.rds")
go_brm_null_zinb_new <- readRDS("brms_models/go_brm_null_zinb_new.rds")

loo(go_brm_null_zinb)

###* K-fold model to get rid of pareto warnings
go_brm_null_zinb_new_k <- kfold(go_brm_null_zinb_new, K = 2, cores = 5)
go_brm_null_zinb_new_k
#saveRDS(go_brm_null_zinb_k, "brms_models/go_brm_null_zinb_k.rds")
go_brm_null_zinb_new_k <- readRDS("brms_models/go_brm_null_zinb_new_k.rds")

loo_compare(go_brm_null_zinb_new_k, go_brm_model_zinb_new_k)
###* model better than null
###* 
###* 