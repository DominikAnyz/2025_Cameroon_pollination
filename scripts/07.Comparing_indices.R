###* This script serves as running and interpreting the results from the best 
###* models created with the data which we have.
pacman::p_load(brms, dplyr, ggplot2, scales, patchwork, loo)

options(brms.backend = "cmdstanr")

source("scripts/06. Setup for index comparison.R")

###* Renaming for easier interpretation
c.pl.final.table.4 <- c.pl.final.table.4 %>%
  rename(
    z_vis_mean = mean.visitors.scaled,
    z_vis_sd   = sd.visitors.scaled,
    z_flow_mean = mean.visited.flowers.scaled,
    z_flow_sd   = sd.visited.flowers.scaled,
    z_morpho_mean = mean.morpho.scaled,
    z_morpho_sd   = sd.morpho.scaled,
    z_func_mean   = mean.func.scaled,
    z_func_sd     = sd.func.scaled
  )

###* tiny floor to avoid zero SDs
eps <- 1e-8 

###* Renaming to make as clean as possible
c.pl.final.table.4 <- c.pl.final.table.4 %>%
  mutate(
    # --- means on the z-scale (already computed) ---
    x = z_flow_mean,   # visited flowers per minute (scaled) -- or use z_vis_mean if that's your x
    y = z_morpho_mean,
    z = z_func_mean,
    # --- SEs of the means on the same scale (SD / sqrt(n)) ---
    sx = pmax(z_flow_sd   / sqrt(pmax(n_reps,  1L)), eps),
    sy = pmax(z_morpho_sd  / sqrt(pmax(n_reps, 1L)), eps),
    sz = pmax(z_func_sd   / sqrt(pmax(n_reps,  1L)), eps)
    
  )


View(c.pl.final.table.4)

###* Generating model to analyze whether visitation has an effect on pollination
###* indices
###* 
###* Generating formula
###* 

form_all <- bf(
  mean_seedset_round | weights(seedset_weight_12) ~
    me(x, sx) +
    me(y, sy) +
    me(z, sz) +
    #1 +
    (1|species),
  zi ~ 1
)

###* Create a weakly informative prior for the intercept:
###* For the ZINB model, the mean parameter is on the log scale.
###* We initialize the intercept prior around log(mean(seed set)), using the
###* overall mean of the response (plus a small constant to avoid log(0)).
mu0 <- log(mean(c.pl.final.table.4$mean_seedset_round, na.rm = TRUE) + 1e-8)
mu0

pri_all <- c(
  prior(normal(3.688879, 1.5), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(student_t(3, 0, 2.5), class = "sd", group = "species")
)

bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt <- brm(
  form_all, 
  data = c.pl.final.table.4,
  family = zero_inflated_negbinomial(), 
  prior = pri_all,
  chains = 5, cores = 5, 
  iter = 2000, warmup = 1000,
  threads = threading(4),
  save_pars = save_pars(all = TRUE, latent = TRUE),
  control = list(adapt_delta = 0.95, max_treedepth = 20),
  seed = 1234
)

###* 
###* 1 is visitation frequency, 3 is morphospecies richness and 5 is functional
###* group richness
###* 
###* 
###* Null model 
saveRDS(bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt)
# needs kfold due to pareto
kfold_bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt <- kfold(bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt, K = 5, cores = 5)
saveRDS(kfold_bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt.rds")
print(kfold_bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt)


###* Visitation + morphospecies model
saveRDS(bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt)
# needs kfold due to pareto
kfold_bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt <- kfold(bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt, K = 5, cores = 5)
saveRDS(kfold_bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt.rds")
print(kfold_bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt)


###* Visitation only model
saveRDS(bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt)
# needs kfold due to pareto
kfold_bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt <- kfold(bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt, K = 5, cores = 5)
saveRDS(kfold_bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt.rds")
print(kfold_bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt)


###* Morphsopecies only model
saveRDS(bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt)
# needs kfold due to pareto
kfold_bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt <- kfold(bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt, K = 5, cores = 5)
saveRDS(kfold_bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt.rds")
print(kfold_bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt)


###* Functional group richness only model
saveRDS(bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt)
# needs kfold due to pareto
kfold_bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt <- kfold(bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt, K = 5, cores = 5)
saveRDS(kfold_bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt.rds")
print(kfold_bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt)


###* Morpshoepcies and functional group richness model
saveRDS(bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt)
# needs kfold due to pareto
kfold_bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt <- kfold(bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt, K = 5, cores = 5)
saveRDS(kfold_bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt.rds")
print(kfold_bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt)


###* Visitation, morphospecies and functional group richness model
saveRDS(bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt)
# needs kfold due to pareto
kfold_bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt <- kfold(bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt, K = 5, cores = 5)
saveRDS(kfold_bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt.rds")
print(kfold_bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt)


###* Visitation and functional group richness model
saveRDS(bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt)
# needs kfold due to pareto
kfold_bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt <- kfold(bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt, K = 5, cores = 5)
saveRDS(kfold_bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt.rds")
print(kfold_bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt)


###* Loading all kfolded opbject to not have to do it individually
kfold_bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt.rds")

###* Comparing existing kfolded object to see, which model is best
loo_compare(kfold_bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt,
            kfold_bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt,
            kfold_bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt,
            kfold_bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt,
            kfold_bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt,
            kfold_bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt,
            kfold_bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt,
            kfold_bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt)
###* Four models are indistinguishable. We will not use model 135 or 35, both 
###* include both correlated richness indices

loo_compare(kfold_bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt,
            kfold_bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt)
pp_check(bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt, type = "dens_overlay")
pp_check(bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt, type = "stat", stat = "mean")
pp_check(bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt, type = "stat", stat = "sd")


loo_compare(kfold_bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt,
            kfold_bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt)
pp_check(bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt, type = "dens_overlay")
pp_check(bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt, type = "stat", stat = "mean")
pp_check(bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt, type = "stat", stat = "sd")
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
###* POLLEN LIMITATION
###* 
###* 

# 1) Filter out the undefined PL case (your rule)
c.pl.final.table.5 <- c.pl.final.table.4 %>%
  filter(species != "Hypericum r" | elevation != 4000)

# 2) Build the predictors + SEs (use SE = SD/sqrt(n) just like you did for seedset)
eps <- 1e-8
c.pl.final.table.5 <- c.pl.final.table.5 %>%
  mutate(
    # the means you want to use as predictors
    vpm    = mean.visited.flowers.per.minute,
    morpho = mean.morpho,
    func   = mean.func,
    
    # SEs of those means: SD / sqrt(n_reps)
    vpm_se    = pmax(sd.visited.flowers.per.minute / sqrt(pmax(n_reps, 1L)), eps),
    morpho_se = pmax(sd.morpho                     / sqrt(pmax(n_reps, 1L)), eps),
    func_se   = pmax(sd.func                       / sqrt(pmax(n_reps, 1L)), eps),
    
    # --- means on the z-scale (already computed) ---
    x = z_flow_mean,   # visited flowers per minute (scaled) -- or use z_vis_mean if that's your x
    y = z_morpho_mean,
    z = z_func_mean,
    
    # --- SEs of the means on the same scale (SD / sqrt(n)) ---
    sx = pmax(z_flow_sd   / sqrt(pmax(n_reps,  1L)), eps),
    sy = pmax(z_morpho_sd  / sqrt(pmax(n_reps, 1L)), eps),
    sz = pmax(z_func_sd   / sqrt(pmax(n_reps,  1L)), eps)
  )

###* POLLEN LIMITATION 
View(c.pl.final.table.5)

bf_all <- bf(
  mean_PL_index | weights(PL_index_weight_12) ~
    me(x, sx) +
    me(y, sy) +
    me(z, sz) +
    #1 +
    (1 | species),
  phi ~ 1,
  zoi ~ 1,
  coi ~ 1
)


###* Create a weakly informative prior for the intercept:
###* For the ZINB model, the mean parameter is on the log scale.
###* We initialize the intercept prior around log(mean(seed set)), using the
###* overall mean of the response (plus a small constant to avoid log(0)).
w <- c.pl.final.table.5$PL_index_weight_12
pl <- c.pl.final.table.5$mean_PL_index

mu_bar <- weighted.mean(pl[pl > 0 & pl < 1], w[pl > 0 & pl < 1], na.rm = TRUE)
# fall back safely if needed
mu_bar <- ifelse(is.finite(mu_bar), mu_bar, weighted.mean(pl, w, na.rm = TRUE))
logit_mu0 <- qlogis(pmin(pmax(mu_bar, 1e-6), 1 - 1e-6))
logit_mu0
-0.1101463


pri_pl <- c(
  prior(normal(-0.1101463, 1.5), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(student_t(3, 0, 2.5), class = "sd", group = "species"),
  prior(normal(0, 1), class = "Intercept", dpar = "phi"),
  prior(normal(0, 1), class = "Intercept", dpar = "zoi"),
  prior(normal(0, 1), class = "Intercept", dpar = "coi")
)


bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2 <- brm(
  formula = bf_all,
  data    = c.pl.final.table.5,
  family  = zero_one_inflated_beta(),
  prior   = pri_pl,
  chains = 5, cores = 5, 
  iter = 2000, warmup = 1000,
  threads = threading(4),
  save_pars = save_pars(all = TRUE, latent = TRUE),
  control = list(adapt_delta = 0.95, max_treedepth = 20),
  seed = 1234
)

###* Visitation and morphospecies model
saveRDS(bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2, file = "../brms_models/bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2.rds")
bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2 <- readRDS("../brms_models/bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2.rds")
summary(bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2)
loo(bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2)
# needs kfold due to pareto
kfold_bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2 <- kfold(bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2, K = 5, cores = 5)
saveRDS(kfold_bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2, "../brms_models/kfold_bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2.rds")
kfold_bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2 <- readRDS("../brms_models/kfold_bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2.rds")
print(kfold_bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2)


###* Null model
saveRDS(bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2, file = "../brms_models/bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2.rds")
bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2 <- readRDS("../brms_models/bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2.rds")
summary(bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2)
loo(bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2)
# needs kfold due to pareto
kfold_bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2 <- kfold(bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2, K = 5, cores = 5)
saveRDS(kfold_bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2, "../brms_models/kfold_bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2.rds")
kfold_bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2 <- readRDS("../brms_models/kfold_bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2.rds")
print(kfold_bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2)


###* Visitation, morphsopecies and functional group model
saveRDS(bayesian_pl_mean_corrected_scaled_135_zoib_k_opt_2, file = "../brms_models/bayesian_pl_mean_corrected_scaled_135_zoib_k_opt_2.rds")
bayesian_pl_mean_corrected_scaled_135_zoib_k_opt_2 <- readRDS("../brms_models/bayesian_pl_mean_corrected_scaled_135_zoib_k_opt_2.rds")
summary(bayesian_pl_mean_corrected_scaled_135_zoib_k_opt_2)
loo(bayesian_pl_mean_corrected_scaled_135_zoib_k_opt_2)
# needs kfold due to pareto
pl_k5_scaled_135_2 <- kfold(bayesian_pl_mean_corrected_scaled_135_zoib_k_opt_2, K = 5, cores = 5)
saveRDS(pl_k5_scaled_135_2, "../brms_models/kfold_bayesian_pl_mean_corrected_scaled_135_zoib_k_opt_2.rds")
kfold_bayesian_pl_mean_corrected_scaled_135_zoib_k_opt_2 <- readRDS("../brms_models/kfold_bayesian_pl_mean_corrected_scaled_135_zoib_k_opt_2.rds")
print(kfold_bayesian_pl_mean_corrected_scaled_135_zoib_k_opt_2)


###* Visitation only model
saveRDS(bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2, file = "../brms_models/bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2.rds")
bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2 <- readRDS("../brms_models/bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2.rds")
summary(bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2)
loo(bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2)
# needs kfold due to pareto
pl_k5_scaled_1_2 <- kfold(bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2, K = 5, cores = 5)
saveRDS(pl_k5_scaled_1_2, "../brms_models/kfold_bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2.rds")
kfold_bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2 <- readRDS("../brms_models/kfold_bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2.rds")
print(kfold_bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2)


###* Morphospecies only model
saveRDS(bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2, file = "../brms_models/bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2.rds")
bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2 <- readRDS("../brms_models/bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2.rds")
summary(bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2)
loo(bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2)
# needs kfold due to pareto
pl_k5_scaled_3_2 <- kfold(bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2, K = 5, cores = 5)
saveRDS(pl_k5_scaled_3_2, "../brms_models/kfold_bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2.rds")
kfold_bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2 <- readRDS("../brms_models/kfold_bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2.rds")
print(kfold_bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2)


###* Functional group richness only model
saveRDS(bayesian_pl_mean_corrected_scaled_5_zoib_k_opt_2, file = "../brms_models/bayesian_pl_mean_corrected_scaled_5_zoib_k_opt_2.rds")
bayesian_pl_mean_corrected_scaled_5_zoib_k_opt_2 <- readRDS("../brms_models/bayesian_pl_mean_corrected_scaled_5_zoib_k_opt_2.rds")
summary(bayesian_pl_mean_corrected_scaled_5_zoib_k_opt_2)
loo(bayesian_pl_mean_corrected_scaled_5_zoib_k_opt_2)
# needs kfold due to pareto
pl_k5_scaled_5_2 <- kfold(bayesian_pl_mean_corrected_scaled_5_zoib_k_opt_2, K = 5, cores = 5)
saveRDS(pl_k5_scaled_5_2, "../brms_models/kfold_bayesian_pl_mean_corrected_scaled_5_zoib_k_opt_2.rds")
kfold_bayesian_pl_mean_corrected_scaled_5_zoib_k_opt_2 <- readRDS("../brms_models/kfold_bayesian_pl_mean_corrected_scaled_5_zoib_k_opt_2.rds")
print(kfold_bayesian_pl_mean_corrected_scaled_5_zoib_k_opt_2)


###* Morphospecies and functional group richness model
saveRDS(bayesian_pl_mean_corrected_scaled_35_zoib_k_opt_2, file = "../brms_models/bayesian_pl_mean_corrected_scaled_35_zoib_k_opt_2.rds")
bayesian_pl_mean_corrected_scaled_35_zoib_k_opt_2 <- readRDS("../brms_models/bayesian_pl_mean_corrected_scaled_35_zoib_k_opt_2.rds")
summary(bayesian_pl_mean_corrected_scaled_35_zoib_k_opt_2)
loo(bayesian_pl_mean_corrected_scaled_35_zoib_k_opt_2)
# needs kfold due to pareto
pl_k5_scaled_35_2 <- kfold(bayesian_pl_mean_corrected_scaled_35_zoib_k_opt_2, K = 5, cores = 5)
saveRDS(pl_k5_scaled_35_2, "../brms_models/kfold_bayesian_pl_mean_corrected_scaled_35_zoib_k_opt_2.rds")
kfold_bayesian_pl_mean_corrected_scaled_35_zoib_k_opt_2 <- readRDS("../brms_models/kfold_bayesian_pl_mean_corrected_scaled_35_zoib_k_opt_2.rds")
print(kfold_bayesian_pl_mean_corrected_scaled_35_zoib_k_opt_2)


###* Visitation and functional group richness model
saveRDS(bayesian_pl_mean_corrected_scaled_15_zoib_k_opt_2, file = "../brms_models/bayesian_pl_mean_corrected_scaled_15_zoib_k_opt_2.rds")
bayesian_pl_mean_corrected_scaled_15_zoib_k_opt_2 <- readRDS("../brms_models/bayesian_pl_mean_corrected_scaled_15_zoib_k_opt_2.rds")
summary(bayesian_pl_mean_corrected_scaled_15_zoib_k_opt_2)
loo(bayesian_pl_mean_corrected_scaled_15_zoib_k_opt_2)
# needs kfold due to pareto
pl_k5_scaled_15_2 <- kfold(bayesian_pl_mean_corrected_scaled_15_zoib_k_opt_2, K = 5, cores = 5)
saveRDS(pl_k5_scaled_15_2, "../brms_models/kfold_bayesian_pl_mean_corrected_scaled_15_zoib_k_opt_2.rds")
kfold_bayesian_pl_mean_corrected_scaled_15_zoib_k_opt_2 <- readRDS("../brms_models/kfold_bayesian_pl_mean_corrected_scaled_15_zoib_k_opt_2.rds")
print(kfold_bayesian_pl_mean_corrected_scaled_15_zoib_k_opt_2)


###* COomparison of individual kfolded objects
loo_compare(kfold_bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2,
            kfold_bayesian_pl_mean_corrected_scaled_135_zoib_k_opt_2,
            kfold_bayesian_pl_mean_corrected_scaled_15_zoib_k_opt_2,
            kfold_bayesian_pl_mean_corrected_scaled_35_zoib_k_opt_2,
            kfold_bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2,
            kfold_bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2,
            kfold_bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2,
            kfold_bayesian_pl_mean_corrected_scaled_5_zoib_k_opt_2)
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
###* 
###* 
###* 
###* AO INDEX
ao.final.table <- ao.final.table %>%
  rename(
    z_vis_mean = mean.visitors.scaled,
    z_vis_sd   = sd.visitors.scaled,
    z_flow_mean = mean.visited.flowers.scaled,
    z_flow_sd   = sd.visited.flowers.scaled,
    z_morpho_mean = mean.morpho.scaled,
    z_morpho_sd   = sd.morpho.scaled,
    z_func_mean   = mean.func.scaled,
    z_func_sd     = sd.func.scaled
  )

View(ao.final.table)

eps <- 1e-8  # tiny floor to avoid zero SDs

ao.final.table <- ao.final.table %>%
  mutate(
    # --- means on the z-scale (already computed) ---
    x = z_flow_mean,   # visited flowers per minute (scaled) -- or use z_vis_mean if that's your x
    y = z_morpho_mean,
    z = z_func_mean,
    
    # --- SEs of the means on the same scale (SD / sqrt(n)) ---
    sx = pmax(z_flow_sd   / sqrt(pmax(n_reps,  1L)), eps),
    sy = pmax(z_morpho_sd  / sqrt(pmax(n_reps, 1L)), eps),
    sz = pmax(z_func_sd   / sqrt(pmax(n_reps,  1L)), eps)
  )

form_all <- bf(
  mean_ao_index | weights(ao_index_weight_12) ~
    me(x, sx) +
    me(y, sy) +
    me(z, sz) +
    #1 +
    (1|species)
)

aow <- ao.final.table$ao_index_weight_12
ao <- ao.final.table$mean_ao_index

mu_bar_ao <- weighted.mean(ao[ao > 0 & ao < 1], aow[ao > 0 & ao < 1], na.rm = TRUE)
# fall back safely if needed
mu_bar_ao <- ifelse(is.finite(mu_bar_ao), mu_bar_ao, weighted.mean(ao, aow, na.rm = TRUE))
logit_ao_mu0 <- qlogis(pmin(pmax(mu_bar_ao, 1e-6), 1 - 1e-6))
logit_ao_mu0
-0.594663

# Suppose the names come back as: mexsx, mex2sx2, meys y..., mezsz, mez2sz2 etc.
pri_all <- c(
  # Fixed effects (logit link for mu):
  prior(normal(0, 1),    class = "b"),                       # all slopes by default
  prior(normal(-0.594663, 1.5), class = "Intercept"),              # mu-intercept on logit scale
  prior(student_t(3, 0, 2.5), class = "sd", group = "species")                    
)

View(ao.final.table)

bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt <- brm(
  form_all, 
  data = ao.final.table, 
  family = zero_inflated_beta(),
  prior = pri_all, 
  chains = 5, cores = 5, 
  iter = 2000, warmup = 1000,
  threads = threading(4),
  save_pars = save_pars(all = TRUE, latent = TRUE),
  control = list(adapt_delta = 0.95, max_treedepth = 20),
  seed = 1234
)

###* Visitation and morphospecies model
saveRDS(bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt, file = "../brms_models/bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt.rds")
bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt <- readRDS("../brms_models/bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt.rds")
summary(bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt)
loo(bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt)
# needs kfold due to pareto
kfold_bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt <- kfold(bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt, K = 5, cores = 5)
saveRDS(ao_k5_scaled_13_2, "../brms_models/kfold_bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt.rds")
print(kfold_bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt)


###* Null model
saveRDS(bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt, file = "../brms_models/bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt.rds")
bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt <- readRDS("../brms_models/bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt.rds")
summary(bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt)
loo(bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt)
# needs kfold due to pareto
kfold_bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt <- kfold(bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt, K = 5, cores = 5)
saveRDS(kfold_bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt, "../brms_models/kfold_bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt.rds")
print(kfold_bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt)


###* Visitation only model
saveRDS(bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt, file = "../brms_models/bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt.rds")
bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt <- readRDS("../brms_models/bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt.rds")
summary(bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt)
loo(bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt)
# needs kfold due to pareto
kfold_bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt <- kfold(bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt, K = 5, cores = 5)
saveRDS(kfold_bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt, "../brms_models/kfold_bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt.rds")
print(kfold_bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt)


###* Morphsopecies only model
saveRDS(bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt, file = "../brms_models/bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt.rds")
bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt <- readRDS("../brms_models/bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt.rds")
summary(bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt)
loo(bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt)
# needs kfold due to pareto
kfold_bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt <- kfold(bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt, K = 5, cores = 5)
saveRDS(kfold_bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt, "../brms_models/kfold_bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt.rds")
print(kfold_bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt)


###* Functional group only model
saveRDS(bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt, file = "../brms_models/bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt.rds")
bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt <- readRDS("../brms_models/bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt.rds")
summary(bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt)
loo(bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt)
# needs kfold due to pareto
kfold_bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt <- kfold(bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt, K = 5, cores = 5)
saveRDS(kfold_bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt, "../brms_models/kfold_bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt.rds")
print(kfold_bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt)


###* Visitation and functional group richness model
saveRDS(bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt, file = "../brms_models/bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt.rds")
bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt <- readRDS("../brms_models/bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt.rds")
summary(bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt)
loo(bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt)
# needs kfold due to pareto
kfold_bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt <- kfold(bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt, K = 5, cores = 5)
saveRDS(kfold_bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt, "../brms_models/kfold_bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt.rds")
print(kfold_bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt)


###* Morphospecies richness and functional group richness model
saveRDS(bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt, file = "../brms_models/bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt.rds")
bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt <- readRDS("../brms_models/bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt.rds")
summary(bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt)
loo(bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt)
# needs kfold due to pareto
kfold_bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt <- kfold(bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt, K = 5, cores = 5)
saveRDS(kfold_bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt, "../brms_models/kfold_bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt.rds")
print(kfold_bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt)


# Visitation, morphospecies and functional group richness model
saveRDS(bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt, file = "../brms_models/bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt.rds")
bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt <- readRDS("../brms_models/bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt.rds")
summary(bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt)
loo(bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt)
# needs kfold due to pareto
kfold_bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt <- kfold(bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt, K = 5, cores = 5)
saveRDS(kfold_bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt, "../brms_models/kfold_bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt.rds")
print(kfold_bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt)


kfold_bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt.rds")


loo_compare(kfold_bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt,
            kfold_bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt,
            kfold_bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt,
            kfold_bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt,
            kfold_bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt,
            kfold_bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt,
            kfold_bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt,
            kfold_bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt)

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
###*
###* GO INDEX
go.final.table <- go.final.table %>%
  rename(
    z_vis_mean = mean.visitors.scaled,
    z_vis_sd   = sd.visitors.scaled,
    z_flow_mean = mean.visited.flowers.scaled,
    z_flow_sd   = sd.visited.flowers.scaled,
    z_morpho_mean = mean.morpho.scaled,
    z_morpho_sd   = sd.morpho.scaled,
    z_func_mean   = mean.func.scaled,
    z_func_sd     = sd.func.scaled
  )

View(go.final.table)

eps <- 1e-8  # tiny floor to avoid zero SDs

go.final.table <- go.final.table %>%
  mutate(
    # --- means on the z-scale (already computed) ---
    x = z_flow_mean,   # visited flowers per minute (scaled) -- or use z_vis_mean if that's your x
    y = z_morpho_mean,
    z = z_func_mean,
    
    # --- SEs of the means on the same scale (SD / sqrt(n)) ---
    sx = pmax(z_flow_sd   / sqrt(pmax(n_reps,  1L)), eps),
    sy = pmax(z_morpho_sd  / sqrt(pmax(n_reps, 1L)), eps),
    sz = pmax(z_func_sd   / sqrt(pmax(n_reps,  1L)), eps)
  )


###* Create a weakly informative prior for the intercept:
###* For the ZINB model, the mean parameter is on the log scale.
###* We initialize the intercept prior around log(mean(seed set)), using the
###* overall mean of the response (plus a small constant to avoid log(0)).
gow <- go.final.table$go_index_weight_12
go <- go.final.table$mean_go_index

mu_bar_go <- weighted.mean(go[go > 0 & go < 1], gow[go > 0 & go < 1], na.rm = TRUE)
mu_bar_go <- ifelse(is.finite(mu_bar_go), mu_bar_go, weighted.mean(go, gow, na.rm = TRUE))
logit_go_mu0 <- qlogis(pmin(pmax(mu_bar_go, 1e-6), 1 - 1e-6))
logit_go_mu0

form_all <- bf(
  mean_go_index | weights(go_index_weight_12) ~
    #me(x, sx) +
    me(y, sy) +
    me(z, sz) + 
    1 +
    (1|species),
  phi ~ 1,
  zoi ~1,
  coi ~1
)

# Suppose the names come back as: mexsx, mex2sx2, meys y..., mezsz, mez2sz2 etc.
pri_all <- c(
  # Fixed effects (logit link for mu):
  prior(normal(0, 1),    class = "b"),                      # all slopes by default
  prior(normal(-0.5758654, 1.5),    class = "Intercept"),   # mu-intercept on logit scale
  prior(constant(-6, 2), class="Intercept", dpar="zoi"),    # median ≈ 0.0025, but flexible
  prior(constant(-2.6, 1), class="Intercept", dpar="coi"),  # median ≈ 0.07, but flexible
  prior(student_t(3, 0, 2.5), class = "sd", group = "species")
)


View(go.final.table)

bayesian_go_mean_corrected_scaled_35_2_zoib_k_opt <- brm(
  form_all, 
  data = go.final.table, 
  family = zero_one_inflated_beta(),
  prior = pri_all, 
  chains = 5, cores = 5, 
  iter = 2000, warmup = 1000,
  threads = threading(4),
  save_pars = save_pars(all = TRUE, latent = TRUE),
  control = list(adapt_delta = 0.95, max_treedepth = 20),
  seed = 123
)

###* Visitation only model
saveRDS(bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt.rds")
bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt)
loo(bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt)

kfold_bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt <- kfold(bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt, K = 5, cores = 5)
saveRDS(kfold_bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt)


###* Visitation and morphospecies richness model
saveRDS(bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt.rds")
bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt)
loo(bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt)

kfold_bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt <- kfold(bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt, K = 5, cores = 5)
saveRDS(kfold_bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt)


###* Null model
saveRDS(bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt.rds")
bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt)
loo(bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt)

kfold_bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt <- kfold(bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt, K = 5, cores = 5)
saveRDS(kfold_bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt)

###* RUNRUNRUNRUNR
###* Visitation and functional group richness model
saveRDS(bayesian_go_mean_corrected_scaled_15_2_zoib_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_15_2_zoib_k_opt.rds")
bayesian_go_mean_corrected_scaled_15_2_zoib_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_15_2_zoib_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_15_2_zoib_k_opt)
loo(bayesian_go_mean_corrected_scaled_15_2_zoib_k_opt)

kfold_bayesian_go_mean_corrected_scaled_15_2_zoib_k_opt <- kfold(bayesian_go_mean_corrected_scaled_15_2_zoib_k_opt, K = 5, cores = 5)
saveRDS(kfold_bayesian_go_mean_corrected_scaled_15_2_zoib_k_opt, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_15_2_zoib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_15_2_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_15_2_zoib_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_15_2_zoib_k_opt)


###* Morphospecies richness model
saveRDS(bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt.rds")
bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt)
loo(bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt)

kfold_bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt <- kfold(bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt, K = 5, cores = 5)
saveRDS(kfold_bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt)


###* Morphospecies richness and functional gropup richness model
saveRDS(bayesian_go_mean_corrected_scaled_35_2_zoib_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_35_2_zoib_k_opt.rds")
bayesian_go_mean_corrected_scaled_35_2_zoib_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_35_2_zoib_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_35_2_zoib_k_opt)
loo(bayesian_go_mean_corrected_scaled_35_2_zoib_k_opt)

kfold_bayesian_go_mean_corrected_scaled_35_2_zoib_k_opt <- kfold(bayesian_go_mean_corrected_scaled_35_2_zoib_k_opt, K = 5, cores = 5)
saveRDS(kfold_bayesian_go_mean_corrected_scaled_35_2_zoib_k_opt, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_35_2_zoib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_35_2_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_35_2_zoib_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_35_2_zoib_k_opt)


###* Visitation, morphospecies and functional group richness
saveRDS(bayesian_go_mean_corrected_scaled_135_2_zoib_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_135_2_zoib_k_opt.rds")
bayesian_go_mean_corrected_scaled_135_2_zoib_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_135_2_zoib_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_135_2_zoib_k_opt)
loo(bayesian_go_mean_corrected_scaled_135_2_zoib_k_opt)

kfold_bayesian_go_mean_corrected_scaled_135_2_zoib_k_opt <- kfold(bayesian_go_mean_corrected_scaled_135_2_zoib_k_opt, K = 5, cores = 5)
saveRDS(kfold_bayesian_go_mean_corrected_scaled_135_2_zoib_k_opt, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_135_2_zoib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_135_2_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_135_2_zoib_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_135_2_zoib_k_opt)


###* Visitation, morphospecies and functional group richness
saveRDS(bayesian_go_mean_corrected_scaled_5_2_zoib_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_5_2_zoib_k_opt.rds")
bayesian_go_mean_corrected_scaled_5_2_zoib_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_5_2_zoib_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_5_2_zoib_k_opt)
loo(bayesian_go_mean_corrected_scaled_5_2_zoib_k_opt)

kfold_bayesian_go_mean_corrected_scaled_5_2_zoib_k_opt <- kfold(bayesian_go_mean_corrected_scaled_5_2_zoib_k_opt, K = 5, cores = 5)
saveRDS(kfold_bayesian_go_mean_corrected_scaled_5_2_zoib_k_opt, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_5_2_zoib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_5_2_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_5_2_zoib_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_5_2_zoib_k_opt)





