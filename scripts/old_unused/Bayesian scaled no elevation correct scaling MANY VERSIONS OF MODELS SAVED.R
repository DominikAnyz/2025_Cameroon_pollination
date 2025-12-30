###* Bayesian scaled June
###* 
###* 
###* 
###* 
###* This script serves as running and interpreting the results from the best 
###* models created with the data which we have.
source("scripts/Bayesian setup June.R")

library(brms)
options(brms.backend = "cmdstanr")

# View(c.pl.final.table.4)
# 
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

eps <- 1e-8  # tiny floor to avoid zero SDs

c.pl.final.table.4 <- c.pl.final.table.4 %>%
  mutate(
    # --- means on the z-scale (already computed) ---
    x = z_flow_mean,   # visited flowers per minute (scaled) -- or use z_vis_mean if that's your x
    y = z_morpho_mean,
    z = z_func_mean,
    # --- SEs of the means on the same scale (SD / sqrt(n)) ---
    sx = pmax(z_flow_sd   / sqrt(pmax(n_reps,  1L)), eps),
    sy = pmax(z_morpho_sd  / sqrt(pmax(n_reps, 1L)), eps),
    sz = pmax(z_func_sd   / sqrt(pmax(n_reps,  1L)), eps),
    
    # --- quadratic terms and their delta-method SEs ---
    x2  = x^2,  sx2 = pmax(2 * abs(x) * sx, eps),
    y2  = y^2,  sy2 = pmax(2 * abs(y) * sy, eps),
    z2  = z^2,  sz2 = pmax(2 * abs(z) * sz, eps)
  )

cor.test(c.pl.final.table.4$x, c.pl.final.table.4$y)
cor.test(c.pl.final.table.4$y, c.pl.final.table.4$z)

View(c.pl.final.table.4)

form_all <- bf(
  mean_seedset_round | weights(seedset_weight_12) ~
    me(x, sx) +
    #me(x2, sx2) +
    me(y, sy) +
    #me(y2, sy2) +
    #me(z, sz) + 
    #1 +
    (1|species),
    zi ~ 1
)


pri_all <- c(
  prior(normal(0, 0.5), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(student_t(3, 0, 2.5), class = "sd", group = "species")
  #prior(normal(0, 0.5), class = "b", coef = "mex2sx2"),
  #prior(normal(0, 0.5), class = "b", coef = "mey2sy2")
)

bayesian_seedset_mean_corrected_scaled_1234_nb_k_opt <- brm(
  form_all, 
  data = c.pl.final.table.4,
  family = negbinomial(), 
  prior = pri_all,
  chains = 5, cores = 5, 
  iter = 2000, warmup = 1000,
  threads = threading(4),
  save_pars = save_pars(all = TRUE, latent = TRUE),
  control = list(adapt_delta = 0.95, max_treedepth = 20),
  seed = 1234
)

bayesian_seedset_mean_corrected_scaled_13_zinb_k_opt <- brm(
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

saveRDS(bayesian_seedset_mean_corrected_scaled_13_nb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_13_nb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_13_nb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_13_nb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_13_nb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_13_nb_k_opt)

seedset_scaled_k5_13 <- kfold(bayesian_seedset_mean_corrected_scaled_13_nb_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
seedset_scaled_k5_13                   # prints the elpd_kfold etc.
print(seedset_scaled_k5_13)
summary(seedset_scaled_k5_13)          # shows estimates and pointwise columns
seedset_scaled_k5_13$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(seedset_scaled_k5_13$pointwise)   # fold-wise contributions per observation
saveRDS(seedset_scaled_k5_13, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_13_nb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_13_nb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_13_nb_k_opt.rds")
print(seedset_scaled_k5_13)


saveRDS(bayesian_seedset_mean_corrected_scaled_null_nb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_null_nb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_null_nb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_null_nb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_null_nb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_null_nb_k_opt)

seedset_scaled_k5_null <- kfold(bayesian_seedset_mean_corrected_scaled_null_nb_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
seedset_scaled_k5_null                   # prints the elpd_kfold etc.
print(seedset_scaled_k5_null)
summary(seedset_scaled_k5_null)          # shows estimates and pointwise columns
seedset_scaled_k5_null$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(seedset_scaled_k5_null$pointwise)   # fold-wise contributions per observation
saveRDS(seedset_scaled_k5_null, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_null_nb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_null_nb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_null_nb_k_opt.rds")
print(seedset_scaled_k5_null)


saveRDS(bayesian_seedset_mean_corrected_scaled_1_nb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_1_nb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_1_nb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_1_nb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_1_nb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_1_nb_k_opt)

seedset_scaled_k5_1 <- kfold(bayesian_seedset_mean_corrected_scaled_1_nb_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
seedset_scaled_k5_1                   # prints the elpd_kfold etc.
print(seedset_scaled_k5_1)
summary(seedset_scaled_k5_1)          # shows estimates and pointwise columns
seedset_scaled_k5_1$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(seedset_scaled_k5_1$pointwise)   # fold-wise contributions per observation
saveRDS(seedset_scaled_k5_1, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_1_nb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_1_nb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_1_nb_k_opt.rds")
print(seedset_scaled_k5_1)


saveRDS(bayesian_seedset_mean_corrected_scaled_3_nb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_3_nb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_3_nb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_3_nb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_3_nb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_3_nb_k_opt)

seedset_scaled_k5_3 <- kfold(bayesian_seedset_mean_corrected_scaled_3_nb_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
seedset_scaled_k5_3                   # prints the elpd_kfold etc.
print(seedset_scaled_k5_3)
summary(seedset_scaled_k5_3)          # shows estimates and pointwise columns
seedset_scaled_k5_3$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(seedset_scaled_k5_3$pointwise)   # fold-wise contributions per observation
saveRDS(seedset_scaled_k5_3, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_3_nb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_3_nb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_3_nb_k_opt.rds")
print(seedset_scaled_k5_3)


saveRDS(bayesian_seedset_mean_corrected_scaled_12_nb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_12_nb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_12_nb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_12_nb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_12_nb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_12_nb_k_opt)

seedset_scaled_k5_12 <- kfold(bayesian_seedset_mean_corrected_scaled_12_nb_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
seedset_scaled_k5_12                   # prints the elpd_kfold etc.
print(seedset_scaled_k5_12)
summary(seedset_scaled_k5_12)          # shows estimates and pointwise columns
seedset_scaled_k5_12$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(seedset_scaled_k5_12$pointwise)   # fold-wise contributions per observation
saveRDS(seedset_scaled_k5_12, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_12_nb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_12_nb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_12_nb_k_opt.rds")
print(seedset_scaled_k5_12)


saveRDS(bayesian_seedset_mean_corrected_scaled_34_nb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_34_nb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_34_nb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_34_nb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_34_nb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_34_nb_k_opt)

seedset_scaled_k5_34 <- kfold(bayesian_seedset_mean_corrected_scaled_34_nb_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
seedset_scaled_k5_34                   # prints the elpd_kfold etc.
print(seedset_scaled_k5_34)
summary(seedset_scaled_k5_34)          # shows estimates and pointwise columns
seedset_scaled_k5_34$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(seedset_scaled_k5_34$pointwise)   # fold-wise contributions per observation
saveRDS(seedset_scaled_k5_34, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_34_nb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_34_nb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_34_nb_k_opt.rds")
print(seedset_scaled_k5_34)


saveRDS(bayesian_seedset_mean_corrected_scaled_1234_nb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_1234_nb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_1234_nb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_1234_nb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_1234_nb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_1234_nb_k_opt)

seedset_scaled_k5_1234 <- kfold(bayesian_seedset_mean_corrected_scaled_1234_nb_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
seedset_scaled_k5_1234                   # prints the elpd_kfold etc.
print(seedset_scaled_k5_1234)
summary(seedset_scaled_k5_1234)          # shows estimates and pointwise columns
seedset_scaled_k5_1234$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(seedset_scaled_k5_1234$pointwise)   # fold-wise contributions per observation
saveRDS(seedset_scaled_k5_1234, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_1234_nb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_1234_nb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_1234_nb_k_opt.rds")
print(seedset_scaled_k5_1234)


loo_compare(seedset_scaled_k5_13, seedset_scaled_k5_null, seedset_scaled_k5_1, seedset_scaled_k5_3, seedset_scaled_k5_12, seedset_scaled_k5_34, seedset_scaled_k5_1234)





saveRDS(bayesian_seedset_mean_corrected_scaled_13_zip_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_13_zip_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_13_zip_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_13_zip_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_13_zip_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_13_zip_k_opt)

seedset_scaled_k5_13_zip <- kfold(bayesian_seedset_mean_corrected_scaled_13_zip_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
seedset_scaled_k5_13_zip                   # prints the elpd_kfold etc.
print(seedset_scaled_k5_13_zip)
summary(seedset_scaled_k5_13_zip)          # shows estimates and pointwise columns
seedset_scaled_k5_13_zip$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(seedset_scaled_k5_13_zip$pointwise)   # fold-wise contributions per observation
saveRDS(seedset_scaled_k5_13_zip, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_13_zip_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_13_zip_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_13_zip_k_opt.rds")
print(seedset_scaled_k5_13_zip)








form_all <- bf(
  mean_seedset_round | weights(seedset_weight_12) ~
    #me(x, sx) +
    #me(x2, sx2) +
    #me(y, sy) +
    #me(y2, sy2) +
    #me(z, sz) + 
    1 +
    (1|species),
  zi ~ 1
)


pri_all <- c(
  prior(normal(0, 0.5), class = "Intercept"),
  #prior(normal(0, 1), class = "b"),
  prior(student_t(3, 0, 2.5), class = "sd", group = "species")
  #prior(normal(0, 0.5), class = "b", coef = "mex2sx2"),
  #prior(normal(0, 0.5), class = "b", coef = "mey2sy2")
)


bayesian_seedset_mean_corrected_scaled_null_zinb_k_opt <- brm(
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


saveRDS(bayesian_seedset_mean_corrected_scaled_13_zinb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_13_zinb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_13_zinb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_13_zinb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_13_zinb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_13_zinb_k_opt)

seedset_scaled_k5_13_zinb <- kfold(bayesian_seedset_mean_corrected_scaled_13_zinb_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
seedset_scaled_k5_13_zinb                   # prints the elpd_kfold etc.
print(seedset_scaled_k5_13_zinb)
summary(seedset_scaled_k5_13_zinb)          # shows estimates and pointwise columns
seedset_scaled_k5_13_zinb$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(seedset_scaled_k5_13_zinb$pointwise)   # fold-wise contributions per observation
saveRDS(seedset_scaled_k5_13_zinb, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_13_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_13_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_13_zinb_k_opt.rds")
print(seedset_scaled_k5_13_zinb)


saveRDS(bayesian_seedset_mean_corrected_scaled_1_zinb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_1_zinb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_1_zinb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_1_zinb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_1_zinb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_1_zinb_k_opt)

seedset_scaled_k5_1_zinb <- kfold(bayesian_seedset_mean_corrected_scaled_1_zinb_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
seedset_scaled_k5_1_zinb                   # prints the elpd_kfold etc.
print(seedset_scaled_k5_1_zinb)
summary(seedset_scaled_k5_1_zinb)          # shows estimates and pointwise columns
seedset_scaled_k5_1_zinb$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(seedset_scaled_k5_1_zinb$pointwise)   # fold-wise contributions per observation
saveRDS(seedset_scaled_k5_1_zinb, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_1_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_1_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_1_zinb_k_opt.rds")
print(seedset_scaled_k5_1_zinb)


saveRDS(bayesian_seedset_mean_corrected_scaled_3_zinb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_3_zinb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_3_zinb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_3_zinb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_3_zinb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_3_zinb_k_opt)

seedset_scaled_k5_3_zinb <- kfold(bayesian_seedset_mean_corrected_scaled_3_zinb_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
seedset_scaled_k5_3_zinb                   # prints the elpd_kfold etc.
print(seedset_scaled_k5_3_zinb)
summary(seedset_scaled_k5_3_zinb)          # shows estimates and pointwise columns
seedset_scaled_k5_3_zinb$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(seedset_scaled_k5_3_zinb$pointwise)   # fold-wise contributions per observation
saveRDS(seedset_scaled_k5_3_zinb, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_3_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_3_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_3_zinb_k_opt.rds")
print(seedset_scaled_k5_3_zinb)


saveRDS(bayesian_seedset_mean_corrected_scaled_null_zinb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_null_zinb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_null_zinb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_null_zinb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_null_zinb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_null_zinb_k_opt)

seedset_scaled_k5_null_zinb <- kfold(bayesian_seedset_mean_corrected_scaled_null_zinb_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
seedset_scaled_k5_null_zinb                   # prints the elpd_kfold etc.
print(seedset_scaled_k5_null_zinb)
summary(seedset_scaled_k5_null_zinb)          # shows estimates and pointwise columns
seedset_scaled_k5_null_zinb$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(seedset_scaled_k5_null_zinb$pointwise)   # fold-wise contributions per observation
saveRDS(seedset_scaled_k5_null_zinb, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_null_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_null_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_null_zinb_k_opt.rds")
print(seedset_scaled_k5_null_zinb)

loo_compare(seedset_scaled_k5_13, seedset_scaled_k5_13_zip, seedset_scaled_k5_13_zinb)

loo_compare(seedset_scaled_k5_13_zinb, seedset_scaled_k5_1_zinb, seedset_scaled_k5_3_zinb, seedset_scaled_k5_null_zinb)



















library(dplyr)
library(ggplot2)
m <- bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt

# 1) Set "typical" values for the nuisance pieces
sx_med <- median(c.pl.final.table.4$sx, na.rm = TRUE)   # ME SD for x
sy_med <- median(c.pl.final.table.4$sy, na.rm = TRUE)   # ME SD for y

# 2) Build a grid over the observed x range (z-scaled visitation)
x_seq <- seq(min(c.pl.final.table.4$x, na.rm = TRUE),
             max(c.pl.final.table.4$x, na.rm = TRUE),
             length.out = 200)

newdat <- tibble(
  x  = x_seq,
  sx = sx_med,        # required by me(x, sx)
  y  = 0,             # hold morpho at its z-mean
  sy = sy_med,        # required by me(y, sy)
  # species is ignored when we drop REs below
  species = NA
)

# 3) Get posterior fitted means (includes zero-inflation) without random effects
pred <- fitted(
  m, newdata = newdat,
  re_formula = NA,          # population-level (no species RE)
  summary    = TRUE,
  ndraws     = 2000
) %>% as.data.frame()

plot_df <- bind_cols(newdat, pred)  # columns: Estimate, Q2.5, Q97.5, etc.

# 4) Plot: model curve + 95% CI ribbon + observed points
ggplot(plot_df, aes(x = x, y = Estimate)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.2) +
  geom_line(size = 1) +
  geom_point(data = c.pl.final.table.4,
             aes(x = x, y = mean_seedset_round, size = seedset_weight_12),
             inherit.aes = FALSE, alpha = 0.7) +
  scale_size_continuous(name = "Weight") +
  labs(x = "Visitation (z-scaled)", y = "Expected seed set",
       title = "Seed set vs visitation (ZINB, population-level)") +
  theme_minimal(base_size = 12)











ce <- conditional_effects(
  m, effects = "x",
  re_formula = NA,                # population-level
  method = "posterior_epred",     # expected counts incl. zero inflation
  ndraws = 2000
)
plot(ce, points = TRUE)















###* CORRECT MODEL
###* CORRECT MODEL
###* CORRECT MODEL
###* CORRECT MODEL
###* CORRECT MODEL
###* THIS IS THE CORRECT MODEL WHICH WE ARE USING:
###* 

View(c.pl.final.table.4)

form_all <- bf(
  mean_seedset_round | weights(seedset_weight_12) ~
    me(x, sx) +
    #me(x2, sx2) +
    me(y, sy) +
    #me(y2, sy2) +
    #me(z, sz) + 
    #1 +
    (1|species),
  zi ~ 1
)


mu0 <- log(mean(c.pl.final.table.4$mean_seedset_round, na.rm = TRUE) + 1e-8)

mu0

pri_all <- c(
  prior(normal(3.688879, 1.5), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(student_t(3, 0, 2.5), class = "sd", group = "species")
)

#prior(normal(0, 0.5), class = "b", coef = "mex2sx2"),
  #prior(normal(0, 0.5), class = "b", coef = "mey2sy2")

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


saveRDS(bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt)

seedset_scaled_k5_null2_zinb <- kfold(bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
seedset_scaled_k5_null2_zinb                   # prints the elpd_kfold etc.
print(seedset_scaled_k5_null2_zinb)
summary(seedset_scaled_k5_null2_zinb)          # shows estimates and pointwise columns
seedset_scaled_k5_null2_zinb$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(seedset_scaled_k5_null2_zinb$pointwise)   # fold-wise contributions per observation
saveRDS(seedset_scaled_k5_null2_zinb, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt.rds")
print(seedset_scaled_k5_null2_zinb)


saveRDS(bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt)

seedset_scaled_k5_13_2_zinb <- kfold(bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
seedset_scaled_k5_13_2_zinb                   # prints the elpd_kfold etc.
print(seedset_scaled_k5_13_2_zinb)
summary(seedset_scaled_k5_13_2_zinb)          # shows estimates and pointwise columns
seedset_scaled_k5_13_2_zinb$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(seedset_scaled_k5_13_2_zinb$pointwise)   # fold-wise contributions per observation
saveRDS(seedset_scaled_k5_13_2_zinb, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt.rds")
print(seedset_scaled_k5_13_2_zinb)

pp_check(bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt, type = "dens_overlay")
pp_check(bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt, type = "stat", stat = "mean")
pp_check(bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt, type = "stat", stat = "sd")




saveRDS(bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt)

seedset_scaled_k5_1_2_zinb <- kfold(bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
seedset_scaled_k5_1_2_zinb                   # prints the elpd_kfold etc.
print(seedset_scaled_k5_1_2_zinb)
summary(seedset_scaled_k5_1_2_zinb)          # shows estimates and pointwise columns
seedset_scaled_k5_1_2_zinb$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(seedset_scaled_k5_1_2_zinb$pointwise)   # fold-wise contributions per observation
saveRDS(seedset_scaled_k5_1_2_zinb, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt.rds")
print(seedset_scaled_k5_1_2_zinb)

loo_compare(seedset_scaled_k5_13_2_zinb, seedset_scaled_k5_null2_zinb, seedset_scaled_k5_1_2_zinb)


loo_compare(kfold_bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt,
            kfold_bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt,
            kfold_bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt)




saveRDS(bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt)

seedset_scaled_k5_3_2_zinb <- kfold(bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
seedset_scaled_k5_3_2_zinb                   # prints the elpd_kfold etc.
print(seedset_scaled_k5_3_2_zinb)
summary(seedset_scaled_k5_3_2_zinb)          # shows estimates and pointwise columns
seedset_scaled_k5_3_2_zinb$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(seedset_scaled_k5_3_2_zinb$pointwise)   # fold-wise contributions per observation
saveRDS(seedset_scaled_k5_3_2_zinb, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt.rds")
print(seedset_scaled_k5_3_2_zinb)


saveRDS(bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt)

seedset_scaled_k5_5_2_zinb <- kfold(bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
seedset_scaled_k5_5_2_zinb                   # prints the elpd_kfold etc.
print(seedset_scaled_k5_5_2_zinb)
summary(seedset_scaled_k5_5_2_zinb)          # shows estimates and pointwise columns
seedset_scaled_k5_5_2_zinb$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(seedset_scaled_k5_5_2_zinb$pointwise)   # fold-wise contributions per observation
saveRDS(seedset_scaled_k5_5_2_zinb, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt.rds")
print(seedset_scaled_k5_5_2_zinb)


saveRDS(bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt)

seedset_scaled_k5_35_2_zinb <- kfold(bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
seedset_scaled_k5_35_2_zinb                   # prints the elpd_kfold etc.
print(seedset_scaled_k5_35_2_zinb)
summary(seedset_scaled_k5_35_2_zinb)          # shows estimates and pointwise columns
seedset_scaled_k5_35_2_zinb$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(seedset_scaled_k5_35_2_zinb$pointwise)   # fold-wise contributions per observation
saveRDS(seedset_scaled_k5_35_2_zinb, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt.rds")
print(seedset_scaled_k5_35_2_zinb)



saveRDS(bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt)

seedset_scaled_k5_135_2_zinb <- kfold(bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
seedset_scaled_k5_135_2_zinb                   # prints the elpd_kfold etc.
print(seedset_scaled_k5_135_2_zinb)
summary(seedset_scaled_k5_135_2_zinb)          # shows estimates and pointwise columns
seedset_scaled_k5_135_2_zinb$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(seedset_scaled_k5_135_2_zinb$pointwise)   # fold-wise contributions per observation
saveRDS(seedset_scaled_k5_135_2_zinb, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt.rds")
print(seedset_scaled_k5_135_2_zinb)


saveRDS(bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt)

seedset_scaled_k5_15_2_zinb <- kfold(bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
seedset_scaled_k5_15_2_zinb                   # prints the elpd_kfold etc.
print(seedset_scaled_k5_15_2_zinb)
summary(seedset_scaled_k5_15_2_zinb)          # shows estimates and pointwise columns
seedset_scaled_k5_15_2_zinb$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(seedset_scaled_k5_15_2_zinb$pointwise)   # fold-wise contributions per observation
saveRDS(seedset_scaled_k5_15_2_zinb, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt.rds")
print(seedset_scaled_k5_15_2_zinb)



loo_compare(kfold_bayesian_seedset_mean_corrected_scaled_3_2_zinb_k_opt,
            kfold_bayesian_seedset_mean_corrected_scaled_5_2_zinb_k_opt,
            kfold_bayesian_seedset_mean_corrected_scaled_35_2_zinb_k_opt,
            kfold_bayesian_seedset_mean_corrected_scaled_135_2_zinb_k_opt,
            kfold_bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt,
            kfold_bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt,
            kfold_bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt,
            kfold_bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt)


loo_compare(kfold_bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt,
            kfold_bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt)












library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)

View(flowering.visited)

# your existing code
m   <- bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt
dat <- c.pl.final.table.4

vis_mean <- mean(flowering.visited$visitors.per.minute, na.rm = TRUE)
vis_sd   <- sd(flowering.visited$visitors.per.minute,   na.rm = TRUE)

sx_med <- median(dat$sx, na.rm = TRUE)
sz_med <- median(dat$sz, na.rm = TRUE)

xr    <- range(dat$x, na.rm = TRUE)
x_seq <- seq(xr[1], xr[2], length.out = 300)

newdat_vis <- tibble(
  x  = x_seq,    # z-scored visitation (predictor in model)
  sx = sx_med,
  z  = 0,        # functional richness at its mean (z-scale)
  sz = sz_med,
  species = NA
)

pred_vis <- fitted(
  m,
  newdata    = newdat_vis,
  re_formula = NA,
  summary    = TRUE
) |> as.data.frame()

plot_df_vis <- bind_cols(newdat_vis, pred_vis) %>%
  mutate(
    x_visitors_per_min = x * vis_sd + vis_mean   # back-transform x
  )

dat_plot_vis <- c.pl.final.table.4 %>%
  mutate(
    x_visitors_per_min = mean.visitors.per.minute
  )

vis_plot <- ggplot(plot_df_vis, aes(x = x_visitors_per_min, y = Estimate)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5),
              fill = "#d4e9ff", alpha = 0.7, colour = NA) +
  geom_line(size = 1.1, colour = "black") +
  geom_point(data = dat_plot_vis,
             aes(x = x_visitors_per_min, y = mean_seedset_round),
             inherit.aes = FALSE,
             size = 2.2, alpha = 0.85) +
  scale_x_continuous(
    name   = "Visitors per minute",
    breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    trans  = log1p_trans(),
    breaks = c(0, 1000, 5000, 10000, 20000, 40000),
    labels = label_comma(big.mark = " "),
    name   = "Expected seed set (log1p)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90"),
    axis.title         = element_text(face = "bold", size = 16),
    axis.text          = element_text(size = 13),
    plot.tag           = element_text(size = 18, face = "bold")
  )

vis_plot

















# from flowering.visited, before summarising to final.table
vis_func_mean <- mean(flowering.visited$total.func, na.rm = TRUE)
vis_func_sd   <- sd(flowering.visited$total.func,   na.rm = TRUE)

m   <- bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt
dat <- c.pl.final.table.4

zr    <- range(dat$z, na.rm = TRUE)
z_seq <- seq(zr[1], zr[2], length.out = 300)

newdat_func <- tibble(
  x  = 0,        # visitation at its mean (z-scale)
  sx = sx_med,
  z  = z_seq,    # vary functional richness
  sz = sz_med,
  species = NA
)

pred_func <- fitted(
  m,
  newdata    = newdat_func,
  re_formula = NA,
  summary    = TRUE
) |> as.data.frame()

plot_df_func <- bind_cols(newdat_func, pred_func) %>%
  mutate(
    func_richness = z * vis_func_sd + vis_func_mean
  )

dat_plot_func <- c.pl.final.table.4 %>%
  mutate(
    func_richness = mean.func   # original functional richness per species×elevation
  )

func_plot <- ggplot(plot_df_func, aes(x = func_richness, y = Estimate)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5),
              fill = "#d4e9ff", alpha = 0.7, colour = NA) +
  geom_line(size = 1.1, colour = "black") +
  geom_point(data = dat_plot_func,
             aes(x = func_richness, y = mean_seedset_round),
             inherit.aes = FALSE,
             size = 2.2, alpha = 0.85) +
  scale_x_continuous(
    name   = "Functional-group richness",
    breaks = pretty(dat_plot_func$func_richness),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    name   = "Expected seed set",
    breaks = c(0, 50, 100, 150, 200, 250, 300, 350),
    limits = c(0, NA)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90"),
    axis.title         = element_text(face = "bold", size = 16),
    axis.text          = element_text(size = 13),
    plot.tag           = element_text(size = 18, face = "bold")
  )

func_plot




# 1) Add panel tags
vis_plot_tagged <- vis_plot +
  labs(tag = "A)") +
  theme(plot.tag = element_text(size = 18, face = "bold"))

func_plot_tagged <- func_plot +
  labs(tag = "B)") +
  theme(plot.tag = element_text(size = 18, face = "bold"))

# 2) Combine into one figure
fig5_combined <- vis_plot_tagged + func_plot_tagged +
  plot_layout(nrow = 1)

# Optional: look at it in RStudio
print(fig5_combined)

# 3) Save as a single PDF
ggsave(
  filename = "figs/fig5_combined.pdf",
  plot     = fig5_combined,
  width    = 10,   # tweak as you like
  height   = 5
)
























  ###* IS the mean seedset round actually so 
View(c.pl.final.table.4)

form_all <- bf(
  mean_seedset_round | weights(seedset_weight_12) ~
    #me(x, sx) +
    #me(x2, sx2) +
    #me(y, sy) +
    #me(y2, sy2) +
    #me(z, sz) + 
    1 +
    (1|species),
  zi ~ 1
)


mu0 <- log(mean(c.pl.final.table.4$mean_seedset_round, na.rm = TRUE) + 1e-8)

mu0

pri_all <- c(
  prior(normal(3.688879, 1), class = "Intercept"),
  #prior(normal(0, 1), class = "b"),
  prior(student_t(3, 0, 2.5), class = "sd", group = "species")
  #prior(normal(0, 0.5), class = "b", coef = "mex2sx2"),
  #prior(normal(0, 0.5), class = "b", coef = "mey2sy2")
)


bayesian_seedset_mean_corrected_scaled_null_3_zinb_k_opt <- brm(
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


saveRDS(bayesian_seedset_mean_corrected_scaled_null_3_zinb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_null_3_zinb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_null_3_zinb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_null_3_zinb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_null_3_zinb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_null_3_zinb_k_opt)

seedset_scaled_k5_null_3_zinb <- kfold(bayesian_seedset_mean_corrected_scaled_null_3_zinb_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
seedset_scaled_k5_null_3_zinb                   # prints the elpd_kfold etc.
print(seedset_scaled_k5_null_3_zinb)
summary(seedset_scaled_k5_null_3_zinb)          # shows estimates and pointwise columns
seedset_scaled_k5_null_3_zinb$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(seedset_scaled_k5_null_3_zinb$pointwise)   # fold-wise contributions per observation
saveRDS(seedset_scaled_k5_null_3_zinb, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_null_3_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_null_3_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_null_3_zinb_k_opt.rds")
print(seedset_scaled_k5_null_3_zinb)


saveRDS(bayesian_seedset_mean_corrected_scaled_13_3_zinb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_13_3_zinb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_13_3_zinb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_13_3_zinb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_13_3_zinb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_13_3_zinb_k_opt)

seedset_scaled_k5_13_3_zinb <- kfold(bayesian_seedset_mean_corrected_scaled_13_3_zinb_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
seedset_scaled_k5_13_3_zinb                   # prints the elpd_kfold etc.
print(seedset_scaled_k5_13_3_zinb)
summary(seedset_scaled_k5_13_3_zinb)          # shows estimates and pointwise columns
seedset_scaled_k5_13_3_zinb$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(seedset_scaled_k5_13_3_zinb$pointwise)   # fold-wise contributions per observation
saveRDS(seedset_scaled_k5_13_3_zinb, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_13_3_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_13_3_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_13_3_zinb_k_opt.rds")
print(seedset_scaled_k5_13_3_zinb)

pp_check(bayesian_seedset_mean_corrected_scaled_13_3_zinb_k_opt, type = "dens_overlay")
pp_check(bayesian_seedset_mean_corrected_scaled_13_3_zinb_k_opt, type = "stat", stat = "mean")
pp_check(bayesian_seedset_mean_corrected_scaled_13_3_zinb_k_opt, type = "stat", stat = "sd")




saveRDS(bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt.rds")
bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt.rds")
summary(bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt)
loo(bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt)

seedset_scaled_k5_1_2_zinb <- kfold(bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
seedset_scaled_k5_1_2_zinb                   # prints the elpd_kfold etc.
print(seedset_scaled_k5_1_2_zinb)
summary(seedset_scaled_k5_1_2_zinb)          # shows estimates and pointwise columns
seedset_scaled_k5_1_2_zinb$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(seedset_scaled_k5_1_2_zinb$pointwise)   # fold-wise contributions per observation
saveRDS(seedset_scaled_k5_1_2_zinb, "../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_scaled_1_2_zinb_k_opt.rds")
print(seedset_scaled_k5_1_2_zinb)

loo_compare(seedset_scaled_k5_13_3_zinb, seedset_scaled_k5_null_3_zinb)
, seedset_scaled_k5_1_2_zinb)



m <- bayesian_seedset_mean_corrected_scaled_13_3_zinb_k_opt

ce <- conditional_effects(
  m, effects = "x",
  re_formula = NA,                # population-level
  method = "posterior_epred",     # expected counts incl. zero inflation
  ndraws = 2000
)
plot(ce, points = TRUE)




















# 
# View(c.pl.final.table.4)
# hist(c.pl.final.table.4$mean_seedset_round)
# 
# fit_all_scaled_log_seedset_gaus <- brm(
#   form_all, 
#   data = c.pl.final.table.4, 
#   family = gaussian(),
#   prior = pri_all, 
#   chains = 4, 
#   cores = 4,
#   iter = 2000,
#   warmup = 1000,
#   control = list(max_treedepth = 20, adapt_delta = 0.99),
#   save_pars = save_pars(all = TRUE, latent = TRUE),
#   seed = 1234
# )
# 
# saveRDS(fit_all_zib, file = "../brms_models/fit_all_zib.rds")
# fit_all_zib <- readRDS("../brms_models/fit_all_zib.rds")
# summary(fit_all_zib)
# 
# 
# 
# 
# form_nb <- bf(
#   mean_seedset_round | weights(seedset_weight_12) ~
#     me(x,  sx) + #1
#     #me(x2, sx2) + #2
#     me(y,  sy) + #3
#     #me(y2, sy2) + #4
#     me(z,  sz) + #5
#     #me(z2, sz2) + #6
#     (1 | species)
# )
# 
# 
# pri_nb_simple <- c(
#   prior(normal(0, 5),  class = "Intercept"),   # log-mean count; very weak
#   prior(normal(0, 1),  class = "b"),           # all slopes (z-scaled inputs)
#   prior(student_t(3, 0, 2.5), class = "sd", group = "species"),
#   #prior(normal(0, 0.5), class = "b", coef = "mex2sx2"),
#   #prior(normal(0, 0.5), class = "b", coef = "mey2sy2"),
#   #prior(normal(0, 0.5), class = "b", coef = "mez2sz2"),
#   prior(normal(0, 1),  class = "shape")        # log(shape)
# )
# 
# 
# bayesian_seedset_mean_corrected_135_nb <- brm(
#   form_nb, data = c.pl.final.table.4,
#   family = negbinomial(), 
#   prior = pri_nb_simple,
#   chains = 4, cores = 4, 
#   iter = 8000, warmup = 3000,
#   control = list(adapt_delta = 0.99, max_treedepth = 20),
#   seed = 1234
# )
# 
# saveRDS(fit_nb, file = "../brms_models/bayesian_seedset_mean_corrected_all_nb.rds")
# bayesian_seedset_mean_corrected_all_nb <- readRDS("../brms_models/bayesian_seedset_mean_corrected_all_nb.rds")
# summary(bayesian_seedset_mean_corrected_all_nb)
# 
# saveRDS(bayesian_seedset_mean_corrected_135_nb, file = "../brms_models/bayesian_seedset_mean_corrected_135_nb.rds")
# bayesian_seedset_mean_corrected_135_nb <- readRDS("../brms_models/bayesian_seedset_mean_corrected_135_nb.rds")
# summary(bayesian_seedset_mean_corrected_135_nb)
# 
# loo_compare(bayesian_seedset_mean_corrected_all_nb, bayesian_seedset_mean_corrected_135_nb)
# 
# 
# View(c.pl.final.table.4)

eps <- 1e-8
c.pl.final.table.4 <- c.pl.final.table.4 %>%
  mutate(
    # the means you want to use as predictors
    vpm    = mean.visited.flowers.per.minute,
    morpho = mean.morpho,
    func   = mean.func,
    
    # SEs of those means: SD / sqrt(n_reps)
    vpm_se    = pmax(sd.visited.flowers.per.minute / sqrt(pmax(n_reps, 1L)), eps),
    morpho_se = pmax(sd.morpho                     / sqrt(pmax(n_reps, 1L)), eps),
    func_se   = pmax(sd.func                       / sqrt(pmax(n_reps, 1L)), eps)
  )

bf_all <- bf(
  mean_seedset_round | weights(seedset_weight_12) ~
    me(vpm,    vpm_se) +
    me(morpho, morpho_se) +
    me(func,   func_se) +
    (1 | species)
)
# 
# pri_nb_simple <- c(
#   prior(normal(0, 5),  class = "Intercept"),   # log-mean count; very weak
#   prior(normal(0, 1),  class = "b"),           # all slopes (z-scaled inputs)
#   prior(student_t(3, 0, 2.5), class = "sd", group = "species"),
#   #prior(normal(0, 0.5), class = "b", coef = "mex2sx2"),
#   #prior(normal(0, 0.5), class = "b", coef = "mey2sy2"),
#   #prior(normal(0, 0.5), class = "b", coef = "mez2sz2"),
#   prior(normal(0, 1),  class = "shape")        # log(shape)
# )
# 
# bayesian_seedset_mean_corrected_orig_135_nb <- brm(
#   bf_all, 
#   data = c.pl.final.table.4,
#   family = negbinomial(), 
#   prior = pri_nb_simple,
#   chains = 4, cores = 4, 
#   iter = 6000, warmup = 2000,
#   #moment_match = TRUE,
#   save_pars = save_pars(all = TRUE, latent = TRUE),
#   control = list(adapt_delta = 0.99, max_treedepth = 20),
#   seed = 1234
# )
# 
# View(c.pl.final.table.4)
# 
# saveRDS(bayesian_seedset_mean_corrected_orig_135_nb, file = "../brms_models/bayesian_seedset_mean_corrected_orig_135_nb.rds")
# bayesian_seedset_mean_corrected_orig_135_nb <- readRDS("../brms_models/bayesian_seedset_mean_corrected_orig_135_nb.rds")
# summary(bayesian_seedset_mean_corrected_orig_135_nb)
# loo(bayesian_seedset_mean_corrected_orig_135_nb)
# 
# saveRDS(bayesian_seedset_mean_corrected_orig_1_nb, file = "../brms_models/bayesian_seedset_mean_corrected_orig_1_nb.rds")
# bayesian_seedset_mean_corrected_orig_1_nb <- readRDS("../brms_models/bayesian_seedset_mean_corrected_orig_1_nb.rds")
# summary(bayesian_seedset_mean_corrected_orig_1_nb)
# loo(bayesian_seedset_mean_corrected_orig_1_nb)
# 
# saveRDS(bayesian_seedset_mean_corrected_orig_3_nb, file = "../brms_models/bayesian_seedset_mean_corrected_orig_3_nb.rds")
# bayesian_seedset_mean_corrected_orig_3_nb <- readRDS("../brms_models/bayesian_seedset_mean_corrected_orig_3_nb.rds")
# summary(bayesian_seedset_mean_corrected_orig_3_nb)
# loo(bayesian_seedset_mean_corrected_orig_3_nb)
# 
# saveRDS(bayesian_seedset_mean_corrected_orig_13_nb, file = "../brms_models/bayesian_seedset_mean_corrected_orig_13_nb.rds")
# bayesian_seedset_mean_corrected_orig_13_nb <- readRDS("../brms_models/bayesian_seedset_mean_corrected_orig_13_nb.rds")
# summary(bayesian_seedset_mean_corrected_orig_13_nb)
# loo(bayesian_seedset_mean_corrected_orig_13_nb)
# 
# saveRDS(bayesian_seedset_mean_corrected_orig_5_nb, file = "../brms_models/bayesian_seedset_mean_corrected_orig_5_nb.rds")
# bayesian_seedset_mean_corrected_orig_5_nb <- readRDS("../brms_models/bayesian_seedset_mean_corrected_orig_5_nb.rds")
# summary(bayesian_seedset_mean_corrected_orig_5_nb)
# loo(bayesian_seedset_mean_corrected_orig_5_nb)
# 
# ###* save next
# 
# loo_M135_mm <- loo(bayesian_seedset_mean_corrected_orig_135_nb, moment_match = TRUE)
# 
# loo_M135 <- loo(bayesian_seedset_mean_corrected_orig_135_nb, reloo = TRUE)
# 
# ###* when   chains = 4, cores = 4, iter = 6000, warmup = 2000,
# ###* then k fold (K = 10) last for more than hour and a half
# loo_M135_k <- kfold(bayesian_seedset_mean_corrected_orig_135_nb, K = 10)
# 
# ###* How come kfold takes so long? Would it be better to increase chains and cores 
# ###* while decreasing iter, or decrease chanins and cores while increasing iter?
# 
# loo_M1   <- loo(bayesian_seedset_mean_corrected_orig_1_nb,   reloo = TRUE)
# loo_compare(loo_M135, loo_M1)



# ###* Make K better
# ###* 
# options(brms.backend = "cmdstanr")
# 
# bf_all <- bf(
#   mean_seedset_round | weights(seedset_weight_12) ~
#     #me(vpm,    vpm_se) +
#     #me(morpho, morpho_se) +
#     #me(func,   func_se) +
#     1 +
#     (1 | species)
# )
# 
# pri_nb_simple <- c(
#   prior(normal(0, 5),  class = "Intercept"),   # log-mean count; very weak
#   prior(normal(0, 1),  class = "b"),           # all slopes (z-scaled inputs)
#   prior(student_t(3, 0, 2.5), class = "sd", group = "species"),
#   #prior(normal(0, 0.5), class = "b", coef = "mex2sx2"),
#   #prior(normal(0, 0.5), class = "b", coef = "mey2sy2"),
#   #prior(normal(0, 0.5), class = "b", coef = "mez2sz2"),
#   prior(normal(0, 1),  class = "shape")        # log(shape)
# )
# 
# bayesian_seedset_mean_corrected_orig_null_nb_k_opt <- brm(
#   bf_all, 
#   data = c.pl.final.table.4,
#   family = negbinomial(), 
#   prior = pri_nb_simple,
#   chains = 5, cores = 5, 
#   iter = 2000, warmup = 1000,
#   threads = threading(4),
#   save_pars = save_pars(all = TRUE, latent = TRUE),
#   control = list(adapt_delta = 0.95, max_treedepth = 20),
#   seed = 1234
# )
# 
# saveRDS(bayesian_seedset_mean_corrected_orig_135_nb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_orig_135_nb_k_opt.rds")
# bayesian_seedset_mean_corrected_orig_135_nb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_orig_135_nb_k_opt.rds")
# summary(bayesian_seedset_mean_corrected_orig_135_nb_k_opt)
# loo(bayesian_seedset_mean_corrected_orig_135_nb_k_opt)
# 
# k5 <- kfold(bayesian_seedset_mean_corrected_orig_135_nb_k_opt, K = 5, cores = 5)
# ###* this lasts more than 15-20 minutes
# k5                   # prints the elpd_kfold etc.
# print(k5)
# summary(k5)          # shows estimates and pointwise columns
# k5$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
# head(k5$pointwise)   # fold-wise contributions per observation
# saveRDS(k5, "../brms_models/kfold_bayesian_seedset_mean_corrected_orig_135_nb_k_opt.rds")
# # ...later...
# kfold_bayesian_seedset_mean_corrected_orig_135_nb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_orig_135_nb_k_opt.rds")
# print(k5)
# 
# saveRDS(bayesian_seedset_mean_corrected_orig_1_nb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_orig_1_nb_k_opt.rds")
# bayesian_seedset_mean_corrected_orig_1_nb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_orig_1_nb_k_opt.rds")
# summary(bayesian_seedset_mean_corrected_orig_1_nb_k_opt)
# loo(bayesian_seedset_mean_corrected_orig_1_nb_k_opt)
# 
# seedset_k5_1 <- kfold(bayesian_seedset_mean_corrected_orig_1_nb_k_opt, K = 5, cores = 5)
# ###* this lasts more than 15-20 minutes
# seedset_k5_1                  # prints the elpd_kfold etc.
# print(seedset_k5_1)
# summary(seedset_k5_1)          # shows estimates and pointwise columns
# seedset_k5_1$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
# head(seedset_k5_1$pointwise)   # fold-wise contributions per observation
# saveRDS(seedset_k5_1, "../brms_models/kfold_bayesian_seedset_mean_corrected_orig_1_nb_k_opt.rds")
# kfold_bayesian_seedset_mean_corrected_orig_1_nb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_orig_1_nb_k_opt.rds")
# print(kfold_bayesian_seedset_mean_corrected_orig_1_nb_k_opt)
# 
# 
# saveRDS(bayesian_seedset_mean_corrected_orig_3_nb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_orig_3_nb_k_opt.rds")
# bayesian_seedset_mean_corrected_orig_3_nb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_orig_3_nb_k_opt.rds")
# summary(bayesian_seedset_mean_corrected_orig_3_nb_k_opt)
# loo(bayesian_seedset_mean_corrected_orig_3_nb_k_opt)
# 
# seedset_k5_3 <- kfold(bayesian_seedset_mean_corrected_orig_3_nb_k_opt, K = 5, cores = 5)
# ###* this lasts more than 15-20 minutes
# seedset_k5_3                  # prints the elpd_kfold etc.
# print(seedset_k5_3)
# summary(seedset_k5_3)          # shows estimates and pointwise columns
# seedset_k5_3$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
# head(seedset_k5_3$pointwise)   # fold-wise contributions per observation
# saveRDS(seedset_k5_3, "../brms_models/kfold_bayesian_seedset_mean_corrected_orig_3_nb_k_opt.rds")
# kfold_bayesian_seedset_mean_corrected_orig_3_nb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_orig_3_nb_k_opt.rds")
# print(kfold_bayesian_seedset_mean_corrected_orig_3_nb_k_opt)
# 
# 
# saveRDS(bayesian_seedset_mean_corrected_orig_5_nb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_orig_5_nb_k_opt.rds")
# bayesian_seedset_mean_corrected_orig_5_nb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_orig_5_nb_k_opt.rds")
# summary(bayesian_seedset_mean_corrected_orig_5_nb_k_opt)
# loo(bayesian_seedset_mean_corrected_orig_5_nb_k_opt)
# 
# seedset_k5_5 <- kfold(bayesian_seedset_mean_corrected_orig_5_nb_k_opt, K = 5, cores = 5)
# ###* this lasts more than 15-20 minutes
# seedset_k5_5                  # prints the elpd_kfold etc.
# print(seedset_k5_5)
# summary(seedset_k5_5)          # shows estimates and pointwise columns
# seedset_k5_5$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
# head(seedset_k5_5$pointwise)   # fold-wise contributions per observation
# saveRDS(seedset_k5_5, "../brms_models/kfold_bayesian_seedset_mean_corrected_orig_5_nb_k_opt.rds")
# kfold_bayesian_seedset_mean_corrected_orig_5_nb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_orig_5_nb_k_opt.rds")
# print(kfold_bayesian_seedset_mean_corrected_orig_5_nb_k_opt)
# 
# 
# saveRDS(bayesian_seedset_mean_corrected_orig_35_nb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_orig_35_nb_k_opt.rds")
# bayesian_seedset_mean_corrected_orig_35_nb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_orig_35_nb_k_opt.rds")
# summary(bayesian_seedset_mean_corrected_orig_35_nb_k_opt)
# loo(bayesian_seedset_mean_corrected_orig_35_nb_k_opt)
# 
# seedset_k5_35 <- kfold(bayesian_seedset_mean_corrected_orig_35_nb_k_opt, K = 5, cores = 5)
# ###* this lasts more than 15-20 minutes
# seedset_k5_35                  # prints the elpd_kfold etc.
# print(seedset_k5_35)
# summary(seedset_k5_35)          # shows estimates and pointwise columns
# seedset_k5_35$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
# head(seedset_k5_35$pointwise)   # fold-wise contributions per observation
# saveRDS(seedset_k5_35, "../brms_models/kfold_bayesian_seedset_mean_corrected_orig_35_nb_k_opt.rds")
# kfold_bayesian_seedset_mean_corrected_orig_35_nb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_orig_35_nb_k_opt.rds")
# print(kfold_bayesian_seedset_mean_corrected_orig_35_nb_k_opt)
# 
# 
# saveRDS(bayesian_seedset_mean_corrected_orig_13_nb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_orig_13_nb_k_opt.rds")
# bayesian_seedset_mean_corrected_orig_13_nb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_orig_13_nb_k_opt.rds")
# summary(bayesian_seedset_mean_corrected_orig_13_nb_k_opt)
# loo(bayesian_seedset_mean_corrected_orig_13_nb_k_opt)
# 
# seedset_k5_13 <- kfold(bayesian_seedset_mean_corrected_orig_13_nb_k_opt, K = 5, cores = 5)
# ###* this lasts more than 15-20 minutes
# seedset_k5_13                  # prints the elpd_kfold etc.
# print(seedset_k5_13)
# summary(seedset_k5_13)          # shows estimates and pointwise columns
# seedset_k5_13$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
# head(seedset_k5_13$pointwise)   # fold-wise contributions per observation
# saveRDS(seedset_k5_13, "../brms_models/kfold_bayesian_seedset_mean_corrected_orig_13_nb_k_opt.rds")
# kfold_bayesian_seedset_mean_corrected_orig_13_nb_k_opt <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_orig_13_nb_k_opt.rds")
# print(kfold_bayesian_seedset_mean_corrected_orig_13_nb_k_opt)
# 
# 
# saveRDS(bayesian_seedset_mean_corrected_orig_null_nb_k_opt, file = "../brms_models/bayesian_seedset_mean_corrected_orig_null_nb_k_opt.rds")
# bayesian_seedset_mean_corrected_orig_null_nb_k_opt <- readRDS("../brms_models/bayesian_seedset_mean_corrected_orig_null_nb_k_opt.rds")
# summary(bayesian_seedset_mean_corrected_orig_null_nb_k_opt)
# loo(bayesian_seedset_mean_corrected_orig_null_nb_k_opt)
# 
# seedset_k5_null <- kfold(bayesian_seedset_mean_corrected_orig_null_nb_k_opt, K = 5, cores = 5)
# ###* this lasts more than 15-20 minutes
# seedset_k5_null                   # prints the elpd_kfold etc.
# print(seedset_k5_null)
# summary(seedset_k5_null)          # shows estimates and pointwise columns
# seedset_k5_null$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
# head(seedset_k5_null$pointwise)   # fold-wise contributions per observation
# saveRDS(seedset_k5_null, "../brms_models/kfold_bayesian_seedset_mean_corrected_orig_135_nb_k_opt.rds")
# # ...later...
# seedset_k5_null <- readRDS("../brms_models/kfold_bayesian_seedset_mean_corrected_orig_null_nb_k_opt.rds")
# print(seedset_k5_null )
# 
# 
# 
# 
# 
# loo_compare(k5, seedset_k5_1, seedset_k5_3, seedset_k5_5, seedset_k5_35, seedset_k5_13, seedset_k5_null)
# ###*
# ###*
# ###* RESULTS FOR SEEDSET



# c.pl.final.table.4 <- c.pl.final.table.4 %>%
#   mutate(y_log1p = log1p(mean_seedset_round))
# 
# form_gaus <- bf(
#   y_log1p | weights(seedset_weight_12) ~
#     me(vpm, vpm_se) + me(morpho, morpho_se) + me(func, func_se) +
#     (1|species)
# )
# 
# hist(c.pl.final.table.4$y_log1p)
# 
# pri_stud <- c(
#   prior(normal(0, 5), class="Intercept"),
#   prior(normal(0, 1), class="b"),
#   prior(student_t(3,0,2.5), class="sd", group="species"),
#   prior(student_t(3,0,2.5), class="sigma")
# )
# 
# bayesian_seedset_mean_corrected_log_135_gaus <- brm(
#   form_gaus, 
#   data = c.pl.final.table.4,
#   family=gaussian(), 
#   prior = pri_stud,
#   chains= 4, cores = 4,
#   iter = 8000, warmup = 3000,
#   save_pars = save_pars(all = TRUE, latent = TRUE),
#   control = list(adapt_delta = 0.99, max_treedepth = 20),
#   seed = 1234
#   )
# 
# saveRDS(bayesian_seedset_mean_corrected_log_135_stud, file = "../brms_models/bayesian_seedset_mean_corrected_log_135_stud.rds")
# bayesian_seedset_mean_corrected_log_135_stud <- readRDS("../brms_models/bayesian_seedset_mean_corrected_log_135_stud.rds")
# summary(bayesian_seedset_mean_corrected_log_135_stud)
# 
# pp_check(bayesian_seedset_mean_corrected_log_135_stud)      # posterior predictive fit
# bayes_R2(bayesian_seedset_mean_corrected_log_135_stud)
# loo(bayesian_seedset_mean_corrected_log_135_stud)           # compare to your other seeds models
# summary(bayesian_seedset_mean_corrected_log_135_stud)


library(ggplot2)

ggplot(c.pl.final.table.4,
       aes(x = mean.visited.flowers.per.minute,
           y = mean_seedset_round)) +
  geom_point()

ggplot(c.pl.final.table.4,
       aes(x = mean.morpho,
           y = mean_seedset_round)) +
  geom_point()

ggplot(c.pl.final.table.4,
       aes(x = mean.func,
           y = mean_seedset_round)) +
  geom_point()

ggplot(c.pl.final.table.4,
       aes(x = mean.visited.flowers.per.minute,
           y = scaled_seedset)) +
  geom_point()

ggplot(c.pl.final.table.4,
       aes(x = mean.morpho,
           y = scaled_seedset)) +
  geom_point()

ggplot(c.pl.final.table.4,
       aes(x = mean.func,
           y = scaled_seedset)) +
  geom_point()

View(c.pl.final.table.4)



























###* POLLEN LIMITATION
###* 
###* HAVE TO RERUN EVERYTHHING, THE MOTHERFUCKER WAS SCALED - CHNAGE
###* THE PREDICTORS FROM SCALED TO NOT SCALED
###* 
library(brms)
options(brms.backend = "cmdstanr")

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
    sz = pmax(z_func_sd   / sqrt(pmax(n_reps,  1L)), eps),
    
    # --- quadratic terms and their delta-method SEs ---
    x2  = x^2,  sx2 = pmax(2 * abs(x) * sx, eps),
    y2  = y^2,  sy2 = pmax(2 * abs(y) * sy, eps),
    z2  = z^2,  sz2 = pmax(2 * abs(z) * sz, eps)
  )

###* cOLENIARITY?
View(c.pl.final.table.5)

# # 1) Simple diagnostics
# with(c.pl.final.table.5, {
#   cat("n =", nrow(c.pl.final.table.5), "\n")
#   print(summary(y <- mean.morpho.scaled))
#   print(summary(z <- mean.func.scaled))
# })
# 
# # 2) Correlations (scaling doesn’t affect r)
# with(c.pl.final.table.5, {
#   print(cor.test(mean.morpho.scaled, mean.func.scaled, method = "pearson"))
#   print(cor.test(mean.morpho.scaled, mean.func.scaled, method = "spearman"))
# })
# 
# # 3) VIF / condition index using a dummy Gaussian model
# m_dummy <- lm(mean_PL_index ~ mean.morpho.scaled + mean.func.scaled, data = c.pl.final.table.5)
# check_collinearity(m_dummy)   # shows VIF and condition number
# 
# # 4) Quick plot
# library(ggplot2)
# ggplot(c.pl.final.table.5, aes(mean.morpho.scaled, mean.func.scaled)) +
#   geom_point(size = 3) +
#   geom_smooth(method = "lm", se = FALSE) +
#   labs(x = "Morphospecies richness (scaled)", y = "Functional-group richness (scaled)")



# 3) Formula with measurement error and a species random intercept.
#    (zoib handles exact 0s and 1s; you can keep it simple with intercepts for phi/zi)
bf_all <- bf(
  mean_seedset_round | weights(seedset_weight_12) ~
    #me(vpm,    vpm_se) +
    #me(morpho, morpho_se) +
    #me(func,   func_se) +
    me(x, sx) +
    me(x2, sx2) +
    me(y, sy) +
    me(y2, sy2) +
    1 +
    (1 | species),
  phi ~ 1,
  zoi ~ 1,
  coi ~ 1
)

pri_pl <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 1), class = "Intercept"),
  prior(student_t(3, 0, 2.5), class = "sd", group = "species"),
  prior(normal(0, 1), class = "Intercept", dpar = "phi"),
  prior(normal(0, 1), class = "Intercept", dpar = "zoi"),
  prior(normal(0, 1), class = "Intercept", dpar = "coi")
)

bayesian_pl_mean_corrected_orig_null_zoib_k_opt <- brm(
  formula = form_pl,
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

saveRDS(bayesian_pl_mean_corrected_orig_13_zoib_k_opt, file = "../brms_models/bayesian_pl_mean_corrected_orig_13_zoib_k_opt.rds")
bayesian_pl_mean_corrected_orig_13_zoib_k_opt <- readRDS("../brms_models/bayesian_pl_mean_corrected_orig_13_zoib_k_opt.rds")
summary(bayesian_pl_mean_corrected_orig_13_zoib_k_opt)
loo(bayesian_pl_mean_corrected_orig_13_zoib_k_opt)

pl_k5_13 <- kfold(bayesian_pl_mean_corrected_orig_13_zoib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
pl_k5_13                  # prints the elpd_kfold etc.
print(pl_k5_13)
summary(pl_k5_13)          # shows estimates and pointwise columns
pl_k5_13$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(pl_k5_13$pointwise)   # fold-wise contributions per observation
saveRDS(pl_k5_13, "../brms_models/kfold_bayesian_pl_mean_corrected_orig_13_zoib_k_opt.rds")
kfold_bayesian_pl_mean_corrected_orig_13_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_pl_mean_corrected_orig_13_zoib_k_opt.rds")
print(kfold_bayesian_pl_mean_corrected_orig_13_zoib_k_opt)


saveRDS(bayesian_pl_mean_corrected_orig_1_zoib_k_opt, file = "../brms_models/bayesian_pl_mean_corrected_orig_1_zoib_k_opt.rds")
bayesian_pl_mean_corrected_orig_1_zoib_k_opt <- readRDS("../brms_models/bayesian_pl_mean_corrected_orig_1_zoib_k_opt.rds")
summary(bayesian_pl_mean_corrected_orig_1_zoib_k_opt)
loo(bayesian_pl_mean_corrected_orig_1_zoib_k_opt)

pl_k5_1 <- kfold(bayesian_pl_mean_corrected_orig_1_zoib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
pl_k5_1                  # prints the elpd_kfold etc.
print(pl_k5_1)
summary(pl_k5_1)          # shows estimates and pointwise columns
pl_k5_1$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(pl_k5_1$pointwise)   # fold-wise contributions per observation
saveRDS(pl_k5_1, "../brms_models/kfold_bayesian_pl_mean_corrected_orig_1_zoib_k_opt.rds")
kfold_bayesian_pl_mean_corrected_orig_1_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_pl_mean_corrected_orig_1_zoib_k_opt.rds")
print(kfold_bayesian_pl_mean_corrected_orig_1_zoib_k_opt)


saveRDS(bayesian_pl_mean_corrected_orig_3_zoib_k_opt, file = "../brms_models/bayesian_pl_mean_corrected_orig_3_zoib_k_opt.rds")
bayesian_pl_mean_corrected_orig_3_zoib_k_opt <- readRDS("../brms_models/bayesian_pl_mean_corrected_orig_3_zoib_k_opt.rds")
summary(bayesian_pl_mean_corrected_orig_3_zoib_k_opt)
loo(bayesian_pl_mean_corrected_orig_3_zoib_k_opt)

pl_k5_3 <- kfold(bayesian_pl_mean_corrected_orig_3_zoib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
pl_k5_3                  # prints the elpd_kfold etc.
print(pl_k5_3)
summary(pl_k5_3)          # shows estimates and pointwise columns
pl_k5_3$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(pl_k5_3$pointwise)   # fold-wise contributions per observation
saveRDS(pl_k5_3, "../brms_models/kfold_bayesian_pl_mean_corrected_orig_3_zoib_k_opt.rds")
kfold_bayesian_pl_mean_corrected_orig_3_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_pl_mean_corrected_orig_3_zoib_k_opt.rds")
print(kfold_bayesian_pl_mean_corrected_orig_3_zoib_k_opt)


saveRDS(bayesian_pl_mean_corrected_orig_null_zoib_k_opt, file = "../brms_models/bayesian_pl_mean_corrected_orig_null_zoib_k_opt.rds")
bayesian_pl_mean_corrected_orig_null_zoib_k_opt <- readRDS("../brms_models/bayesian_pl_mean_corrected_orig_null_zoib_k_opt.rds")
summary(bayesian_pl_mean_corrected_orig_null_zoib_k_opt)
loo(bayesian_pl_mean_corrected_orig_null_zoib_k_opt)

pl_k5_null <- kfold(bayesian_pl_mean_corrected_orig_null_zoib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
pl_k5_null                  # prints the elpd_kfold etc.
print(pl_k5_null)
summary(pl_k5_null)          # shows estimates and pointwise columns
pl_k5_null$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(pl_k5_null$pointwise)   # fold-wise contributions per observation
saveRDS(pl_k5_null, "../brms_models/kfold_bayesian_pl_mean_corrected_orig_null_zoib_k_opt.rds")
kfold_bayesian_pl_mean_corrected_orig_null_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_pl_mean_corrected_orig_null_zoib_k_opt.rds")
print(kfold_bayesian_pl_mean_corrected_orig_null_zoib_k_opt)

loo_compare(pl_k5_13, pl_k5_1, pl_k5_3, pl_k5_null)
















###* POLLEN LIMITATION SCALED
View(c.pl.final.table.5)

bf_all <- bf(
  mean_PL_index | weights(PL_index_weight_12) ~
    me(x, sx) +
    #me(x2, sx2) +
    me(y, sy) +
    #(y2, sy2) +
    #1 +
    (1 | species),
  phi ~ 1,
  zoi ~ 1,
  coi ~ 1
)

pri_pl <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 1), class = "Intercept"),
  prior(student_t(3, 0, 2.5), class = "sd", group = "species"),
  prior(normal(0, 1), class = "Intercept", dpar = "phi"),
  prior(normal(0, 1), class = "Intercept", dpar = "zoi"),
  prior(normal(0, 1), class = "Intercept", dpar = "coi"),
  #prior(normal(0, 0.5), class = "b", coef = "mex2sx2"),
  #prior(normal(0, 0.5), class = "b", coef = "mey2sy2")
)


bayesian_pl_mean_corrected_scaled_1234_zoib_k_opt <- brm(
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

saveRDS(bayesian_pl_mean_corrected_scaled_13_zoib_k_opt, file = "../brms_models/bayesian_pl_mean_corrected_scaled_13_zoib_k_opt.rds")
bayesian_pl_mean_corrected_scaled_13_zoib_k_opt <- readRDS("../brms_models/bayesian_pl_mean_corrected_scaled_13_zoib_k_opt.rds")
summary(bayesian_pl_mean_corrected_scaled_13_zoib_k_opt)
loo(bayesian_pl_mean_corrected_scaled_13_zoib_k_opt)

pl_k5_scaled_13 <- kfold(bayesian_pl_mean_corrected_scaled_13_zoib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
pl_k5_scaled_13                  # prints the elpd_kfold etc.
print(pl_k5_scaled_13)
summary(pl_k5_scaled_13)          # shows estimates and pointwise columns
pl_k5_scaled_13$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(pl_k5_scaled_13$pointwise)   # fold-wise contributions per observation
saveRDS(pl_k5_scaled_13, "../brms_models/kfold_bayesian_pl_mean_corrected_scaled_13_zoib_k_opt.rds")
kfold_bayesian_pl_mean_corrected_scaled_13_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_pl_mean_corrected_scaled_13_zoib_k_opt.rds")
print(kfold_bayesian_pl_mean_corrected_scaled_13_zoib_k_opt)


saveRDS(bayesian_pl_mean_corrected_scaled_null_zoib_k_opt, file = "../brms_models/bayesian_pl_mean_corrected_scaled_null_zoib_k_opt.rds")
bayesian_pl_mean_corrected_scaled_null_zoib_k_opt <- readRDS("../brms_models/bayesian_pl_mean_corrected_scaled_null_zoib_k_opt.rds")
summary(bayesian_pl_mean_corrected_scaled_null_zoib_k_opt)
loo(bayesian_pl_mean_corrected_scaled_null_zoib_k_opt)

pl_k5_scaled_null <- kfold(bayesian_pl_mean_corrected_scaled_null_zoib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
pl_k5_scaled_null                  # prints the elpd_kfold etc.
print(pl_k5_scaled_null)
summary(pl_k5_scaled_null)          # shows estimates and pointwise columns
pl_k5_scaled_null$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(pl_k5_scaled_null$pointwise)   # fold-wise contributions per observation
saveRDS(pl_k5_scaled_null, "../brms_models/kfold_bayesian_pl_mean_corrected_scaled_null_zoib_k_opt.rds")
kfold_bayesian_pl_mean_corrected_scaled_null_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_pl_mean_corrected_scaled_null_zoib_k_opt.rds")
print(kfold_bayesian_pl_mean_corrected_scaled_null_zoib_k_opt)


saveRDS(bayesian_pl_mean_corrected_scaled_1_zoib_k_opt, file = "../brms_models/bayesian_pl_mean_corrected_scaled_1_zoib_k_opt.rds")
bayesian_pl_mean_corrected_scaled_1_zoib_k_opt <- readRDS("../brms_models/bayesian_pl_mean_corrected_scaled_1_zoib_k_opt.rds")
summary(bayesian_pl_mean_corrected_scaled_1_zoib_k_opt)
loo(bayesian_pl_mean_corrected_scaled_1_zoib_k_opt)

pl_k5_scaled_1 <- kfold(bayesian_pl_mean_corrected_scaled_1_zoib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
pl_k5_scaled_1                  # prints the elpd_kfold etc.
print(pl_k5_scaled_1)
summary(pl_k5_scaled_1)          # shows estimates and pointwise columns
pl_k5_scaled_1$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(pl_k5_scaled_1$pointwise)   # fold-wise contributions per observation
saveRDS(pl_k5_scaled_1, "../brms_models/kfold_bayesian_pl_mean_corrected_scaled_1_zoib_k_opt.rds")
kfold_bayesian_pl_mean_corrected_scaled_1_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_pl_mean_corrected_scaled_1_zoib_k_opt.rds")
print(kfold_bayesian_pl_mean_corrected_scaled_1_zoib_k_opt)


saveRDS(bayesian_pl_mean_corrected_scaled_3_zoib_k_opt, file = "../brms_models/bayesian_pl_mean_corrected_scaled_3_zoib_k_opt.rds")
bayesian_pl_mean_corrected_scaled_3_zoib_k_opt <- readRDS("../brms_models/bayesian_pl_mean_corrected_scaled_3_zoib_k_opt.rds")
summary(bayesian_pl_mean_corrected_scaled_3_zoib_k_opt)
loo(bayesian_pl_mean_corrected_scaled_3_zoib_k_opt)

pl_k5_scaled_3 <- kfold(bayesian_pl_mean_corrected_scaled_3_zoib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
pl_k5_scaled_3                  # prints the elpd_kfold etc.
print(pl_k5_scaled_3)
summary(pl_k5_scaled_3)          # shows estimates and pointwise columns
pl_k5_scaled_3$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(pl_k5_scaled_3$pointwise)   # fold-wise contributions per observation
saveRDS(pl_k5_scaled_3, "../brms_models/kfold_bayesian_pl_mean_corrected_scaled_3_zoib_k_opt.rds")
kfold_bayesian_pl_mean_corrected_scaled_3_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_pl_mean_corrected_scaled_3_zoib_k_opt.rds")
print(kfold_bayesian_pl_mean_corrected_scaled_3_zoib_k_opt)


saveRDS(bayesian_pl_mean_corrected_scaled_1234_zoib_k_opt, file = "../brms_models/bayesian_pl_mean_corrected_scaled_1234_zoib_k_opt.rds")
bayesian_pl_mean_corrected_scaled_1234_zoib_k_opt <- readRDS("../brms_models/bayesian_pl_mean_corrected_scaled_1234_zoib_k_opt.rds")
summary(bayesian_pl_mean_corrected_scaled_1234_zoib_k_opt)
loo(bayesian_pl_mean_corrected_scaled_1234_zoib_k_opt)

pl_k5_scaled_1234 <- kfold(bayesian_pl_mean_corrected_scaled_1234_zoib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
pl_k5_scaled_1234                  # prints the elpd_kfold etc.
print(pl_k5_scaled_1234)
summary(pl_k5_scaled_1234)          # shows estimates and pointwise columns
pl_k5_scaled_1234$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(pl_k5_scaled_1234$pointwise)   # fold-wise contributions per observation
saveRDS(pl_k5_scaled_1234, "../brms_models/kfold_bayesian_pl_mean_corrected_scaled_1234_zoib_k_opt.rds")
kfold_bayesian_pl_mean_corrected_scaled_1234_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_pl_mean_corrected_scaled_1234_zoib_k_opt.rds")
print(kfold_bayesian_pl_mean_corrected_scaled_1234_zoib_k_opt)


saveRDS(bayesian_pl_mean_corrected_scaled_12_zoib_k_opt, file = "../brms_models/bayesian_pl_mean_corrected_scaled_12_zoib_k_opt.rds")
bayesian_pl_mean_corrected_scaled_12_zoib_k_opt <- readRDS("../brms_models/bayesian_pl_mean_corrected_scaled_12_zoib_k_opt.rds")
summary(bayesian_pl_mean_corrected_scaled_12_zoib_k_opt)
loo(bayesian_pl_mean_corrected_scaled_12_zoib_k_opt)

pl_k5_scaled_12 <- kfold(bayesian_pl_mean_corrected_scaled_12_zoib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
pl_k5_scaled_12                  # prints the elpd_kfold etc.
print(pl_k5_scaled_12)
summary(pl_k5_scaled_12)          # shows estimates and pointwise columns
pl_k5_scaled_12$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(pl_k5_scaled_12$pointwise)   # fold-wise contributions per observation
saveRDS(pl_k5_scaled_12, "../brms_models/kfold_bayesian_pl_mean_corrected_scaled_12_zoib_k_opt.rds")
kfold_bayesian_pl_mean_corrected_scaled_12_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_pl_mean_corrected_scaled_12_zoib_k_opt.rds")
print(kfold_bayesian_pl_mean_corrected_scaled_12_zoib_k_opt)


saveRDS(bayesian_pl_mean_corrected_scaled_34_zoib_k_opt, file = "../brms_models/bayesian_pl_mean_corrected_scaled_34_zoib_k_opt.rds")
bayesian_pl_mean_corrected_scaled_34_zoib_k_opt <- readRDS("../brms_models/bayesian_pl_mean_corrected_scaled_34_zoib_k_opt.rds")
summary(bayesian_pl_mean_corrected_scaled_34_zoib_k_opt)
loo(bayesian_pl_mean_corrected_scaled_34_zoib_k_opt)

pl_k5_scaled_34 <- kfold(bayesian_pl_mean_corrected_scaled_34_zoib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
pl_k5_scaled_34                  # prints the elpd_kfold etc.
print(pl_k5_scaled_34)
summary(pl_k5_scaled_34)          # shows estimates and pointwise columns
pl_k5_scaled_34$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(pl_k5_scaled_34$pointwise)   # fold-wise contributions per observation
saveRDS(pl_k5_scaled_34, "../brms_models/kfold_bayesian_pl_mean_corrected_scaled_34_zoib_k_opt.rds")
kfold_bayesian_pl_mean_corrected_scaled_34_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_pl_mean_corrected_scaled_34_zoib_k_opt.rds")
print(kfold_bayesian_pl_mean_corrected_scaled_34_zoib_k_opt)

loo_compare(pl_k5_scaled_13, pl_k5_scaled_1, pl_k5_scaled_3, pl_k5_scaled_null, pl_k5_scaled_1234, , pl_k5_scaled_12, pl_k5_scaled_34)

library(loo)
set.seed(1234)
folds <- kfold_split_grouped(K = 5, x = c.pl.final.table.5$species)

k5_null <- kfold(bayesian_pl_mean_corrected_scaled_null_zoib_k_opt, K = 5, folds = folds)
k5_1    <- kfold(bayesian_pl_mean_corrected_scaled_1_zoib_k_opt,    K = 5, folds = folds)
k5_3    <- kfold(bayesian_pl_mean_corrected_scaled_3_zoib_k_opt,    K = 5, folds = folds)
k5_13   <- kfold(bayesian_pl_mean_corrected_scaled_13_zoib_k_opt,   K = 5, folds = folds)

loo_compare(k5_1, k5_13, k5_3, k5_null)























###* POLLEN LIMITATION fixed priors
View(c.pl.final.table.5)

bf_all <- bf(
  mean_PL_index | weights(PL_index_weight_12) ~
    #me(x, sx) +
    #me(x2, sx2) +
    #me(y, sy) +
    #me(y2, sy2) +
    #me(z, sz) +
    #me(z2, sz2) +
    1 +
    (1 | species),
  phi ~ 1,
  zoi ~ 1,
  coi ~ 1
)

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
  #prior(normal(0, 1), class = "b"),
  prior(student_t(3, 0, 2.5), class = "sd", group = "species"),
  prior(normal(0, 1), class = "Intercept", dpar = "phi"),
  prior(normal(0, 1), class = "Intercept", dpar = "zoi"),
  prior(normal(0, 1), class = "Intercept", dpar = "coi")
  #prior(normal(0, 0.5), class = "b", coef = "mex2sx2"),
  #prior(normal(0, 0.5), class = "b", coef = "mey2sy2")
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

saveRDS(bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2, file = "../brms_models/bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2.rds")
bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2 <- readRDS("../brms_models/bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2.rds")
summary(bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2)
loo(bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2)

pl_k5_scaled_13_2 <- kfold(bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
pl_k5_scaled_13_2                  # prints the elpd_kfold etc.
print(pl_k5_scaled_13_2)
summary(pl_k5_scaled_13_2)          # shows estimates and pointwise columns
pl_k5_scaled_13_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(pl_k5_scaled_13_2$pointwise)   # fold-wise contributions per observation
saveRDS(pl_k5_scaled_13_2, "../brms_models/kfold_bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2.rds")
kfold_bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2 <- readRDS("../brms_models/kfold_bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2.rds")
print(kfold_bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2)


saveRDS(bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2, file = "../brms_models/bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2.rds")
bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2 <- readRDS("../brms_models/bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2.rds")
summary(bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2)
loo(bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2)

pl_k5_scaled_null_2 <- kfold(bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
pl_k5_scaled_null_2                  # prints the elpd_kfold etc.
print(pl_k5_scaled_null_2)
summary(pl_k5_scaled_null_2)          # shows estimates and pointwise columns
pl_k5_scaled_null_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(pl_k5_scaled_null_2$pointwise)   # fold-wise contributions per observation
saveRDS(pl_k5_scaled_null_2, "../brms_models/kfold_bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2.rds")
kfold_bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2 <- readRDS("../brms_models/kfold_bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2.rds")
print(kfold_bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2)



saveRDS(bayesian_pl_mean_corrected_scaled_135_zoib_k_opt_2, file = "../brms_models/bayesian_pl_mean_corrected_scaled_135_zoib_k_opt_2.rds")
bayesian_pl_mean_corrected_scaled_135_zoib_k_opt_2 <- readRDS("../brms_models/bayesian_pl_mean_corrected_scaled_135_zoib_k_opt_2.rds")
summary(bayesian_pl_mean_corrected_scaled_135_zoib_k_opt_2)
loo(bayesian_pl_mean_corrected_scaled_135_zoib_k_opt_2)

pl_k5_scaled_135_2 <- kfold(bayesian_pl_mean_corrected_scaled_135_zoib_k_opt_2, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
pl_k5_scaled_135_2                  # prints the elpd_kfold etc.
print(pl_k5_scaled_135_2)
summary(pl_k5_scaled_135_2)          # shows estimates and pointwise columns
pl_k5_scaled_135_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(pl_k5_scaled_135_2$pointwise)   # fold-wise contributions per observation
saveRDS(pl_k5_scaled_135_2, "../brms_models/kfold_bayesian_pl_mean_corrected_scaled_135_zoib_k_opt_2.rds")
kfold_bayesian_pl_mean_corrected_scaled_135_zoib_k_opt_2 <- readRDS("../brms_models/kfold_bayesian_pl_mean_corrected_scaled_135_zoib_k_opt_2.rds")
print(kfold_bayesian_pl_mean_corrected_scaled_135_zoib_k_opt_2)


saveRDS(bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2, file = "../brms_models/bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2.rds")
bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2 <- readRDS("../brms_models/bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2.rds")
summary(bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2)
loo(bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2)

pl_k5_scaled_1_2 <- kfold(bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
pl_k5_scaled_1_2                  # prints the elpd_kfold etc.
print(pl_k5_scaled_1_2)
summary(pl_k5_scaled_1_2)          # shows estimates and pointwise columns
pl_k5_scaled_1_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(pl_k5_scaled_1_2$pointwise)   # fold-wise contributions per observation
saveRDS(pl_k5_scaled_1_2, "../brms_models/kfold_bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2.rds")
kfold_bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2 <- readRDS("../brms_models/kfold_bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2.rds")
print(kfold_bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2)


saveRDS(bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2, file = "../brms_models/bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2.rds")
bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2 <- readRDS("../brms_models/bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2.rds")
summary(bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2)
loo(bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2)

pl_k5_scaled_3_2 <- kfold(bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
pl_k5_scaled_3_2                  # prints the elpd_kfold etc.
print(pl_k5_scaled_3_2)
summary(pl_k5_scaled_3_2)          # shows estimates and pointwise columns
pl_k5_scaled_3_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(pl_k5_scaled_3_2$pointwise)   # fold-wise contributions per observation
saveRDS(pl_k5_scaled_3_2, "../brms_models/kfold_bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2.rds")
kfold_bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2 <- readRDS("../brms_models/kfold_bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2.rds")
print(kfold_bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2)



saveRDS(bayesian_pl_mean_corrected_scaled_5_zoib_k_opt_2, file = "../brms_models/bayesian_pl_mean_corrected_scaled_5_zoib_k_opt_2.rds")
bayesian_pl_mean_corrected_scaled_5_zoib_k_opt_2 <- readRDS("../brms_models/bayesian_pl_mean_corrected_scaled_5_zoib_k_opt_2.rds")
summary(bayesian_pl_mean_corrected_scaled_5_zoib_k_opt_2)
loo(bayesian_pl_mean_corrected_scaled_5_zoib_k_opt_2)

pl_k5_scaled_5_2 <- kfold(bayesian_pl_mean_corrected_scaled_5_zoib_k_opt_2, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
pl_k5_scaled_5_2                  # prints the elpd_kfold etc.
print(pl_k5_scaled_5_2)
summary(pl_k5_scaled_5_2)          # shows estimates and pointwise columns
pl_k5_scaled_5_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(pl_k5_scaled_5_2$pointwise)   # fold-wise contributions per observation
saveRDS(pl_k5_scaled_5_2, "../brms_models/kfold_bayesian_pl_mean_corrected_scaled_5_zoib_k_opt_2.rds")
kfold_bayesian_pl_mean_corrected_scaled_5_zoib_k_opt_2 <- readRDS("../brms_models/kfold_bayesian_pl_mean_corrected_scaled_5_zoib_k_opt_2.rds")
print(kfold_bayesian_pl_mean_corrected_scaled_5_zoib_k_opt_2)



saveRDS(bayesian_pl_mean_corrected_scaled_35_zoib_k_opt_2, file = "../brms_models/bayesian_pl_mean_corrected_scaled_35_zoib_k_opt_2.rds")
bayesian_pl_mean_corrected_scaled_35_zoib_k_opt_2 <- readRDS("../brms_models/bayesian_pl_mean_corrected_scaled_35_zoib_k_opt_2.rds")
summary(bayesian_pl_mean_corrected_scaled_35_zoib_k_opt_2)
loo(bayesian_pl_mean_corrected_scaled_35_zoib_k_opt_2)

pl_k5_scaled_35_2 <- kfold(bayesian_pl_mean_corrected_scaled_35_zoib_k_opt_2, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
pl_k5_scaled_35_2                  # prints the elpd_kfold etc.
print(pl_k5_scaled_35_2)
summary(pl_k5_scaled_35_2)          # shows estimates and pointwise columns
pl_k5_scaled_35_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(pl_k5_scaled_35_2$pointwise)   # fold-wise contributions per observation
saveRDS(pl_k5_scaled_35_2, "../brms_models/kfold_bayesian_pl_mean_corrected_scaled_35_zoib_k_opt_2.rds")
kfold_bayesian_pl_mean_corrected_scaled_35_zoib_k_opt_2 <- readRDS("../brms_models/kfold_bayesian_pl_mean_corrected_scaled_35_zoib_k_opt_2.rds")
print(kfold_bayesian_pl_mean_corrected_scaled_35_zoib_k_opt_2)


#DOING
saveRDS(bayesian_pl_mean_corrected_scaled_15_zoib_k_opt_2, file = "../brms_models/bayesian_pl_mean_corrected_scaled_15_zoib_k_opt_2.rds")
bayesian_pl_mean_corrected_scaled_15_zoib_k_opt_2 <- readRDS("../brms_models/bayesian_pl_mean_corrected_scaled_15_zoib_k_opt_2.rds")
summary(bayesian_pl_mean_corrected_scaled_15_zoib_k_opt_2)
loo(bayesian_pl_mean_corrected_scaled_15_zoib_k_opt_2)

pl_k5_scaled_15_2 <- kfold(bayesian_pl_mean_corrected_scaled_15_zoib_k_opt_2, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
pl_k5_scaled_15_2                  # prints the elpd_kfold etc.
print(pl_k5_scaled_15_2)
summary(pl_k5_scaled_15_2)          # shows estimates and pointwise columns
pl_k5_scaled_15_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(pl_k5_scaled_15_2$pointwise)   # fold-wise contributions per observation
saveRDS(pl_k5_scaled_15_2, "../brms_models/kfold_bayesian_pl_mean_corrected_scaled_15_zoib_k_opt_2.rds")
kfold_bayesian_pl_mean_corrected_scaled_15_zoib_k_opt_2 <- readRDS("../brms_models/kfold_bayesian_pl_mean_corrected_scaled_15_zoib_k_opt_2.rds")
print(kfold_bayesian_pl_mean_corrected_scaled_15_zoib_k_opt_2)



loo_compare(kfold_bayesian_pl_mean_corrected_scaled_13_zoib_k_opt_2,
            kfold_bayesian_pl_mean_corrected_scaled_135_zoib_k_opt_2,
            kfold_bayesian_pl_mean_corrected_scaled_15_zoib_k_opt_2,
            kfold_bayesian_pl_mean_corrected_scaled_35_zoib_k_opt_2,
            kfold_bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2,
            kfold_bayesian_pl_mean_corrected_scaled_1_zoib_k_opt_2,
            kfold_bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2,
            kfold_bayesian_pl_mean_corrected_scaled_5_zoib_k_opt_2)











loo_compare(pl_k5_scaled_13_2, pl_k5_scaled_null_2)



















































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
    sz = pmax(z_func_sd   / sqrt(pmax(n_reps,  1L)), eps),
    
    # --- quadratic terms and their delta-method SEs ---
    x2  = x^2,  sx2 = pmax(2 * abs(x) * sx, eps),
    y2  = y^2,  sy2 = pmax(2 * abs(y) * sy, eps),
    z2  = z^2,  sz2 = pmax(2 * abs(z) * sz, eps)
  )

summary(ao.final.table$sx)   # typically small, e.g. 0.05–0.4 depending on n
summary(ao.final.table$sy)
summary(ao.final.table$sz)

summary(ao.final.table$sx2)  # same order of magnitude as sx, not gigantic
summary(ao.final.table$sy2)
summary(ao.final.table$sz2)

ggplot(ao.final.table,
       aes(x = z_flow_mean,
           y = mean_ao_index)) +
  geom_point()

ggplot(ao.final.table,
       aes(x = z_morpho_mean,
           y = mean_ao_index)) +
  geom_point()





form_all <- bf(
  mean_ao_index | weights(ao_index_weight_12) ~
    #me(x, sx) +
    #me(x2, sx2) +
    me(y, sy) +
    me(y2, sy2) + 
    #1 +
    (1|species)
)


# Suppose the names come back as: mexsx, mex2sx2, meys y..., mezsz, mez2sz2 etc.
pri_all <- c(
  # Fixed effects (logit link for mu):
  prior(normal(0, 1),    class = "b"),                       # all slopes by default
  prior(normal(0, 2),    class = "Intercept"),               # mu-intercept on logit scale
  #prior(normal(0, 0.5), class = "b", coef = "mex2sx2"),
  #prior(normal(0, 0.5), class = "b", coef = "mey2sy2"),
  # Random intercept SD (half-Student-t):
  prior(student_t(3, 0, 2.5), class = "sd", group = "species")                    
)

View(ao.final.table)

bayesian_ao_mean_corrected_scaled_34_zib_k_opt <- brm(
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


saveRDS(bayesian_ao_mean_corrected_scaled_13_zib_k_opt, file = "../brms_models/bayesian_ao_mean_corrected_scaled_13_zib_k_opt.rds")
bayesian_ao_mean_corrected_scaled_13_zib_k_opt <- readRDS("../brms_models/bayesian_ao_mean_corrected_scaled_13_zib_k_opt.rds")
summary(bayesian_ao_mean_corrected_scaled_13_zib_k_opt)
loo(bayesian_ao_mean_corrected_scaled_13_zib_k_opt)

ao_k5_scaled_13_2 <- kfold(bayesian_ao_mean_corrected_scaled_13_zib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
ao_k5_scaled_13_2                  # prints the elpd_kfold etc.
print(ao_k5_scaled_13_2)
summary(ao_k5_scaled_13_2)          # shows estimates and pointwise columns
ao_k5_scaled_13_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(ao_k5_scaled_13_2$pointwise)   # fold-wise contributions per observation
saveRDS(ao_k5_scaled_13_2, "../brms_models/kfold_bayesian_ao_mean_corrected_scaled_13_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_13_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_13_zib_k_opt.rds")
print(kfold_bayesian_ao_mean_corrected_scaled_13_zib_k_opt)


saveRDS(bayesian_ao_mean_corrected_scaled_null_zib_k_opt, file = "../brms_models/bayesian_ao_mean_corrected_scaled_null_zib_k_opt.rds")
bayesian_ao_mean_corrected_scaled_null_zib_k_opt <- readRDS("../brms_models/bayesian_ao_mean_corrected_scaled_null_zib_k_opt.rds")
summary(bayesian_ao_mean_corrected_scaled_null_zib_k_opt)
loo(bayesian_ao_mean_corrected_scaled_null_zib_k_opt)

ao_k5_scaled_null <- kfold(bayesian_ao_mean_corrected_scaled_null_zib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
ao_k5_scaled_null                  # prints the elpd_kfold etc.
print(ao_k5_scaled_null)
summary(ao_k5_scaled_null)          # shows estimates and pointwise columns
ao_k5_scaled_null$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(ao_k5_scaled_null$pointwise)   # fold-wise contributions per observation
saveRDS(ao_k5_scaled_null, "../brms_models/kfold_bayesian_ao_mean_corrected_scaled_null_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_null_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_null_zib_k_opt.rds")
print(kfold_bayesian_ao_mean_corrected_scaled_null_zib_k_opt)


saveRDS(bayesian_ao_mean_corrected_scaled_1_zib_k_opt, file = "../brms_models/bayesian_ao_mean_corrected_scaled_1_zib_k_opt.rds")
bayesian_ao_mean_corrected_scaled_1_zib_k_opt <- readRDS("../brms_models/bayesian_ao_mean_corrected_scaled_1_zib_k_opt.rds")
summary(bayesian_ao_mean_corrected_scaled_1_zib_k_opt)
loo(bayesian_ao_mean_corrected_scaled_1_zib_k_opt)

ao_k5_scaled_1 <- kfold(bayesian_ao_mean_corrected_scaled_1_zib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
ao_k5_scaled_1                  # prints the elpd_kfold etc.
print(ao_k5_scaled_1)
summary(ao_k5_scaled_1)          # shows estimates and pointwise columns
ao_k5_scaled_1$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(ao_k5_scaled_1$pointwise)   # fold-wise contributions per observation
saveRDS(ao_k5_scaled_1, "../brms_models/kfold_bayesian_ao_mean_corrected_scaled_1_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_1_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_1_zib_k_opt.rds")
print(kfold_bayesian_ao_mean_corrected_scaled_1_zib_k_opt)


saveRDS(bayesian_ao_mean_corrected_scaled_1_zib_k_opt_tighter, file = "../brms_models/bayesian_ao_mean_corrected_scaled_1_zib_k_opt_tighter.rds")
bayesian_ao_mean_corrected_scaled_1_zib_k_opt_tighter <- readRDS("../brms_models/bayesian_ao_mean_corrected_scaled_1_zib_k_opt_tighter.rds")
summary(bayesian_ao_mean_corrected_scaled_1_zib_k_opt_tighter)
loo(bayesian_ao_mean_corrected_scaled_1_zib_k_opt_tighter)

ao_k5_scaled_1_tighter <- kfold(bayesian_ao_mean_corrected_scaled_1_zib_k_opt_tighter, K = 5, cores = 5)
ao_k5_scaled_1_tighter                  # prints the elpd_kfold etc.
print(ao_k5_scaled_1_tighter)
summary(ao_k5_scaled_1_tighter)          # shows estimates and pointwise columns
ao_k5_scaled_1_tighter$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(ao_k5_scaled_1_tighter$pointwise)   # fold-wise contributions per observation
saveRDS(ao_k5_scaled_1_tighter, "../brms_models/kfold_bayesian_ao_mean_corrected_scaled_1_zib_k_opt_tighter.rds")
kfold_bayesian_ao_mean_corrected_scaled_1_zib_k_opt_tighter <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_1_zib_k_opt_tighter.rds")
print(kfold_bayesian_ao_mean_corrected_scaled_1_zib_k_opt_tighter)


saveRDS(bayesian_ao_mean_corrected_scaled_3_zib_k_opt, file = "../brms_models/bayesian_ao_mean_corrected_scaled_3_zib_k_opt.rds")
bayesian_ao_mean_corrected_scaled_3_zib_k_opt <- readRDS("../brms_models/bayesian_ao_mean_corrected_scaled_3_zib_k_opt.rds")
summary(bayesian_ao_mean_corrected_scaled_3_zib_k_opt)
loo(bayesian_ao_mean_corrected_scaled_3_zib_k_opt)

ao_k5_scaled_3 <- kfold(bayesian_ao_mean_corrected_scaled_3_zib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
ao_k5_scaled_3                  # prints the elpd_kfold etc.
print(ao_k5_scaled_3)
summary(ao_k5_scaled_3)          # shows estimates and pointwise columns
ao_k5_scaled_3$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(ao_k5_scaled_3$pointwise)   # fold-wise contributions per observation
saveRDS(ao_k5_scaled_3, "../brms_models/kfold_bayesian_ao_mean_corrected_scaled_3_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_3_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_3_zib_k_opt.rds")
print(kfold_bayesian_ao_mean_corrected_scaled_3_zib_k_opt)


saveRDS(bayesian_ao_mean_corrected_scaled_12_zib_k_opt, file = "../brms_models/bayesian_ao_mean_corrected_scaled_12_zib_k_opt.rds")
bayesian_ao_mean_corrected_scaled_12_zib_k_opt <- readRDS("../brms_models/bayesian_ao_mean_corrected_scaled_12_zib_k_opt.rds")
summary(bayesian_ao_mean_corrected_scaled_12_zib_k_opt)
loo(bayesian_ao_mean_corrected_scaled_12_zib_k_opt)

ao_k5_scaled_12 <- kfold(bayesian_ao_mean_corrected_scaled_12_zib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
ao_k5_scaled_12                  # prints the elpd_kfold etc.
print(ao_k5_scaled_12)
summary(ao_k5_scaled_12)          # shows estimates and pointwise columns
ao_k5_scaled_12$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(ao_k5_scaled_12$pointwise)   # fold-wise contributions per observation
saveRDS(ao_k5_scaled_12, "../brms_models/kfold_bayesian_ao_mean_corrected_scaled_12_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_12_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_12_zib_k_opt.rds")
print(kfold_bayesian_ao_mean_corrected_scaled_12_zib_k_opt)


saveRDS(bayesian_ao_mean_corrected_scaled_34_zib_k_opt, file = "../brms_models/bayesian_ao_mean_corrected_scaled_34_zib_k_opt.rds")
bayesian_ao_mean_corrected_scaled_34_zib_k_opt <- readRDS("../brms_models/bayesian_ao_mean_corrected_scaled_34_zib_k_opt.rds")
summary(bayesian_ao_mean_corrected_scaled_34_zib_k_opt)
loo(bayesian_ao_mean_corrected_scaled_34_zib_k_opt)

ao_k5_scaled_34 <- kfold(bayesian_ao_mean_corrected_scaled_34_zib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
ao_k5_scaled_34                  # prints the elpd_kfold etc.
print(ao_k5_scaled_34)
summary(ao_k5_scaled_34)          # shows estimates and pointwise columns
ao_k5_scaled_34$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(ao_k5_scaled_34$pointwise)   # fold-wise contributions per observation
saveRDS(ao_k5_scaled_34, "../brms_models/kfold_bayesian_ao_mean_corrected_scaled_34_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_34_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_34_zib_k_opt.rds")
print(kfold_bayesian_ao_mean_corrected_scaled_34_zib_k_opt)

loo_compare(ao_k5_scaled_13_2, ao_k5_scaled_null, ao_k5_scaled_1, ao_k5_scaled_1_tighter, ao_k5_scaled_12, ao_k5_scaled_34)












ao.final.table <- ao.final.table %>%
  mutate(
    # --- means on the z-scale (already computed) ---
    x = z_flow_mean,   # visited flowers per minute (scaled) -- or use z_vis_mean if that's your x
    y = z_morpho_mean,
    z = z_func_mean,
    
    # --- SEs of the means on the same scale (SD / sqrt(n)) ---
    sx = pmax(z_flow_sd   / sqrt(pmax(n_reps,  1L)), eps),
    sy = pmax(z_morpho_sd  / sqrt(pmax(n_reps, 1L)), eps),
    sz = pmax(z_func_sd   / sqrt(pmax(n_reps,  1L)), eps),
    
    # --- quadratic terms and their delta-method SEs ---
    x2  = x^2,  sx2 = pmax(2 * abs(x) * sx, eps),
    y2  = y^2,  sy2 = pmax(2 * abs(y) * sy, eps),
    z2  = z^2,  sz2 = pmax(2 * abs(z) * sz, eps)
  )

# main model: nonlinear effects via smooths of the latent predictors
bf_main <- bf(mean_ao_index | weights(ao_index_weight_12) ~
                s(z_flow_mean) +
                s(z_morpho_mean) +
                s(z_func_mean) +
                (1|species))


# measurement models: link observed predictors to their latent true values
bf_me1 <- bf(x | mi(sx) ~ 1)
bf_me2 <- bf(y | mi(sy) ~ 1)
bf_me3 <- bf(z | mi(sz) ~ 1)

# combine into a multivariate model
fit <- brm(bf_main + bf_me1 + bf_me2 + bf_me3 + set_rescor(FALSE))

# Suppose the names come back as: mexsx, mex2sx2, meys y..., mezsz, mez2sz2 etc.
pri_all <- c(
  # Fixed effects (logit link for mu):
  #prior(normal(0, 1),    class = "b"),                       # all slopes by default
  prior(normal(0, 2),    class = "Intercept")               # mu-intercept on logit scale
  #prior(normal(0, 0.5), class = "b", coef = "mex2sx2"),
  #prior(normal(0, 0.5), class = "b", coef = "mey2sy2"),
  # Random intercept SD (half-Student-t):
  #prior(student_t(3, 0, 2.5), class = "sd", group = "species")                    
)

View(ao.final.table)

bayesian_ao_complex_1_side <- brm(
  bf_me1, 
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



saveRDS(bayesian_ao_complex_135_main file = "../brms_models/bayesian_ao_complex_135_main.rds")

###


























form_all <- bf(
  mean_ao_index | weights(ao_index_weight_12) ~
    me(x, sx) +
    #me(x2, sx2) +
    me(y, sy) +
    #me(y2, sy2) +
    #me(z, sz) +
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
  #prior(normal(0, 0.5), class = "b", coef = "mex2sx2"),
  #prior(normal(0, 0.5), class = "b", coef = "mey2sy2"),
  # Random intercept SD (half-Student-t):
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

#DOING
saveRDS(bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt, file = "../brms_models/bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt.rds")
bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt <- readRDS("../brms_models/bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt.rds")
summary(bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt)
loo(bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt)

ao_k5_scaled_13_2 <- kfold(bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
ao_k5_scaled_13_2                  # prints the elpd_kfold etc.
print(ao_k5_scaled_13_2)
summary(ao_k5_scaled_13_2)          # shows estimates and pointwise columns
ao_k5_scaled_13_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(ao_k5_scaled_13_2$pointwise)   # fold-wise contributions per observation
saveRDS(ao_k5_scaled_13_2, "../brms_models/kfold_bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt.rds")
print(kfold_bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt)


saveRDS(bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt, file = "../brms_models/bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt.rds")
bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt <- readRDS("../brms_models/bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt.rds")
summary(bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt)
loo(bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt)

ao_k5_scaled_null_2 <- kfold(bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
ao_k5_scaled_null_2                  # prints the elpd_kfold etc.
print(ao_k5_scaled_null_2)
summary(ao_k5_scaled_null_2)          # shows estimates and pointwise columns
ao_k5_scaled_null_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(ao_k5_scaled_null_2$pointwise)   # fold-wise contributions per observation
saveRDS(ao_k5_scaled_null_2, "../brms_models/kfold_bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt.rds")
print(kfold_bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt)


# saveRDS(bayesian_ao_mean_corrected_scaled_1234_2_zib_k_opt, file = "../brms_models/bayesian_ao_mean_corrected_scaled_1234_2_zib_k_opt.rds")
# bayesian_ao_mean_corrected_scaled_1234_2_zib_k_opt <- readRDS("../brms_models/bayesian_ao_mean_corrected_scaled_1234_2_zib_k_opt.rds")
# summary(bayesian_ao_mean_corrected_scaled_1234_2_zib_k_opt)
# loo(bayesian_ao_mean_corrected_scaled_1234_2_zib_k_opt)
# 
# ao_k5_scaled_2 <- kfold(bayesian_ao_mean_corrected_scaled_1234_2_zib_k_opt, K = 5, cores = 5)
# ###* this lasts more than 15-20 minutes
# ao_k5_scaled_2                  # prints the elpd_kfold etc.
# print(ao_k5_scaled_2)
# summary(ao_k5_scaled_2)          # shows estimates and pointwise columns
# ao_k5_scaled_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
# head(ao_k5_scaled_2$pointwise)   # fold-wise contributions per observation
# saveRDS(ao_k5_scaled_2, "../brms_models/kfold_bayesian_ao_mean_corrected_scaled_1234_2_zib_k_opt.rds")
# kfold_bayesian_ao_mean_corrected_scaled_1234_2_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_1234_2_zib_k_opt.rds")
# print(kfold_bayesian_ao_mean_corrected_scaled_1234_2_zib_k_opt)


saveRDS(bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt, file = "../brms_models/bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt.rds")
bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt <- readRDS("../brms_models/bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt.rds")
summary(bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt)
loo(bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt)

ao_k5_scaled_1_2 <- kfold(bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
ao_k5_scaled_1_2                  # prints the elpd_kfold etc.
print(ao_k5_scaled_1_2)
summary(ao_k5_scaled_1_2)          # shows estimates and pointwise columns
ao_k5_scaled_1_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(ao_k5_scaled_1_2$pointwise)   # fold-wise contributions per observation
saveRDS(ao_k5_scaled_1_2, "../brms_models/kfold_bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt.rds")
print(kfold_bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt)


# saveRDS(bayesian_ao_mean_corrected_scaled_12_2_zib_k_opt, file = "../brms_models/bayesian_ao_mean_corrected_scaled_12_2_zib_k_opt.rds")
# bayesian_ao_mean_corrected_scaled_12_2_zib_k_opt <- readRDS("../brms_models/bayesian_ao_mean_corrected_scaled_12_2_zib_k_opt.rds")
# summary(bayesian_ao_mean_corrected_scaled_12_2_zib_k_opt)
# loo(bayesian_ao_mean_corrected_scaled_12_2_zib_k_opt)
# 
# ao_k5_scaled_12_2 <- kfold(bayesian_ao_mean_corrected_scaled_12_2_zib_k_opt, K = 5, cores = 5)
# ###* this lasts more than 15-20 minutes
# ao_k5_scaled_12_2                  # prints the elpd_kfold etc.
# print(ao_k5_scaled_12_2)
# summary(ao_k5_scaled_12_2)          # shows estimates and pointwise columns
# ao_k5_scaled_12_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
# head(ao_k5_scaled_12_2$pointwise)   # fold-wise contributions per observation
# saveRDS(ao_k5_scaled_12_2, "../brms_models/kfold_bayesian_ao_mean_corrected_scaled_12_2_zib_k_opt.rds")
# kfold_bayesian_ao_mean_corrected_scaled_12_2_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_12_2_zib_k_opt.rds")
# print(kfold_bayesian_ao_mean_corrected_scaled_12_2_zib_k_opt)


# saveRDS(bayesian_ao_mean_corrected_scaled_34_2_zib_k_opt, file = "../brms_models/bayesian_ao_mean_corrected_scaled_34_2_zib_k_opt.rds")
# bayesian_ao_mean_corrected_scaled_34_2_zib_k_opt <- readRDS("../brms_models/bayesian_ao_mean_corrected_scaled_34_2_zib_k_opt.rds")
# summary(bayesian_ao_mean_corrected_scaled_34_2_zib_k_opt)
# loo(bayesian_ao_mean_corrected_scaled_34_2_zib_k_opt)
# 
# ao_k5_scaled_34_2 <- kfold(bayesian_ao_mean_corrected_scaled_34_2_zib_k_opt, K = 5, cores = 5)
# ###* this lasts more than 15-20 minutes
# ao_k5_scaled_34_2                  # prints the elpd_kfold etc.
# print(ao_k5_scaled_34_2)
# summary(ao_k5_scaled_34_2)          # shows estimates and pointwise columns
# ao_k5_scaled_34_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
# head(ao_k5_scaled_34_2$pointwise)   # fold-wise contributions per observation
# saveRDS(ao_k5_scaled_34_2, "../brms_models/kfold_bayesian_ao_mean_corrected_scaled_34_2_zib_k_opt.rds")
# kfold_bayesian_ao_mean_corrected_scaled_34_2_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_34_2_zib_k_opt.rds")
# print(kfold_bayesian_ao_mean_corrected_scaled_34_2_zib_k_opt)



saveRDS(bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt, file = "../brms_models/bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt.rds")
bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt <- readRDS("../brms_models/bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt.rds")
summary(bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt)
loo(bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt)

ao_k5_scaled_3_2 <- kfold(bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
ao_k5_scaled_3_2                  # prints the elpd_kfold etc.
print(ao_k5_scaled_3_2)
summary(ao_k5_scaled_3_2)          # shows estimates and pointwise columns
ao_k5_scaled_3_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(ao_k5_scaled_3_2$pointwise)   # fold-wise contributions per observation
saveRDS(ao_k5_scaled_3_2, "../brms_models/kfold_bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt.rds")
print(kfold_bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt)


saveRDS(bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt, file = "../brms_models/bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt.rds")
bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt <- readRDS("../brms_models/bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt.rds")
summary(bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt)
loo(bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt)

ao_k5_scaled_5_2 <- kfold(bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
ao_k5_scaled_5_2                  # prints the elpd_kfold etc.
print(ao_k5_scaled_5_2)
summary(ao_k5_scaled_5_2)          # shows estimates and pointwise columns
ao_k5_scaled_5_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(ao_k5_scaled_5_2$pointwise)   # fold-wise contributions per observation
saveRDS(ao_k5_scaled_5_2, "../brms_models/kfold_bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt.rds")
print(kfold_bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt)



saveRDS(bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt, file = "../brms_models/bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt.rds")
bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt <- readRDS("../brms_models/bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt.rds")
summary(bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt)
loo(bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt)

ao_k5_scaled_15_2 <- kfold(bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
ao_k5_scaled_15_2                  # prints the elpd_kfold etc.
print(ao_k5_scaled_15_2)
summary(ao_k5_scaled_15_2)          # shows estimates and pointwise columns
ao_k5_scaled_15_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(ao_k5_scaled_15_2$pointwise)   # fold-wise contributions per observation
saveRDS(ao_k5_scaled_15_2, "../brms_models/kfold_bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt.rds")
print(kfold_bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt)



saveRDS(bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt, file = "../brms_models/bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt.rds")
bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt <- readRDS("../brms_models/bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt.rds")
summary(bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt)
loo(bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt)

ao_k5_scaled_35_2 <- kfold(bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
ao_k5_scaled_35_2                  # prints the elpd_kfold etc.
print(ao_k5_scaled_35_2)
summary(ao_k5_scaled_35_2)          # shows estimates and pointwise columns
ao_k5_scaled_35_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(ao_k5_scaled_35_2$pointwise)   # fold-wise contributions per observation
saveRDS(ao_k5_scaled_35_2, "../brms_models/kfold_bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt.rds")
print(kfold_bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt)


#DOing
saveRDS(bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt, file = "../brms_models/bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt.rds")
bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt <- readRDS("../brms_models/bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt.rds")
summary(bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt)
loo(bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt)

ao_k5_scaled_135_2 <- kfold(bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
ao_k5_scaled_135_2                  # prints the elpd_kfold etc.
print(ao_k5_scaled_135_2)
summary(ao_k5_scaled_135_2)          # shows estimates and pointwise columns
ao_k5_scaled_135_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(ao_k5_scaled_135_2$pointwise)   # fold-wise contributions per observation
saveRDS(ao_k5_scaled_135_2, "../brms_models/kfold_bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt.rds")
print(kfold_bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt)






loo_compare(kfold_bayesian_ao_mean_corrected_scaled_135_2_zib_k_opt,
            kfold_bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt,
            kfold_bayesian_ao_mean_corrected_scaled_3_2_zib_k_opt,
            kfold_bayesian_ao_mean_corrected_scaled_5_2_zib_k_opt,
            kfold_bayesian_ao_mean_corrected_scaled_13_2_zib_k_opt,
            kfold_bayesian_ao_mean_corrected_scaled_15_2_zib_k_opt,
            kfold_bayesian_ao_mean_corrected_scaled_35_2_zib_k_opt,
            kfold_bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt)















































# ###* GO INDEX
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
    sz = pmax(z_func_sd   / sqrt(pmax(n_reps,  1L)), eps),

    # --- quadratic terms and their delta-method SEs ---
    x2  = x^2,  sx2 = pmax(2 * abs(x) * sx, eps),
    y2  = y^2,  sy2 = pmax(2 * abs(y) * sy, eps),
    z2  = z^2,  sz2 = pmax(2 * abs(z) * sz, eps)
  )

summary(go.final.table$sx)   # typically small, e.g. 0.05–0.4 depending on n
summary(go.final.table$sy)
summary(go.final.table$sz)

summary(go.final.table$sx2)  # same order of magnitude as sx, not gigantic
summary(go.final.table$sy2)
summary(go.final.table$sz2)

ggplot(go.final.table,
       aes(x = z_flow_mean,
           y = mean_go_index)) +
  geom_point()

ggplot(go.final.table,
       aes(x = z_morpho_mean,
           y = mean_go_index)) +
  geom_point()





form_all <- bf(
  mean_go_index | weights(go_index_weight_12) ~
    #me(x, sx) +
    #me(x2, sx2) +
    me(y, sy) +
    #me(y2, sy2) +
    #1 +
    (1|species),
  phi ~ 1,
  zoi ~ 1,
  coi ~ 1
)


# Suppose the names come back as: mexsx, mex2sx2, meys y..., mezsz, mez2sz2 etc.
pri_all <- c(
  # Fixed effects (logit link for mu):
  prior(normal(0, 1),    class = "b"),                       # all slopes by default
  prior(normal(0, 2),    class = "Intercept"),               # mu-intercept on logit scale
  #prior(normal(0, 0.5), class = "b", coef = "mex2sx2"),
  #prior(normal(0, 0.5), class = "b", coef = "mey2sy2"),
  # Random intercept SD (half-Student-t):
  prior(student_t(3, 0, 2.5), class = "sd", group = "species"),
  prior(normal(0, 1), class = "Intercept", dpar = "phi"),
  prior(normal(0, 1), class = "Intercept", dpar = "zoi"),
  prior(normal(0, 1), class = "Intercept", dpar = "coi")
)

View(go.final.table)

bayesian_go_mean_corrected_scaled_3_zib_k_opt <- brm(
  form_all,
  data = go.final.table,
  family = zero_one_inflated_beta(),
  prior = pri_all,
  chains = 5, cores = 5,
  iter = 2000, warmup = 1000,
  threads = threading(4),
  save_pars = save_pars(all = TRUE, latent = TRUE),
  control = list(adapt_delta = 0.95, max_treedepth = 20),
  seed = 1234
)


saveRDS(bayesian_go_mean_corrected_scaled_13_zib_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_13_zib_k_opt.rds")
bayesian_go_mean_corrected_scaled_13_zib_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_13_zib_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_13_zib_k_opt)
loo(bayesian_go_mean_corrected_scaled_13_zib_k_opt)

go_k5_scaled_13 <- kfold(bayesian_go_mean_corrected_scaled_13_zib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
go_k5_scaled_13                  # prints the elpd_kfold etc.
print(go_k5_scaled_13)
summary(go_k5_scaled_13)          # shows estimates and pointwise columns
go_k5_scaled_13$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(go_k5_scaled_13$pointwise)   # fold-wise contributions per observation
saveRDS(go_k5_scaled_13, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_13_zib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_13_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_13_zib_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_13_zib_k_opt)


saveRDS(bayesian_go_mean_corrected_scaled_null_zib_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_null_zib_k_opt.rds")
bayesian_go_mean_corrected_scaled_null_zib_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_null_zib_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_null_zib_k_opt)
loo(bayesian_go_mean_corrected_scaled_null_zib_k_opt)

go_k5_scaled_null <- kfold(bayesian_go_mean_corrected_scaled_null_zib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
go_k5_scaled_null                  # prints the elpd_kfold etc.
print(go_k5_scaled_null)
summary(go_k5_scaled_null)          # shows estimates and pointwise columns
go_k5_scaled_null$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(go_k5_scaled_null$pointwise)   # fold-wise contributions per observation
saveRDS(go_k5_scaled_null, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_null_zib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_null_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_null_zib_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_null_zib_k_opt)


saveRDS(bayesian_go_mean_corrected_scaled_1_zib_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_1_zib_k_opt.rds")
bayesian_go_mean_corrected_scaled_1_zib_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_1_zib_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_1_zib_k_opt)
loo(bayesian_go_mean_corrected_scaled_1_zib_k_opt)

go_k5_scaled_1 <- kfold(bayesian_go_mean_corrected_scaled_1_zib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
go_k5_scaled_1                  # prints the elpd_kfold etc.
print(go_k5_scaled_1)
summary(go_k5_scaled_1)          # shows estimates and pointwise columns
go_k5_scaled_1$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(go_k5_scaled_1$pointwise)   # fold-wise contributions per observation
saveRDS(go_k5_scaled_1, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_1_zib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_1_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_1_zib_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_1_zib_k_opt)


saveRDS(bayesian_go_mean_corrected_scaled_3_zib_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_3_zib_k_opt.rds")
bayesian_go_mean_corrected_scaled_3_zib_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_3_zib_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_3_zib_k_opt)
loo(bayesian_go_mean_corrected_scaled_3_zib_k_opt)

go_k5_scaled_3 <- kfold(bayesian_go_mean_corrected_scaled_3_zib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
go_k5_scaled_3                  # prints the elpd_kfold etc.
print(go_k5_scaled_3)
summary(go_k5_scaled_3)          # shows estimates and pointwise columns
go_k5_scaled_3$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(go_k5_scaled_3$pointwise)   # fold-wise contributions per observation
saveRDS(go_k5_scaled_3, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_3_zib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_3_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_3_zib_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_3_zib_k_opt)


saveRDS(bayesian_go_mean_corrected_scaled_12_zib_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_12_zib_k_opt.rds")
bayesian_go_mean_corrected_scaled_12_zib_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_12_zib_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_12_zib_k_opt)
loo(bayesian_go_mean_corrected_scaled_12_zib_k_opt)

go_k5_scaled_12 <- kfold(bayesian_go_mean_corrected_scaled_12_zib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
go_k5_scaled_12                  # prints the elpd_kfold etc.
print(go_k5_scaled_12)
summary(go_k5_scaled_12)          # shows estimates and pointwise columns
go_k5_scaled_12$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(go_k5_scaled_12$pointwise)   # fold-wise contributions per observation
saveRDS(go_k5_scaled_12, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_12_zib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_12_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_12_zib_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_12_zib_k_opt)


saveRDS(bayesian_go_mean_corrected_scaled_34_zib_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_34_zib_k_opt.rds")
bayesian_go_mean_corrected_scaled_34_zib_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_34_zib_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_34_zib_k_opt)
loo(bayesian_go_mean_corrected_scaled_34_zib_k_opt)

go_k5_scaled_34 <- kfold(bayesian_go_mean_corrected_scaled_34_zib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
go_k5_scaled_34                  # prints the elpd_kfold etc.
print(go_k5_scaled_34)
summary(go_k5_scaled_34)          # shows estimates and pointwise columns
go_k5_scaled_34$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(go_k5_scaled_34$pointwise)   # fold-wise contributions per observation
saveRDS(go_k5_scaled_34, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_34_zib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_34_zib_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_34_zib_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_34_zib_k_opt)

loo_compare(go_k5_scaled_13, go_k5_scaled_null, go_k5_scaled_1, go_k5_scaled_1_tighter, go_k5_scaled_12, go_k5_scaled_34)




















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
    sz = pmax(z_func_sd   / sqrt(pmax(n_reps,  1L)), eps),

    # --- quadratic terms and their delta-method SEs ---
    x2  = x^2,  sx2 = pmax(2 * abs(x) * sx, eps),
    y2  = y^2,  sy2 = pmax(2 * abs(y) * sy, eps),
    z2  = z^2,  sz2 = pmax(2 * abs(z) * sz, eps)
  )

summary(go.final.table$sx)   # typically small, e.g. 0.05–0.4 depending on n
summary(go.final.table$sy)
summary(go.final.table$sz)

summary(go.final.table$sx2)  # same order of magnitude as sx, not gigantic
summary(go.final.table$sy2)
summary(go.final.table$sz2)

ggplot(go.final.table,
       aes(x = z_flow_mean,
           y = mean_go_index)) +
  geom_point()

ggplot(go.final.table,
       aes(x = z_morpho_mean,
           y = mean_go_index)) +
  geom_point()


# # quick boundary check
table(go.final.table$mean_go_index == 0, useNA="ifany")
table(go.final.table$mean_go_index == 1, useNA="ifany")


# # move exact 0/1 slightly into (0,1)
k <- 100  # 50–200 is typical; k=100 makes tiny shifts (0.005)
go.final.table <- go.final.table |>
  dplyr::mutate(
    mean_go_index_beta = (mean_go_index * (k - 1) + 0.5) / k
  )



form_all <- bf(
  mean_go_index | weights(go_index_weight_12) ~
    #me(x, sx) +
    #me(x2, sx2) +
    me(y, sy) +
    me(y2, sy2) +
    #1 +
    (1|species),
  phi ~ 1
)

gow <- go.final.table$go_index_weight_12
go <- go.final.table$mean_go_index

mu_bar_go <- weighted.mean(go[go > 0 & go < 1], gow[go > 0 & go < 1], na.rm = TRUE)
# fall back safely if needed
mu_bar_go <- ifelse(is.finite(mu_bar_go), mu_bar_go, weighted.mean(go, gow, na.rm = TRUE))
logit_go_mu0 <- qlogis(pmin(pmax(mu_bar_go, 1e-6), 1 - 1e-6))
logit_go_mu0
-0.3329464

# Suppose the names come back as: mexsx, mex2sx2, meys y..., mezsz, mez2sz2 etc.
pri_all <- c(
  # Fixed effects (logit link for mu):
  prior(normal(0, 1),    class = "b"),                       # all slopes by default
  prior(normal(-0.5868932, 1.5),    class = "Intercept"),               # mu-intercept on logit scale
  #prior(normal(0, 0.5), class = "b", coef = "mex2sx2"),
  #prior(normal(0, 0.5), class = "b", coef = "mey2sy2"),
  prior(student_t(3, 0, 2.5), class = "sd", group = "species")
)

View(go.final.table)

bayesian_go_mean_corrected_scaled_34_2_zoib_k_opt <- brm(
  form_all,
  data = go.final.table,
  family = zero_one_inflated_beta(),
  prior = pri_all,
  chains = 5, cores = 5,
  iter = 2000, warmup = 1000,
  threads = threading(4),
  save_pars = save_pars(all = TRUE, latent = TRUE),
  control = list(adapt_delta = 0.95, max_treedepth = 20),
  seed = 1234
)


saveRDS(bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt.rds")
bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt)
loo(bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt)

go_k5_scaled_13_2 <- kfold(bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
go_k5_scaled_13_2                  # prints the elpd_kfold etc.
print(go_k5_scaled_13_2)
summary(go_k5_scaled_13_2)          # shows estimates and pointwise columns
go_k5_scaled_13_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(go_k5_scaled_13_2$pointwise)   # fold-wise contributions per observation
saveRDS(go_k5_scaled_13_2, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt)


saveRDS(bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt.rds")
bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt)
loo(bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt)

go_k5_scaled_null_2 <- kfold(bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
go_k5_scaled_null_2                  # prints the elpd_kfold etc.
print(go_k5_scaled_null_2)
summary(go_k5_scaled_null_2)          # shows estimates and pointwise columns
go_k5_scaled_null_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(go_k5_scaled_null_2$pointwise)   # fold-wise contributions per observation
saveRDS(go_k5_scaled_null_2, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt)


saveRDS(bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt.rds")
bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt)
loo(bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt)

go_k5_scaled_1_2 <- kfold(bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
go_k5_scaled_1_2                  # prints the elpd_kfold etc.
print(go_k5_scaled_1_2)
summary(go_k5_scaled_1_2)          # shows estimates and pointwise columns
go_k5_scaled_1_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(go_k5_scaled_1_2$pointwise)   # fold-wise contributions per observation
saveRDS(go_k5_scaled_1_2, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt)


saveRDS(bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt.rds")
bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt)
loo(bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt)

go_k5_scaled_3_2 <- kfold(bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
go_k5_scaled_3_2                  # prints the elpd_kfold etc.
print(go_k5_scaled_3_2)
summary(go_k5_scaled_3_2)          # shows estimates and pointwise columns
go_k5_scaled_3_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(go_k5_scaled_3_2$pointwise)   # fold-wise contributions per observation
saveRDS(go_k5_scaled_3_2, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt)


saveRDS(bayesian_go_mean_corrected_scaled_12_2_beta_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_12_2_zoib_k_opt.rds")
bayesian_go_mean_corrected_scaled_12_2_zoib_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_12_2_zoib_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_12_2_zoib_k_opt)
loo(bayesian_go_mean_corrected_scaled_12_2_zoib_k_opt)

go_k5_scaled_12_2 <- kfold(bayesian_go_mean_corrected_scaled_12_2_zoib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
go_k5_scaled_12_2                  # prints the elpd_kfold etc.
print(go_k5_scaled_12_2)
summary(go_k5_scaled_12_2)          # shows estimates and pointwise columns
go_k5_scaled_12_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(go_k5_scaled_12_2$pointwise)   # fold-wise contributions per observation
saveRDS(go_k5_scaled_12_2, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_12_2_zoib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_12_2_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_12_2_zoib_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_12_2_zoib_k_opt)


saveRDS(bayesian_go_mean_corrected_scaled_34_2_zoib_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_34_2_zoib_k_opt.rds")
bayesian_go_mean_corrected_scaled_34_2_zoib_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_34_2_zoib_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_34_2_zoib_k_opt)
loo(bayesian_go_mean_corrected_scaled_34_2_zoib_k_opt)

go_k5_scaled_34_2 <- kfold(bayesian_go_mean_corrected_scaled_34_2_zoib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
go_k5_scaled_34_2                  # prints the elpd_kfold etc.
print(go_k5_scaled_34_2)
summary(go_k5_scaled_34_2)          # shows estimates and pointwise columns
go_k5_scaled_34_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(go_k5_scaled_34_2$pointwise)   # fold-wise contributions per observation
saveRDS(go_k5_scaled_34_2, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_34_2_zoib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_34_2_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_34_2_zoib_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_34_2_zoib_k_opt)


saveRDS(bayesian_go_mean_corrected_scaled_1234_2_zoib_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_1234_2_zoib_k_opt.rds")
bayesian_go_mean_corrected_scaled_1234_2_zoib_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_1234_2_zoib_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_1234_2_zoib_k_opt)
loo(bayesian_go_mean_corrected_scaled_1234_2_zoib_k_opt)

go_k5_scaled_1234_2 <- kfold(bayesian_go_mean_corrected_scaled_1234_2_zoib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
go_k5_scaled_1234_2                  # prints the elpd_kfold etc.
print(go_k5_scaled_1234_2)
summary(go_k5_scaled_1234_2)          # shows estimates and pointwise columns
go_k5_scaled_1234_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(go_k5_scaled_1234_2$pointwise)   # fold-wise contributions per observation
saveRDS(go_k5_scaled_1234_2, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_1234_2_zoib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_1234_2_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_1234_2_zoib_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_1234_2_zoib_k_opt)




loo_compare(go_k5_scaled_null_2, bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt)

loo_compare(go_k5_scaled_13_2, go_k5_scaled_null_2)

, go_k5_scaled_1, go_k5_scaled_1_tighter, go_k5_scaled_12, go_k5_scaled_34)




go_k5_scaled_13_2 <- kfold(bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt, K = 5, cores = 5)
saveRDS(go_k5_scaled_13_2, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt.rds")
go_k5_scaled_null_2 <- kfold(bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt, K = 5, cores = 5)
saveRDS(go_k5_scaled_null_2, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt.rds")
go_k5_scaled_1_2 <- kfold(bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt, K = 5, cores = 5)
saveRDS(go_k5_scaled_1_2, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt.rds")
go_k5_scaled_3_2 <- kfold(bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt, K = 5, cores = 5)
saveRDS(go_k5_scaled_3_2, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt.rds")
go_k5_scaled_12_2 <- kfold(bayesian_go_mean_corrected_scaled_12_2_zoib_k_opt, K = 5, cores = 5)
saveRDS(go_k5_scaled_12_2, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_12_2_zoib_k_opt.rds")
go_k5_scaled_34_2 <- kfold(bayesian_go_mean_corrected_scaled_34_2_zoib_k_opt, K = 5, cores = 5)
saveRDS(go_k5_scaled_34_2, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_34_2_zoib_k_opt.rds")
go_k5_scaled_1234_2 <- kfold(bayesian_go_mean_corrected_scaled_1234_2_zoib_k_opt, K = 5, cores = 5)
saveRDS(go_k5_scaled_1234_2, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_1234_2_zoib_k_opt.rds")

















fit_go <- bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt  # your already-fitted model

# # Use the same folds (prefer grouped by species if that’s what you used elsewhere)
folds <- loo::kfold_split_grouped(K = 5, x = go.final.table$species)

go_k5_scaled_13_2 <- kfold(
  fit_go,
  K = 5,
  folds = folds,
  chains = 2, cores = 2,         # fewer chains/cores for CV
  iter = 1500, warmup = 750,     # shorter chains work fine for ELPD ranking
  control = list(adapt_delta = 0.95)
)




fit_go_null <- bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt  # your already-fitted model

# # Use the same folds (prefer grouped by species if that’s what you used elsewhere)
folds <- loo::kfold_split_grouped(K = 5, x = go.final.table$species)

go_k5_scaled_null_2 <- kfold(
  fit_go_null,
  K = 5,
  folds = folds,
  chains = 2, cores = 2,         # fewer chains/cores for CV
  iter = 1500, warmup = 750,     # shorter chains work fine for ELPD ranking
  control = list(adapt_delta = 0.95)
)



# pointwise contributions for each fold
pw_full <- attr(go_k5_scaled_13_2, "pointwise")  # matrix: obs × components
pw_null <- attr(go_k5_scaled_null_2, "pointwise")

# quick look at elpd_kfold per obs
full_elpd <- pw_full[, "elpd_kfold"]
null_elpd <- pw_null[, "elpd_kfold"]

which.min(full_elpd)               # index of worst obs (held-out) for the full model
cbind(obs = 1:length(full_elpd),
      full = round(full_elpd, 1),
      null = round(null_elpd, 1))[order(full_elpd), ][1:5, ]






























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
    sz = pmax(z_func_sd   / sqrt(pmax(n_reps,  1L)), eps),
    
    # --- quadratic terms and their delta-method SEs ---
    x2  = x^2,  sx2 = pmax(2 * abs(x) * sx, eps),
    y2  = y^2,  sy2 = pmax(2 * abs(y) * sy, eps),
    z2  = z^2,  sz2 = pmax(2 * abs(z) * sz, eps)
  )

summary(go.final.table$sx)   # typically small, e.g. 0.05–0.4 depending on n
summary(go.final.table$sy)
summary(go.final.table$sz)

summary(go.final.table$sx2)  # same order of magnitude as sx, not gigantic
summary(go.final.table$sy2)
summary(go.final.table$sz2)

ggplot(go.final.table,
       aes(x = z_flow_mean,
           y = mean_go_index)) +
  geom_point()

ggplot(go.final.table,
       aes(x = z_morpho_mean,
           y = mean_go_index)) +
  geom_point()


# quick boundary check
table(go.final.table$mean_go_index == 0, useNA="ifany")
table(go.final.table$mean_go_index == 1, useNA="ifany")


# move exact 0/1 slightly into (0,1)
k <- 100  # 50–200 is typical; k=100 makes tiny shifts (0.005)
go.final.table <- go.final.table |>
  dplyr::mutate(
    mean_go_index_beta = (mean_go_index * (k - 1) + 0.5) / k
  )



form_all <- bf(
  mean_go_index | weights(go_index_weight_12) ~
    me(x, sx) +
    #me(x2, sx2) +
    me(y, sy) +
    #me(y2, sy2) + 
    #1 +
    (1|species),
  phi ~ 1,
  zoi ~1,
  coi ~1
)

gow <- go.final.table$go_index_weight_12
go <- go.final.table$mean_go_index_beta

mu_bar_go <- weighted.mean(go[go > 0 & go < 1], gow[go > 0 & go < 1], na.rm = TRUE)
# fall back safely if needed
mu_bar_go <- ifelse(is.finite(mu_bar_go), mu_bar_go, weighted.mean(go, gow, na.rm = TRUE))
logit_go_mu0 <- qlogis(pmin(pmax(mu_bar_go, 1e-6), 1 - 1e-6))
logit_go_mu0
-0.3329464

# Suppose the names come back as: mexsx, mex2sx2, meys y..., mezsz, mez2sz2 etc.
pri_all <- c(
  # Fixed effects (logit link for mu):
  prior(normal(0, 1),    class = "b"),                       # all slopes by default
  prior(normal(-0.3329464, 1.5),    class = "Intercept"),               # mu-intercept on logit scale
  #prior(normal(0, 0.5), class = "b", coef = "mex2sx2"),
  #prior(normal(0, 0.5), class = "b", coef = "mey2sy2"),
  prior(constant(-10), class="Intercept", dpar="zoi"),    # ≈ 0
  prior(constant(-2.58), class="Intercept", dpar="coi"),  # ≈ 0.07
  prior(student_t(3, 0, 2.5), class = "sd", group = "species")
)

View(go.final.table)

bayesian_go_mean_corrected_scaled_13_2_k_opt <- brm(
  form_all, 
  data = go.final.table, 
  family = zero_one_inflated_beta(),
  prior = pri_all, 
  chains = 5, cores = 5, 
  iter = 2000, warmup = 1000,
  threads = threading(4),
  save_pars = save_pars(all = TRUE, latent = TRUE),
  control = list(adapt_delta = 0.95, max_treedepth = 20),
  seed = 12345
)


loo_compare()


loo_go <- loo(bayesian_go_mean_corrected_scaled_13_2_k_opt,
              moment_match = TRUE,  # fix bad Pareto-k without full refits
              reloo        = TRUE,  # if still bad, refit only those points
              cores        = 4)
print(loo_go)

###Last all the time in the world
saveRDS(bayesian_go_mean_corrected_scaled_12_2_beta_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_12_2_beta_k_opt.rds")
bayesian_go_mean_corrected_scaled_12_2_beta_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_12_2_beta_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_12_2_beta_k_opt)
loo(bayesian_go_mean_corrected_scaled_12_2_beta_k_opt)

go_k5_scaled_12_2 <- kfold(bayesian_go_mean_corrected_scaled_12_2_beta_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
go_k5_scaled_12_2                  # prints the elpd_kfold etc.
print(go_k5_scaled_12_2)
summary(go_k5_scaled_12_2)          # shows estimates and pointwise columns
go_k5_scaled_12_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(go_k5_scaled_12_2$pointwise)   # fold-wise contributions per observation
saveRDS(go_k5_scaled_12_2, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_12_2_beta_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_12_2_beta_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_12_2_beta_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_12_2_beta_k_opt)




saveRDS(bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt.rds")
bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt)
loo(bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt)

go_k5_scaled_13_2 <- kfold(bayesian_go_mean_corrected_scaled_13_2_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
go_k5_scaled_13_2                  # prints the elpd_kfold etc.
print(go_k5_scaled_13_2)
summary(go_k5_scaled_13_2)          # shows estimates and pointwise columns
go_k5_scaled_13_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(go_k5_scaled_13_2$pointwise)   # fold-wise contributions per observation
saveRDS(go_k5_scaled_13_2, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt)


saveRDS(bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt.rds")
bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt)
loo(bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt)

go_k5_scaled_null_2 <- kfold(bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
go_k5_scaled_null_2                  # prints the elpd_kfold etc.
print(go_k5_scaled_null_2)
summary(go_k5_scaled_null_2)          # shows estimates and pointwise columns
go_k5_scaled_null_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(go_k5_scaled_null_2$pointwise)   # fold-wise contributions per observation
saveRDS(go_k5_scaled_null_2, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt)


saveRDS(bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt.rds")
bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt)
loo(bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt)

go_k5_scaled_1_2 <- kfold(bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
go_k5_scaled_1_2                  # prints the elpd_kfold etc.
print(go_k5_scaled_1_2)
summary(go_k5_scaled_1_2)          # shows estimates and pointwise columns
go_k5_scaled_1_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(go_k5_scaled_1_2$pointwise)   # fold-wise contributions per observation
saveRDS(go_k5_scaled_1_2, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_1_2_zoib_k_opt)


saveRDS(bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt.rds")
bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt)
loo(bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt)

go_k5_scaled_3_2 <- kfold(bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
go_k5_scaled_3_2                  # prints the elpd_kfold etc.
print(go_k5_scaled_3_2)
summary(go_k5_scaled_3_2)          # shows estimates and pointwise columns
go_k5_scaled_3_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(go_k5_scaled_3_2$pointwise)   # fold-wise contributions per observation
saveRDS(go_k5_scaled_3_2, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_3_2_zoib_k_opt)




#DOING
saveRDS(bayesian_go_mean_corrected_scaled_34_2_zoib_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_34_2_zoib_k_opt.rds")
bayesian_go_mean_corrected_scaled_34_2_zoib_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_34_2_zoib_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_34_2_zoib_k_opt)
loo(bayesian_go_mean_corrected_scaled_34_2_zoib_k_opt)

go_k5_scaled_34_2 <- kfold(bayesian_go_mean_corrected_scaled_34_2_zoib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
go_k5_scaled_34_2                  # prints the elpd_kfold etc.
print(go_k5_scaled_34_2)
summary(go_k5_scaled_34_2)          # shows estimates and pointwise columns
go_k5_scaled_34_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(go_k5_scaled_34_2$pointwise)   # fold-wise contributions per observation
saveRDS(go_k5_scaled_34_2, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_34_2_zoib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_34_2_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_34_2_zoib_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_34_2_zoib_k_opt)


saveRDS(bayesian_go_mean_corrected_scaled_1234_2_zoib_k_opt, file = "../brms_models/bayesian_go_mean_corrected_scaled_1234_2_zoib_k_opt.rds")
bayesian_go_mean_corrected_scaled_1234_2_zoib_k_opt <- readRDS("../brms_models/bayesian_go_mean_corrected_scaled_1234_2_zoib_k_opt.rds")
summary(bayesian_go_mean_corrected_scaled_1234_2_zoib_k_opt)
loo(bayesian_go_mean_corrected_scaled_1234_2_zoib_k_opt)

go_k5_scaled_1234_2 <- kfold(bayesian_go_mean_corrected_scaled_1234_2_zoib_k_opt, K = 5, cores = 5)
###* this lasts more than 15-20 minutes
go_k5_scaled_1234_2                  # prints the elpd_kfold etc.
print(go_k5_scaled_1234_2)
summary(go_k5_scaled_1234_2)          # shows estimates and pointwise columns
go_k5_scaled_1234_2$estimates         # data frame with elpd_kfold, se, p_kfold, kfoldic
head(go_k5_scaled_1234_2$pointwise)   # fold-wise contributions per observation
saveRDS(go_k5_scaled_1234_2, "../brms_models/kfold_bayesian_go_mean_corrected_scaled_1234_2_zoib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_1234_2_zoib_k_opt <- readRDS("../brms_models/kfold_bayesian_go_mean_corrected_scaled_1234_2_zoib_k_opt.rds")
print(kfold_bayesian_go_mean_corrected_scaled_1234_2_zoib_k_opt)




loo_compare(go_k5_scaled_null_2, bayesian_go_mean_corrected_scaled_13_2_zoib_k_opt)

loo_compare(go_k5_scaled_13_2, go_k5_scaled_null_2)

, go_k5_scaled_1, go_k5_scaled_1_tighter, go_k5_scaled_12, go_k5_scaled_34)






