###* Bayesian scaled June
###* 
###* 
###* 
###* 
###* This script serves as running and interpreting the results from the best 
###* models created with the data which we have.
source("scripts/Bayesian setup June.R")

View(ao.final.table)
###* SEEDSET INDEX 
###* 
###* even though we are using Bayesian statistics, our dataset is
###* very small and so not we cannot get meaningfull results from it. This is 
###* because Hypericum produces hundreds of seeds, wheras the other plants 
###* produce significantly less

View(c.pl.final.table.4)
hist(log(c.pl.final.table.4$mean_seedset_round))

formula <- bf(mean_seedset_round | weights (seedset_weight_12) ~ 
                #1+
                me(mean.visited.flowers.scaled, sd.visited.flowers.scaled) +
                me(mean.morpho.scaled, sd.morpho.scaled) +
                me(mean.func.scaled, sd.func.scaled) + 
                (1|species))
###* Setting priors
prior <- c(
  set_prior("normal(0, 1)", class = "b"), 
  set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))

bayesian_seedset_noel_w12_p1_vmf_h4 <- brm(
  formula = formula,
  data = c.pl.final.table.4,
  chains = 10,
  cores = 10,
  iter = 12000,
  warmup = 5000,
  family = poisson(),
  prior = prior,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234 # Addition to code
)

###* Trying version h4 version for seedset, which has info about Hypericum from
###* elevation 4000 / no produces seeds. Only usable in seedset though

saveRDS(bayesian_seedset_noel_w12_p1_v_h4, file = "../brms_models/bayesian_seedset_noel_w12_p1_v_h4.rds")
saveRDS(bayesian_seedset_noel_w12_p1_m_h4, file = "../brms_models/bayesian_seedset_noel_w12_p1_v_h4.rds")
saveRDS(bayesian_seedset_noel_w12_p1_f_h4, file = "../brms_models/bayesian_seedset_noel_w12_p1_v_h4.rds")
saveRDS(bayesian_seedset_noel_w12_p1_null_h4, file = "../brms_models/bayesian_seedset_noel_w12_p1_null_h4.rds")
saveRDS(bayesian_seedset_noel_w12_p1_vmf_h4, file = "../brms_models/bayesian_seedset_noel_w12_p1_vmf_h4.rds")



#not finished
saveRDS(bayesian_seedset_noel_w12_p1_v_h4, file = "../brms_models/bayesian_seedset_noel_w12_p1_v_h4.rds")



bayesian_seedset_w12_p1_null_h4 <- readRDS("../brms_models/bayesian_seedset_w12_p1_null_h4.rds")
bayesian_seedset_w12_p0.5_h4 <- readRDS("../brms_models/bayesian_seedset_w12_p0.5_h4.rds")
bayesian_seedset_w12_p1_h4 <- readRDS("../brms_models/bayesian_seedset_w12_p1_h4.rds")
bayesian_seedset_w12_p2.5_h4 <- readRDS("../brms_models/bayesian_seedset_w12_p2.5_h4.rds")
bayesian_seedset_w12_pnone_h4 <- readRDS("../brms_models/bayesian_seedset_w12_pnone_h4.rds")


fit <- bayesian_seedset_w12_p0.5_h4
fit <- bayesian_seedset_w12_p1_h4
fit <- bayesian_seedset_w12_p2.5_h4
fit <- bayesian_seedset_w12_pnone_h4

summary(fit)
pp_check(fit, type = "dens_overlay")
pp_check(fit, type = "stat")
bayes_R2(fit)
plot(fit)
fitted_values <- fitted(fit)
residuals <- residuals(fit)
plot(fitted_values, residuals)
abline(h = 0, col = "red")

# Extract predictor range from model data
x_range <- insight::get_data(fit) %>%
  pull(mean.visited.flowers.scaled)
# Generate marginal predictions across the predictor range
data_pred_indiv <- predictions(
  fit,
  newdata = datagrid(
    mean.visited.flowers.scaled = seq(from = min(x_range), to = max(x_range), length.out = 100)
  ))
# Extract the original observed data
observed_data <- insight::get_data(fit)
# Create ggplot
ggplot(data = as.data.frame(data_pred_indiv), 
       aes(x = mean.visited.flowers.scaled, y = estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "steelblue") + # Uncertainty ribbon
  geom_line(color = "steelblue", size = 1) + # Fitted line
  geom_point(data = observed_data, # Observed data points
             aes(x = mean.visited.flowers.scaled, y = mean_seedset_round), 
             color = "black", alpha = 0.6, size = 2) +
  labs(
    title = "Predicted seedset vs. Visitation Rate",
    x = "Mean visitors per minute (scaled)",
    y = "Predicted seedset (mean_seedset_round)"
  ) +
  theme_minimal()

###* As mentioned, we do not see any clear pattern, with the seedset
###* 


formula <- bf(mean_PL_index | weights (PL_index_weight_12) ~ 
                me(mean.visited.flowers.scaled, sd.visited.flowers.scaled) + 
                me(mean.morpho, sd.morpho) +
                (1|species))
priors <- get_prior(
  formula = formula,
  data = c.pl.final.table.5,
  family = zero_one_inflated_beta()
)

# View the priors
print(priors)
###* 
###* POLLEN LIMITATION INDEX
c.pl.final.table.5 <- c.pl.final.table.4 %>%
  filter(species != "Hypericum r" | elevation != 4000)

View(c.pl.final.table.5)

bayesian_pl_noel_w12_p1_vmf <- brm(
  formula = bf(
    mean_PL_index | weights(PL_index_weight_12) ~ 
      #1 +
      me(mean.visited.flowers.scaled, sd.visited.flowers.scaled) +
      me(mean.morpho.scaled, sd.morpho.scaled) +
      me(mean.func.scaled, sd.func.scaled) + 
      (1|species)
  ),
  data = c.pl.final.table.5,
  family = zero_one_inflated_beta(),
  prior = c(
    set_prior("normal(0, 1)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 12000,
  warmup = 5000,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234
)


# saveRDS(bayesian_pl_w12_pnone_12k, file = "../brms_models/bayesian_pl_w12_pnone_12k.rds")
# saveRDS(bayesian_pl_w12_p10_12k, file = "../brms_models/bayesian_pl_w12_p10_12k.rds")
# saveRDS(bayesian_pl_w12_p5_12k, file = "../brms_models/bayesian_pl_w12_p5_12k.rds")
# saveRDS(bayesian_pl_w12_p2.5_12k, file = "../brms_models/bayesian_pl_w12_p2.5_12k.rds")

#saved
saveRDS(bayesian_pl_noel_w12_p1_v_flat, file = "../brms_models/bayesian_pl_noel_w12_p1_v_flat.rds")
saveRDS(bayesian_pl_noel_w12_p1_null_flat, file = "../brms_models/bayesian_pl_noel_w12_p1_null_flat.rds")
saveRDS(bayesian_pl_noel_w12_p1_v, file = "../brms_models/bayesian_pl_noel_w12_p1_v.rds")
saveRDS(bayesian_pl_noel_w12_p1_vm, file = "../brms_models/bayesian_pl_noel_w12_p1_vm.rds")
saveRDS(bayesian_pl_noel_w12_p1_vmf, file = "../brms_models/bayesian_pl_noel_w12_p1_vmf.rds")
saveRDS(bayesian_pl_noel_w12_p1_null, file = "../brms_models/bayesian_pl_noel_w12_p1_null.rds")
saveRDS(bayesian_pl_noel_w12_p1_mf, file = "../brms_models/bayesian_pl_noel_w12_p1_mf.rds")
saveRDS(bayesian_pl_noel_w12_p0.5_vmf, file = "../brms_models/bayesian_pl_noel_w12_p0.5_vmf.rds")
saveRDS(bayesian_pl_noel_w12_p1_m, file = "../brms_models/bayesian_pl_noel_w12_p1_m.rds")
saveRDS(bayesian_pl_noel_w12_p1_f, file = "../brms_models/bayesian_pl_noel_w12_p1_f.rds")


#not finised

    #neXt one with priors for species = "normal (0, 1)"
saveRDS(bayesian_pl_noel_w12_p1_v_2, file = "../brms_models/bayesian_pl_noel_w12_p1_v_2.rds")
saveRDS(bayesian_pl_noel_w12_p1_vm_flat, file = "../brms_models/bayesian_pl_noel_w12_p1_vm_flat.rds")
saveRDS(bayesian_pl_noel_w12_p1_vf, file = "../brms_models/bayesian_pl_noel_w12_p1_vf.rds")
    # THE FOLLOWING DOES NOT CONVERGE PROPERLY



# saveRDS(bayesian_pl_noel_w12_p1_null, file = "../brms_models/bayesian_pl_noel_w12_p1_null.rds")
# saveRDS(bayesian_pl_w12_p0.75_12k, file = "../brms_models/bayesian_pl_w12_p0.75_12k.rds")
# saveRDS(bayesian_pl_w12_p0.5, file = "../brms_models/bayesian_pl_w12_p0.5.rds")
# bayesian_pl_w12_pnone_12k <- readRDS("../brms_models/bayesian_pl_w12_pnone_12k.rds")
# bayesian_pl_w12_p10_12k <- readRDS("../brms_models/bayesian_pl_w12_p10_12k.rds")
# bayesian_pl_w12_p5_12k <- readRDS("../brms_models/bayesian_pl_w12_p5_12k.rds")
# bayesian_pl_w12_p2.5_12k <- readRDS("../brms_models/bayesian_pl_w12_p2.5_12k.rds")
bayesian_pl_noel_w12_p1_null <- readRDS("../brms_models/bayesian_pl_noel_w12_p1_null.rds")
bayesian_pl_noel_w12_p1_v <- readRDS("../brms_models/bayesian_pl_noel_w12_p1_v.rds")
bayesian_pl_noel_w12_p1_m <- readRDS("../brms_models/bayesian_pl_noel_w12_p1_m.rds")
bayesian_pl_noel_w12_p1_f <- readRDS("../brms_models/bayesian_pl_noel_w12_p1_f.rds")
bayesian_pl_noel_w12_p1_vmf <- readRDS("../brms_models/bayesian_pl_noel_w12_p1_vmf.rds")

bayesian_pl_noel_w12_p1_mf <- readRDS("../brms_models/bayesian_pl_noel_w12_p1_mf.rds")

bayesian_pl_noel_w12_p1_vm <- readRDS("../brms_models/bayesian_pl_noel_w12_p1_vm.rds")
bayesian_pl_noel_w12_p1_vf <- readRDS("../brms_models/bayesian_pl_noel_w12_p1_vf.rds")

bayesian_pl_noel_w12_p1_m_flat <- readRDS("../brms_models/bayesian_pl_noel_w12_p1_m_flat.rds")
bayesian_pl_noel_w12_p1_vm_flat <- readRDS("../brms_models/bayesian_pl_noel_w12_p1_vm_flat.rds")
bayesian_pl_noel_w12_p1_v_2 <- readRDS("../brms_models/bayesian_pl_noel_w12_p1_v_2.rds")

###* Cannot compare with loo
###* We can compare bayes_R2
bayes_R2(bayesian_pl_noel_w12_p1_null)
bayes_R2(bayesian_pl_noel_w12_p1_v)
bayes_R2(bayesian_pl_noel_w12_p1_m)
bayes_R2(bayesian_pl_noel_w12_p1_f)
bayes_R2(bayesian_pl_noel_w12_p1_vmf)



bayes_R2(bayesian_pl_noel_w12_p1_v_2)
bayes_R2(bayesian_pl_noel_w12_p1_vm)
bayes_R2(bayesian_pl_noel_w12_p1_vf)
bayes_R2(bayesian_pl_noel_w12_p1_mf)
bayes_R2(bayesian_pl_noel_w12_p0.5_vmf)

bayes_R2(bayesian_pl_noel_w12_p1_m_flat)
bayes_R2(bayesian_pl_noel_w12_p1_vm_flat)
bayes_R2(bayesian_pl_noel_w12_p1_v_flat)
bayes_R2(bayesian_pl_noel_w12_p1_null_flat)


pp_check(bayesian_pl_noel_w12_p1_null, type = "dens_overlay", nsamples = 100)
pp_check(bayesian_pl_noel_w12_p1_v, type = "dens_overlay", nsamples = 100)
pp_check(bayesian_pl_noel_w12_p1_v_2, type = "dens_overlay", nsamples = 100)
pp_check(bayesian_pl_noel_w12_p1_vm, type = "dens_overlay", nsamples = 100)
pp_check(bayesian_pl_noel_w12_p1_vmf, type = "dens_overlay", nsamples = 100)
pp_check(bayesian_pl_noel_w12_p1_m, type = "dens_overlay", nsamples = 100)
pp_check(bayesian_pl_noel_w12_p1_vf, type = "dens_overlay", nsamples = 100)
pp_check(bayesian_pl_noel_w12_p1_f, type = "dens_overlay", nsamples = 100)
pp_check(bayesian_pl_noel_w12_p1_mf, type = "dens_overlay", nsamples = 100)
pp_check(bayesian_pl_noel_w12_p1_m_flat, type = "dens_overlay", nsamples = 100)
pp_check(bayesian_pl_noel_w12_p1_vm_flat, type = "dens_overlay", nsamples = 100)

pp_check(bayesian_pl_noel_w12_p1_null, type = "stat")
pp_check(bayesian_pl_noel_w12_p1_v, type = "stat")
pp_check(bayesian_pl_noel_w12_p1_v_2, type = "stat")
pp_check(bayesian_pl_noel_w12_p1_vm, type = "stat")
pp_check(bayesian_pl_noel_w12_p1_vmf, type = "stat")
pp_check(bayesian_pl_noel_w12_p1_m, type = "stat")
pp_check(bayesian_pl_noel_w12_p1_vf, type = "stat")
pp_check(bayesian_pl_noel_w12_p1_f, type = "stat")
pp_check(bayesian_pl_noel_w12_p1_mf, type = "stat")
pp_check(bayesian_pl_noel_w12_p1_m_flat, type = "stat")
pp_check(bayesian_pl_noel_w12_p1_vm_flat, type = "stat")
library(dplyr)

pp_check(bayesian_pl_noel_w12_p1_vm, , type = "scatter_avg")


# Get unique species and assign them to K folds
set.seed(1234)  # For reproducibility
K <- 5
species_folds <- tibble(
  species = unique(c.pl.final.table.5$species),
  fold = sample(rep(1:K, length.out = n_distinct(c.pl.final.table.5$species)))
)

# Join back to the full dataset to assign folds
folds <- c.pl.final.table.5 %>%
  left_join(species_folds, by = "species") %>%
  pull(fold)

options(future.globals.maxSize = 2 * 1024^3)  # 2 GB

options(mc.cores = 10)
kfold_v <- kfold(bayesian_pl_noel_w12_p1_v_2, folds = folds)
kfold_null <- kfold(bayesian_pl_noel_w12_p1_null, folds = folds)


# fit <- bayesian_pl_w12_pnone_12k
# fit <- bayesian_pl_w12_p10_12k
# fit <- bayesian_pl_w12_p5_12k
# fit <- bayesian_pl_w12_p2.5_12k
fit <- bayesian_pl_noel_w12_p1_null
fit <- bayesian_pl_noel_w12_p1_v
fit <- bayesian_pl_noel_w12_p1_vm
fit <- bayesian_pl_noel_w12_p1_vmf
fit <- bayesian_pl_noel_w12_p1_m
fit <- bayesian_pl_noel_w12_p1_m_flat
fit <- bayesian_pl_noel_w12_p1_vm_flat
fit <- bayesian_pl_noel_w12_p1_v_flat
fit <- bayesian_pl_noel_w12_p0.5_vmf
# fit <- bayesian_pl_w12_p0.75_12k
# fit <- bayesian_pl_w12_p0.5




library(loo)
waic_result <- waic(bayesian_pl_noel_w12_p1_vmf)
print(waic_result)

library(marginaleffects)


summary(fit)
pp_check(fit, type = "dens_overlay", nsamples = 100)
pp_check(fit, type = "stat")
bayes_R2(fit)
plot(fit)
fitted_values <- fitted(fit)
residuals <- residuals(fit)
plot(fitted_values, residuals)
abline(h = 0, col = "red")

# Extract predictor range from model data
x_range <- insight::get_data(fit) %>%
  pull(mean.visited.flowers.scaled)
# Generate marginal predictions across the predictor range
data_pred_indiv <- predictions(
  fit,
  newdata = datagrid(
    mean.visited.flowers.scaled = seq(from = min(x_range), to = max(x_range), length.out = 100)
  ))
# Extract the original observed data
observed_data <- insight::get_data(fit)
# Create ggplot
ggplot(data = as.data.frame(data_pred_indiv), 
       aes(x = mean.visited.flowers.scaled, y = estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "steelblue") + # Uncertainty ribbon
  geom_line(color = "steelblue", size = 1) + # Fitted line
  geom_point(data = observed_data, # Observed data points
             aes(x = mean.visited.flowers.scaled, y = mean_PL_index), 
             color = "black", alpha = 0.6, size = 2) +
  labs(
    title = "Predicted PL index vs. Visitation Rate",
    x = "Mean visitors per minute (scaled)",
    y = "Predicted PL index (mean_PL_index)"
  ) +
  theme_minimal()


# Extract predictor range from model data
x_range <- insight::get_data(fit) %>%
  pull(mean.morpho.scaled)
# Generate marginal predictions across the predictor range
data_pred_indiv <- predictions(
  fit,
  newdata = datagrid(
    mean.morpho.scaled = seq(from = min(x_range), to = max(x_range), length.out = 100)
  ))
# Extract the original observed data
observed_data <- insight::get_data(fit)
# Create ggplot
ggplot(data = as.data.frame(data_pred_indiv), 
       aes(x = mean.morpho.scaled, y = estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "steelblue") + # Uncertainty ribbon
  geom_line(color = "steelblue", size = 1) + # Fitted line
  geom_point(data = observed_data, # Observed data points
             aes(x = mean.morpho.scaled, y = mean_PL_index), 
             color = "black", alpha = 0.6, size = 2) +
  labs(
    title = "Predicted PL index vs. Mean morpho",
    x = "Mean morpho (scaled)",
    y = "Predicted PL index (mean_PL_index)"
  ) +
  theme_minimal()















View(c.pl.final.table.5)







View(ao.final.table)

formula <- bf(mean_ao_index | weights(ao_index_weight_12) ~ me(mean.visited.flowers.scaled, sd.visited.flowers.scaled) + (1|species) + (1|elevation))

priors <- get_prior(
  formula = formula,
  data = ao.final.table,
  family = zero_inflated_beta()
)

# View the priors
print(priors)

formula <- bf(mean_ao_index | weights(ao_index_weight_12) ~ 
                #1 +
                #me(mean.visited.flowers.scaled, sd.visited.flowers.scaled) + 
                me(mean.morpho.scaled, sd.morpho.scaled) +
                me(mean.func.scaled, sd.func.scaled) +
                (1|species))

bayesian_ao_noel_w12_p1_mf <- brm(
  formula = formula,
  data = ao.final.table,
  family = zero_inflated_beta(),
  prior = c(
    set_prior("normal(0, 1)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 12000,
  warmup = 5000,
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234
)

saveRDS(bayesian_ao_noel_w12_p1_v, file = "../brms_models/bayesian_ao_noel_w12_p1_v.rds")
saveRDS(bayesian_ao_noel_w12_p1_null, file = "../brms_models/bayesian_ao_noel_w12_p1_null.rds")
saveRDS(bayesian_ao_noel_w12_p1_m, file = "../brms_models/bayesian_ao_noel_w12_p1_m.rds")
saveRDS(bayesian_ao_noel_w12_p1_f, file = "../brms_models/bayesian_ao_noel_w12_p1_f.rds")
saveRDS(bayesian_ao_noel_w12_p1_vmf, file = "../brms_models/bayesian_ao_noel_w12_p1_vmf.rds")
saveRDS(bayesian_ao_noel_w12_p1_mf, file = "../brms_models/bayesian_ao_noel_w12_p1_mf.rds")

#not finished




bayesian_ao_noel_w12_p1_v <- readRDS("../brms_models/bayesian_ao_noel_w12_p1_v.rds")
bayesian_ao_noel_w12_p1_mf <- readRDS("../brms_models/bayesian_ao_noel_w12_p1_mf.rds")
bayesian_ao_noel_w12_p1_vmf <- readRDS("../brms_models/bayesian_ao_noel_w12_p1_vmf.rds")
bayesian_ao_noel_w12_p1_null <- readRDS("../brms_models/bayesian_ao_noel_w12_p1_null.rds")
bayesian_ao_noel_w12_p1_m <- readRDS("../brms_models/bayesian_ao_noel_w12_p1_m.rds")
bayesian_ao_noel_w12_p1_f <- readRDS("../brms_models/bayesian_ao_noel_w12_p1_f.rds")

bayes_R2(bayesian_ao_noel_w12_p1_v)
bayes_R2(bayesian_ao_noel_w12_p1_mf)
bayes_R2(bayesian_ao_noel_w12_p1_null)
bayes_R2(bayesian_ao_noel_w12_p1_vmf)
bayes_R2(bayesian_ao_noel_w12_p1_m)
bayes_R2(bayesian_ao_noel_w12_p1_f)


###* Somehow it seems that with rising functional group, there is a rise in AO index.. 
###* wtf. Not with morpho

fit <- bayesian_ao_noel_w12_p1_v
fit <- bayesian_ao_noel_w12_p1_mf
fit <- bayesian_ao_noel_w12_p1_vmf
fit <- bayesian_ao_noel_w12_p1_null
fit <- bayesian_ao_noel_w12_p1_m
fit <- bayesian_ao_noel_w12_p1_f

summary(fit)
pp_check(fit, type = "dens_overlay", nsamples = 100)
pp_check(fit, type = "stat")
bayes_R2(fit)
plot(fit)
fitted_values <- fitted(fit)
residuals <- residuals(fit)
plot(fitted_values, residuals)
abline(h = 0, col = "red")


# Extract predictor range from model data
x_range <- insight::get_data(fit) %>%
  pull(mean.func.scaled)
# Generate marginal predictions across the predictor range
data_pred_indiv <- predictions(
  fit,
  newdata = datagrid(
    mean.func.scaled = seq(from = min(x_range), to = max(x_range), length.out = 100)
  ))
# Extract the original observed data
observed_data <- insight::get_data(fit)
# Create ggplot
ggplot(data = as.data.frame(data_pred_indiv), 
       aes(x = mean.func.scaled, y = estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "steelblue") + # Uncertainty ribbon
  geom_line(color = "steelblue", size = 1) + # Fitted line
  geom_point(data = observed_data, # Observed data points
             aes(x = mean.func.scaled, y = mean_ao_index), 
             color = "black", alpha = 0.6, size = 2) +
  labs(
    title = "Predicted AO index vs. func",
    x = "Number of func (scaled)",
    y = "Predicted AO index (mean_AO_index)"
  ) +
  theme_minimal()




















###* 
View(go.final.table)

formula <- bf(mean_go_index | weights(go_index_weight_12) ~ me(mean.visited.flowers.scaled, sd.visited.flowers.scaled) + (1|species) + (1|elevation))

priors <- get_prior(
  formula = formula,
  data = go.final.table,
  family = zero_inflated_beta()
)

# View the priors
print(priors)

formula <- bf(mean_go_index | weights(go_index_weight_12) ~ 
                1 +
                #me(mean.visited.flowers.scaled, sd.visited.flowers.scaled) + 
                #me(mean.morpho.scaled, sd.morpho.scaled) +
                #me(mean.func.scaled, sd.func.scaled) +
                (1|species))

bayesian_go_noel_w12_p1_null <- brm(
  formula = formula,
  data = go.final.table,
  family = zero_one_inflated_beta(),
  prior = c(
    #set_prior("normal(0, 1)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 12000,
  warmup = 5000,
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234
)

saveRDS(bayesian_go_noel_w12_p1_v, file = "../brms_models/bayesian_go_noel_w12_p1_v.rds")
saveRDS(bayesian_go_noel_w12_p1_m, file = "../brms_models/bayesian_go_noel_w12_p1_m.rds")
saveRDS(bayesian_go_noel_w12_p1_f, file = "../brms_models/bayesian_go_noel_w12_p1_f.rds")
saveRDS(bayesian_go_noel_w12_p1_vmf, file = "../brms_models/bayesian_go_noel_w12_p1_vmf.rds")
saveRDS(bayesian_go_noel_w12_p1_null, file = "../brms_models/bayesian_go_noel_w12_p1_null.rds")


#not finished
saveRDS(bayesian_go_noel_w12_p1_mf, file = "../brms_models/bayesian_go_noel_w12_p1_mf.rds")


bayesian_go_noel_w12_p1_v <- readRDS("../brms_models/bayesian_go_noel_w12_p1_v.rds")
bayesian_go_noel_w12_p1_mf <- readRDS("../brms_models/bayesian_go_noel_w12_p1_mf.rds")
bayesian_go_noel_w12_p1_vmf <- readRDS("../brms_models/bayesian_go_noel_w12_p1_vmf.rds")
bayesian_go_noel_w12_p1_null <- readRDS("../brms_models/bayesian_go_noel_w12_p1_null.rds")
bayesian_go_noel_w12_p1_m <- readRDS("../brms_models/bayesian_go_noel_w12_p1_m.rds")
bayesian_go_noel_w12_p1_f <- readRDS("../brms_models/bayesian_go_noel_w12_p1_f.rds")

bayes_R2(bayesian_go_noel_w12_p1_v)
bayes_R2(bayesian_go_noel_w12_p1_mf)
bayes_R2(bayesian_go_noel_w12_p1_null)
bayes_R2(bayesian_go_noel_w12_p1_vmf)
bayes_R2(bayesian_go_noel_w12_p1_m)
bayes_R2(bayesian_go_noel_w12_p1_f)


fit <- bayesian_go_noel_w12_p1_v
fit <- bayesian_go_noel_w12_p1_mf
fit <- bayesian_go_noel_w12_p1_vmf
fit <- bayesian_go_noel_w12_p1_null
fit <- bayesian_go_noel_w12_p1_m
fit <- bayesian_go_noel_w12_p1_f


summary(fit)
pp_check(fit, type = "dens_overlay", nsamples = 100)
pp_check(fit, type = "stat")
bayes_R2(fit)
plot(fit)
fitted_values <- fitted(fit)
residuals <- residuals(fit)
plot(fitted_values, residuals)
abline(h = 0, col = "red")


# Extract predictor range from model data
x_range <- insight::get_data(fit) %>%
  pull(mean.visited.flowers.scaled)
# Generate marginal predictions across the predictor range
data_pred_indiv <- predictions(
  fit,
  newdata = datagrid(
    mean.visited.flowers.scaled = seq(from = min(x_range), to = max(x_range), length.out = 100)
  ))
# Extract the original observed data
observed_data <- insight::get_data(fit)
# Create ggplot
ggplot(data = as.data.frame(data_pred_indiv), 
       aes(x = mean.visited.flowers.scaled, y = estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "steelblue") + # Uncertainty ribbon
  geom_line(color = "steelblue", size = 1) + # Fitted line
  geom_point(data = observed_data, # Observed data points
             aes(x = mean.visited.flowers.scaled, y = mean_go_index), 
             color = "black", alpha = 0.6, size = 2) +
  labs(
    title = "Predicted GO index vs. Visitation Rate",
    x = "Mean visitors per minute (scaled)",
    y = "Predicted GO index (mean_GO_index)"
  ) +
  theme_minimal()
