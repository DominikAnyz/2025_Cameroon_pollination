source("scripts/Bayesian setup June no Hyp.R")

View(c.pl.final.table.4)

formula <- bf(mean_seedset_round | weights (seedset_weight_12) ~ me(mean.visited.flowers.scaled, sd.visited.flowers.scaled) + (1|species) + (1|elevation))
###* Setting priors
prior <- c(
  #set_prior("normal(0, 1)", class = "b"), 
  set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"), 
  set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))

fit_seedset_w12_pnone_nh <- brm(
  formula = formula,
  data = c.pl.final.table.4,
  chains = 10,
  cores = 10,
  iter = 12000,
  warmup = 5000,
  family = poisson(),
  prior = prior,
  save_pars = save_pars(all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 123 # Addition to code
)

saveRDS(fit_seedset_w12_p0.5_nh, file = "../brms_models/fit_seedset_w12_p0.5_nh.rds")
saveRDS(fit_seedset_w12_p1_nh, file = "../brms_models/fit_seedset_w12_p1_nh.rds")
saveRDS(fit_seedset_w12_pnone_nh, file = "../brms_models/fit_seedset_w12_pnone_nh.rds")
fit_seedset_w12_p0.5_nh <- readRDS("../brms_models/fit_seedset_w12_p0.5_nh.rds")
fit_seedset_w12_p1_nh <- readRDS("../brms_models/fit_seedset_w12_p1_nh.rds")
fit_seedset_w12_pnone_nh <- readRDS("../brms_models/fit_seedset_w12_pnone_nh.rds")

fit <- fit_seedset_w12_p0.5_nh
fit <- fit_seedset_w12_p1_nh
fit <- fit_seedset_w12_pnone_nh

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
             aes(x = mean.visited.flowers.scaled, y = mean_seedset_round), 
             color = "black", alpha = 0.6, size = 2) +
  labs(
    title = "Predicted seedset vs. Visitation Rate",
    x = "Mean visitors per minute (scaled)",
    y = "Predicted seedset (mean_seedset_round)"
  ) +
  theme_minimal()


###* POLLEN LIMITATION INDEX

bayesian_pl_noel_w12_p0.5_vm_nh <- brm(
  formula = bf(
    mean_PL_index | weights(PL_index_weight_12) ~ 
      me(mean.visited.flowers.scaled, sd.visited.flowers.scaled) + 
      me(mean.morpho.scaled, sd.morpho.scaled) +
      #me(mean.func.scaled, sd.func.scaled) +
      (1|species)
  ),
  data = c.pl.final.table.4,
  family = zero_one_inflated_beta(),
  prior = c(
    set_prior("normal(0, 0.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 12000,
  warmup = 5000,
  save_pars = save_pars(all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234
)

saveRDS(bayesian_pl_noel_w12_p1_v_nh, file = "../brms_models/bayesian_pl_noel_w12_p1_v_nh.rds")
saveRDS(bayesian_pl_noel_w12_p1_vm_nh, file = "../brms_models/bayesian_pl_noel_w12_p1_vm_nh.rds")
saveRDS(bayesian_pl_noel_w12_p0.5_vm_nh, file = "../brms_models/bayesian_pl_noel_w12_p1_vm_nh.rds")

#not finished
saveRDS(bayesian_pl_noel_w12_p1_vmf_nh, file = "../brms_models/bayesian_pl_noel_w12_p1_vmf_nh.rds")
saveRDS(bayesian_pl_noel_w12_p1_null_nh, file = "../brms_models/bayesian_pl_noel_w12_p1_null_nh.rds")
saveRDS(bayesian_pl_noel_w12_p0.5_v_nh, file = "../brms_models/bayesian_pl_noel_w12_p0.5_v_nh.rds")

bayesian_pl_noel_w12_p1_v_nh <- readRDS("../brms_models/bayesian_pl_noel_w12_p1_v_nh.rds")
bayesian_pl_noel_w12_p1_vm_nh <- readRDS("../brms_models/bayesian_pl_noel_w12_p1_vm_nh.rds")
bayesian_pl_noel_w12_p1_vmf_nh <- readRDS("../brms_modelsbayesian_pl_noel_w12_p1_vmf_nh.rds")
bayesian_pl_noel_w12_p1_null_nh <- readRDS("../brms_models/bayesian_pl_noel_w12_p1_null_nh.rds")
bayesian_pl_noel_w12_p0.5_v_nh <- readRDS("../brms_models/bayesian_pl_noel_w12_p0.5_v_nh.rds")

bayes_R2(bayesian_pl_noel_w12_p1_null_nh)
bayes_R2(bayesian_pl_noel_w12_p1_v_nh)
bayes_R2(bayesian_pl_noel_w12_p1_vm_nh)
bayes_R2(bayesian_pl_noel_w12_p1_vmf_nh)
bayes_R2(bayesian_pl_noel_w12_p1_mf_nh)
bayes_R2(bayesian_pl_noel_w12_p0.5_v_nh)



fit <- bayesian_pl_noel_w12_p0.5_vm_nh
fit <- bayesian_pl_noel_w12_p1_v_nh
fit <- bayesian_pl_noel_w12_pnone_nh
fit <- bayesian_pl_noel_noel_w12_p1_vm_nh
fit <- bayesian_pl_noel_w12_pnone_nh
fit <- bayesian_pl_noel_w12_pnone_nh
fit <- bayesian_pl_noel_w12_pnone_nh




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



###* AUTOGAMY INDEX
###* 
###* 

formula <- bf(mean_ao_index | weights(ao_index_weight_12) ~ me(mean.visited.flowers.scaled, sd.visited.flowers.scaled) + (1|species) + (1|elevation))

bayesian_ao_w12_p1_nh <- brm(
  formula = formula,
  data = ao.final.table,
  family = zero_inflated_beta(),
  prior = c(
    set_prior("normal(0, 1)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 12000,
  warmup = 5000,
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234
)

saveRDS(bayesian_ao_w12_pnone_nh, file = "../brms_models/bayesian_ao_w12_pnone_nh.rds")
saveRDS(bayesian_ao_w12_p2.5_nh, file = "../brms_models/bayesian_ao_w12_p2.5_nh.rds")
saveRDS(bayesian_ao_w12_p1_nh, file = "../brms_models/bayesian_ao_w12_p1_nh.rds")
saveRDS(bayesian_ao_w12_p0.5_nh, file = "../brms_models/bayesian_ao_w12_p0.5_nh.rds")
bayesian_ao_w12_pnone_nh <- readRDS("../brms_models/bayesian_ao_w12_pnone_nh.rds")
bayesian_ao_w12_p2.5_nh <- readRDS("../brms_models/bayesian_ao_w12_p2.5_nh.rds")
bayesian_ao_w12_p1_nh <- readRDS("../brms_models/bayesian_ao_w12_p1_nh.rds")
bayesian_ao_w12_p0.5_nh <- readRDS("../brms_models/bayesian_ao_w12_p0.5_nh.rds")


fit <- bayesian_ao_w12_pnone_nh
fit <- bayesian_ao_w12_p2.5_nh
fit <- bayesian_ao_w12_p1_nh
fit <- bayesian_ao_w12_p0.5_nh


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
             aes(x = mean.visited.flowers.scaled, y = mean_ao_index), 
             color = "black", alpha = 0.6, size = 2) +
  labs(
    title = "Predicted AO index vs. Visitation Rate",
    x = "Mean visitors per minute (scaled)",
    y = "Predicted AO index (mean_AO_index)"
  ) +
  theme_minimal()





formula <- bf(mean_go_index | weights(go_index_weight_12) ~ me(mean.visited.flowers.scaled, sd.visited.flowers.scaled) + (1|species) + (1|elevation))

bayesian_go_w12_p0.5_12k_nh <- brm(
  formula = formula,
  data = go.final.table,
  family = zero_one_inflated_beta(),
  prior = c(
    set_prior("normal(0, 0.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 12000,
  warmup = 5000,
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234
)

saveRDS(bayesian_go_w12_p1_12k_nh, file = "../brms_models/bayesian_go_w12_p1_12k_nh.rds")
saveRDS(bayesian_go_w12_p0.5_12k_nh, file = "../brms_models/bayesian_go_w12_p0.5_12k_nh.rds")saveRDS(bayesian_go_w12_p2.5_12k_nh, file = "../brms_models/bayesian_go_w12_p2.5_12k_nh.rds")
saveRDS(bayesian_go_w12_pnone_12k_nh, file = "../brms_models/bayesian_go_w12_pnone_12k_nh.rds")
bayesian_go_w12_p1_12k_nh <- readRDS("../brms_models/bayesian_go_w12_p1_12k_nh.rds")
bayesian_go_w12_p0.5_12k_nh <- readRDS("../brms_models/bayesian_go_w12_p0.5_12k_nh.rds")
bayesian_go_w12_p2.5_12k_nh <- readRDS("../brms_models/bayesian_go_w12_p2.5_12k_nh.rds")
bayesian_go_w12_pnone_12k_nh <- readRDS("../brms_models/bayesian_go_w12_pnone_12k_nh.rds")

fit <- bayesian_go_w12_p0.5_12k_nh
fit <- bayesian_go_w12_p1_12k_nh
fit <- bayesian_go_w12_p2.5_12k_nh
fit <- bayesian_go_w12_pnone_12k_nh

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
