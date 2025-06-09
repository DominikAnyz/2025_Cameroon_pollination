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

formula <- bf(mean_seedset_round | weights (seedset_weight_12) ~ me(mean.visited.flowers.scaled, sd.visited.flowers.scaled) + (1|species) + (1|elevation))
###* Setting priors
prior <- c(
  set_prior("normal(0, 0.5)", class = "b"), 
  set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"), 
  set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))

bayesian_seedset_w12_p0.5_h4 <- brm(
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
  seed = 1234 # Addition to code
)

###* Trying version h4 version for seedset, which has info about Hypericum from
###* elevation 4000 / no produces seeds. Only usable in seedset though

saveRDS(bayesian_seedset_w12_p0.5_h4, file = "../brms_models/bayesian_seedset_w12_p0.5_h4.rds")
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


formula <- bf(mean_PL_index | weights (PL_index_weight_12) ~ me(mean.visited.flowers.scaled, sd.visited.flowers.scaled) + (1|species) + (1|elevation))
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

bayesian_pl_w12_p0.75_12k <- brm(
  formula = bf(
    mean_PL_index | weights(PL_index_weight_12) ~ me(mean.visited.flowers.scaled, sd.visited.flowers.scaled) + (1|species) + (1|elevation)
  ),
  data = c.pl.final.table.5,
  family = zero_one_inflated_beta(),
  prior = c(
    set_prior("normal(0, 0.75)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
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

saveRDS(bayesian_pl_w12_pnone_12k, file = "../brms_models/bayesian_pl_w12_pnone_12k.rds")
saveRDS(bayesian_pl_w12_p10_12k, file = "../brms_models/bayesian_pl_w12_p10_12k.rds")
saveRDS(bayesian_pl_w12_p5_12k, file = "../brms_models/bayesian_pl_w12_p5_12k.rds")
saveRDS(bayesian_pl_w12_p2.5_12k, file = "../brms_models/bayesian_pl_w12_p2.5_12k.rds")
saveRDS(bayesian_pl_w12_p1_12k, file = "../brms_models/bayesian_pl_w12_p1_12k.rds")
saveRDS(bayesian_pl_w12_p0.75_12k, file = "../brms_models/bayesian_pl_w12_p0.75_12k.rds")
saveRDS(bayesian_pl_w12_p0.5, file = "../brms_models/bayesian_pl_w12_p0.5.rds")
bayesian_pl_w12_pnone_12k <- readRDS("../brms_models/bayesian_pl_w12_pnone_12k.rds")
bayesian_pl_w12_p10_12k <- readRDS("../brms_models/bayesian_pl_w12_p10_12k.rds")
bayesian_pl_w12_p5_12k <- readRDS("../brms_models/bayesian_pl_w12_p5_12k.rds")
bayesian_pl_w12_p2.5_12k <- readRDS("../brms_models/bayesian_pl_w12_p2.5_12k.rds")
bayesian_pl_w12_p1_12k <- readRDS("../brms_models/bayesian_pl_w12_p1_12k.rds")
bayesian_pl_w12_p0.75_12k <- readRDS("../brms_models/bayesian_pl_w12_p0.75_12k.rds")
bayesian_pl_w12_p0.5 <- readRDS("../brms_models/bayesian_pl_w12_p0.5.rds")

fit <- bayesian_pl_w12_pnone_12k
fit <- bayesian_pl_w12_p10_12k
fit <- bayesian_pl_w12_p5_12k
fit <- bayesian_pl_w12_p2.5_12k
fit <- bayesian_pl_w12_p1_12k
fit <- bayesian_pl_w12_p0.75_12k
fit <- bayesian_pl_w12_p0.5

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

formula <- bf(mean_ao_index | weights(ao_index_weight_12) ~ me(mean.visited.flowers.scaled, sd.visited.flowers.scaled) + (1|species) + (1|elevation))

bayesian_ao_w12_p2.5_12k <- brm(
  formula = formula,
  data = ao.final.table,
  family = zero_inflated_beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
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

saveRDS(bayesian_ao_w12_p1, file = "../brms_models/bayesian_ao_w12_p1.rds")
saveRDS(bayesian_ao_w12_p1_12k, file = "../brms_models/bayesian_ao_w12_p1_12k.rds")
saveRDS(bayesian_ao_w12_p0.5_12k, file = "../brms_models/bayesian_ao_w12_p0.5_12k.rds")
saveRDS(bayesian_ao_w12_p2.5_12k, file = "../brms_models/bayesian_ao_w12_p2.5_12k.rds")
saveRDS(bayesian_ao_w12_pnone_12k, file = "../brms_models/bayesian_ao_w12_pnone_12k.rds")
bayesian_ao_w12_p1 <- readRDS("../brms_models/bayesian_ao_w12_p1.rds")
bayesian_ao_w12_p1_12k <- readRDS("../brms_models/bayesian_ao_w12_p1_12k.rds")
bayesian_ao_w12_p0.5_12k <- readRDS("../brms_models/bayesian_ao_w12_p0.5_12k.rds")
bayesian_ao_w12_p2.5_12k <- readRDS("../brms_models/bayesian_ao_w12_p2.5_12k.rds")
bayesian_ao_w12_pnone_12k <- readRDS("../brms_models/bayesian_ao_w12_pnone_12k.rds")

fit <- bayesian_ao_w12_p0.5_12k
fit <- bayesian_ao_w12_p1_12k
fit <- bayesian_ao_w12_p2.5_12k
fit <- bayesian_ao_w12_pnone_12k

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

formula <- bf(mean_go_index | weights(go_index_weight_12) ~ me(mean.visited.flowers.scaled, sd.visited.flowers.scaled) + (1|species) + (1|elevation))

bayesian_go_w12_p0.5_12k <- brm(
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

saveRDS(bayesian_go_w12_p1_12k, file = "../brms_models/bayesian_go_w12_p1_12k.rds")
saveRDS(bayesian_go_w12_p0.5_12k, file = "../brms_models/bayesian_go_w12_p0.5_12k.rds")
saveRDS(bayesian_go_w12_p2.5_12k, file = "../brms_models/bayesian_go_w12_p2.5_12k.rds")
saveRDS(bayesian_go_w12_pnone_12k, file = "../brms_models/bayesian_go_w12_pnone_12k.rds")
bayesian_go_w12_p1_12k <- readRDS("../brms_models/bayesian_go_w12_p1_12k.rds")
bayesian_go_w12_p0.5_12k <- readRDS("../brms_models/bayesian_go_w12_p0.5_12k.rds")
bayesian_go_w12_p2.5_12k <- readRDS("../brms_models/bayesian_go_w12_p2.5_12k.rds")
bayesian_go_w12_pnone_12k <- readRDS("../brms_models/bayesian_go_w12_pnone_12k.rds")

fit <- bayesian_go_w12_p0.5_12k
fit <- bayesian_go_w12_p1_12k
fit <- bayesian_go_w12_p2.5_12k
fit <- bayesian_go_w12_pnone_12k

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
