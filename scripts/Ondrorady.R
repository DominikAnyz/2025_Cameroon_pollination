1 + 1/(1+sd)
rescale(-x, to = c(1,2))


###* Plotting the trend
data_pred <-
  marginaleffects::avg_predictions(
    fit_seedset_r_p_11sd,
    by = "mean.visitors.scaled"
  ) %>%
  as.data.frame()

View(data_pred)

nejaky_vektor <- insight::get_data(fit_seedset_r_p_11sd) %>%
  dplyr::pull(mean.visitors.scaled)

nejaky_vektor

data_pred_indiv <-
  marginaleffects::predictions(
    fit_seedset_r_p_11sd,
    newdata = datagrid(
      mean.visitors.scaled =seq(from = min(nejaky_vektor),to = max(nejaky_vektor), length.out = 100)
    )
  )

data_pred_indiv

ggplot(data = as.data.frame(data_pred_indiv),
       mapping = aes(x = mean.visitors.scaled, y = estimate)
       ) +
  geom_ribbon(mapping = aes(ymin = conf.low, ymax = conf.high)) +
  geom_line(
  ) 

View(data_pred_indiv)

colnames(data_pred_indiv)

library(marginaleffects)
library(ggeffects)

ggeffects::predict_response(
  fit_seedset_r_p_11sd, terms = "mean.visitors.scaled",
  interval = "prediction",
  bias_correction = TRUE
)

###* Fitting null model
formula <- bf(mean_seedset_round | weights (seedset_weight) ~ 1 + (1|species) + (1|elevation))
###* We can use it later to asses whether our model is better than the null model





library(marginaleffects)
library(ggplot2)
library(dplyr)
library(insight)



###* FIT SEEDSET

###* fit_seedset_12 
###* prior <- c(set_prior("normal(0, 2.5)", class = "b"), set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"), set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))
fit_seedset_12 <- readRDS("../brms_models/fit_seedset_12.rds")
###* - RHAT = 1.05

###* fit_seedset_12_2
###* prior <- c(set_prior("normal(0, 1)", class = "b"), set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"), set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))
fit_seedset_12_2 <- readRDS("../brms_models/fit_seedset_12_2.rds")

###* fit_seedset_12_3 
###* prior <- c(set_prior("normal(0, 0.5)", class = "b"), set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"), set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))
fit_seedset_12_3 <- readRDS("../brms_models/fit_seedset_12_3.rds")
###* - 1 DIVERGENT

###* fit_seedset_11sd
###* family = poisson
###* prior <- c(set_prior("normal(0, 2.5)", class = "b"), set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"), set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))
fit_seedset_11sd <- readRDS("../brms_models/fit_seedset_11sd.rds")

###* fit_seedset_11sd_2
###* family = negbinomial()
fit_seedset_11sd_2 <- readRDS("../brms_models/fit_seedset_11sd_2.rds")

###* fit_seedset_11sd_3
###* family = poisson
###* prior <- c(set_prior("normal(0, 1)", class = "b"), set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"), set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))
fit_seedset_11sd_3 <- readRDS("../brms_models/fit_seedset_11sd_3.rds")

###* fit_seedset_11sd_4
###* family = poisson
###* prior <- c(set_prior("normal(0, 0.5)", class = "b"), set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"), set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))
fit_seedset_11sd_4 <- readRDS("../brms_models/fit_seedset_11sd_4.rds")

###* fit_seedset_12_new_3 
###* prior <- c(set_prior("normal(0, 0.5)", class = "b"), set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"), set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))
fit_seedset_12_new_3 <- readRDS("../brms_models/fit_seedset_12_new_3.rds")


###* fit_seedset_12_new_nh_3 
###* prior <- c(set_prior("normal(0, 0.5)", class = "b"), set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"), set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))
fit_seedset_12_new_nh_3 <- readRDS("../brms_models/fit_seedset_12_new_nh_3.rds")

fit <- fit_seedset_12_3
fit <- fit_seedset_12_new_3
fit <- fit_seedset_12_new_nh_3

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
  pull(mean.visitors.scaled)
# Generate marginal predictions across the predictor range
data_pred_indiv <- predictions(
  fit,
  newdata = datagrid(
    mean.visitors.scaled = seq(from = min(x_range), to = max(x_range), length.out = 100)
  ))
# Extract the original observed data
observed_data <- insight::get_data(fit)
# Create ggplot
ggplot(data = as.data.frame(data_pred_indiv), 
       aes(x = mean.visitors.scaled, y = estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "steelblue") + # Uncertainty ribbon
  geom_line(color = "steelblue", size = 1) + # Fitted line
  geom_point(data = observed_data, # Observed data points
             aes(x = mean.visitors.scaled, y = mean_seedset_round), 
             color = "black", alpha = 0.6, size = 2) +
  labs(
    title = "Predicted seedset vs. Visitation Rate",
    x = "Mean visitors per minute (scaled)",
    y = "Predicted seedset (mean_seedset_round)"
  ) +
  theme_minimal()


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


###* FIT PL INDEX




###* fit_PL_12
###* prior = c(set_prior("normal(0, 2)", class = "b")
fit_PL_12 <- readRDS("../brms_models/fit_PL_12.rds")
###* RHAT = 1.02

###* fit_PL_12_2
###* prior = c(set_prior("normal(0, 1)", class = "b")
fit_PL_12_2 <- readRDS("../brms_models/fit_PL_12_2.rds")

###* fit_PL_12_3
###* prior = c(set_prior("normal(0, 0.5)", class = "b")
fit_PL_12_3 <- readRDS("../brms_models/fit_PL_12_3.rds")

# ###* fit_PL_11sd
# ###* set_prior("normal(0, 2.5)", class = "b")
# fit_PL_11sd <- readRDS("../brms_models/fit_PL_11sd.rds")
# ###* RHAT = 1.42

# ###* fit_PL_11sd_2
# ###* set_prior("normal(0, 1)", class = "b")
# fit_PL_11sd_2 <- readRDS("../brms_models/fit_PL_11sd_2.rds")
# ###* RHAT = 1.41

# ###* fit_PL_11sd_3
# ###* set_prior("normal(0, 0.5)", class = "b")
# fit_PL_11sd_3 <- readRDS("../brms_models/fit_PL_11sd_3.rds")
# ###* RHAT = 1.16

###* fit_PL_12_new_3
###* prior = c(set_prior("normal(0, 0.5)", class = "b")
fit_PL_12_new_3 <- readRDS("../brms_models/fit_PL_12_new_3.rds")

###* fit_PL_12_new_nh_3
###* prior = c(set_prior("normal(0, 0.5)", class = "b")
fit_PL_12_new_nh_3 <- readRDS("../brms_models/fit_PL_12_new_nh_3.rds")

fit <- fit_PL_12_3
fit <- fit_PL_12_new_3
fit <- fit_PL_12_new_nh_3

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
  pull(mean.visitors.scaled)
# Generate marginal predictions across the predictor range
data_pred_indiv <- predictions(
  fit,
  newdata = datagrid(
    mean.visitors.scaled = seq(from = min(x_range), to = max(x_range), length.out = 100)
  ))
# Extract the original observed data
observed_data <- insight::get_data(fit)
# Create ggplot
ggplot(data = as.data.frame(data_pred_indiv), 
       aes(x = mean.visitors.scaled, y = estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "steelblue") + # Uncertainty ribbon
  geom_line(color = "steelblue", size = 1) + # Fitted line
  geom_point(data = observed_data, # Observed data points
             aes(x = mean.visitors.scaled, y = mean_PL_index), 
             color = "black", alpha = 0.6, size = 2) +
  labs(
    title = "Predicted PL index vs. Visitation Rate",
    x = "Mean visitors per minute (scaled)",
    y = "Predicted PL index (mean_PL_index)"
  ) +
  theme_minimal()


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




###* FIT AO INDEX

###* fit_ao_12
###* set_prior("normal(0, 2.5)", class = "b")
fit_ao_12 <- readRDS("../brms_models/fit_ao_12.rds")

# ###* fit_ao_12_2
# ###* set_prior("normal(0, 1)", class = "b")
# ###* iter = 10000, warmup = 4000
# fit_ao_12_2 <- readRDS("../brms_models/fit_ao_12_2.rds")
# ###* RHAT = 1.12

###* fit_ao_12_3
###* set_prior("normal(0, 0.5)", class = "b")
###* iter = 10000, warmup = 4000
fit_ao_12_3 <- readRDS("../brms_models/fit_ao_12_3.rds")

# ###* fit_ao_11sd
# ###* set_prior("normal(0, 2.5)", class = "b")
# fit_ao_11sd <- readRDS("../brms_models/fit_ao_11sd.rds")
# ###* RHAT = 1.65

# ###* fit_ao_11sd_2
# ###* set_prior("normal(0, 1)", class = "b")
# fit_ao_11sd_2 <- readRDS("../brms_models/fit_ao_11sd_2.rds")
# ###* RHAT = 1.25

###* fit_ao_11sd_3
###* set_prior("normal(0, 0.5)", class = "b")
fit_ao_11sd_3 <- readRDS("../brms_models/fit_ao_11sd_3.rds")

###* fit_ao_12_new_3
###* set_prior("normal(0, 0.5)", class = "b")
###* iter = 10000, warmup = 4000
fit_ao_12_new_3 <- readRDS("../brms_models/fit_ao_12_new_3.rds")

###* fit_ao_12_new_nh_3
###* set_prior("normal(0, 0.5)", class = "b")
###* iter = 10000, warmup = 4000
fit_ao_12_new_nh_3 <- readRDS("../brms_models/fit_ao_12_new_nh_3.rds")

###* fit_ao_12_new_nh
###* set_prior("normal(0, 0.5)", class = "b")
###* iter = 10000, warmup = 4000
fit_ao_12_new_nh <- readRDS("../brms_models/fit_ao_12_new_nh.rds")

fit <- fit_ao_12_3
fit <- fit_ao_12_new_3
fit <- fit_ao_12_new_nh_3
fit <- fit_ao_12_new_nh

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
  pull(mean.visitors.scaled)
# Generate marginal predictions across the predictor range
data_pred_indiv <- predictions(
  fit,
  newdata = datagrid(
    mean.visitors.scaled = seq(from = min(x_range), to = max(x_range), length.out = 100)
  ))
# Extract the original observed data
observed_data <- insight::get_data(fit)
# Create ggplot
ggplot(data = as.data.frame(data_pred_indiv), 
       aes(x = mean.visitors.scaled, y = estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "steelblue") + # Uncertainty ribbon
  geom_line(color = "steelblue", size = 1) + # Fitted line
  geom_point(data = observed_data, # Observed data points
             aes(x = mean.visitors.scaled, y = mean_ao_index), 
             color = "black", alpha = 0.6, size = 2) +
  labs(
    title = "Predicted AO index vs. Visitation Rate",
    x = "Mean visitors per minute (scaled)",
    y = "Predicted AO index (mean_AO_index)"
  ) +
  theme_minimal()




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
    x = "Mean visited flowers per minute (scaled)",
    y = "Predicted AO index (mean_AO_index)"
  ) +
  theme_minimal()







###* FIT GO INDEX

###* fit_go_12
###* set_prior("normal(0, 2.5)", class = "b")
fit_go_12 <- readRDS("../brms_models/fit_go_12.rds")

###* fit_go_12_2
###* set_prior("normal(0, 1)", class = "b")
###* iter = 10000, warmup = 4000
fit_go_12_2 <- readRDS("../brms_models/fit_go_12_2.rds")

###* fit_go_12_3
###* set_prior("normal(0, 0.5)", class = "b")
###* iter = 10000, warmup = 4000
fit_go_12_3 <- readRDS("../brms_models/fit_go_12_3.rds")

###* fit_go_11sd
###* set_prior("normal(0, 2.5)", class = "b")
fit_go_11sd <- readRDS("../brms_models/fit_go_11sd.rds")

###* fit_go_11sd_2
###* set_prior("normal(0, 1)", class = "b")
fit_go_11sd_2 <- readRDS("../brms_models/fit_go_11sd_2.rds")

###* fit_go_11sd_3
###* set_prior("normal(0, 0.5)", class = "b")
fit_go_11sd_3 <- readRDS("../brms_models/fit_go_11sd_3.rds")


fit <- fit_go_11sd_3

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
  pull(mean.visitors.scaled)
# Generate marginal predictions across the predictor range
data_pred_indiv <- predictions(
  fit,
  newdata = datagrid(
    mean.visitors.scaled = seq(from = min(x_range), to = max(x_range), length.out = 100)
  ))
# Extract the original observed data
observed_data <- insight::get_data(fit)
# Create ggplot
ggplot(data = as.data.frame(data_pred_indiv), 
       aes(x = mean.visitors.scaled, y = estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "steelblue") + # Uncertainty ribbon
  geom_line(color = "steelblue", size = 1) + # Fitted line
  geom_point(data = observed_data, # Observed data points
             aes(x = mean.visitors.scaled, y = mean_go_index), 
             color = "black", alpha = 0.6, size = 2) +
  labs(
    title = "Predicted GO index vs. Visitation Rate",
    x = "Mean visitors per minute (scaled)",
    y = "Predicted GO index (mean_GO_index)"
  ) +
  theme_minimal()
