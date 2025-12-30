source("scripts/setup.R")
source("scripts/Bayesian setup.R")

c.pl.final.table$elevation <- as.factor(c.pl.final.table$elevation)
c.pl.final.table$species <- as.factor(c.pl.final.table$species)
c.pl.final.table <- na.omit(c.pl.final.table)

View(c.pl.final.table)

# c.pl.final.table.2 <- c.pl.final.table %>%
#   mutate(mean_seedset = ifelse(mean_seedset == 0, 0.001, mean_seedset),
#          mean_PL_index = ifelse(mean_PL_index == 1, 0.999, mean_PL_index))

c.pl.final.table <- c.pl.final.table %>%
  mutate(mean_seedset_round = round(mean_seedset))

c.pl.final.table <- c.pl.final.table %>%
  mutate(log_mean_seedset = log(mean_seedset),
         se_log_mean_seedset = se_seedset / mean_seedset,
         log_visitors = log(mean.visitors.per.minute),
         sd_log_visitors = sd.visitors.per.minute / mean.visitors.per.minute,
         mean.scaled.plus.one = mean.visitors.scaled + 1,
         sd.scaled.plus.one = sd.visitors.scaled + 1,
         log.m.s.p.o = log(mean.scaled.plus.one),
         log.s.s.p.o = log(sd.scaled.plus.one)
  )

c.pl.final.table.2 <- c.pl.final.table

c.pl.final.table <- na.omit(c.pl.final.table)

library(scales)  # for rescale()

c.pl.final.table.4 <- c.pl.final.table.2 %>%
  mutate(
    # 1. Inverse-transformed weight: 1 / (1 + original)
    seedset_weight_11sd = 1 / (1 + sd_seedset),
    # 2. Rescaled negative weights between 1 and 2
    seedset_weight_12 = rescale(-sd_seedset, to = c(1, 2))
  )

View(c.pl.final.table.4)

# mean_real_weight_seedset <- mean(
#   c.pl.final.table$seedset_weight[c.pl.final.table$seedset_weight != 1],
#   na.rm = TRUE
# )
# 
# mean_real_weight_PL <- mean(
#   c.pl.final.table$PL_index_weight[c.pl.final.table$PL_index_weight != 1],
#   na.rm = TRUE
# )
# 
# c.pl.final.table.3 <- c.pl.final.table.2 %>%
#   mutate(
#     seedset_weight = ifelse(seedset_weight == 1, mean_real_weight_seedset, seedset_weight),
#     PL_index_weight = ifelse(PL_index_weight == 1, mean_real_weight_PL, PL_index_weight)
#   )

View(c.pl.final.table.3)

ao.final.table <- na.omit(ao.final.table)
ao.final.table$elevation <- as.factor(ao.final.table$elevation)
ao.final.table$species <- as.factor(ao.final.table$species)

ao.final.table <- ao.final.table %>%
  mutate(
    log_visitors = log(mean.visitors.per.minute),
    sd_log_visitors = sd.visitors.per.minute / mean.visitors.per.minute,
    mean.scaled.plus.one = mean.visitors.scaled + 1,
    sd.scaled.plus.one = sd.visitors.scaled + 1,
    log.m.s.p.o = log(mean.scaled.plus.one),
    log.s.s.p.o = log(sd.scaled.plus.one)
  )

mean_real_weight_ao <- mean(
  ao.final.table$ao_index_weight[ao.final.table$ao_index_weight != 1],
  na.rm = TRUE
)

ao.final.table.2 <- ao.final.table %>%
  mutate(
    ao_index_weight = ifelse(ao_index_weight == 1, mean_real_weight_ao, ao_index_weight)
  )

View(ao.final.table.2)

mean_real_weight_go <- mean(
  go.final.table$go_index_weight[go.final.table$go_index_weight != 1 & go.final.table$go_index_weight != 0],
  na.rm = TRUE
)

go.final.table.2 <- ao.pl.final.table %>%
  mutate(
    go_index_weight = ifelse(go_index_weight %in% c(0,1), mean_real_weight_go, go_index_weight)
  )
















###* SEEDSET FULL








formula <- bf(mean_seedset_round | weights (seedset_weight) ~ me(log_visitors, sd_log_visitors) + (1|species) + (1|elevation))
###* Setting priors
prior <- c(set_prior("normal(0, 2.5)", class = "b"), set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"), set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))

fit_seedset_nr_p_m <- brm(
  formula = formula,
  data = c.pl.final.table.3,
  chains = 10,
  cores = 10,
  iter = 8000,
  warmup = 3000,
  family = poisson(),
  prior = prior,
  save_pars = save_pars(all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234 # Addition to code
)

saveRDS(fit_seedset_nr_p_m, file = "../brms_models/fit_seedset_nr_p_m.rds")
fit_seedset_nr_p_m <- readRDS("../brms_models/fit_seedset_nr_p_m.rds")
summary(fit_seedset_nr_p_m)
pp_check(fit_seedset_nr_p_m)
pp_check(fit_seedset_nr_p_m, type = "stat")
bayes_R2(fit_seedset_nr_p_m)
plot(fit_seedset_nr_p_m)
fitted_values <- fitted(fit_seedset_nr_p_m)
residuals <- residuals(fit_seedset_nr_p_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")

formula <- bf(mean_seedset_round | weights (seedset_weight) ~ me(mean.visitors.per.minute, sd.visitors.per.minute) + (1|species) + (1|elevation))
###* Setting priors
prior <- c(set_prior("normal(0, 2.5)", class = "b"), set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"), set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))

fit_seedset_nr_p_m <- brm(
  formula = formula,
  data = c.pl.final.table.3,
  chains = 10,
  cores = 10,
  iter = 8000,
  warmup = 3000,
  family = poisson(),
  prior = prior,
  save_pars = save_pars(all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234 # Addition to code
)

saveRDS(fit_seedset_nr_p_m, file = "../brms_models/fit_seedset_nr_p_m.rds")
fit_seedset_nr_p_m <- readRDS("../brms_models/fit_seedset_nr_p_m.rds")
summary(fit_seedset_nr_p_m)
pp_check(fit_seedset_nr_p_m)
pp_check(fit_seedset_nr_p_m, type = "stat")
bayes_R2(fit_seedset_nr_p_m)
plot(fit_seedset_nr_p_m)
fitted_values <- fitted(fit_seedset_nr_p_m)
residuals <- residuals(fit_seedset_nr_p_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")

formula <- bf(mean_seedset_round | weights (seedset_weight) ~ me(mean.visitors.scaled, sd.visitors.scaled) + (1|species) + (1|elevation))

###* Setting priors
prior <- c(set_prior("normal(0, 2.5)", class = "b"), set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"), set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))

fit_seedset_r_p_m <- brm(
  formula = formula,
  data = c.pl.final.table.3,
  chains = 10,
  cores = 10,
  iter = 8000,
  warmup = 3000,
  family = poisson(),
  prior = prior,
  save_pars = save_pars(all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234 # Addition to code
)

saveRDS(fit_seedset_r_p_m, file = "../brms_models/fit_seedset_r_p_m.rds")
fit_seedset_r_p_m <- readRDS("../brms_models/fit_seedset_r_p_m.rds")
summary(fit_seedset_r_p_m)
pp_check(fit_seedset_r_p_m)
pp_check(fit_seedset_r_p_m, type = "stat")
bayes_R2(fit_seedset_r_p_m)
plot(fit_seedset_r_p_m)
fitted_values <- fitted(fit_seedset_r_p_m)
residuals <- residuals(fit_seedset_r_p_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")

formula <- bf(mean_seedset_round | weights (seedset_weight) ~ me(log.m.s.p.o, log.s.s.p.o) + (1|species) + (1|elevation))
###* Setting priors
prior <- c(set_prior("normal(0, 2.5)", class = "b"), set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"), set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))

fit_seedset_r2_p_m <- brm(
  formula = formula,
  data = c.pl.final.table.3,
  chains = 10,
  cores = 10,
  iter = 8000,
  warmup = 3000,
  family = poisson(),
  prior = prior,
  save_pars = save_pars(all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234 # Addition to code
)

saveRDS(fit_seedset_r2_p_m, file = "../brms_models/fit_seedset_r2_p_m.rds")
fit_seedset_r2_p_m <- readRDS("../brms_models/fit_seedset_r2_p_m.rds")
summary(fit_seedset_r2_p_m)
pp_check(fit_seedset_r2_p_m)
pp_check(fit_seedset_r2_p_m, type = "stat")
bayes_R2(fit_seedset_r2_p_m)
plot(fit_seedset_r2_p_m)
fitted_values <- fitted(fit_seedset_r2_p_m)
residuals <- residuals(fit_seedset_r2_p_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")














###* SEEDSET PARTIAL













formula <- bf(mean_seedset_round | weights (seedset_weight) ~ mean.visitors.per.minute + (1|species) + (1|elevation))
###* Setting priors
prior <- c(set_prior("normal(0, 2.5)", class = "b"), set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"), set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))

fit_seedset_nr_0_partial_o_m <- brm(
  formula = formula,
  data = c.pl.final.table.3,
  chains = 10,
  cores = 10,
  iter = 6000,
  warmup = 2000,
  family = poisson(),
  prior = prior,
  save_pars = save_pars(all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234 # Addition to code
)

saveRDS(fit_seedset_nr_0_partial_o_m, file = "../brms_models/fit_seedset_nr_p_partial_o_m.rds")
fit_seedset_nr_p_partial_o_m <- readRDS("../brms_models/fit_seedset_nr_p_partial_o_m.rds")
summary(fit_seedset_nr_p_partial_o_m)
pp_check(fit_seedset_nr_p_partial_o_m)
pp_check(fit_seedset_nr_p_partial_o_m, type = "stat")
bayes_R2(fit_seedset_nr_p_partial_o_m)
plot(fit_seedset_nr_p_partial_o_m)
fitted_values <- fitted(fit_seedset_nr_p_partial_o_m)
residuals <- residuals(fit_seedset_nr_p_partial_o_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")


formula <- bf(mean_seedset_round | weights (seedset_weight) ~ log_visitors + (1|species) + (1|elevation))
###* Setting priors
prior <- c(set_prior("normal(0, 2.5)", class = "b"), set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"), set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))

fit_seedset_nr_0_partial_m <- brm(
  formula = formula,
  data = c.pl.final.table.3,
  chains = 10,
  cores = 10,
  iter = 6000,
  warmup = 2000,
  family = poisson(),
  prior = prior,
  save_pars = save_pars(all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234 # Addition to code
)

saveRDS(fit_seedset_nr_0_partial_m, file = "../brms_models/fit_seedset_nr_p_partial_m.rds")
fit_seedset_nr_p_partial_m <- readRDS("../brms_models/fit_seedset_nr_p_partial_m.rds")
summary(fit_seedset_nr_p_partial_m)
pp_check(fit_seedset_nr_p_partial_m)
pp_check(fit_seedset_nr_p_partial_m, type = "stat")
bayes_R2(fit_seedset_nr_p_partial_m)
plot(fit_seedset_nr_p_partial_m)
fitted_values <- fitted(fit_seedset_nr_p_partial_m)
residuals <- residuals(fit_seedset_nr_p_partial_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")

formula <- bf(mean_seedset_round | weights (seedset_weight) ~ mean.visitors.scaled + (1|species) + (1|elevation))
###* Setting priors
prior <- c(set_prior("normal(0, 2.5)", class = "b"), set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"), set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))

fit_seedset_r_p_partial_m <- brm(
  formula = formula,
  data = c.pl.final.table.3,
  chains = 10,
  cores = 10,
  iter = 8000,
  warmup = 3000,
  family = poisson(),
  prior = prior,
  save_pars = save_pars(all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234 # Addition to code
)

saveRDS(fit_seedset_r_p_partial_m, file = "../brms_models/fit_seedset_r_p_partial_m.rds")
fit_seedset_r_p_partial_m <- readRDS("../brms_models/fit_seedset_r_p_partial_m.rds")
summary(fit_seedset_r_p_partial_m)
pp_check(fit_seedset_r_p_partial_m)
pp_check(fit_seedset_r_p_partial_m, type = "stat")
bayes_R2(fit_seedset_r_p_partial_m)
plot(fit_seedset_r_p_partial_m)
fitted_values <- fitted(fit_seedset_r_p_partial_m)
residuals <- residuals(fit_seedset_r_p_partial_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")

formula <- bf(mean_seedset_round | weights (seedset_weight) ~ log.m.s.p.o + (1|species) + (1|elevation))
###* Setting priors
prior <- c(set_prior("normal(0, 2.5)", class = "b"), set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"), set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))

fit_seedset_r2_p_partial_m <- brm(
  formula = formula,
  data = c.pl.final.table.3,
  chains = 10,
  cores = 10,
  iter = 8000,
  warmup = 3000,
  family = poisson(),
  prior = prior,
  save_pars = save_pars(all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 123 # Addition to code
)

saveRDS(fit_seedset_r2_p_partial_m, file = "../brms_models/fit_seedset_r2_p_partial_m.rds")
fit_seedset_r2_p_partial_m <- readRDS("../brms_models/fit_seedset_r2_p_partial_m.rds")
summary(fit_seedset_r2_p_partial_m)
pp_check(fit_seedset_r2_p_partial_m)
pp_check(fit_seedset_r2_p_partial_m, type = "stat")
bayes_R2(fit_seedset_r2_p_partial_m)
plot(fit_seedset_r2_p_partial_m)
fitted_values <- fitted(fit_seedset_r2_p_partial_m)
residuals <- residuals(fit_seedset_r2_p_partial_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")










###* PL INDEX FULL









formula <- bf(mean_PL_index | weights(PL_index_weight) ~ me(mean.visitors.per.minute, sd.visitors.per.minute) + (1|species) + (1|elevation)
)
fit_PL_nr_y_o_m <- brm(
  formula = formula,
  data = c.pl.final.table.3,
  family = zero_one_inflated_beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 8000,
  warmup = 3000,
  save_pars = save_pars(all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234
)
saveRDS(fit_PL_nr_y_o_m, file = "../brms_models/fit_PL_nr_y_o_m.rds")
fit_PL_nr_y_o_m <- readRDS("../brms_models/fit_PL_nr_y_o_m.rds")
summary(fit_PL_nr_y_o_m)
pp_check(fit_PL_nr_y_o_m, type = "dens_overlay")
pp_check(fit_PL_nr_y_o_m, type = "stat")
bayes_R2(fit_PL_nr_y_o_m)
plot(fit_PL_nr_y_o_m)
fitted_values <- fitted(fit_PL_nr_y_o_m)
residuals <- residuals(fit_PL_nr_y_o_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")


formula <- bf(mean_PL_index | weights(PL_index_weight) ~ me(log_visitors, sd_log_visitors) + (1|species) + (1|elevation)
)

fit_PL_nr_y_l_m <- brm(
  formula = formula,
  data = c.pl.final.table.3,
  family = zero_one_inflated_beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 6000,
  warmup = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
saveRDS(fit_PL_nr_y_l_m, file = "../brms_models/fit_PL_nr_y_l_m.rds")
fit_PL_nr_y_l_m <- readRDS("../brms_models/fit_PL_nr_y_l_m.rds")
summary(fit_PL_nr_y_l_m)
pp_check(fit_PL_nr_y_l_m, type = "dens_overlay")
pp_check(fit_PL_nr_y_l_m, type = "stat")
bayes_R2(fit_PL_nr_y_l_m)
plot(fit_PL_nr_y_l_m)
fitted_values <- fitted(fit_PL_nr_y_l_m)
residuals <- residuals(fit_PL_nr_y_l_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")


formula <- bf(mean_PL_index | weights(PL_index_weight) ~ me(mean.visitors.scaled, sd.visitors.scaled) + (1|species) + (1|elevation))

fit_PL_r_y_m <- brm(
  formula = formula,
  data = c.pl.final.table,
  family = zero_one_inflated_beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 8000,
  warmup = 3000,
  save_pars = save_pars(all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234
)

saveRDS(fit_PL_r_y_m, file = "../brms_models/fit_PL_r_y_m.rds")
fit_PL_r_y_m <- readRDS("../brms_models/fit_PL_r_y_m.rds")
summary(fit_PL_r_y_m)
pp_check(fit_PL_r_y_m, type = "dens_overlay")
pp_check(fit_PL_r_y_m, type = "stat")
bayes_R2(fit_PL_r_y_m)
plot(fit_PL_r_y_m)
fitted_values <- fitted(fit_PL_r_y_m)
residuals <- residuals(fit_PL_r_y_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")



formula <- bf(mean_PL_index | weights(PL_index_weight) ~ me(log.m.s.p.o, log.s.s.p.o) + (1|species) + (1|elevation))

fit_PL_r2_y_m <- brm(
  formula = formula,
  data = c.pl.final.table.3,
  family = zero_one_inflated_beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 8000,
  warmup = 3000,
  save_pars = save_pars(all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234
)
saveRDS(fit_PL_r2_y_m, file = "../brms_models/fit_PL_r2_y_m.rds")
fit_PL_r2_y_m <- readRDS("../brms_models/fit_PL_r2_y_m.rds")
summary(fit_PL_r2_y_m)
pp_check(fit_PL_r2_y_m, type = "dens_overlay")
pp_check(fit_PL_r2_y_m, type = "stat")
bayes_R2(fit_PL_r2_y_m)
plot(fit_PL_r2_y_m)
fitted_values <- fitted(fit_PL_r2_y_m)
residuals <- residuals(fit_PL_r2_y_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")







###* PL INDEX PARTIAL








formula <- bf(mean_PL_index | weights(PL_index_weight) ~ mean.visitors.per.minute + (1|species) + (1|elevation)
)
fit_PL_nr_y_o_partial_m <- brm(
  formula = formula,
  data = c.pl.final.table.3,
  family = zero_one_inflated_beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 7000,
  warmup = 2500,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)

saveRDS(fit_PL_nr_y_o_partial_m, file = "../brms_models/fit_PL_nr_y_o_partial_m.rds")
fit_PL_nr_y_o_partial_m <- readRDS("../brms_models/fit_PL_nr_y_o_partial_m.rds")
summary(fit_PL_nr_y_o_partial_m)
pp_check(fit_PL_nr_y_o_partial_m, type = "dens_overlay")
pp_check(fit_PL_nr_y_o_partial_m, type = "stat")
bayes_R2(fit_PL_nr_y_o_partial_m)
plot(fit_PL_nr_y_o_partial_m)
fitted_values <- fitted(fit_PL_nr_y_o_partial_m)
residuals <- residuals(fit_PL_nr_y_o_partial_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")


formula <- bf(mean_PL_index | weights(PL_index_weight) ~ log_visitors + (1|species) + (1|elevation)
)

fit_PL_nr_y_l_partial_m <- brm(
  formula = formula,
  data = c.pl.final.table.3,
  family = zero_one_inflated_beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 8000,
  warmup = 3000,
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234
)

saveRDS(fit_PL_nr_y_l_partial_m, file = "../brms_models/fit_PL_nr_y_l_partial_m.rds")
fit_PL_nr_y_l_partial_m <- readRDS("../brms_models/fit_PL_nr_y_l_partial_m.rds")
summary(fit_PL_nr_y_l_partial_m)
pp_check(fit_PL_nr_y_l_partial_m, type = "dens_overlay")
pp_check(fit_PL_nr_y_l_partial_m, type = "stat")
bayes_R2(fit_PL_nr_y_l_partial_m)
plot(fit_PL_nr_y_l_partial_m)
fitted_values <- fitted(fit_PL_nr_y_l_partial_m)
residuals <- residuals(fit_PL_nr_y_l_partial_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")


formula <- bf(mean_PL_index | weights(PL_index_weight) ~ mean.visitors.scaled + (1|species) + (1|elevation))

fit_PL_r2_1_partial_m <- brm(
  formula = formula,
  data = c.pl.final.table.3,
  family = zero_one_inflated_beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 8000,
  warmup = 3000,
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234
)

saveRDS(fit_PL_r2_1_partial_m, file = "../brms_models/fit_PL_r2_1_partial_m.rds")
fit_PL_r2_1_partial_m <- readRDS("../brms_models/fit_PL_r2_1_partial_m.rds")
summary(fit_PL_r2_1_partial_m)
pp_check(fit_PL_r2_1_partial_m, type = "dens_overlay")
pp_check(fit_PL_r2_1_partial_m, type = "stat")
bayes_R2(fit_PL_r2_1_partial_m)
plot(fit_PL_r2_1_partial_m)
fitted_values <- fitted(fit_PL_r2_1_partial_m)
residuals <- residuals(fit_PL_r2_1_partial_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")


formula <- bf(mean_PL_index | weights(PL_index_weight) ~ log.m.s.p.o + (1|species) + (1|elevation))

fit_PL_r2_y_partial_m <- brm(
  formula = formula,
  data = c.pl.final.table.3,
  family = zero_one_inflated_beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 6000,
  warmup = 2000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)

saveRDS(fit_PL_r2_y_partial_m, file = "../brms_models/fit_PL_r2_y_partial_m.rds")
fit_PL_r2_y_partial_m <- readRDS("../brms_models/fit_PL_r2_y_partial_m.rds")
summary(fit_PL_r2_y_partial_m)
pp_check(fit_PL_r2_y_partial_m, type = "dens_overlay")
pp_check(fit_PL_r2_y_partial_m, type = "stat")
bayes_R2(fit_PL_r2_y_partial_m)
plot(fit_PL_r2_y_partial_m)
fitted_values <- fitted(fit_PL_r2_y_partial_m)
residuals <- residuals(fit_PL_r2_y_partial_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")














###* AUTOGAMY FULL

















formula <- bf(mean_ao_index | weights(ao_index_weight) ~ me(mean.visitors.per.minute, sd.visitors.per.minute) + (1|species) + (1|elevation))

fit_ao_nr_0_m <- brm(
  formula = formula,
  data = ao.final.table.2,
  family = zero_inflated_beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 6000,
  warmup = 2000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 123
)
saveRDS(fit_ao_nr_0_m, file = "../brms_models/fit_ao_nr_0_m.rds")
fit_ao_nr_0_m <- readRDS("../brms_models/fit_ao_nr_0_m.rds")
summary(fit_ao_nr_0_m)
pp_check(fit_ao_nr_0_m, type = "dens_overlay")
pp_check(fit_ao_nr_0_m, type = "stat")
bayes_R2(fit_ao_nr_0_n)
plot(fit_ao_nr_0_n)
fitted_values <- fitted(fit_ao_nr_0_n)
residuals <- residuals(fit_ao_nr_0_n)
plot(fitted_values, residuals)
abline(h = 0, col = "red")




formula <- bf(mean_ao_index | weights(ao_index_weight) ~ me(log_visitors, sd_log_visitors) + (1|species) + (1|elevation))

fit_ao_nr_l_m <- brm(
  formula = formula,
  data = ao.final.table.2,
  family = zero_inflated_beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 8000,
  warmup = 3000,
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 123
)
saveRDS(fit_ao_nr_l_m, file = "../brms_models/fit_ao_nr_l_m.rds")
fit_ao_nr_l_m <- readRDS("../brms_models/fit_ao_nr_l_m.rds")
summary(fit_ao_nr_l_m)
pp_check(fit_ao_nr_l_m, type = "dens_overlay")
pp_check(fit_ao_nr_l_m, type = "stat")
bayes_R2(fit_ao_nr_l_m)
plot(fit_ao_nr_l_m)
fitted_values <- fitted(fit_ao_nr_l_m)
residuals <- residuals(fit_ao_nr_l_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")



formula <- bf(mean_ao_index | weights(ao_index_weight) ~ me(mean.visitors.scaled, sd.visitors.scaled) + (1|species) + (1|elevation))

fit_ao_r_o_m <- brm(
  formula = formula,
  data = ao.final.table.2,
  family = zero_inflated_beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 8000,
  warmup = 3000,
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234
)
saveRDS(fit_ao_r_o_m, file = "../brms_models/fit_ao_r_o_m.rds")
fit_ao_r_o_m <- readRDS("../brms_models/fit_ao_r_o_m.rds")
summary(fit_ao_r_o_m)
pp_check(fit_ao_r_o_m, type = "dens_overlay")
pp_check(fit_ao_r_o_m, type = "stat")
bayes_R2(fit_ao_r_o_m)
plot(fit_ao_r_o_m)
fitted_values <- fitted(fit_ao_r_o_m)
residuals <- residuals(fit_ao_r_o_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")



formula <- bf(mean_ao_index | weights(ao_index_weight) ~ me(log.m.s.p.o, log.s.s.p.o) + (1|species) + (1|elevation))

fit_ao_r2_o_m <- brm(
  formula = formula,
  data = ao.final.table.2,
  family = zero_inflated_beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 6000,
  warmup = 2000,
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234
)

saveRDS(fit_ao_r2_o_m, file = "../brms_models/fit_ao_r2_o_m.rds")
fit_ao_r2_o_m <- readRDS("../brms_models/fit_ao_r2_o_m.rds")
summary(fit_ao_r2_o_m)
pp_check(fit_ao_r2_o_m, type = "dens_overlay")
pp_check(fit_ao_r2_o_m, type = "stat")
bayes_R2(fit_ao_r2_o_m)
plot(fit_ao_r2_o_m)
fitted_values <- fitted(fit_ao_r2_o_m)
residuals <- residuals(fit_ao_r2_o_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")












###* AUTOGAMY PARTIAL















formula <- bf(mean_ao_index | weights(ao_index_weight) ~ mean.visitors.per.minute + (1|species) + (1|elevation))
fit_ao_nr_0_partial_m <- brm(
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
  iter = 8000,
  warmup = 3000,
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234
)

saveRDS(fit_ao_nr_0_partial_m, file = "../brms_models/fit_ao_nr_0_partial_m.rds")
fit_ao_nr_0_partial_m <- readRDS("../brms_models/fit_ao_nr_0_partial_m.rds")
summary(fit_ao_nr_0_partial_m)
pp_check(fit_ao_nr_0_partial_m, type = "dens_overlay")
pp_check(fit_ao_nr_0_partial_m, type = "stat")
bayes_R2(fit_ao_nr_0_partial_m)
plot(fit_ao_nr_0_partial_m)
fitted_values <- fitted(fit_ao_nr_0_partial_m)
residuals <- residuals(fit_ao_nr_0_partial_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")



formula <- bf(mean_ao_index | weights(ao_index_weight) ~ log_visitors + (1|species) + (1|elevation))
fit_ao_nr_l_partial_m <- brm(
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
  iter = 6000,
  warmup = 2000,
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234
)
saveRDS(fit_ao_nr_l_partial_m, file = "../brms_models/fit_ao_nr_l_partial_m.rds")
fit_ao_nr_l_partial_m <- readRDS("../brms_models/fit_ao_nr_l_partial_m.rds")
summary(fit_ao_nr_l_partial_m)
pp_check(fit_ao_nr_l_partial_m, type = "dens_overlay")
pp_check(fit_ao_nr_l_partial_m, type = "stat")
bayes_R2(fit_ao_nr_l_partial_m)
plot(fit_ao_nr_l_partial_m)
fitted_values <- fitted(fit_ao_nr_l_partial_m)
residuals <- residuals(fit_ao_nr_l_partial_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")




formula <- bf(mean_ao_index | weights(ao_index_weight) ~ mean.visitors.scaled + (1|species) + (1|elevation))
fit_ao_r_o_partial_m <- brm(
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
  iter = 8000,
  warmup = 3000,
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234
)
saveRDS(fit_ao_r_o_partial_m, file = "../brms_models/fit_ao_r_o_partial_m.rds")
fit_ao_r_o_partial_m <- readRDS("../brms_models/fit_ao_r_o_partial_m.rds")
summary(fit_ao_r_o_partial_m)
pp_check(fit_ao_r_o_partial_m, type = "dens_overlay")
pp_check(fit_ao_r_o_partial_m, type = "stat")
bayes_R2(fit_ao_r_o_partial_m)
plot(fit_ao_r_o_partial_m)
fitted_values <- fitted(fit_ao_r_o_partial_m)
residuals <- residuals(fit_ao_r_o_partial_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")



formula <- bf(mean_ao_index | weights(ao_index_weight) ~ log.m.s.p.o + (1|species) + (1|elevation))
fit_ao_r2_o_partial_m <- brm(
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
  iter = 7000,
  warmup = 2500,
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234
)
saveRDS(fit_ao_r2_o_partial_m, file = "../brms_models/fit_ao_r2_o_partial_m.rds")
fit_ao_r2_o_partial_m <- readRDS("../brms_models/fit_ao_r2_o_partial_m.rds")
summary(fit_ao_r2_o_partial_m)
pp_check(fit_ao_r2_o_partial_m, type = "dens_overlay")
pp_check(fit_ao_r2_o_partial_m, type = "stat")
bayes_R2(fit_ao_r2_o_partial_m)
plot(fit_ao_r2_o_partial_m)
fitted_values <- fitted(fit_ao_r2_o_partial_m)
residuals <- residuals(fit_ao_r2_o_partial_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")




















###* GEITONOGAMY FULL
















formula <- bf(mean_go_index | weights(go_index_weight) ~ me(mean.visitors.per.minute, sd.visitors.per.minute) + (1|species) + (1|elevation))
fit_go_nr_0_m <- brm(
  formula = formula,
  data = go.final.table,
  family = zero_one_inflated_beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 8,
  cores = 8,
  iter = 4000,
  warmup = 1000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
saveRDS(fit_go_nr_0_m, file = "../brms_models/fit_go_nr_0_m.rds")
fit_go_nr_0_m <- readRDS("../brms_models/fit_go_nr_0_m.rds")
summary(fit_go_nr_0_m)
pp_check(fit_go_nr_0_m, type = "dens_overlay")
pp_check(fit_go_nr_0_m, type = "stat")
bayes_R2(fit_go_nr_0_m)
plot(fit_go_nr_0_m)
fitted_values <- fitted(fit_go_nr_0_m)
residuals <- residuals(fit_go_nr_0_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")


formula <- bf(mean_go_index | weights(go_index_weight) ~ me(log_visitors, sd_log_visitors) + (1|species) + (1|elevation))
fit_go_nr_l_m <- brm(
  formula = formula,
  data = go.final.table,
  family = zero_inflated_beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 6000,
  warmup = 2000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
saveRDS(fit_go_nr_l_m, file = "../brms_models/fit_go_nr_l_m.rds")
fit_go_nr_l_m <- readRDS("../brms_models/fit_go_nr_l_m.rds")
summary(fit_go_nr_l_m)
pp_check(fit_go_nr_l_m, type = "dens_overlay")
pp_check(fit_go_nr_l_m, type = "stat")
bayes_R2(fit_go_nr_l_m)
plot(fit_go_nr_l_m)
fitted_values <- fitted(fit_go_nr_l_m)
residuals <- residuals(fit_go_nr_l_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")



formula <- bf(mean_go_index | weights(go_index_weight) ~ me(mean.visitors.scaled, sd.visitors.scaled) + (1|species) + (1|elevation))
fit_go_r_o_m <- brm(
  formula = formula,
  data = go.final.table,
  family = zero_inflated_beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 6000,
  warmup = 2000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
saveRDS(fit_go_r_o_m, file = "../brms_models/fit_go_r_o_m.rds")
fit_go_r_o_m <- readRDS("../brms_models/fit_go_r_o_m.rds")
summary(fit_go_r_o_m)
pp_check(fit_go_r_o_m, type = "dens_overlay")
pp_check(fit_go_r_o_m, type = "stat")
bayes_R2(fit_go_r_o_m)
plot(fit_go_r_o_m)
fitted_values <- fitted(fit_go_r_o_m)
residuals <- residuals(fit_go_r_o_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")



formula <- bf(mean_go_index | weights(go_index_weight) ~ me(log.m.s.p.o, log.s.s.p.o) + (1|species) + (1|elevation))
fit_go_r2_o_m <- brm(
  formula = formula,
  data = go.final.table,
  family = zero_inflated_beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 6000,
  warmup = 2000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
saveRDS(fit_go_r2_o_m, file = "../brms_models/fit_go_r2_o_m.rds")
fit_go_r2_o_m <- readRDS("../brms_models/fit_go_r2_o_m.rds")
summary(fit_go_r2_o_m)
pp_check(fit_go_r2_o_m, type = "dens_overlay")
pp_check(fit_go_r2_o_m, type = "stat")
bayes_R2(fit_go_r2_o_m)
plot(fit_go_r2_o_m)
fitted_values <- fitted(fit_go_r2_o_m)
residuals <- residuals(fit_go_r2_o_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")





















###* GEITONOGAMY PARTIAL




















formula <- bf(mean_go_index | weights(go_index_weight) ~ me(mean.visitors.per.minute, sd.visitors.per.minute) + (1|species) + (1|elevation))
fit_go_nr_p_m <- brm(
  formula = formula,
  data = go.final.table,
  family = zero_one_inflated_beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 8,
  cores = 8,
  iter = 4000,
  warmup = 1000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
saveRDS(fit_go_nr_p_m, file = "../brms_models/fit_go_nr_p_m.rds")
fit_go_nr_p_m <- readRDS("../brms_models/fit_go_nr_p_m.rds")
summary(fit_go_nr_p_m)
pp_check(fit_go_nr_p_m, type = "dens_overlay")
pp_check(fit_go_nr_p_m, type = "stat")
bayes_R2(fit_go_nr_p_m)
plot(fit_go_nr_p_m)
fitted_values <- fitted(fit_go_nr_p_m)
residuals <- residuals(fit_go_nr_p_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")


formula <- bf(mean_go_index | weights(go_index_weight) ~ me(log_visitors, sd_log_visitors) + (1|species) + (1|elevation))
fit_go_nr_p_m <- brm(
  formula = formula,
  data = go.final.table,
  family = zero_inflated_beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 6000,
  warmup = 2000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
saveRDS(fit_go_nr_p_m, file = "../brms_models/fit_go_nr_p_m.rds")
fit_go_nr_p_m <- readRDS("../brms_models/fit_go_nr_p_m.rds")
summary(fit_go_nr_p_m)
pp_check(fit_go_nr_p_m, type = "dens_overlay")
pp_check(fit_go_nr_p_m, type = "stat")
bayes_R2(fit_go_nr_p_m)
plot(fit_go_nr_p_m)
fitted_values <- fitted(fit_go_nr_p_m)
residuals <- residuals(fit_go_nr_p_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")



formula <- bf(mean_go_index | weights(go_index_weight) ~ me(mean.visitors.scaled, sd.visitors.scaled) + (1|species) + (1|elevation))
fit_go_r_p_m <- brm(
  formula = formula,
  data = go.final.table,
  family = zero_inflated_beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 6000,
  warmup = 2000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
saveRDS(fit_go_r_p_m, file = "../brms_models/fit_go_r_p_m.rds")
fit_go_r_p_m <- readRDS("../brms_models/fit_go_r_p_m.rds")
summary(fit_go_r_p_m)
pp_check(fit_go_r_p_m, type = "dens_overlay")
pp_check(fit_go_r_p_m, type = "stat")
bayes_R2(fit_go_r_p_m)
plot(fit_go_r_p_m)
fitted_values <- fitted(fit_go_r_p_m)
residuals <- residuals(fit_go_r_p_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")



formula <- bf(mean_go_index | weights(go_index_weight) ~ me(log.m.s.p.o, log.s.s.p.o) + (1|species) + (1|elevation))
fit_go_r2_p_m <- brm(
  formula = formula,
  data = go.final.table,
  family = zero_inflated_beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 6000,
  warmup = 2000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
saveRDS(fit_go_r2_p_m, file = "../brms_models/fit_go_r2_p_m.rds")
fit_go_r2_p_m <- readRDS("../brms_models/fit_go_r2_p_m.rds")
summary(fit_go_r2_p_m)
pp_check(fit_go_r2_p_m, type = "dens_overlay")
pp_check(fit_go_r2_p_m, type = "stat")
bayes_R2(fit_go_r2_p_m)
plot(fit_go_r2_p_m)
fitted_values <- fitted(fit_go_r2_p_m)
residuals <- residuals(fit_go_r2_p_m)
plot(fitted_values, residuals)
abline(h = 0, col = "red")


