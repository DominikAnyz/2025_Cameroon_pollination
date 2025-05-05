source("scripts/setup.R")
source("scripts/Bayesian setup.R")

###* Many options have to be tested:
###* FOR SEEDSET, fit models:
###* assuming a normal distribution. Maybe thing about zero-inflated data? Or any 
###* other distribution?
###*
###* - not rescaled
###* --- with mean seedset = 0 (will use poisson and rounded seedset) ---> c.pl.final.table.2
###* ----- formula <- bf(mean_seedset_round | weights (seedset_weight) ~ me(log_visitors, sd_log_visitors) + (1|species) + (1|elevation))
###*
###* - rescaled with scale function (scale)
###* --- with mean seedset = 0 will use poisson and rounded seedset) ---> c.pl.final.table.2
###* ----- formula <- bf(mean_seedset_round | weights (seedset_weight) ~ me(mean.visitors.scaled, sd.visitors.sclaed) + (1|species) + (1|elevation))
###*
###* - rescaled with scale function (log + 1)
###* --- with mean seedset = 0 (will use poisson and rounded seedset) ---> c.pl.final.table.2
###* ----- formula <- bf(mean_seedset_round | weights (seedset_weight) ~ me(log.m.s.p.o, log.s.s.p.o) + (1|species) + (1|elevation))
###* 
###* All the results from these should be compared to find the best fitting model
###* Then we can talk about the option of doing ordbeta, but I do not think that
###* it is possible with if just one extreme is present
###* 
c.pl.final.table$elevation <- as.factor(c.pl.final.table$elevation)
c.pl.final.table$species <- as.factor(c.pl.final.table$species)
c.pl.final.table <- na.omit(c.pl.final.table)

# c.pl.final.table.2 <- c.pl.final.table %>%
#   mutate(mean_seedset = ifelse(mean_seedset == 0, 0.001, mean_seedset),
#          mean_PL_index = ifelse(mean_PL_index == 1, 0.999, mean_PL_index))

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

c.pl.final.table.2 <- c.pl.final.table %>%
  mutate(log_mean_seedset = log(mean_seedset),
         se_log_mean_seedset = se_seedset / mean_seedset,
         log_visitors = log(mean.visitors.per.minute),
         sd_log_visitors = sd.visitors.per.minute / mean.visitors.per.minute,
         mean.scaled.plus.one = mean.visitors.scaled + 1,
         sd.scaled.plus.one = sd.visitors.scaled + 1,
         log.m.s.p.o = log(mean.scaled.plus.one),
         log.s.s.p.o = log(sd.scaled.plus.one)
  )

c.pl.final.table <- c.pl.final.table %>%
  mutate(mean_seedset_round = round(mean_seedset))

c.pl.final.table.2 <- c.pl.final.table.2 %>%
  mutate(mean_seedset_round = round(mean_seedset))

c.pl.final.table <- na.omit(c.pl.final.table)

hist(c.pl.final.table$mean_seedset)

View(c.pl.final.table)
View(c.pl.final.table.2)
###*
###*
###*
###*
###* - not rescaled
###* --- with mean_seedset = 0 ---> use poisson
formula <- bf(mean_seedset_round | weights (n_replicates) ~ me(log_visitors, sd_log_visitors) + (1|species) + (1|elevation))
###* Setting priors
prior <- c(set_prior("normal(0, 2.5)", class = "b"), set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"), set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))

fit_seedset_nr_p <- brm(
  formula = formula,
  data = c.pl.final.table.2,
  chains = 10,
  cores = 10,
  iter = 4000,
  warmup = 1000,
  family = poisson(),
  prior = prior,
  save_pars = save_pars(all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234 # Addition to code
)

saveRDS(fit_seedset_nr_p, file = "../brms_models/fit_seedset_nr_p.rds")
fit_seedset_nr_p <- readRDS("../brms_models/fit_seedset_nr_p.rds")
summary(fit_seedset_nr_p)
pp_check(fit_seedset_nr_p)
bayes_R2(fit_seedset_nr_p)
plot(fit_seedset_nr_p)
fitted_values <- fitted(fit_seedset_nr_p)
residuals <- residuals(fit_seedset_nr_p)
plot(fitted_values, residuals)
abline(h = 0, col = "red")
###*
###*
###*
###*
###*
###*
###* - not rescaled
###* --- with mean_seedset = 0 ---> use poisson
formula <- bf(mean_seedset_round | weights (n_replicates) ~ me(mean.visitors.per.minute, sd.visitors.per.minute) + (1|species) + (1|elevation))
###* Setting priors
prior <- c(set_prior("normal(0, 2.5)", class = "b"), set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"), set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))

fit_seedset_nr_p_n <- brm(
  formula = formula,
  data = c.pl.final.table.2,
  chains = 10,
  cores = 10,
  iter = 4000,
  warmup = 1000,
  family = poisson(),
  prior = prior,
  save_pars = save_pars(all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234 # Addition to code
)

saveRDS(fit_seedset_nr_p_n, file = "../brms_models/fit_seedset_nr_p_n.rds")
fit_seedset_nr_p_n <- readRDS("../brms_models/fit_seedset_nr_p_n.rds")
summary(fit_seedset_nr_p_n)
pp_check(fit_seedset_nr_p_n)
bayes_R2(fit_seedset_nr_p_n)
plot(fit_seedset_nr_p_n)
fitted_values <- fitted(fit_seedset_nr_p_n)
residuals <- residuals(fit_seedset_nr_p_n)
plot(fitted_values, residuals)
abline(h = 0, col = "red")
###*
###*
###*
###*
###*
###*
###* - rescaled with scale function (scale)
###* --- with mean seedset = 0 (will be changed to 0.001, se will stay 0) ---> c.pl.final.table.2
formula <- bf(mean_seedset_round | weights (n_replicates) ~ me(mean.visitors.scaled, sd.visitors.scaled) + (1|species) + (1|elevation))
###* Setting priors
prior <- c(set_prior("normal(0, 2.5)", class = "b"), set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"), set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))

fit_seedset_r_p_n <- brm(
  formula = formula,
  data = c.pl.final.table.2,
  chains = 10,
  cores = 10,
  iter = 4000,
  warmup = 1000,
  family = poisson(),
  prior = prior,
  save_pars = save_pars(all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234 # Addition to code
)

saveRDS(fit_seedset_r_p, file = "../brms_models/fit_seedset_r_p.rds")
fit_seedset_r_p <- readRDS("../brms_models/fit_seedset_r_p.rds")
summary(fit_seedset_r_p)
pp_check(fit_seedset_r_p)
bayes_R2(fit_seedset_r_p)
plot(fit_seedset_r_p)
fitted_values <- fitted(fit_seedset_r_p)
residuals <- residuals(fit_seedset_r_p)
plot(fitted_values, residuals)
abline(h = 0, col = "red")
###*
###*
###*
###*
###*
###*
###*
###* - rescaled with scale function (log + 1)
###* --- with mean seedset = 0 (will be changed to 0.001, se will stay 0) ---> c.pl.final.table.2
formula <- bf(mean_seedset_round | weights (n_replicates) ~ me(log.m.s.p.o, log.s.s.p.o) + (1|species) + (1|elevation))
###* Setting priors
prior <- c(set_prior("normal(0, 2.5)", class = "b"), set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"), set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))

fit_seedset_r2_p_n <- brm(
  formula = formula,
  data = c.pl.final.table.2,
  chains = 10,
  cores = 10,
  iter = 6000,
  warmup = 2000,
  family = poisson(),
  prior = prior,
  save_pars = save_pars(all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234 # Addition to code
)

saveRDS(fit_seedset_r2_p_n, file = "../brms_models/fit_seedset_r2_p_n.rds")
fit_seedset_r2_p_n <- readRDS("../brms_models/fit_seedset_r2_p_n.rds")
summary(fit_seedset_r2_p_n)
pp_check(fit_seedset_r2_p_n)
bayes_R2(fit_seedset_r2_p_n)
plot(fit_seedset_r2_p_n)
fitted_values <- fitted(fit_seedset_r2_p_n)
residuals <- residuals(fit_seedset_r2_p_n)
plot(fitted_values, residuals)
abline(h = 0, col = "red")
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
###* - not rescaled, loged
###* --- with mean_seedset = 0 ---> use poisson
formula <- bf(mean_seedset_round | weights (n_replicates) ~ mean.visitors.per.minute + (1|species) + (1|elevation))
###* Setting priors
prior <- c(set_prior("normal(0, 2.5)", class = "b"), set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"), set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))

fit_seedset_nr_0_partial_o <- brm(
  formula = formula,
  data = c.pl.final.table.2,
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

saveRDS(fit_seedset_nr_p_partial_o, file = "../brms_models/fit_seedset_nr_p_partial_o.rds")
fit_seedset_nr_p_partial_o <- readRDS("../brms_models/fit_seedset_nr_p_partial_o.rds")
summary(fit_seedset_nr_p_partial_o)
pp_check(fit_seedset_nr_p_partial_o)
bayes_R2(fit_seedset_nr_p_partial_o)
plot(fit_seedset_nr_p_partial_o)
fitted_values <- fitted(fit_seedset_nr_p_partial_o)
residuals <- residuals(fit_seedset_nr_p_partial_o)
plot(fitted_values, residuals)
abline(h = 0, col = "red")
###* 
###* 
###* 
###*
###*
###*
###*
###* - not rescaled, loged
###* --- with mean_seedset = 0 ---> use poisson
formula <- bf(mean_seedset_round | weights (n_replicates) ~ log_visitors + (1|species) + (1|elevation))
###* Setting priors
prior <- c(set_prior("normal(0, 2.5)", class = "b"), set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"), set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))

fit_seedset_nr_0_partial_n <- brm(
  formula = formula,
  data = c.pl.final.table.2,
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

saveRDS(fit_seedset_nr_p_partial_n, file = "../brms_models/fit_seedset_nr_p_partial_n.rds")
fit_seedset_nr_p_partial_n <- readRDS("../brms_models/fit_seedset_nr_p_partial_n.rds")
summary(fit_seedset_nr_p_partial_n)
pp_check(fit_seedset_nr_p_partial_n)
bayes_R2(fit_seedset_nr_p_partial_n)
plot(fit_seedset_nr_p_partial_n)
fitted_values <- fitted(fit_seedset_nr_p_partial_n)
residuals <- residuals(fit_seedset_nr_p_partial_n)
plot(fitted_values, residuals)
abline(h = 0, col = "red")
###*
###*
###*
###*
###*
###* - rescaled with scale function (scale)
###* --- with mean seedset = 0 (will be changed to 0.001, se will stay 0) ---> c.pl.final.table.2
formula <- bf(mean_seedset_round | weights (n_replicates) ~ mean.visitors.scaled + (1|species) + (1|elevation))
###* Setting priors
prior <- c(set_prior("normal(0, 2.5)", class = "b"), set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"), set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))

fit_seedset_r_p_partial_n <- brm(
  formula = formula,
  data = c.pl.final.table.2,
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

saveRDS(fit_seedset_r_p_partial_n, file = "../brms_models/fit_seedset_r_p_partial_n.rds")
fit_seedset_r_p_partial_n <- readRDS("../brms_models/fit_seedset_r_p_partial_n.rds")
summary(fit_seedset_r_p_partial_n)
pp_check(fit_seedset_r_p_partial_n)
bayes_R2(fit_seedset_r_p_partial_n)
plot(fit_seedset_r_p_partial_n)
fitted_values <- fitted(fit_seedset_r_p_partial_n)
residuals <- residuals(fit_seedset_r_p_partial_n)
plot(fitted_values, residuals)
abline(h = 0, col = "red")
###*
###*
###*
###*
###*
###*
###*
###* - rescaled with scale function (log + 1)
###* --- with mean seedset = 0 (will be changed to 0.001, se will stay 0) ---> c.pl.final.table.2
formula <- bf(mean_seedset_round | weights (n_replicates) ~ log.m.s.p.o + (1|species) + (1|elevation))
###* Setting priors
prior <- c(set_prior("normal(0, 2.5)", class = "b"), set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"), set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))

fit_seedset_r2_p_partial_n <- brm(
  formula = formula,
  data = c.pl.final.table.2,
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

saveRDS(fit_seedset_r2_p_partial_n, file = "../brms_models/fit_seedset_r2_p_partial_n.rds")
fit_seedset_r2_p_partial_n <- readRDS("../brms_models/fit_seedset_r2_p_partial_n.rds")
summary(fit_seedset_r2_p_partial_n)
pp_check(fit_seedset_r2_p_partial_n)
bayes_R2(fit_seedset_r2_p_partial_n)
plot(fit_seedset_r2_p_partial_n)
fitted_values <- fitted(fit_seedset_r2_p_partial_n)
residuals <- residuals(fit_seedset_r2_p_partial_n)
plot(fitted_values, residuals)
abline(h = 0, col = "red")
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
###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* FOR POLLEN LIMITATION, fit models:
###* all with new weights. Instead of number of replicates, using PL_index_weights
###* 
###* - not rescaled
###* --- without mean_PL_index = 1 ---> c.pl.final.table
###* ----- use actual values in visitation
###* ------- formula <- bf(mean_PL_index | weights(PL_index_weights) ~ me(mean.visitors.per.minute, sd.visitors.per.minute) + (1|species) + (1|elevation))
###* ----- use logged values in visitation
###* ------- formula <- bf(mean_PL_index | weights(PL_index_weights) ~ me(log_visitors, sd_log_visitors) + (1|species) + (1|elevation))
###* --- with mean_PL_index = 1 (will be changed to 0.9999) ---> c.pl.final.table.2
###* ----- use actual values in visitation
###* ------- formula <- bf(mean_PL_index | weights(PL_index_weights) ~ me(mean.visitors.per.minute, sd.visitors.per.minute) + (1|species) + (1|elevation))
###* ----- use logged values in visitation 
###* ------- formula <- bf(mean_PL_index | weights(PL_index_weights) ~ me(log_visitors, sd_log_visitors) + (1|species) + (1|elevation))
###* 
###* - rescaled with scale function (scale)
###* --- without mean_PL_index = 1 ---> c.pl.final.table
###* ----- formula <- bf(mean_PL_index | weights(PL_index_weights) ~ me(mean.visitors.scaled, sd.visitors.sclaed) + (1|species) + (1|elevation))
###* --- with mean_PL_index = 1 (will be changed to 0.9999) ---> c.pl.final.table.2
###* ----- formula <- bf(mean_PL_index | weights(PL_index_weights) ~ me(mean.visitors.scaled, sd.visitors.sclaed) + (1|species) + (1|elevation))
###* 
###* - rescaled with scale function (log + 1)
###* --- without mean_PL_index = 1 ---> c.pl.final.table
###* ----- formula <- bf(mean_PL_index | weights(PL_index_weights) ~ me(log.m.s.p.o, log.s.s.p.o) + (1|species) + (1|elevation))
###* --- with mean_PL_index = 1 (will be changed to 0.9999) ---> c.pl.final.table.2
###* ----- formula <- bf(mean_PL_index | weights(PL_index_weights) ~ me(log.m.s.p.o, log.s.s.p.o) + (1|species) + (1|elevation))
###* 
###* VARIATIONS 
###* - not rescaled
###* --- without mean_PL_index = 1 ---> c.pl.final.table
###* ----- use actual values in visitation
formula <- bf(mean_PL_index | weights(PL_index_weight) ~ me(mean.visitors.per.minute, sd.visitors.per.minute) + (1|species) + (1|elevation))

fit_PL_nr_1_o <- brm(
  formula = formula,
  data = c.pl.final.table,
  family = Beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 8000,
  warmup = 2000,
  save_pars = save_pars(all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
summary(fit_PL_nr_1_o)
saveRDS(fit_PL_nr_1_o, file = "../brms_models/fit_PL_nr_1_o.rds")
fit_PL_nr_1_o <- readRDS("../brms_models/fit_PL_nr_1_o.rds")
summary(fit_seedset_r2_p_partial)
pp_check(fit_seedset_r2_p_partial, type = "dens_overlay")
pp_check(fit_seedset_r2_p_partial, type = "stat")
bayes_R2(fit_seedset_r2_p_partial)
plot(fit_seedset_r2_p_partial)
fitted_values <- fitted(fit_seedset_r2_p_partial)
residuals <- residuals(fit_seedset_r2_p_partial)
plot(fitted_values, residuals)
abline(h = 0, col = "red")
###*
###*
###*
###*
###* - not rescaled
###* --- with mean_PL_index = 1 ---> c.pl.final.table 
###* ----- use logged values in visitation
formula <- bf(mean_PL_index | weights(PL_index_weight) ~ me(log_visitors, sd_log_visitors) + (1|species) + (1|elevation)
)

fit_PL_nr_1_l <- brm(
  formula = formula,
  data = c.pl.final.table,
  family = Beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 4000,
  warmup = 1000,
  save_pars = save_pars(all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
summary(fit_PL_nr_1_l)
saveRDS(fit_PL_nr_1_l, file = "../brms_models/fit_PL_nr_1_l.rds")
fit_PL_nr_1_l <- readRDS("../brms_models/fit_PL_nr_1_l.rds")
summary(fit_PL_nr_1_l)
pp_check(fit_PL_nr_1_l, type = "dens_overlay")
pp_check(fit_PL_nr_1_l, type = "stat")
bayes_R2(fit_PL_nr_1_l)
plot(fit_PL_nr_1_l)
fitted_values <- fitted(fit_PL_nr_1_l)
residuals <- residuals(fit_PL_nr_1_l)
plot(fitted_values, residuals)
abline(h = 0, col = "red")
###* 
###* 
###* 
###* 
###* - not rescaled
###* --- with mean_PL_index = 1 (will use zero-one-inflated-beta) ---> c.pl.final.table.2
###* ----- use actual values in visitation
formula <- bf(mean_PL_index | weights(PL_index_weight) ~ me(mean.visitors.per.minute, sd.visitors.per.minute) + (1|species) + (1|elevation)
)
fit_PL_nr_y_o <- brm(
  formula = formula,
  data = c.pl.final.table.2,
  family = zero_one_inflated_beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 4000,
  warmup = 1000,
  save_pars = save_pars(all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
summary(fit_PL_nr_y_o)
saveRDS(fit_PL_nr_y_o, file = "../brms_models/fit_PL_nr_y_o.rds")
fit_PL_nr_y_o <- readRDS("../brms_models/fit_PL_nr_y_o.rds")
summary(fit_PL_nr_y_o)
pp_check(fit_PL_nr_y_o, type = "dens_overlay")
pp_check(fit_PL_nr_y_o, type = "stat")
bayes_R2(fit_PL_nr_y_o)
plot(fit_PL_nr_y_o)
fitted_values <- fitted(fit_PL_nr_y_o)
residuals <- residuals(fit_PL_nr_y_o)
plot(fitted_values, residuals)
abline(h = 0, col = "red")
###*
###*
###*
###*
###*
###* - not rescaled
###* --- with mean_PL_index = 1 (will use zero-one-inflated-beta) ---> c.pl.final.table.2
###* ----- use loged values in visitation
formula <- bf(mean_PL_index | weights(PL_index_weight) ~ me(log_visitors, sd_log_visitors) + (1|species) + (1|elevation)
)

fit_PL_nr_y_l <- brm(
  formula = formula,
  data = c.pl.final.table.2,
  family = zero_one_inflated_beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 4000,
  warmup = 1000,
  save_pars = save_pars(all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
summary(fit_PL_nr_y_l)
saveRDS(fit_PL_nr_y_l, file = "../brms_models/fit_PL_nr_y_l.rds")
fit_PL_nr_y_l <- readRDS("../brms_models/fit_PL_nr_y_l.rds")
summary(fit_PL_nr_1_l)
pp_check(fit_PL_nr_1_l, type = "dens_overlay")
pp_check(fit_PL_nr_1_l, type = "stat")
bayes_R2(fit_PL_nr_1_l)
plot(fit_PL_nr_1_l)
fitted_values <- fitted(fit_PL_nr_1_l)
residuals <- residuals(fit_PL_nr_1_l)
plot(fitted_values, residuals)
abline(h = 0, col = "red")
###* 
###* 
###* 
###* 
###* 
###* - rescaled with scale function (scale)
###* --- without mean_PL_index = 1 ---> c.pl.final.table
formula <- bf(mean_PL_index | weights(PL_index_weight) ~ me(mean.visitors.scaled, sd.visitors.scaled) + (1|species) + (1|elevation))

fit_PL_r_1 <- brm(
  formula = formula,
  data = c.pl.final.table,
  family = Beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 4000,
  warmup = 1000,
  save_pars = save_pars(all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
summary(fit_PL_r_1)
saveRDS(fit_PL_r_1, file = "../brms_models/fit_PL_r_1.rds")
fit_PL_r_1 <- readRDS("../brms_models/fit_PL_r_1.rds")
pp_check(fit_PL_r_1, type = "dens_overlay")
pp_check(fit_PL_r_1, type = "stat")
bayes_R2(fit_PL_r_1)
plot(fit_PL_r_1)
fitted_values <- fitted(fit_PL_r_1)
residuals <- residuals(fit_PL_r_1)
plot(fitted_values, residuals)
abline(h = 0, col = "red")
###*
###*
###*
###*
###*
###* - rescaled with scale function (scale)
###* --- with mean_PL_index = 1 (will use zero-one-inflated-beta) ---> c.pl.final.table.2
formula <- bf(mean_PL_index | weights(PL_index_weight) ~ me(mean.visitors.scaled, sd.visitors.scaled) + (1|species) + (1|elevation))

fit_PL_r_y <- brm(
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
  iter = 4000,
  warmup = 1000,
  save_pars = save_pars(all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
summary(fit_PL_r_y)
saveRDS(fit_PL_r_y, file = "../brms_models/fit_PL_r_y.rds")
#fit_PL_r_y <- readRDS("../brms_models/fit_PL_r_y.rds")
###* 
###* 
###* 
###* 
###*
###*
###* - rescaled with scale function (log + 1)
###* --- without mean_PL_index = 1 ---> c.pl.final.table
formula <- bf(mean_PL_index | weights(PL_index_weight) ~ me(log.m.s.p.o, log.s.s.p.o) + (1|species) + (1|elevation))

fit_PL_r2_1 <- brm(
  formula = formula,
  data = c.pl.final.table,
  family = Beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 4000,
  warmup = 1000,
  save_pars = save_pars(all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
summary(fit_PL_r2_1)
saveRDS(fit_PL_r2_1, file = "../brms_models/fit_PL_r2_1.rds")
#fit_PL_r2_1 <- readRDS("../brms_models/fit_PL_r2_1.rds")
###*
###*
###* - rescaled with scale function (log + 1)
###* --- with mean_PL_index = 1 (will be changed to 0.9999) ---> c.pl.final.table.2
formula <- bf(mean_PL_index | weights(PL_index_weight) ~ me(log.m.s.p.o, log.s.s.p.o) + (1|species) + (1|elevation))

fit_PL_r2_y <- brm(
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
  iter = 4000,
  warmup = 1000,
  save_pars = save_pars(all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
summary(fit_PL_r2_y)
saveRDS(fit_PL_r2_y, file = "../brms_models/fit_PL_r2_y.rds")
#fit_PL_r2_y <- readRDS("../brms_models/fit_PL_r2_y.rds")
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
###* VARIATIONS 
###* - not rescaled
###* --- without mean_PL_index = 1 ---> c.pl.final.table
###* ----- use actual values in visitation
formula <- bf(mean_PL_index | weights(PL_index_weight) ~ mean.visitors.per.minute + (1|species) + (1|elevation))

fit_PL_nr_1_o_partial <- brm(
  formula = formula,
  data = c.pl.final.table,
  family = Beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 8000,
  warmup = 2000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
summary(fit_PL_nr_1_o_partial)
saveRDS(fit_PL_nr_1_o_partial, file = "../brms_models/fit_PL_nr_1_o_partial.rds")
#fit_PL_nr_1_o_partial <- readRDS("../brms_models/fit_PL_nr_1_o_partial.rds")
###*
###*
###*
###*
###* - not rescaled
###* --- with mean_PL_index = 1 ---> c.pl.final.table 
###* ----- use logged values in visitation
formula <- bf(mean_PL_index | weights(PL_index_weight) ~ log_visitors + (1|species) + (1|elevation)
)

fit_PL_nr_1_l_partial <- brm(
  formula = formula,
  data = c.pl.final.table,
  family = Beta(),
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
summary(fit_PL_nr_1_l_partial)
saveRDS(fit_PL_nr_1_l_partial, file = "../brms_models/fit_PL_nr_1_l_partial.rds")
#fit_PL_nr_1_l_partial <- readRDS("../brms_models/fit_PL_nr_1_l_partial.rds")
###* 
###* 
###* 
###* 
###* - not rescaled
###* --- with mean_PL_index = 1 (will use zero-one-inflated-beta) ---> c.pl.final.table.2
###* ----- use actual values in visitation
formula <- bf(mean_PL_index | weights(PL_index_weight) ~ mean.visitors.per.minute + (1|species) + (1|elevation)
)
fit_PL_nr_y_o_partial <- brm(
  formula = formula,
  data = c.pl.final.table.2,
  family = zero_one_inflated_beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 4000,
  warmup = 1000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
summary(fit_PL_nr_y_o_partial)
saveRDS(fit_PL_nr_y_o_partial, file = "../brms_models/fit_PL_nr_y_o_partial.rds")
#fit_PL_nr_y_o_partial <- readRDS("../brms_models/fit_PL_nr_y_o_partial.rds")
###*
###*
###*
###*
###*
###* - not rescaled
###* --- with mean_PL_index = 1 (will use zero-one-inflated-beta) ---> c.pl.final.table.2
###* ----- use loged values in visitation
formula <- bf(mean_PL_index | weights(PL_index_weight) ~ log_visitors + (1|species) + (1|elevation)
)

fit_PL_nr_y_l_partial <- brm(
  formula = formula,
  data = c.pl.final.table.2,
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
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234
)
summary(fit_PL_nr_y_l_partial)
saveRDS(fit_PL_nr_y_l_partial, file = "../brms_models/fit_PL_nr_y_l_partial.rds")
#fit_PL_nr_y_l_partial <- readRDS("../brms_models/fit_PL_nr_y_l_partial.rds")
###* 
###* 
###* 
###* 
###* 
###* - rescaled with scale function (scale)
###* --- without mean_PL_index = 1 ---> c.pl.final.table
formula <- bf(mean_PL_index | weights(PL_index_weight) ~ mean.visitors.scaled + (1|species) + (1|elevation))

fit_PL_r_1_partial <- brm(
  formula = formula,
  data = c.pl.final.table,
  family = Beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 4000,
  warmup = 1000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
summary(fit_PL_r_1_partial)
saveRDS(fit_PL_r_1_partial, file = "../brms_models/fit_PL_r_1_partial.rds")
#fit_PL_r_1_partial <- readRDS("../brms_models/fit_PL_r_1_partial.rds")
###*
###*
###*
###*
###*
###* - rescaled with scale function (scale)
###* --- with mean_PL_index = 1 (will use zero-one-inflated-beta) ---> c.pl.final.table.2
formula <- bf(mean_PL_index | weights(PL_index_weight) ~ mean.visitors.scaled + (1|species) + (1|elevation))

fit_PL_r_y_partial <- brm(
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
  iter = 4000,
  warmup = 1000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
summary(fit_PL_r_y_partial)
saveRDS(fit_PL_r_y_partial, file = "../brms_models/fit_PL_r_y_partial.rds")
#fit_PL_r_y_partial <- readRDS("../brms_models/fit_PL_r_y_partial.rds")
###* 
###* 
###* 
###* 
###*
###*
###* - rescaled with scale function (log + 1)
###* --- without mean_PL_index = 1 ---> c.pl.final.table
formula <- bf(mean_PL_index | weights(PL_index_weight) ~ log.m.s.p.o + (1|species) + (1|elevation))

fit_PL_r2_1_partial <- brm(
  formula = formula,
  data = c.pl.final.table,
  family = Beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 4000,
  warmup = 1000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
summary(fit_PL_r2_1_partial)
saveRDS(fit_PL_r2_1_partial, file = "../brms_models/fit_PL_r2_1_partial.rds")
#fit_PL_r2_1_partial <- readRDS("../brms_models/fit_PL_r2_1_partial.rds")
###*
###*
###* - rescaled with scale function (log + 1)
###* --- with mean_PL_index = 1 (will be changed to 0.9999) ---> c.pl.final.table.2
formula <- bf(mean_PL_index | weights(PL_index_weight) ~ log.m.s.p.o + (1|species) + (1|elevation))

fit_PL_r2_y_partial <- brm(
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
  iter = 4000,
  warmup = 1000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
summary(fit_PL_r2_y_partial)
saveRDS(fit_PL_r2_y_partial, file = "../brms_models/fit_PL_r2_y_partial.rds")
#fit_PL_r2_y_partial <- readRDS("../brms_models/fit_PL_r2_y_partial.rds")
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
###* 
###* 
###* FOR AUTOGAMY, fit models:
###* all with new weights. Instead of number of replicates, using AO_index_weight
###* Also, we cannot even try doing 
###* 
###* - not rescaled
###* --- with mean_ao_index = 0 (will use zero-inflated-beta) 
###* ----- use actual values in visitation
###* ------- formula <- bf(mean_ao_index | weights(ao_index_weight) ~ me(mean.visitors.per.minute, sd.visitors.per.minute) + (1|species) + (1|elevation))
###* ----- use logged values in visitation 
###* ------- formula <- bf(mean_ao_index | weights(ao_index_weight) ~ me(log_visitors, sd_log_visitors) + (1|species) + (1|elevation))
###* 
###* - rescaled with scale function (scale)
###* --- with mean_ao_index = 0 (will use zero-inflated-beta) 
###* ----- formula <- bf(mean_PL_index | weights(PL_index_weight) ~ me(mean.visitors.scaled, sd.visitors.sclaed) + (1|species) + (1|elevation))
###* 
###* - rescaled with scale function (log + 1)
###* --- with mean_ao_index = 0 (will use zero-inflated-beta)
###* ----- formula <- bf(mean_ao_index | weights(ao_index_weight) ~ me(log.m.s.p.o, log.s.s.p.o) + (1|species) + (1|elevation))
View(ao.final.table)

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




###* 
###* 
###*
###*
###*
###* - not rescaled
###* --- with mean_ao_index = 0 (will use zero-inflated-beta) 
###* ----- use actual values in visitation
formula <- bf(mean_ao_index | weights(ao_index_weight) ~ me(mean.visitors.per.minute, sd.visitors.per.minute) + (1|species) + (1|elevation))

fit_ao_nr_0 <- brm(
  formula = formula,
  data = ao.final.table,
  family = zero_inflated_beta(),
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
saveRDS(fit_ao_nr_0, file = "../brms_models/fit_ao_nr_0.rds")
#fit_ao_nr_0 <- readRDS("../brms_models/fit_ao_nr_0.rds")
###*
###*
###*
###*
###*
###* - not rescaled
###* --- with mean_ao_index = 0 (will use zero-inflated-beta) 
###* ----- use logged values in visitation 
formula <- bf(mean_ao_index | weights(ao_index_weight) ~ me(log_visitors, sd_log_visitors) + (1|species) + (1|elevation))

fit_ao_nr_l <- brm(
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
  iter = 4000,
  warmup = 1000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
saveRDS(fit_ao_nr_l, file = "../brms_models/fit_ao_nr_l.rds")
#fit_ao_nr_l <- readRDS("../brms_models/fit_ao_nr_l.rds")
###*
###*
###*
###*
###*
###* - rescaled with scale function (scale)
###* --- with mean_ao_index = 0 (will use zero-inflated-beta) 
formula <- bf(mean_ao_index | weights(ao_index_weight) ~ me(mean.visitors.scaled, sd.visitors.scaled) + (1|species) + (1|elevation))

fit_ao_r_o <- brm(
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
  iter = 4000,
  warmup = 1000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
saveRDS(fit_ao_r_o, file = "../brms_models/fit_ao_r_o.rds")
#fit_ao_r_o <- readRDS("../brms_models/fit_ao_r_o.rds")
###*
###*
###*
###*
###*
###* - rescaled with scale function (log + 1)
###* --- with mean_ao_index = 0 (will use zero-inflated-beta)
formula <- bf(mean_ao_index | weights(ao_index_weight) ~ me(log.m.s.p.o, log.s.s.p.o) + (1|species) + (1|elevation))

fit_ao_r2_o <- brm(
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
  iter = 4000,
  warmup = 1000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)

summary(fit_ao_r2_o)
saveRDS(fit_ao_r2_o, file = "../brms_models/fit_ao_r2_o.rds")
#fit_ao_r2_o <- readRDS("../brms_models/fit_ao_r2_o.rds")
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
###* - not rescaled
###* --- with mean_ao_index = 0 (will use zero-inflated-beta) 
###* ----- use actual values in visitation
formula <- bf(mean_ao_index | weights(ao_index_weight) ~ mean.visitors.per.minute + (1|species) + (1|elevation))

fit_ao_nr_0_partial <- brm(
  formula = formula,
  data = ao.final.table,
  family = zero_inflated_beta(),
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
summary(fit_ao_nr_0_partial)
saveRDS(fit_ao_nr_0_partial, file = "../brms_models/fit_ao_nr_0_partial.rds")
#fit_ao_nr_0_partial <- readRDS("../brms_models/fit_ao_nr_0_partial.rds")
###*
###*
###*
###*
###*
###* - not rescaled
###* --- with mean_ao_index = 0 (will use zero-inflated-beta) 
###* ----- use logged values in visitation 
formula <- bf(mean_ao_index | weights(ao_index_weight) ~ log_visitors + (1|species) + (1|elevation))

fit_ao_nr_l_partial <- brm(
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
summary(fit_ao_nr_l_partial)
saveRDS(fit_ao_nr_l_partial, file = "../brms_models/fit_ao_nr_l_partial.rds")
#fit_ao_nr_l_partial <- readRDS("../brms_models/fit_ao_nr_l_partial.rds")
###*
###*
###*
###*
###*
###* - rescaled with scale function (scale)
###* --- with mean_ao_index = 0 (will use zero-inflated-beta) 
formula <- bf(mean_ao_index | weights(ao_index_weight) ~ mean.visitors.scaled + (1|species) + (1|elevation))

fit_ao_r_o_partial <- brm(
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
  iter = 4000,
  warmup = 1000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
summary(fit_ao_r_o_partial)
saveRDS(fit_ao_r_o_partial, file = "../brms_models/fit_ao_r_o_partial.rds")
#fit_ao_r_o_partial <- readRDS("../brms_models/fit_ao_r_o_partial.rds")
###*
###*
###*
###*
###*
###* - rescaled with scale function (log + 1)
###* --- with mean_ao_index = 0 (will use zero-inflated-beta)
formula <- bf(mean_ao_index | weights(ao_index_weight) ~ log.m.s.p.o + (1|species) + (1|elevation))

fit_ao_r2_o_partial <- brm(
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
  iter = 4000,
  warmup = 1000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)

summary(fit_ao_r2_o_partial)
saveRDS(fit_ao_r2_o_partial, file = "../brms_models/fit_ao_r2_o_partial.rds")
#fit_ao_r2_o_partial <- readRDS("../brms_models/fit_ao_r2_o_partial.rds")
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
###*
###*
###*
###*
###* FOR GEITONOGAMY, fit models:
###* all with new weights. Instead of number of replicates, using AO_index_weight
###* Also, we cannot even try doing 
###* 
###* - not rescaled
###* --- with mean_go_index = 0 and 1 (will use zero-one-inflated-beta) 
###* ----- use actual values in visitation
###* ------- formula <- bf(mean_go_index | weights(go_index_weight) ~ me(mean.visitors.per.minute, sd.visitors.per.minute) + (1|species) + (1|elevation))
###* ----- use logged values in visitation 
###* ------- formula <- bf(mean_go_index | weights(go_index_weight) ~ me(log_visitors, sd_log_visitors) + (1|species) + (1|elevation))
###* 
###* - rescaled with scale function (scale)
###* --- with mean_go_index = 0 and 1(will use zero-one-inflated-beta) 
###* ----- formula <- bf(mean_go_index | weights(go_index_weight) ~ me(mean.visitors.scaled, sd.visitors.sclaed) + (1|species) + (1|elevation))
###* 
###* - rescaled with scale function (log + 1)
###* --- with mean_go_index = 0 and 1 (will use zero-one-inflated-beta)
###* ----- formula <- bf(mean_go_index | weights(go_index_weight) ~ me(log.m.s.p.o, log.s.s.p.o) + (1|species) + (1|elevation))
View(go.final.table)

go.final.table <- na.omit(go.final.table)
go.final.table$elevation <- as.factor(go.final.table$elevation)
go.final.table$species <- as.factor(go.final.table$species)

go.final.table <- go.final.table %>%
  mutate(
    log_visitors = log(mean.visitors.per.minute),
    sd_log_visitors = sd.visitors.per.minute / mean.visitors.per.minute,
    mean.scaled.plus.one = mean.visitors.scaled + 1,
    sd.scaled.plus.one = sd.visitors.scaled + 1,
    log.m.s.p.o = log(mean.scaled.plus.one),
    log.s.s.p.o = log(sd.scaled.plus.one)
  )

go.final.table <- go.final.table %>%
  mutate(mean_go_index = ifelse(mean_go_index > 1, 1, mean_go_index))
###*
###*
###*
###*
###*
###* - not rescaled
###* --- with mean_go_index = 0 and 1 (will use zero-one-inflated-beta) 
###* ----- use actual values in visitation
formula <- bf(mean_go_index | weights(go_index_weight) ~ me(mean.visitors.per.minute, sd.visitors.per.minute) + (1|species) + (1|elevation))

fit_go_nr_0 <- brm(
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
saveRDS(fit_go_nr_0, file = "../brms_models/fit_go_nr_0.rds")
#fit_go_nr_0 <- readRDS("../brms_models/fit_go_nr_0.rds")
###*
###*
###*
###*
###*
###* - not rescaled
###* --- with mean_go_index = 0 and 1 (will use zero-one-inflated-beta)
###* ----- use logged values in visitation 
formula <- bf(mean_go_index | weights(go_index_weight) ~ me(log_visitors, sd_log_visitors) + (1|species) + (1|elevation))

fit_go_nr_l <- brm(
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
  iter = 4000,
  warmup = 1000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
saveRDS(fit_go_nr_l, file = "../brms_models/fit_go_nr_l.rds")
#fit_go_nr_l <- readRDS("../brms_models/fit_go_nr_l.rds")
###*
###*
###*
###*
###*
###* - rescaled with scale function (scale)
###* --- with mean_go_index = 0 and 1(will use zero-one-inflated-beta) 
formula <- bf(mean_go_index | weights(go_index_weight) ~ me(mean.visitors.scaled, sd.visitors.scaled) + (1|species) + (1|elevation))

fit_go_r_o <- brm(
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
  iter = 4000,
  warmup = 1000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
saveRDS(fit_go_r_o, file = "../brms_models/fit_go_r_o.rds")
#fit_go_r_o <- readRDS("../brms_models/fit_go_r_o.rds")
###*
###*
###*
###*
###*
###* - rescaled with scale function (log + 1)
###* --- with mean_go_index = 0 and 1 (will use zero-one-inflated-beta)
formula <- bf(mean_go_index | weights(go_index_weight) ~ me(log.m.s.p.o, log.s.s.p.o) + (1|species) + (1|elevation))

fit_go_r2_o <- brm(
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
  iter = 4000,
  warmup = 1000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)
saveRDS(fit_go_r2_o, file = "../brms_models/fit_go_r2_o.rds")
#fit_go_r2_o <- readRDS("../brms_models/fit_go_r2_o.rds")
###*
###*
###*
###*
###*
###*









