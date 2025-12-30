###* Loading the necessary packages
pacman::p_load(tidyverse, glmmTMB, DHARMa, emmeans)

select <- dplyr::select

set.seed(123)
source("scripts/Bayesian setup June short.R")


###* This script is to test the effect of the predictors "total.morphospecies", 
###* "total.functional.groups" and "visited.flowers.per.minute"
###* 
###* 

#View(final.table)

final.table$elevation <- as.factor(final.table$elevation)

# Recode incorrect elevation values
final.table <- final.table %>%
  mutate(elevation = recode(elevation,
                            `3500` = "3400",
                            `4000` = "3800"))

hist(final.table$total.morpho)
hist(log(final.table$total.morpho))
hist(final.table$total.func)
hist(log(final.table$total.func))
hist(final.table$visited.flowers.per.minute)
hist(log(final.table$visited.flowers.per.minute))
hist(final.table$scaled.visited.flowers)
hist(log(final.table$scaled.visited.flowers))

library(brms)








m_brm_model_po <- brm(
  formula = bf(total.morpho ~ elevation  + (1|plant.species)),
  data = final.table,
  family = poisson(),
  prior = c(
    set_prior("normal(0, 2)", class = "b"),          # elevation effect
    set_prior("student_t(3, 0, 5)", class = "sd")   # random effects
  ),
  chains = 5, cores = 5,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  iter = 5000, warmup = 2000,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)

saveRDS(m_brm_model_po, file = "../brms_models/m_brm_model_po.rds")
m_brm_model_po <- readRDS("../brms_models/m_brm_model_po.rds")
summary(m_brm_model_po)
loo(m_brm_model, m_brm_model_po)
bayes_R2(m_brm_model)
bayes_R2(m_brm_model_po)


m_brm_model_zipo <- brm(
  formula = bf(total.morpho ~ elevation  + (1|plant.species),
               zi ~ elevation),
  data = final.table,
  family = zero_inflated_poisson(),
  prior = c(
    set_prior("normal(0, 2)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd"),
    set_prior("normal(0, 2)", class = "b", dpar = "zi")
  ),
  chains = 5, cores = 5,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  iter = 5000, warmup = 2000,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)

saveRDS(m_brm_model_zipo, file = "../brms_models/m_brm_model_zipo.rds")
m_brm_model_zipo <- readRDS("../brms_models/m_brm_model_zipo.rds")
summary(m_brm_model_zipo)
loo(m_brm_model_po, m_brm_model_zipo)
bayes_R2(m_brm_model_po)
bayes_R2(m_brm_model_zipo)


m_brm_model_nb <- brm(
  formula = bf(total.morpho ~ elevation  + (1|plant.species)),
  data = final.table,
  family = negbinomial(),
  prior = c(
    set_prior("normal(0, 2)", class = "b"),          # elevation effect
    set_prior("student_t(3, 0, 5)", class = "sd")   # random effects
  ),
  chains = 6, cores = 6,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  iter = 8000, warmup = 3000,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)

saveRDS(m_brm_model_nb, file = "../brms_models/m_brm_model_nb.rds")
m_brm_model_nb <- readRDS("../brms_models/m_brm_model_nb.rds")
summary(m_brm_model_nb)
loo(m_brm_model_pozi, m_brm_model_nb)
bayes_R2(m_brm_model_pozi)
bayes_R2(m_brm_model_nb)



m_brm_model_zinb <- brm(
  formula = bf(total.morpho ~ elevation  + (1|plant.species),
               zi ~ elevation),
  data = final.table,
  family = zero_inflated_negbinomial(),
  prior = c(
    set_prior("normal(0, 2)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd"),
    set_prior("normal(0, 2)", class = "b", dpar = "zi")
  ),
  chains = 5, cores = 5,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  iter = 5000, warmup = 2000,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)

saveRDS(m_brm_model_zinb, file = "../brms_models/m_brm_model_zinb.rds")
m_brm_model_zinb <- readRDS("../brms_models/m_brm_model_zinb.rds")
summary(m_brm_model_zinb)
loo(m_brm_model_po, m_brm_model_zipo, m_brm_model_nb, m_brm_model_zinb)
bayes_R2(m_brm_model_po)
bayes_R2(m_brm_model_zipo)
bayes_R2(m_brm_model_nb)
bayes_R2(m_brm_model_zinb)
















###* Try models to see the effect of elevation on the the number of morphospecies
mod.morpho <- glmmTMB(total.morpho ~ elevation  + (1|plant.species),
                      data = final.table,
                      family = poisson
)

mod.morpho.null <- glmmTMB(total.morpho ~ 1  + (1|plant.species),
                           data = final.table,
                           family = poisson
)


summary(mod.morpho.null)
saveRDS(mod.morpho, "glm_outputs/m_model.rds")
saveRDS(mod.morpho.null, "glm_outputs/m_null.rds")

AIC(mod.morpho, mod.morpho.null)
emm <- emmeans(mod.morpho, ~elevation)
pairs(emm)
pairs(emm, type = "response")

summary(mod.morpho)
simulationOutput3 <- simulateResiduals(fittedModel = mod.morpho, plot = F)
testDispersion(simulationOutput3)
testZeroInflation(simulationOutput3)
testResiduals(simulationOutput3)
plot(simulationOutput3)
deviance(mod.morpho)



















f_brm_model_po <- brm(
  formula = bf(total.func ~ elevation  + (1|plant.species)),
  data = final.table,
  family = poisson(),
  prior = c(
    set_prior("normal(0, 2)", class = "b"),          # elevation effect
    set_prior("student_t(3, 0, 5)", class = "sd")   # random effects
  ),
  chains = 5, cores = 5,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  iter = 5000, warmup = 2000,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)

saveRDS(f_brm_model_po, file = "../brms_models/f_brm_model_po.rds")
f_brm_model_po <- readRDS("../brms_models/f_brm_model_po.rds")
summary(f_brm_model_po)
loo(f_brm_model, f_brm_model_po)
bayes_R2(f_brm_model)
bayes_R2(f_brm_model_po)


f_brm_model_zipo <- brm(
  formula = bf(total.func ~ elevation  + (1|plant.species),
               zi ~ elevation),
  data = final.table,
  family = zero_inflated_poisson(),
  prior = c(
    set_prior("normal(0, 2)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd"),
    set_prior("normal(0, 2)", class = "b", dpar = "zi")
  ),
  chains = 5, cores = 5,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  iter = 5000, warmup = 2000,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)

saveRDS(f_brm_model_zipo, file = "../brms_models/f_brm_model_zipo.rds")
f_brm_model_zipo <- readRDS("../brms_models/f_brm_model_zipo.rds")
summary(f_brm_model_zipo)
loo(f_brm_model_po, f_brm_model_zipo)
bayes_R2(f_brm_model_po)
bayes_R2(f_brm_model_zipo)


f_brm_model_nb <- brm(
  formula = bf(total.func ~ elevation  + (1|plant.species)),
  data = final.table,
  family = negbinomial(),
  prior = c(
    set_prior("normal(0, 2)", class = "b"),          # elevation effect
    set_prior("student_t(3, 0, 5)", class = "sd")   # random effects
  ),
  chains = 3, cores = 3,
  #save_pars = save_pars(latent = TRUE, all = TRUE),
  iter = 16000, warmup = 5000,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)

saveRDS(f_brm_model_nb, file = "../brms_models/f_brm_model_nb.rds")
f_brm_model_nb <- readRDS("../brms_models/f_brm_model_nb.rds")
summary(f_brm_model_nb)
loo(f_brm_model_pozi, f_brm_model_nb)
bayes_R2(f_brm_model_pozi)
bayes_R2(f_brm_model_nb)



f_brm_model_zinb <- brm(
  formula = bf(total.func ~ elevation  + (1|plant.species),
               zi ~ elevation),
  data = final.table,
  family = zero_inflated_negbinomial(),
  prior = c(
    set_prior("normal(0, 2)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd"),
    set_prior("normal(0, 2)", class = "b", dpar = "zi")
  ),
  chains = 5, cores = 5,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  iter = 5000, warmup = 2000,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)

saveRDS(f_brm_model_zinb, file = "../brms_models/f_brm_model_zinb.rds")
f_brm_model_zinb <- readRDS("../brms_models/f_brm_model_zinb.rds")
summary(f_brm_model_zinb)
loo(f_brm_model_po, f_brm_model_zipo,  f_brm_model_zinb)
bayes_R2(f_brm_model_po)
bayes_R2(f_brm_model_zipo)
bayes_R2(f_brm_model_nb)
bayes_R2(f_brm_model_zinb)



















###* Try models to see the effect of elevation on the total number of functional groups

mod.func <- glmmTMB(total.func ~ elevation  + (1|plant.species),
                    data = final.table,
                    family = poisson
)

mod.func.null <- glmmTMB(total.func ~ 1  + (1|plant.species),
                         data = final.table,
                         family = poisson
)

summary(mod.func.null)
saveRDS(mod.func, "glm_outputs/f_model.rds")
saveRDS(mod.func.null, "glm_outputs/f_null.rds")

AIC(mod.func, mod.func.null)
emm <- emmeans(mod.func, ~elevation)
pairs(emm)
pairs(emm, type = "response")

summary(mod.func)
simulationOutput3 <- simulateResiduals(fittedModel = mod.func, plot = F)
testDispersion(simulationOutput3)
testZeroInflation(simulationOutput3)
testResiduals(simulationOutput3)
plot(simulationOutput3)
deviance(mod.morpho)























library(scales)

final.table <- final.table %>%
  ungroup() %>%  # Ensure no grouping is active
  mutate(visited.0.1.scaled = rescale(visited.flowers.per.minute, to = c(0, 1), na.rm = TRUE))







library(ordbetareg)

vis_ord_model <- ordbetareg(
  formula = visited.0.1.scaled ~ elevation  + (1|plant.species),
  data = final.table,
  chains = 5, iter = 5000, warmup = 2000, cores = 5,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)

saveRDS(vis_ord_model, file = "../brms_models/vis_ord_model.rds")
vis_ord_model <- readRDS("../brms_models/vis_ord_model.rds")


summary(vis_ord_model)


vis_ord_null <- ordbetareg(
  formula = visited.0.1.scaled ~ 1 + (1|plant.species),
  data = final.table,
  coef_prior_mean = 0,
  coef_prior_SD = 0.5,
  chains = 5, iter = 5000, warmup = 2000, cores = 5,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)


# Fit a simple model with a predictor, just to extract priors
tmp_full <- ordbetareg(
  formula = visited.0.1.scaled ~ elevation + (1|plant.species),
  data = final.table,
  chains = 1, iter = 500, warmup = 250, cores = 1
)

# Extract its default priors
priors <- default_prior(tmp_full)

# Remove the slope prior ("b"), keep only priors that exist for the null
priors <- priors[names(priors) != "b"]

# Now run the null model using the cleaned priors
vis_ord_null <- ordbetareg(
  formula = visited.0.1.scaled ~ 1 + (1|plant.species),
  data = final.table,
  chains = 5, iter = 5000, warmup = 2000, cores = 5,
  priors = priors,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)



### Check if different seed has different outcome
vis_ord_model_2 <- ordbetareg(
  formula = visited.0.1.scaled ~ elevation  + (1|plant.species),
  data = final.table,
  chains = 5, iter = 5000, warmup = 2000, cores = 5,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 123
)

saveRDS(vis_ord_model_2, file = "../brms_models/vis_ord_model_2.rds")
vis_ord_model_2 <- readRDS("../brms_models/vis_ord_model_2.rds")

summary(vis_ord_model_2)









View(final.table)

###* Try models to see the effect of elevation on the visited flowers per minute
mod.visited <- glmmTMB(visited.0.1.scaled ~ elevation  + (1|plant.species),
                       data = final.table,
                       family = ordbeta
)

mod.visited.null <- glmmTMB(visited.0.1.scaled ~ 1  + (1|plant.species),
                            data = final.table,
                            family = ordbeta
)

summary(mod.visited.null)
saveRDS(mod.visited, "glm_outputs/v_model.rds")
saveRDS(mod.visited.null, "glm_outputs/v_null.rds")


AIC(mod.visited, mod.visited.null)
emm <- emmeans(mod.visited, ~elevation)
pairs(emm)
pairs(emm, type = "response")

summary(mod.visited)
simulationOutput3 <- simulateResiduals(fittedModel = mod.visited, plot = F)
testDispersion(simulationOutput3)
testZeroInflation(simulationOutput3)
testResiduals(simulationOutput3)
plot(simulationOutput3)
deviance(mod.visited)


###* Try models to see the effect of elevation on the visited flowers per minute
mod.visited.twed <- glmmTMB(visited.flowers.per.minute ~ elevation  + (1|plant.species),
                            data = final.table,
                            family = tweedie(link="log")
)

mod.visited.twed.null <- glmmTMB(visited.flowers.per.minute ~ 1  + (1|plant.species),
                                 data = final.table,
                                 family = tweedie(link="log")
)

AIC(mod.visited.twed, mod.visited.twed.null)
emm <- emmeans(mod.visited.twed, ~elevation)
pairs(emm)
pairs(emm, type = "response")

summary(mod.visited.twed)
simulationOutput3 <- simulateResiduals(fittedModel = mod.visited.twed, plot = F)
testDispersion(simulationOutput3)
testZeroInflation(simulationOutput3)
testResiduals(simulationOutput3)
plot(simulationOutput3)
deviance(mod.visited)
