###* I would like to try and carry out the same models, which I carried out in 
###* "2025_glmmTMB.R", except with Bayesian statistics
source("scripts/setup.R")
source("scripts/Bayesian setup.R")

###* Seedset looked like this
C.glmer3 <- glmmTMB(seedset ~ elevation  + (1|species) + (1|plant.id),
                    ziformula=~elevation,
                    data = c.index,
                    #family = poisson)
                    family = nbinom2)
# 

formula_seedset <- bf(seedset ~ elevation + (1|species) + (1|plant.id))

glmm_seedset <- brm(
  formula = formula_seedset,
  data = c.index,
  family = gaussian(),
  chains = 10,
  cores = 10,
  iter = 4000,
  warmup = 1000,
  control = list(adapt_delta = 0.999, max_treedepth = 20),
  seed = 1234
)

summary(glmm_seedset)
saveRDS(glmm_seedset, file = "brms_glmm_seedset.rds")
glmm_seedset <- readRDS("brms_glmm_seedset.rds")
###* Very high uncertainty

###* poisson
formula_seedset_p <- bf(seedset ~ elevation + (1|species) + (1|plant.id))

glmm_seedset_p <- brm(
  formula = formula_seedset,
  data = c.index,
  family = poisson(),
  chains = 10,
  cores = 10,
  iter = 4000,
  warmup = 1000,
  control = list(adapt_delta = 0.999, max_treedepth = 20),
  seed = 1234
)

summary(glmm_seedset_p)
saveRDS(glmm_seedset_p, file = "brms_glmm_seedset_p.rds")
glmm_seedset_p <- readRDS("brms_glmm_seedset_p.rds")

###* zero inflated poisson
formula_seedset_zip <- bf(seedset ~ elevation + (1|species) + (1|plant.id),
                          zi ~ elevation)

glmm_seedset_zip <- brm(
  formula = formula_seedset_zip,
  data = c.index,
  family = zero_inflated_poisson(),
  chains = 10,
  cores = 10,
  iter = 4000,
  warmup = 1000,
  control = list(adapt_delta = 0.999, max_treedepth = 20),
  seed = 1234
)

summary(glmm_seedset_zip)
saveRDS(glmm_seedset_zip, file = "brms_glmm_seedset_zip.rds")
glmm_seedset_zip <- readRDS("brms_glmm_seedset_zip.rds")

###* negative binomial
formula_seedset_nb <- bf(seedset ~ elevation + (1|species) + (1|plant.id))

glmm_seedset_nb <- brm(
  formula = formula_seedset_nb,
  data = c.index,
  family = negbinomial(),
  chains = 10,
  cores = 10,
  iter = 4000,
  warmup = 1000,
  control = list(adapt_delta = 0.999, max_treedepth = 20),
  seed = 1234
)

summary(glmm_seedset_nb)
saveRDS(glmm_seedset_nb, file = "brms_glmm_seedset_nb.rds")
glmm_seedset_nb <- readRDS("brms_glmm_seedset_nb.rds")


###* zero inflated negative binomial
formula_seedset_zinb <- bf(seedset ~ elevation + (1|species) + (1|plant.id),
                           zi ~ elevation)

glmm_seedset_zinb <- brm(
  formula = formula_seedset_zinb,
  data = c.index,
  family = zero_inflated_negbinomial(),
  chains = 10,
  cores = 10,
  iter = 4000,
  warmup = 1000,
  control = list(adapt_delta = 0.999, max_treedepth = 20),
  seed = 1234
)

summary(glmm_seedset_zinb)
saveRDS(glmm_seedset_zinb, file = "brms_glmm_seedset_zinb.rds")
glmm_seedset_zinb <- readRDS("brms_glmm_seedset_zinb.rds")


loo_gaussian <- loo(glmm_seedset)
loo_poisson <- loo(glmm_seedset_p)
loo_zip <- loo(glmm_seedset_zip)
loo_nb <- loo(glmm_seedset_nb)
loo_zinb <- loo(glmm_seedset_zinb)

loo_compare(loo_gaussian, loo_poisson, loo_zip, loo_nb, loo_zinb)

View(ao.final.table)




###* PL index looked liek this
# glm_model3 <- glmmTMB(PL.index ~ elevation + (1|species) +(1|plant.id), 
#                       data = c.index,
#                       family = ordbeta())
library(ordbetareg)
formula_PL_ordbeta <- bf(PL.index ~ elevation + (1|species) + (1|plant.id))

glmm_PL_ordbeta <- ordbetareg(
  formula = formula_PL_ordbeta,
  data = c.index,
  chains = 10,
  cores = 10,
  iter = 4000,
  warmup = 1000,
  control = list(adapt_delta = 0.999, max_treedepth = 20),
  seed = 1234
)

summary(glmm_PL_ordbeta)
saveRDS(glmm_PL_ordbeta, file = "brms_glmm_PL_ordbeta.rds")
glmm_PL_ordbeta <- readRDS("brms_glmm_PL_ordbeta.rds")



###* Autogamy index













