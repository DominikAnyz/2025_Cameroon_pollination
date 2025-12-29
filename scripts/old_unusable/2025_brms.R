###* In this script are all of the scripts for the glmmTMB models which I fit for
###* individual treatments/measurements and explaining how I got to them
###* 
###* Loading the necessary packages
pacman::p_load(tidyverse, glmmTMB, DHARMa, emmeans)

select <- dplyr::select

set.seed(123)

###* 
###* Load data----

seed.data<- read.delim("data/clean_seeds.txt", na = c("na"))

seed.data$elevation<-as.factor(seed.data$elevation)

seed.data <- seed.data %>%
  mutate(elevation = recode(elevation,
                            `2300` = 2300,
                            `2800` = 2800,
                            `3500` = 3400,
                            `4000` = 3800))

View(seed.data)

###* 
###* Create new datset "seed.indices" with some configuraions
###* The following code creates several new columns

seed.indices <- 
  seed.data %>% 
  group_by(elevation,plant_number, species) %>% 
  #* O.mean.plantnumber will return the mean seedset from the outcrossing per
  #* elevation, plant.number and species 
  mutate(O.mean.plantnumber = mean(seedset[treatment == "outcrossing"], na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(elevation, species) %>% 
  #* O.mean.elevation will return the mean seedset from the outcrossing treatment
  #* for a given elevation and species
  mutate(O.mean.elevation = mean(seedset[treatment == "outcrossing"], na.rm = TRUE)) %>%
  #* O.mean will combine O.mean.plantnumber and O.mean.elevation, specifically 
  #* when it is possible to use O.mean.plantnumber, it is used, however in cases
  #* when a given plant did not have any outcrossing results, the mean for the 
  #* elevation will be used
  mutate(O.mean = case_when(is.nan(O.mean.plantnumber) ~ O.mean.elevation,
                            .default = O.mean.plantnumber)) %>% 
  #* create a ne column named index, which calculats the index
  mutate(index = seedset/O.mean) %>%
  #* round the seedset, if it is not rounded for some reason
  mutate(seedset = if_else(seedset != round(seedset), round(seedset), seedset))

###* If the seedset is NA or Inf, then change to 0

seed.indices <- seed.indices %>% 
  filter(!(index %in% c("NA", "Inf")|is.na(seedset))) %>% 
  mutate(index = case_when(
    is.finite (index) ~index,
    .default = 0
  ))

###* Elevation and species as factor

seed.indices$elevation<-as.factor(seed.indices$elevation)
seed.indices$species<-as.factor(seed.indices$species)

###* Create a dataset specifically for the control treatment

c.index <- seed.indices %>%
  filter(treatment == "control") %>%
  select(seedset, index, elevation, species, plant_number) %>%
  mutate(seedset = round (seedset)) %>%
  #* calculate the pollen limitation index
  mutate(PL.index = replace(1 - index, 1 - index < 0, 0)) %>%
  mutate(elevation = factor(elevation),
         species = factor(species)) %>%
  tibble::rowid_to_column("ID") %>%
  #* create column plant.id which will be plant specific, since it will have all
  #* elevation, species and plant.number in the name
  mutate(plant.id = as.factor(paste0(elevation,species,plant_number))) %>%
  #* filter out Hypericum in elevation 4000, since the plant didn't produce any
  #* seeds in the highest elevation and this could scew the results
  #filter(species != "Hypericum r" | elevation != 3800) %>%
  mutate(flower.id = as.factor(paste0(elevation, species, plant_number,"_", ID)))










###* CONTROL SEEDSET
###* 
###* I have created four models based on the distiburtion of my response, seedset
###* 
###* ZINB
c_brm_model_nb_new <- brm(
  formula = bf(seedset ~ elevation + (1|species) + (1|plant.id)#,
               #zi ~ elevation
               ),
  data = c.index,
  family = negbinomial(),
  prior = c(
    set_prior("normal(0, 1.5)", class = "b"),
    #set_prior("normal(0, 2)", class = "b", dpar = "zi"),
    set_prior("student_t(3, 0, 5)", class = "sd")
  ),
  chains = 5, cores = 5,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  iter = 2000, warmup = 1000,
  control = list(max_treedepth = 15, adapt_delta = 0.95),
  seed = 1234
)

###* Save models with different distributions
# saveRDS(c_brm_model_zinb_new, "../brms_models/c_brm_model_zinb_new.rds")
# saveRDS(c_brm_model_po_new, "../brms_models/c_brm_model_po_new.rds")
# saveRDS(c_brm_model_zipo_new, "../brms_models/c_brm_model_zipo_new.rds")
# saveRDS(c_brm_model_nb_new, "../brms_models/c_brm_model_nb_new.rds")

###* Load models with different distributions
c_brm_model_po_new <- readRDS("../brms_models/c_brm_model_po_new.rds")
c_brm_model_zipo_new <- readRDS("../brms_models/c_brm_model_zipo_new.rds")
c_brm_model_nb_new <- readRDS("../brms_models/c_brm_model_nb_new.rds")
c_brm_model_zinb_new <- readRDS("../brms_models/c_brm_model_zinb_new.rds")

###* Check loo's of models
loo(c_brm_model_po_new, c_brm_model_zipo_new, c_brm_model_nb_new, c_brm_model_zinb_new)
loo(c_brm_model_zinb_new)
###* zinb seems to be best, but there are pareto warings

###* Getting rid of pareto warnings by using kfold
###* 
# c_brm_model_zinb_new_k <- kfold(c_brm_model_zinb_new, K = 5, cores = 5)
# saveRDS(c_brm_model_zinb_new_k, "../brms_models/c_brm_model_zinb_new_k.rds")
# c_brm_model_po_new_k <- kfold(c_brm_model_po_new, K = 5, cores = 5)
# saveRDS(c_brm_model_po_new_k, "../brms_models/c_brm_model_po_new_k.rds")
# c_brm_model_zipo_new_k <- kfold(c_brm_model_zipo_new, K = 5, cores = 5)
# saveRDS(c_brm_model_zipo_new_k, "../brms_models/c_brm_model_zipo_new_k.rds")
# c_brm_model_nb_new_k <- kfold(c_brm_model_nb_new, K = 5, cores = 5)
# saveRDS(c_brm_model_nb_new_k, "../brms_models/c_brm_model_nb_new_k.rds")

###* Prepared loo images for models
c_brm_model_zinb_new_k <- readRDS("../brms_models/c_brm_model_zinb_new_k.rds")
c_brm_model_po_new_k <- readRDS("../brms_models/c_brm_model_po_new_k.rds")
c_brm_model_zipo_new_k <- readRDS("../brms_models/c_brm_model_zipo_new_k.rds")
c_brm_model_nb_new_k <- readRDS("../brms_models/c_brm_model_nb_new_k.rds")

###* Comparing models after kfolding
loo_compare(c_brm_model_po_new_k, c_brm_model_zipo_new_k, c_brm_model_nb_new_k, c_brm_model_zinb_new_k)
###* Zinb model comes out on top

###* Creating null model to my zinb model to see, whether full model shows
###* better support
###* NULL
c_brm_null_zinb_new <- brm(
  formula = bf(seedset ~ 1 + (1|species) + (1|plant.id),
               zi ~ elevation),
  data = c.index,
  family = zero_inflated_negbinomial(),
  prior = c(
    set_prior("student_t(3, 0, 5)", class = "sd"),
    set_prior("normal(0, 2)", class = "b", dpar = "zi")
  ),
  chains = 5, cores = 5,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  iter = 2000, warmup = 1000,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)

###* Save model
saveRDS(c_brm_null_zinb_new, file = "../brms_models/c_brm_null_zinb_new.rds")
c_brm_null_zinb_new <- readRDS("../brms_models/c_brm_null_zinb_new.rds")
loo(c_brm_null_zinb_new)
###* Also problem with pareto

###* Fix pareto problem in null
c_brm_null_zinb_new_k <- kfold(c_brm_null_zinb_new, K = 5, cores = 5)
saveRDS(c_brm_null_zinb_new_k, "../brms_models/c_brm_null_zinb_new_k.rds")

###* Fit kfold fixed model
c_brm_null_zinb_new_k <- readRDS("../brms_models/c_brm_null_zinb_new_k.rds")

c_brm_model_zinb_new_k
c_brm_null_zinb_new_k

loo_compare(c_brm_null_zinb_new_k, c_brm_model_zinb_new_k)

###* Clearly, the model explains more than the null model
###* For clarification, how does loo work?
###* 
###* For every observation, loo asks "how well can the model predict this point,
###* if it had been left out when fitting"?
###* 
###* The results is expected log predictive density (elpd_loo)
###* this is the sum of log predictive densities across all points under loo CV
###* SE is the uncertainty for the estinate
###* -> the higher the eldp_loo, the better predictive accuracy
###* 
###* p_loo / effective number of parameters (model complexity penalty)
###* 
###* elpd_diff, se_diff
###* if Δelpd > ~2 × se_diff, then evidence for the better model is strong
###* 
###* pareto k diagnostics
###* Each observation has a “k” diagnostic telling whether the LOO estimate is reliable
###* k < 0.7 = good; 0.7–1 = problematic; >1 = unreliable.
###* 
###* So now that we know which model is the best and that it is more reliable
###* than the null, can it be trusted?
###* 
c_brm_model_zinb_new <- readRDS("../brms_models/c_brm_model_zinb_new.rds")
summary(c_brm_model_zinb_new)

pp_check(c_brm_model_zinb_new, type = "dens_overlay")
pp_check(c_brm_model_zinb_new, type = "hist")
pp_check(c_brm_model_zinb_new, type = "stat", stat = "mean")
pp_check(c_brm_model_zinb_new, type = "stat", stat = "sd")

###* compute residuals and check patterns unexplained by model
resids <- residuals(c_brm_model_zinb_new, type = "pearson")
plot(resids)      # residual distribution
hist(resids)

bayes_R2(c_brm_model_zinb_new)

























###POLLEN LIMITATION
c.index <- c.index %>%
  filter(species != "Hypericum r" | elevation != 3800)

library(brms)
library(ordbetareg)

# pl_ord_model <- ordbetareg(
#   formula = PL.index ~ elevation + (1|species) + (1|plant.id),
#   data = c.index,
#   chains = 10, iter = 12000, warmup = 5000, cores = 10,
#   control = list(adapt_delta = 0.99),
#   seed = 1234
# )
# 
# saveRDS(pl_ord_model, file = "../brms_models/pl_ord_model.rds")
# pl_ord_model <- readRDS("../brms_models/pl_ord_model.rds")
# 
# c.index$intercept_only <- 1
# 
# pl_ord_null <- ordbetareg(
#   formula = PL.index ~ intercept_only + (1|species) + (1|plant.id),
#   data = c.index,
#   chains = 10, iter = 12000, warmup = 5000, cores = 10,
#   control = list(adapt_delta = 0.99),
#   seed = 1234
# )
# 
# saveRDS(pl_ord_null, file = "../brms_models/pl_ord_null.rds")
# pl_ord_null <- readRDS("../brms_models/pl_ord_null.rds")


pl_ord_model_new <- ordbetareg(
  formula = PL.index ~ elevation + (1|species) + (1|plant.id),
  data = c.index,
  chains = 5, cores = 5,
  iter = 2000, warmup = 1000,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  control = list(adapt_delta = 0.99),
  seed = 1234
)

saveRDS(pl_ord_model_new, file = "../brms_models/pl_ord_model_new.rds")
pl_ord_model_new <- readRDS("../brms_models/pl_ord_model_new.rds")
loo(pl_ord_model_new)

pl_ord_model_new_k <- kfold(pl_ord_model_new, K = 5, cores = 5)
saveRDS(pl_ord_model_new_k, "../brms_models/pl_ord_model_new_k.rds")
pl_ord_model_new_k <- readRDS("../brms_models/pl_ord_model_new_k.rds")

c.index$intercept_only <- 1

pl_ord_null_new <- ordbetareg(
  formula = PL.index ~ intercept_only + (1|species) + (1|plant.id),
  data = c.index,
  chains = 5, cores = 5,
  iter = 2000, warmup = 1000, 
  save_pars = save_pars(latent = TRUE, all = TRUE), 
  control = list(adapt_delta = 0.99),
  seed = 1234
)

saveRDS(pl_ord_null_new, file = "../brms_models/pl_ord_null_new.rds")
pl_ord_null_new <- readRDS("../brms_models/pl_ord_null_new.rds")

pl_ord_null_new_k <- kfold(pl_ord_null_new, K = 5, cores = 5)
saveRDS(pl_ord_null_new_k, "../brms_models/pl_ord_null_new_k.rds")
pl_ord_null_new_k <- readRDS("../brms_models/pl_ord_null_new_k.rds")

loo_compare(pl_ord_null_new_k, pl_ord_model_new_k)

###* yes, model is better than null
###* MODEL OUTPUT
pl_ord_model_new <- readRDS("../brms_models/pl_ord_model_new.rds")
summary(pl_ord_model_new)

pp_check(pl_ord_model_new, type = "dens_overlay")
pp_check(pl_ord_model_new, type = "hist")
pp_check(pl_ord_model_new, type = "stat", stat = "mean")
pp_check(pl_ord_model_new, type = "stat", stat = "sd")

###* compute residuals and check patterns unexplained by model
resids <- residuals(pl_ord_model_new, type = "pearson")
plot(resids)      # residual distribution
hist(resids)

bayes_R2(pl_ord_model_new)


















###* AO index
ao.index <- seed.indices %>%
  filter(treatment == "autogamy") %>%
  select(index, elevation, species, plant_number, plant_number) %>%
  mutate(elevation = factor(elevation),
         species = factor(species)) %>%
  mutate(index = round(index * 100)) %>%
  mutate(plant.number.el = as.factor(paste0(elevation,plant_number))) %>%
  mutate(plant.id = as.factor(paste0(elevation,species,plant_number))) %>%
  filter(species != "Hypericum r" | elevation != 3800) %>%
  mutate(species.sp = as.factor(paste0(elevation,species)))

ao_brm_model_zipo_new <- brm(
  formula = bf(index ~ elevation + (1|species) + (1|plant.id),
               zi ~ elevation
               ),
  data = ao.index,
  family = poisson(),
  prior = c(
    set_prior("normal(0, 1.5)", class = "b"),
    set_prior("normal(0, 2)", class = "b", dpar = "zi"),
    set_prior("student_t(3, 0, 5)", class = "sd")
  ),
  chains = 5, cores = 5,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  iter = 2000, warmup = 1000,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)



###* Save models with different distributions
saveRDS(ao_brm_model_po_new, "../brms_models/ao_brm_model_po_new.rds")
saveRDS(ao_brm_model_zipo_new, "../brms_models/ao_brm_model_zipo_new.rds")
saveRDS(ao_brm_model_nb_new, "../brms_models/ao_brm_model_nb_new.rds")
saveRDS(ao_brm_model_zinb_new, "../brms_models/ao_brm_model_zinb_new.rds")

###* Load models with different distributions
ao_brm_model_po_new <- readRDS("../brms_models/ao_brm_model_po_new.rds")
ao_brm_model_zipo_new <- readRDS("../brms_models/ao_brm_model_zipo_new.rds")
ao_brm_model_nb_new <- readRDS("../brms_models/ao_brm_model_nb_new.rds")
ao_brm_model_zinb_new <- readRDS("../brms_models/ao_brm_model_zinb_new.rds")

###* Check loo's of models
loo(ao_brm_model_po_new, ao_brm_model_zipo_new, ao_brm_model_nb_new, ao_brm_model_zinb_new)
loo(ao_brm_model_zinb_new)
###* zinb seems to be best, but there are pareto warings

###* Getting rid of pareto warnings by using kfold
###* 
# ao_brm_model_po_new_k <- kfold(ao_brm_model_po_new, K = 2, cores = 5)
# ao_brm_model_po_new_k
# saveRDS(ao_brm_model_po_new_k, "../brms_models/ao_brm_model_po_new_k.rds")
# ao_brm_model_zipo_new_k <- kfold(ao_brm_model_zipo_new, K = 3, cores = 5)
# ao_brm_model_zipo_new_k
# saveRDS(ao_brm_model_zipo_new_k, "../brms_models/ao_brm_model_zipo_new_k.rds")
# ao_brm_model_nb_new_k <- kfold(ao_brm_model_nb_new, K = 3, cores = 5)
# ao_brm_model_nb_new_k
# saveRDS(ao_brm_model_nb_new_k, "../brms_models/ao_brm_model_nb_new_k.rds")
# ao_brm_model_zinb_new_k <- kfold(ao_brm_model_zinb_new, K = 3, cores = 5)
# ao_brm_model_zinb_new_k
# saveRDS(ao_brm_model_zinb_new_k, "../brms_models/ao_brm_model_zinb_new_k.rds")


###* Prepared loo images for models
ao_brm_model_zinb_new_k <- readRDS("../brms_models/ao_brm_model_zinb_new_k.rds")
ao_brm_model_po_new_k <- readRDS("../brms_models/ao_brm_model_po_new_k.rds")
ao_brm_model_zipo_new_k <- readRDS("../brms_models/ao_brm_model_zipo_new_k.rds")
ao_brm_model_nb_new_k <- readRDS("../brms_models/ao_brm_model_nb_new_k.rds")


###* null
ao_brm_null_zinb <- brm(
  formula = bf(index ~ 1 + (1|species) + (1|plant.id),
               zi ~ elevation),
  data = ao.index,
  family = zero_inflated_negbinomial(),
  prior = c(
    set_prior("student_t(3, 0, 5)", class = "sd"),
    set_prior("normal(0, 2)", class = "b", dpar = "zi")
  ),
  chains = 5, cores = 5,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  iter = 5000, warmup = 2000,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)

saveRDS(ao_brm_null_zinb, file = "../brms_models/ao_brm_null_zinb.rds")
ao_brm_null_zinb <- readRDS("../brms_models/ao_brm_null_zinb.rds")

loo(ao_brm_null_zinb, ao_brm_model_zinb)
###* model betterthan null
###* 
###* 
###* MODEL OUTPUT
ao_brm_model_zinb <- readRDS("../brms_models/ao_brm_model_zinb.rds")
summary(ao_brm_model_zinb)

pp_check(ao_brm_model_zinb, type = "dens_overlay")
pp_check(ao_brm_model_zinb, type = "hist")
pp_check(ao_brm_model_zinb, type = "stat", stat = "mean")
pp_check(ao_brm_model_zinb, type = "stat", stat = "sd")

###* compute residuals and check patterns unexplained by model
resids <- residuals(ao_brm_model_zinb, type = "pearson")
plot(resids)      # residual distribution
hist(resids)

bayes_R2(ao_brm_model_zinb)








###GEITONOGAMY
go.index <- seed.indices %>%
  filter(treatment == "geitonogamy") %>%
  select(index, elevation, species, plant_number) %>%
  mutate(elevation = factor(elevation),
         species = factor(species)) %>%
  mutate(index100 = round(index * 100)) %>%
  mutate(plant.id = as.factor(paste0(elevation,species,plant_number)))%>%
  filter(species != "Lactuca i") %>%
  filter(species != "Hypericum r" | elevation != 3800)

go_brm_model_zinb_new <- brm(
  formula = bf(index100 ~ elevation + (1|species) + (1|plant.id),
               zi ~ elevation
  ),
  data = go.index,
  family = negbinomial(),
  prior = c(
    set_prior("normal(0, 1.5)", class = "b"),
    set_prior("normal(0, 2)", class = "b", dpar = "zi"),
    set_prior("student_t(3, 0, 5)", class = "sd")
  ),
  chains = 5, cores = 5,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  iter = 2000, warmup = 1000,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)



###* Save models with different distributions
saveRDS(go_brm_model_po_new, "../brms_models/go_brm_model_po_new.rds")
saveRDS(go_brm_model_zipo_new, "../brms_models/go_brm_model_zipo_new.rds")
saveRDS(go_brm_model_nb_new, "../brms_models/go_brm_model_nb_new.rds")
saveRDS(go_brm_model_zinb_new, "../brms_models/go_brm_model_zinb_new.rds")

###* Load models with different distributions
go_brm_model_po_new <- readRDS("../brms_models/go_brm_model_po_new.rds")
go_brm_model_zipo_new <- readRDS("../brms_models/go_brm_model_zipo_new.rds")
go_brm_model_nb_new <- readRDS("../brms_models/go_brm_model_nb_new.rds")
go_brm_model_zinb_new <- readRDS("../brms_models/go_brm_model_zinb_new.rds")

###* Check loo's of models
loo(go_brm_model_po_new, go_brm_model_zipo_new, go_brm_model_nb_new, go_brm_model_zinb_new)
loo(go_brm_model_zinb_new)
###* zinb seems to be best, but there are pareto warings

###* Getting rid of pareto warnings by using kfold
###* 
# go_brm_model_po_new_k <- kfold(go_brm_model_po_new, K = 2, cores = 5)
# go_brm_model_po_new_k
# saveRDS(go_brm_model_po_new_k, "../brms_models/go_brm_model_po_new_k.rds")
# go_brm_model_zipo_new_k <- kfold(go_brm_model_zipo_new, K = 3, cores = 5)
# go_brm_model_zipo_new_k
# saveRDS(go_brm_model_zipo_new_k, "../brms_models/go_brm_model_zipo_new_k.rds")
# go_brm_model_nb_new_k <- kfold(go_brm_model_nb_new, K = 3, cores = 5)
# go_brm_model_nb_new_k
# saveRDS(go_brm_model_nb_new_k, "../brms_models/go_brm_model_nb_new_k.rds")
# go_brm_model_zinb_new_k <- kfold(go_brm_model_zinb_new, K = 3, cores = 5)
# go_brm_model_zinb_new_k
# saveRDS(go_brm_model_zinb_new_k, "../brms_models/go_brm_model_zinb_new_k.rds")


###* Prepared loo images for models
go_brm_model_zinb_new_k <- readRDS("../brms_models/go_brm_model_zinb_new_k.rds")
go_brm_model_po_new_k <- readRDS("../brms_models/go_brm_model_po_new_k.rds")
go_brm_model_zipo_new_k <- readRDS("../brms_models/go_brm_model_zipo_new_k.rds")
go_brm_model_nb_new_k <- readRDS("../brms_models/go_brm_model_nb_new_k.rds")

loo_compare(go_brm_model_po_new_k, go_brm_model_zipo_new_k, go_brm_model_nb_new_k, go_brm_model_zinb_new_k)

###* null
go_brm_null_zinb <- brm(
  formula = bf(index ~ 1 + (1|species) + (1|plant.id),
               zi ~ elevation),
  data = go.index,
  family = zero_inflated_negbinomial(),
  prior = c(
    set_prior("student_t(3, 0, 5)", class = "sd"),
    set_prior("normal(0, 2)", class = "b", dpar = "zi")
  ),
  chains = 5, cores = 5,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  iter = 2000, warmup = 1000,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)

saveRDS(go_brm_null_zinb, file = "../brms_models/go_brm_null_zinb.rds")
go_brm_null_zinb <- readRDS("../brms_models/go_brm_null_zinb.rds")

loo(go_brm_null_zinb, go_brm_model_zinb)
###* model better than null
###* 
###* 



































library(brms)
library(dplyr)

c_brm_model_zinb <- readRDS("../brms_models/c_brm_model_zinb.rds")
c_brm_null_zinb <- readRDS("../brms_models/c_brm_null_zinb.rds")
summary(c_brm_model_zinb)
summary(c_brm_model_zinb_new)
loo(c_brm_null_zinb, c_brm_model_zinb)
pl_ord_model <- readRDS("../brms_models/pl_ord_model.rds")
pl_ord_null <- readRDS("../brms_models/pl_ord_null.rds")
loo(pl_ord_null, pl_ord_model)
ao_brm_model_zinb <- readRDS("../brms_models/ao_brm_model_zinb.rds")
ao_brm_null_zinb <- readRDS("../brms_models/ao_brm_null_zinb.rds")
loo(ao_brm_null_zinb, ao_brm_model_zinb)
go_brm_model_zinb <- readRDS("../brms_models/go_brm_model_zinb.rds")
go_brm_null_zinb <- readRDS("../brms_models/go_brm_null_zinb.rds")
loo(go_brm_null_zinb, go_brm_model_zinb)
m_brm_model_nb <- readRDS("../brms_models/m_brm_model_nb.rds")
m_brm_null_nb <- readRDS("../brms_models/m_brm_null_nb.rds")
loo(m_brm_null_nb, m_brm_model_nb)
f_brm_model_po <- readRDS("../brms_models/f_brm_model_po.rds")
f_brm_null_po <- readRDS("../brms_models/f_brm_null_po.rds")
loo(f_brm_null_po, f_brm_model_po)
vis_ord_model <- readRDS("../brms_models/vis_ord_model.rds")
vis_ord_null <- readRDS("../brms_models/vis_ord_null.rds")
loo(vis_ord_model, vis_ord_null)

# Helper: extract info from loo comparison
extract_loo <- function(null_model, full_model, dist, zi = FALSE, metric_name) {
  cmp <- loo(null_model, full_model)
  diff <- cmp$diffs[2, "elpd_diff"]
  se <- cmp$diffs[2, "se_diff"]
  z <- diff / se
  tibble(
    Metric = metric_name,
    Distribution = dist,
    Zero_inflation = ifelse(zi, "Yes", "No"),
    Delta_LOO = round(diff, 1),
    SE = round(se, 1),
    z_ELPD = round(z, 2),
    `Include elevation effect?` = ifelse(diff > 2*se, "Yes", "No")
  )
}

# Build table row by row
results <- bind_rows(
  extract_loo(c_brm_null_zinb, c_brm_model_zinb, "NegBinomial", TRUE,  "Natural seed-set"),
  extract_loo(pl_ord_null,     pl_ord_model,     "Ordered Beta", FALSE, "Pollen limitation"),
  extract_loo(ao_brm_null_zinb, ao_brm_model_zinb, "NegBinomial", TRUE, "Autogamy index"),
  extract_loo(go_brm_null_zinb, go_brm_model_zinb, "NegBinomial", TRUE, "Geitonogamy index"),
  extract_loo(m_brm_null_nb,   m_brm_model_nb,   "NegBinomial", FALSE, "Morphospecies richness"),
  extract_loo(f_brm_null_po,   f_brm_model_po,   "Poisson",     FALSE, "Functional group richness"),
  extract_loo(vis_ord_null,    vis_ord_model,    "Ordered Beta", FALSE, "Visitation frequency")
)

print(results)

readr::write_csv(results, "tables/new_bayes.csv")

summary(c_brm_model_zinb)






library(report)
r <- report(c_brm_model_zinb, verbose = TRUE)
r



library(emmeans)
library(dplyr)

# Example: contrasts for Natural seed-set model
c_emm <- emmeans(c_brm_model_zinb, ~ elevation)
c_contr <- contrast(c_emm, method = "pairwise") %>%
  as.data.frame()

# Format columns
c_contr <- c_contr %>%
  mutate(
    Metric = "Natural seed-set",
    Estimate = round(estimate, 2),
    l95 = round(lower.HPD, 2),
    u95 = round(upper.HPD, 2)
  ) %>%
  select(Metric, contrast, Estimate, l95, u95)

# Pollen limitation
pl_emm <- emmeans(pl_ord_model, ~ elevation)
pl_contr <- contrast(pl_emm, method = "pairwise") %>%
  as.data.frame() %>%
  mutate(Metric = "Pollen limitation",
         Estimate = round(estimate, 2),
         l95 = round(lower.HPD, 2),
         u95 = round(upper.HPD, 2)) %>%
  select(Metric, contrast, Estimate, l95, u95)

# ... repeat for AO, GO, M, F, Visitation

# Autogamy index
ao_emm <- emmeans(ao_brm_model_zinb, ~ elevation)
ao_contr <- contrast(ao_emm, method = "pairwise") %>%
  as.data.frame() %>%
  mutate(Metric = "Autogamy",
         Estimate = round(estimate, 2),
         l95 = round(lower.HPD, 2),
         u95 = round(upper.HPD, 2)) %>%
  select(Metric, contrast, Estimate, l95, u95)

# Geitonogamy index
summary(go_brm_model_zinb_new)

go_emm <- emmeans(go_brm_model_zinb_new, ~ elevation)
go_contr <- contrast(go_emm, method = "pairwise") %>%
  as.data.frame() %>%
  mutate(Metric = "Geitonogamy",
         Estimate = round(estimate, 2),
         l95 = round(lower.HPD, 2),
         u95 = round(upper.HPD, 2)) %>%
  select(Metric, contrast, Estimate, l95, u95)

go_contr
# Morpho richness
m_emm <- emmeans(m_brm_model_nb, ~ elevation)
m_contr <- contrast(m_emm, method = "pairwise") %>%
  as.data.frame() %>%
  mutate(Metric = "Morphospecies richness",
         Estimate = round(estimate, 2),
         l95 = round(lower.HPD, 2),
         u95 = round(upper.HPD, 2)) %>%
  select(Metric, contrast, Estimate, l95, u95)

# Functional richness
f_emm <- emmeans(f_brm_model_po, ~ elevation)
f_contr <- contrast(f_emm, method = "pairwise") %>%
  as.data.frame() %>%
  mutate(Metric = "Functional group richness",
         Estimate = round(estimate, 2),
         l95 = round(lower.HPD, 2),
         u95 = round(upper.HPD, 2)) %>%
  select(Metric, contrast, Estimate, l95, u95)

# Visitation frequency
v_emm <- emmeans(vis_ord_model, ~ elevation)
v_contr <- contrast(v_emm, method = "pairwise") %>%
  as.data.frame() %>%
  mutate(Metric = "Visitation frequency",
         Estimate = round(estimate, 2),
         l95 = round(lower.HPD, 2),
         u95 = round(upper.HPD, 2)) %>%
  select(Metric, contrast, Estimate, l95, u95)

all_contrasts <- bind_rows(c_contr, pl_contr, ao_contr, m_contr, f_contr, v_contr)

all_contrasts <- all_contrasts %>%
  mutate(
    CI_95 = paste0("[", l95, ", ", u95, "]")
  ) %>%
  select(Metric, contrast, Estimate, CI_95)

print(all_contrasts)

readr::write_csv(all_contrasts, "tables/new_bayes_post.csv")

###* Natural seedset 

























