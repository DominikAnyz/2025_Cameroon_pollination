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

#View(c.index)

###* Fit a ZINB model to the data for the control treatment. This was after
###* trying various variations and finding out that the data did not fit well
###* with any other distributions. The other tested distributions were poisson,
###* zero inflated poisson and  negative binomial distribution.



###Control

library(brms)


c_brm_model_po <- brm(
  formula = bf(seedset ~ elevation + (1|species) + (1|plant.id)),
  data = c.index,
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

saveRDS(c_brm_model_po, file = "../brms_models/c_brm_model_po.rds")
summary(c_brm_model_po)
loo(c_brm_model, c_brm_model_po, moment_match = TRUE)
bayes_R2(c_brm_model)
bayes_R2(c_brm_model_po)


c_brm_model_zipo <- brm(
  formula = bf(seedset ~ elevation + (1|species) + (1|plant.id),
               zi ~ elevation),
  data = c.index,
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

saveRDS(c_brm_model_zipo, file = "../brms_models/c_brm_model_zipo.rds")
summary(c_brm_model_zipo)
loo(c_brm_model_po, c_brm_model_zipo, moment_match = TRUE)
bayes_R2(c_brm_model_po)
bayes_R2(c_brm_model_zipo)


c_brm_model_nb <- brm(
  formula = bf(seedset ~ elevation + (1|species) + (1|plant.id)),
  data = c.index,
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

saveRDS(c_brm_model_nb, file = "../brms_models/c_brm_model_nb.rds")
summary(c_brm_model_nb)
loo(c_brm_model_pozi, c_brm_model_nb, moment_match = TRUE)
bayes_R2(c_brm_model_pozi)
bayes_R2(c_brm_model_nb)



c_brm_model_zinb <- brm(
  formula = bf(seedset ~ elevation + (1|species) + (1|plant.id),
               zi ~ elevation),
  data = c.index,
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

saveRDS(c_brm_model_zinb, file = "../brms_models/c_brm_model_zinb.rds")
summary(c_brm_model_zinb)
loo(c_brm_model_po, c_brm_model_zipo, c_brm_model_nb, c_brm_model_zinb)
bayes_R2(c_brm_model_po)
bayes_R2(c_brm_model_zipo)
bayes_R2(c_brm_model_nb)
bayes_R2(c_brm_model_zinb)


summary(c_brm_model_nb)
summary(c_brm_model_zinb)

































library(brms)

# Fit model
c_brm_model <- brm(
  formula = bf(seedset ~ elevation + (1|species) + (1|plant.id)),
  data = c.index,
  family = zero_inflated_negbinomial(),
  prior = c(
    set_prior("normal(0, 2)", class = "b"),  # elevation effect (main part)
    set_prior("student_t(3, 0, 5)", class = "sd")  # random effects
  ),
  chains = 10, cores = 10,
  iter = 12000, warmup = 5000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)

saveRDS(c_brm_model, file = "../brms_models/c_brm_model.rds")
c_brm_model <- readRDS("../brms_models/c_brm_model.rds")

# Null model (no elevation effect)
c_brm_null <- brm(
  formula = bf(seedset ~ 1 + (1|species) + (1|plant.id)),
  data = c.index,
  family = zero_inflated_negbinomial(),
  prior = c(
    set_prior("student_t(3, 0, 5)", class = "sd")  # random effects
  ),
  chains = 10, cores = 10,
  iter = 8000, warmup = 3000,
  control = list(adapt_delta = 0.99),
  seed = 1234
)

saveRDS(c_brm_null, file = "../brms_models/c_brm_null.rds")
c_brm_null <- readRDS("../brms_models/c_brm_null.rds")

# Pairwise comparisons
pairs(emm)

emm <- emmeans(C.glmer3, ~elevation)
pairs(emm)

# Model comparison
bayes_R2(c_brm_model, c_brm_null)
loo(c_brm_model, c_brm_null)









fit <- c_brm_model

summary(fit)
pp_check(fit, type = "dens_overlay")
pp_check(fit, type = "stat")
bayes_R2(fit)
plot(fit)
fitted_values <- fitted(fit)
residuals <- residuals(fit)
plot(fitted_values, residuals)
abline(h = 0, col = "red")

library(emmeans)

# Estimated marginal means (posterior marginal means) by elevation
c_emm_brm <- emmeans(c_brm_model, ~ elevation)

pairs(c_emm_brm)


# Compact letter display
cld_c_brm <- multcomp::cld(
  c_emm_brm,
  adjust = "none",     # Tukey-style adjustment is not meaningful in Bayesian
  Letters = letters,
  alpha = 0.05         # 95% HPDIs used to decide overlap
) %>%
  as.data.frame() %>%
  mutate(
    elevation = factor(elevation, levels = c("2300", "2800", "3400", "3800")),
    .group = str_trim(.group)
  )

# Max y-value for positioning letters
y_max_c <- max(c.index$seedset, na.rm = TRUE)

# Add letters to your violin plot
control_combined_plot_with_letters_brm <- control_combined_plot +
  geom_text(
    data = cld_c_brm,
    aes(x = elevation, y = y_max_c + 40, label = .group),
    size = 7,
    fontface = "bold"
  )

control_combined_plot_with_letters_brm









library(emmeans)
library(coda)

# Estimated marginal means for the brms model
c_emm_brm <- emmeans(c_brm_model, ~ elevation)

# Pairwise contrasts
contr_brm <- contrast(c_emm_brm, method = "pairwise")

# Extract posterior draws (mcmc.list)
contr_mcmc <- as.mcmc(contr_brm)

# Combine chains into one matrix
contr_draws <- do.call(rbind, lapply(contr_mcmc, as.matrix))

# Check structure
dim(contr_draws)
head(colnames(contr_draws))





# Posterior probability that difference > 0
post_probs <- apply(contr_draws, 2, function(x) mean(x > 0))

# Format results
post_probs_df <- data.frame(
  contrast = colnames(contr_draws),
  P_greater0 = round(post_probs, 3),
  P_less0 = round(1 - post_probs, 3)
)

post_probs_df

# Posterior probability that difference > 0
post_probs <- apply(contr_draws, 2, function(x) mean(x > 0))

# Format results
post_probs_df <- data.frame(
  contrast = colnames(contr_draws),
  P_greater0 = round(post_probs, 3),
  P_less0 = round(1 - post_probs, 3)
)

post_probs_df




library(multcompView)

library(dplyr)

library(multcompView)

make_bayes_cld <- function(contr_draws, alpha = 0.95) {
  # Posterior probabilities
  probs <- apply(contr_draws, 2, function(x) mean(x > 0))
  
  # Extract elevation names
  contrasts <- colnames(contr_draws)
  lvls <- unique(unlist(strsplit(contrasts, " - ")))
  
  # Build empty matrix
  comp_mat <- matrix(FALSE, nrow = length(lvls), ncol = length(lvls),
                     dimnames = list(lvls, lvls))
  
  # Fill matrix: TRUE if credibly different
  for (c in contrasts) {
    sp <- strsplit(c, " - ")[[1]]
    p <- probs[c]
    
    if (p > alpha | p < (1 - alpha)) {
      comp_mat[sp[1], sp[2]] <- TRUE
      comp_mat[sp[2], sp[1]] <- TRUE
    }
  }
  
  # Generate letters
  letters <- multcompView::multcompLetters(!comp_mat)$Letters
  return(letters)
}

# Example usage
cld_letters <- make_bayes_cld(contr_draws, alpha = 0.95)

cld_df <- data.frame(
  elevation = names(cld_letters),
  .group = cld_letters,
  stringsAsFactors = FALSE
)

cld_df


cld_letters


# Run it
cld_letters <- make_bayes_cld(contr_draws, alpha = 0.95)
cld_letters



# R

# Run on your contrast draws
cld_letters <- make_bayes_cld(contr_draws, alpha = 0.95)

# Make a tidy df for ggplot
cld_df <- data.frame(
  elevation = names(cld_letters),
  .group = cld_letters,
  stringsAsFactors = FALSE
)


cld_letters























###POLLEN LIMITATION
c.index <- c.index %>%
  filter(species != "Hypericum r" | elevation != 3800)

# Full model
pl_brm_model <- brm(
  formula = bf(PL.index ~ elevation + (1|species) + (1|plant.id)),
  data = c.index,
  family = zero_one_inflated_beta(),
  prior = c(
    prior(normal(0, 2), class = "b"),
    prior(student_t(3, 0, 5), class = "sd")
  ),
  chains = 10, cores = 10,
  iter = 12000, warmup = 5000,
  control = list(adapt_delta = 0.99),
  seed = 1234
)

saveRDS(pl_brm_model, file = "../brms_models/pl_brm_model.rds")
pl_brm_model <- readRDS("../brms_models/pl_brm_model.rds")

summary(pl_brm_model)








library(ordbetareg)

pl_ord_model <- ordbetareg(
  formula = PL.index ~ elevation + (1|species) + (1|plant.id),
  data = c.index,
  chains = 10, iter = 12000, warmup = 5000, cores = 10,
  control = list(adapt_delta = 0.99),
  seed = 1234
)

saveRDS(pl_ord_model, file = "../brms_models/pl_ord_model.rds")
pl_ord_model <- readRDS("../brms_models/pl_ord_model.rds")

summary(pl_ord_model)
pp_check(pl_ord_model, type = "dens_overlay")
pp_check(pl_ord_model, type = "stat")







pr <- default_prior(
  PL.index ~ 1 + (1|species) + (1|plant.id),
  data = c.index
)

pr <- subset(pr, class != "b") 

# Now run the null model without the "b" prior
pl_ord_null <- ordbetareg(
  formula = PL.index ~ 1 + (1|species) + (1|plant.id),
  data = c.index,
  prior = pr,
  coef_prior_mean = NULL,
  coef_prior_SD = NULL,
  chains = 10, iter = 12000, warmup = 5000, cores = 10,
  control = list(adapt_delta = 0.99),
  seed = 1234
)




pl_ord_null <- ordbetareg(
  formula = PL.index ~ 1 +(1|species) + (1|plant.id),
  data = c.index,
  intercept_prior_mean = NULL,
  intercept_prior_SD = NULL,
  chains = 10, iter = 12000, warmup = 5000, cores = 10,
  control = list(adapt_delta = 0.99),
  seed = 1234
)

samanual_prior = saveRDS(pl_ord_null, file = "../brms_models/pl_ord_null.rds")
pl_ord_null <- readRDS("../brms_models/pl_ord_null.rds")






# Null model (no elevation effect)
pl_brm_null <- brm(
  formula = bf(PL.index ~ 1 + (1|species) + (1|plant.id)),
  data = c.index,
  family = beta_or(),
  prior = c(
    set_prior("student_t(3, 0, 5)", class = "sd") # random effects
  ),
  chains = 10, cores = 10,
  iter = 12000, warmup = 5000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)

saveRDS(pl_brm_null, file = "../brms_models/pl_brm_null.rds")
pl_brm_null <- readRDS("../brms_models/pl_brm_null.rds")































library(brms)
###AUTOGAMY
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


ao_brm_model <- brm(
  formula = bf(index ~ elevation + (1|species) + (1|plant.id),
               zi ~ elevation),
  data = ao.index,
  family = zero_inflated_poisson(),
  prior = c(
    set_prior("normal(0, 2)", class = "b"),          # elevation effect
    set_prior("student_t(3, 0, 5)", class = "sd"),   # random effects
    set_prior("normal(0, 2)", class = "b", dpar = "zi")  # zero-inflation
  ),
  chains = 10, cores = 10,
  iter = 7000, warmup = 2000,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)

saveRDS(ao_brm_model, file = "../brms_models/ao_brm_model.rds")
saveRDS(ao_brm_model, file = "../brms_models/ao_brm_pozi_model.rds")
ao_brm_model <- readRDS("../brms_models/ao_brm_model.rds")


summary(ao_brm_model)
pp_check(ao_brm_model, type = "dens_overlay")   # overlay histograms/densities
pp_check(ao_brm_model, type = "hist")           # distribution of simulated data
pp_check(ao_brm_model, type = "stat", stat = "mean") # means
pp_check(ao_brm_model, type = "stat", stat = "sd")   # variances
ppc_residuals(ao_brm_model)   # if you have bayesplot loaded
bayes_R2(ao_brm_model)

posterior_summary(ao_brm_model)

ao_brm_model_po <- brm(
  formula = bf(index ~ elevation + (1|species) + (1|plant.id)),
  data = ao.index,
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

saveRDS(ao_brm_model_po, file = "../brms_models/ao_brm_model_po.rds")
summary(ao_brm_model_po)
loo(ao_brm_model, ao_brm_model_po, moment_match = TRUE)
bayes_R2(ao_brm_model)
bayes_R2(ao_brm_model_po)

ao_brm_model_nb <- brm(
  formula = bf(index ~ elevation + (1|species) + (1|plant.id)),
  data = ao.index,
  family = negbinomial(),
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

saveRDS(ao_brm_model_nb, file = "../brms_models/ao_brm_model_nb.rds")

summary(ao_brm_model_nb)

pp_check(ao_brm_model_nb, type = "dens_overlay")
pp_check(ao_brm_model_nb, type = "hist")
pp_check(ao_brm_model_nb, type = "stat", stat = "mean")
pp_check(ao_brm_model_nb, type = "stat", stat = "sd")


loo(ao_brm_model, ao_brm_model_nb)



ao_brm_model_zinb <- brm(
  formula = bf(index ~ elevation + (1|species) + (1|plant.id),
               zi ~ elevation),
  data = ao.index,
  family = zero_inflated_negbinomial(),
  prior = c(
    set_prior("normal(0, 2)", class = "b"),          # elevation effect
    set_prior("student_t(3, 0, 5)", class = "sd"),
    set_prior("normal(0, 2)", class = "b", dpar = "zi")
  ),
  chains = 5, cores = 5,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  iter = 5000, warmup = 2000,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)

saveRDS(ao_brm_model_zinb, file = "../brms_models/ao_brm_model_zinb.rds")

summary(ao_brm_model_zinb)

pp_check(ao_brm_model_zinb, type = "dens_overlay")
pp_check(ao_brm_model_zinb, type = "hist")
pp_check(ao_brm_model_zinb, type = "stat", stat = "mean")
pp_check(ao_brm_model_zinb, type = "stat", stat = "sd")

loo(ao_brm_model_nb, ao_brm_model_zinb)

# Null model (no elevation effect on mean, still allow zi ~ elevation)
ao_brm_null_zinb <- brm(
  formula = bf(index ~ 1 + (1|species) + (1|plant.id),
               zi ~ elevation),
  data = ao.index,
  family = zero_inflated_poisson(),
  prior = c(
    set_prior("student_t(3, 0, 5)", class = "sd"),   # random effects
    set_prior("normal(0, 2)", class = "b", dpar = "zi")  # zero-inflation
  ),
  chains = 5, cores = 5,
  iter = 5000, warmup = 2000,
  control = list(max_treedepth = 15, adapt_delta = 0.99),
  seed = 1234
)

saveRDS(ao_brm_null_zinb, file = "../brms_models/ao_brm_null_zinb.rds")
ao_brm_null <- readRDS("../brms_models/ao_brm_null.rds")





















###GEITONOGAMY
go.index <- seed.indices %>%
  filter(treatment == "geitonogamy") %>%
  select(index, elevation, species, plant_number) %>%
  mutate(elevation = factor(elevation),
         species = factor(species)) %>%
  mutate(index = round(index * 100)) %>%
  mutate(plant.id = as.factor(paste0(elevation,species,plant_number)))%>%
  filter(species != "Lactuca i") %>%
  filter(species != "Hypericum r" | elevation != 3800)


go_brm_model_po <- brm(
  formula = bf(index ~ elevation + (1|species) + (1|plant.id)),
  data = go.index,
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

saveRDS(go_brm_model_po, file = "../brms_models/go_brm_model_po.rds")
summary(go_brm_model_po)
loo(go_brm_model, go_brm_model_po, moment_match = TRUE)
bayes_R2(go_brm_model)
bayes_R2(go_brm_model_po)


go_brm_model_zipo <- brm(
  formula = bf(index ~ elevation + (1|species) + (1|plant.id),
               zi ~ elevation),
  data = go.index,
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

saveRDS(go_brm_model_zipo, file = "../brms_models/go_brm_model_zipo.rds")
summary(go_brm_model_zipo)
loo(go_brm_model_po, go_brm_model_zipo, moment_match = TRUE)
bayes_R2(go_brm_model_po)
bayes_R2(go_brm_model_zipo)


go_brm_model_nb <- brm(
  formula = bf(index ~ elevation + (1|species) + (1|plant.id)),
  data = go.index,
  family = negbinomial(),
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

saveRDS(go_brm_model_nb, file = "../brms_models/go_brm_model_nb.rds")
summary(go_brm_model_nb)
loo(go_brm_model_pozi, go_brm_model_nb, moment_match = TRUE)
bayes_R2(go_brm_model_pozi)
bayes_R2(go_brm_model_nb)



go_brm_model_zinb <- brm(
  formula = bf(index ~ elevation + (1|species) + (1|plant.id),
               zi ~ elevation),
  data = go.index,
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

saveRDS(go_brm_model_zinb, file = "../brms_models/go_brm_model_zinb.rds")
summary(go_brm_model_zinb)
loo(go_brm_model_po, go_brm_model_zipo, go_brm_model_nb, go_brm_model_zinb)
bayes_R2(go_brm_model_po)
bayes_R2(go_brm_model_zipo)
bayes_R2(go_brm_model_nb)
bayes_R2(go_brm_model_zinb)


summary(go_brm_model_nb)
summary(go_brm_model_zinb)












###* For the pollen limitation data, I originally though that I would have to 
###* calculate it using bayesian statistics, however glmmTMB has a package 
###* "ordbetareg" which can simulate ordered beta regression without the need
###* for bayesian statistics
c.index <- c.index %>%
  filter(species != "Hypericum r" | elevation != 3800)
View(c.index)

glm_model3 <- glmmTMB(PL.index ~ elevation + (1|species) +(1|plant.id), 
                      data = c.index,
                      family = ordbeta())

glm_model3.null <- glmmTMB(PL.index ~ 1 + (1|species) +(1|plant.id), 
                           data = c.index,
                           family = ordbeta())

summary(glm_model3.null)
AIC(glm_model3, glm_model3.null)
emm <- emmeans(glm_model3, ~elevation)
pairs(emm)

saveRDS(glm_model3, "glm_outputs/pl_model.rds")
saveRDS(glm_model3.null, "glm_outputs/pl_null.rds")


summary(glm_model3)
simulationOutput <- simulateResiduals(fittedModel = glm_model3, plot = F)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)
testResiduals(simulationOutput)
plot(simulationOutput)

###* Create a dataset just for treatment "autogamy"

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

###* The model which best fits our data is a negative binomial with zero inflation

ao.po.zi <- glmmTMB(index ~ elevation + (1| species) + (1|plant.id),
                    ziformula=~elevation,
                    data = ao.index,
                    family = poisson)

ao.po.zi.null <- glmmTMB(index ~ 1 + (1| species) + (1|plant.id),
                         ziformula=~elevation,
                         data = ao.index,
                         family = poisson)

AIC(ao.po.zi, ao.po.zi.null)
emm <- emmeans(ao.po.zi, ~elevation)
pairs(emm)

saveRDS(ao.po.zi, "glm_outputs/ao_model.rds")
saveRDS(ao.po.zi.null, "glm_outputs/ao_null.rds")


summary(ao.nb.zi)
simulationOutputnbzi <- simulateResiduals(fittedModel = ao.nb.zi, plot = F)
testDispersion(simulationOutputnbzi)
testZeroInflation(simulationOutputnbzi)
testResiduals(simulationOutputnbzi)
plot(simulationOutputnbzi)

###* Create a dataset just for treatment "geitonogamy"

go.index <- seed.indices %>%
  filter(treatment == "geitonogamy") %>%
  select(index, elevation, species, plant_number) %>%
  mutate(elevation = factor(elevation),
         species = factor(species)) %>%
  mutate(index = round(index * 100)) %>%
  mutate(plant.id = as.factor(paste0(elevation,species,plant_number)))%>%
  filter(species != "Lactuca i") %>%
  filter(species != "Hypericum r" | elevation != 3800)

###* The model which best fits our data is a poisson with zero inflation

go.po.zi <- glmmTMB(index ~ elevation + (1| species) + (1|plant.id),
                    ziformula=~elevation,
                    data = go.index,
                    family = poisson)

go.po.zi.null <- glmmTMB(index ~ 1 + (1| species) + (1|plant.id),
                         ziformula=~elevation,
                         data = go.index,
                         family = poisson)

AIC(go.po.zi, go.po.zi.null)
emm <- emmeans(go.po.zi, ~elevation)
pairs(emm)

saveRDS(go.po.zi, "glm_outputs/go_model.rds")
saveRDS(go.po.zi.null, "glm_outputs/go_null.rds")

summary(go.po.zi)
simulationOutputpozi <- simulateResiduals(fittedModel = go.po.zi, plot = F)
testDispersion(simulationOutputpozi)
testZeroInflation(simulationOutputpozi)
testResiduals(simulationOutputpozi)
plot(simulationOutputpozi)

###*
###* PLOTTING THE DATA

#*THIS IS THE CORRECT CODE!!! BUT ONLY WHEN ZOOMED IN
summary_data <- ao.index %>%
  group_by(species, elevation) %>%
  summarise(count = n()) %>%
  ungroup()

View(seed.indices.autogamy)

# Plot the data
ao.index %>%
  ggplot(aes(x = as.factor(elevation), y = index, fill = species)) +
  geom_violin(trim = TRUE, scale = "width", width = 0.7,) +
  geom_boxplot(alpha = 0.1,
               outlier.color = "red",
               outlier.size = 2,
               outlier.shape = 8) +
  facet_wrap(. ~ species, nrow = 2, scales = "free_y",
             labeller = label_wrap_gen(width = 50)) +
  labs(x = "Elevation", y = "A/O Index") +
  theme_minimal() +
  geom_text(data = summary_data, 
            aes(x = as.factor(elevation), y = Inf, label = paste("n =", count)),
            position = position_dodge(width = 0.3),
            vjust = 0.9, hjust = 0.5, size = 5, inherit.aes = FALSE) +
  theme(axis.text = element_text( size = 15 ),
        axis.title = element_text( size = 15, face = "bold" ),
        strip.text = element_text(size = 15))

#*CONTROL
#*THIS IS THE CORRECT CODE!!! BUT ONLY WHEN ZOOMED IN
summary_data <- c.index %>%
  group_by(species, elevation) %>%
  summarise(count = n()) %>%
  ungroup()

# Plot the data
c.index %>%
  ggplot(aes(x = as.factor(elevation), y = seedset, fill = species)) +
  geom_violin(trim = TRUE, scale = "width", width = 0.7,) +
  geom_boxplot(alpha = 0.1,
               outlier.color = "red",
               outlier.size = 2,
               outlier.shape = 8) +
  facet_wrap(. ~ species, nrow = 2, scales = "free_y",
             labeller = label_wrap_gen(width = 50)) +
  labs(x = "Elevation", y = "Control index") +
  theme_minimal() +
  geom_text(data = summary_data, 
            aes(x = as.factor(elevation), y = Inf, label = paste("n =", count)),
            position = position_dodge(width = 0.3),
            vjust = 1, hjust = 0.5, size = 5, inherit.aes = FALSE) +
  theme(axis.text = element_text( size = 20 ),
        axis.title = element_text( size = 20, face = "bold" ),
        strip.text = element_text(size = 20))

#PL.INDEX
#*CONTROL
#*THIS IS THE CORRECT CODE!!! BUT ONLY WHEN ZOOMED IN
summary_data <- seed.indices.control %>%
  group_by(species, elevation) %>%
  summarise(count = n()) %>%
  ungroup()

# Plot the data
seed.indices.control %>%
  ggplot(aes(x = as.factor(elevation), y = PL.index, fill = species)) +
  geom_violin(trim = TRUE, scale = "width", width = 0.7,) +
  geom_boxplot(alpha = 0.1,
               outlier.color = "red",
               outlier.size = 2,
               outlier.shape = 8) +
  facet_wrap(. ~ species, nrow = 2, scales = "free_y",
             labeller = label_wrap_gen(width = 50)) +
  labs(x = "Elevation", y = "PL.index") +
  theme_minimal() +
  geom_text(data = summary_data, 
            aes(x = as.factor(elevation), y = Inf, label = paste("n =", count)),
            position = position_dodge(width = 0.3),
            vjust = 0.85, hjust = 0.5, size = 5, inherit.aes = FALSE) +
  theme(axis.text = element_text( size = 15 ),
        axis.title = element_text( size = 15, face = "bold" ),
        strip.text = element_text(size = 15))

#GEITONOGAMY
summary_data <- seed.indices.geitonogamy %>%
  group_by(species, elevation) %>%
  summarise(count = n()) %>%
  ungroup()


# Plot the data
seed.indices.geitonogamy %>%
  ggplot(aes(x = as.factor(elevation), y = index, fill = species)) +
  geom_violin(trim = TRUE, scale = "width", width = 0.7,) +
  geom_boxplot(alpha = 0.1,
               outlier.color = "red",
               outlier.size = 2,
               outlier.shape = 8) +
  facet_wrap(. ~ species, nrow = 2, scales = "free_y",
             labeller = label_wrap_gen(width = 50)) +
  labs(x = "Elevation", y = "G/O index") +
  theme_minimal() +
  geom_text(data = summary_data, 
            aes(x = as.factor(elevation), y = Inf, label = paste("n =", count)),
            position = position_dodge(width = 0.3),
            vjust = 1, hjust = 0.5, size = 5, inherit.aes = FALSE) +
  theme(axis.text = element_text( size = 15 ),
        axis.title = element_text( size = 15, face = "bold" ),
        strip.text = element_text(size = 15))



###*
###*
###*
###* 
###* 
###* 