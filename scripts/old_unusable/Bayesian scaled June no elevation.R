###* Bayesian scaled June
###* 
###* 
###* 
###* 
###* This script serves as running and interpreting the results from the best 
###* models created with the data which we have.
source("scripts/Bayesian setup June.R")

library(brms)

View(ao.final.table)
###* SEEDSET INDEX 
###* 
###* even though we are using Bayesian statistics, our dataset is
###* very small and so not we cannot get meaningfull results from it. This is 
###* because Hypericum produces hundreds of seeds, wheras the other plants 
###* produce significantly less

View(c.pl.final.table.4)
hist(log(c.pl.final.table.4$mean_seedset_round))
hist(c.pl.final.table.4$mean_seedset_round)

formula <- bf(mean_seedset_round | weights (seedset_weight_12) ~ 
                #1+
                me(mean.visited.flowers.scaled, sd.visited.flowers.scaled) +
                I(mean.visited.flowers.scaled^2) + 
                #me(mean.morpho.scaled, sd.morpho.scaled) +
                #me(mean.func.scaled, sd.func.scaled) + 
                (1|species))
###* Setting priors
prior <- c(
  set_prior("normal(0, 2)", class = "b"), 
  set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))

poly_bayesian_seedset_noel_w12_p2_v_h4_nb <- brm(
  formula = formula,
  data = c.pl.final.table.4,
  chains = 10,
  cores = 10,
  iter = 12000,
  warmup = 5000,
  family = negbinomial(),
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
















hist(log(c.pl.final.table.4$scaled_seedset))
hist(c.pl.final.table.4$scaled_log_seedset)
View(c.pl.final.table.4)

formula <- bf(scaled_seedset | weights (seedset_weight_12) ~ 
                #1+
                me(mean.visited.flowers.scaled, sd.visited.flowers.scaled) +
                #me(mean.morpho.scaled, sd.morpho.scaled) +
                #me(mean.func.scaled, sd.func.scaled) + 
                (1|species))
###* Setting priors
prior <- c(
  set_prior("normal(0, 2)", class = "b"), 
  set_prior("student_t(3, 0, 5)", class = "sd", group = "species"))

scaled_seedset_noel_w12_p2_v_h4 <- brm(
  formula = formula,
  data = c.pl.final.table.4,
  chains = 10,
  cores = 10,
  iter = 12000,
  warmup = 5000,
  family = student(),
  prior = prior,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  control = list(max_treedepth = 15, adapt_delta = 0.95),
  seed = 1234 # Addition to code
)


###* Trying version h4 version for seedset, which has info about Hypericum from
###* elevation 4000 / no produces seeds. Only usable in seedset though

saveRDS(bayesian_seedset_noel_w12_p1_m_h4, file = "../brms_models/bayesian_seedset_noel_w12_p1_v_h4.rds")
saveRDS(bayesian_seedset_noel_w12_p1_f_h4, file = "../brms_models/bayesian_seedset_noel_w12_p1_v_h4.rds")
saveRDS(bayesian_seedset_noel_w12_p1_null_h4, file = "../brms_models/bayesian_seedset_noel_w12_p1_null_h4.rds")
saveRDS(bayesian_seedset_noel_w12_p1_vmf_h4, file = "../brms_models/bayesian_seedset_noel_w12_p1_vmf_h4.rds")



saveRDS(scaled_log_seedset_noel_w12_p2_v_h4, file = "../brms_models/scaled_seedset_noel_w12_p2_v_h4.rds")



saveRDS(scaled_seedset_noel_w12_p2_v_h4, file = "../brms_models/scaled_seedset_noel_w12_p2_v_h4.rds")

bayesian_pl_noel_w12_p1_mf <- readRDS("../brms_models/bayesian_pl_noel_w12_p1_mf.rds")


















#me(mean.morpho.scaled, sd.morpho.scaled) +
#me(mean.func.scaled, sd.func.scaled) + 


# x = predictor with ME; sx = its row-wise SD
c.pl.final.table.4 <- c.pl.final.table.4 %>%
  mutate(
    x  = mean.visited.flowers.scaled,
    sx = sd.visited.flowers.scaled,
    x2 = x^2,
    sx2 = 2 * abs(x) * sx, # delta-method SD for x^2
    y  = mean.morpho.scaled,
    sy = sd.morpho.scaled,
    y2 = y^2,
    sy2 = 2 * abs(y) * sy,
    z  = mean.func.scaled,
    sz = sd.func.scaled,
    z2 = z^2,
    sz2 = 2 * abs(z) * sz
  )

form_quad_me <- bf(
  mean_seedset_round | weights(seedset_weight_12) ~
    me(x,  sx) +             # linear ME term
    me(x2, sx2) +            # quadratic ME term with propagated SD
    me(y,  sy) +             
    me(y2, sy2) +  
    me(z,  sz) +             
    me(z2, sz2) +  
    (1|species)
)

prior <- c(
  set_prior("normal(0, 2)", class = "b"),
  set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
)

fit_quad_me_po_vmf <- brm(
  form_quad_me,
  data = c.pl.final.table.4,
  family = poisson(),
  prior  = prior,
  chains = 10, cores = 10, iter = 12000, warmup = 5000,
  control = list(adapt_delta = 0.999, max_treedepth = 20),
  seed = 1234
)


summary(fit_quad_me)
summary(fit_quad_me_nb)
summary(fit_quad_me_po_vmf)

bayes_R2(fit_quad_me)
bayes_R2(fit_quad_me_nb)

loo(fit_quad_me, fit_quad_me_nb)

saveRDS(fit_quad_me, file = "../brms_models/fit_quad_me_po.rds")
saveRDS(fit_quad_me_nb, file = "../brms_models/fit_quad_me_nb.rds")
saveRDS(fit_quad_me_po_vmf, file = "../brms_models/fit_quad_me_po_vmf.rds")

# New dataset for prediction
newdat <- data.frame(
  x = seq(min(c.pl.final.table.4$mean.visited.flowers.scaled, na.rm = TRUE),
          max(c.pl.final.table.4$mean.visited.flowers.scaled, na.rm = TRUE),
          length.out = 100)
)
newdat$x2 <- newdat$x^2
newdat$sx <- mean(c.pl.final.table.4$sd.visited.flowers.scaled, na.rm = TRUE)  # use typical SD
newdat$sx2 <- 2 * abs(newdat$x) * newdat$sx   # delta method again
newdat$species <- NA  # placeholder

# Posterior predictions
pp <- posterior_epred(fit_quad_me, newdata = newdat, re_formula = NA)

# Summarise posterior
preds <- apply(pp, 2, function(x) c(fit_mean=mean(x), 
                                    lwr=quantile(x,0.025), 
                                    upr=quantile(x,0.975)))
preds <- as.data.frame(t(preds))

newdat <- bind_cols(newdat, preds)

names(newdat)[names(newdat) == "lwr.2.5%"] <- "lwr"
names(newdat)[names(newdat) == "upr.97.5%"] <- "upr"

# Plot
ggplot(newdat, aes(x, fit_mean)) +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_point(data = c.pl.final.table.4,
             aes(x = mean.visited.flowers.scaled,
                 y = mean_seedset_round),
             inherit.aes = FALSE,
             color = "black") +
  labs(x = "Scaled visitation", y = "Predicted seed set") +
  theme_minimal()








c.pl.final.table.4 <- c.pl.final.table.4 %>%
  mutate(
    x   = mean.visited.flowers.scaled,
    sx  = sd.visited.flowers.scaled,
    x2  = x^2,
    sx2 = 2 * abs(x) * sx   # delta method SD for x^2
  )

form_quad_log <- bf(
  scaled_log_seedset | weights(seedset_weight_12) ~ 
    me(x, sx) + 
    me(x2, sx2) + 
    (1|species)
)

priors <- c(
  set_prior("normal(0, 2)", class = "b"),
  set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
)

fit_quad_log_student <- brm(
  formula = form_quad_log,
  data = c.pl.final.table.4,
  family = student(),
  prior = priors,
  chains = 6, cores = 6, iter = 6000, warmup = 2000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 1234
)

summary(fit_quad_log)

saveRDS(fit_quad_log, file = "../brms_models/fit_quad_log.rds")
saveRDS(fit_quad_log_student, file = "../brms_models/fit_quad_log_student.rds")







































































###* New 6.9.2025
###* Trying to fit model - it doesn't make sense for visitation and seedset to 
###* not be correlated
saveRDS(bayesian_seedset_noel_w12_p2_v_h4, file = "../brms_models/bayesian_seedset_noel_w12_p2_v_h4.rds")
summary(bayesian_seedset_noel_w12_p2_v_h4)
# RHAT very high, Blk and Tail ESS both low.
# Does more chains even make sense? THis is more than normal..

# Try negative binomial, amybe that could help
saveRDS(bayesian_seedset_noel_w12_p2_v_h4_nb, file = "../brms_models/bayesian_seedset_noel_w12_p2_v_h4_nb.rds")
summary(bayesian_seedset_noel_w12_p2_v_h4_nb)

###* Or I do not know, if I have tried to log the response. I will try to do this next.
###* It would make sense, since all the other variables are somehow transformed




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








library(ordbetareg)

c.pl.final.table.5$intercept_only <- 1

ordbeta_pl_noel_w12_p1_null <- ordbetareg(
  formula = bf(
    mean_PL_index | weights(PL_index_weight_12) ~ 
      intercept_only +
      #me(mean.visited.flowers.scaled, sd.visited.flowers.scaled) +
      #me(mean.morpho.scaled, sd.morpho.scaled) +
      #me(mean.func.scaled, sd.func.scaled) + 
      (1|species)
  ),
  data = c.pl.final.table.5,
  chains = 10,
  cores = 10,
  iter = 12000,
  warmup = 5000,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)





#save




#saved
saveRDS(ordbeta_pl_noel_w12_p1_vmf, file = "../brms_models/ordbeta_pl_noel_w12_p1_vmf.rds")
saveRDS(ordbeta_pl_noel_w12_p1_v, file = "../brms_models/ordbeta_pl_noel_w12_p1_v.rds")
saveRDS(ordbeta_pl_noel_w12_p1_vm, file = "../brms_models/ordbeta_pl_noel_w12_p1_vm.rds")
saveRDS(ordbeta_pl_noel_w12_p1_m, file = "../brms_models/ordbeta_pl_noel_w12_p1_m.rds")
saveRDS(ordbeta_pl_noel_w12_p1_f, file = "../brms_models/ordbeta_pl_noel_w12_p1_f.rds")
saveRDS(ordbeta_pl_noel_w12_p1_mf, file = "../brms_models/ordbeta_pl_noel_w12_p1_mf.rds")
saveRDS(ordbeta_pl_noel_w12_p1_null, file = "../brms_models/ordbeta_pl_noel_w12_p1_null.rds")

loo(ordbeta_pl_noel_w12_p1_vmf, ordbeta_pl_noel_w12_p1_v, ordbeta_pl_noel_w12_p1_vm,
    ordbeta_pl_noel_w12_p1_m, ordbeta_pl_noel_w12_p1_f, ordbeta_pl_noel_w12_p1_mf,
    ordbeta_pl_noel_w12_p1_null)

bayes_R2(ordbeta_pl_noel_w12_p1_null)
bayes_R2(ordbeta_pl_noel_w12_p1_v)
bayes_R2(ordbeta_pl_noel_w12_p1_m)
bayes_R2(ordbeta_pl_noel_w12_p1_f)
bayes_R2(ordbeta_pl_noel_w12_p1_vmf)
bayes_R2(ordbeta_pl_noel_w12_p1_vm)
bayes_R2(ordbeta_pl_noel_w12_p1_mf)

summary(ordbeta_pl_noel_w12_p1_vmf)

View(c.pl.final.table.5)




























library(ordbetareg)

###* POLLEN LIMITATION INDEX
c.pl.final.table.5 <- c.pl.final.table.4 %>%
  filter(species != "Hypericum r" | elevation != 4000)

# x = predictor with ME; sx = its row-wise SD
c.pl.final.table.5 <- c.pl.final.table.4 %>%
  mutate(
    x  = mean.visited.flowers.scaled,
    sx = sd.visited.flowers.scaled,
    x2 = x^2,
    sx2 = 2 * abs(x) * sx, # delta-method SD for x^2
    y  = mean.morpho.scaled,
    sy = sd.morpho.scaled,
    y2 = y^2,
    sy2 = 2 * abs(y) * sy,
    z  = mean.func.scaled,
    sz = sd.func.scaled,
    z2 = z^2,
    sz2 = 2 * abs(z) * sz
  )

form_quad_me_ord <- bf(
  mean_PL_index | weights(PL_index_weight_12) ~
    me(x,  sx) +             # linear ME term
    me(x2, sx2) +            # quadratic ME term with propagated SD
    me(y,  sy) +             
    me(y2, sy2) +  
    me(z,  sz) +             
    me(z2, sz2) +  
    (1|species)
)

ordbeta_pl_noel_w12_p1_vmf_quad <- ordbetareg(
  formula = form_quad_me_ord,
  data = c.pl.final.table.5,
  chains = 10,
  cores = 10,
  iter = 12000,
  warmup = 5000,
  save_pars = save_pars(latent = TRUE, all = TRUE),
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)

summary(ordbeta_pl_noel_w12_p1_vmf_quad)
saveRDS(ordbeta_pl_noel_w12_p1_vmf_quad, file = "../brms_models/ordbeta_pl_noel_w12_p1_vmf_quad.rds")


























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
                #me(mean.morpho.scaled, sd.morpho.scaled) +
                me(mean.func.scaled, sd.func.scaled) +
                (1|species))

bayesian_ao_noel_w12_p2_f <- brm(
  formula = formula,
  data = ao.final.table,
  family = zero_inflated_beta(),
  prior = c(
    set_prior("normal(0, 2)", class = "b"),
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








###* What happens, when priors are set more flat?
saveRDS(bayesian_ao_noel_w12_p2_f, file = "../brms_models/bayesian_ao_noel_w12_p2_f.rds")
summary(bayesian_ao_noel_w12_p2_f)
fit <- bayesian_ao_noel_w12_p2_f

















# x = predictor with ME; sx = its row-wise SD
ao.final.table <- ao.final.table %>%
  mutate(
    x  = mean.visited.flowers.scaled,
    sx = sd.visited.flowers.scaled,
    x2 = x^2,
    sx2 = 2 * abs(x) * sx, # delta-method SD for x^2
    y  = mean.morpho.scaled,
    sy = sd.morpho.scaled,
    y2 = y^2,
    sy2 = 2 * abs(y) * sy,
    z  = mean.func.scaled,
    sz = sd.func.scaled,
    z2 = z^2,
    sz2 = 2 * abs(z) * sz
  )

formula <- bf(
  mean_ao_index | weights(ao_index_weight_12) ~
    #me(x,  sx) +             # linear ME term
    #me(x2, sx2) +            # quadratic ME term with propagated SD
    #me(y,  sy) +             
    #me(y2, sy2) +  
    me(z,  sz) +             
    #me(z2, sz2) +  
    (1|species)
)

View(ao.final.table)

bayesian_ao_noel_w12_p2_f_only_quad <- brm(
  formula = formula,
  data = ao.final.table,
  family = zero_inflated_beta(),
  prior = c(
    set_prior("normal(0, 2)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 10,
  cores = 10,
  iter = 12000,
  warmup = 5000,
  control = list(max_treedepth = 20, adapt_delta = 0.999),
  seed = 1234
)

saveRDS(bayesian_ao_noel_w12_p2_vmf_quad, file = "../brms_models/bayesian_ao_noel_w12_p2_vmf_quad.rds")
bayesian_ao_noel_w12_p2_vmf_quad <- readRDS("../brms_models/bayesian_ao_noel_w12_p2_vmf_quad.rds")
summary(bayesian_ao_noel_w12_p2_vmf_quad)

#without quadratic term
saveRDS(bayesian_ao_noel_w12_p2_f_only_quad, file = "../brms_models/bayesian_ao_noel_w12_p2_f_only_quad.rds")
bayesian_ao_noel_w12_p2_f_only_quad <- readRDS("../brms_models/bayesian_ao_noel_w12_p2_f_only_quad.rds")
summary(bayesian_ao_noel_w12_p2_f_only_quad)

#with quadratic term
saveRDS(bayesian_ao_noel_w12_p2_f_quad, file = "../brms_models/bayesian_ao_noel_w12_p2_f_quad.rds")
bayesian_ao_noel_w12_p2_f_quad <- readRDS("../brms_models/bayesian_ao_noel_w12_p2_f_quad.rds")
summary(bayesian_ao_noel_w12_p2_f_quad)
plot(conditional_effects(bayesian_ao_noel_w12_p2_f_quad), points = TRUE)

loo(bayesian_ao_noel_w12_p2_f, bayesian_ao_noel_w12_p2_f_quad)






library(dplyr)
library(ggplot2)
library(brms)

## 1) Range of the predictor actually used in the fit (z)
# If you still have 'z' in ao.final.table:
rng <- range(ao.final.table$z, na.rm = TRUE)

# If you don't have 'z' saved, use the original column directly:
# rng <- range(ao.final.table$mean.func.scaled, na.rm = TRUE)

## 2) Build newdata with the SAME variable names the models expect
newdat <- tibble(
  z  = seq(rng[1], rng[2], length.out = 200)
) %>%
  mutate(
    z2 = z^2,                                         # exactly as in the model
    sz = median(ao.final.table$sz, na.rm = TRUE),     # typical ME SD for z
    sz2 = 2 * abs(z) * sz,                            # delta-method SD for z^2
    species = NA                                      # population-level curve
  )

## 3) Population-level predictions (no random effects)
pp_lin  <- posterior_epred(bayesian_ao_noel_w12_p2_f_only_quad,
                           newdata = newdat, re_formula = NA)
pp_quad <- posterior_epred(bayesian_ao_noel_w12_p2_f_quad,
                           newdata = newdat, re_formula = NA)

## 4) Summarize posterior predictions
summ_from_pp <- function(pp) {
  s <- apply(pp, 2, function(x) c(
    fit_mean = mean(x),
    lwr = quantile(x, 0.025),
    upr = quantile(x, 0.975)
  ))
  as.data.frame(t(s))
}

pred_lin  <- bind_cols(newdat, summ_from_pp(pp_lin))  %>% mutate(model = "linear")
pred_quad <- bind_cols(newdat, summ_from_pp(pp_quad)) %>% mutate(model = "quadratic")
pred_all  <- bind_rows(pred_lin, pred_quad)

## 5) Raw data to overlay (use the same predictor z used in model)
raw_df <- ao.final.table %>%
  transmute(z = z %||% mean.func.scaled,  # if z doesn't exist, fall back
            mean_ao_index = mean_ao_index)

## 6) Plot combined effects (mean curve + 95% band)
ggplot(pred_all, aes(z, fit_mean, color = model, fill = model)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, color = NA) +
  geom_line(size = 1) +
  geom_point(data = raw_df,
             aes(x = z, y = mean_ao_index),
             inherit.aes = FALSE, color = "black", alpha = 0.6) +
  labs(x = "Functional groups (scaled z)",
       y = "AO index (population-level fit)",
       title = "Combined effect of functional groups on AO index") +
  theme_minimal() +
  theme(legend.position = "top")


head(pred_all)

pred_all <- pred_all %>%
  mutate(
    lwr = `lwr.2.5%`,
    upr = `upr.97.5%`
  )

loo_lin  <- loo(bayesian_ao_noel_w12_p2_f_only_quad)
loo_quad <- loo(bayesian_ao_noel_w12_p2_f_quad)
loo_compare(loo_lin, loo_quad)

k10_lin  <- kfold(bayesian_ao_noel_w12_p2_f_only_quad, K = 2)
k10_quad <- kfold(bayesian_ao_noel_w12_p2_f_quad, K = 10)
kfold_compare(k10_lin, k10_quad)


















# Which entries of the pointwise log-likelihood are non-finite?
ll_quad <- log_lik(bayesian_ao_noel_w12_p2_f_quad, pointwise = TRUE)
bad_draw  <- !is.finite(ll_quad) & !is.infinite(ll_quad)  # TRUE where NaN or +Inf
any(bad_draw)  # should be TRUE given your error

# Count how many bad entries per observation (columns are obs)
colSums(bad_draw)

# Check if your mu’s ever go numerically to 0 or 1 for those obs
mu_quad <- posterior_linpred(bayesian_ao_noel_w12_p2_f_quad,
                             transform = TRUE, re_formula = NA)
range(mu_quad)                  # should be strictly (0,1); if not, that’s a clue


summary(ao.final.table$sz)
summary(ao.final.table$sz2)
any(!is.finite(ao.final.table$sz2))























priors_quad <- c(
  set_prior("normal(0, 1)",   class = "b", coef = "mezsz"),     # linear
  #set_prior("normal(0, 0.5)", class = "b", coef = "mez2sz2"),   # quadratic tighter
  set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
)

eps <- 1e-8
ao.final.table <- ao.final.table %>%
  mutate(
    sz  = pmax(sz,  eps),
    sz2 = pmax(2 * abs(z) * sz, eps)   # delta-method with tiny floor
  )

bayesian_ao_noel_w12_p1_f_only_quad <- brm(
  mean_ao_index | weights(ao_index_weight_12) ~ 
    me(z, sz) + 
    #me(z2, sz2) + 
    (1|species),
  data = ao.final.table,
  family = zero_inflated_beta(),                    # keep family (no 1s in data)
  prior  = priors_quad,
  save_pars = save_pars(all = TRUE, latent = TRUE), # <-- critical for LOO with me()
  control   = list(adapt_delta = 0.995),            # a bit higher to avoid extremes
  chains = 8, cores = 8, iter = 6000, warmup = 2500, seed = 1234
)

saveRDS(bayesian_ao_noel_w12_p1_f_quad, file = "../brms_models/bayesian_ao_noel_w12_p1_f_quad.rds")
bayesian_ao_noel_w12_p1_f_quad <- readRDS("../brms_models/bayesian_ao_noel_w12_p1_f_quad.rds")
summary(bayesian_ao_noel_w12_p1_f_quad)


saveRDS(bayesian_ao_noel_w12_p1_f_only_quad, file = "../brms_models/bayesian_ao_noel_w12_p1_f_only_quad.rds")
bayesian_ao_noel_w12_p1_f_only_quad <- readRDS("../brms_models/bayesian_ao_noel_w12_p1_f_only_quad.rds")
summary(bayesian_ao_noel_w12_p1_f_only_quad)

loo_f <- loo(bayesian_ao_noel_w12_p1_f_only_quad, moment_match = TRUE)
loo_f <- loo(bayesian_ao_noel_w12_p1_f_only_quad)
loo_quad <- loo(bayesian_ao_noel_w12_p1_f_quad, moment_match = TRUE)
loo_quad <- loo(bayesian_ao_noel_w12_p1_f_quad)

loo_compare(loo_f, loo_quad)













































ao.final.table <- ao.final.table %>%
  mutate(
    x  = mean.visited.flowers.scaled,
    sx = sd.visited.flowers.scaled,
    x2 = x^2,
    sx2 = 2 * abs(x) * sx, # delta-method SD for x^2
    y  = mean.morpho.scaled,
    sy = sd.morpho.scaled,
    y2 = y^2,
    sy2 = 2 * abs(y) * sy,
    z  = mean.func.scaled,
    sz = sd.func.scaled,
    z2 = z^2,
    sz2 = 2 * abs(z) * sz
  )

summary(ao.final.table$x)
summary(ao.final.table$sx)
summary(ao.final.table$sx2)
summary(ao.final.table$sy)
summary(ao.final.table$sy2)
summary(ao.final.table$sz)
summary(ao.final.table$sz2)


priors_quad <- c(
  set_prior("normal(0, 1)",   class = "b", coef = "mezsz"),     
  set_prior("normal(0, 0.5)", class = "b", coef = "mez2sz2"),   
  set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
)

eps <- 1e-8
ao.final.table <- ao.final.table %>%
  mutate(
    sz  = pmax(sz,  eps),
    sz2 = pmax(2 * abs(z) * sz, eps)   # delta-method with tiny floor
  )


bayesian_ao_noel_w12_p1_vmf_quad <- brm(
  mean_ao_index | weights(ao_index_weight_12) ~ 
    me(x,  sx) + 
    me(x2, sx2) + 
    me(y,  sy) +             
    me(y2, sy2) +  
    me(z,  sz) +             
    me(z2, sz2) +  
    (1|species),
  data = ao.final.table,
  family = zero_inflated_beta(),                    # keep family (no 1s in data)
  prior  = priors_quad,
  save_pars = save_pars(all = TRUE, latent = TRUE), # <-- critical for LOO with me()
  control   = list(adapt_delta = 0.995),            # a bit higher to avoid extremes
  chains = 8, cores = 8, iter = 6000, warmup = 2500, seed = 1234
)





























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
