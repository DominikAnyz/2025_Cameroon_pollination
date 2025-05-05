
source("scripts/setup.R")
source("scripts/Bayesian setup.R")

###* The brms models saves are too big for uplaoding to github, they can be found here:
###* https://drive.google.com/drive/folders/1edMD7uW1J__-lk-HkRRdAHg0w6QoqzVO?usp=sharing

set.seed(1234)

Sys.setenv(TZ = "UTC")

###* ---CONTROL SEEDSET---
###* Rob sent me this formula from Ondra
###* formula =  bf(R_mean | se(R_sd) ~ Elevation + me(A_mean, A_sd) + me(B_mean, B_sd)+ me(C_mean, C_sd)
###* 
###* The formula by Ondra would let us measure the uncertainty both in our response,
###* the index, and in our predictor, visitation rate.
###* 
###* The formula can only be used for the gaussian distribution, so perhaps for 
###* seedset, but for the other indexes, we will have to find a workaround.
###* 
###* For the seedset, in this formula, the response is expected to have a normal 
###* distribution. The "R_mean | se(R_sd)" takes into account not the standard 
###* deviation as I originally thought, but the standard error. That is why,
###* SE was added in script "Bayesian setup.R" on line 470.
###* 
###* For the response, in the fomula, me(A_mean, A_sd), me means measurement 
###* error and the error now really does use the standard deviation.
###* column sd.visitors.per.minute was created in "Bayesian setup" on line 345
###* 
###* We so not have standard error for the elevation or species.
###* 
###* What should our model look like? We are asking what is the effect of 
###* the index on visitation. In the glm's we were asking if the indexes were being 
###* affected by elevation, but now we are asking something different. In this case,
###* I think that both the elevation and species should act as a random effect. 
###* And since we are using the means of most of the other columns, it does not
###* make sense to me to add any other predictor than visitation. 
###* 
###* We cannot add most common functional group, since it is mostly "Other fly". 
###* This makes sense for it to be the most common. We cannot measure species
###* richness, since not all species are identified. We could theoretically look
###* at species richness of species identified to morphospecies level.. We could 
###* theoretically also look at the effect of total amount of functionalo groups.
###* 
###* First however, I would likt to try a simple formula and see how it goes:
###* formula <- bf(mean_seedset | se(se_mean_seedset) ~ me(visitors, sd_visitors) + (1|species) + (1|elevation))
###* 
###* The first step: visualize data:
View(c.pl.final.table)

###* Delete rows with missing values (some species were treated for the pollination
###* experiment, but not filmed), make sure that elevation and species are factors
c.pl.final.table <- na.omit(c.pl.final.table)
c.pl.final.table$elevation <- as.factor(c.pl.final.table$elevation)
c.pl.final.table$species <- as.factor(c.pl.final.table$species)

###* Histograms
hist(c.pl.final.table$mean_seedset)
hist(log(c.pl.final.table$mean_seedset))
hist(c.pl.final.table$mean.visitors.per.minute)
hist(log(c.pl.final.table$mean.visitors.per.minute))

###* Both the historgam of mean seedset and visitation or skewed and look to have
###* a more normal distribution when loged
###*
###* We have to alter all: mean seedset, standard error of the seedset, \
###* mean visitors per minute and the standard error of visitors per minute
###* 
###* We will makes logs of the seedset and visitors, we will use the delta method 
###* on the se of seedset and sd visitors
c.pl.final.table <- c.pl.final.table %>%
  mutate(log_mean_seedset = log(mean_seedset),
         se_log_mean_seedset = se_seedset / mean_seedset,
         log_visitors = log(mean.visitors.per.minute),
         sd_log_visitors = sd.visitors.per.minute / mean.visitors.per.minute
  )

###* Check
hist(c.pl.final.table$log_mean_seedset)
hist(c.pl.final.table$log_visitors)

###* Set up formula
formula <- bf(log_mean_seedset | se(se_log_mean_seedset) ~ me(log_visitors, sd_log_visitors) + (1|species) + (1|elevation))

###* Setting priors
prior <- c(
  set_prior("normal(0, 1)", class = "b"),
  set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
  set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
)

###* Run model
fit <- brm(
  formula = formula,
  data = c.pl.final.table,
  chains = 5,
  cores = 10,
  iter = 6000,
  prior = prior,
  control = list(max_treedepth = 30, adapt_delta = 0.999)
)
#fit <- readRDS("../brms_models/brms_seedset.rds")

summary(fit)
###* The model summary tells us:
###* -visitation rate is a significant predictor of seedset (loglog scale)
###* -both species and elevation play a role in the seedset, but species more
###* -number of observations is low (17), which is why we are using SE and so
###* the sigma is 0
###* -Rhat=1, meaning good convergance
###* We can took at the relationship of visitors on seedset
conditional_effects(fit)

###* What do posterior predictions tell us?
pp_check(fit)
###* It seems that predicted vs observed data seem pretty similar, however the
###* tails are a bit misalligned -> check residuals
###* 
###* Extract fitted values and residuals, plot them
fitted_values <- fitted(fit)
residuals <- residuals(fit)
plot(fitted_values, residuals)
abline(h = 0, col = "red")
###* It seems that the residuals could be better. However, would we not expect 
###* such residuals with our model? Can we do anything about it? Yes, we can make
###* a simpler model, but since we cannot get a simpler model, we cannot do much
###* about it.
###* Histogram of residuals
hist(residuals, main = "Residuals Histogram", xlab = "Residuals")
###* Althought the residuals could be nicer, they have a normal distribution.
###* I would consider this as a well fit model
###* IS THIS CORRECT?
saveRDS(fit, file = "../brms_models/brms_seedset.rds")










###* --- POLLEN LIMITATION---
###* So we have some results for seedset, how about pollen limitation (PL)?
###*
###* The first problem with pollen limitation is that we are now looking at a
###* beta distribution
View(c.pl.final.table)
###* Fortunately, it seems that no values are exactly 0 or 1, so we are actually
###* looking at a pretty normal beta distribution. OR?
hist(c.pl.final.table$mean_PL_index)
###* Could be better, could be worse. However, since we will be using a the beta 
###* distibution, there is no longer an assumption of a normal distribution.
###* 
###* Unfortunately, there is also the ability to add uncertainty into the response, 
###* since it is modeled differently than the gaussian distribution. The documentation
###* of brms tells us this about beta distribution:
###* 
###* It seems that we can try a few variations: The beta family does support 
###* measurement error in the predictor, but not in the response. This is 
###* because the variance is implicitly tied to the mean via its' formula. So we
###* can only model the mean and the preciscion (phi)
###* 
###* At any rate, in PL we will not be using the SE or SD of the mean of our response
###* 
###* We can try the mean PL index, without the "measurement error". For better 
###* visualization, we will keep the data about the visitors in the log scale
formula <- bf(mean_PL_index ~ me(log_visitors, sd_log_visitors) + (1|species) + (1|elevation))

###* Setting priors
prior <- c(
  set_prior("normal(0, 2.5)", class = "b"),
  set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
  set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
)

fit_PL_orig <- brm(
  formula = formula,
  data = c.pl.final.table,
  chains = 6,
  family = Beta(),
  cores = 10,
  iter = 10000,
  warmup = 5000,
  prior = prior,
  control = list(max_treedepth = 20, adapt_delta = 0.9999)
)
saveRDS(fit_PL_orig, file = "../brms_models/brms_PL.rds")
#fit_PL <- readRDS("../brms_models/brms_PL.rds")

###* When rerunning the code, I realized that I was getting different results, so
###* I tried this experiment, fitting multiple same models to see results.
###* In the middle, I thought that maybe I had not set the seed properly.
#fit_PL_orig - 10k iter, 5k warm, 0 divergent, not significant
#fit_PL_orig_2 -10k iter, 5 k warm, 1 divergent, significant
#fit_PL_orig_3 - 12k iter, 6k warm, 3 divergent, significant
#fit_PL_orig_4 - 12k iter, 6k warm, 2 divergent, significant
#set.seed(1234)
#fit_PL_orig_5 - 12k iter, 6k warm, 0 divergent, significant
#fit_PL_orig_6 - 12k iter, 6k warm, 1 divergent, significant
#fit_PL_orig_7 - 12k iter, 6k warm, 1 divergent, significant
###* 
###* Apparantly there is a problem with setting seed, which I did not realize 
###* before. However, even with identical models, I am unable to get the same 
###* result multiple times
###* 
###* So apparently seed is set in the model in brms...
fit_PL_orig_8 <- brm(
  formula = formula,
  data = c.pl.final.table,
  chains = 6,
  family = Beta(),
  cores = 10,
  iter = 12000,
  warmup = 6000,
  prior = prior,
  control = list(max_treedepth = 20, adapt_delta = 0.9999),
  seed = 1234 # Addition to code
)
saveRDS(fit_PL_orig_8, file = "../brms_models/brms_PL_orig_8.rds")
#fit_PL_orig_8 <- readRDS("../brms_models/brms_PL_orig_8.rds")

###* When I set the seed I was originally using, I get convergence and 
###* significance. However, now that I have tried many models I think I should 
###* be careful with these results.
###* 
###* However, we still have to analyze them:
summary(fit_PL_orig_8)
###* The model summary tells us:
###* -There is a negative effect of log_visitors on mean_PL_index - as (log) 
###* visitation rate increases, pollen limitation decreases. This would make sense.
###* However, as we have seen, not all models which have the same parameters get 
###* to the same results
###* -The intercept's confidence interval contains 0, suggesting uncertainty in 
###* the baseline estimate.
###* -Both the species and elevation play a role in the variation in the PL index,
###* as expected
###* -ESS is high, so posteriors have been explored well
conditional_effects(fit_PL_orig_8)
pairs(fit_PL_orig_8)

###* What can the pp tell us?
pp_check(fit_PL_orig_8, nsamples = 100)
###* We can see a lot of noise even in the middle, where we would expect more 
###* similar patterns
###* 
###* Check residuals
fitted_values <- fitted(fit_PL_orig_8)
residuals <- residuals(fit_PL_orig_8)
plot(fitted_values, residuals)
abline(h = 0, col = "red")
###* Oof, that does not look too good. But again, there is not much we can do 
###* about our residuals. I do not know if we can call this a good model. On the
###* other hand, would we not expect the residauls to be heteroscandic? 
kfold(fit_PL_orig_8)



###* Are there other ways which we could use to try and measure the data? I would
###* like to try and find something, which would take into account the measuremement
###* error, same as with the seedset.
###*
###* A way to try it would be to try a weighted beta regression, to assign 
###* more weight to groups with more pollen limitation replicates. 
###* 
###* For this, line 475 was added to "Bayesian setup.R" to create colum "n_replicates"
formula_weighted <- bf(mean_PL_index | weights(n_replicates) ~ me(log_visitors, sd_log_visitors) + (1|species) + (1|elevation))

fit_weighted_full <- brm(
  formula = formula_weighted,
  data = c.pl.final.table,
  family = Beta(), 
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 6,
  cores = 10,
  iter = 4000,
  warmup = 1000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)

saveRDS(fit_weighted_full, file = "../brms_models/brms_PL_weights.rds")
#fit_weighted_full <- readRDS("../brms_models/brms_PL_weights.rds")

summary(fit_weighted_full)
###* - This model tells me something different than what I would expect: it tells 
###* me that with increasing visitation, the pollen limitation index also increases.
###* However, the CI spans from -1.18 to 1.26, meaning that this cannot be viewed
###* as statistically significant.
###* - The model suggests that elevation and species have a much higher effect on 
###* the PL_index, both with CI not including zero, meaning statistical significance
###* - Rhat values are higher than 1, suggesting convergance issues
pp_check(fit_weighted_full)
###* 
###* 
###* 
###*
###* I would like to try the model without log transformed visitors, although
###* it seems to me, taht it make more sense to have the model use the log 
###* transformed version
formula_weighted_2 <- bf(mean_PL_index | weights(n_replicates) ~ me(mean.visitors.per.minute, sd.visitors.per.minute) + (1|species) + (1|elevation))

fit_weighted_full_2 <- brm(
  formula = formula_weighted_2,
  data = c.pl.final.table,
  family = Beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 6,
  cores = 10,
  iter = 4000,
  warmup = 1000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)

saveRDS(fit_weighted_full_2, file = "../brms_models/brms_PL_weights_2.rds")
#fit_weighted_full_2 <- readRDS("../brms_models/brms_PL_weights_2.rds")

summary(fit_weighted_full_2)
###* The results of this model are similar to those of the first weighted model
###* - Again, we have the results that with increasing visitation comes increasing
###* pollen limitation, but again it is not statistically significant
###* - Again, elevation and species have a statistically significant effect on the 
###* PL index.
###* Again, Rhat is higher than 1, (although less than previous model), suggesting
###* convergence issues
###* 
###* 
###* Without writting the whole model code, I would like to try different things:
###* running the models without the sd, of visitors both normal and loged, just to 
###* see what this would do to results
formula_fit_weighted_partial <- bf(mean_PL_index | weights(n_replicates) ~ log_visitors + (1 | species) + (1 | elevation))
saveRDS(fit_weighted_partial, file = "../brms_models/brms_PL_weighted_partial.rds")
fit_weighted_partial <- readRDS("../brms_models/brms_PL_weighted_partial.rds")
summary(fit_weighted_partial)

formula_fit_weighted_partial_2 <- bf(mean_PL_index | weights(n_replicates) ~ mean.visitors.per.minute + (1|species) + (1|elevation))
saveRDS(fit_weighted_partial, file = "../brms_models/brms_PL_weighted_partial_2.rds")
fit_weighted_partial_2 <- readRDS("../brms_models/brms_PL_weighted_partial_2.rds")
summary(fit_weighted_partial_2)

###* Both models are better converged than the "full" versions, but show an opposite
###* trend. The log visitors show a very slight decline in PL with increasing 
###* visitation, wheras the non-transformed means should and incline. However, 
###* none of these are statistically significant
###* 
###* 
###* AT ANY RATE, it is correct to use weights? Are they not affecting the the
###* predictor also? 



###* However, in my opinion, not using the log
###* scale for the predictor when possible is a shame, since the relationshiop 
###* between pollination and seedset doesn't necessarilly need to be linear.




###* Trying to increase iterations to lower Rhat
formula_weighted <- bf(mean_PL_index | weights(n_replicates) ~ me(log_visitors, sd_log_visitors) + (1|species) + (1|elevation))

fit_weighted_full_more <- brm(
  formula = formula_weighted,
  data = c.pl.final.table,
  family = Beta(), 
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 6,
  cores = 10,
  iter = 8000,
  warmup = 3000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)

summary(fit_weighted_full_more)
saveRDS(fit_weighted_full_more, file = "../brms_models/brms_PL_weighted_full_more.rds")
fit_weighted_full_more <- readRDS("../brms_models/brms_PL_weighted_full_more.rds")

formula_weighted_2 <- bf(mean_PL_index | weights(n_replicates) ~ me(mean.visitors.per.minute, sd.visitors.per.minute) + (1|species) + (1|elevation))

fit_weighted_full_2 <- brm(
  formula = formula_weighted_2,
  data = c.pl.final.table,
  family = Beta(),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 6,
  cores = 10,
  iter = 8000,
  warmup = 3000,
  control = list(max_treedepth = 20, adapt_delta = 0.99),
  seed = 1234
)


summary(fit_weighted_full_2)
saveRDS(fit_weighted_full_2, file = "../brms_models/brms_PL_weighted_full_more_2.rds")
fit_weighted_full_more_2 <- readRDS("../brms_models/brms_PL_weighted_full_more_2.rds")







###* Ondra recomendations:
###* 
###* RESPONSE
###* 1) Instead of n_replicates, created index_weights: before making means of index,
###* do 1 / sd
###* 
###* PRESICTOR
###* 1) Try scaling the predictor. Function scale
###* 2) Try log + 1 on the predictor
###* 
###* 
###* 
###* How to test validity of models?
###* 
###* summary()
###* pp_check()
###* bayes_R2()
###* plot()
###* ppc_residuals()
###* 
###* loo() - propbably will not work due to low repetitions
###* fitted_values()
###* 


summary()
pp_check()
bayes_R2()
plot()
ppc_residuals()









formula <- bf(mean_PL_index ~ me(mean.visitors.per.minute, sd.visitors.per.minute) +  (1|species) + (1|elevation)) 
#              phi ~  (1|species) + (1|elevation))






###*
###*
###*
###*
###*
###*
###*
###* 
###* --- AUTOGAMY ---
View(ao.final.table)
hist(ao.final.table$mean_ao_index)
hist(log(ao.final.table$mean_ao_index_trans))

ao.final.table <- na.omit(ao.final.table)
ao.final.table$elevation <- as.factor(ao.final.table$elevation)
ao.final.table$species <- as.factor(ao.final.table$species)


###* Preparation:
###* We have a few ways in which we could approach our autogamy data. Unfortunately,
###* the data seems to have plenty of zero's and values exceed 1, which makes it
###* very compicated to try and fix into the beta distribution, as we did with 
###* PL.
###* 
###* To fit A into a beta, we would need to use a zero_inflated beta, but also alter
###* the data, so that no value was above 1.
###* 
###* This would mean to alter at least 5 of our 17 observations, which is a lot
###* 
###* Another way we could transform the data is to multiply the index and use 
###* the full numbers to try and fit a gaussian distribution. This would seem the most
###* appropriate way to me at first. Let us try that
###* 


# ao.final.table$mean_ao_index <- pmin(pmax(ao.final.table$mean_ao_index, 0.0000001), 0.9999999)

ao.final.table <- ao.final.table %>%
  mutate(
    log_visitors = log(mean.visitors.per.minute),
    sd_log_visitors = sd.visitors.per.minute / mean.visitors.per.minute
  )


###*
formula_ao <- bf(mean_ao_index ~ me(log_visitors, sd_log_visitors) + (1|species) + (1|elevation))

fit_ao <- brm(
  formula = formula_ao,
  data = ao.final.table,
  family = zero_inflated_beta(),
  chains = 6,
  cores = 10,
  iter = 4000,
  warmup = 1000,
  control = list(adapt_delta = 0.99, max_treedepth = 20),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  )
)


summary(fit_ao)
saveRDS(fit_ao, file = "../brms_models/brms_ao_1.rds")
fit_ao <- readRDS("../brms_models/brms_ao_1.rds")
###* The results from the summary tell us that with higher visitation, higher 
###* autogamy would be expected. However, this is not statistically significant!
###* (credible interval includes 0)
###* 
###* This makes very little sense, I would expect the opposite. So, am I setting 
###* the model correctly?
###* 
###* I will try to set a model in which I use weights, ie the amount of replicates 
###* used to measure the mean. For this, row 511 and 517 were added in "Bayesian setup.R".
###* The model should take into account how strongly each value from mean_ao_index
###* should be accounted for
formula_ao_weighted <- bf(mean_ao_index | weights(n_replicates) ~ me(log_visitors, sd_log_visitors) + (1|species) + (1|elevation))

fit_ao_weighted <- brm(
  formula = formula_ao_weighted,
  data = ao.final.table,
  family = zero_inflated_beta(),
  chains = 10,
  cores = 10,
  iter = 4000,
  warmup = 1000,
  control = list(adapt_delta = 0.999, max_treedepth = 20),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  seed = 1234
)

saveRDS(fit_ao_weighted, file = "../brms_models/brms_ao_weighted.rds")
fit_ao_weighted <- readRDS("../brms_models/brms_ao_weighted.rds")
summary(fit_ao_weighted)
###* These results tell us that with more visitors, there will be higher autogamy.
###* However, the model has very high Rhat and simultaniuosly very low values of
###* Bulk_ESS and Tail_ESS, both pointing for a bad model fit. It would seem that
###* the model structure might be too complex
###* 
###* I will try to exclude measurement error in the predictor, since it makes more
###* sense to me than getting rid of the weigths
formula_ao_weighted_2 <- bf(mean_ao_index | weights(n_replicates) ~ log_visitors + (1|species) + (1|elevation))

saveRDS(fit_ao_weighted_2, file = "../brms_models/brms_ao_weighted_2.rds")
fit_ao_weighted_2 <- readRDS("../brms_models/brms_ao_weighted_2.rds")
summary(fit_ao_weighted_2)
###* These results tell us, that with higher visitation, less autogamy would be 
###* expected. I would expect these results. Rhat is 1 and both ESS are reasonably
###* high. This migth be a good model. The result are statistically significant
###* 
###* However, we were using the log of visitors per minute. What is we used the
###* original value?
formula_ao_weighted_3 <- bf(mean_ao_index | weights(n_replicates) ~ mean.visitors.per.minute + (1|species) + (1|elevation))

saveRDS(fit_ao_weighted_3, file = "../brms_models/brms_ao_weighted_3.rds")
fit_ao_weighted_3 <- readRDS("../brms_models/brms_ao_weighted_3.rds")
summary(fit_ao_weighted_3)
###* These results tell us the opposite of the model with log transformed values
###* of visitation, while also being statistically significant. Also, Rhat is 1
###* and ESS are fine
###*
###*Which of the models should be believed? 
pp_check(fit_ao_weighted_2, type = "dens_overlay")
pp_check(fit_ao_weighted_3, type = "dens_overlay")
###* both have a reasonable pp
###* 
###* What about loo?
loo2 <- loo(fit_ao_weighted_2)
loo3 <- loo(fit_ao_weighted_3)
loo_compare(loo2, loo3)
###* The models are indistinguishable in therms of predictive perfomrmance
###* 
###* Both have many problematic data points (14 observations with paraketo_k > 0.7),
###* which is most of our observations. 
###*
###* Would a model with all info but original data istead of logged be better? 
formula_ao_weighted_4 <- bf(mean_ao_index | weights(n_replicates) ~ me(mean.visitors.per.minute, sd.visitors.per.minute) + (1|species) + (1|elevation))

saveRDS(fit_ao_weighted_4, file = "../brms_models/brms_ao_weighted_4.rds")
fit_ao_weighted_4 <- readRDS("../brms_models/brms_ao_weighted_4.rds")
summary(fit_ao_weighted_4)
###* Nope, same problem as the fit_ao_weighted
###* 
###* So, we are looking at weighted_2 and weighted_3. I do not know, which woul be 
###* better. I would say that the loged verion, as I would not expect autogamy 
###* to increase lineraly and normal data should be better for interpretation,
###* however this might be just me talking, the person who want results which make
###* sense





loo2 <- loo(fit_ao_weighted_2, moment_match = TRUE)
loo3 <- loo(fit_ao_weighted_3, moment_match = TRUE)
loo_compare(loo2, loo3)








View(ao.final.table)

formula_ao_weighted <- bf(mean_ao_index | weights(n_replicates) ~ me(log_visitors, sd_log_visitors) + (1|species) + (1|elevation))

ao_ordbetareg <- ordbetareg(
  formula = formula_ao_weighted,
  data = ao.final.table,
  chains = 10,
  cores = 10,
  iter = 4000,
  warmup = 1000,
  control = list(adapt_delta = 0.999, max_treedepth = 20),
  seed = 1234
)

saveRDS(ao_ordbetareg, file = "../brms_models/brms_ao_ordbetareg.rds")
ao_ordbetareg <- readRDS("../brms_models/brms_ao_ordbetareg.rds")
summary(ao_ordbetareg)


formula_ao <- bf(mean_ao_index ~ me(log_visitors, sd_log_visitors) + (1|ID1|species) + (1|ID2|elevation),
                 zi ~ log_visitors + (1|ID1|species) + (1|ID2|elevation))

fit_ao_2 <- brm(
  formula = formula_ao,
  data = ao.final.table,
  family = beta,
  chains = 7,
  cores = 10,
  iter = 12000,
  warmup = 6000,
  control = list(adapt_delta = 0.999, max_treedepth = 20),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("normal(0, 1)", class = "b", dpar = "zi"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  )
)

pairs(fit_ao_2)
rhat(fit_ao_2)
pp_check(fit_ao_2)


summary(fit_ao_2)

###* 
###* 
###* 
###* 
###* Now let's try the go_index
###* 
View(go.final.table)

go.final.table <- na.omit(go.final.table)
go.final.table$elevation <- as.factor(go.final.table$elevation)
go.final.table$species <- as.factor(go.final.table$species)

go.final.table <- go.final.table %>%
  mutate(
    log_visitors = log(mean.visitors.per.minute),
    sd_log_visitors = sd.visitors.per.minute / mean.visitors.per.minute
  )

formula_go <- bf(mean_go_index ~ me(log_visitors, sd_log_visitors) + (1|species) + (1|elevation),
                 zi ~ log_visitors + (1|species) + (1|elevation))

fit_go <- brm(
  formula = formula_go,
  data = go.final.table,
  family = zero_one_inflated_beta(),
  chains = 6,
  cores = 10,
  iter = 4000,
  warmup = 1000,
  control = list(adapt_delta = 0.99, max_treedepth = 20),
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  )
)

summary(fit_ao)
saveRDS(fit_ao, file = "../brms_models/brms_ao_1.rds")
