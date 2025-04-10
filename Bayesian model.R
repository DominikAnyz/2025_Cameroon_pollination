source("setup.R")
source("Bayesian setup.R")

set.seed(1234)

Sys.setenv(TZ = "UTC")

View(c.pl.final.table)

c.pl.final.table <- na.omit(c.pl.final.table)
ao.final.table <- na.omit(ao.final.table)
go.final.table <- na.omit(go.final.table)

c.pl.final.table$elevation <- as.factor(c.pl.final.table$elevation)
c.pl.final.table$species <- as.factor(c.pl.final.table$species)
ao.final.table$elevation <- as.factor(ao.final.table$elevation)
ao.final.table$species <- as.factor(ao.final.table$species)
go.final.table$elevation <- as.factor(go.final.table$elevation)
go.final.table$species <- as.factor(go.final.table$species)

summary(c.pl.final.table$mean_seedset)
hist(c.pl.final.table$mean_seedset)
hist(log(c.pl.final.table$mean_seedset))

###* 
###* 
###* ---CONTROL SEEDSET---
###* Thought: both the response and the predictor have to have a normal distribution,
###* specifically they are assumed to have a normal distribution when being used 
###* in the way which we are using them in
###* 
###* We have to look at the the distributions and alter them accordingly 
###* 
hist(c.pl.final.table$mean_seedset)
hist(log(c.pl.final.table$mean_seedset))
hist(c.pl.final.table$mean.visitors.per.minute)
hist(log(c.pl.final.table$mean.visitors.per.minute))

###* We have to alter all: mean seedset, standard error of the seedset, \
###* mean visitors per minute and the standard error of visitors per minute
###* 
###* We will makes logs of the seedset and visitors, we will use the delta method 
###* on the se of seedset and sd visitors
c.pl.final.table <- c.pl.final.table %>%
  mutate(log_mean_seedset = log(mean_seedset),
         se_log_mean_seedset = se_seedset / mean_seedset
         )

c.pl.final.table <- c.pl.final.table %>%
  mutate(
    log_visitors = log(mean.visitors.per.minute),
    sd_log_visitors = sd.visitors.per.minute / mean.visitors.per.minute
        )

###* Check
hist(c.pl.final.table$log_mean_seedset)
hist(c.pl.final.table$log_visitors)

###* Set up formula
formula <- bf(log_mean_seedset | se(se_log_mean_seedset) ~ me(log_visitors, sd_log_visitors) + (1|species) + (1|elevation))

#* Setting priors
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

summary(fit)
pairs(fit)
pp_check(fit)
conditional_effects(fit)
saveRDS(fit, file = "brms_seedset.rds")

###* The model summary tells us:
###* -visitation rate is a significant predictor of seedset (loglog scale)
###* -both species and elevation play a role in the seedset, but species more
###* -number of observations is low (17), which is why we are using SE and so
###* the sigma is 0
###*
###*
###*
###*
###*
###*
###*
###* --- POLLEN LIMITATION---
###* So we have some results for seedset, how about pollen limitation (PL)?
###* 
# Transform mean_PL_index and sd_PL_index to avoid exact 0s and 1s
c.pl.final.table$mean_PL_index <- pmin(pmax(c.pl.final.table$mean_PL_index, 0.0000001), 0.9999999)

hist(c.pl.final.table$mean_PL_index)

###* Look into documentation of brms and check what can work for the Beta distribution 
###* instead of standard error
###* 
###* It seems that we can try a few variations: The beta family does support 
###* measurement error in the predictor, but not in the response. This is 
###* because the variance is implicitly tied to the mean via its' formula. So we
###* can only model the mean and the preciscion (phi)
###* 
###* At any rate, in PL we will not be using the SE or SD of the mean
###* 
###* We can try the mean PL index, without the "measurement error". For better 
###* visualization, we will keep the data about the visitors in the log scale
formula <- bf(mean_PL_index ~ me(log_visitors, sd_log_visitors) + (1|species) + (1|elevation))

#* Setting priors
prior <- c(
  set_prior("normal(0, 2.5)", class = "b"),
  set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
  set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
)

fit2 <- brm(
  formula = formula,
  data = c.pl.final.table,
  chains = 6,
  family = Beta(),
  cores = 10,
  iter = 8000,
  warmup = 4000,
  prior = prior,
  control = list(max_treedepth = 20, adapt_delta = 0.999)
)

summary(fit2)
pairs(fit2)
saveRDS(fit2, file = "brms_PL_2.rds")
###* The model summary tells us:
###* There is a negative effect of log_visitors on mean_PL_index - as (log) 
###* visitation rate increases, pollen limitation decreases
###*  
###*
###*
###* Another way to try is would be to try a weighted beta regression, to assign 
###* more weight to groups with more pollen limitation replicates
formula <- bf(mean_PL_index | weights(n_replicates) ~ log_visitors + (1|species) + (1|elevation))

fit_weighted <- brm(
  formula = formula,
  data = c.pl.final.table,
  family = Beta(),
#  weights = "n_replicates", 
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 6,
  cores = 10,
  iter = 4000,
  warmup = 1000,
  control = list(max_treedepth = 20, adapt_delta = 0.99)
)

summary(fit_weighted)
saveRDS(fit_weighted, file = "brms_PL_3.rds")
###* This model tells me what I would expect: that with increased visitation, 
###* pollen limitation decreases. It is borderline significant though
###*
###*
###* I would like to try the moel without log transformed visitors, although
###* it seems to me, taht it make more sense to have the model use the log 
###* transformed version
formula <- bf(mean_PL_index | weights(n_replicates) ~ mean.visitors.per.minute + (1|species) + (1|elevation))

fit_weighted_2 <- brm(
  formula = formula,
  data = c.pl.final.table,
  family = Beta(),
  #  weights = "n_replicates", 
  prior = c(
    set_prior("normal(0, 2.5)", class = "b"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
    set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
  ),
  chains = 6,
  cores = 10,
  iter = 4000,
  warmup = 1000,
  control = list(max_treedepth = 20, adapt_delta = 0.99)
)

summary(fit_weighted_2)
saveRDS(fit_weighted_2, file = "brms_PL_4.rds")
###* Tjis model shows us a statistically insignificant trend in the opposite
###* direction than would be expected. However, in my opinion, not using the log
###* scale for the predictor when possible is a shame, since the relationshiop 
###* between pollination and seedset doesn't necessarilly need to be linear.
###* 
###* 
###* 

formula <- bf(mean_PL_index ~ me(mean.visitors.per.minute, sd.visitors.per.minute) +  (1|species) + (1|elevation)) 
#              phi ~  (1|species) + (1|elevation))

###* 
###* 
###* --- AUTOGAMY ---
View(ao.final.table)
hist(ao.final.table$mean_ao_index)

###* Preparation
ao.final.table$mean_ao_index <- pmin(pmax(ao.final.table$mean_ao_index, 0.0000001), 0.9999999)

ao.final.table <- na.omit(ao.final.table)

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
  family = Beta(),
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
saveRDS(fit_ao, file = "brms_ao_1.rds")
###* The results from the summary tell us that with higher visitation, higher 
###* autogamy would be expected. However, this is not statistically significant!
###* (credible interval includes 0)
###* 
###* This makes very little sense, I would expect the opposit. So, am I setting 
###* the model correctly?
###* 
ao.final.table <- left_join(final.table, ao.index2, by = "elevation.species")

ao.final.table <- na.omit(ao.final.table)

ao.final.table <- ao.final.table %>%
  mutate(
    log_visitors = log(mean.visitors.per.minute),
    sd_log_visitors = sd.visitors.per.minute / mean.visitors.per.minute
  )

ao.final.table$elevation <- as.factor(ao.final.table$elevation)
ao.final.table$species <- as.factor(ao.final.table$species)

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
saveRDS(fit_ao, file = "brms_ao_1.rds")















