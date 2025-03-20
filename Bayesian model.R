source("setup.R")
source("Bayesian approach.R")

Sys.setenv(TZ = "UTC")

View(c.pl.final.table)
hist(log(c.pl.final.table$mean_seedset))

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

# Fit the model with optimized computation
fit <- brm(
  formula = formula,
  data = c.pl.final.table,
  chains = 10,   # More chains = better convergence
  cores = 10,    # Maximize CPU usage
  iter = 4000,   # Default is 2000; more for better precision
  warmup = 1000, # Standard warmup period
  control = list(adapt_delta = 0.95),  # Prevents divergent transitions
  prior = prior
)  



#* Specify the formula which I want to use
formula <- bf(mean_seedset | se(sd_seedset) ~ elevation + species + me(mean.visitors.per.minute, sd.visitors.per.minute))

#* Setting priors
prior <- c(
  set_prior("normal(0, 5)", class = "b"),  # Elevation effect
  set_prior("normal(40, 30)", class = "Intercept")  # Centered around observed mean (42)
)

fit <- brm(
  formula = formula,
  data = c.pl.final.table,
  chains = 5,  # Debugging mode
  cores = 10,
  iter = 4000,
  prior = prior,
  control = list(max_treedepth = 20, adapt_delta = 0.99)
)

summary(fit)
fixef(fit)
plot(fit)
saveRDS(fit, file = "brms_elevation_species.rds")

# Check for zero variance or constant variables
apply(c.pl.final.table, 2, function(x) var(x, na.rm = TRUE))

summary(c.pl.final.table$sd_seedset)

c.pl.final.table$sd_seedset <- ifelse(c.pl.final.table$sd_seedset == 0, 0.001, c.pl.final.table$sd_seedset)

#* Put species are rf or ff?

formula <- bf(mean_seedset | se(sd_seedset) ~ elevation + me(mean.visitors.per.minute, sd.visitors.per.minute) + (1|species))

#* Setting priors
prior <- c(
  set_prior("normal(0, 5)", class = "b"),  # Elevation effect
  set_prior("normal(40, 30)", class = "Intercept"),
  set_prior("lognormal(0, 1)", class = "sd", group = "species")
)

fit2 <- brm(
  formula = formula,
  data = c.pl.final.table,
  chains = 5,  # Debugging mode
  cores = 15,
  iter = 4000,
  prior = prior,
  control = list(max_treedepth = 20, adapt_delta = 0.99)
)

summary(fit2)
saveRDS(fit2, file = "brms_elevation_(1|species).rds")

###*
###*
###*

#* Put species are rf or ff?

formula <- bf(mean_seedset | se(sd_seedset) ~ elevation + me(mean.visitors.per.minute, sd.visitors.per.minute) + (1|species))

#* Setting priors
prior <- c(
  set_prior("normal(0, 5)", class = "b"),  # Elevation effect
  set_prior("normal(40, 30)", class = "Intercept"),
  set_prior("lognormal(0, 1)", class = "sd", group = "species")
)

fit3 <- brm(
  formula = formula,
  data = c.pl.final.table,
  chains = 5,
  cores = 25,
  iter = 6000,
  prior = prior,
  control = list(max_treedepth = 30, adapt_delta = 0.999)
)

summary(fit3)
pairs(fit3)
saveRDS(fit3, file = "brms_elevation_(1|species)_2.rds")

#*How is the histogram of seedset?
#*What about the other indexes?
#*Is the data in tune with the priors
#*
#*Am I now actually asking what is the effect of visitors per minute on the 
#*index, with elevation and species as a random effect?

###*

#* Put species are rf or ff?
hist(c.pl.final.table$mean_seedset)
formula <- bf(log1p(mean_seedset) | se(sd_seedset) ~ me(mean.visitors.per.minute, sd.visitors.per.minute) + (1|species) + (1|elevation))

#* Setting priors
prior <- c(
  set_prior("normal(0, 1)", class = "b"),
  set_prior("student_t(3, 0, 10)", class = "sd", group = "elevation"),
  set_prior("student_t(3, 0, 10)", class = "sd", group = "species")
)

fit4 <- brm(
  formula = formula,
  data = c.pl.final.table,
  chains = 5,
  cores = 25,
  iter = 6000,
  prior = prior,
  control = list(max_treedepth = 30, adapt_delta = 0.999)
)

summary(fit4)
pairs(fit4)
saveRDS(fit4, file = "brms_(1|elevation)_(1|species).rds")

###*
###*
###*

#* Put species are rf or ff?
hist(c.pl.final.table$mean_seedset)
formula <- bf(log1p(mean_seedset) | se(sd_seedset) ~ me(mean.visitors.per.minute, sd.visitors.per.minute) + (1|species) + (1|elevation))

hist(log(c.pl.final.table$mean.visitors.per.minute))

#* Setting priors
prior <- c(
  set_prior("normal(0, 1)", class = "b"),
  set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
  set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
)

fit5 <- brm(
  formula = formula,
  data = c.pl.final.table,
  chains = 5,
  cores = 25,
  iter = 6000,
  prior = prior,
  control = list(max_treedepth = 30, adapt_delta = 0.999)
)

summary(fit5)
pairs(fit5)
saveRDS(fit5, file = "brms_(1|elevation)_(1|species).rds")

pp_check(fit5)

###*
###* SO we have some rsuts for seedset, how about PL?
###*

#* Put species are rf or ff?
#* 
# Transform mean_PL_index and sd_PL_index to avoid exact 0s and 1s
c.pl.final.table$mean_PL_index <- pmin(pmax(c.pl.final.table$mean_PL_index, 0.0000001), 0.9999999)
c.pl.final.table$sd_PL_index   <- pmin(pmax(c.pl.final.table$sd_PL_index, 0.0000001), 0.9999999)

hist(logit(c.pl.final.table$mean_PL_index))

# Apply the logit transformation
c.pl.final.table$mean_PL_index_logit <- log(c.pl.final.table$mean_PL_index / (1 - c.pl.final.table$mean_PL_index))
View(c.pl.final.table)
hist(c.pl.final.table$mean_PL_index_tr)

# Create new columns by multiplying the existing ones by 10,000
c.pl.final.table$mean_PL_index_tr <- c.pl.final.table$mean_PL_index * 10000
c.pl.final.table$sd_PL_index_tr <- c.pl.final.table$sd_PL_index * 10000

###* add species richness to model as fixed factor
###* 
###* look into documentation of brms and check what can work for the Beta distribution 
###* instead of standard deviation
formula <- bf(mean_PL_index ~ me(mean.visitors.per.minute, sd.visitors.per.minute) +  (1|species) + (1|elevation), 
              phi ~  (1|species) + (1|elevation))


#* Setting priors
prior <- c(
  set_prior("normal(0, 2.5)", class = "b"),
  set_prior("student_t(3, 0, 5)", class = "sd", group = "elevation"),
  set_prior("student_t(3, 0, 5)", class = "sd", group = "species")
)

fit6 <- brm(
  formula = formula,
  data = c.pl.final.table,
  chains = 6,
  family = Beta(),
  cores = 25,
  iter = 6000,
  warmup = 3000,
  prior = prior,
  control = list(max_treedepth = 20, adapt_delta = 0.999)
)

summary(c.pl.final.table)

summary(fit6)
pairs(fit6)
saveRDS(fit6, file = "brms_(1|elevation)_(1|species).rds")

###*
###*
###* And now for ao index
View(ao.final.table)
summary(ao.final.table$mean_ao_index)
hist(log(ao.final.table$mean_ao_index))

###*
###*
###*And now for the go index
summary(go.final.table$mean_seedset)
hist(go.final.table$mean_seedset)