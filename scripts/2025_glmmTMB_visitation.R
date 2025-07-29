###* Loading the necessary packages
pacman::p_load(tidyverse, glmmTMB, DHARMa, emmeans)

select <- dplyr::select

set.seed(123)
source("scripts/Bayesian setup June short.R")


###* This script is to test the effect of the predictors "total.morphospecies", 
###* "total.functional.groups" and "visited.flowers.per.minute"
###* 
###* 

View(final.table)

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
