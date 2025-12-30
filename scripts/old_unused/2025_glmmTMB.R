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


C.glmer3 <- glmmTMB(seedset ~ elevation  + (1|species) + (1|plant.id),
                    ziformula=~elevation,
                    data = c.index,
                    family = nbinom2)


C.glmer3.null <- glmmTMB(seedset ~ 1  + (1|species) + (1|plant.id),
                    ziformula=~elevation,
                    data = c.index,
                    family = nbinom2)

summary(C.glmer3.null)
AIC(C.glmer3, C.glmer3.null)
emm <- emmeans(C.glmer3, ~elevation)
pairs(emm)

saveRDS(C.glmer3, "glm_outputs/c_model.rds")
saveRDS(C.glmer3.null, "glm_outputs/c_null.rds")


###* Check the summary and the redults with package DHARMA
summary(C.glmer3)
simulationOutput3 <- simulateResiduals(fittedModel = C.glmer3, plot = F)
testDispersion(simulationOutput3)
testZeroInflation(simulationOutput3)
testResiduals(simulationOutput3)
plot(simulationOutput3)
deviance(C.glmer3)

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