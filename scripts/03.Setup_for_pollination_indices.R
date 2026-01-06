###* Setup for pollination indices 
###* This script contains all of the code used to datasets to be used in the 
###* following script, 04.GLMM_pollination_indices
###* 
###* Loading the necessary packages
pacman::p_load(tidyverse)

set.seed(1234)

###* 
###* Load data---
seed.data<- read.delim("data/clean_seeds.txt", na = c("na"))

###* elevation as factor
seed.data$elevation<-as.factor(seed.data$elevation)

###* Change elevations to correct names. The values 3500 and 4000 were original
###* names, however they are changed in accordance with actual elevation
seed.data <- seed.data %>%
  mutate(elevation = recode(elevation,
                            `2300` = 2300,
                            `2800` = 2800,
                            `3500` = 3400,
                            `4000` = 3800))
###* Check dataset
#View(seed.data)

###* Create new dataset "seed.indices" with some configurations
###* The following code creates several new columns
seed.indices <- 
  seed.data %>% 
  group_by(elevation, plant_number, species) %>% 
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

###* Keep rows where seedset is not NA, change any values in column index to 0,
###* if they are NA or Inf. This happened when both a treatment "x" and the
###* outcrossing treatment produced 0 seeds. 
seed.indices <- seed.indices %>% 
  filter(!(index %in% c("NA", "Inf")|is.na(seedset))) %>% 
  mutate(index = case_when(
    is.finite (index) ~index,
    .default = 0
  ))

###* Elevation and species as factor
seed.indices$elevation<-as.factor(seed.indices$elevation)
seed.indices$species<-as.factor(seed.indices$species)

#Additional column "elevation.species"
seed.indices <- seed.indices %>%
  mutate(elevation.species = paste0(elevation, species))

###* Create a dataset specifically for the control treatment
c.index <- seed.indices %>%
  filter(treatment == "control") %>%
  select(seedset, index, elevation, species, plant_number, elevation.species) %>%
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
  mutate(flower.id = as.factor(paste0(elevation, species, plant_number,"_", ID)))

###* Create dataset for autogamy index 
ao.index <- seed.indices %>%
  filter(treatment == "autogamy") %>%
  select(index, elevation, species, plant_number, elevation.species) %>%
  mutate(elevation = factor(elevation),
         species = factor(species)) %>%
  # convert proportional index to 0–100 pseudo-count to model with count families
  mutate(index_trans = round(index * 100)) %>%
  mutate(plant.number.el = as.factor(paste0(elevation,plant_number))) %>%
  mutate(plant.id = as.factor(paste0(elevation,species,plant_number))) %>%
  #filter out Hypericum, as it produced no seeds in the highest elevation
  filter(species != "Hypericum r" | elevation != 3800) %>%
  mutate(species.sp = as.factor(paste0(elevation,species)))

###* Create dataset for geitonogamy index
go.index <- seed.indices %>%
  filter(treatment == "geitonogamy") %>%
  select(index, elevation, species, plant_number,elevation.species) %>%
  mutate(elevation = factor(elevation),
         species = factor(species)) %>%
  mutate(index_trans = round(index * 100)) %>%
  # convert proportional index to 0–100 pseudo-count to model with count families
  mutate(plant.id = as.factor(paste0(elevation,species,plant_number)))%>%
  # filter out Lactuca, as it does not have enough replications of the geitonogamy treatment
  filter(species != "Lactuca i") %>%
  # filter out HYpericum in highest elevation, as it produced not seeds there
  filter(species != "Hypericum r" | elevation != 3800)











