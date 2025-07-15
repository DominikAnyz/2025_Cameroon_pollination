###* Here we will try a Bayesian approach to our data
###* Everything prepared should be here
###* 
###* Comparing seed data with visitation data
###* 
###* The first part is loading the data as we have already done in 2025_glmmTMB.
###* Individual steps for this part are explained in that script, here all the 
###* comments are deleted to free up as much space as possible.
###* 
###* The only column difference between the datasets created here and the ones
###* created in 2025_glmmTMB is that these have an additional column 
###* "elevation.species", which is a combination of values of elevation and 
###* species for a given entry. This column is used for merging the visitor and
###* seedset data
pacman::p_load(tidyverse, glmmTMB, DHARMa)

select <- dplyr::select

set.seed(1234)

seed.data<- read.delim("data/clean_seeds.txt", na = c("na"))

seed.indices <- 
  seed.data %>% 
  group_by(elevation,plant_number, species) %>% 
  mutate(O.mean.plantnumber = mean(seedset[treatment == "outcrossing"], na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(elevation, species) %>% 
  mutate(O.mean.elevation = mean(seedset[treatment == "outcrossing"], na.rm = TRUE)) %>%
  mutate(O.mean = case_when(is.nan(O.mean.plantnumber) ~ O.mean.elevation,
                            .default = O.mean.plantnumber)) %>% 
  mutate(index = seedset/O.mean) %>%
  mutate(seedset = if_else(seedset != round(seedset), round(seedset), seedset))

seed.indices <- seed.indices %>% 
  filter(!(index %in% c("NA", "Inf")|is.na(seedset))) %>% 
  mutate(index = case_when(
    is.finite (index) ~index,
    .default = 0
  ))

seed.indices$elevation<-as.factor(seed.indices$elevation)
seed.indices$species<-as.factor(seed.indices$species)

#Additional column "elevation.species"
seed.indices <- seed.indices %>%
  mutate(elevation.species = paste0(elevation, species))

c.index <- seed.indices %>%
  filter(treatment == "control") %>%
  select(seedset, index, elevation, species, plant_number, elevation.species) %>%
  mutate(seedset = round (seedset)) %>%
  mutate(PL.index = replace(1 - index, 1 - index < 0, 0)) %>%
  mutate(elevation = factor(elevation),
         species = factor(species)) %>%
  tibble::rowid_to_column("ID") %>%
  mutate(plant.id = as.factor(paste0(elevation,species,plant_number))) %>%
  #filter(species != "Hypericum r" | elevation != 4000) %>%
  mutate(flower.id = as.factor(paste0(elevation, species, plant_number,"_", ID)))

# str(c.index)
# summary(c.index)
# any(is.na(c.index$seedset))
# any(is.na(c.index$elevation))

ao.index <- seed.indices %>%
  filter(treatment == "autogamy") %>%
  select(index, elevation, species, plant_number, plant_number, elevation.species) %>%
  mutate(elevation = factor(elevation),
         species = factor(species)) %>%
  mutate(index = index) %>%
  mutate(index_trans = round(index * 100)) %>%
  mutate(plant.number.el = as.factor(paste0(elevation,plant_number))) %>%
  mutate(plant.id = as.factor(paste0(elevation,species,plant_number))) %>%
  filter(species != "Hypericum r" | elevation != 4000) %>%
  mutate(species.sp = as.factor(paste0(elevation,species)))

go.index <- seed.indices %>%
  filter(treatment == "geitonogamy") %>%
  select(index, elevation, species, plant_number, elevation.species) %>%
  mutate(elevation = factor(elevation),
         species = factor(species)) %>%
  mutate(index = index) %>%
  mutate(index_trans = round(index * 100)) %>%
  mutate(plant.id = as.factor(paste0(elevation,species,plant_number)))%>%
  filter(species != "Lactuca i") %>%
  filter(species != "Hypericum r" | elevation != 4000)
###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* Next I would like to try two things getting info from visitation data by
###* indivdual and then by elevation combined with species
###*
###* So how do we extract this info from the visitation data?

###* Loading visitor data, which has all the information of the video recordings 
###* and the visitors from the videos
visitors<- read.delim("Visitors/Visitors2.txt")
#View(functional)

###* Loading .txt "functional", which contains the functional groups of all of 
###* visitors from the table "Visitors2"
functional <- read.delim("Visitors/functional.txt")
#View (visitors)

###* Creating a dataframe "functional_groups", which will contain only the 8
###* functional groups, which are important for us
functional_groups <- c("Hoverfly", "Bee", "Wasp", "Bird", "Beetle", "Butterfly", "Moth", "Other fly")

###* Create column with duplicate values of "minutes" and then count minutes in 
###* recording and flowering minutes.
visitors <- visitors %>% 
  mutate(duplicate.minutes = if_else(min == lag(min), number.of.observed.flowers, NA_real_)) %>%
  group_by(plant.code) %>%
  mutate(minutes.in.recording = n() - sum(!is.na(duplicate.minutes))) %>%
  mutate(flowering.minutes = sum(number.of.observed.flowers, na.rm = TRUE) -
           sum(duplicate.minutes, na.rm = TRUE)) %>%
  mutate(insect.order = ifelse(insect.order == "Syrphidae", "Diptera", insect.order)) %>%
  ungroup()

#View(visitors)

###* Merging the "visitors" and "functional" datasets into "visitors" based on
###* column "SD.s.ID" (shich stands for Sylvain Delabye's ID, who is the entomologist
###* in charge of identification)
visitors <- merge(visitors, functional, by = "SD.s.ID", all.x = TRUE)

#View(visitors)

###* Counting the flowering minutes, ie the amount of flowers present for visitors
###* to visit in a given video
flowering <- visitors %>%
  group_by(plant.code, plant.species, elevation) %>%
  summarize(
    flowering.minutes = sum(number.of.observed.flowers, na.rm = TRUE) - sum(duplicate.minutes, na.rm = TRUE)
  )

#View(flowering)

###* Created dataset "visited", in which only visitors from the functional groups
###* which we are interested in are present
###* 
###* For each "plant.code", summarize the total number of visited flowers, the 
###* total visitor count and the visitor count only if we were able to identify 
###* the visitor into morphospecies level
visited <- visitors %>%
  filter(functional.group %in% functional_groups) %>%
  group_by(plant.code) %>%
  summarize(
    number.of.visited.flowers = sum(number.of.visited.flowers, na.rm = TRUE),
    total.visitor.count = n(),
    morpho.visitor.count = sum(morphospecies == 1, na.rm = TRUE),
    .groups = 'drop'
  )

#View(visited)

###* Join the tables "flowering" and "visited" by "plant.code"
###* 
###* Replace "na"'s in columns "number.of.visited.flowers", "total.visitor.count" 
###* and "morpho.visitor.count" with values "0"
flowering.visited <- left_join(flowering, visited, by = "plant.code") %>%
  mutate(number.of.visited.flowers = replace_na(number.of.visited.flowers, 0),
         total.visitor.count = replace_na(total.visitor.count, 0),
         morpho.visitor.count = replace_na(morpho.visitor.count, 0))


###* Genearte a value of "visitors.per.minute", "morpho.per.minute" and "scaled.visitors"
###* for every "plant.code"
flowering.visited <- flowering.visited %>%
  mutate(
    visited.flowers.per.minute = ifelse(flowering.minutes > 0, number.of.visited.flowers / flowering.minutes, NA),
    visitors.per.minute = ifelse(flowering.minutes > 0, total.visitor.count / flowering.minutes, NA),
    morpho.per.minute = ifelse(flowering.minutes > 0, morpho.visitor.count / flowering.minutes, NA)
  )

#View(flowering.visited)

vis_mean <- mean(flowering.visited$visitors.per.minute, na.rm = TRUE)
vis_sd   <- sd(flowering.visited$visitors.per.minute, na.rm = TRUE)

vis_flow_mean <- mean(flowering.visited$visited.flowers.per.minute, na.rm = TRUE)
vis_flow_sd   <- sd(flowering.visited$visited.flowers.per.minute, na.rm = TRUE)

flowering.visited <- flowering.visited %>%
  mutate(
    scaled.visitors = if (!is.na(vis_sd) && vis_sd > 0) {
      (visitors.per.minute - vis_mean) / vis_sd
    } else {
      visitors.per.minute  # fallback
    },
    scaled.visited.flowers = if (!is.na(vis_flow_sd) && vis_flow_sd > 0) {
      (visited.flowers.per.minute - vis_flow_mean) / vis_flow_sd
    } else {
      visited.flowers.per.minute  # fallback
    }
  )

#View(flowering.visited)

###* Create a mode function to be able to identify the most common visitor in 
###* our dataset
mode_function <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_character_)
  }
  tbl <- table(x)
  max_freq <- max(tbl)
  most_common <- as.character(names(tbl)[which(tbl == max_freq)])
  if (length(most_common) > 1) {
    return(paste(most_common, collapse = ", "))
  } else {
    return(most_common)
  }
}

###* Generate a dataset, where we get the following information about every 
###* "plant.code":
###* - "total.morpho": number of individual visitor "species" (which were 
###* identified into morphospecies level)
###* - "total.func": number of functional groups (provided that they were 
###* identified into morphospecies level)
###* - "most.common.func": most common functional group
###* - "most.common.func.morpho": most common functional group (provided that 
###* they were identified into morphospecies level)
###* - "most.common.morpho": most common species (provided that they were 
###* identified into morphospecies level)
###* - "most.common.morpho.count": the amount of times the most common species
###* was present
###* - "most.common.func.count": number of functional groups identified
###* - "most.common.morpho.func.count": number of functional groups identified 
###* (provided that they were identified into morphospecies level)
flowering.visited <- visitors %>%
  filter(functional.group %in% functional_groups) %>%
  group_by(plant.code) %>%
  summarize(
    total.morpho = n_distinct(SD.s.ID[morphospecies == 1]),
    total.func = n_distinct(functional.group),
    #most.common.func = mode_function(functional.group),
    #most.common.func.morpho = mode_function(functional.group[morphospecies == 1]),
    #most.common.morpho = mode_function(SD.s.ID[morphospecies == 1]),
    #most.common.morpho.count = sum(SD.s.ID[morphospecies == 1] == most.common.morpho, na.rm = TRUE),
    #most.common.func.count = sum(functional.group == most.common.func),
    #most.common.morpho.func.count = sum(functional.group == most.common.func.morpho),
    .groups = 'drop'
  ) %>%
###* Make sure that the functional groups, which were not identified into
###* morphospecies level are taken into account. SO, if total.func (total 
###* amount of functional groups) is higher than the toatal.morpho (total amount
###* of morphospecies), the number from total.func will be added to total.morpho
###* instead of the original number
  mutate(
    total.morpho = ifelse(total.func > total.morpho, total.func, total.morpho)
  ) %>%
  right_join(flowering.visited, by = "plant.code") %>%
  mutate(
    total.morpho = replace_na(total.morpho, 0),
    total.func = replace_na(total.func, 0)
  )

#View(flowering.visited)


vis_morpho_mean <- mean(flowering.visited$total.morpho, na.rm = TRUE)
vis_morpho_sd   <- sd(flowering.visited$total.morpho, na.rm = TRUE)
vis_func_mean <- mean(flowering.visited$total.func, na.rm = TRUE)
vis_func_sd   <- sd(flowering.visited$total.func, na.rm = TRUE)

flowering.visited <- flowering.visited %>%
  mutate(
    scaled.total.morpho = if (!is.na(vis_morpho_sd) && vis_morpho_sd > 0) {
      (total.morpho - vis_morpho_mean) / vis_morpho_sd
    } else {
      total.morpho  # fallback
    },
    scaled.total.func = if (!is.na(vis_func_sd) && vis_func_sd > 0) {
      (total.func - vis_func_mean) / vis_func_sd
    } else {
      total.func  # fallback
    }
  )


#View(flowering.visited)

###* Generate a new table, where instead of "plant.code", we will be grouping by
###* "elevation.species"
final.table <- flowering.visited %>%
  mutate(elevation.species = paste0(elevation, plant.species))

#View(final.table)

# ###* If the column "most.common.func.morpho" is empty, replace it with the value
# ###* from column "most.common.func"
# ###* 
# ###* If the column "most.common.func.morpho.count" is empty, replace it with the 
# ###* value from column "most.common.func.count"
# ###* 
# ###* Create column "proportion.most.common.morpho" from the number of visits by
# ###* the most common species identified into morphospecies level and the total 
# ###* number of visits on a given plant
# final.table <- final.table %>%
#   mutate(
#     most.common.func.morpho = ifelse(is.na(most.common.func.morpho), most.common.func, most.common.func.morpho),
#     most.common.morpho.func.count = ifelse(is.na(most.common.morpho.func.count), most.common.func.count, most.common.morpho.func.count),
#     proportion.most.common.visitor = most.common.morpho.count / total.visitor.count,
#     proportion.most.common.morpho = most.common.morpho.count / morpho.visitor.count
#   ) %>%
#   select(-most.common.func)

#View(final.table)

###* Now I would like to summarize what I can so that it is based on elevation.species
###* In this way, I will get the standard deviances from the 7 recordings per species
###* per elevation, which I can then compare with the indexes and their standard deviances
final.table <- final.table %>%
  group_by(elevation.species, elevation, plant.species) %>%
  summarize(
    flowering.minutes = sum(flowering.minutes, na.rm = TRUE),
    number.of.visited.flowers = sum(number.of.visited.flowers, na.rm = TRUE),
    total.visitor.count = sum(total.visitor.count, na.rm = TRUE),
    mean.visited.flowers.per.minute = mean(visited.flowers.per.minute, na.rm = TRUE),
    sd.visited.flowers.per.minute = sd(visited.flowers.per.minute, na.rm = TRUE),    
    mean.visitors.per.minute = mean(visitors.per.minute, na.rm = TRUE),
    sd.visitors.per.minute = sd(visitors.per.minute, na.rm = TRUE),
    mean.visitors.scaled = mean(scaled.visitors, na.rm = TRUE),
    sd.visitors.scaled = sd(scaled.visitors, na.rm = TRUE) / vis_sd,
    mean.visited.flowers.scaled = mean(scaled.visited.flowers, na.rm = TRUE),
    sd.visited.flowers.scaled = sd(scaled.visited.flowers, na.rm = TRUE) / vis_flow_sd,
    mean.morpho = mean(total.morpho, na.rm = TRUE),
    sd.morpho = sd(total.morpho, na.rm = TRUE),
    mean.morpho.scaled = mean(scaled.total.morpho, na.rm = TRUE),
    sd.morpho.scaled = sd(scaled.total.morpho, na.rm = TRUE) / vis_morpho_sd,
    mean.func = mean(total.func, na.rm = TRUE),
    sd.func = sd(total.func, na.rm = TRUE),
    mean.func.scaled = mean(scaled.total.func, na.rm = TRUE),
    sd.func.scaled = sd(scaled.total.func, na.rm = TRUE) / vis_func_sd,
    .groups = 'drop'
  )

#hist(final.table$mean.visitors.per.minute)

###*
###*
###*
###*
###* Now I need to get information, which should be grouped by elevation.species
###* in order to be correct and them merge the two tables. This includes the most 
###* common visitors and their Genuses
###* 
###* 
###* 
###* 
###* 
###* Generate a new table, where instead of "plant.code", we will be grouping by
###* "elevation.species"
###*
###*
###* 
# visitors2<- read.delim("Visitors/Visitors2.txt")
# visitors2 <- visitors2 %>%
#   mutate(elevation.species = paste0(elevation, plant.species))
# 
# ###* Create column with duplicate values of "minutes" and then count minutes in 
# ###* recording and flowering minutes
# visitors2 <- visitors2 %>% 
#   mutate(duplicate.minutes = if_else(min == lag(min), number.of.observed.flowers, NA_real_)) %>%
#   group_by(elevation.species) %>%
#   mutate(minutes.in.recording = n() - sum(!is.na(duplicate.minutes))) %>%
#   mutate(flowering.minutes = sum(number.of.observed.flowers, na.rm = TRUE) -
#            sum(duplicate.minutes, na.rm = TRUE)) %>%
#   mutate(insect.order = ifelse(insect.order == "Syrphidae", "Diptera", insect.order)) %>%
#   ungroup()
# 
# #View(visitors2)
# 
# ###* Merging the "visitors" and "functional" datasets into "visitors" based on
# ###* column "SD.s.ID"
# visitors2 <- merge(visitors2, functional, by = "SD.s.ID", all.x = TRUE)
# 
# ###* Generate a dataset, where we get the following information about every 
# ###* "elevation.species":
# ###* - "total.visitors": number of individual visitor "species" 
# ###* - "total.func": number of functional groups (from 8 groups)
# ###* - "most.common.func": most common functional group
# ###* - "most.common.morpho": most common species (provided that they were 
# ###* identified into morphospecies level)
# ###* - "most.common.morpho.count": the amount of times the most common species
# ###* was present
# ###* - "most.common.func.count": number of functional groups identified
# 
# result.es <- visitors2 %>%
#   filter(functional.group %in% functional_groups) %>%
#   group_by(elevation.species) %>%
#   summarize(
#     sp.richness = n_distinct(SD.s.ID),
#     sp.richness.morpho = n_distinct(SD.s.ID[morphospecies == 1]),    
#     total.func = n_distinct(functional.group),
#     most.common.func = mode_function(functional.group),
#     most.common.func.count = sum(functional.group == most.common.func),    
#     most.common.morpho = mode_function(SD.s.ID[morphospecies == 1]),
#     most.common.morpho.count = sum(SD.s.ID[morphospecies == 1] == most.common.morpho, na.rm = TRUE),
#     .groups = 'drop'
#   )
# 
# #View(result.es)
# 
# ###* Now I want to join the two tables together
# final.table <- left_join(final.table, result.es, by = "elevation.species")

#View(final.table)


###* Change the names of the plant species to match other analyses
final.table <- final.table %>%
  rename(species = plant.species) %>%
  mutate(species = case_when(
    species == "arabicum" ~ "Geranium a",
    species == "simensis" ~ "Clematis s",
    species == "revolutum" ~ "Hypericum r",
    species == "inermis" ~ "Lactuca i",
    species == "purpureus" ~ "Senecio p",
    species == "hypochoeridea" ~ "Crepis h",
    species == "burtonii" ~ "Senecio b",
    TRUE ~ species  # keep the original value if no match is found
  ))

# ###* Calculate the proportions of the most common morphospecies from the 
# ###* "most.common.morpho.count" and the "total.visitor.count"
# ###* 
# ###* Calculate the proportions of the most common functional group from the 
# ###* "most.common.morpho.func.count" and the "total.visitor.count"
# final.table <- final.table %>%
#   mutate(
#     proportion.most.common.morpho = most.common.morpho.count / total.visitor.count,
#     proportion.most.common.func = most.common.morpho.func.count / total.visitor.count
#   )
#
#final.table$most.common.func.morpho <- gsub(",.*", "", final.table$most.common.func.morpho)

# View(final.table)

###* Delete unecessary columns
final.table <- final.table %>%
  select(-elevation.species, -flowering.minutes, -number.of.visited.flowers)

final.table <- final.table %>%
  mutate(elevation.species = paste0(elevation, species)) %>%
  select(elevation.species, everything())


###* 
###* 
###* 
###* 
###* 
###* 
###* In the next part, seedset data is prepared for merging with visitor data
###* visitation rate mean table

###* C/PL mean and sd table

c.index2 <- c.index %>%
  group_by(elevation.species) %>%
  summarise(
    mean_seedset = round(mean(seedset, na.rm = TRUE), 3),
    #se_seedset = round(sd(seedset, na.rm = TRUE) / sqrt(sum(!is.na(seedset))), 3),
    sd_seedset = round(sd(seedset, na.rm = TRUE), 3),
    mean_PL_index = round(mean(PL.index, na.rm = TRUE), 3),
    #se_PL_index = round(sd(PL.index, na.rm = TRUE) / sqrt(sum(!is.na(PL.index))), 3),
    sd_PL_index = sd(PL.index, na.rm = TRUE)#,
    #PL_index_weight = ifelse(sd_PL_index == 0 | is.na(sd_PL_index), 1, 1 / sd_PL_index),
    #seedset_weight = ifelse(sd_seedset == 0 | is.na(sd_seedset), 1, 1 / sd_seedset)
  )

###* Generate number of replicates per eelvation.species
replicates <- seed.indices %>%
  filter(treatment == "control") %>%
  group_by(elevation.species) %>%
  summarise(n_replicates = n())

###* Join tables
c.pl.final.table <- left_join(final.table, c.index2, by = "elevation.species")
c.pl.final.table <- left_join(c.pl.final.table, replicates, by = "elevation.species")

###* Make sure everything is as factor, without NA's
c.pl.final.table$elevation <- as.factor(c.pl.final.table$elevation)
c.pl.final.table$species <- as.factor(c.pl.final.table$species)
c.pl.final.table <- na.omit(c.pl.final.table)

#View(c.pl.final.table)

c.pl.final.table <- c.pl.final.table %>%
  mutate(mean_seedset_round = round(mean_seedset))

c.pl.final.table.2 <- c.pl.final.table

c.pl.final.table <- na.omit(c.pl.final.table)

library(scales)  # for rescale()

c.pl.final.table.4 <- c.pl.final.table.2 %>%
  mutate(
    seedset_weight_11sd = 1 + 1 / (1 + sd_seedset),
    seedset_weight_12 = scales::rescale(-sd_seedset, to = c(1, 2)),
    PL_index_weight_11sd = 1 + 1 / (1 + sd_PL_index),
    PL_index_weight_12 = scales::rescale(-sd_PL_index, to = c(1, 2)),    
  )


#View(c.pl.final.table.4)
###* AO mean and sd table

#View(ao.index)

ao.index2 <- ao.index %>%
  group_by(elevation.species) %>%
  summarise(
    mean_ao_index = round(mean(index, na.rm = TRUE),3),
    #se_ao_index = round(sd(index, na.rm = TRUE) / sqrt(sum(!is.na(index)))),
    sd_ao_index = sd(index, na.rm = TRUE),
    mean_ao_index_trans = round(mean(index_trans, na.rm = TRUE),2),
    #se_ao_index_trans = round(sd(index_trans, na.rm = TRUE) / sqrt(sum(!is.na(index_trans)))),
    ao_index_weight = ifelse(sd_ao_index == 0 | is.na(sd_ao_index), 1, 1 / sd_ao_index)
  )

replicates_ao <- seed.indices %>%
  filter(treatment == "autogamy") %>%
  group_by(elevation.species) %>%
  summarise(n_replicates = n())

ao.final.table <- left_join(final.table, ao.index2, by = "elevation.species")
ao.final.table <- left_join(ao.final.table, replicates_ao, by = "elevation.species")

ao.final.table <- na.omit(ao.final.table)
ao.final.table$elevation <- as.factor(ao.final.table$elevation)
ao.final.table$species <- as.factor(ao.final.table$species)

ao.final.table <- ao.final.table %>%
  mutate(
    ao_index_weight_11sd = 1 + 1 / (1 + sd_ao_index),
    ao_index_weight_12 = scales::rescale(-sd_ao_index, to = c(1, 2)),    
  )

###* GO mean and sd table

go.index2 <- go.index %>%
  group_by(elevation.species) %>%
  summarise(
    mean_go_index = round(mean(index, na.rm = TRUE),3),
    #se_go_index = round(sd(index, na.rm = TRUE) / sqrt(sum(!is.na(index)))),
    sd_go_index = sd(index, na.rm = TRUE),
    mean_go_index_trans = round(mean(index_trans, na.rm = TRUE),3),
    #se_go_index_trans = round(sd(index_trans, na.rm = TRUE) / sqrt(sum(!is.na(index_trans)))),
    go_index_weight = ifelse(sd_go_index == 0 | is.na(sd_go_index), 1, 1 / sd_go_index)
  )

go.final.table <- left_join(final.table, go.index2, by = "elevation.species")

go.final.table <- na.omit(go.final.table)
go.final.table$elevation <- as.factor(go.final.table$elevation)
go.final.table$species <- as.factor(go.final.table$species)

go.final.table <- go.final.table %>%
  mutate(mean_go_index = ifelse(mean_go_index > 1, 1, mean_go_index))

go.final.table <- go.final.table %>%
  mutate(
    go_index_weight_11sd = 1 + 1 / (1 + sd_go_index),
    go_index_weight_12 = scales::rescale(-sd_go_index, to = c(1, 2)),    
  )


