###* This script is for setting up the data for the GLMM for analyzing 
###* visitation indices based on elevation
###* 
###* The difference between this setup and the next is that for the setup of 
###* analysis for visitation indices based on elevation, we need to have 
###* everything agregated at the observation level. For the analyses of 
###* pollination indices based on visitation indices, we need to have 
###* everything agregated at the elevation and species levels. Although it sounds 
###* like a small fix, I belive it is best to run these scripts separately in 
###* order to avoid complications.
###* 
pacman::p_load(tidyverse)

select <- dplyr::select

set.seed(1234)

seed.data<- read.delim("data/clean_seeds.txt", na = c("na"))

seed.data$elevation<-as.factor(seed.data$elevation)

seed.data <- seed.data %>%
  mutate(elevation = recode(elevation,
                            `2300` = 2300,
                            `2800` = 2800,
                            `3500` = 3400,
                            `4000` = 3800))

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
  #filter(species != "Hypericum r" | elevation != 3800) %>%
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
  filter(species != "Hypericum r" | elevation != 3800) %>%
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
  filter(species != "Hypericum r" | elevation != 3800)
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
visitors<- read.delim("data/visitors.txt")
#View(functional)

###* Loading .txt "functional", which contains the functional groups of all of 
###* visitors from the table "Visitors2"
functional <- read.delim("data/functional.txt")
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

# Count the total number of unique morphospecies
total.morphospecies <- visitors %>%
  filter(functional.group %in% functional_groups, morphospecies == 1) %>%
  distinct(SD.s.ID) %>%
  nrow()

total.morphospecies


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
result <- visitors %>%
  filter(functional.group %in% functional_groups) %>%
  group_by(plant.code) %>%
  summarize(
    total.morpho = n_distinct(SD.s.ID[morphospecies == 1]),
    total.func = n_distinct(functional.group),
    most.common.func = mode_function(functional.group),
    most.common.func.morpho = mode_function(functional.group[morphospecies == 1]),
    most.common.morpho = mode_function(SD.s.ID[morphospecies == 1]),
    most.common.morpho.count = sum(SD.s.ID[morphospecies == 1] == most.common.morpho, na.rm = TRUE),
    most.common.func.count = sum(functional.group == most.common.func),
    most.common.morpho.func.count = sum(functional.group == most.common.func.morpho),
    .groups = 'drop'
  )

#View(result)

###* Join the tables "flowering.visited" and "results" by "plant.code"
final.table <- left_join(flowering.visited, result, by = "plant.code")

###* Replace all "NA" values in columns "total.func" and "total.morpho" with "0"
final.table <- final.table %>%
  mutate(
    total.func = ifelse(is.na(total.func), 0, total.func),
    total.morpho = ifelse(is.na(total.morpho), 0, total.morpho)
  )

#View(final.table)

###* If the column "most.common.func.morpho" is empty, replace it with the value
###* from column "most.common.func"
###* 
###* If the column "most.common.func.morpho.count" is empty, replace it with the 
###* value from column "most.common.func.count"
###* 
###* Create column "proportion.most.common.morpho" from the number of visits by
###* the most common species identified into morphospecies level and the total 
###* number of visits on a given plant
final.table <- final.table %>%
  mutate(
    most.common.func.morpho = ifelse(is.na(most.common.func.morpho), most.common.func, most.common.func.morpho),
    most.common.morpho.func.count = ifelse(is.na(most.common.morpho.func.count), most.common.func.count, most.common.morpho.func.count),
    proportion.most.common.visitor = most.common.morpho.count / total.visitor.count,
    proportion.most.common.morpho = most.common.morpho.count / morpho.visitor.count
  ) %>%
  select(-most.common.func)

#View(final.table)