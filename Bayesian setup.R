###* Here we will try a Bayesian approach to our data
###* Everything prepared should be here
###* 
###* Comparing seed data with visitation data
###* 
###* The first part is loading the data as we have already done in 2025_glmmTMB.
###* Individual steps for this part are explained in that script, here all the 
###* comments are deleteed to free up as much space as possible.
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
  filter(species != "Hypericum r" | elevation != 4000) %>%
  mutate(flower.id = as.factor(paste0(elevation, species, plant_number,"_", ID)))

# ###* Code to send data to Rob
# View(c.index)
# ###* We only want important columns, so delete "plant.id", "index" and "flower.id"
# c.index.rob <- c.index %>% select(-plant.id, -index, -flower.id)
# ###* Rearange columns
# c.index.rob <- c.index.rob %>% select(ID, elevation, species, plant_number, seedset, PL.index)
# ###* Rename columns
# c.index.rob <- c.index.rob %>% rename(plant.number = plant_number, pl.index = PL.index)
# ###* View
# View(c.index.rob)
# ###* Write table as csv to send to rob
# write.csv(c.index.rob, file = "send to rob/c.pl.index.rob.xlsx", row.names = FALSE)

str(c.index)
summary(c.index)
any(is.na(c.index$seedset))
any(is.na(c.index$elevation))


# C.glmer3 <- glmmTMB(seedset ~ elevation  + (1|species) + (1|plant.id),
#                     ziformula=~elevation,
#                     data = c.index,
#                     #family = poisson)
#                     family = nbinom2)
# 
# glm_model3 <- glmmTMB(PL.index ~ elevation + (1|species) +(1|plant.id), 
#                       data = c.index,
#                       family = ordbeta())

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

# ###* Code to send data to Rob
# ###* NECESSARY TO CHANGE INDEX.TRANS IN CODE WHEN ROUNDING
# View(ao.index)
# ###* We only want important columns, so delete "plant.id", "index" and "flower.id"
# ao.index.rob <- ao.index %>% select(-plant.id, -species.sp, -plant.number.el)
# ###* Add column ID
# ao.index.rob <- ao.index.rob %>% mutate(ID = 1:n())
# ###* Rearange columns
# ao.index.rob <- ao.index.rob %>% select(ID, elevation, species, plant_number, index, index.trans)
# ###* Rename ceratin columns
# ao.index.rob <- ao.index.rob %>% rename(plant.number = plant_number, ao.index = index, ao.index.trans = index.trans)
# ###* View table
# View(ao.index.rob)
# ###* Write table as csv to send to rob
# write.csv(ao.index.rob, file = "send to rob/ao.index.rob.csv", row.names = FALSE)


# ao.nb.zi <- glmmTMB(index ~ elevation + (1| species) + (1|plant.id),
#                     ziformula=~elevation,
#                     data = ao.index,
#                     #family = poisson)
#                     family = nbinom2)

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

# go.po.zi <- glmmTMB(index ~ elevation + (1| species) + (1|plant.id),
#                     ziformula=~elevation,
#                     data = go.index,
#                     family = poisson)

# ###* Code to send data to Rob
# ###* NECESSARY TO CHANGE INDEX.TRANS IN CODE WHEN ROUNDING
# View(go.index)
# ###* We only want important columns, so delete "plant.id", "index" and "flower.id"
# go.index.rob <- go.index %>% select(-plant.id, -species.sp, -plant.number.el)
# ###* Add column ID
# go.index.rob <- go.index.rob %>% mutate(ID = 1:n())
# ###* Rearange columns
# go.index.rob <- go.index.rob %>% select(ID, elevation, species, plant_number, index, index.trans)
# ###* Rename ceratin columns
# go.index.rob <- go.index.rob %>% rename(plant.number = plant_number, ao.index = index, ao.index.trans = index.trans)
# ###* View table
# View(go.index.rob)
# ###* Write table as csv to send to rob
# write.csv(go.index.rob, file = "send to rob/go.index.rob.csv", row.names = FALSE)

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
###* indivudal and then by elevation combined with species
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

###* Created datset "visited", in which only visitors from the functional groups
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

###* Genearte a value of "visitors.per.minute" and "morpho.per.minute" for every
###* "plant.code"
flowering.visited <- flowering.visited %>%
  mutate(visitors.per.minute = ifelse(flowering.minutes > 0, total.visitor.count / flowering.minutes, NA)) %>%
  mutate(morpho.per.minute = ifelse(flowering.minutes > 0, morpho.visitor.count / flowering.minutes, NA))

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

###* Generate a new table, where instead of "plant.code", we will be grouping by
###* "elevation.species"
final.table <- final.table %>%
  mutate(elevation.species = paste0(elevation, plant.species))

###* Now I would like to summarize what I can so that it is based on elevation.species
###* In this way, I will get the standard deviances from the 7 recordings per species
###* per elevation, which I can then compare with the indexes and their standard deviances
final.table <- final.table %>%
  group_by(elevation.species, elevation, plant.species) %>%
  summarize(
    flowering.minutes = sum(flowering.minutes, na.rm = TRUE),
    number.of.visited.flowers = sum(number.of.visited.flowers, na.rm = TRUE),
    total.visitor.count = sum(total.visitor.count, na.rm = TRUE),
    mean.visitors.per.minute = mean(visitors.per.minute, na.rm = TRUE),
    sd.visitors.per.minute = sd(visitors.per.minute, na.rm = TRUE),
    .groups = 'drop'
  )

hist(final.table$mean.visitors.per.minute)

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
visitors2<- read.delim("Visitors/Visitors2.txt")
visitors2 <- visitors2 %>%
  mutate(elevation.species = paste0(elevation, plant.species))

###* Create column with duplicate values of "minutes" and then count minutes in 
###* recording and flowering minutes
visitors2 <- visitors2 %>% 
  mutate(duplicate.minutes = if_else(min == lag(min), number.of.observed.flowers, NA_real_)) %>%
  group_by(elevation.species) %>%
  mutate(minutes.in.recording = n() - sum(!is.na(duplicate.minutes))) %>%
  mutate(flowering.minutes = sum(number.of.observed.flowers, na.rm = TRUE) -
           sum(duplicate.minutes, na.rm = TRUE)) %>%
  mutate(insect.order = ifelse(insect.order == "Syrphidae", "Diptera", insect.order)) %>%
  ungroup()

#View(visitors2)

###* Merging the "visitors" and "functional" datasets into "visitors" based on
###* column "SD.s.ID"
visitors2 <- merge(visitors2, functional, by = "SD.s.ID", all.x = TRUE)

###* Generate a dataset, where we get the following information about every 
###* "elevation.species":
###* - "total.visitors": number of individual visitor "species" 
###* - "total.func": number of functional groups (from 8 groups)
###* - "most.common.func": most common functional group
###* - "most.common.morpho": most common species (provided that they were 
###* identified into morphospecies level)
###* - "most.common.morpho.count": the amount of times the most common species
###* was present
###* - "most.common.func.count": number of functional groups identified

result.es <- visitors2 %>%
  filter(functional.group %in% functional_groups) %>%
  group_by(elevation.species) %>%
  summarize(
    sp.richness = n_distinct(SD.s.ID),
    sp.richness.morpho = n_distinct(SD.s.ID[morphospecies == 1]),    
    total.func = n_distinct(functional.group),
    most.common.func = mode_function(functional.group),
    most.common.func.count = sum(functional.group == most.common.func),    
    most.common.morpho = mode_function(SD.s.ID[morphospecies == 1]),
    most.common.morpho.count = sum(SD.s.ID[morphospecies == 1] == most.common.morpho, na.rm = TRUE),
    .groups = 'drop'
  )

#View(result.es)

###* Now I want to join the two tables together
final.table <- left_join(final.table, result.es, by = "elevation.species")

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

#View(final.table)

###* Export csv for Rob
###* Delete unecessary columns
final.table <- final.table %>%
  select(-elevation.species, -flowering.minutes, -number.of.visited.flowers)

final.table <- final.table %>%
  mutate(elevation.species = paste0(elevation, species)) %>%
  select(elevation.species, everything())

# ###* write csv for Rob
# write.csv(final.table, file = "send to rob/visitation.rob.csv", row.names = FALSE)

###* In the next part, seedset data is prepared for merging with visitor data
###* visitation rate mean table

###* C/PL mean and sd table

c.index2 <- c.index %>%
  group_by(elevation.species) %>%
  summarise(
    mean_seedset = round(mean(seedset, na.rm = TRUE), 3),
    se_seedset = round(sd(seedset, na.rm = TRUE) / sqrt(sum(!is.na(seedset))), 3),
    mean_PL_index = round(mean(PL.index, na.rm = TRUE), 3),
    se_PL_index = round(sd(PL.index, na.rm = TRUE) / sqrt(sum(!is.na(PL.index))), 3)
  )

replicates <- seed.indices %>%
  filter(treatment == "control") %>%
  group_by(elevation.species) %>%
  summarise(n_replicates = n())

#View(c.index)

# c.index %>%
#   ggplot(aes(x = log(seedset))) + 
#   geom_histogram(bins = 20, color = "black", fill = "skyblue", alpha = 0.7) +
#   facet_wrap(~ elevation.species) +   # Create separate plots for each species
#   theme_minimal() +                   # Use a minimal theme
#   labs(title = "Histograms of 'index' by Species",
#        x = "Index",
#        y = "Frequency") +
#   theme(strip.text = element_text(size = 10, face = "bold"))  # Customize facet labels

#hist(ao.index2$mean_ao_index)

c.pl.final.table <- left_join(final.table, c.index2, by = "elevation.species")
c.pl.final.table <- left_join(c.pl.final.table, replicates, by = "elevation.species")

#View(c.pl.final.table)
###* AO mean and sd table

# View(ao.index2)

ao.index2 <- ao.index %>%
  group_by(elevation.species) %>%
  summarise(
    mean_ao_index = mean(index, na.rm = TRUE),
    se_ao_index = round(sd(index, na.rm = TRUE) / sqrt(sum(!is.na(index))))
  )

ao.final.table <- left_join(final.table, ao.index2, by = "elevation.species")
#View(ao.index2)

###* GO mean and sd table

go.index2 <- go.index %>%
  group_by(elevation.species) %>%
  summarise(
    mean_go_index = mean(index, na.rm = TRUE),
    se_go_index = round(sd(index, na.rm = TRUE) / sqrt(sum(!is.na(index))))
  )

go.final.table <- left_join(final.table, go.index2, by = "elevation.species")
#View(go.index2)

# write.csv(merged_table, "visitors.rob.csv")

