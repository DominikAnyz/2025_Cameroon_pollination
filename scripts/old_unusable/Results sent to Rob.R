#*What do we have and what do we want?
#*
#*What data do we have?
#*
#*We currently have two datasets from which we have to extract information and
#*combine it.
#*
#*We have data of various indexes based on elevation, with plant species and 
#*plant number (plant ID) as random effects

pacman::p_load(tidyverse, glmmTMB, DHARMa)

#*
# Load data----
seed.data<- read.delim("indexes untangled 2.txt", na = c("na"))

# Give each observation a number
seed.data <- seed.data %>%  
  tibble::rowid_to_column("flowerID")

#* 1) Create mean plant number per plant species, per individual, per elevation
#* 2) Create mean plant number per plant species, per elevation
#* 3) Create column for the proper mean to use. If a plant individual has an
#* outcrossing tretment, use that. If no outcrossing treatment, use mean for 
#* elevation per species
#* 4) Create indexes
#* 5) Create column with "seedset", which is rounded (some Hypericum were not 
#* rounded)
seed.indices <- 
  seed.data %>% 
  group_by(elevation,plant.number, species) %>% 
  mutate(O.mean.plantnumber = mean(seedset[treatment == "outcrossing"], na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(elevation, species) %>% 
  mutate(O.mean.elevation = mean(seedset[treatment == "outcrossing"], na.rm = TRUE)) %>%
  mutate(O.mean = case_when(is.nan(O.mean.plantnumber) ~ O.mean.elevation,
                           .default = O.mean.plantnumber)) %>% 
  mutate(index = seedset/O.mean.elevation) %>%
  mutate(seedset = if_else(seedset != round(seedset), round(seedset), seedset))

#View(seed.indices)

#If the seedset is NA or Inf, then change to 0
seed.indices <- seed.indices %>% 
  filter(!(index %in% c("NA", "Inf")|is.na(seedset))) %>% 
  mutate(index = case_when(
    is.finite (index) ~index,
    .default = 0
  ))

#Elevation and species as factor
seed.indices$elevation<-as.factor(seed.indices$elevation)
seed.indices$species<-as.factor(seed.indices$species)

#*Control index and PL.index 
c.index <- seed.indices %>%
  filter(treatment == "control") %>%
  select(seedset, index, elevation, species, plant.number, O.mean.elevation) %>%
  mutate(seedset = round (seedset)) %>%
  mutate(PL.index = replace(1 - index, 1 - index < 0, 0)) %>%
  mutate(elevation = factor(elevation),
         species = factor(species)) %>%
  tibble::rowid_to_column("ID") %>%
  mutate(plant.id = as.factor(paste0(elevation,species,plant.number))) %>%
  filter(species != "Hypericum r" | elevation != 4000) %>%
  mutate(flower.id = as.factor(paste0(elevation, species, plant.number,"_", ID))) 

#View(c.index)

#*A/O index
ao.index <- seed.indices %>%
  filter(treatment == "autogamy") %>%
  select(index, elevation, species, plant.number, plant.number) %>%
  mutate(elevation = factor(elevation),
         species = factor(species)) %>%
  mutate(index = round(index * 1000)) %>%
  mutate(plant.number.el = as.factor(paste0(elevation,plant.number))) %>%
  mutate(plant.id = as.factor(paste0(elevation,species,plant.number))) %>%
  filter(species != "Hypericum r" | elevation != 4000) %>%
  mutate(species.sp = as.factor(paste0(elevation,species)))

#View(ao.index)

#*G/O index
go.index <- seed.indices %>%
  filter(treatment == "geitonogamy") %>%
  select(index, elevation, species, plant.number) %>%
  mutate(elevation = factor(elevation),
         species = factor(species)) %>%
  mutate(index = round(index * 1000)) %>%
  mutate(plant.id = as.factor(paste0(elevation,species,plant.number)))%>%
  filter(species != "Lactuca i") %>%
  filter(species != "Hypericum r" | elevation != 4000)

#View(go.index)

#*
#*We have results from these analyses fit into glmmTMB, however we will probably
#*have to change the model fit for the SEM, since we will carry out an SEM for
#*every index
#*
#*First we have to focus on what we will want to say/find out from our data
#*We want to know if elevation (E) has an effect on the visitation rate (VR) and 
#*visitor species richness (VSR).
#*
#*We also then want to know, if E has an effect on the following indexes:
#*Autogamy (A), Geitonogamy (G), Outrcossing (O), Control (C), and Pollen
#*limitation (PL).
#*
#*We also want to know, whether VR and VSR have an effect on the indexes
#*   
#*So we must comprise a dataset by index, which will have a value for for every 
#*plant individual, for every species of these various elevation. To this we will
#*add columns specific from visitation data in different elevations. The data 
#*which is necessary to be exctracted is as follows:
#*   -visitation rate (VR): the mean value of visitation frequencies
#*   -visitation rate standard deviation (VRSD): the mean of the sd from the 
#*   mean of visitarion rate
#*   -visitor species richness (VSR): the sum of unique species richness
#*   - primary pollinator richnes (PPR): the proprotion of visits by the most
#*   active pollinator
#*   - primary functional group abundance:
#*   - amount of time spent by pollinators on flowers
#*   
#*So how do we ectract this info from the visitation data?
#load .txt file
visitors<- read.delim("C:/Users/domin/Desktop/Cameroon pollination/All additional data/Visitors2.txt")
functional <- read.delim("C:/Users/domin/Desktop/Cameroon pollination/All additional data/functional.txt")
View(functional)
View (visitors)


functional_groups <- c("Hoverfly", "Bee", "Wasp", "Bird", "Beetle", "Butterfly", "Moth", "Other fly")

View(functional_groups)

#Create column with duplicate values of "minutes" and then count minutes in 
#recording and flowering minutes.
visitors <- visitors %>% 
  mutate(duplicate.minutes = if_else(min == lag(min), number.of.observed.flowers, NA_real_)) %>%
  group_by(plant.code) %>%
  mutate(minutes.in.recording = n() - sum(!is.na(duplicate.minutes))) %>%
  mutate(flowering.minutes = sum(number.of.observed.flowers, na.rm = TRUE) -
           sum(duplicate.minutes, na.rm = TRUE)) %>%
  mutate(insect.order = ifelse(insect.order == "Syrphidae", "Diptera", insect.order)) %>%
  ungroup()

View(visitors)

visitors <- merge(visitors, functional, by = "SD.s.ID", all.x = TRUE)

flowering <- visitors %>%
  group_by(plant.code, plant.species, elevation) %>%
  summarize(
    flowering.minutes = sum(number.of.observed.flowers, na.rm = TRUE) - sum(duplicate.minutes, na.rm = TRUE)
  )
#View(flowering)

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

flowering.visited <- left_join(flowering, visited, by = "plant.code") %>%
  mutate(number.of.visited.flowers = replace_na(number.of.visited.flowers, 0),
         total.visitor.count = replace_na(total.visitor.count, 0),
         morpho.visitor.count = replace_na(morpho.visitor.count, 0))

flowering.visited <- flowering.visited %>%
  mutate(visitors.per.minute = ifelse(flowering.minutes > 0, total.visitor.count / flowering.minutes, NA))

View(flowering.visited)

mode_function <- function(x) {
  if (length(x) == 0) {
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

View(result)

final.table <- left_join(flowering.visited, result, by = "plant.code")

final.table <- final.table %>%
  mutate(
    total.func = ifelse(is.na(total.func), 0, total.func),
    total.morpho = ifelse(is.na(total.morpho), 0, total.morpho)
  )

final.table <- final.table %>%
  mutate(
    most.common.func.morpho = ifelse(is.na(most.common.func.morpho), most.common.func, most.common.func.morpho),
    proportion.most.common.morpho = most.common.morpho.count / total.visitor.count
  ) %>%
  select(-most.common.func)

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

final.table <- final.table %>%
  mutate(
    most.common.morpho.func.count = ifelse(is.na(most.common.morpho.func.count), most.common.func.count, most.common.morpho.func.count)
  )

final.table <- final.table %>%
  mutate(
    proportion.most.common.morpho = most.common.morpho.count / total.visitor.count,
    proportion.most.common.func = most.common.morpho.func.count / total.visitor.count
  )

final.table <- final.table %>%
  select(-flowering.minutes, -number.of.visited.flowers, -most.common.morpho,
         -most.common.morpho.count, - most.common.func.count, 
         -most.common.morpho.func.count)

final.table$most.common.func.morpho <- gsub(",.*", "", final.table$most.common.func.morpho)

View(final.table)

#write.csv(final.table, "visitors.rob.csv")
###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* I would like to tryt the whole thing in the same way, expect grouped not by 
###* plant.code, but by plant_elevation. This would give us mean for every plant
###* species in every elevation

visitors2 <- visitors %>%
  mutate(elevation_species = paste0(elevation, plant.species))

visitors2 <- visitors2 %>% 
  mutate(duplicate.minutes = if_else(min == lag(min), number.of.observed.flowers, NA_real_)) %>%
  group_by(elevation_species) %>%
  mutate(minutes.in.recording = n() - sum(!is.na(duplicate.minutes))) %>%
  mutate(flowering.minutes = sum(number.of.observed.flowers, na.rm = TRUE) -
           sum(duplicate.minutes, na.rm = TRUE)) %>%
  mutate(insect.order = ifelse(insect.order == "Syrphidae", "Diptera", insect.order)) %>%
  ungroup()

View(visitors2)

visitors2 <- merge(visitors2, functional, by = "SD.s.ID", all.x = TRUE)

flowering2 <- visitors2 %>%
  group_by(elevation_species, plant.species, elevation) %>%
  summarize(
    flowering.minutes = sum(number.of.observed.flowers, na.rm = TRUE) - sum(duplicate.minutes, na.rm = TRUE)
  )
#View(flowering)

visited2 <- visitors2 %>%
  filter(functional.group %in% functional_groups) %>%
  group_by(elevation_species) %>%
  summarize(
    number.of.visited.flowers = sum(number.of.visited.flowers, na.rm = TRUE),
    total.visitor.count = n(),
    morpho.visitor.count = sum(morphospecies == 1, na.rm = TRUE),
    .groups = 'drop'
  )
#View(visited)

flowering.visited2 <- left_join(flowering2, visited2, by = "elevation_species") %>%
  mutate(number.of.visited.flowers = replace_na(number.of.visited.flowers, 0),
         total.visitor.count = replace_na(total.visitor.count, 0),
         morpho.visitor.count = replace_na(morpho.visitor.count, 0))

flowering.visited2 <- flowering.visited2 %>%
  mutate(visitors.per.minute = ifelse(flowering.minutes > 0, total.visitor.count / flowering.minutes, NA))

View(flowering.visited2)

mode_function <- function(x) {
  if (length(x) == 0) {
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


result2 <- visitors2 %>%
  filter(functional.group %in% functional_groups) %>%
  group_by(elevation_species) %>%
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

View(result2)

final.table2 <- left_join(flowering.visited2, result2, by = "elevation_species")

final.table2 <- final.table2 %>%
  mutate(
    total.func = ifelse(is.na(total.func), 0, total.func),
    total.morpho = ifelse(is.na(total.morpho), 0, total.morpho)
  )

final.table2 <- final.table2 %>%
  mutate(
    most.common.func.morpho = ifelse(is.na(most.common.func.morpho), most.common.func, most.common.func.morpho),
    proportion.most.common.morpho = most.common.morpho.count / total.visitor.count
  ) %>%
  select(-most.common.func)

final.table2 <- final.table2 %>%
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

final.table2 <- final.table2 %>%
  mutate(
    most.common.morpho.func.count = ifelse(is.na(most.common.morpho.func.count), most.common.func.count, most.common.morpho.func.count)
  )

final.table2 <- final.table2 %>%
  mutate(
    proportion.most.common.morpho = most.common.morpho.count / total.visitor.count,
    proportion.most.common.func = most.common.morpho.func.count / total.visitor.count
  )

final.table2 <- final.table2 %>%
  select(-flowering.minutes, -number.of.visited.flowers, -most.common.morpho,
         -most.common.morpho.count, - most.common.func.count, 
         -most.common.morpho.func.count)

final.table2$most.common.func.morpho <- gsub(",.*", "", final.table2$most.common.func.morpho)

# Modify final.table2 as described
final.table2 <- final.table2 %>%
  select(-elevation_species) %>%  # Remove the original column
  mutate(elevation_species = paste0(elevation, species))  # Create the new column

View(final.table2)
###* visitation rate mean table
###* 

seed.data<- read.delim("data/clean_seeds.txt", na = c("na"))

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

View(seed.indices)

seed.indices <- seed.indices %>%
  mutate(elevation_species = paste0(elevation, species))



###* Create a dataset specifically for the control and PL treatment

c.index <- seed.indices %>%
  filter(treatment == "control") %>%
  select(seedset, index, elevation, species, plant_number, elevation_species) %>%
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
  filter(species != "Hypericum r" | elevation != 4000) %>%
  mutate(flower.id = as.factor(paste0(elevation, species, plant_number,"_", ID)))

View(c.index)

c.index2 <- c.index %>%
  group_by(elevation_species) %>%
  summarise(
    mean_seedset = round(mean(seedset, na.rm = TRUE) / n(), 3),
    mean_PL_index = mean(PL.index, na.rm = TRUE) / n()
  )

View(c.index2)

###* AO mean table
###* 

ao.index <- seed.indices %>%
  filter(treatment == "autogamy") %>%
  select(index, elevation, species, plant_number, plant_number, elevation_species) %>%
  mutate(elevation = factor(elevation),
         species = factor(species)) %>%
  #mutate(index = round(index * 1000)) %>%
  mutate(plant.number.el = as.factor(paste0(elevation,plant_number))) %>%
  mutate(plant.id = as.factor(paste0(elevation,species,plant_number))) %>%
  filter(species != "Hypericum r" | elevation != 4000) %>%
  mutate(species.sp = as.factor(paste0(elevation,species)))

ao.index2 <- ao.index %>%
  group_by(elevation_species) %>%
  summarise(
    mean_ao_index = round(mean(index, na.rm = TRUE) / n(), 3)
  )

View(ao.index2)

###* 
###* 
###* GO mean table
###* 

go.index <- seed.indices %>%
  filter(treatment == "geitonogamy") %>%
  select(index, elevation, species, plant_number, elevation_species) %>%
  mutate(elevation = factor(elevation),
         species = factor(species)) %>%
  #mutate(index = round(index * 1000)) %>%
  mutate(plant.id = as.factor(paste0(elevation,species,plant_number)))%>%
  filter(species != "Lactuca i") %>%
  filter(species != "Hypericum r" | elevation != 4000)

go.index2 <- go.index %>%
  group_by(elevation_species) %>%
  summarise(
    mean_go_index = round(mean(index, na.rm = TRUE) / n(), 3)
  )

View(go.index2)


merged_table <- final.table2 %>%
  left_join(c.index2, by = "elevation_species") %>%
  left_join(ao.index2, by = "elevation_species")

merged_table2 <- final.table2 %>%
  left_join(go.index2, by = "elevation_species")

View(merged_table)

merged_table_clean <- merged_table[!is.na(merged_table$mean_seedset), ]

View(merged_table_clean)

# Scatterplots
pairs(merged_table_clean[, c("mean_seedset", "mean_PL_index", "mean_ao_index", "visitors.per.minute")])

# Correlation analysis for control seedset
cor.test(merged_table_clean$mean_seedset, merged_table_clean$visitors.per.minute, method = "pearson")
cor.test(merged_table_clean$mean_seedset, merged_table_clean$visitors.per.minute, method = "spearman")

# Correlation analysis for pl index
cor.test(merged_table_clean$mean_PL_index, merged_table_clean$visitors.per.minute, method = "pearson")
cor.test(merged_table_clean$mean_PL_index, merged_table_clean$visitors.per.minute, method = "spearman")

# Correlation analysis for ao seedset
cor.test(merged_table_clean$mean_ao_index, merged_table_clean$visitors.per.minute, method = "pearson")
cor.test(merged_table_clean$mean_ao_index, merged_table_clean$visitors.per.minute, method = "spearman")






