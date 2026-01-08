fg_by_elev <- visitors %>%
  filter(functional.group %in% functional_groups) %>%   # keep only your 8 focal groups
  distinct(elevation, functional.group) %>%
  arrange(as.numeric(as.character(elevation)), functional.group)

fg_by_elev


fg_counts <- visitors %>%
  filter(functional.group %in% functional_groups) %>%
  count(elevation, functional.group, name = "n_records") %>%
  arrange(as.numeric(as.character(elevation)), desc(n_records))

fg_counts


max_elev <- max(as.numeric(as.character(visitors$elevation)), na.rm = TRUE)

bees_highest <- visitors %>%
  filter(functional.group == "Bee",
         as.numeric(as.character(elevation)) == max_elev)

nrow(bees_highest)          # >0 means yes, there are bees
bees_highest %>% distinct(SD.s.ID)  # which IDs (morphospecies IDs)


bees_highest_morpho <- visitors %>%
  filter(functional.group == "Bee",
         morphospecies == 1,
         as.numeric(as.character(elevation)) == max_elev)

nrow(bees_highest_morpho)
bees_highest_morpho %>% distinct(SD.s.ID)


max_elev <- max(as.numeric(as.character(visitors$elevation)), na.rm = TRUE)

bee_plants_highest <- visitors %>%
  filter(functional.group == "Bee",
         as.numeric(as.character(elevation)) == max_elev) %>%
  distinct(plant.code, plant.species, elevation) %>%
  arrange(plant.species, plant.code)

bee_plants_highest

bee_visits_per_plant_highest <- visitors %>%
  filter(functional.group == "Bee",
         as.numeric(as.character(elevation)) == max_elev) %>%
  count(plant.code, plant.species, name = "bee_records") %>%
  arrange(desc(bee_records))

bee_visits_per_plant_highest





