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

visitors_morpho <- visitors %>%
  filter(morphospecies == 1, functional.group %in% functional_groups)

fg_elev_summary <- visitors_morpho %>%
  group_by(elevation, functional.group) %>%
  summarise(visits = n(), .groups = "drop")

fg_elev_table <- fg_elev_summary %>%
  pivot_wider(names_from = functional.group, values_from = visits, values_fill = 0)

fg_elev_table <- fg_elev_table %>%
  rowwise() %>%
  mutate(richness = sum(c_across(all_of(functional_groups)) > 0)) %>%
  ungroup()

ggplot(fg_elev_summary, aes(x = functional.group, y = visits, fill = functional.group)) +
  geom_col() +
  facet_wrap(~ elevation, ncol = 2) +
  labs(x = "Functional group", y = "Number of visits (morphospecies only)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

View(fg_elev_summary)







fg_elev_summary <- visitors %>%
  filter(morphospecies == 1, functional.group %in% functional_groups) %>%
  group_by(elevation, functional.group) %>%
  summarise(visits = n(), .groups = "drop")

func_richness <- fg_elev_summary %>%
  group_by(elevation) %>%
  summarise(number_of_func_total = sum(visits > 0), .groups = "drop")

most_common_func <- fg_elev_summary %>%
  group_by(elevation) %>%
  slice_max(order_by = visits, n = 1, with_ties = FALSE) 

fg_summary_table <- left_join(func_richness, most_common_func, by = "elevation")

View(fg_summary_table)




















