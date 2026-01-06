###* This script is for setting up the data for the GLMM for analyzing 
###* visitation indices based on elevation
###* 
###* The difference between this setup and the next (06.) is that for the setup 
###* of analysis for visitation indices based on elevation, we need to have 
###* everything agregated at the observation level. For the analyses of 
###* pollination indices based on visitation indices, we need to have 
###* everything agregated at the elevation and species levels. Although it sounds 
###* like a small fix, I belive it is best to run these scripts separately in 
###* order to avoid complications.
###* 
pacman::p_load(tidyverse)

###* Loading visitor data, which has all the information of the video recordings 
###* and the visitors from the videos
visitors2<- read.delim("data/visitors.txt")
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
###*
###*
###*
###* FIXED
View(visitors2)

visitors2 <- visitors2 %>%
  arrange(plant.code, day, hour, min) %>%        # include hour if you still have it
  group_by(plant.code) %>%
  mutate(
    duplicate.minutes = if_else(hour == lag(hour) & min == lag(min),
                                number.of.observed.flowers, NA_real_),
    minutes.in.recording = n() - sum(!is.na(duplicate.minutes)),
    flowering.minutes = sum(number.of.observed.flowers, na.rm = TRUE) -
      sum(duplicate.minutes, na.rm = TRUE)
  ) %>%
  ungroup()

sum(!is.na(visitors2$duplicate.minutes))
sum(visitors2$duplicate.minutes, na.rm = TRUE)

#View(visitors)

###* Merging the "visitors" and "functional" datasets into "visitors" based on
###* column "SD.s.ID" (shich stands for Sylvain Delabye's ID, who is the entomologist
###* in charge of identification)
any(duplicated(functional$SD.s.ID))   # should be FALSE

functional %>%
  count(SD.s.ID, sort = TRUE) %>%
  filter(n > 1)

functional %>%
  filter(SD.s.ID %in% c("Adelidaesp.01","Anthomyiasp.02","BlackandyellowEchthromorpha",
                        "Dolichopodidaesp.01","Micromothsp.04","Tephritoideasp.3dots","leafhopper")) %>%
  arrange(SD.s.ID)

functional_unique <- functional %>% distinct(SD.s.ID, .keep_all = TRUE)

visitors2 <- visitors2 %>% left_join(functional_unique, by = "SD.s.ID")

View(visitors2)





View(functional)

#View(visitors2)

###* Counting the flowering minutes, ie the amount of flowers present for visitors
###* to visit in a given video
flowering <- visitors2 %>%
  group_by(plant.code, plant.species, elevation) %>%
  summarize(
    flowering.minutes = sum(number.of.observed.flowers, na.rm = TRUE) - sum(duplicate.minutes, na.rm = TRUE)
  )

#View(flowering)

###* Created dataset "visited", in which only visitors2 from the functional groups
###* which we are interested in are present
###* 
###* For each "plant.code", summarize the total number of visited flowers, the 
###* total visitor count and the visitor count only if we were able to identify 
###* the visitor into morphospecies level
visited <- visitors2 %>%
  filter(functional.group %in% functional_groups) %>%
  group_by(plant.code) %>%
  summarize(
    number.of.visited.flowers = sum(number.of.visited.flowers, na.rm = TRUE),
    total.visitor.count = n(),
    morpho.visitor.count = sum(morphospecies == 1, na.rm = TRUE),
    .groups = 'drop'
  )

visited %>% summarize(min_total = min(total.visitor.count), max_total = max(total.visitor.count))

functional %>% count(SD.s.ID) %>% filter(n > 1) %>% arrange(desc(n))


# Count the total number of unique morphospecies
total.morphospecies <- visitors2 %>%
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
###* - "total.func": number of functional groups 
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