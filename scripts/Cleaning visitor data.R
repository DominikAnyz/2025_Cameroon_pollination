# Visitors bud in the correct folder

###* visitors
library(tidyverse)

###* load .txt file
visitors<- read.delim("C:/Users/domin/Desktop/Cameroon pollination/All additional data/Visitors2.txt")

###* View
View (visitors)

###* Create column with duplicate values of "minutes" and then count minutes in 
###* recording and flowering minutes.
visitors <- visitors %>% 
  mutate(duplicate.minutes= if_else(min == lag(min), number.of.observed.flowers, NA_real_)) %>%
  group_by(plant.code) %>%
  mutate(minutes.in.recording = n() - sum(!is.na(duplicate.minutes))) %>%
  mutate(flowering.minutes = sum(number.of.observed.flowers, na.rm = TRUE) -
           sum(duplicate.minutes, na.rm = TRUE)) %>%
  mutate(insect.order = ifelse(insect.order == "Syrphidae", "Diptera", insect.order)) %>%
  ungroup()

###* Create pivot to see sucessful calculation of flowering minutes per plant 
###* individual per elevation
visitors %>%
  group_by(plant.code) %>%
  summarize(flowering.minutes = sum(number.of.observed.flowers, na.rm = TRUE) -
              sum(duplicate.minutes., na.rm = TRUE)) %>%
  ungroup()

###* Pivot for flowering minutes per plant species per elevation
print(visitors %>%
        group_by(plant.species, elevation) %>%
        summarize(flowering.minutes = sum(number.of.observed.flowers, na.rm = TRUE) -
                    sum(duplicate.minutes, na.rm = TRUE)),
      n=30) %>%
  ungroup()

###* Visitors in each elevation per plant species and those visitors classified
###* as "Morphospicies = 1" in each elevation
print(visitors %>%
        group_by(plant.species, elevation) %>%
        summarize(SD.s.ID = sum(nzchar(SD.s.ID)),
                  morphospecies = sum(morphospecies == 1, na.rm = TRUE)),
      #morphospecies = sum(morphospecies == 0, na.rm = TRUE),
      n=30)









