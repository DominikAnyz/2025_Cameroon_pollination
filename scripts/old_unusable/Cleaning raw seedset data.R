#* This R script serves as the first part in analyzing Mt. Cameroon pollination
#* data. I would like all of my steps to be properly written down, so that
#* my data is reproducible

### Load the necessary libraries
library(here)
library(readxl)
library(tidyverse)
library(stringr)
library(writexl)

### Define the file path using here()
file_path <- here("data", "Mt. Cameroon seed count raw.xlsx")

### Load the "All" sheet from the Excel file
data <- read_excel(file_path, sheet = "All")

### View the first few rows of the data
head(data)

### Change names of elevations to elevations a.s.l.
data <- data %>%
  mutate(
    Elevation = case_when(
      Elevation == "Mann's Spring" ~ "2300",
      Elevation == "Hut 2" ~ "2800",
      Elevation == "Wevondi cave" ~ "3500",
      Elevation == "Hut 3" ~ "4000",
      TRUE ~ Elevation # Keep other values unchanged, if any
    )
  )

### Modify the "Treatment" column
data <- data %>%
  mutate(
    Treatment = case_when(
      Treatment == "II" ~ "geitonogamy",
      Treatment == "I"  ~ "autogamy",
      Treatment == "O"  ~ "control",
      Treatment == "+"  ~ "outcrossing",
      TRUE ~ Treatment # Keep other values unchanged, if any
    )
  )

### Create the "species" column
data <- data %>%
  mutate(
    species = paste(Genus, str_sub(Species, 1, 1))
  )

### Make seedsets numeric
data$`Seed count` <- as.numeric(data$`Seed count`)
data$`Hypericum seedset` <- as.numeric(data$`Hypericum seedset`)

### Update Seed count by the actual seedset for Hypericum
### Hypericum pods consist of 5 chambers, however in our experiment, at least one
### chamber was often infected. The value in "Hypericum seedset" is the total number
### of seeds divided by the number of uninfetced chambers
data <- data %>%
  mutate(
    `Seed count` = if_else(Genus == "Hypericum" & !is.na(`Hypericum seedset`), `Hypericum seedset`, `Seed count`)
  )

### Remove unwanted columns columns
data <- data %>%
  select(-Genus, -Species, -`chambers`, -`Hypericum seedset`)

### Rename and reposition the column
data <- data %>%
  relocate(species, .after = 2) # Move "species" to the third position (after the 2nd column)

### Rename columns 
data <- data %>%
  rename(elevation = Elevation,
         plot = Plot,
         plant_number = 'Plant number',
         treatment = Treatment,
         seedset = 'Seed count',
         comment = Comments)

### Make seedsets numeric
data$seedset <- as.numeric(data$seedset)

### Keep only species which have enough data for reproduction
data <- data %>%
  filter(species %in% c("Senecio p", "Crepis h", "Clematis s", 
                        "Geranium a", "Hypericum r", "Lactuca i", "Senecio b"))

### Reorganize the data the way I want them
data <- data %>%
  mutate(species = factor(species, levels = c("Senecio p", "Crepis h", "Clematis s", 
                                              "Geranium a", "Hypericum r", "Lactuca i", "Senecio b"))) %>%
  arrange(species, elevation, plant_number, seedset)

#View(data)

### Delete any comments, which would suggest unusable data
data <- data %>%
  filter(!str_detect(comment, regex("^absent$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^all chambers infected$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^all chambers infected, tiny eaten$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^all chambers infected, tiny eaten, with insect!$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^broke from plant$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^broken from main plant, al chambers eaten$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^broken from main plant, didnt have time to mature, 5 chambers$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^brown, infected with insect, check bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^but infected$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^but infected, would have been more$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^but open bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^but some black, infected$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^clearly not finished, bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^dead plant$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^dead side$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^dried$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^fell off$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^fell off plant$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^flower not present$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^II's probably not viable$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%  
  filter(!str_detect(comment, regex("^inf$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infeted$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%  
  filter(!str_detect(comment, regex("^infcted, black balls$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected - check bag / not usefull$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected a bit$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected check bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected with 2 diptera!$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%  
  filter(!str_detect(comment, regex("^infected with bug$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected with hymenoptera, black balls, check bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected with insect$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected with insect, black balls, cgheck bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected with insect, check bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected with insect, fell off plant$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected with insect, which check bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected with insect, which one, check bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected with insect, which, check bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%  
  filter(!str_detect(comment, regex("^infected with insects$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%  
  filter(!str_detect(comment, regex("^infected with multiple insects, check bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infecte, with small grains - is it eaten seeds - check bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected, 10 tiny$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected, 30 small green$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected, again fly in bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected, almost, check bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected, balls, check bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%  
  filter(!str_detect(comment, regex("^infected, black balls$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>% 
  filter(!str_detect(comment, regex("^infected, black balls, check bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected, black sell$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected, black shell$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%  
  filter(!str_detect(comment, regex("^infected, black shell, many insects, check bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected, but by what? check bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected, but looks like there were small black$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected, but rest look small and non mature$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%  
  filter(!str_detect(comment, regex("^infected, but what to do with them, heck bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected, catepillar indise, black balls$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected, check bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%  
  filter(!str_detect(comment, regex("^infected, cool kukla, check bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%  
  filter(!str_detect(comment, regex("^infected, fallen off of plant$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected, hymenoptera$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected, might have been good, check bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected, not sure!!! was on table$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected, one seed with hole, black balls, check bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected, small balls$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected, some almost$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected, sticky$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^infected, sticky, check bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^is infected$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%  
  filter(!str_detect(comment, regex("^li, with insect, check bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>% 
  filter(!str_detect(comment, regex("^looks infected$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%   
  filter(!str_detect(comment, regex("^looks infected, check$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%   
  filter(!str_detect(comment, regex("^looks infected, check bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%    
  filter(!str_detect(comment, regex("^more buds in bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%    
  filter(!str_detect(comment, regex("^more buds inside, keep on side$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%  
  filter(!str_detect(comment, regex("^multiple buds", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^nm!$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^no bud$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%  
  filter(!str_detect(comment, regex("^no flower$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^no flower in bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%    
  filter(!str_detect(comment, regex("^nothing was there$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%  
  filter(!str_detect(comment, regex("^recount$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%  
  filter(!str_detect(comment, regex("^small$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%  
  filter(!str_detect(comment, regex("^tied badly$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%    
  filter(!str_detect(comment, regex("^tiny$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%  
  filter(!str_detect(comment, regex("^tiny, unfinished, infected$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^tiny, unfinished, infected with live insect$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%
  filter(!str_detect(comment, regex("^too many flowers$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%  
  filter(!str_detect(comment, regex("^two inflorescenses in bag$", ignore_case = TRUE)) | is.na(comment) | comment == "") %>%  
  filter(!str_detect(comment, regex("^infecte, with small grains - is it eaten seeds? - check bag$", ignore_case = TRUE)) | is.na(comment) | comment == "")

data %>%
  filter(!(`seedset` == 0))

data <- data %>%
  filter(!(species == "Crepis h" & comment == "li" & elevation == "4000")) 

data <- data %>%
  filter(!(species == "Crepis h" & comment == "li, check bag" & elevation == "4000"))

# Round up all values in the "seedset" column
data <- data %>%
  mutate(seedset = ceiling(seedset))

#View(data)

write.table(data, file = "data/clean_seeds.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE, 
            na = "") # Replace NA with an empty string


