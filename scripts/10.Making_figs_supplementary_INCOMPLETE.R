pacman::p_load(tidyverse, glmmTMB, DHARMa, emmeans, patchwork, dplyr, ggplot2, 
               corrplot, stringr)

source("scripts/04. Setup for visitation indices.R")

vars <- c("visited.flowers.per.minute", "total.morpho", "total.func")
R <- cor(final.table[, vars], use = "pairwise.complete.obs", method = "pearson")

corrplot::corrplot(R, method = "color", type = "upper", addCoef.col = "black")

vars
R








###* 
###* 
###* 
###* 
###* FIG S1 A - D ARE HERE SO THAT THE READER CAN LOOK INTO THE 
###* INTRASPECIFIC DIFFERENCES
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


select <- dplyr::select

set.seed(123)

###* 
###* Load data----

seed.data<- read.delim("data/clean_seeds.txt", na = c("na"))

seed.data$elevation<-as.factor(seed.data$elevation)

seed.data <- seed.data %>%
  mutate(elevation = recode(elevation,
                            `2300` = 2300,
                            `2800` = 2800,
                            `3500` = 3400,
                            `4000` = 3800))

#View(seed.data)

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

###* Create a dataset specifically for the control treatment

c.index <- seed.indices %>%
  filter(treatment == "control") %>%
  select(seedset, index, elevation, species, plant_number) %>%
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
  #filter(species != "Hypericum r" | elevation != 3800) %>%
  mutate(flower.id = as.factor(paste0(elevation, species, plant_number,"_", ID)))

c.index.2 <- c.index %>%
  filter(species != "Hypericum r" | elevation != 3800)

ao.index <- seed.indices %>%
  filter(treatment == "autogamy") %>%
  select(index, elevation, species, plant_number, plant_number) %>%
  mutate(elevation = factor(elevation),
         species = factor(species)) %>%
  #mutate(index = round(index * 100)) %>%
  mutate(plant.number.el = as.factor(paste0(elevation,plant_number))) %>%
  mutate(plant.id = as.factor(paste0(elevation,species,plant_number))) %>%
  filter(species != "Hypericum r" | elevation != 3800) %>%
  mutate(species.sp = as.factor(paste0(elevation,species)))

go.index <- seed.indices %>%
  filter(treatment == "geitonogamy") %>%
  select(index, elevation, species, plant_number) %>%
  mutate(elevation = factor(elevation),
         species = factor(species)) %>%
  #mutate(index = round(index * 100)) %>%
  mutate(plant.id = as.factor(paste0(elevation,species,plant_number)))%>%
  filter(species != "Lactuca i") %>%
  filter(species != "Hypericum r" | elevation != 3800)

library(dplyr)

species_full <- c(
  "Clematis s" = "Clematis simensis",
  "Geranium a" = "Geranium arabicum",
  "Hypericum r" = "Hypericum revolutum",
  "Lactuca i"  = "Lactuca inermis",
  "Senecio p"  = "Senecio purpureus",
  "Crepis h"   = "Crepis hypochoeridea",
  "Senecio b"  = "Senecio burtonii"
)

relabel_species <- function(df) {
  df %>%
    mutate(species = recode(species, !!!species_full)) %>%
    mutate(species = factor(species))   # clean factor levels
}

c.index   <- relabel_species(c.index)
c.index.2 <- relabel_species(c.index.2)
ao.index  <- relabel_species(ao.index)
go.index  <- relabel_species(go.index)










elev_colors <- c("2300"="#9bcd9b","2800"="#d4af37","3400"="#d4c7bd","3800"="#fdfbf2")

# Ensure elevation is an ordered factor in each table we touch
fix_elev <- function(df) {
  df %>% mutate(elevation = factor(as.character(elevation),
                                   levels = c("2300","2800","3400","3800")))
}

# ---- Generic plotting function (linear scale, species-specific breaks)
plot_index_by_species_linear <- function(df_species, value_col, y_lab, title_stub) {
  # value_col is a string; we use aes_string for simplicity
  sp_name <- unique(df_species$species)
  
  y_vals  <- df_species[[value_col]]
  y_min   <- min(y_vals, na.rm = TRUE)
  y_max   <- max(y_vals, na.rm = TRUE)
  
  # Keep zero in view (useful if indices can be < 0)
  span    <- c(min(0, y_min), y_max)
  y_breaks <- pretty(span, n = 8)
  y_lower  <- min(y_breaks)
  y_upper  <- max(y_breaks)
  
  ggplot(df_species, aes_string(x = "elevation", y = value_col, fill = "elevation")) +
    geom_violin(trim = TRUE, scale = "width", width = 0.7) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
    stat_summary(fun = median, geom = "point", shape = 18, size = 4, color = "red") +
    stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "red", fatten = 2) +
    labs(title = paste0(sp_name, " — ", title_stub, " by Elevation"),
         x = "Elevation", y = y_lab) +
    scale_y_continuous(breaks = y_breaks, limits = c(y_lower, y_upper),
                       expand = expansion(mult = c(0, 0.03))) +
    geom_hline(yintercept = y_breaks, linetype = "dotted", color = "grey70") +
    scale_fill_manual(values = elev_colors, drop = FALSE) +
    scale_x_discrete(labels = function(x) scales::comma(as.numeric(as.character(x)))) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.text  = element_text(size = 15),
      axis.title = element_text(size = 20, face = "bold"),
      legend.position = "none",
      plot.title = element_text(size = 16, face = "bold")
    )
}

# ---- Helper to generate & save species plots from a table
make_species_plots <- function(df, value_col, y_lab, title_stub, file_prefix, out_dir = "figs") {
  df <- fix_elev(df)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  species_list <- sort(unique(df$species))
  plot_list <- list()
  
  for (sp in species_list) {
    df_sp <- df %>% filter(species == sp)
    p <- plot_index_by_species_linear(df_sp, value_col, y_lab, title_stub)
    plot_list[[sp]] <- p
    
    safe_name <- str_to_lower(str_replace_all(sp, "[^A-Za-z0-9]+", "_"))
    ggsave(file.path(out_dir, paste0(file_prefix, "_", safe_name, ".pdf")),
           p, width = 7.5, height = 6)
  }
  
  # Multi-page PDF (one species per page)
  pdf(file.path(out_dir, paste0(file_prefix, "_all_species_multipage.pdf")),
      width = 7.5, height = 6)
  for (sp in species_list) print(plot_list[[sp]])
  dev.off()
  
  # Faceted quick overview with free y for easy per-species visibility
  facet_plot <- ggplot(df, aes_string(x = "elevation", y = value_col, fill = "elevation")) +
    geom_violin(trim = TRUE, scale = "width", width = 0.7) +
    geom_jitter(width = 0.15, size = 1.2, alpha = 0.6, color = "black") +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3.5, color = "red") +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "red", fatten = 2) +
    scale_fill_manual(values = elev_colors, drop = FALSE) +
    labs(x = "Elevation (m a.s.l.)", y = y_lab) +
    facet_wrap(~ species, ncol = 3, scales = "free_y") +
    scale_x_discrete(labels = function(x) scales::comma(as.numeric(as.character(x)))) +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 16, face = "bold"),
      legend.position = "none",
      strip.text = element_text(size = 12, face = "bold.italic")
    )
  
  ggsave(file.path(out_dir, paste0(file_prefix, "_faceted_freeY.pdf")),
         facet_plot, width = 16, height = 12)
  
  invisible(list(plots = plot_list, facet = facet_plot))
}


c_out <- make_species_plots(
  df         = c.index,
  value_col  = "seedset",
  y_lab      = "Natural seed set",
  title_stub = "Seedset",
  file_prefix= "supp_seedset"
)

c_out$facet

pl_out <- make_species_plots(
  df         = c.index.2,
  value_col  = "PL.index",
  y_lab      = "Pollen limitation",
  title_stub = "Pollen limitation",
  file_prefix= "supp_PLindex"
)

pl_out$facet

ao_out <- make_species_plots(
  df         = ao.index,
  value_col  = "index",
  y_lab      = "Autogamy",
  title_stub = "Autogamy",
  file_prefix= "supp_AOindex"
)

ao_out$facet

go_out <- make_species_plots(
  df         = go.index,
  value_col  = "index",
  y_lab      = "Geitonogamy",
  title_stub = "Geitonogamy",
  file_prefix= "supp_GOindex"
)

go_out$facet

# Save the faceted overviews as separate PDFs (10 x 7)

ggsave(
  filename = "figs/supp_seedset_facet.pdf",
  plot     = c_out$facet,
  width    = 10,
  height   = 6
)

ggsave(
  filename = "figs/supp_PLindex_facet.pdf",
  plot     = pl_out$facet,
  width    = 10,
  height   = 6
)

ggsave(
  filename = "figs/supp_AOindex_facet.pdf",
  plot     = ao_out$facet,
  width    = 10,
  height   = 6
)

ggsave(
  filename = "figs/supp_GOindex_facet.pdf",
  plot     = go_out$facet,
  width    = 10,
  height   = 6
)
