pacman::p_load(tidyverse, glmmTMB, DHARMa, emmeans, paletteer, emmeans, ggpubr,
               multcompView, sf, ggplot2, ggspatial, elevatr, dplyr, raster,
               paletteer, marginaleffects, scales, multcomp, readr, stringr,
               tidybayes, brms, ordbetareg, patchwork)

select <- dplyr::select

set.seed(123)

seed.data<- read.delim("data/clean_seeds.txt", na = c("na"))

seed.data$elevation<-as.factor(seed.data$elevation)

# Recode elevation values globally before factor conversion
seed.data <- seed.data %>%
  mutate(elevation = recode(elevation,
                            `2300` = 2300,
                            `2800` = 2800,
                            `3500` = 3400,
                            `4000` = 3800))

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

#View(seed.indices)








###* FIGS 3A - control seedset
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
  filter(species != "Hypericum r" | elevation != 3800) %>%
  mutate(flower.id = as.factor(paste0(elevation, species, plant_number,"_", ID)))




















library(tidyverse)
library(emmeans)
library(ggplot2)

## Elevation gradient colours (same order as elsewhere)
elev_gradient_cols <- c(
  "#9bcd9b", # 2300
  "#d4af37", # 2800
  "#d4c7bd", # 3400
  "#fdfbf2"  # 3800
)

## Generic function to make the 2-sided gradient “violin”
## emm_df must have columns: elevation (factor/character), and a response column
make_gradient_violin <- function(emm_df,
                                 response_col = "mean",
                                 xlab = "",
                                 file_name = NULL,
                                 width = 5,
                                 height = 6) {
  df <- emm_df %>%
    transmute(
      elevation      = elevation,
      elevation_num  = as.numeric(as.character(elevation)),
      mean_val       = .data[[response_col]]
    ) %>%
    arrange(elevation_num) %>%
    filter(!is.na(mean_val))
  
  # Elevation sequence + linear interpolation between means
  min_e <- min(df$elevation_num)
  max_e <- max(df$elevation_num)
  
  elev_seq <- seq(min_e, max_e, length.out = 200)
  
  boundary <- approx(
    x    = df$elevation_num,
    y    = df$mean_val,
    xout = elev_seq,
    rule = 2
  )$y
  
  bdf <- tibble(
    elev_num = elev_seq,
    boundary = boundary
  ) %>%
    filter(!is.na(boundary))
  
  dy <- diff(bdf$elev_num)[1]
  
  # Right side rectangles
  shape_df <- tibble(
    xmin     = 0,
    xmax     = bdf$boundary,
    ymin     = bdf$elev_num - dy / 2,
    ymax     = bdf$elev_num + dy / 2,
    elev_num = bdf$elev_num
  )
  
  # Mirror to left side
  shape_df_both <- bind_rows(
    shape_df,
    shape_df %>%
      mutate(
        xmin = -xmax,
        xmax = -xmin
      )
  )
  
  boundary_df        <- tibble(x =  bdf$boundary, y = bdf$elev_num)
  boundary_df_mirror <- tibble(x = -bdf$boundary, y = bdf$elev_num)
  
  max_val <- max(abs(bdf$boundary), na.rm = TRUE)
  
  base_breaks <- pretty(c(0, max_val), n = 4)
  base_breaks <- unique(base_breaks)
  x_breaks <- c(-rev(base_breaks[base_breaks > 0]), base_breaks)
  
  p <- ggplot() +
    geom_rect(
      data = shape_df_both,
      aes(xmin = xmin, xmax = xmax,
          ymin = ymin, ymax = ymax,
          fill = elev_num),
      colour = NA
    ) +
    geom_path(
      data = boundary_df,
      aes(x = x, y = y),
      colour = "black",
      linewidth = 0.7
    ) +
    geom_path(
      data = boundary_df_mirror,
      aes(x = x, y = y),
      colour = "black",
      linewidth = 0.7
    ) +
    geom_point(
      data = df,
      aes(x = mean_val, y = elevation_num),
      colour = "black",
      size   = 3
    ) +
    geom_point(
      data = df,
      aes(x = -mean_val, y = elevation_num),
      colour = "black",
      size   = 3
    ) +
    scale_fill_gradientn(
      colours = elev_gradient_cols,
      name    = "Elevation (m)"
    ) +
    scale_x_continuous(
      limits = c(-max_val * 1.05, max_val * 1.05),
      breaks = x_breaks,
      labels = function(x) abs(x),
      name   = xlab
    ) +
    scale_y_continuous(
      breaks = df$elevation_num,
      labels = as.character(df$elevation),
      name   = "Elevation (m a.s.l.)"
    ) +
    coord_cartesian(expand = FALSE) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold", size = 16)
    )
  
  if (!is.null(file_name)) {
    ggsave(file_name, plot = p, width = width, height = height)
  }
  
  p
}





## Emmeans for control seed-set model
emm_c <- emmeans(
  c_brm_model_zinb_new,
  ~ elevation,
  type       = "response",
  re_formula = NA
)

emm_c_df <- as.data.frame(emm_c)

# Ensure we have a "mean" column
if (!"mean" %in% names(emm_c_df)) {
  if ("prob" %in% names(emm_c_df))      emm_c_df$mean <- emm_c_df$prob
  else if ("response" %in% names(emm_c_df)) emm_c_df$mean <- emm_c_df$response
  else if ("emmean" %in% names(emm_c_df))   emm_c_df$mean <- emm_c_df$emmean
}

seedset_violin <- make_gradient_violin(
  emm_c_df,
  response_col = "mean",
  xlab         = "Natural seed-set",
  file_name    = "figs/pres_seedset_violin_model.pdf"
)

seedset_violin

###*PL
pl_ord_model_new <- readRDS("../brms_models/pl_ord_model_new.rds")

# Get emmeans; ensure we are on the response scale.
# Some ordbetareg methods ignore type="response"; regrid() guarantees back-transform.
# 1) Emmeans on the LINK scale (this works)
emm_pl_link <- emmeans(pl_ord_model_new, ~ elevation)  # returns logit-scale summaries

pl_means <- emm_pl_df %>%
  transmute(
    elevation      = elevation,
    elevation_num  = as.numeric(as.character(elevation)),
    mean_PL        = mean
  ) %>%
  arrange(elevation_num) %>%
  filter(!is.na(mean_PL))          # just in case

pl_means

pl_violin <- make_gradient_violin(
  emm_pl_df,
  response_col = "mean",
  xlab         = "Pollen limitation index",
  file_name    = "figs/pres_PL_violin_model.pdf"
)

pl_violin

###*ao
## AUTOGAMY INDEX – build emmeans and gradient violin

ao_ord_model_new <- readRDS("../brms_models/ao_brm_model_zinb_new.rds")  # adjust name if needed

emm_ao_link <- emmeans(ao_ord_model_new, ~ elevation)  # on link scale
df_ao_link  <- as.data.frame(emm_ao_link)

# Identify link-scale columns robustly
mu_col_ao <- intersect(c("response","emmean","prob"), names(df_ao_link))[1]
lo_col_ao <- intersect(c("lower.HPD","lower.CL","asymp.LCL"), names(df_ao_link))[1]
hi_col_ao <- intersect(c("upper.HPD","upper.CL","asymp.UCL"), names(df_ao_link))[1]

# Back-transform via inverse-logit (if logit link; adjust if different)
emm_ao_df <- df_ao_link %>%
  transmute(
    elevation = elevation,
    mean = plogis(.data[[mu_col_ao]]),
    LCL  = plogis(.data[[lo_col_ao]]),
    UCL  = plogis(.data[[hi_col_ao]])
  )

ao_violin <- make_gradient_violin(
  emm_ao_df,
  response_col = "mean",
  xlab         = "Autogamy index",
  file_name    = "figs/pres_AO_violin_model.pdf"
)

ao_violin


###* GO
#run go up until getting emm_go_diff
# Make sure we have 'mean' column
if (!"mean" %in% names(emm_go_df)) {
  if ("prob" %in% names(emm_go_df)) emm_go_df$mean <- emm_go_df$prob
}

go_violin <- make_gradient_violin(
  emm_go_df,
  response_col = "mean",
  xlab         = "Geitonogamy index",
  file_name    = "figs/pres_GO_violin_model.pdf"
)

go_violin


###* Vis
#run vis up until getting emm_vis_diff
vis_violin <- make_gradient_violin(
  emm_vis_df,
  response_col = "mean",
  xlab         = "Visitation frequency (scaled)",
  file_name    = "figs/pres_VIS_violin_model.pdf"
)

vis_violin


###* Morph
#run morph up until getting emm_morpho_diff
# Ensure 'mean' column
if (!"mean" %in% names(emm_morpho_df)) {
  if ("prob" %in% names(emm_morpho_df)) emm_morpho_df$mean <- emm_morpho_df$prob
  else if ("response" %in% names(emm_morpho_df)) emm_morpho_df$mean <- emm_morpho_df$response
  else if ("emmean" %in% names(emm_morpho_df)) emm_morpho_df$mean <- emm_morpho_df$emmean
}

morpho_violin <- make_gradient_violin(
  emm_morpho_df,
  response_col = "mean",
  xlab         = "Morphospecies richness",
  file_name    = "figs/pres_MORPHO_violin_model.pdf"
)

morpho_violin


###* Func
#run func up until getting emm_func_df
# Ensure 'mean' column
if (!"mean" %in% names(emm_func_df)) {
  if ("rate" %in% names(emm_func_df)) emm_func_df$mean <- emm_func_df$rate
  else if ("response" %in% names(emm_func_df)) emm_func_df$mean <- emm_func_df$response
  else if ("emmean" %in% names(emm_func_df)) emm_func_df$mean <- emm_func_df$emmean
}

func_violin <- make_gradient_violin(
  emm_func_df,
  response_col = "mean",
  xlab         = "Functional-group richness",
  file_name    = "figs/pres_FUNC_violin_model.pdf"
)

func_violin





















