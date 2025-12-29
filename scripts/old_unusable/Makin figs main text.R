###* SCRIPT FOR THE CURRENT FIGURES APPEARING IN THE MAIN TEXT OF MANUSCRIPT
###* 
###* FIG 1 MAP IS GENERATED HERE, REST WAS DONE IN INKSCAPE
###* 
###* 
###* 
###* 
###* 
###* 
###* Making graphs
###* library(ggplot2)
pacman::p_load(tidyverse, glmmTMB, DHARMa, emmeans, paletteer, emmeans, ggpubr,
               multcompView, sf, ggplot2, ggspatial, elevatr, dplyr, raster,
               paletteer, marginaleffects, scales, multcomp, readr, stringr,
               rnaturalearth, cowplot)

select <- dplyr::select

set.seed(123)

###* FIGURE 1
# 1. Define sites and convert to sf
sites <- data.frame(
  elevation = c("2300 m", "2800 m", "3400 m", "3800 m"),
  lat = c(4 + 8.67 / 60, 4 + 11.64 / 60, 4 + 12.15 / 60, 4 + 12.56 / 60),
  lon = c(9 + 7.15 / 60, 9 + 11.86 / 60, 9 + 11.27 / 60, 9 + 10.80 / 60)
) %>%
  mutate(label = as.character(elevation))  

# Add Buea location
buea <- data.frame(
  name = "Buea",
  lon = 9.2306592,
  lat = 4.1568167
)

buea_sf <- st_as_sf(buea, coords = c("lon", "lat"), crs = 4326)

sites_sf <- st_as_sf(sites, coords = c("lon", "lat"), crs = 4326)

# 2. Create expanded bounding box polygon
bbox <- st_bbox(sites_sf) + c(-0.04, -0.04, 0.04, 0.05)  # Zoomed-in
bbox_poly <- st_as_sfc(st_bbox(bbox, crs = 4326))
bbox_sf <- st_sf(geometry = bbox_poly)

cameroon <- ne_countries(country = "Cameroon", returnclass = "sf")

# 3. Download elevation data using the polygon (not bbox directly)
elev_raster <- get_elev_raster(locations = bbox_sf, z = 10, clip = "bbox", prj = "+proj=longlat +ellps=WGS84")

# 4. Convert raster to dataframe
elev_df <- as.data.frame(rasterToPoints(elev_raster))
colnames(elev_df) <- c("Longitude", "Latitude", "elevation")

contour_lines <- rasterToContour(elev_raster, levels = seq(0, 4000, by = 200))
contour_sf <- st_as_sf(contour_lines)

# # 5. Color palette (same as your plots)
# elev_colors <- paletteer_c("grDevices::Green-Yellow", n = 8)

# 5. Custom color palette for elevation
elev_colors <- c("#123f3e", "#9bcd9b", "#d4af37", "#d4c7bd", "#fdfbf2")

elev_breaks <- c(0, 2300, 2800, 3400, 4000)

# Rescale values to [0, 1] range (as required by scale_fill_gradientn)
elev_values <- scales::rescale(elev_breaks, to = c(0, 1))

# Apply in the scale
scale_fill_gradientn(
  colors = elev_colors,
  values = elev_values,
  name = "Elevation (m)"
)

map <- ggplot() +
  # DEM raster
  geom_raster(data = elev_df, aes(x = Longitude, y = Latitude, fill = elevation)) +
  # Contours
  geom_sf(data = contour_sf, color = "grey30", size = 0.2, alpha = 0.6) +
  # MCNP boundary
  #geom_sf(data = mcnp_boundary, fill = NA, color = "black", size = 0.8) +
  # Sampling sites
  geom_sf(data = sites_sf, size = 4, shape = 21, fill = "white", color = "black") +
  # Add Buea point and label
  geom_sf(data = buea_sf, shape = 21, size = 4, fill = "red", color = "white") +
  geom_text(data = buea, aes(x = lon, y = lat, label = "Buea"), inherit.aes = FALSE, 
            hjust = -0.1, vjust = 1.5, size = 5,fontface = "bold") +
  
  geom_text(
    data = sites,
    aes(x = lon + 0.012, y = lat, label = elevation),
    size = 5, fontface = "bold"
  ) +
  # Aesthetics
  scale_fill_gradientn(
    colors = c("#123f3e", "#9bcd9b", "#d4af37", "#d4c7bd", "#fdfbf2"),
    values = scales::rescale(c(0, 2300, 2800, 3400, 4000)),
    name = "Elevation (m)"
  ) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering()) +
  coord_sf(
    xlim = c(9.10, 9.25),
    ylim = c(4.13, 4.23),
    expand = FALSE
  ) +
  scale_x_continuous(
    breaks = c(9.11, 9.24),
    labels = c("9.11°E", "9.24°E")
  ) +
  scale_y_continuous(
    breaks = c(4.14, 4.22),
    labels = c("4.14°N", "4.22°N")
  ) +
  labs(title = "Study Sites in Mount Cameroon National Park") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )

map

ggsave("figs/fig1_final.pdf", plot = map, width = 8, height = 6)
###*
###*
###*
###*
###*
###*
###*
###*GENERATING MAP TO SHOW WHICH PART OF CAMERRON THE EXPERIMENT IS TAKING PLACE IN
cameroon <- ne_countries(country = "Cameroon", returnclass = "sf")

bbox <- st_bbox(sites_sf) + c(-0.04, -0.04, 0.04, 0.05)
bbox_poly <- st_as_sfc(st_bbox(bbox, crs = 4326))
bbox_sf <- st_sf(geometry = bbox_poly)

cameroon_map <- ggplot() +
  geom_sf(data = cameroon, fill = "white", color = "black", size = 0.3) +
  geom_sf(data = bbox_sf, fill = "red", color = "red", alpha = 1, size = 1.2)+  # outline your study area
  coord_sf(xlim = c(8.4, 17), ylim = c(1.5, 13), expand = FALSE) +
  theme_void() +
  #labs(title = "Location of Study Area in Cameroon") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

cameroon_map

ggsave("figs/cameroon_outline_map.pdf", cameroon_map, width = 6, height = 8)
###
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
###* FIG 2 IS CREATED IN INKSCAPE
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
###* 
###* 
###* 
###* FIG 3 (A-F):
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

# Prepare summary counts for annotation
summary_data_control <- c.index %>%
  group_by(elevation) %>%
  summarise(count = n())

# Define consistent elevation colors
elev_colors <- c(
  "2300" = "#9bcd9b",
  "2800" = "#d4af37",
  "3400" = "#d4c7bd",
  "3800" = "#fdfbf2"
)

c.index$elevation <- factor(c.index$elevation, levels = c("2300", "2800", "3400", "3800"))
names(elev_colors) <- levels(c.index$elevation)

control_combined_plot <- ggplot(c.index, aes(x = elevation, y = seedset, fill = elevation)) +
  geom_violin(trim = TRUE, scale = "width", width = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
  labs(x = "Elevation", y = "Natural seed set") +
  scale_y_continuous(trans = pseudo_log_trans(base = 10), 
                     breaks = c(0, 5, 10, 25, 50, 100, 200, 400, 600)) +
  geom_hline(yintercept = c(0, 5, 10, 25, 50, 100, 200, 400, 600),
             linetype = "dotted", color = "grey70") +
  scale_fill_manual(values = elev_colors) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 15, face = "bold"),
    legend.position = "none"
  )


# Load Bayesian model
c_brm_model_zinb_new <- readRDS("../brms_models/c_brm_model_zinb_new.rds")

# posterior draws for pairwise elevation contrasts
contr_draws <- emmeans(c_brm_model_zinb_new, ~ elevation) %>%
  contrast(method = "pairwise") %>%
  gather_emmeans_draws()

# Summarise contrasts
post_probs <- contr_draws %>%
  group_by(contrast) %>%
  summarise(
    Estimate   = mean(.value),
    l95        = quantile(.value, 0.025),
    u95        = quantile(.value, 0.975),
    P_greater0 = mean(.value > 0),
    P_less0    = mean(.value < 0),
    .groups = "drop"
  ) %>%
  mutate(contrast_clean = gsub("elevation", "", contrast))

post_probs

manual_letters <- tibble(
  elevation = c("2300", "2800", "3400", "3800"),
  .group    = c("b", "c", "b", "a")
)

# Manually assign letters (lowest = a, mid = b, highest = c)
cld_c <- tibble(
  elevation = factor(c("2300", "2800", "3400", "3800"),
                     levels = c("2300", "2800", "3400", "3800")),
  .group = c("b", "c", "b", "a")  # your chosen letters
)

# Get max y-value for positioning letters
y_max_c <- max(c.index$seedset, na.rm = TRUE)

# Plot with letters
control_combined_plot_with_letters <- control_combined_plot +
  #coord_cartesian(ylim = c(0, 700)) +
  geom_text(data = cld_c,
            aes(x = elevation, y = 1100, label = .group),
            #aes(x = elevation, y = y_max_c * 1.6, label = .group), # 20% above max
            size = 5,
            fontface = "bold") #+
  # stat_summary(fun = median, geom = "point",
  #              shape = 18, size = 4, color = "red") +
  # stat_summary(fun = median, geom = "crossbar",
  #              width = 0.5, color = "red", fatten = 2)

control_combined_plot_with_letters

ggsave("figs/fig_C_bayes.pdf", 
       plot = control_combined_plot_with_letters, width = 8, height = 6)

summary(c_brm_model_zinb_new)


library(emmeans)
library(dplyr)

# make sure intervals are included and on the response scale
emm <- emmeans(c_brm_model_zinb_new, ~ elevation,
               type = "response", re_formula = NA)

emm_df <- as.data.frame(emm)

# Robustly standardize column names: mean, LCL, UCL
if (!"mean" %in% names(emm_df)) {
  if ("response" %in% names(emm_df)) emm_df$mean <- emm_df$response
  else if ("emmean" %in% names(emm_df)) emm_df$mean <- emm_df$emmean
}

if (!"LCL" %in% names(emm_df)) {
  if ("asymp.LCL" %in% names(emm_df)) emm_df$LCL <- emm_df$asymp.LCL
  else if ("lower.CL" %in% names(emm_df)) emm_df$LCL <- emm_df$lower.CL
  else if ("lower.HPD" %in% names(emm_df)) emm_df$LCL <- emm_df$lower.HPD
}

if (!"UCL" %in% names(emm_df)) {
  if ("asymp.UCL" %in% names(emm_df)) emm_df$UCL <- emm_df$asymp.UCL
  else if ("upper.CL" %in% names(emm_df)) emm_df$UCL <- emm_df$upper.CL
  else if ("upper.HPD" %in% names(emm_df)) emm_df$UCL <- emm_df$upper.HPD
}

# sanity check
names(emm_df); head(emm_df)

control_combined_plot_with_letters <- control_combined_plot_with_letters +
  geom_pointrange(data = emm_df,
                  aes(x = elevation, y = prob, ymin = LCL, ymax = UCL),
                  inherit.aes = FALSE, color = "red", size = 0.6)

control_combined_plot_with_letters <- control_combined_plot_with_letters +
  scale_x_discrete(labels = c(
    "2300" = "2,300",
    "2800" = "2,800",
    "3400" = "3,400",
    "3800" = "3,800"
  ))

control_combined_plot_with_letters
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
###* FIG 3B - POLLEN LIMITATION
# Colors + factor order (match your other figs)
elev_colors <- c("2300"="#9bcd9b","2800"="#d4af37","3400"="#d4c7bd","3800"="#fdfbf2")
c.index$elevation <- factor(c.index$elevation, levels = c("2300","2800","3400","3800"))
names(elev_colors) <- levels(c.index$elevation)

# Base violin + points
pl_plot <- ggplot(c.index, aes(x = elevation, y = PL.index, fill = elevation)) +
  geom_violin(trim = TRUE, scale = "width", width = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
  labs(x = "Elevation", y = "Pollen limitation") +
  scale_fill_manual(values = elev_colors) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1.1)) +
  geom_hline(yintercept = seq(0, 1, 0.25), linetype = "dotted", color = "grey70") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 15, face = "bold"),
    legend.position = "none"
  )

# --- MODEL-BASED marginal means on RESPONSE scale (+ 95% CrI) ----------------
# Load your Bayesian ordbetareg model (already fitted)
pl_ord_model_new <- readRDS("../brms_models/pl_ord_model_new.rds")

# Get emmeans; ensure we are on the response scale.
# Some ordbetareg methods ignore type="response"; regrid() guarantees back-transform.
# 1) Emmeans on the LINK scale (this works)
emm_pl_link <- emmeans(pl_ord_model_new, ~ elevation)  # returns logit-scale summaries

# 2) Manually back-transform via inverse-logit (plogis)
df_link <- as.data.frame(emm_pl_link)

# Identify the column names that exist in your object
mu_col <- intersect(c("response","emmean","prob"), names(df_link))[1]         # "emmean" in your printout
lo_col <- intersect(c("lower.HPD","lower.CL","asymp.LCL"), names(df_link))[1] # "lower.HPD"
hi_col <- intersect(c("upper.HPD","upper.CL","asymp.UCL"), names(df_link))[1] # "upper.HPD"

# Build response-scale summary in [0,1]
emm_pl_df <- df_link %>%
  transmute(
    elevation = elevation,
    mean = plogis(.data[[mu_col]]),
    LCL  = plogis(.data[[lo_col]]),
    UCL  = plogis(.data[[hi_col]])
  )

# 3) Overlay on your existing violin+points plot (pl_plot)
# (Assumes you've already created pl_plot and cld_pl as in your code)
y_max_pl <- max(c.index$PL.index, na.rm = TRUE)

cld_pl <- tibble(
  elevation = factor(c("2300", "2800", "3400", "3800"),
                     levels = levels(c.index$elevation)),
  .group    = c("a", "a", "a", "b")
)

pl_plot_with_letters <- pl_plot +
  geom_text(data = cld_pl,
            aes(x = elevation, y = y_max_pl + 0.1, label = .group),
            size = 5, fontface = "bold") +
  # model-based marginal mean (dot) + 95% CrI (vertical line)
  geom_pointrange(data = emm_pl_df,
                  aes(x = elevation, y = mean, ymin = LCL, ymax = UCL),
                  inherit.aes = FALSE, color = "red", size = 0.6)

pl_plot_with_letters <- pl_plot_with_letters +
  scale_x_discrete(labels = c(
    "2300" = "2,300",
    "2800" = "2,800",
    "3400" = "3,400",
    "3800" = "3,800"
  ))

pl_plot_with_letters

ggsave("figs/fig_PL_bayes.pdf", plot = pl_plot_with_letters, width = 8, height = 6)
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
###* AUTOGAMY
###* FIG 3C - AUTOGAMY
pacman::p_load(tidyverse, emmeans, multcomp, scales, ggplot2)

## 1) (Re)build ao.index in a clean way ---------------------------------------
##    - index_raw  = original G/O (0–4)
##    - index      = index_raw * 100  (model scale)

ao.index <- seed.indices %>%
  filter(treatment == "autogamy") %>%
  select(index, elevation, species, plant_number) %>%
  rename(index_raw = index) %>%                # keep original G/O
  mutate(
    index      = round(index_raw * 100),       # model scale (0–400)
    elevation  = factor(elevation),
    species    = factor(species),
    plant.id   = as.factor(paste0(elevation, species, plant_number))
  ) %>%
  filter(species != "Hypericum r" | elevation != 3800)
View(ao.index)
## If you already have ao.index with index on 0–400 and no index_raw,
## you can instead do:
## ao.index <- ao.index %>% mutate(index_raw = index / 100)

## 2) Elevation factor order + colours ----------------------------------------

elev_colors <- c(
  "2300" = "#9bcd9b",
  "2800" = "#d4af37",
  "3400" = "#d4c7bd",
  "3800" = "#fdfbf2"
)

ao.index <- ao.index %>%
  mutate(elevation = factor(elevation, levels = c("2300","2800","3400","3800")))

names(elev_colors) <- levels(ao.index$elevation)

## 3) Base violin + jitter on ORIGINAL scale (0–4) ----------------------------

ao_plot_base <- ggplot(ao.index,
                       aes(x = elevation, y = index_raw, fill = elevation)) +
  geom_violin(trim = TRUE, scale = "width", width = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
  labs(x = "Elevation", y = "Autogamy") +
  scale_fill_manual(values = elev_colors) +
  scale_y_continuous(breaks = seq(0, 6, 1), 
                     limits = c(0, 6.5)) +
  geom_hline(yintercept = seq(0, 6, 1),
             linetype = "dotted", color = "grey70") +
  theme_minimal() +
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    axis.text          = element_text(size = 10),
    axis.title         = element_text(size = 15, face = "bold"),
    legend.position    = "none"
  )

## 4) MODEL-BASED marginal means (ZINB, response scale, /100) -----------------

# Load the fitted model if needed
ao_brm_model_zinb_new <- readRDS("../brms_models/ao_brm_model_zinb_new.rds")

emm_ao <- emmeans(
  ao_brm_model_zinb_new,
  ~ elevation,
  type       = "response",   # expected index * 100
  re_formula = NA            # population-level effect
)

emm_ao_df <- as.data.frame(emm_ao)

## Standardise names: mean, LCL, UCL
if (!"mean" %in% names(emm_ao_df)) {
  if ("response" %in% names(emm_ao_df))      emm_ao_df$mean <- emm_ao_df$response
  else if ("emmean" %in% names(emm_ao_df))   emm_ao_df$mean <- emm_ao_df$emmean
}

if (!"LCL" %in% names(emm_ao_df)) {
  if ("asymp.LCL" %in% names(emm_ao_df))     emm_ao_df$LCL <- emm_ao_df$asymp.LCL
  else if ("lower.CL" %in% names(emm_ao_df)) emm_ao_df$LCL <- emm_ao_df$lower.CL
  else if ("lower.HPD" %in% names(emm_ao_df))emm_ao_df$LCL <- emm_ao_df$lower.HPD
}

if (!"UCL" %in% names(emm_ao_df)) {
  if ("asymp.UCL" %in% names(emm_ao_df))     emm_ao_df$UCL <- emm_ao_df$asymp.UCL
  else if ("upper.CL" %in% names(emm_ao_df)) emm_ao_df$UCL <- emm_ao_df$upper.CL
  else if ("upper.HPD" %in% names(emm_ao_df))emm_ao_df$UCL <- emm_ao_df$upper.HPD
}

head(emm_ao_df)

## Back-transform from index*100 to original 0–4 scale
emm_ao_df <- emm_ao_df %>%
  mutate(
    prob = prob / 100,
    LCL  = LCL  / 100,
    UCL  = UCL  / 100
  )

## 5) Compact letter display for elevations -----------------------------------

# cld_ao_manual <- tibble(
#   elevation = factor(c("2300", "2800", "3400", "3800"),
#                      levels = levels(ao.index$elevation)),
#   .group    = c("a", "a", "a", "a")
# )

y_max_ao <- max(ao.index$index_raw, na.rm = TRUE)
y_max_ao
## 6) Final plot: violins + points + CLD letters + red CrI --------------------

ao_plot_with_letters <- ao_plot_base +
  # geom_text(data = cld_ao_manual,
  #           aes(x = elevation, y = 6.5, label = .group),
  #           size = 5, fontface = "bold") +
  geom_pointrange(
    data        = emm_ao_df,
    aes(x = elevation, y = prob, ymin = LCL, ymax = UCL),
    inherit.aes = FALSE,
    color       = "red",
    size        = 0.6
  )

ao_plot_with_letters <- ao_plot_with_letters +
  scale_x_discrete(labels = c(
    "2300" = "2,300",
    "2800" = "2,800",
    "3400" = "3,400",
    "3800" = "3,800"
  ))

ao_plot_with_letters
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
###* FIG 3C - GEITONOGAMY
pacman::p_load(tidyverse, emmeans, multcomp, scales, ggplot2)

## 1) (Re)build go.index in a clean way ---------------------------------------
##    - index_raw  = original G/O (0–4)
##    - index      = index_raw * 100  (model scale)

go.index <- seed.indices %>%
  filter(treatment == "geitonogamy") %>%
  select(index, elevation, species, plant_number) %>%
  rename(index_raw = index) %>%                # keep original G/O
  mutate(
    index      = round(index_raw * 100),       # model scale (0–400)
    elevation  = factor(elevation),
    species    = factor(species),
    plant.id   = as.factor(paste0(elevation, species, plant_number))
  ) %>%
  filter(species != "Lactuca i") %>%
  filter(species != "Hypericum r" | elevation != 3800)

## If you already have go.index with index on 0–400 and no index_raw,
## you can instead do:
## go.index <- go.index %>% mutate(index_raw = index / 100)

## 2) Elevation factor order + colours ----------------------------------------

elev_colors <- c(
  "2300" = "#9bcd9b",
  "2800" = "#d4af37",
  "3400" = "#d4c7bd",
  "3800" = "#fdfbf2"
)

go.index <- go.index %>%
  mutate(elevation = factor(elevation, levels = c("2300","2800","3400","3800")))

names(elev_colors) <- levels(go.index$elevation)

## 3) Base violin + jitter on ORIGINAL scale (0–4) ----------------------------

go_plot_base <- ggplot(go.index,
                       aes(x = elevation, y = index_raw, fill = elevation)) +
  geom_violin(trim = TRUE, scale = "width", width = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
  labs(x = "Elevation", y = "Geitonogamy") +
  scale_fill_manual(values = elev_colors) +
  scale_y_continuous(breaks = seq(0, 4.2, 1), 
                     limits = c(0, 4.5)) +
  geom_hline(yintercept = seq(0, 4.2, 1),
             linetype = "dotted", color = "grey70") +
  theme_minimal() +
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    axis.text          = element_text(size = 10),
    axis.title         = element_text(size = 15, face = "bold"),
    legend.position    = "none"
  )

## 4) MODEL-BASED marginal means (ZINB, response scale, /100) -----------------

# Load the fitted model if needed
go_brm_model_zinb_new <- readRDS("../brms_models/go_brm_model_zinb_new.rds")

emm_go <- emmeans(
  go_brm_model_zinb_new,
  ~ elevation,
  type       = "response",   # expected index * 100
  re_formula = NA            # population-level effect
)

emm_go_df <- as.data.frame(emm_go)

## Standardise names: mean, LCL, UCL
if (!"mean" %in% names(emm_go_df)) {
  if ("response" %in% names(emm_go_df))      emm_go_df$mean <- emm_go_df$response
  else if ("emmean" %in% names(emm_go_df))   emm_go_df$mean <- emm_go_df$emmean
}

if (!"LCL" %in% names(emm_go_df)) {
  if ("asymp.LCL" %in% names(emm_go_df))     emm_go_df$LCL <- emm_go_df$asymp.LCL
  else if ("lower.CL" %in% names(emm_go_df)) emm_go_df$LCL <- emm_go_df$lower.CL
  else if ("lower.HPD" %in% names(emm_go_df))emm_go_df$LCL <- emm_go_df$lower.HPD
}

if (!"UCL" %in% names(emm_go_df)) {
  if ("asymp.UCL" %in% names(emm_go_df))     emm_go_df$UCL <- emm_go_df$asymp.UCL
  else if ("upper.CL" %in% names(emm_go_df)) emm_go_df$UCL <- emm_go_df$upper.CL
  else if ("upper.HPD" %in% names(emm_go_df))emm_go_df$UCL <- emm_go_df$upper.HPD
}

head(emm_go_df)

## Back-transform from index*100 to original 0–4 scale
emm_go_df <- emm_go_df %>%
  mutate(
    prob = prob / 100,
    LCL  = LCL  / 100,
    UCL  = UCL  / 100
  )

## 5) Compact letter display for elevations -----------------------------------

cld_go_manual <- tibble(
  elevation = factor(c("2300", "2800", "3400", "3800"),
                     levels = levels(go.index$elevation)),
  .group    = c("a", "ab", "ab", "b")
)

y_max_go <- max(go.index$index_raw, na.rm = TRUE)

## 6) Final plot: violins + points + CLD letters + red CrI --------------------

go_plot_with_letters <- go_plot_base +
  geom_text(data = cld_go_manual,
            aes(x = elevation, y = y_max_go + 0.3, label = .group),
            size = 5, fontface = "bold") +
  geom_pointrange(
    data        = emm_go_df,
    aes(x = elevation, y = prob, ymin = LCL, ymax = UCL),
    inherit.aes = FALSE,
    color       = "red",
    size        = 0.6
  )

go_plot_with_letters <- go_plot_with_letters +
  scale_x_discrete(labels = c(
    "2300" = "2,300",
    "2800" = "2,800",
    "3400" = "3,400",
    "3800" = "3,800"
  ))

go_plot_with_letters

# ggsave("figs/fig_GO_bayes.pdf", plot = go_plot_with_letters,
#        width = 8, height = 6)
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
###* COMBINE FIG A-D
add_star_top_right <- function(p, star = "*") {
  p +
    annotate("text",
             x = Inf, y = Inf, label = star,
             hjust = 1.2, vjust = 0.3,   # nudge left/down from corner
             size = 14, fontface = "bold") +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(t = 10, r = 15, b = 5, l = 5))
}

# Panel A
control_combined_plot_with_letters <- control_combined_plot_with_letters +
  labs(tag = "A)") +
  theme(plot.tag = element_text(size = 20, face = "bold"))

control_combined_plot_with_letters <- add_star_top_right(control_combined_plot_with_letters, "**")

# Panel B
pl_plot_with_letters <- pl_plot_with_letters +
  labs(tag = "B)") +
  theme(plot.tag = element_text(size = 20, face = "bold"))

pl_plot_with_letters <- add_star_top_right(pl_plot_with_letters, "**")

# Panel C
ao_plot_with_letters <- ao_plot_with_letters +
  labs(tag = "C)") +
  theme(plot.tag = element_text(size = 20, face = "bold"))

# Panel D
go_plot_with_letters <- go_plot_with_letters +
  labs(tag = "D)") +
  theme(plot.tag = element_text(size = 20, face = "bold"))

go_plot_with_letters <- add_star_top_right(go_plot_with_letters, "*")

# Combine
combined_plots_3 <- (control_combined_plot_with_letters |
  pl_plot_with_letters) /
  (ao_plot_with_letters |
  go_plot_with_letters)

combined_plots_3

ggsave("figs/fig2_final_bayes_3_with_ao.pdf",
       plot = combined_plots_3,
       width = 10, height = 8)
###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* VISITATION METRICS
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
###* FIG 3E - VISITATION
library(tidyverse)
library(emmeans)
library(scales)
library(patchwork)

source("scripts/Bayesian setup June short.R")

final.table$elevation <- as.factor(final.table$elevation)

final.table <- final.table %>%
  mutate(elevation = recode(elevation,
                            `2300` = 2300,
                            `2800` = 2800,
                            `3500` = 3400,
                            `4000` = 3800))

final.table$elevation <- factor(final.table$elevation,
                                levels = c("2300", "2800", "3400", "3800"))

elev_colors <- c(
  "2300" = "#9bcd9b",
  "2800" = "#d4af37",
  "3400" = "#d4c7bd",
  "3800" = "#fdfbf2"
)
names(elev_colors) <- levels(final.table$elevation)

## range & breaks on ORIGINAL scale
vis_range  <- range(final.table$visited.flowers.per.minute, na.rm = TRUE)
y_breaks   <- pretty(c(0, vis_range[2]), n = 5)
y_max_plot <- max(y_breaks)

# helper to undo 0–1 scaling
unscale01 <- function(x, rng = vis_range) rng[1] + x * diff(rng)

# reusable plotting function (now general in y-range)
plot_visitation_metric <- function(data, yvar, ylabel, breaks, filename = NULL) {
  
  p <- ggplot(data, aes(x = elevation, y = .data[[yvar]], fill = elevation)) +
    geom_violin(trim = TRUE, scale = "width", width = 0.7) +
    geom_jitter(width = 0.1, height = 0, size = 1.5, alpha = 0.6, color = "black") +
    scale_fill_manual(values = elev_colors) +
    scale_y_continuous(breaks = seq(0, 0.8, 0.2), limits = c(0, 0.85)) +
    labs(x = "Elevation",
         y = ylabel) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.text  = element_text(size = 10),
      axis.title = element_text(size = 15, face = "bold"),
      legend.position = "none"
    )
  
  if (!is.null(filename)) {
    ggsave(filename, plot = p, width = 8, height = 6)
  }
  
  p
}

## Base plot: ORIGINAL visitation frequency
vis_freq_plot <- plot_visitation_metric(
  data    = final.table,
  yvar    = "visited.flowers.per.minute",
  ylabel  = "Visitation frequency",
  breaks  = y_breaks
)

## Load model (still fitted on 0–1–scaled response)
vis_ord_model <- readRDS("../brms_models/vis_ord_model.rds")

## Emmeans on link scale, back-transform to original units
emm_vis_link <- emmeans(vis_ord_model, ~ elevation)
df_link      <- as.data.frame(emm_vis_link)

mu_col <- intersect(c("response", "emmean", "prob"), names(df_link))[1]
lo_col <- intersect(c("lower.HPD", "lower.CL", "asymp.LCL"), names(df_link))[1]
hi_col <- intersect(c("upper.HPD", "upper.CL", "asymp.UCL"), names(df_link))[1]

emm_vis_df <- df_link %>%
  transmute(
    elevation = factor(elevation,
                       levels = c("2300", "2800", "3400", "3800")),
    mean = unscale01(plogis(.data[[mu_col]])),
    LCL  = unscale01(plogis(.data[[lo_col]])),
    UCL  = unscale01(plogis(.data[[hi_col]]))
  )

## CLD letters
cld_v <- tibble(
  elevation = factor(c("2300", "2800", "3400", "3800"),
                     levels = c("2300", "2800", "3400", "3800")),
  .group = c("b", "c", "bc", "a")
)

## Final visitation-frequency plot on original scale
vis_freq_plot_with_letters <- vis_freq_plot +
  geom_text(data = cld_v,
            aes(x = elevation, y = 0.83, label = .group),
            size = 5, fontface = "bold") +
  geom_hline(yintercept = y_breaks,
             linetype = "dotted", color = "grey70") +
  geom_pointrange(
    data        = emm_vis_df,
    aes(x = elevation, y = mean, ymin = LCL, ymax = UCL),
    inherit.aes = FALSE,
    color       = "red",
    size        = 0.6
  )

vis_freq_plot_with_letters <- vis_freq_plot_with_letters +
  scale_x_discrete(labels = c(
    "2300" = "2,300",
    "2800" = "2,800",
    "3400" = "3,400",
    "3800" = "3,800"
  ))

vis_freq_plot_with_letters
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
###* FIG 3F - MORPHOSPECIES
## Load model
m_brm_model_nb <- readRDS("../brms_models/m_brm_model_nb.rds")

## Base plot (counts; we override y-scale later)
morpho_plot <- plot_visitation_metric(
  data    = final.table,
  yvar    = "total.morpho",
  ylabel  = "Morphospecies richness",
  filename = NULL
)

## Emmeans on RESPONSE scale (expected counts) -------------------------------
emm_morpho <- emmeans(
  m_brm_model_nb,
  ~ elevation,
  type       = "response",
  re_formula = NA      # population-level prediction
)

emm_morpho_df <- as.data.frame(emm_morpho)

# Standardise column names
if (!"mean" %in% names(emm_morpho_df)) {
  if ("response" %in% names(emm_morpho_df))    emm_morpho_df$mean <- emm_morpho_df$response
  else if ("emmean" %in% names(emm_morpho_df)) emm_morpho_df$mean <- emm_morpho_df$emmean
}
if (!"LCL" %in% names(emm_morpho_df)) {
  if ("asymp.LCL" %in% names(emm_morpho_df))   emm_morpho_df$LCL <- emm_morpho_df$asymp.LCL
  else if ("lower.CL" %in% names(emm_morpho_df)) emm_morpho_df$LCL <- emm_morpho_df$lower.CL
  else if ("lower.HPD" %in% names(emm_morpho_df)) emm_morpho_df$LCL <- emm_morpho_df$lower.HPD
}
if (!"UCL" %in% names(emm_morpho_df)) {
  if ("asymp.UCL" %in% names(emm_morpho_df))   emm_morpho_df$UCL <- emm_morpho_df$asymp.UCL
  else if ("upper.CL" %in% names(emm_morpho_df)) emm_morpho_df$UCL <- emm_morpho_df$upper.CL
  else if ("upper.HPD" %in% names(emm_morpho_df)) emm_morpho_df$UCL <- emm_morpho_df$upper.HPD
}

## Manual letters
cld_m <- tibble(
  elevation = factor(c("2300", "2800", "3400", "3800"),
                     levels = c("2300", "2800", "3400", "3800")),
  .group = c("c", "b", "b", "a")
)

## Final morpho plot
morpho_plot_with_letters <- morpho_plot +
  geom_text(data = cld_m,
            aes(x = elevation, y = 33, label = .group),
            size = 5,
            fontface = "bold") +
  scale_y_continuous(
    breaks = seq(0, 30, 5),
    limits = c(0, 33)
  ) +
  geom_hline(yintercept = seq(0, 30, 5),
             linetype = "dotted", color = "grey70") +
  # MODEL-BASED marginal mean ± 95% CrI
  geom_pointrange(
    data        = emm_morpho_df,
    aes(x = elevation, y = prob, ymin = LCL, ymax = UCL),
    inherit.aes = FALSE,
    color       = "red",
    size        = 0.6
  )

morpho_plot_with_letters

morpho_plot_with_letters <- morpho_plot_with_letters +
  scale_x_discrete(labels = c(
    "2300" = "2,300",
    "2800" = "2,800",
    "3400" = "3,400",
    "3800" = "3,800"
  ))

ggsave("figs/fig3.2_bayes.pdf",
       plot = morpho_plot_with_letters,
       width = 8, height = 6)
###*
###*
###*
###*
###*
###*
###*
###*
###*
###* FIG 3F - FUNCTIONAL GROUP
## Load model
f_brm_model_po <- readRDS("../brms_models/f_brm_model_po.rds")

## Base plot
func_plot <- plot_visitation_metric(
  data    = final.table,
  yvar    = "total.func",
  ylabel  = "Functional-group richness",
  filename = NULL
)

## Emmeans on RESPONSE scale (expected counts) -------------------------------
emm_func <- emmeans(
  f_brm_model_po,
  ~ elevation,
  type       = "response",
  re_formula = NA
)



emm_func_df <- as.data.frame(emm_func)

if (!"mean" %in% names(emm_func_df)) {
  if ("response" %in% names(emm_func_df))    emm_func_df$mean <- emm_func_df$response
  else if ("emmean" %in% names(emm_func_df)) emm_func_df$mean <- emm_func_df$emmean
}
if (!"LCL" %in% names(emm_func_df)) {
  if ("asymp.LCL" %in% names(emm_func_df))   emm_func_df$LCL <- emm_func_df$asymp.LCL
  else if ("lower.CL" %in% names(emm_func_df)) emm_func_df$LCL <- emm_func_df$lower.CL
  else if ("lower.HPD" %in% names(emm_func_df)) emm_func_df$LCL <- emm_func_df$lower.HPD
}
if (!"UCL" %in% names(emm_func_df)) {
  if ("asymp.UCL" %in% names(emm_func_df))   emm_func_df$UCL <- emm_func_df$asymp.UCL
  else if ("upper.CL" %in% names(emm_func_df)) emm_func_df$UCL <- emm_func_df$upper.CL
  else if ("upper.HPD" %in% names(emm_func_df)) emm_func_df$UCL <- emm_func_df$upper.HPD
}

## Manual letters
cld_f <- tibble(
  elevation = factor(c("2300", "2800", "3400", "3800"),
                     levels = c("2300", "2800", "3400", "3800")),
  .group = c("c", "b", "b", "a")
)

## Final functional richness plot
func_plot_with_letters <- func_plot +
  geom_text(data = cld_f,
            aes(x = elevation, y = 8.7, label = .group),
            size = 5,
            fontface = "bold") +
  scale_y_continuous(
    breaks = seq(0, 8, 2),
    limits = c(0, 8.7)
  ) +
  geom_hline(yintercept = seq(0, 8, 2),
             linetype = "dotted", color = "grey70") +
  geom_pointrange(
    data        = emm_func_df,
    aes(x = elevation, y = rate, ymin = LCL, ymax = UCL),
    inherit.aes = FALSE,
    color       = "red",
    size        = 0.6
  )

func_plot_with_letters

func_plot_with_letters <- func_plot_with_letters +
  scale_x_discrete(labels = c(
    "2300" = "2,300",
    "2800" = "2,800",
    "3400" = "3,400",
    "3800" = "3,800"
  ))

ggsave("figs/fig3.3_bayes.pdf",
       plot = func_plot_with_letters,
       width = 8, height = 6)
###*
###*
###*
###*
###*
###*
###*
###*
###*
###* COMBINING PLOTS 3D-F
add_star_top_right <- function(p, star = "*") {
  p +
    annotate("text",
             x = Inf, y = Inf, label = star,
             hjust = 1.2, vjust = 0.3,   # nudge left/down from corner
             size = 14, fontface = "bold") +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(t = 10, r = 15, b = 5, l = 5))
}

# Panel A
vis_freq_plot_with_letters <- vis_freq_plot_with_letters +
  theme(plot.tag = element_text(size = 20, face = "bold")) +
  labs(tag = "E)")

vis_freq_plot_with_letters <- add_star_top_right(vis_freq_plot_with_letters, "**")

# Panel B
morpho_plot_with_letters <- morpho_plot_with_letters +
  theme(plot.tag = element_text(size = 20, face = "bold")) +
  labs(tag = "F)")

morpho_plot_with_letters <- add_star_top_right(morpho_plot_with_letters, "**")

# Panel C
func_plot_with_letters <- func_plot_with_letters +
  theme(plot.tag = element_text(size = 20, face = "bold")) +
  labs(tag = "G)")

func_plot_with_letters <- add_star_top_right(func_plot_with_letters, "**")

# Combine
combined_plot <- vis_freq_plot_with_letters +
  morpho_plot_with_letters +
  func_plot_with_letters +
  plot_layout(nrow = 1)

print(combined_plot)

ggsave("figs/fig3_final_bayes.pdf",
       plot = combined_plot,
       width = 12, height = 5)
###*
###*
###*
###*
###*
###*
###*
###* COMBINING ALL PLOTS TOGETHER
combined_plots_3 <- (control_combined_plot_with_letters |
                     pl_plot_with_letters |
                     ao_plot_with_letters |
                     go_plot_with_letters) /
  (vis_freq_plot_with_letters |
  morpho_plot_with_letters |
  func_plot_with_letters ) 
  

combined_plots_3

ggsave("figs/fig3_final_bayes_ALL.pdf",
       plot = combined_plots_3,
       width = 12, height = 9)
###*
###*
###*
###*
###*
###*
###*
###*
###*
###* FIG 4 - GRAPHS FOR CONTROL SEED-SET MODEL WITH VISITATION AND FUNCTIONAL 
###* GROUP AS PREDICTORS
###* 
###* 
###* 
###* 
###* 
###* FIG 4A - PLOT SEEDSET ON VISITATION
source("scripts/05. Setup for index comparison.R")

library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)

#View(flowering.visited)

# your existing code

c.pl.final.table.4 <- c.pl.final.table.4 %>%
  rename(
    z_vis_mean = mean.visitors.scaled,
    z_vis_sd   = sd.visitors.scaled,
    z_flow_mean = mean.visited.flowers.scaled,
    z_flow_sd   = sd.visited.flowers.scaled,
    z_morpho_mean = mean.morpho.scaled,
    z_morpho_sd   = sd.morpho.scaled,
    z_func_mean   = mean.func.scaled,
    z_func_sd     = sd.func.scaled
  )

eps <- 1e-8  # tiny floor to avoid zero SDs

c.pl.final.table.4 <- c.pl.final.table.4 %>%
  mutate(
    # --- means on the z-scale (already computed) ---
    x = z_flow_mean,   # visited flowers per minute (scaled) -- or use z_vis_mean if that's your x
    y = z_morpho_mean,
    z = z_func_mean,
    # --- SEs of the means on the same scale (SD / sqrt(n)) ---
    sx = pmax(z_flow_sd   / sqrt(pmax(n_reps,  1L)), eps),
    sy = pmax(z_morpho_sd  / sqrt(pmax(n_reps, 1L)), eps),
    sz = pmax(z_func_sd   / sqrt(pmax(n_reps,  1L)), eps),
    
    # --- quadratic terms and their delta-method SEs ---
    x2  = x^2,  sx2 = pmax(2 * abs(x) * sx, eps),
    y2  = y^2,  sy2 = pmax(2 * abs(y) * sy, eps),
    z2  = z^2,  sz2 = pmax(2 * abs(z) * sz, eps)
  )

m   <- bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt
dat <- c.pl.final.table.4


vis_mean <- mean(flowering.visited$visitors.per.minute, na.rm = TRUE)
vis_sd   <- sd(flowering.visited$visitors.per.minute,   na.rm = TRUE)

sx_med <- median(dat$sx, na.rm = TRUE)
sz_med <- median(dat$sz, na.rm = TRUE)

xr    <- range(dat$x, na.rm = TRUE)
x_seq <- seq(xr[1], xr[2], length.out = 300)

newdat_vis <- tibble(
  x  = x_seq,    # z-scored visitation (predictor in model)
  sx = sx_med,
  z  = 0,        # functional richness at its mean (z-scale)
  sz = sz_med,
  species = NA
)

pred_vis <- fitted(
  m,
  newdata    = newdat_vis,
  re_formula = NA,
  summary    = TRUE
) |> as.data.frame()

plot_df_vis <- bind_cols(newdat_vis, pred_vis) %>%
  mutate(
    x_visitors_per_min = x * vis_sd + vis_mean   # back-transform x
  )

dat_plot_vis <- c.pl.final.table.4 %>%
  mutate(
    x_visitors_per_min = mean.visitors.per.minute
  )

vis_plot <- ggplot(plot_df_vis, aes(x = x_visitors_per_min, y = Estimate)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5),
              fill = "#d4e9ff", alpha = 0.7, colour = NA) +
  geom_line(size = 1.1, colour = "black") +
  geom_point(data = dat_plot_vis,
             aes(x = x_visitors_per_min, y = mean_seedset_round),
             inherit.aes = FALSE,
             size = 2.2, alpha = 0.85) +
  scale_x_continuous(
    name   = "Visitation frequency",
    breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    trans  = log1p_trans(),
    breaks = c(0, 10, 50, 250, 1000, 5000, 40000),
    labels = label_comma(big.mark = ","),
    name   = "Natural seed set"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "grey90"),
    panel.grid.major.y = element_line(colour = "grey90"),
    axis.title         = element_text(face = "bold", size = 16),
    axis.text          = element_text(size = 13),
    plot.tag           = element_text(size = 18, face = "bold")
  )

vis_plot
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
###* FIG 4B - PLOT SEEDSET ON FUNCITONAL GROUP RICHNESS
# from flowering.visited, before summarising to final.table
m   <- bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt
dat <- c.pl.final.table.4

func_mean <- mean(dat$mean.func, na.rm = TRUE)
func_sd   <- sd(dat$mean.func,   na.rm = TRUE)

# use 0 as the lower bound for x, plus the observed max
func_range <- range(c(0, dat$mean.func), na.rm = TRUE)
func_seq   <- seq(func_range[1], func_range[2], length.out = 300)

# build newdata on ORIGINAL scale, convert to z for the model
newdat_func <- tibble(
  x  = 0,  # visitation at mean (z-scale)
  sx = sx_med,
  z  = (func_seq - func_mean) / func_sd,  # z-scored functional richness
  sz = sz_med,
  species = NA
)

pred_func <- fitted(
  m,
  newdata    = newdat_func,
  re_formula = NA,
  summary    = TRUE
) |> as.data.frame()

plot_df_func <- bind_cols(newdat_func, pred_func) %>%
  mutate(func_richness = func_seq)   # directly use the grid from 0 to max


dat_plot_func <- c.pl.final.table.4 %>%
  mutate(
    func_richness = mean.func   # original functional richness per species×elevation
  )

func_plot <- ggplot(plot_df_func, aes(x = func_richness, y = Estimate)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5),
              fill = "#d4e9ff", alpha = 0.7, colour = NA) +
  geom_line(size = 1.1, colour = "black") +
  geom_point(data = dat_plot_func,
             aes(x = func_richness, y = mean_seedset_round),
             inherit.aes = FALSE,
             size = 2.2, alpha = 0.85) +
  scale_x_continuous(
    name   = "Functional-group richness",
    breaks = pretty(dat_plot_func$func_richness),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    name   = "Natural seed set",
    breaks = c(0, 50, 100, 150, 200, 250, 300, 350),
    limits = c(0, NA)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor   = element_blank(),
    #panel.grid.major.x = element_blank(),
    panel.grid.major.x = element_line(colour = "grey90"),
    panel.grid.major.y = element_line(colour = "grey90"),
    axis.title         = element_text(face = "bold", size = 16),
    axis.text          = element_text(size = 13),
    plot.tag           = element_text(size = 18, face = "bold")
  )

func_plot
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
###* MERGE AND LABEL PLOTS
# 1) Add panel tags
vis_plot_tagged <- vis_plot +
  labs(tag = "A)") +
  theme(plot.tag = element_text(size = 18, face = "bold"))

func_plot_tagged <- func_plot +
  labs(tag = "B)") +
  theme(plot.tag = element_text(size = 18, face = "bold"))

# 2) Combine into one figure
fig5_combined <- vis_plot_tagged + func_plot_tagged +
  plot_layout(nrow = 1)

# Optional: look at it in RStudio
print(fig5_combined)

# 3) Save as a single PDF
ggsave(
  filename = "figs/fig5_combined.pdf",
  plot     = fig5_combined,
  width    = 10,   # tweak as you like
  height   = 5
)






