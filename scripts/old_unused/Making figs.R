###* Making graphs
###* library(ggplot2)
pacman::p_load(tidyverse, glmmTMB, DHARMa, emmeans, paletteer, emmeans, ggpubr,
               multcompView, sf, ggplot2, ggspatial, elevatr, dplyr, raster,
               paletteer, marginaleffects, scales, multcomp, readr, stringr)

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





library(grid)
library(cowplot)  # or patchwork

# Create the inset map
inset_map <- ggplot() +
  geom_sf(data = cameroon, fill = "white", color = "black", size = 0.3) +
  geom_sf(data = bbox_sf, fill = "red", color = "red", alpha = 1, size = 1) +
  coord_sf(xlim = c(8.4, 17), ylim = c(1.5, 13), expand = FALSE) +
  theme_void()

final_map <- ggdraw() +
  draw_plot(map) +  # your main map object
  draw_plot(inset_map, x = 0.04, y = 0.5, width = 0.3, height = 0.3)  # adjust position and size

final_map

ggsave("figs/final_map_with_Cameroon.pdf", final_map, width = 8, height = 6)







































###* FIGURE 2
###* Load data----

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

View(seed.indices)

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

summary_data <- c.index %>%
  group_by(species, elevation) %>%
  summarise(count = n()) %>%
  ungroup()

# Plot the data
p <- c.index %>%
  ggplot(aes(x = as.factor(elevation), y = seedset, fill = species)) +
  geom_violin(trim = TRUE, scale = "width", width = 0.7,) +
  # geom_boxplot(alpha = 0.1,
  #              outlier.color = "red",
  #              outlier.size = 2,
  #              outlier.shape = 8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
  facet_wrap(. ~ species, nrow = 2, scales = "free_y",
             labeller = label_wrap_gen(width = 50)) +
  labs(x = "Elevation", y = "Control index") +
  theme_minimal() +
  geom_text(data = summary_data, 
            aes(x = as.factor(elevation), y = Inf, label = paste("n =", count)),
            position = position_dodge(width = 0.3),
            vjust = 1, hjust = 0.5, size = 3.5, inherit.aes = FALSE) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.text = element_text( size = 15), #"elevations" size
        axis.title = element_text( size = 20, face = "bold" ), #
        strip.text = element_text(size = 15, face = "bold"))

ggsave("figs/fig2.1.pdf", plot = p, width = 12, height = 7)
















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

# Plot all control seed set indices together
control_combined_plot <- ggplot(c.index, aes(x = elevation, y = seedset, fill = elevation)) +
  geom_violin(trim = TRUE, scale = "width", width = 0.7) +
  # geom_boxplot(alpha = 0.1,
  #              outlier.color = "red",
  #              outlier.size = 2,
  #              outlier.shape = 8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
  # geom_text(data = summary_data_control, 
  #           aes(x = as.factor(elevation), y = y_max_c + 35, label = paste("n =", count)),
  #           vjust = -0.55, hjust = 0.5, size = 5, fontface = "bold", inherit.aes = FALSE) +
  labs(x = "Elevation", y = "Natural Seed-Set") +
  scale_fill_manual(values = elev_colors) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20, face = "bold"),
    legend.position = "none"
  )

#control_combined_plot

c_model <- readRDS("glm_outputs/c_model.rds")

# Extract estimated marginal means and post-hoc comparisons
c_emm <- emmeans(c_model, ~ elevation)
pairs(c_emm)
summary(c_model)
# Generate compact letter display
cld_c <- multcomp::cld(c_emm, adjust = "tukey", Letters = letters) %>%
  as.data.frame() %>%
  mutate(
    elevation = factor(elevation, levels = c("2300", "2800", "3400", "3800")),
    .group = str_trim(.group)
  )

cld_c

# Get max y-value to position letters
y_max_c <- max(c.index$seedset, na.rm = TRUE)

# Plot with letters
control_combined_plot_with_letters <- control_combined_plot +
  geom_text(data = cld_c,
            aes(x = elevation, y = y_max_c + 40, label = .group),
            size = 7,
            fontface = "bold")

control_combined_plot_with_letters

# Save the final figure
ggsave("figs/fig2.2.pdf",
       plot = control_combined_plot_with_letters,
       width = 8, height = 6)

control_combined_plot_with_letters











###* FIGURE 2.3
# Summary n per elevation
summary_data_pl <- c.index %>%
  group_by(elevation) %>%
  summarise(count = n())

elev_colors <- c(
  "2300" = "#9bcd9b",
  "2800" = "#d4af37",
  "3400" = "#d4c7bd",
  "3800" = "#fdfbf2"
)

c.index$elevation <- factor(c.index$elevation, levels = c("2300", "2800", "3400", "3800"))

# Assign the extracted palette
names(elev_colors) <- levels(c.index$elevation)


# Plot
pl_plot <- ggplot(c.index, aes(x = elevation, y = PL.index, fill = elevation)) +
  geom_violin(trim = TRUE, scale = "width", width = 0.7) +
  # geom_boxplot(alpha = 0.1,
  #              outlier.color = "red",
  #              outlier.size = 2,
  #              outlier.shape = 8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
  # geom_text(data = summary_data_pl, 
  #           aes(x = elevation, y = Inf, label = paste("n =", count)),
  #           vjust = 1, hjust = 0.5, size = 5, fontface = "bold", inherit.aes = FALSE) +
  labs(x = "Elevation", y = "Pollen Limitation") +
  scale_fill_manual(values = elev_colors) +
  #scale_fill_brewer(palette = "YlGnBu") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20, face = "bold"),
    strip.text = element_text(size = 15),
    legend.position = "none"
  )

pl_plot

# Load model
pl_model <- readRDS("glm_outputs/pl_model.rds")

# Tukey post-hoc comparisons
pl_emm <- emmeans(pl_model, ~ elevation)

# Generate compact letter display
cld_pl <- multcomp::cld(pl_emm, adjust = "tukey", Letters = letters) %>%
  as.data.frame() %>%
  mutate(elevation = factor(elevation, levels = c("2300", "2800", "3400", "3800")),
         .group = str_trim(.group))

# Get max y value to position letters above violins
y_max_pl <- max(c.index$PL.index, na.rm = TRUE)

# Add letters to the plot
pl_plot_with_letters <- pl_plot +
  geom_text(data = cld_pl,
            aes(x = elevation, y = y_max_pl + 0.07, label = .group),
            size = 7,
            fontface = "bold")

pl_plot_with_letters

# Save the updated plot
ggsave("figs/fig2.3.pdf", plot = pl_plot_with_letters, width = 8, height = 6)










# Tag and style each plot manually
control_combined_plot_with_letters <- control_combined_plot_with_letters +
  theme(plot.tag = element_text(size = 25, face = "bold")) +
  labs(tag = "A)")

pl_plot_with_letters <- pl_plot_with_letters +
  theme(plot.tag = element_text(size = 25, face = "bold")) +
  labs(tag = "B)")

# Combine plots into one figure
combined_plots_2 <- control_combined_plot_with_letters +
  pl_plot_with_letters +
  plot_layout(ncol = 2)  # stacked vertically

# Display in RStudio
print(combined_plots_2)

# Save
ggsave("figs/fig2_final.pdf",
       plot = combined_plots_2,
       width = 12, height = 6)






















###* FIGURE 3
source("scripts/Bayesian setup June short.R")

#View(final.table)

final.table$elevation<-as.factor(final.table$elevation)

final.table <- final.table %>%
  mutate(elevation = recode(elevation,
                            `2300` = 2300,
                            `2800` = 2800,
                            `3500` = 3400,
                            `4000` = 3800))

final.table <- final.table %>%
  ungroup() %>%  # Ensure no grouping is active
  mutate(visited.0.1.scaled = rescale(visited.flowers.per.minute, to = c(0, 1), na.rm = TRUE))

# Factor and color setup
final.table$elevation <- factor(final.table$elevation, levels = c("2300", "2800", "3400", "3800"))
#elev_colors <- paletteer_c("grDevices::Green-Yellow", n = 4)
elev_colors <- c(
  "2300" = "#9bcd9b",
  "2800" = "#d4af37",
  "3400" = "#d4c7bd",
  "3800" = "#fdfbf2"
)
names(elev_colors) <- levels(final.table$elevation)

# Reusable function to generate sample size summary per elevation
get_summary_n <- function(df) {
  df %>%
    group_by(elevation) %>%
    summarise(count = n(), .groups = "drop")
}

plot_visitation_metric <- function(data, yvar, ylabel, filename) {
  summary_data <- get_summary_n(data)
  
  ggplot(data, aes(x = elevation, y = .data[[yvar]], fill = elevation)) +
    geom_violin(trim = TRUE, scale = "width", width = 0.7) +
    # geom_boxplot(alpha = 0.1,
    #              outlier.color = "red",
    #              outlier.size = 2,
    #              outlier.shape = 8) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
    # geom_text(data = summary_data, 
    #           aes(x = elevation, y = Inf, label = paste("n =", count)),
    #           vjust = 1, hjust = 0.5, size = 5, fontface = "bold", inherit.aes = FALSE) +
    scale_fill_manual(values = elev_colors) +
    labs(x = "Elevation", y = ylabel) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 20, face = "bold"),
      legend.position = "none"
    ) -> plot_out
  
  # Save plot
  ggsave(filename, plot = plot_out, width = 8, height = 6)
  
  return(plot_out)
}

# Generate the plot object for visitation frequency
vis_freq_plot <- plot_visitation_metric(
  data = final.table,
  yvar = "visited.0.1.scaled",
  ylabel = "Visitation frequency (scaled)",
  filename = "figs/visitation_frequency.pdf"
)

v_model <- readRDS("glm_outputs/v_model.rds")

# Tukey post-hoc comparisons
v_emm <- emmeans(v_model, ~ elevation)

# Generate compact letter display
cld_letters <- multcomp::cld(v_emm, adjust = "tukey", Letters = letters) %>%
  as.data.frame() %>%
  mutate(
    elevation = as.factor(elevation),
    .group = str_trim(.group)
  )

v_contrasts_all <- as.data.frame(pairs(v_emm, adjust = "tukey")) %>%
  mutate(
    group1 = str_extract(contrast, "^[0-9]+"),
    group2 = str_extract(contrast, "(?<=-)\\d+"),
    label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "n.s."
    ),
    y.position = 1.1 + 0.05 * row_number()
  ) %>%
  select(group1, group2, y.position, label)


y_max <- max(final.table$visited.0.1.scaled, na.rm = TRUE)

vis_freq_plot_with_letters <- vis_freq_plot +
  geom_text(data = cld_letters,
            aes(x = elevation, y = y_max + 0.1, label = .group),
            size = 6,
            fontface = "bold")


vis_freq_plot_with_letters

ggsave("figs/fig3.1.pdf",
       plot = vis_freq_plot_with_letters,
       width = 8, height = 6
)

















###* FIGURE 3.2
# Morphospecies richness
m_model <- readRDS("glm_outputs/m_model.rds")

# Tukey post-hoc comparisons for morphospecies richness
m_emm <- emmeans(m_model, ~ elevation)

# Compact letter display
cld_morpho <- multcomp::cld(m_emm, adjust = "tukey", Letters = letters) %>%
  as.data.frame() %>%
  mutate(
    elevation = as.factor(elevation),
    .group = str_trim(.group)
  )

# Get max y value for label placement
y_max_morpho <- max(final.table$total.morpho, na.rm = TRUE)

# Base plot
morpho_plot <- plot_visitation_metric(
  data = final.table,
  yvar = "total.morpho",
  ylabel = "Morphospecies richness",
  filename = "figs/morphospecies_richness.pdf"  # Optional save of base plot
)

morpho_plot_with_letters <- morpho_plot +
  geom_text(data = cld_morpho,
            aes(x = elevation, y = y_max_morpho + 2, label = .group),
            size = 6,
            fontface = "bold")

morpho_plot_with_letters

ggsave("figs/fig3.2.pdf",
       plot = morpho_plot_with_letters,
       width = 8, height = 6
)











###* FIGURE 3.3
# Functional group richness
f_model <- readRDS("glm_outputs/f_model.rds")

# Tukey post-hoc comparisons for morphospecies richness
f_emm <- emmeans(f_model, ~ elevation)

# Compact letter display
cld_func <- multcomp::cld(f_emm, adjust = "tukey", Letters = letters) %>%
  as.data.frame() %>%
  mutate(
    elevation = as.factor(elevation),
    .group = str_trim(.group)
  )

# Get max y value for label placement
y_max_func <- max(final.table$total.func, na.rm = TRUE)

# Base plot
func_plot <- plot_visitation_metric(
  data = final.table,
  yvar = "total.func",
  ylabel = "Functional group richness",
  filename = "figs/functional_richness.pdf"  # Optional save of base plot
)

func_plot_with_letters <- func_plot +
  geom_text(data = cld_func,
            aes(x = elevation, y = 10, label = .group),
            size = 6,
            fontface = "bold") +
  coord_cartesian(ylim = c(0, 10))  # Set y-axis range

func_plot_with_letters

ggsave("figs/fig3.3.pdf",
       plot = func_plot_with_letters,
       width = 8, height = 6
)




# Load patchwork for combining plots
library(patchwork)

vis_freq_plot_with_letters <- vis_freq_plot_with_letters +
  theme(plot.tag = element_text(size = 25, face = "bold")) +
  labs(tag = "C)")

morpho_plot_with_letters <- morpho_plot_with_letters +
  theme(plot.tag = element_text(size = 25, face = "bold")) +
  labs(tag = "D)")

func_plot_with_letters <- func_plot_with_letters +
  theme(plot.tag = element_text(size = 25, face = "bold")) +
  labs(tag = "E)")

combined_plot <- vis_freq_plot_with_letters +
  morpho_plot_with_letters +
  func_plot_with_letters +
  plot_layout(nrow = 1)

# Display in RStudio
print(combined_plot)

# Save to PDF
ggsave("figs/fig3_final.pdf",
       plot = combined_plot,
       width = 12, height = 6)

combined_plot




# Total number of visits
total_visits <- sum(final.table$total.visitor.count, na.rm = TRUE)

# Number of visits identified to morphospecies level
morpho_visits <- sum(final.table$morpho.visitor.count, na.rm = TRUE)

# Output the results
cat("We recorded", total_visits, "floral visits in total,", 
    morpho_visits, "of which were identified to the morphospecies level and used for analysis.\n")















source("scripts/Bayesian setup June.R")
bayesian_ao_noel_w12_p1_f <- readRDS("../brms_models/bayesian_ao_noel_w12_p1_f.rds")
fit <- bayesian_ao_noel_w12_p1_f

# Extract predictor range from model data
x_range <- insight::get_data(fit) %>%
  pull(mean.func.scaled)
# Generate marginal predictions across the predictor range
data_pred_indiv <- predictions(
  fit,
  newdata = datagrid(
    mean.func.scaled = seq(from = min(x_range), to = max(x_range), length.out = 100)
  ))
# Extract the original observed data
observed_data <- insight::get_data(fit)
# Create ggplot
bayes_ao <- ggplot(data = as.data.frame(data_pred_indiv), 
       aes(x = mean.func.scaled, y = estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "steelblue") + # Uncertainty ribbon
  geom_line(color = "steelblue", size = 1) + # Fitted line
  geom_point(data = observed_data, # Observed data points
             aes(x = mean.func.scaled, y = mean_ao_index), 
             color = "black", alpha = 0.6, size = 3) +
  labs(
    #title = "Predicted AO index vs. func",
    x = "Number of functional groups (scaled)",
    y = "Predicted AO index"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(size = 15),       # Axis tick labels
    axis.title = element_text(size = 20, face = "bold"),  # Axis titles
    legend.position = "none"
  )

# Save to PDF
ggsave("figs/fig4_final.pdf",
       plot = bayes_ao,
       width = 8, height = 6)
