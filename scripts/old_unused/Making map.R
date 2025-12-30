pacman::p_load(sf, ggplot2, ggspatial, elevatr, dplyr, raster, paletteer)

# 1. Define sites and convert to sf
sites <- data.frame(
  elevation = c("2300 m", "2800 m", "3500 m", "4000 m"),
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

contour_lines <- rasterToContour(elev_raster, levels = seq(0, 4000, by = 100))
contour_sf <- st_as_sf(contour_lines)

# 5. Color palette (same as your plots)
elev_colors <- paletteer_c("grDevices::Green-Yellow", n = 8)

ggplot() +
  # DEM raster
  geom_raster(data = elev_df, aes(x = Longitude, y = Latitude, fill = elevation)) +
  # Contours
  geom_sf(data = contour_sf, color = "grey30", size = 0.2, alpha = 0.6) +
  # MCNP boundary
  #geom_sf(data = mcnp_boundary, fill = NA, color = "black", size = 0.8) +
  # Sampling sites
  geom_sf(data = sites_sf, size = 3, shape = 21, fill = "white", color = "black") +
  # Add Buea point and label
  geom_sf(data = buea_sf, shape = 21, size = 3, fill = "red", color = "white") +
  geom_text(data = buea, aes(x = lon, y = lat, label = "Buea"), inherit.aes = FALSE, 
            hjust = -0.1, vjust = 1.5, size = 4) +

  geom_text(
    data = sites,
    aes(x = lon + 0.015, y = lat, label = elevation),
    size = 4, fontface = "bold"
  ) +
  # Aesthetics
  scale_fill_gradientn(colors = elev_colors, name = "Elevation (m)") +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering()) +
  coord_sf(
    xlim = c(9.10, 9.25),
    ylim = c(4.12, 4.25),
    expand = FALSE
  ) +
  labs(title = "Study Sites in Mount Cameroon National Park") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )



















# # 6. Plot
# ggplot() +
#   geom_raster(data = elev_df, aes(x = x, y = y, fill = elevation)) +
#   geom_sf(data = contour_sf, color = "grey40", size = 0.2, alpha = 0.6) +
#   scale_fill_gradientn(colors = elev_colors) +
#   geom_sf(data = sites_sf, size = 3, shape = 21, fill = "black", color = "white") +
#   geom_text(data = sites,
#             aes(x = lon, y = lat, label = label),
#             hjust = -0.2, vjust = -0.5, size = 3.5, fontface = "bold") +
#   annotation_scale(location = "bl", width_hint = 0.3) +
#   annotation_north_arrow(location = "tl", which_north = "true",
#                          style = north_arrow_fancy_orienteering()) +
#   coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
#            ylim = c(bbox["ymin"], bbox["ymax"]),
#            expand = FALSE) +
#   labs(title = "Study Sites on Mount Cameroon",
#        fill = "Elevation (m)") +
#   theme_minimal() +
#   theme(
#     legend.position = "right",
#     plot.title = element_text(size = 16, face = "bold")
#   )


