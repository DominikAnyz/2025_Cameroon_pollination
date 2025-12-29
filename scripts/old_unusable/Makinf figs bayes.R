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





































### FIgs 2
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
  scale_y_log10() +
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






control_combined_plot <- ggplot(c.index, aes(x = elevation, y = seedset, fill = elevation)) +
  geom_violin(trim = TRUE, scale = "width", width = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
  labs(x = "Elevation", y = "Natural Seed-Set") +
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
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20, face = "bold"),
    legend.position = "none"
  )


# Load Bayesian model
c_brm_model_zinb <- readRDS("../brms_models/c_brm_model_zinb.rds")

# posterior draws for pairwise elevation contrasts
contr_draws <- emmeans(c_brm_model_zinb, ~ elevation) %>%
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
            aes(x = elevation, y = 800, label = .group),
            #aes(x = elevation, y = y_max_c * 1.6, label = .group), # 20% above max
            size = 7,
            fontface = "bold") +
  stat_summary(fun = median, geom = "point",
               shape = 18, size = 4, color = "red") +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.5, color = "red", fatten = 2)

# Save the figure
ggsave("figs/fig2.2_bayes.pdf",
       plot = control_combined_plot_with_letters,
       width = 8, height = 6)

control_combined_plot_with_letters


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

control_combined_plot_with_letters +
  geom_pointrange(data = emm_df,
                  aes(x = elevation, y = prob, ymin = LCL, ymax = UCL),
                  inherit.aes = FALSE, color = "red", size = 0.6)






















###* FIGURE 2.3
# Summary n per elevation
###POLLEN LIMITATION
c.index <- c.index %>%
  filter(species != "Hypericum r" | elevation != 3800)

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
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1.1)) +
  geom_hline(yintercept = seq(0, 1, 0.25),
             linetype = "dotted", color = "grey70") +
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
pl_ord_model <- readRDS("../brms_models/pl_ord_model_new.rds")

# posterior draws for pairwise elevation contrasts
contr_draws <- emmeans(pl_ord_model_new, ~ elevation) %>%
  contrast(method = "pairwise") %>%
  gather_emmeans_draws()

# Summarise contrasts
post_probs <- contr_draws %>%
  group_by(contrast) %>%
  summarise(
    Estimate   = median(.value),
    l95        = quantile(.value, 0.025),
    u95        = quantile(.value, 0.975),
    P_greater0 = mean(.value > 0),
    P_less0    = mean(.value < 0),
    .groups = "drop"
  ) %>%
  mutate(contrast_clean = gsub("elevation", "", contrast))

post_probs

# Manually assign letters (lowest = a, mid = b, highest = c)
cld_pl <- tibble(
  elevation = factor(c("2300", "2800", "3400", "3800"),
                     levels = c("2300", "2800", "3400", "3800")),
  .group = c("a", "a", "a", "b")  # your chosen letters
)

# Get max y value to position letters above violins
y_max_pl <- max(c.index$PL.index, na.rm = TRUE)

# Add letters to the plot
pl_plot_with_letters <- pl_plot +
  geom_text(data = cld_pl,
            aes(x = elevation, y = y_max_pl + 0.09, label = .group),
            size = 7,
            fontface = "bold") +
  stat_summary(fun = median, geom = "point", 
               shape = 18, size = 4, color = "red") +
  stat_summary(fun = median, geom = "crossbar", 
               width = 0.5, color = "red", fatten = 2)

pl_plot_with_letters



# Marginal means on the response scale (i.e., 0–1, including any model links)
emm_pl <- emmeans(pl_ord_model_new, ~ elevation, type = "response")

# Full posterior draws of those means
pl_mean_draws <- gather_emmeans_draws(emm_pl)   # columns: elevation, .draw, .value

# Summarise for plotting or tables
pl_mean_summ <- pl_mean_draws %>%
  group_by(elevation) %>%
  summarise(
    mean = mean(.value),
    median = median(.value),
    lo = quantile(.value, 0.025),
    hi = quantile(.value, 0.975),
    .groups = "drop"
  )

pl_diff_draws <- contrast(emm_pl, method = "pairwise", adjust = "none") %>%
  gather_emmeans_draws()                         # columns: contrast, .draw, .value

pl_diff_summ <- pl_diff_draws %>%
  group_by(contrast) %>%
  summarise(
    mean = mean(.value),
    median = median(.value),
    lo = quantile(.value, 0.025),
    hi = quantile(.value, 0.975),
    p_gt0 = mean(.value > 0),
    p_lt0 = mean(.value < 0),
    .groups = "drop"
  )

# get emmeans on the link (default), then regrid to the response scale
emm_pl_link <- emmeans(pl_ord_model_new, ~ elevation)     # link scale
df <- as.data.frame(emm_pl_link)

# figure out the columns present
mu_col  <- intersect(c("response","emmean","prob"), names(df))[1]
l_col   <- intersect(c("lower.HPD","lower.CL","asymp.LCL"), names(df))[1]
u_col   <- intersect(c("upper.HPD","upper.CL","asymp.UCL"), names(df))[1]

emm_pl_df <- transform(df,
                       mean = plogis(df[[mu_col]]),
                       LCL  = plogis(df[[l_col]]),
                       UCL  = plogis(df[[u_col]])
)
emm_pl_df <- emm_pl_df[, c("elevation","mean","LCL","UCL")]

# standardize column names for plotting
nm <- names(emm_pl_df)
if ("response" %in% nm) names(emm_pl_df)[nm=="response"] <- "mean"
if ("emmean"   %in% nm) names(emm_pl_df)[nm=="emmean"]   <- "mean"
if ("prob"     %in% nm) names(emm_pl_df)[nm=="prob"]     <- "mean"
if ("lower.HPD" %in% nm) names(emm_pl_df)[nm=="lower.HPD"] <- "LCL"
if ("upper.HPD" %in% nm) names(emm_pl_df)[nm=="upper.HPD"] <- "UCL"
if ("lower.CL"  %in% nm) names(emm_pl_df)[nm=="lower.CL"]  <- "LCL"
if ("upper.CL"  %in% nm) names(emm_pl_df)[nm=="upper.CL"]  <- "UCL"

# Now emm_pl_df$mean, LCL, UCL are all in [0,1] and safe to plot
# Overlay on your existing pl_plot (remove raw median crossbars if you want a model-only figure)
pl_plot_with_letters +
  geom_pointrange(data = emm_pl_df,
                  aes(x = elevation, y = mean, ymin = LCL, ymax = UCL),
                  inherit.aes = FALSE, color = "red", size = 0.6)


names(emm_pl_df); head(emm_pl_df)

















# --- Pollen Limitation figure (ordbetareg, Bayesian) -------------------------

library(ggplot2)
library(dplyr)
library(emmeans)
library(tidybayes)  # only used if you want draws later

# Colors + factor order (match your other figs)
elev_colors <- c("2300"="#9bcd9b","2800"="#d4af37","3400"="#d4c7bd","3800"="#fdfbf2")
c.index$elevation <- factor(c.index$elevation, levels = c("2300","2800","3400","3800"))
names(elev_colors) <- levels(c.index$elevation)

# Base violin + points
pl_plot <- ggplot(c.index, aes(x = elevation, y = PL.index, fill = elevation)) +
  geom_violin(trim = TRUE, scale = "width", width = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
  labs(x = "Elevation", y = "Pollen Limitation") +
  scale_fill_manual(values = elev_colors) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1.1)) +
  geom_hline(yintercept = seq(0, 1, 0.25), linetype = "dotted", color = "grey70") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20, face = "bold"),
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

pl_plot_with_letters <- pl_plot +
  geom_text(data = cld_pl,
            aes(x = elevation, y = y_max_pl + 0.09, label = .group),
            size = 7, fontface = "bold") +
  # model-based marginal mean (dot) + 95% CrI (vertical line)
  geom_pointrange(data = emm_pl_df,
                  aes(x = elevation, y = mean, ymin = LCL, ymax = UCL),
                  inherit.aes = FALSE, color = "red", size = 0.6)

print(pl_plot_with_letters)

ggsave("figs/fig_PL_bayes.pdf", plot = pl_plot_with_letters, width = 8, height = 6)

# --- (Optional) Posterior draws & pairwise differences on the response scale --
# Full posterior draws of marginal means (useful for ribbons/densities)
pl_draws_resp <- regrid(emm_pl_link, transform = "response") |>
  gather_emmeans_draws()  # columns: elevation, .draw, .value

# Pairwise differences with posterior draws (response scale)
pl_diff_draws <- contrast(regrid(emm_pl_link, transform = "response"),
                          method = "pairwise", adjust = "none") |>
  gather_emmeans_draws()

# Summaries table if you need it:
pl_diff_summ <- pl_diff_draws %>%
  group_by(contrast) %>%
  summarise(
    mean = mean(.value),
    median = median(.value),
    lo = quantile(.value, 0.025),
    hi = quantile(.value, 0.975),
    p_gt0 = mean(.value > 0),
    .groups = "drop"
  )

pl_diff_summ




















###* GEITONOGAMY
go.index <- seed.indices %>%
  filter(treatment == "geitonogamy") %>%
  select(index, elevation, species, plant_number) %>%
  mutate(elevation = factor(elevation),
         species = factor(species)) %>%
  mutate(index = round(index * 100)) %>%
  mutate(plant.id = as.factor(paste0(elevation,species,plant_number)))%>%
  filter(species != "Lactuca i") %>%
  filter(species != "Hypericum r" | elevation != 3800)

View(go.index)

go_brm_model_zinb_new <- readRDS("../brms_models/go_brm_model_zinb_new.rds")

# (Optional) posterior draws for pairwise elevation contrasts
go_contr_draws <- emmeans(go_brm_model_zinb_new, ~ elevation) %>%
  contrast(method = "pairwise") %>%
  gather_emmeans_draws()

go_post_probs <- go_contr_draws %>%
  group_by(contrast) %>%
  summarise(
    Estimate   = median(.value),
    l95        = quantile(.value, 0.025),
    u95        = quantile(.value, 0.975),
    P_greater0 = mean(.value > 0),
    P_less0    = mean(.value < 0),
    .groups = "drop"
  )
go_post_probs  # for Table SXXX1 if needed

# Geitonogamy plot
go_plot <- ggplot(go.index, aes(x = elevation, y = index, fill = elevation)) +
  geom_violin(trim = TRUE, scale = "width", width = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
  labs(x = "Elevation", y = "Geitonogamy index") +
  scale_fill_manual(values = elev_colors) +
  scale_y_continuous(breaks = seq(0, 400, by =100), 
                     limits = c(0, 410)) +
  geom_hline(yintercept = seq(0, 1, 0.25),
             linetype = "dotted", color = "grey70") +
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

# Manual letters for geitonogamy
cld_go <- tibble(
  elevation = factor(c("2300","2800","3400","3800"),
                     levels = c("2300","2800","3400","3800")),
  .group = c("a","ab","ab","b")
)
y_max_go <- max(go.index$index, na.rm = TRUE)

go_plot_with_letters <- go_plot +
  geom_text(data = cld_go,
            aes(x = elevation, y = y_max_go + 0.09, label = .group),
            size = 7, fontface = "bold") +
  stat_summary(fun = mean, geom = "point",
               shape = 18, size = 4, color = "red") +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.5, color = "red", fatten = 2) +
  theme(plot.tag = element_text(size = 25, face = "bold")) +
  labs(tag = "C)")

go_plot_with_letters









library(patchwork)

# Tag and style each plot manually
control_combined_plot_with_letters <- control_combined_plot_with_letters +
  theme(plot.tag = element_text(size = 25, face = "bold")) +
  labs(tag = "A)")

pl_plot_with_letters <- pl_plot_with_letters +
  theme(plot.tag = element_text(size = 25, face = "bold")) +
  labs(tag = "B)")

# Combine plots into one figure (side by side)
combined_plots_2 <- control_combined_plot_with_letters | pl_plot_with_letters

# Display in RStudio
print(combined_plots_2)

# Save
ggsave("figs/fig2_final_bayes_2.pdf",
       plot = combined_plots_2,
       width = 12, height = 6)













library(emmeans)

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
    geom_jitter(width = 0.1, height = 0, size = 1.5, alpha = 0.6, color = "black") +
    # geom_text(data = summary_data, 
    #           aes(x = elevation, y = Inf, label = paste("n =", count)),
    #           vjust = 1, hjust = 0.5, size = 5, fontface = "bold", inherit.aes = FALSE) +
    scale_fill_manual(values = elev_colors) +
    # geom_hline(yintercept = seq(0, 1, 0.25),
    #            linetype = "dotted", color = "grey70") +
    scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1.1)) +
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

vis_ord_model <- readRDS("../brms_models/vis_ord_model.rds")

# posterior draws for pairwise elevation contrasts
contr_draws <- emmeans(vis_ord_model, ~ elevation) %>%
  contrast(method = "pairwise") %>%
  gather_emmeans_draws()

# Summarise contrasts
post_probs <- contr_draws %>%
  group_by(contrast) %>%
  summarise(
    Estimate   = median(.value),
    l95        = quantile(.value, 0.025),
    u95        = quantile(.value, 0.975),
    P_greater0 = mean(.value > 0),
    P_less0    = mean(.value < 0),
    .groups = "drop"
  ) %>%
  mutate(contrast_clean = gsub("elevation", "", contrast))

post_probs

# Manually assign letters (lowest = a, mid = b, highest = c)
cld_v <- tibble(
  elevation = factor(c("2300", "2800", "3400", "3800"),
                     levels = c("2300", "2800", "3400", "3800")),
  .group = c("b", "c", "bc", "a")  # your chosen letters
)

y_max <- max(final.table$visited.0.1.scaled, na.rm = TRUE)

vis_freq_plot_with_letters <- vis_freq_plot +
  geom_text(data = cld_v,
            aes(x = elevation, y = y_max + 0.1, label = .group),
            size = 6,
            fontface = "bold") +
  geom_hline(yintercept = seq(0, 1, 0.25),
             linetype = "dotted", color = "grey70") +
  stat_summary(fun = median, geom = "point", 
               shape = 18, size = 4, color = "red") +
  stat_summary(fun = median, geom = "crossbar", 
               width = 0.5, color = "red", fatten = 2)

vis_freq_plot_with_letters

ggsave("figs/fig3.1.pdf",
       plot = vis_freq_plot_with_letters,
       width = 8, height = 6
)











###* FIGURE 3.2
# Morphospecies richness
m_brm_model_nb <- readRDS("../brms_models/m_brm_model_nb.rds")

contr_draws <- emmeans(m_brm_model_nb, ~ elevation) %>%
  contrast(method = "pairwise") %>%
  gather_emmeans_draws()

# Summarise contrasts
post_probs <- contr_draws %>%
  group_by(contrast) %>%
  summarise(
    Estimate   = median(.value),
    l95        = quantile(.value, 0.025),
    u95        = quantile(.value, 0.975),
    P_greater0 = mean(.value > 0),
    P_less0    = mean(.value < 0),
    .groups = "drop"
  ) %>%
  mutate(contrast_clean = gsub("elevation", "", contrast))

post_probs

# Manually assign letters (lowest = a, mid = b, highest = c)
cld_m <- tibble(
  elevation = factor(c("2300", "2800", "3400", "3800"),
                     levels = c("2300", "2800", "3400", "3800")),
  .group = c("c", "b", "b", "a")  # your chosen letters
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
  geom_text(data = cld_m,
            aes(x = elevation, 
                y = 33, 
                label = .group),
            size = 6,
            fontface = "bold") +
  scale_y_continuous(
    breaks = seq(0, 30, 5),
    limits = c(0, 33)
  ) +
  geom_hline(yintercept = seq(0, 30, 5),
             linetype = "dotted", color = "grey70") +
  stat_summary(fun = median, geom = "point", 
               shape = 18, size = 4, color = "red") +
  stat_summary(fun = median, geom = "crossbar", 
               width = 0.5, color = "red", fatten = 2)

morpho_plot_with_letters

ggsave("figs/fig3.2.pdf",
       plot = morpho_plot_with_letters,
       width = 8, height = 6
)













###* FIGURE 3.3
# Functional group richness
f_brm_model_po <- readRDS("../brms_models/f_brm_model_po.rds")

contr_draws <- emmeans(f_brm_model_po, ~ elevation) %>%
  contrast(method = "pairwise") %>%
  gather_emmeans_draws()

# Summarise contrasts
post_probs <- contr_draws %>%
  group_by(contrast) %>%
  summarise(
    Estimate   = median(.value),
    l95        = quantile(.value, 0.025),
    u95        = quantile(.value, 0.975),
    P_greater0 = mean(.value > 0),
    P_less0    = mean(.value < 0),
    .groups = "drop"
  ) %>%
  mutate(contrast_clean = gsub("elevation", "", contrast))

post_probs

# Manually assign letters (lowest = a, mid = b, highest = c)
cld_f <- tibble(
  elevation = factor(c("2300", "2800", "3400", "3800"),
                     levels = c("2300", "2800", "3400", "3800")),
  .group = c("c", "b", "b", "a")  # your chosen letters
)

# Get max y value for label placement
y_max_func <- max(final.table$total.func, na.rm = TRUE)

# Base plot
func_plot <- plot_visitation_metric(
  data = final.table,
  yvar = "total.func",
  ylabel = "Functional-group richness",
  filename = "figs/functional_richness.pdf"  # Optional save of base plot
)

func_plot_with_letters <- func_plot +
  geom_text(data = cld_f,
            aes(x = elevation,
                y = 8.7,
                label = .group),
            size = 6,
            fontface = "bold") +
  scale_y_continuous(
      breaks = seq(0, 8, 2),
      limits = c(0, 8.7)
    ) +
  geom_hline(yintercept = seq(0, 8, 2),
             linetype = "dotted", color = "grey70") +
  stat_summary(fun = median, geom = "point",
               shape = 18, size = 4, color = "red") +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.5, color = "red", fatten = 2)

func_plot_with_letters

ggsave("figs/fig3.3.pdf",
       plot = func_plot_with_letters,
       width = 8, height = 6
)














vis_freq_plot_with_letters <- vis_freq_plot_with_letters +
  theme(plot.tag = element_text(size = 25, face = "bold")) +
  labs(tag = "D)")

morpho_plot_with_letters <- morpho_plot_with_letters +
  theme(plot.tag = element_text(size = 25, face = "bold")) +
  labs(tag = "E)")

func_plot_with_letters <- func_plot_with_letters +
  theme(plot.tag = element_text(size = 25, face = "bold")) +
  labs(tag = "F)")

combined_plot <- vis_freq_plot_with_letters +
  morpho_plot_with_letters +
  func_plot_with_letters +
  plot_layout(nrow = 1)

# Display in RStudio
print(combined_plot)

# Save to PDF
ggsave("figs/fig3_final_bayes.pdf",
       plot = combined_plot,
       width = 12, height = 6)

combined_plot





















source("scripts/Bayesian setup June.R")
bayesian_ao_noel_w12_p1_f <- readRDS("../brms_models/bayesian_ao_noel_w12_p1_f.rds")
fit <- bayesian_ao_noel_w12_p1_f

# Extract predictor range from model data
x_range <- insight::get_data(fit) %>%
  pull(mean.func.scaled)

View(final.table)

func_mean <- mean(final.table$mean.func, na.rm = TRUE)
func_sd   <- sd(final.table$mean.func, na.rm = TRUE)


# Generate marginal predictions across the predictor range
data_pred_indiv <- predictions(
  fit,
  newdata = datagrid(
    mean.func.scaled = seq(from = min(x_range), to = max(x_range), length.out = 100)
  ))
# Extract the original observed data
observed_data <- insight::get_data(fit)



###* Transform back
data_pred_indiv <- data_pred_indiv %>%
  mutate(mean.func = mean.func.scaled * func_sd + func_mean)

observed_data <- observed_data %>%
  mutate(mean.func = mean.func.scaled * func_sd + func_mean)



# Create ggplot
bayes_ao <- ggplot(data = as.data.frame(data_pred_indiv), 
                   aes(x = mean.func, y = estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "steelblue") +
  geom_line(color = "steelblue", size = 1) +
  geom_point(data = observed_data,
             aes(x = mean.func, y = mean_ao_index),
             color = "black", alpha = 0.6, size = 3) +
  labs(
    x = "Number of functional groups",
    y = "Predicted AO index"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20, face = "bold"),
    legend.position = "none"
  )

bayes_ao

# Save to PDF
ggsave("figs/fig4_final_fixed_x.pdf",
       plot = bayes_ao,
       width = 8, height = 6)




