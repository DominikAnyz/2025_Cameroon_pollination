###* Tables for the manuscript
###* 
###* 
###* 
###* 
###* Table 1 - Plant species summary
###* 
###* 
###* 
###* 
###* I would like to create a table in which I have information about the amounts 
###* of observations from distinct treatmen
# Load packages
pacman::p_load(tidyverse, emmeans, gt)

# Ensure dplyr::select is preferred
select <- dplyr::select

# Set seed for reproducibility
set.seed(1234)

# Load and preprocess data
seed.data <- read.delim("data/clean_seeds.txt", na = c("na"))

# Calculate indices
seed.indices <- 
  seed.data %>%
  group_by(elevation, plant_number, species) %>%
  mutate(O.mean.plantnumber = mean(seedset[treatment == "outcrossing"], na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(elevation, species) %>%
  mutate(O.mean.elevation = mean(seedset[treatment == "outcrossing"], na.rm = TRUE)) %>%
  mutate(O.mean = case_when(
    is.nan(O.mean.plantnumber) ~ O.mean.elevation,
    .default = O.mean.plantnumber
  )) %>%
  mutate(index = seedset / O.mean) %>%
  mutate(seedset = if_else(seedset != round(seedset), round(seedset), seedset))

# Clean index values and structure
seed.indices <- seed.indices %>%
  filter(!(index %in% c("NA", "Inf") | is.na(seedset))) %>%
  mutate(index = ifelse(is.finite(index), index, 0))

# Ensure factor variables
seed.indices$elevation <- as.factor(seed.indices$elevation)
seed.indices$species <- as.factor(seed.indices$species)

# Add elevation.species column
seed.indices <- seed.indices %>%
  mutate(elevation.species = paste0(elevation, species))

# Create replicate summary for table
replicate_summary_2 <- seed.indices %>%
  mutate(treatment = recode(treatment,
                            "control" = "C",
                            "autogamy" = "A",
                            "geitonogamy" = "G",
                            "outcrossing" = "O")) %>%
  filter(treatment %in% c("A", "G", "O", "C")) %>%
  group_by(species, elevation, treatment) %>%
  summarise(n_replicates = n(), 
            mean_seedset = mean(seedset, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(mean_seedset = round(mean_seedset, 1)) %>%
  pivot_wider(
    names_from = treatment,
    values_from = c(n_replicates, mean_seedset),
    names_glue = "{treatment}_{.value}",
    values_fill = 0
  ) %>%
  arrange(species, elevation)
###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* 
###* TABLE 1
# Recode incorrect elevation values
replicate_summary_2 <- replicate_summary_2 %>%
  mutate(elevation = recode(elevation,
                            `3500` = "3400",
                            `4000` = "3800"))

replicate_summary_2$elevation <- factor(replicate_summary_2$elevation,
                                        levels = c("2300", "2800", "3400", "3800"))


# Format species names with italics and family names
rep_table <- replicate_summary_2 %>%
  mutate(Species = case_when(
    species == "Clematis s"  ~ "<i>Clematis simensis</i> (Ranunculaceae)",
    species == "Crepis h"    ~ "<i>Crepis hypochoeridea</i> (Asteraceae)",
    species == "Geranium a"  ~ "<i>Geranium arabicum</i> (Geraniaceae)",
    species == "Hypericum r" ~ "<i>Hypericum revolutum</i> (Hypericaceae)",
    species == "Lactuca i"   ~ "<i>Lactuca inermis</i> (Asteraceae)",
    species == "Senecio b"   ~ "<i>Senecio burtonii</i> (Asteraceae)",
    species == "Senecio p"   ~ "<i>Senecio purpureus</i> (Asteraceae)"
  )) %>%
  arrange(Species, elevation) %>%
  group_by(Species) %>%
  mutate(Species = ifelse(row_number() == 1, Species, "")) %>%
  ungroup() %>%
  select(
    Species,
    Elevation = elevation,
    `Autogamy n`          = A_n_replicates,
    `Autogamy mean`       = A_mean_seedset,
    `Geitonogamy n`       = G_n_replicates,
    `Geitonogamy mean`    = G_mean_seedset,
    `Outcrossing n`       = O_n_replicates,
    `Outcrossing mean`    = O_mean_seedset,
    `Control n`           = C_n_replicates,
    `Control mean`        = C_mean_seedset
  )

rep_table

readr::write_csv(rep_table, "tables/table1_with_seedset.csv")
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
###* TABLE 2
###* 
###* The table is easiest to create in Word with the proper information. The 
###* only in for which is necessary to generate from R is the the elpd and the
###* SE of it
###* 
###* For the reproductive indices, only the k fold is used. Full models are
###* however added to check
###* 
###* Natural seedset
c_brm_model_zinb_new   <- readRDS("brms_models/c_brm_model_zinb_new.rds")
c_brm_model_zinb_new_k <- readRDS("brms_models/c_brm_model_zinb_new_k.rds")
c_brm_null_zinb_new    <- readRDS("brms_models/c_brm_null_zinb_new.rds")
c_brm_null_zinb_new_k  <- readRDS("brms_models/c_brm_null_zinb_new_k.rds")

loo_compare(c_brm_model_zinb_new_k, c_brm_null_zinb_new_k)

###* Pollen limitation
pl_ord_model_new   <- readRDS("brms_models/pl_ord_model_new.rds")
pl_ord_model_new_k <- readRDS("brms_models/pl_ord_model_new_k.rds")
pl_ord_null_new    <- readRDS("brms_models/pl_ord_null_new.rds")
pl_ord_null_new_k  <- readRDS("brms_models/pl_ord_null_new_k.rds")

loo_compare(pl_ord_model_new_k, pl_ord_null_new_k)

###* Autogamy index
ao_brm_model_zinb_new   <- readRDS("brms_models/ao_brm_model_zinb_new.rds")
ao_brm_model_zinb_new_k <- readRDS("brms_models/ao_brm_model_zinb_new_k.rds")
ao_brm_null_zinb        <- readRDS("brms_models/ao_brm_null_zinb_new.rds")
ao_brm_null_zinb_k      <- readRDS("brms_models/ao_brm_null_zinb_new_k.rds")

loo_compare(ao_brm_model_zinb_new_k, ao_brm_null_zinb_new_k)

###* Geitonogamy
go_brm_model_zinb_new   <- readRDS("brms_models/go_brm_model_zinb_new.rds")
go_brm_model_zinb_new_k <- readRDS("brms_models/go_brm_model_zinb_new_k.rds")
go_brm_null_zinb        <- readRDS("brms_models/go_brm_null_zinb_new.rds")
go_brm_null_zinb_k      <- readRDS("brms_models/go_brm_null_zinb_new_k.rds")

loo_compare(go_brm_model_zinb_new_k, go_brm_null_zinb_new_k)

###* Morphospecies richness 
m_brm_model_nb <- readRDS("brms_models/m_brm_model_nb.rds")
m_brm_null_nb  <- readRDS("brms_models/m_brm_null_nb.rds")

loo(m_brm_model_nb, m_brm_null_nb)

###* Functional-group richness
f_brm_model_po <- readRDS("brms_models/f_brm_model_po.rds")
f_brm_null_po  <- readRDS("brms_models/f_brm_null_po.rds")

loo(f_brm_model_po, f_brm_null_po)

###* Visitation frequency
vis_ord_model <- readRDS("brms_models/vis_ord_model.rds")
vis_ord_null  <- readRDS("brms_models/vis_ord_null.rds")

loo(vis_ord_model, vis_ord_null)
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
###* TABLE 3
###* 
###* The table is easiest to create in Word with the proper information. The 
###* only in for which is necessary to generate from R is the the elpd and the
###* SE of it
###* 
###* Natural seedset
###* 
###* This index had two models tied - one with visitation and morphospecies, the
###* other with visitation and functional-group richness
kfold_bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt <- readRDS("brms_models/kfold_bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt <- readRDS("brms_models/kfold_bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt.rds")
kfold_bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt <- readRDS("brms_models/kfold_bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt.rds")

loo_compare(kfold_bayesian_seedset_mean_corrected_scaled_15_2_zinb_k_opt,
            kfold_bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt)

loo_compare(kfold_bayesian_seedset_mean_corrected_scaled_13_2_zinb_k_opt,
            kfold_bayesian_seedset_mean_corrected_scaled_null2_zinb_k_opt)

###* Pollen limitation
kfold_bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2 <- readRDS("brms_models/kfold_bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2.rds")
kfold_bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2 <- readRDS("brms_models/kfold_bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2.rds")

loo_compare(kfold_bayesian_pl_mean_corrected_scaled_null_zoib_k_opt_2,
            kfold_bayesian_pl_mean_corrected_scaled_3_zoib_k_opt_2)

###* Autogamy
kfold_bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt <- readRDS("brms_models/kfold_bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt.rds")
kfold_bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt <- readRDS("brms_models/kfold_bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt.rds")

loo_compare(kfold_bayesian_ao_mean_corrected_scaled_1_2_zib_k_opt,
            kfold_bayesian_ao_mean_corrected_scaled_null_2_zib_k_opt)

###* Geitonogamy
kfold_bayesian_go_mean_corrected_scaled_5_2_zoib_k_opt <- readRDS("brms_models/kfold_bayesian_go_mean_corrected_scaled_5_2_zoib_k_opt.rds")
kfold_bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt <- readRDS("brms_models/kfold_bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt.rds")

loo_compare(kfold_bayesian_go_mean_corrected_scaled_null_2_zoib_k_opt,
            kfold_bayesian_go_mean_corrected_scaled_5_2_zoib_k_opt)