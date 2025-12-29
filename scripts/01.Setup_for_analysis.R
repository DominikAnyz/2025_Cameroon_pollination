# ------------------------------------------------------------
# SETUP: run this first
# ------------------------------------------------------------

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# pacman first (so we can p_load the rest)
install_if_missing("pacman")
library(pacman)

# Packages actually used across scripts 01–11
pacman::p_load(tidyverse, here, readxl, brms, loo, emmeans, ordbetareg, corrplot,
  patchwork, scales, gt, sf, ggspatial, elevatr, raster, rnaturalearth
)

# CmdStan setup (needed for brms backend)
install_if_missing("cmdstanr")
library(cmdstanr)

if (is.null(cmdstanr::cmdstan_version())) {
  cmdstanr::install_cmdstan(cores = 4)
}

# Recommended global option (only if you always want this)
options(brms.backend = "cmdstanr")
