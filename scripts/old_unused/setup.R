# List of required packages
packages <- c("tidyverse", "glmmTMB", "DHARMa", "brms", "ordbetareg", "bayesplot")

# Install missing packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Install all required packages
lapply(packages, install_if_missing)

# Load the packages
lapply(packages, library, character.only = TRUE)

# Ensure RStan and cmdstanr are set up correctly for brms
if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
}

library(cmdstanr)

# Install CmdStan (if not installed)
if (is.null(cmdstan_version())) {
  cmdstanr::install_cmdstan(cores = 4)  # Adjust cores if needed
}

#* Check if CmdStan is installed
cmdstanr::cmdstan_path()

cmdstanr::install_cmdstan(cores = 4)

cmdstanr::cmdstan_version()

library(brms)



