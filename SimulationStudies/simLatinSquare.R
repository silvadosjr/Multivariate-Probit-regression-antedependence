###############################################################################
# Multivariate Probit Regression with Structured Dependence
# Simulation Script
#
# Purpose
# -------
# This script generates simulated datasets for a multivariate probit
# regression model with repeated measures and structured within-subject
# dependence. The correlation structure can be specified as:
#   - AR(1): first-order autoregressive correlation
#   - Toeplitz: lag-dependent correlation with decreasing strength
#
# The script is used to create Monte Carlo simulation scenarios for the paper:
# "Multivariate Probit Regression Using Antedependence Structures".
#
#
# Main Features
# -------------
# - Supports multiple scenario configurations (vT, n, correlation strength)
# - Balanced treatment/control design
# - Flexible correlation structure (AR1 or Toeplitz)
# - Uses Stan via rstansim for data generation
# - Automatically organizes outputs by correlation type and scenario
#
#
# Inputs
# ------
# - corr_type : character
#       Type of correlation structure. Must be either:
#       "AR1" or "Toeplitz"
#
# - Scenario : matrix
#       Each row defines a simulation scenario with:
#       (1) vT        = number of repeated measures
#       (2) n         = sample size
#       (3) rho_start = correlation level (or AR1 rho)
#       (4) rho_end   = final correlation level (Toeplitz only)
#
# - R : integer
#       Number of Monte Carlo replications per scenario
#
#
# Outputs
# -------
# For each scenario and replication, the script saves:
# - Simulated binary responses Y
#
# Files are stored in:
#   LatinSquareScenarios/<corr_type>/
#
# with filenames encoding:
#   - correlation type
#   - number of time points (vT)
#   - sample size (n)
#   - correlation strength label (Null / Moderate / Strong)
#
#
# Dependencies
# ------------
# R packages:
#   - rstan
#   - rstansim
#   - here
#
# Stan models (expected in Programs/):
#   - simMultProbit_AR1.stan
#   - simMultProbit_HT.stan
#
# Auxiliary R functions:
#   - Programs/Aux_Functions.R
#
#
# How to Run
# ----------
# 1. Clone the repository:
#      git clone <repository-url>
#
# 2. Set the working directory to the project root.
#
# 3. Choose the correlation type by setting:
#      corr_type <- "AR1"  or  "Toeplitz"
#
# 4. Run the script. Simulation outputs will be saved automatically.
#
#
# Expected Project Structure
# --------------------------
# Multivariate-Probit-regression-antedependence/
# ├── Programs/
# │   ├── Aux_Functions.R
# │   ├── simMultProbit_AR1.stan
# │   └── simMultProbit_HT.stan
# ├── LatinSquareScenarios/
# │   ├── AR1/
# │   └── Toeplitz/
# └── (analysis / figures / paper files)
#
#
# Notes
# -----
# - The design matrix assumes equal allocation to treatment and control.
# - Time covariates are centered to improve interpretability and stability.
# - Correlation labels (Null / Moderate / Strong) are used for file naming
#   and post-simulation summaries.
#
###############################################################################



# Load packages ------------------------------------------------------------

library(rstan)     # Interface to Stan (compiling and fitting Stan models)
library(rstansim)  # Helper utilities for simulation workflows with Stan

# Set project root directory (GitHub repo) --------------------------------

setwd('~/GitHub/Multivariate-Probit-regression-antedependence/')

library(here)      # Build robust file paths relative to the project root
options(mc.cores = parallel::detectCores())  # Use all available CPU cores for Stan

# Choose the within-subject correlation structure --------------------------
# Options currently supported in this script: "AR1" or "Toeplitz"
corr_type <- "AR1"

# Correlation ranges used in the simulation design
# [.10, .40] Weak; [.40, .70] Moderate; [.70, .90] Strong; or Null (0)

# Define paths for scripts and simulation outputs --------------------------

programas_dir <- here("Programs")                    # Folder with Stan/R helper scripts
pathSim <- here("LatinSquareScenarios", corr_type)   # Output folder depends on corr_type

# Load auxiliary functions (e.g., simulate_data()) -------------------------

source(file.path(programas_dir, "Aux_Functions.R"))

# Define simulation scenarios ---------------------------------------------
# Each row corresponds to one scenario:
# (1) vT        = number of time points / repeated measures
# (2) n         = sample size
# (3) rho_start = starting correlation level (or AR1 correlation)
# (4) rho_end   = ending correlation level (used only for Toeplitz)

Scenario <- rbind(
  c(4,  50, .9, .7), c(8,  50, 0,  0), c(12,  50, .7, .4),
  c(4, 100, .7, .4), c(8, 100, .9, .7), c(12, 100, 0,  0),
  c(4, 200, 0,  0),  c(8, 200, .7, .4), c(12, 200, .9, .7)
)

R <- 20  # Number of Monte Carlo replications per scenario (nsim in rstansim)

# Loop over scenarios ------------------------------------------------------

for (j in seq_len(nrow(Scenario))) {
  
  # Extract scenario settings
  vT <- Scenario[j, 1]         # number of repeated measures
  n  <- Scenario[j, 2]         # sample size
  rho_start <- Scenario[j, 3]  # start correlation (or AR1 rho)
  rho_end   <- Scenario[j, 4]  # end correlation (Toeplitz only)
  
  # Build the design matrices ---------------------------------------------
  # Time index centered around zero (helps interpret intercept/slope)
  t_age <- 1:vT - mean(1:vT)
  
  # X_s  : treated group design matrix (includes treatment effects)
  # X_ns : control group design matrix (treatment terms set to zero)
  #
  # Model form implied here:
  #   eta_it = beta1 + beta2 * t_age + (treat) * (beta3 + beta4 * t_age)
  #
  X_s  <- cbind(1, t_age, 1, t_age)
  X_ns <- cbind(1, t_age, 0, 0)
  
  # Split sample equally into control/treatment groups
  treat_counts <- c(n / 2, n / 2)
  
  p <- 4  # number of regression coefficients
  
  # Build a 3D array of covariates X_array with dimensions:
  #   [subject, time, covariate]
  X_array <- array(NA, dim = c(length(t_age), p, sum(treat_counts)))
  
  # Fill subjects 1:(n/2) with control design (X_ns)
  X_array[, , 1:treat_counts[1]] <- replicate(treat_counts[1], X_ns, simplify = "array")
  
  # Fill subjects (n/2+1):n with treated design (X_s)
  X_array[, , (treat_counts[1] + 1):sum(treat_counts)] <- replicate(treat_counts[2], X_s, simplify = "array")
  
  # Reorder dimensions to match typical Stan expectation:
  #   X[subject, time, covariate]
  X_array <- aperm(X_array, c(3, 1, 2))
  
  # True regression parameters used to simulate data ----------------------
  
  beta_verd <- c(1, .5, .8, .6)
  
  # Define the true correlation parameter(s) ------------------------------
  # AR1: single rho
  # Toeplitz: vector of length (vT - 1) defining lag-specific correlations
  if (corr_type == "AR1") {
    rho_verd <- rho_start
  } else if (corr_type == "Toeplitz") {
    rho_verd <- round(seq(rho_start, to = rho_end, length.out = vT - 1), 2)
  } else {
    stop("corr_type must be either 'AR1' or 'Toeplitz'")
  }
  
  # Label the scenario by correlation strength for file naming ------------
  corr_label <- if (rho_start == 0) {
    "Null"
  } else if (rho_start == 0.7) {
    "Moderate"
  } else {
    "Strong"
  }
  
  # Prepare inputs for the Stan-based simulator ---------------------------
  # data_simAux  : list passed as 'data' to the Stan program
  # param_sim1   : list of true parameters used by rstansim for simulation
  data_simAux <- list(vT = vT, p = length(beta_verd), n = n, X = X_array)
  param_sim1  <- list(beta = beta_verd, rho = rho_verd)
  
  # Run simulation --------------------------------------------------------
  # simulate_data() is assumed to be defined in Aux_Functions.R.
  # It compiles the Stan file (if needed), simulates nsim datasets,
  # and saves the requested variables (here only "Y") into pathSim.
  
  begin <- Sys.time()
  print(begin)
  
  sim_out <- simulate_data(
    file = paste0(
      programas_dir,
      ifelse(corr_type == "AR1", "/simMultProbit_AR1.stan", "/simMultProbit_HT.stan")
    ),
    data_name = paste0("SimMultProbit_", corr_type, "_", vT, "_", n, "_", corr_label),
    input_data = data_simAux,
    param_values = param_sim1,
    vars = c("Y"),     # Only store the simulated binary outcomes
    nsim = R,          # Number of replications
    path = pathSim     # Output folder
  )
  
  end <- Sys.time()
  print(end - begin)   # Report elapsed time for this scenario
}








