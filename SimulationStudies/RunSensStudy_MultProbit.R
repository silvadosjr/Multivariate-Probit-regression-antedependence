###############################################################################
# Multivariate Probit Regression with Structured Dependence
# Sensitivity Study Script (Priors)
#
# Purpose
# -------
# This script fits a Bayesian multivariate probit regression model to simulated
# datasets under different prior settings. It loops over:
#   (i) correlation scenarios (vT, n, rho strength),
#   (ii) Monte Carlo replications (R),
#   (iii) prior configurations for (beta, rho).
#
# The simulated datasets are assumed to have been generated previously and
# saved as .rds files. For each dataset, we run Stan (rstan) to obtain posterior
# summaries, then save the summary table to disk.
#
# Correlation structures supported by the simulation design:
#   - AR1      : single correlation parameter rho
#   - Toeplitz : vector of lag correlations (not used in this fitting script;
#               this script calls FitMultProbit_AR1.stan)
#
# Outputs are written as one summary .rds per replication/scenario/prior.
#
###############################################################################

# Packages -----------------------------------------------------------------

library(rstan)     # Stan interface for Bayesian inference
library(rstansim)  # Simulation utilities (not directly used below, but kept for workflow consistency)

# Project root -------------------------------------------------------------

setwd('~/GitHub/Multivariate-Probit-regression-antedependence/')
library(here)  # Robust path construction relative to the project root

# Use all available CPU cores for Stan
options(mc.cores = parallel::detectCores())

# Correlation structure ----------------------------------------------------
# Choose between "AR1" and "Toeplitz".
# NOTE: This fitting script calls FitMultProbit_AR1.stan (AR1 fitting model).
corr_type <- "AR1"  # or "Toeplitz"

# Correlation ranges (used in the scenario grid)
# [.40, .70] Moderate; [.70, .90] Strong; or Null (0)

# Paths -------------------------------------------------------------------
# Organized folder structure using `here()`.
programas_dir <- here("Programs")

# Folder containing the previously simulated datasets (RDS files)
pathDataSets <- here("SimulationStudies", "LatinSquareScenarios", corr_type)

# Folder where posterior summaries (per replication) will be saved
pathResults <- here("SimulationStudies", "SensStudy")

# Load auxiliary functions -------------------------------------------------
# Assumes helper utilities (if any) are defined in Programs/Aux_Functions.R.
source(file.path(programas_dir, "Aux_Functions.R"))

# Prior settings -----------------------------------------------------------
# Each row defines a prior configuration indexed by Prior.
# sigma_beta : prior scale for beta coefficients
# sigma_rho  : prior scale for rho
PriorSettings <- data.frame(
  Prior = paste0("Prior", 1:4),
  sigma_beta = c(1, 7, 10, 1),
  sigma_rho  = c(1, 1, 10, .5)
)

# Simulation scenarios -----------------------------------------------------
# Each row defines a scenario:
#   vT        : number of repeated measures (time points)
#   n         : sample size
#   rho_start : correlation strength (AR1 rho, or Toeplitz initial rho)
#   rho_end   : Toeplitz end rho (ignored for AR1; kept for naming consistency)
Scenario <- rbind(
  c(4,   50, .9, .7), c(8,   50, 0,  0), c(12,  50, .7, .4),
  c(4,  100, .7, .4), c(8,  100, .9, .7), c(12, 100, 0,  0),
  c(4,  200, 0,  0),  c(8,  200, .7, .4), c(12, 200, .9, .7)
)

R <- 20  # number of Monte Carlo replications per scenario

# Main loops ---------------------------------------------------------------
# k : index of prior setting (here, Prior2, Prior3, Prior4)
# j : index of simulation scenario
# r : replication index within scenario
for (k in c(2, 3, 4)) {
  
  for (j in seq_len(nrow(Scenario))) {
    
    # Extract scenario parameters
    vT        <- Scenario[j, 1]  # number of time points
    n         <- Scenario[j, 2]  # sample size
    rho_start <- Scenario[j, 3]  # correlation level for naming / AR1 rho
    rho_end   <- Scenario[j, 4]  # Toeplitz end rho (if applicable)
    
    # Build centered time covariate (helps interpretability and stability)
    t_age <- 1:vT - mean(1:vT)
    
    # Design matrices:
    # - X_ns: control group (no treatment effect)
    # - X_s : treated group (adds treatment intercept and slope)
    X_s  <- cbind(1, t_age, 1, t_age)
    X_ns <- cbind(1, t_age, 0, 0)
    
    # Balanced allocation to control/treatment groups
    treat_counts <- rep(n / 2, 2)
    p <- 4  # number of regression coefficients
    
    # Build 3D array X with dimensions [subject, time, covariate]
    X_array <- array(NA, dim = c(length(t_age), p, n))
    X_array[, , 1:treat_counts[1]] <- replicate(treat_counts[1], X_ns, simplify = "array")
    X_array[, , (treat_counts[1] + 1):n] <- replicate(treat_counts[2], X_s, simplify = "array")
    X_array <- aperm(X_array, c(3, 1, 2))
    
    # True regression coefficients used in simulation (kept here for reference)
    beta_verd <- c(1, .5, .8, .6)
    
    # True correlation parameter(s) used in simulation (kept here for reference)
    # AR1: scalar rho; Toeplitz: vector of lag correlations
    rho_verd <- if (corr_type == "AR1") {
      rho_start
    } else if (corr_type == "Toeplitz") {
      round(seq(rho_start, to = rho_end, length.out = vT - 1), 2)
    } else {
      stop("corr_type must be either 'AR1' or 'Toeplitz'")
    }
    
    # Label used for file naming / scenario grouping
    corr_label <- if (rho_start == 0) {
      "Null"
    } else if (rho_start == 0.7) {
      "Moderate"
    } else {
      "Strong"
    }
    
    # Replication loop -----------------------------------------------------
    for (r in seq_len(R)) {
      
      # Load the r-th simulated dataset for this scenario
      rds_path <- file.path(
        pathDataSets,
        paste0("SimMultProbit_", corr_type, "_", vT, "_", n, "_", corr_label, "_", r, ".rds")
      )
      sim_AR1 <- readRDS(rds_path)
      
      # Prepare the data list for Stan
      # NOTE:
      #  - N is used here as the sample size in Stan
      #  - sigma_beta and sigma_rho correspond to the chosen prior setting (k)
      data_list <- list(
        vT = vT,
        p  = p,
        N  = n,
        Y  = sim_AR1$Y,
        X  = X_array,
        sigma_beta = PriorSettings$sigma_beta[k],
        sigma_rho  = PriorSettings$sigma_rho[k]
      )
      
      # Initial values for Stan parameters
      # For AR1, rho is scalar; for Toeplitz models, this would need to be a vector.
      ini <- function() list(beta = rep(0.1, p), rho = 0.1)
      
      # Parameters to monitor
      params_AR1 <- c("beta", "rho")
      
      # MCMC settings
      nChains       <- 1
      burnInSteps   <- 1000
      thinSteps     <- 15
      numSavedSteps <- 1000
      
      # Total iterations: warmup + thinning-adjusted draws
      nIter <- ceiling(burnInSteps + (numSavedSteps * thinSteps) / nChains)
      
      # Fit model in Stan ---------------------------------------------------
      begin <- Sys.time()
      print(begin)
      
      samp_AR1 <- stan(
        data   = data_list,
        file   = file.path(programas_dir, "FitMultProbit_AR1.stan"),
        init   = ini,
        chains = nChains,
        pars   = params_AR1,
        iter   = nIter,
        warmup = burnInSteps,
        thin   = thinSteps,
        control = list(adapt_delta = 0.8, max_treedepth = 10),
        save_dso = TRUE
      )
      
      end <- Sys.time()
      print(end - begin)
      
      # Extract posterior summaries (mean, sd, quantiles, Rhat, n_eff, etc.)
      est_AR1 <- summary(samp_AR1)
      
      # Save summary table for downstream aggregation
      out_file <- file.path(pathResults, paste0("summaryMult_AR1_R", r, "_S", j, "_P", k, ".rds"))
      saveRDS(est_AR1$summary, out_file)
      
      message("Replication ", r, " completed.")
    }
    
    message("Scenario ", j, " completed.")
  }
  
  message("Prior ", k, " completed.")
}
