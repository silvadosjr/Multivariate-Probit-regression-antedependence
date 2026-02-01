###############################################################################
# Multivariate Probit Regression with Structured Dependence
# Robust Simulation Study: Fit Competing Correlation Models
#
# Purpose
# -------
# This script fits one or more Bayesian multivariate probit regression models
# to simulated datasets generated under a chosen correlation structure
# (AR1 or Toeplitz). The goal is to evaluate robustness / model misspecification
# by fitting alternative correlation models to the same data, e.g.:
#
#   - "Lcorr": an unstructured correlation model (e.g., LKJ/Cholesky-based)
#   - "AR1"  : first-order autoregressive correlation model
#   - "Ind"  : conditional independence (no residual correlation)
#
# The set of fitted models is controlled by the character vector `fit_models`.
# For each scenario and replication, the script:
#   1) Loads a simulated dataset (.rds) from disk
#   2) Builds the design array X for treated/control groups
#   3) Fits the requested model(s) in Stan via rstan
#   4) Saves posterior summary tables (est$summary) to disk
#
# Inputs
# ------
# - Simulated datasets saved as .rds files, located at:
#     SimulationStudies/LatinSquareScenarios/<corr_type>/
#
#   Expected naming convention:
#     SimMultProbit_<corr_type>_<vT>_<n>_<corr_label>_<r>.rds
#
# - Stan model files expected under Programs/:
#     FitMultProbit_Lcorr.stan
#     FitMultProbit_AR1.stan
#     FitProbit_Ind.stan
#
# - Auxiliary functions (if any) expected in:
#     Programs/Aux_Functions.R
#
# Outputs
# -------
# For each fitted model, scenario index j, and replication r, a summary table
# is saved as an .rds file in:
#     SimulationStudies/ResultsSimRobus/
#
# Naming convention:
#   summaryMult_<Model>_R<r>_S<j>.rds
#
# Notes / Assumptions
# -------------------
# - This script uses 1 chain by default (nChains = 1). With 1 chain, Rhat may
#   be undefined or not meaningful; multiple chains are recommended if you
#   intend to use Rhat as a convergence diagnostic.
# - The initialization function `ini()` includes parameters for all models.
#   Stan ignores any extra parameters not declared in the model.
# - For Toeplitz-generated data, the AR1 fitting model is a deliberate
#   misspecification used for robustness assessment.
#
###############################################################################

# Packages -----------------------------------------------------------------

library(rstan)     # Stan interface for Bayesian inference
library(rstansim)  # Simulation utilities (kept for workflow consistency)

# Use all available CPU cores for Stan
options(mc.cores = parallel::detectCores())

# Project root -------------------------------------------------------------

setwd('~/GitHub/Multivariate-Probit-regression-antedependence/')
library(here)

# Correlation structure used to generate the datasets ----------------------
# Choose between "Toeplitz" and "AR1"
corr_type <- "Toeplitz"  # or "AR1"

# Model selection: which models to fit ------------------------------------
# Examples:
#   fit_models <- c("Lcorr")
#   fit_models <- c("AR1", "Ind")
#   fit_models <- c("Lcorr", "AR1", "Ind")
fit_models <- c("AR1", "Lcorr")

# Paths --------------------------------------------------------------------

programas_dir <- here("Programs")
pathDataSets  <- here("SimulationStudies", "LatinSquareScenarios", corr_type)
pathResults   <- here("SimulationStudies", "ResultsSimRobus")

# Load auxiliary functions (if used) --------------------------------------
source(file.path(programas_dir, "Aux_Functions.R"))

# Correlation ranges (for scenario labeling) ------------------------------
# [.10, .40] Weak; [.40, .70] Moderate; [.70, .90] Strong

# Scenario grid ------------------------------------------------------------
# Each row: (vT, n, rho_start, rho_end)
Scenario <- rbind(
  c(4,   50, .9, .7),
  c(8,   50, .4, .1),
  c(12,  50, .7, .4),
  c(4,  100, .7, .4),
  c(8,  100, .9, .7),
  c(12, 100, .4, .1),
  c(4,  200, .4, .1),
  c(8,  200, .7, .4),
  c(12, 200, .9, .7)
)

R <- 20  # number of replications per scenario

# Main loop over scenarios -------------------------------------------------
for (j in 1:nrow(Scenario)) {
  
  # Extract scenario parameters
  vT        <- Scenario[j, 1]  # number of repeated measures
  n         <- Scenario[j, 2]  # sample size
  rho_start <- Scenario[j, 3]  # correlation strength (or start of Toeplitz)
  rho_end   <- Scenario[j, 4]  # end of Toeplitz correlation sequence
  
  # Build time covariate (centered) ---------------------------------------
  t_age <- 1:vT - mean(1:vT)
  
  # Design matrices for treatment/control groups
  X_s  <- cbind(1, t_age, 1, t_age)  # treated
  X_ns <- cbind(1, t_age, 0, 0)      # control
  
  # Balanced allocation to treatment and control groups
  treat_counts <- rep(n / 2, 2)
  
  p <- 4  # number of regression coefficients
  
  # Build X as a 3D array with dimensions [subject, time, covariate]
  X_array <- array(NA, dim = c(length(t_age), p, n))
  X_array[, , 1:treat_counts[1]] <- replicate(treat_counts[1], X_ns, simplify = "array")
  X_array[, , (treat_counts[1] + 1):n] <- replicate(treat_counts[2], X_s, simplify = "array")
  X_array <- aperm(X_array, c(3, 1, 2))
  
  # True parameters used to simulate data (kept for reference) ------------
  beta_verd <- c(1, .5, .8, .6)
  
  # True correlation parameter(s) used during simulation ------------------
  rho_verd <- if (corr_type == "AR1") {
    rho_start
  } else if (corr_type == "Toeplitz") {
    round(seq(rho_start, to = rho_end, length.out = vT - 1), 2)
  } else {
    stop("corr_type must be either 'AR1' or 'Toeplitz'")
  }
  
  # Correlation strength label (used in dataset file naming) --------------
  corr_label <- if (rho_start == 0.4) {
    "Weak"
  } else if (rho_start == 0.7) {
    "Moderate"
  } else {
    "Strong"
  }
  
  # Parameters monitored in Stan outputs ----------------------------------
  params_Lcorr <- c("beta", "C")    # unstructured correlation parameters (example: Cholesky factor C)
  params_AR1   <- c("beta", "rho")  # AR1 correlation parameter
  params_Ind   <- c("beta")         # independence model has no correlation parameters
  
  # MCMC settings ----------------------------------------------------------
  nChains        <- 1
  burnInSteps    <- 1000
  thinSteps      <- 15
  numSavedSteps  <- 1000  # total saved draws across chains
  nIter <- ceiling(burnInSteps + (numSavedSteps * thinSteps) / nChains)
  
  # Generic initialization function ---------------------------------------
  # Includes parameters used by different models; unused entries are ignored by Stan.
  ini <- function() {
    list(
      beta = rep(.1, p),
      rho  = rep(.1, max(vT - 1, 1)),  # AR1 expects scalar; Toeplitz may expect a vector
      C    = diag(vT)                  # used if the Lcorr model defines a matrix parameter C
    )
  }
  
  # Replication loop -------------------------------------------------------
  for (r in 1:R) {
    
    # Load simulated dataset ----------------------------------------------
    rds_path <- file.path(
      pathDataSets,
      paste0("SimMultProbit_", corr_type, "_", vT, "_", n, "_", corr_label, "_", r, ".rds")
    )
    sim_Toeplitz <- readRDS(rds_path)
    
    # Data list passed to Stan --------------------------------------------
    data_list <- list(
      vT         = vT,
      p          = p,
      N          = n,
      Y          = sim_Toeplitz$Y,
      X          = X_array,
      sigma_beta = 1,
      sigma_rho  = 1
    )
    
    ##============================== Fit requested models ==============================##
    
    # --- Lcorr model (unstructured correlation) --------------------------
    if ("Lcorr" %in% fit_models) {
      
      begin <- Sys.time()
      print(begin)
      
      samp_Lcorr <- stan(
        data     = data_list,
        file     = file.path(programas_dir, "FitMultProbit_Lcorr.stan"),
        init     = ini,
        chains   = nChains,
        pars     = params_Lcorr,
        iter     = nIter,
        warmup   = burnInSteps,
        thin     = thinSteps,
        control  = list(adapt_delta = 0.8, max_treedepth = 10),
        save_dso = TRUE
      )
      
      end <- Sys.time()
      print(end - begin)
      
      est_Lcorr <- summary(samp_Lcorr)
      
      Lcorr_file <- file.path(
        pathResults,
        paste0("summaryMult_Lcorr_R", r, "_S", j, ".rds")
      )
      saveRDS(est_Lcorr$summary, Lcorr_file)
    }
    
    # --- AR1 model -------------------------------------------------------
    if ("AR1" %in% fit_models) {
      
      begin <- Sys.time()
      print(begin)
      
      samp_AR1 <- stan(
        data     = data_list,
        file     = file.path(programas_dir, "FitMultProbit_AR1.stan"),
        init     = ini,
        chains   = nChains,
        pars     = params_AR1,
        iter     = nIter,
        warmup   = burnInSteps,
        thin     = thinSteps,
        control  = list(adapt_delta = 0.8, max_treedepth = 10),
        save_dso = TRUE
      )
      
      end <- Sys.time()
      print(end - begin)
      
      est_AR1 <- summary(samp_AR1)
      
      AR1_file <- file.path(
        pathResults,
        paste0("summaryMult_AR1_R", r, "_S", j, ".rds")
      )
      saveRDS(est_AR1$summary, AR1_file)
    }
    
    # --- Independence model ----------------------------------------------
    if ("Ind" %in% fit_models) {
      
      begin <- Sys.time()
      print(begin)
      
      # Independence model: beta only (no correlation parameters)
      samp_Ind <- stan(
        data     = data_list,
        file     = file.path(programas_dir, "FitProbit_Ind.stan"),
        init     = ini,
        chains   = nChains,
        pars     = params_Ind,
        iter     = nIter,
        warmup   = burnInSteps,
        thin     = thinSteps,
        control  = list(adapt_delta = 0.8, max_treedepth = 10),
        save_dso = TRUE
      )
      
      end <- Sys.time()
      print(end - begin)
      
      est_Ind <- summary(samp_Ind)
      
      Ind_file <- file.path(
        pathResults,
        paste0("summaryMult_Ind_R", r, "_S", j, ".rds")
      )
      saveRDS(est_Ind$summary, Ind_file)
    }
    
    message("Replication ", r, " completed.")
  }
  
  message("Scenario ", j, " completed.")
}
