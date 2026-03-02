# Packages -----------------------------------------------------------------
library(rstan)
library(dplyr)
library(tidyr)
library(loo)


options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

setwd('~/GitHub/Multivariate-Probit-regression-antedependence/')
library(here)
source(here("Programs", "Aux_Functions.R"))

# Directory to store fitted models
pathFit <- here('RealDataAnalysis','Fit')
dir.create(pathFit, recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# Data preparation: Six Cities data
# ===========================================================================
six_cities <- read.csv(file.path('RealDataAnalysis','Data','six_cities_expanded.csv'))

# Number of repeated measures (time points)
vT <- 4

# Centered time covariate
t_age <- 1:vT - mean(1:vT)

# Design matrices for treated and non-treated groups
X_s  <- cbind(1, t_age, 1, t_age)
X_ns <- cbind(1, t_age, 0, 0)

# Number of subjects per group
treat_counts <- table(six_cities$maternal_smoking)
n <- sum(treat_counts)

# Number of regression coefficients
p <- 4

# Build array of design matrices: dimensions (N x vT x p)
X_array <- array(NA, dim = c(length(t_age), p, n))
X_array[, , 1:treat_counts[1]] <- replicate(treat_counts[1], X_ns, simplify = "array")
X_array[, , (treat_counts[1] + 1):n] <- replicate(treat_counts[2], X_s, simplify = "array")
X_array <- aperm(X_array, c(3, 1, 2))

# Binary response matrix (N x vT)
Y <- six_cities[, -5]

# ===========================================================================
# MCMC configuration
# ===========================================================================
nChains       <- 3
burnInSteps   <- 1000
thinSteps     <- 10
numSavedSteps <- 1000

nIter <- ceiling(burnInSteps + (numSavedSteps * thinSteps) / nChains)

control_list <- list(adapt_delta = 0.9, max_treedepth = 15)

# ===========================================================================
# Stan model files
# ===========================================================================
stan_files <- list(
  structured   = file.path('Programs', "FitMultProbit_Structured_Diag.stan"),
  unstructured = file.path('Programs', "FitMultProbit_Lcorr_Diag.stan"),
  independent  = file.path('Programs', "FitProbit_Ind_Diag.stan")
)

# Base data shared by all models
base_data <- list(
  vT = vT,
  p  = p,
  N  = n,
  Y  = Y,
  X  = X_array,
  sigma_beta = 10,
  sigma_rho  = 1
)

# ===========================================================================
# Helper functions
# ===========================================================================


# Data list by model
make_data_list <- function(model_name) {
  if (model_name %in% c("Toeplitz","AR1","AD","ARMA11")) {
    cor_type <- switch(model_name,
                       "Toeplitz" = 1L,
                       "AR1"      = 2L,
                       "AD"       = 3L,
                       "ARMA11"   = 4L)
    return(c(base_data, list(cor_type = cor_type)))
  }
  base_data
}

# Select Stan file by model
pick_stan_file <- function(model_name) {
  if (model_name %in% c("Toeplitz","AR1","AD","ARMA11")) return(stan_files$structured)
  if (model_name == "Unstructured") return(stan_files$unstructured)
  if (model_name == "Independent")  return(stan_files$independent)
  stop("Unknown model in pick_stan_file().")
}

# ===========================================================================
# Fit one model and compute LOO / WAIC
# ===========================================================================
fit_one_model <- function(model_name) {
  
  data_list <- make_data_list(model_name)
  pars      <- make_pars(model_name)
  init_fun  <- make_inits(model_name, p = p, vT = vT)
  stan_file <- pick_stan_file(model_name)
  
  cat("\n=================================\n")
  cat("Fitting model:", model_name, "\n")
  cat("Stan file:", stan_file, "\n")
  cat("=================================\n")
  
  begin <- Sys.time()
  fit <- stan(
    data   = data_list,
    file   = stan_file,
    init   = init_fun,
    chains = nChains,
    pars   = pars,
    iter   = nIter,
    warmup = burnInSteps,
    thin   = thinSteps,
    control = control_list,
    save_dso = TRUE
  )
  end <- Sys.time()
  cat("Elapsed time:", end - begin, "\n")
  
  # Save fitted object
  out_rds <- file.path(pathFit, paste0("ResultsFit_", model_name, "_Real.rds"))
  saveRDS(fit, out_rds)
  
  # Extract pointwise log-likelihood (S x N)
  log_lik <- loo::extract_log_lik(fit, parameter_name = "log_lik",
                                  merge_chains = TRUE)
  
  # Information criteria
  loo_res  <- loo::loo(log_lik)
  waic_res <- loo::waic(log_lik)
  
  list(
    model = model_name,
    fit   = fit,
    loo   = loo_res,
    waic  = waic_res,
    log_lik = log_lik,
    rds_path = out_rds
  )
}

# ===========================================================================
# Fit all models and compare
# ===========================================================================
models <- c("Independent", "Toeplitz", "AR1", "AD", "ARMA11", "Unstructured")

results <- lapply(models, fit_one_model)

# -----------------------------------------------------------------------
# Load fitted models and compute LOO / WAIC
# -----------------------------------------------------------------------
results <- lapply(models, function(m) {
  
  cat("Loading model:", m, "\n")
  
  fit <- readRDS(file.path(pathFit, paste0("ResultsFit_", m, "_Real.rds")))
  
  # Extract pointwise log-likelihood (S x N)
  log_lik <- loo::extract_log_lik(
    fit,
    parameter_name = "log_lik",
    merge_chains   = TRUE
  )
  
  list(
    model   = m,
    fit     = fit,
    log_lik = log_lik,
    loo     = loo::loo(log_lik),
    waic    = loo::waic(log_lik)
  )
})
names(results) <- models

# LOO comparison
loo_list <- lapply(results, `[[`, "loo")
loo_comp <- loo::loo_compare(loo_list)
print(loo_comp,simplify = F)

# WAIC comparison (less stable, but reported)
waic_list <- lapply(results, `[[`, "waic")
waic_comp <- loo::loo_compare(waic_list)
print(waic_comp)

# Summary table
summary_table <- tibble(
  model = models,
  elpd_loo  = sapply(loo_list,  function(x) x$estimates["elpd_loo","Estimate"]),
  se_elpd   = sapply(loo_list,  function(x) x$estimates["elpd_loo","SE"]),
  looic     = sapply(loo_list,  function(x) x$estimates["looic","Estimate"]),
  p_loo     = sapply(loo_list,  function(x) x$estimates["p_loo","Estimate"]),
  elpd_waic = sapply(waic_list, function(x) x$estimates["elpd_waic","Estimate"]),
  se_waic   = sapply(waic_list, function(x) x$estimates["elpd_waic","SE"]),
  waic      = sapply(waic_list, function(x) x$estimates["waic","Estimate"]),
  p_waic    = sapply(waic_list, function(x) x$estimates["p_waic","Estimate"]),
  rds_path  = sapply(results, `[[`, "rds_path")
) |>
  arrange(desc(elpd_loo))

print(summary_table)
