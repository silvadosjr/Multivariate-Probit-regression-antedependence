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
pathFit <- here('RealDataAnalysis','Madras','Fit')
dir.create(pathFit, recursive = TRUE, showWarnings = FALSE)


load(file.path('RealDataAnalysis','Data','madras.rda'))

# Convert long data to wide format (one row per subject)
madras_wide <- madras %>%
  arrange(id, month) %>%
  pivot_wider(
    id_cols = c(id, age, gender),
    names_from = month,
    values_from = y,
    names_prefix = "month_"
  )

# Remove subjects with any missing repeated measurements
#madras_wide <- na.omit(madras_wide)


## Design matrix for: y ~ month + age + gender + monthXage + monthXgender
## (Intercept is included automatically)

# Number of repeated measures (time points)
vT <- 12

# Centered time covariate (improves numerical stability and interpretation)
t_month <- 0:(vT-1) - mean(0:(vT-1))

# Group-specific design blocks for (age, gender) ∈ {0,1} × {0,1}
# Columns: (Intercept, month, age, gender, month*age, month*gender)

X_00 <- cbind(1, t_month, 0, 0, 0,        0)        # age=0, gender=0
X_01 <- cbind(1, t_month, 0, 1, 0,        t_month)  # age=0, gender=1
X_10 <- cbind(1, t_month, 1, 0, t_month,  0)        # age=1, gender=0
X_11 <- cbind(1, t_month, 1, 1, t_month,  t_month)  # age=1, gender=1

# Count number of subjects in each (age, gender) group
# Order enforced as: (0,0), (0,1), (1,0), (1,1)
tab <- table(madras_wide$age, madras_wide$gender)
treat_counts <- c(tab["0","0"], tab["0","1"], tab["1","0"], tab["1","1"])
n <- sum(treat_counts)

# Number of regression coefficients
p <- ncol(X_00)

# Build design array with dimensions (vT x p x N)
X_array <- array(NA, dim = c(vT, p, n))

# Fill array by contiguous group blocks
i1 <- 1
i2 <- treat_counts[1]
X_array[,, i1:i2] <- replicate(treat_counts[1], X_00, simplify = "array")

i1 <- i2 + 1
i2 <- i2 + treat_counts[2]
X_array[,, i1:i2] <- replicate(treat_counts[2], X_01, simplify = "array")

i1 <- i2 + 1
i2 <- i2 + treat_counts[3]
X_array[,, i1:i2] <- replicate(treat_counts[3], X_10, simplify = "array")

i1 <- i2 + 1
i2 <- i2 + treat_counts[4]
X_array[,, i1:i2] <- replicate(treat_counts[4], X_11, simplify = "array")

# Rearrange dimensions to (N x vT x p), matching typical Stan structure
X_array <- aperm(X_array, c(3, 1, 2))

# Reorder madras_wide to match the block structure used in X_array
madras_wide <- madras_wide %>%
  mutate(age = as.integer(age), gender = as.integer(gender)) %>%
  arrange(age, gender, id)

# Binary response matrix (N x vT)
# Removes id, age, and gender columns
Y <- as.matrix(madras_wide[,-c(1,2,3)])


t_index <- t(apply(Y, 1, function(x) {
  idx <- which(!is.na(x))
  c(idx, rep(0L, ncol(Y) - length(idx)))
}))

n_obs <- rowSums(t_index > 0)


Y[is.na(Y)]<- -1L

storage.mode(Y) <- "integer"



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
  structured   = file.path('Programs', "FitMultProbit_Structured_Dropout_Diag.stan"),
  unstructured = file.path('Programs', "FitMultProbit_Lcorr_Missing_Diag.stan"),
  independent  = file.path('Programs', "FitProbit_Ind_Missing_Diag.stan")
)

# Base data shared by all models
base_data <- list(
  vT = vT,
  p  = p,
  N  = n,
  Y  = Y,
  X  = X_array,
  sigma_beta = 10,
  sigma_rho  = 1,
  n_obs=n_obs,
  t_index=t_index
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
  pars      <- make_pars(model_name,log_lik = T)
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
models <- c("Independent","Unstructured",'AD')

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
print(loo_comp)

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











