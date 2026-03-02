# Packages -----------------------------------------------------------------
library(rstan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

setwd('~/GitHub/Multivariate-Probit-regression-antedependence/')
library(here)
source(here("Programs", "Aux_Functions.R"))

# Directory to store fitted models
pathFit <- here('RealDataAnalysis','Fit')
dir.create(pathFit, recursive = TRUE, showWarnings = FALSE)

pathFig<-here('Figures')

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
Y <- as.matrix(six_cities[, -5])

# ===========================================================================
# Train/test split for HPC (subject-level split)
# ===========================================================================
test_frac <- 0.2

split <- make_train_test_split(
  Y = Y,
  X = X_array,
  test_frac = test_frac,
  seed = 258
)

# Training data
Y_train <- split$Y_train
X_train <- split$X_train
N_train <- dim(Y_train)[1]

# Test data
Y_test <- split$Y_test
X_test <- split$X_test
N_test <- dim(Y_test)[1]

# Sanity checks
dim(Y_train)  # N_train x vT
dim(Y_test)   # N_test  x vT
dim(X_train)  # N_train x vT x p
dim(X_test)   # N_test  x vT x p

# ===========================================================================
# Stan model files (as specified)
# ===========================================================================
stan_files <- list(
  structured   = file.path('Programs', "FitMultProbit_Structured.stan"),
  unstructured = file.path('Programs', "FitMultProbit_Lcorr.stan"),
  independent  = file.path('Programs', "FitProbit_Ind.stan")
)

# ===========================================================================
# Base training data (shared)
# ===========================================================================
base_train_data <- list(
  vT = vT,
  p  = p,
  N  = N_train,
  Y  = Y_train,
  X  = X_train,
  sigma_beta = 10,
  sigma_rho  = 1
)

# ===========================================================================
# MCMC configuration
# ===========================================================================
nChains       <- 1
burnInSteps   <- 1000
thinSteps     <- 10
numSavedSteps <- 1000

nIter <- ceiling(burnInSteps + (numSavedSteps * thinSteps) / nChains)

control_list <- list(adapt_delta = 0.9, max_treedepth = 15)

# ===========================================================================
# Model registry (what to fit, with which Stan file and (if needed) cor_type)
# ---------------------------------------------------------------------------
# Structured cor_type:
#   1 = Toeplitz
#   2 = AR(1)
#   3 = AD(1)
#   4 = ARMA(1,1)
# ===========================================================================
model_grid <- tibble::tibble(
  model_name = c("Independent", "Toeplitz", "AR1", "AD", "ARMA11", "Unstructured"),
  stan_key   = c("independent", "structured", "structured", "structured", "structured", "unstructured"),
  cor_type   = c(NA_integer_, 1L, 2L, 3L, 4L, NA_integer_)
)

# ===========================================================================
# Fit a single model on training data (optionally re-use saved fit)
# ===========================================================================
fit_one_model_train <- function(model_name, stan_key, cor_type,
                                base_train_data, pathFit,
                                nChains, nIter, burnInSteps, thinSteps,
                                control_list,
                                overwrite = FALSE,
                                log_lik = FALSE) {
  
  # Output path (consistent naming)
  out_rds <- file.path(pathFit, paste0("ResultsFit_", model_name, "_train.rds"))
  
  # Reuse if already saved (unless overwrite = TRUE)
  if (file.exists(out_rds) && !overwrite) {
    cat("Loading saved fit:", model_name, "\n")
    return(readRDS(out_rds))
  }
  
  # Build data list
  data_list <- base_train_data
  if (!is.na(cor_type)) data_list$cor_type <- as.integer(cor_type)
  
  # Select Stan file
  stan_file <- stan_files[[stan_key]]
  
  # Parameters to monitor (log_lik typically FALSE for HPC fitting)
  pars <- make_pars(model_name, log_lik = log_lik)
  
  # Initial values
  init_fun <- make_inits(model_name, p = base_train_data$p, vT = base_train_data$vT)
  
  cat("\n=================================\n")
  cat("Fitting model:", model_name, "\n")
  cat("Stan file:", stan_file, "\n")
  if (!is.na(cor_type)) cat("cor_type:", cor_type, "\n")
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
  
  saveRDS(fit, out_rds)
  fit
}

# ===========================================================================
# Run HPC (Q3 adjusted) for a fitted model
# ===========================================================================
hpc_one_model <- function(fit, model_name, cor_type,
                          X_test, Y_test, vT,
                          S = 400, seed = 258,
                          global_type = "max_abs") {
  
  hpc_q3_adjusted(
    fit = fit,
    model_name = model_name,
    X_test = X_test,
    Y_test = Y_test,
    vT = vT,
    S = S,
    seed = seed,
    cor_type = cor_type,
    global_type = global_type
  )
}

# ===========================================================================
# Fit all models (training) + run HPC in batch + summarize results
# ===========================================================================
overwrite_fits <- FALSE   # set TRUE to refit everything
S_hpc <- 400
seed_hpc <- 258
global_type <- "max_abs"

# Fit (or load) all models on training data
fits_train <- lapply(seq_len(nrow(model_grid)), function(i) {
  row <- model_grid[i, ]
  fit_one_model_train(
    model_name = row$model_name,
    stan_key   = row$stan_key,
    cor_type   = row$cor_type,
    base_train_data = base_train_data,
    pathFit = pathFit,
    nChains = nChains,
    nIter = nIter,
    burnInSteps = burnInSteps,
    thinSteps = thinSteps,
    control_list = control_list,
    overwrite = overwrite_fits,
    log_lik = FALSE
  )
})
names(fits_train) <- model_grid$model_name

# Run HPC for each model
hpc_results <- lapply(seq_len(nrow(model_grid)), function(i) {
  row <- model_grid[i, ]
  cat("Running HPC (Q3) for:", row$model_name, "\n")
  hpc_one_model(
    fit = fits_train[[row$model_name]],
    model_name = row$model_name,
    cor_type = row$cor_type,
    X_test = X_test,
    Y_test = Y_test,
    vT = vT,
    S = S_hpc,
    seed = seed_hpc,
    global_type = global_type
  )
})
names(hpc_results) <- model_grid$model_name

# ===========================================================================
# Summaries: global p-value + a couple of handy Q3 summaries
# ===========================================================================
hpc_summary <- tibble::tibble(
  model = model_grid$model_name,
  ppp_global = sapply(hpc_results, function(res) res$ppp_global),
  # Magnitude of the (posterior-mean) "observed" Q3 on the holdout
  Q3_new_max_abs  = sapply(hpc_results, function(res) max(abs(res$Q3_new_mean), na.rm = TRUE)),
  Q3_new_mean_abs = sapply(hpc_results, function(res) {
    off <- res$Q3_new_mean[upper.tri(res$Q3_new_mean)]
    mean(abs(off), na.rm = TRUE)
  })
) |>
  arrange(desc(ppp_global))

print(hpc_summary)

# Example: inspect pairwise adjusted p-values for one model
# hpc_results$Toeplitz$ppp_pair
# hpc_results$Toeplitz$Q3_new_mean



# Compare Independent vs AR(1)
out <- plot_Q3_discrepancy_independent(
  hpc_independent = hpc_results$Independent,
  hpc_model       = hpc_results$AR1,
  model_label     = "AR(1)",
  common_scale    = TRUE,
  exp_fun         = "median",
  show_labels     = TRUE
)

out$plot_independent
out$plot_model




panel_ar1 <- plot_Q3_discrepancy_panel(
  hpc_independent = hpc_results$Independent,
  hpc_model       = hpc_results$AR1,
  model_label     = "AR(1)",
  exp_fun         = "median",
  show_labels     = TRUE
)


panel_ar1



ggsave(
  file.path(pathFig, "Heatmap_ar1_SixCities.eps"),
  panel_ar1,
  device = cairo_ps,
  width  = 6,
  height = 5,
  units  = "in"
)








