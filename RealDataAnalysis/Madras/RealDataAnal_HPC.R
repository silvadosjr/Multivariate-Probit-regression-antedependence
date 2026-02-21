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

Y[is.na(Y)]<- -1L

storage.mode(Y) <- "integer"


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

t_index_train <- t(apply(Y_train, 1, function(x) {
  idx <- which(!is.na(x))
  c(idx, rep(0L, ncol(Y_train) - length(idx)))
}))

n_obs_train <- rowSums(t_index_train > 0)


# Test data
Y_test <- split$Y_test
X_test <- split$X_test
N_test <- dim(Y_test)[1]

t_index_test <- t(apply(Y_test, 1, function(x) {
  idx <- which(!is.na(x))
  c(idx, rep(0L, ncol(Y_test) - length(idx)))
}))

n_obs_test <- rowSums(t_index_test > 0)



# Sanity checks
dim(Y_train)  # N_train x vT
dim(Y_test)   # N_test  x vT
dim(X_train)  # N_train x vT x p
dim(X_test)   # N_test  x vT x p

# ===========================================================================
# Stan model files (as specified)
# ===========================================================================
stan_files <- list(
  structured   = file.path('Programs', "FitMultProbit_Structured_Dropout.stan"),
  unstructured = file.path('Programs', "FitMultProbit_Lcorr_Missing.stan"),
  independent  = file.path('Programs', "FitProbit_Ind_Missing.stan")
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
  sigma_rho  = 1,
  n_obs=n_obs_train,
  t_index=t_index_train
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
  model_name = c("Independent", "Unstructured",'AD'),
  stan_key   = c("independent","unstructured","structured"),
  cor_type   = c(NA_integer_,NA_integer_,3L)
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



# Compare Independent vs AD
out <- plot_Q3_discrepancy_independent(
  hpc_independent = hpc_results$Independent,
  hpc_model       = hpc_results$AD,
  model_label     = "AD",
  common_scale    = TRUE,
  exp_fun         = "median",
  show_labels     = TRUE
)

out$plot_independent
out$plot_model




panel_AD <- plot_Q3_discrepancy_panel(
  hpc_independent = hpc_results$Independent,
  hpc_model       = hpc_results$AD,
  model_label     = "AD",
  exp_fun         = "median",
  show_labels     = TRUE
)


panel_AD










