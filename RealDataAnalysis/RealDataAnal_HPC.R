# Packages -----------------------------------------------------------------
library(rstan)
library(dplyr)
library(tidyr)
library(here)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

setwd('~/GitHub/Multivariate-Probit-regression-antedependence/')
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


# Fraction of data held out for HPC
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


##================================= Fitting models ========================================##

# Escolha da estrutura de correlação
# 1 = Toeplitz
# 2 = AR(1)
# 3 = AD(1)
# 4 = ARMA(1,1)

data_train <- list(
  vT = vT,
  p  = p,
  N  = N_train,
  Y  = Y_train,
  X  = X_train,
  sigma_beta = 10,
  sigma_rho  = 1,
  cor_type   = 2   # example: AR(1)
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

