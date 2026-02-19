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



