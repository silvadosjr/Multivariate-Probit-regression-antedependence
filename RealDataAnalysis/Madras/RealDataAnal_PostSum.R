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

t_index <- t(apply(Y, 1, function(x) {
  idx <- which(!is.na(x))
  c(idx, rep(0L, ncol(Y) - length(idx)))
}))

n_obs <- rowSums(t_index > 0)

Y[is.na(Y)]<- -1L

storage.mode(Y) <- "integer"



##================================= Fitting models ========================================##

# Escolha da estrutura de correlação
# 1 = Toeplitz
# 2 = AR(1)
# 3 = AD(1)
# 4 = ARMA(1,1)

data_list=list(vT=vT,p=p,N=n,Y=Y,X=X_array,sigma_beta=10,sigma_rho=1,cor_type=3,n_obs=n_obs,t_index=t_index) 


model_name<-'AD'

params <- make_pars(model_name,log_lik = F)
nChains = 1
burnInSteps = 1000
thinSteps = 5
numSavedSteps = 1000  #across all chains
nIter = ceiling(burnInSteps + (numSavedSteps * thinSteps)/nChains)
nIter

#nIter=11000


ini<-make_inits(model_name, p = p, vT = vT)

begin<-Sys.time()
print(begin)


samp <- stan(data = data_list, file =file.path('Programs', "FitMultProbit_Structured_Dropout.stan"),init = ini,
              chains = nChains, pars = params, iter = nIter,
              warmup = burnInSteps, thin = thinSteps,control = list(adapt_delta = 0.9,max_treedepth=15),
              save_dso = T)


end<-Sys.time()
print(end-begin)


saveRDS(samp,paste0(pathFit,'/Samples_',model_name,'_madras.rds'))


samp<-readRDS(paste0(pathFit,'/Samples_',model_name,'_madras.rds'))


#launch_shinystan(samp)


fit <- rstan::extract(samp,permuted=T,inc_warmup=F)

#-----------------------------#
# Helper: posterior summaries #
#-----------------------------#
post_summ <- function(draws, par_names = NULL, prob = c(0.025, 0.5, 0.975)) {
  # draws: matrix (S x K) OR vector (S)
  if (is.vector(draws)) draws <- matrix(draws, ncol = 1)
  
  qs <- apply(draws, 2, quantile, probs = prob)
  out <- data.frame(
    mean   = colMeans(draws),
    sd     = apply(draws, 2, sd),
    q2.5   = qs[1, ],
    median = qs[2, ],
    q97.5  = qs[3, ],
    pr_gt0 = colMeans(draws > 0)
  )
  
  if (!is.null(par_names)) rownames(out) <- par_names
  out
}

#-----------------------------#
# 1) Betas                    #
#-----------------------------#
# Expected structure: fit$beta is (S x p)
stopifnot(!is.null(fit$beta))
S <- nrow(fit$beta)
p <- ncol(fit$beta)

beta_names <- c("(Intercept)", "month", "age", "gender", "monthXage", "monthXgender")
if (length(beta_names) != p) beta_names <- paste0("beta[", 1:p, "]")

beta_sum <- post_summ(fit$beta, par_names = beta_names)
beta_sum

#-----------------------------#
# 2) Rhos (AD)                #
#-----------------------------#
# Common naming in your Stan workflow: rho is (S x (vT-1)) for AD(1),
# where rho[t] is the lag-1 partial correlation parameter at transition t -> t+1.
stopifnot(!is.null(fit$rho))
K <- ncol(fit$rho)

rho_names <- paste0("rho[", 1:K, "]")   # or use rho_0,...,rho_{K-1} if you prefer
rho_sum <- post_summ(fit$rho, par_names = rho_names)
rho_sum

#-----------------------------#
# 3) Combine in one table      #
#-----------------------------#
beta_tbl <- tibble::rownames_to_column(beta_sum, var = "parameter") %>%
  mutate(block = "beta")

rho_tbl  <- tibble::rownames_to_column(rho_sum,  var = "parameter") %>%
  mutate(block = "rho")

post_table <- bind_rows(beta_tbl, rho_tbl) %>%
  select(block, parameter, mean, sd, median, q2.5, q97.5, pr_gt0)

post_table


# rho[1] corresponds to Corr(Z_0, Z_1) in the AD(1) lag-1 band
rho1  <- fit$rho[, 1]
rho_rest <- fit$rho[, 2:K, drop = FALSE]

delta <- rho1 - rowMeans(rho_rest)

c(
  mean_delta   = mean(delta),
  median_delta = median(delta),
  q2.5         = quantile(delta, 0.025),
  q97.5        = quantile(delta, 0.975),
  pr_delta_gt0 = mean(delta > 0)
)


#-----------------------------#
# Helper: posterior summaries #
#-----------------------------#
post_ci_df <- function(draws, par_names = NULL, level = 0.95, point = c("mean","median")) {
  point <- match.arg(point)
  alpha <- (1 - level) / 2
  probs <- c(alpha, 0.5, 1 - alpha)
  
  if (is.vector(draws)) draws <- matrix(draws, ncol = 1)
  
  qs <- apply(draws, 2, quantile, probs = probs)
  est <- if (point == "mean") colMeans(draws) else qs[2, ]
  
  out <- tibble::tibble(
    parameter = if (!is.null(par_names)) par_names else paste0("par[", seq_len(ncol(draws)), "]"),
    estimate  = as.numeric(est),
    lower     = as.numeric(qs[1, ]),
    upper     = as.numeric(qs[3, ])
  )
  out
}

#-----------------------------#
# 1) Beta plot                #
#-----------------------------#
beta_names <- c("(Intercept)", "month", "age", "gender", "monthXage", "monthXgender")
if (is.null(fit$beta)) stop("fit$beta not found.")
if (ncol(fit$beta) != length(beta_names)) beta_names <- paste0("beta[", 1:ncol(fit$beta), "]")

df_beta <- post_ci_df(fit$beta, par_names = beta_names, level = 0.95, point = "mean") %>%
  mutate(parameter = factor(parameter, levels = rev(parameter)))  # reverse for top-to-bottom order

p_beta <- ggplot(df_beta, aes(x = parameter, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  geom_point(size = 2) +
  coord_flip() +
  labs(x = NULL, y = "Estimate", title = "Posterior estimates (mean) and 95% CrI: Betas") +
  theme_minimal(base_size = 12)

p_beta


#-----------------------------#
# 2) Rho plot                 #
#-----------------------------#
if (is.null(fit$rho)) stop("fit$rho not found.")
rho_names <- paste0("rho[", 1:ncol(fit$rho), "]")

df_rho <- post_ci_df(fit$rho, par_names = rho_names, level = 0.95, point = "mean") %>%
  mutate(parameter = factor(parameter, levels = rev(parameter)))

p_rho <- ggplot(df_rho, aes(x = parameter, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  geom_point(size = 2) +
  coord_flip() +
  labs(x = NULL, y = "Estimate", title = "Posterior estimates (mean) and 95% CrI: Rhos (AD)") +
  theme_minimal(base_size = 12)

p_rho



# fit$rho: matrix (S x (vT-1))
rho_draws <- fit$rho
S  <- nrow(rho_draws)
vT <- ncol(rho_draws) + 1L

# Build R for each draw: array (S x vT x vT)
R_draws <- array(NA_real_, dim = c(S, vT, vT))
for (s in 1:S) {
  R_draws[s, , ] <- AD_matrix(rho_draws[s, ])
}

# Point estimate: posterior mean matrix
R_mean <- apply(R_draws, c(2,3), mean)

# Optional: elementwise 95% CrI
R_lo <- apply(R_draws, c(2,3), quantile, probs = 0.025)
R_hi <- apply(R_draws, c(2,3), quantile, probs = 0.975)




df_Rmean <- as.data.frame(R_mean) %>%
  mutate(row = row_number() - 1L) %>%   # time index starting at 0
  pivot_longer(-row, names_to = "col", values_to = "value") %>%
  mutate(
    col = as.integer(gsub("V", "", col)) - 1L
  )

p_Rmean_gray <- ggplot(df_Rmean, aes(x = col, y = row, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
  coord_equal() +
  scale_y_reverse(breaks = 0:(vT-1)) +
  scale_x_continuous(breaks = 0:(vT-1)) +
  scale_fill_gradient(
    low = "white",
    high = "black",
    limits = c(0, 1),
    name = "Corr"
  ) +
  labs(
    x = "Time",
    y = "Time",
    title = "Estimated correlation matrix (posterior mean)"
  ) +
  theme_minimal(base_size = 12)

p_Rmean_gray

R_width <- R_hi - R_lo

df_Rw <- as.data.frame(R_width) %>%
  mutate(row = row_number() - 1L) %>%
  pivot_longer(-row, names_to = "col", values_to = "value") %>%
  mutate(
    col = as.integer(gsub("V", "", col)) - 1L
  )

p_Rw_gray <- ggplot(df_Rw, aes(x = col, y = row, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
  coord_equal() +
  scale_y_reverse(breaks = 0:(vT-1)) +
  scale_x_continuous(breaks = 0:(vT-1)) +
  scale_fill_gradient(
    low = "white",
    high = "black",
    name = "CrI width"
  ) +
  labs(
    x = "Time",
    y = "Time",
    title = "Uncertainty in correlation (95% CrI width)"
  ) +
  theme_minimal(base_size = 12)

p_Rw_gray



df_Rmean_low <- as.data.frame(R_mean) %>%
  mutate(row = row_number() - 1L) %>%                 # time index 0,...,vT-1
  pivot_longer(-row, names_to = "col", values_to = "value") %>%
  mutate(
    col = as.integer(gsub("V", "", col)) - 1L
  ) %>%
  filter(row >= col)                                  # keep lower triangle

p_Rmean_low <- ggplot(df_Rmean_low, aes(x = col, y = row, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
  coord_equal() +
  scale_y_reverse(breaks = 0:(vT-1)) +
  scale_x_continuous(breaks = 0:(vT-1)) +
  scale_fill_gradient(
    low = "white",
    high = "black",
    limits = c(0, 1),
    name = "Corr"
  ) +
  labs(
    x = "Time",
    y = "Time",
    title = "Estimated correlation matrix (posterior mean) — lower triangle"
  ) +
  theme_minimal(base_size = 12)

p_Rmean_low














