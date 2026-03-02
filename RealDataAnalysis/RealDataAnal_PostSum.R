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
Y <- six_cities[, -5]



##================================= Fitting models ========================================##

# Escolha da estrutura de correlação
# 1 = Toeplitz
# 2 = AR(1)
# 3 = AD(1)
# 4 = ARMA(1,1)

data_list=list(vT=vT,p=p,N=n,Y=Y,X=X_array,sigma_beta=10,sigma_rho=1,cor_type=2) 


model_name<-'AR1'

params <- make_pars(model_name,log_lik = F)
nChains = 1
burnInSteps = 1000
thinSteps = 10
numSavedSteps = 1000  #across all chains
nIter = ceiling(burnInSteps + (numSavedSteps * thinSteps)/nChains)
nIter

#nIter=11000


ini<-make_inits(model_name, p = p, vT = vT)

begin<-Sys.time()
print(begin)


samp <- stan(data = data_list, file =file.path('Programs', "FitMultProbit_Structured.stan"),init = ini,
             chains = nChains, pars = params, iter = nIter,
             warmup = burnInSteps, thin = thinSteps,control = list(adapt_delta = 0.9,max_treedepth=15),
             save_dso = T)


end<-Sys.time()
print(end-begin)


saveRDS(samp,paste0(pathFit,'/Samples_',model_name,'_SixCities.rds'))


samp<-readRDS(paste0(pathFit,'/Samples_',model_name,'_SixCities.rds'))


#launch_shinystan(samp)


fit <- rstan::extract(samp,permuted=T,inc_warmup=F)

#-----------------------------#
# Helper: posterior summaries #
#-----------------------------#
post_summ <- function(draws, par_names = NULL,
                      prob = c(0.025, 0.5, 0.975),
                      drop = FALSE) {
  # draws can be:
  #  - vector (S)                            -> single parameter (e.g., rho_ar1)
  #  - matrix (S x K)                        -> multiple parameters
  #  - array  (S x K x ...) (e.g., Lcorr)    -> will be flattened to S x K
  
  # 1) Coerce to S x K matrix
  if (is.vector(draws)) {
    draws_mat <- matrix(draws, ncol = 1)
  } else if (is.matrix(draws)) {
    draws_mat <- draws
  } else if (is.array(draws)) {
    S <- dim(draws)[1]
    draws_mat <- matrix(draws, nrow = S)  # flatten remaining dims into columns
  } else {
    stop("'draws' must be a vector, matrix, or array.")
  }
  
  # 2) Default parameter names
  K <- ncol(draws_mat)
  if (is.null(par_names)) {
    par_names <- if (K == 1) "param" else paste0("param[", seq_len(K), "]")
  } else {
    if (length(par_names) != K) {
      stop("Length of 'par_names' must match number of parameters (ncol(draws)).")
    }
  }
  
  # 3) Summaries
  qs <- apply(draws_mat, 2, quantile, probs = prob, na.rm = TRUE)
  
  out <- data.frame(
    mean   = colMeans(draws_mat, na.rm = TRUE),
    sd     = apply(draws_mat, 2, sd, na.rm = TRUE),
    q2.5   = qs[1, ],
    median = qs[2, ],
    q97.5  = qs[3, ],
    pr_gt0 = colMeans(draws_mat > 0, na.rm = TRUE),
    row.names = par_names,
    check.names = FALSE
  )
  
  if (K == 1 && drop) return(out[1, , drop = FALSE])
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
# 2) Rhos (AR1)                #
#-----------------------------#
# Common naming in your Stan workflow: rho is (S x (vT-1)) for AD(1),
# where rho[t] is the lag-1 partial correlation parameter at transition t -> t+1.
stopifnot(!is.null(fit$rho_ar1))
#K <- ncol(fit$rho_ar1)

#rho_names <- paste0("rho[", 1:K, "]")   # or use rho_0,...,rho_{K-1} if you prefer
rho_sum <- post_summ(fit$rho_ar1,par_names = 'rho')
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




#-----------------------------#
# Helper: posterior summaries #
#-----------------------------#
post_ci_df <- function(draws,
                       par_names = NULL,
                       level = 0.95,
                       point = c("mean","median")) {
  point <- match.arg(point)
  alpha <- (1 - level) / 2
  probs <- c(alpha, 0.5, 1 - alpha)
  
  # Coerce draws to an S x K matrix
  if (is.null(dim(draws))) {
    # vector case (e.g., rho_ar1)
    draws <- matrix(draws, ncol = 1)
  } else if (is.matrix(draws)) {
    # ok
  } else if (is.array(draws)) {
    # flatten arrays (e.g., if ever used)
    S <- dim(draws)[1]
    draws <- matrix(draws, nrow = S)
  } else {
    stop("'draws' must be a vector, matrix, or array.")
  }
  
  K <- ncol(draws)
  
  # Handle parameter names
  if (is.null(par_names)) {
    par_names <- if (K == 1) "par" else paste0("par[", seq_len(K), "]")
  } else {
    if (length(par_names) != K) {
      stop("Length of 'par_names' must match ncol(draws).")
    }
  }
  
  qs <- apply(draws, 2, quantile, probs = probs, na.rm = TRUE)
  est <- if (point == "mean") colMeans(draws, na.rm = TRUE) else qs[2, ]
  
  tibble::tibble(
    parameter = par_names,
    estimate  = as.numeric(est),
    lower     = as.numeric(qs[1, ]),
    upper     = as.numeric(qs[3, ])
  )
}

#-----------------------------------------#
# Helper: choose correlation draws by model
#-----------------------------------------#
get_cor_draws <- function(fit, model_name, vT = NULL) {
  
  if (model_name %in% c("Toeplitz", "AD")) {
    if (is.null(fit$rho_vec)) stop("fit$rho_vec not found.")
    return(list(draws = fit$rho_vec, names = paste0("rho[", 1:ncol(fit$rho_vec), "]")))
  }
  
  if (model_name == "AR1") {
    if (is.null(fit$rho_ar1)) stop("fit$rho_ar1 not found.")
    return(list(draws = fit$rho_ar1, names = "rho"))
  }
  
  if (model_name == "ARMA11") {
    if (is.null(fit$phi) || is.null(fit$theta)) stop("fit$phi and/or fit$theta not found.")
    draws <- cbind(phi = fit$phi, theta = fit$theta)
    return(list(draws = draws, names = c("phi", "theta")))
  }
  
  if (model_name == "Independent") {
    # No correlation parameters
    return(list(draws = NULL, names = NULL))
  }
  
  if (model_name == "Unstructured") {
    if (is.null(fit$Lcorr)) stop("fit$Lcorr not found.")
    if (is.null(vT)) stop("For Unstructured, provide vT.")
    
    # Convert Lcorr draws (S x vT x vT) to correlation entries (S x K) for i<j
    Larr <- fit$Lcorr
    S <- dim(Larr)[1]
    if (length(dim(Larr)) != 3 || dim(Larr)[2] != vT || dim(Larr)[3] != vT) {
      stop("fit$Lcorr must have dimension S x vT x vT.")
    }
    
    pairs <- which(upper.tri(matrix(1, vT, vT)), arr.ind = TRUE)
    K <- nrow(pairs)
    cor_draws <- matrix(NA_real_, nrow = S, ncol = K)
    cor_names <- character(K)
    
    for (k in 1:K) {
      i <- pairs[k, 1]; j <- pairs[k, 2]
      # For each draw: C = L %*% t(L)
      # Compute C_ij efficiently per draw
      # (still looped, but vT is small)
      cor_draws[, k] <- vapply(1:S, function(s) {
        L <- Larr[s, , ]
        C <- L %*% t(L)
        C[i, j]
      }, numeric(1))
      cor_names[k] <- paste0("Corr[", i, ",", j, "]")
    }
    
    return(list(draws = cor_draws, names = cor_names))
  }
  
  stop("Unknown model_name in get_cor_draws().")
}


#-----------------------------#
# One plot: betas + AR(1) rho  #
#-----------------------------#

beta_names <- c("(Intercept)", "age", "mother's smoking", "Interaction")
if (is.null(fit$beta)) stop("fit$beta not found.")
if (ncol(fit$beta) != length(beta_names)) beta_names <- paste0("beta[", 1:ncol(fit$beta), "]")

# Betas
df_beta <- post_ci_df(fit$beta, par_names = beta_names, level = 0.95, point = "mean") %>%
  mutate(block = "Regression coefficients")

# AR(1) correlation parameter
if (is.null(fit$rho_ar1)) stop("fit$rho_ar1 not found (AR1).")
df_rho <- post_ci_df(fit$rho_ar1, par_names = "rho", level = 0.95, point = "mean") %>%
  mutate(block = "Dependence (AR1)")

# Combine and order (betas first, rho last)
df_all <- bind_rows(df_beta, df_rho) %>%
  mutate(parameter = factor(parameter, levels = rev(parameter)))  # top-to-bottom order

p_all <- ggplot(df_all, aes(x = parameter, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  geom_point(size = 2) +
  coord_flip() +
  facet_wrap(~ block, scales = "free_y", ncol = 1) +
  labs(x = NULL, y = "Estimate",
       title = "") +
  theme_minimal(base_size = 12)

p_all


ggsave(
  file.path(pathFig, "Est_ar1_SixCities.eps"),
  p_all,
  device = cairo_ps,
  width  = 6,
  height = 5,
  units  = "in"
)




#------------------------------------------------------------#
# Heatmap: estimated AR(1) correlation matrix (posterior mean)
#------------------------------------------------------------#

# rho_ar1 draws: vector length S
rho_draws <- fit$rho_ar1
S <- length(rho_draws)

# vT must be known (recommended). If not, set it explicitly:
# vT <- 4
if (!exists("vT")) stop("Please define vT (number of time points).")

# Build R for each draw: array (S x vT x vT)
R_draws <- array(NA_real_, dim = c(S, vT, vT))
for (s in 1:S) {
  R_draws[s, , ] <- AR1_matrix_R(vT, rho_draws[s])
}

# Posterior mean correlation matrix
R_mean <- apply(R_draws, c(2, 3), mean, na.rm = TRUE)

# (Optional) elementwise 95% CrI
R_lo <- apply(R_draws, c(2, 3), quantile, probs = 0.025, na.rm = TRUE)
R_hi <- apply(R_draws, c(2, 3), quantile, probs = 0.975, na.rm = TRUE)

# Long DF: lower triangle only (including diagonal)
df_Rmean_low <- as.data.frame(R_mean) %>%
  mutate(row = row_number() - 1L) %>%                 # time index 0,...,vT-1
  pivot_longer(-row, names_to = "col", values_to = "value") %>%
  mutate(
    col = as.integer(gsub("^V", "", col)) - 1L
  ) %>%
  filter(row >= col)

# Heatmap (grayscale) with values
p_Rmean_low <- ggplot(df_Rmean_low, aes(x = col, y = row, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
  coord_equal() +
  scale_y_reverse(breaks = 0:(vT - 1)) +
  scale_x_continuous(breaks = 0:(vT - 1)) +
  scale_fill_gradient(
    low = "white",
    high = "black",
    limits = c(0, 1),
    name = "Corr"
  ) +
  labs(
    x = "Time",
    y = "Time",
    title = ""
  ) +
  theme_minimal(base_size = 12)

p_Rmean_low



ggsave(
  file.path(pathFig, "CorMatrix_ar1_SixCities.eps"),
  p_Rmean_low,
  device = cairo_ps,
  width  = 6,
  height = 5,
  units  = "in"
)










