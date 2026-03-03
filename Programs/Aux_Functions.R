
## parâmetros globais

# Y: matriz de dados (n x T) em que "n" é o total de observações por tempo e "T" é o número de tempos (serve para dados balanceados).

# sigma2: vetor de variâncias marginais

# mu: vetor médias marginais

# Sigma: matriz de covariância



## Para rescalonar os dados simulados com vetor de médias "mu" e matriz de covariância "Sigma". Retorna a matriz de dados "Y", reescalonada.

#Y<-Y_C[,-1];mu<-mu_verd1;Sigma<-Sigma_verd

autocorr_lag <- function(stanfit, param, lag = 1, chain = NULL) {
  # Extrai amostras brutas (não permutadas)
  samples <- rstan::extract(stanfit, pars = param, permuted = FALSE)  # array: iter x chains x parameters
  
  # Seleciona dimensão correta: supõe parâmetro escalar
  if (is.null(chain)) {
    # Agrega todas as cadeias
    series <- as.vector(samples[, , 1])
  } else {
    if (!chain %in% 1:dim(samples)[2]) stop("Número de cadeia inválido.")
    series <- samples[, chain, 1]
  }
  
  # Verifica se lag é válido
  if (lag >= length(series)) stop("Lag maior ou igual ao tamanho da cadeia.")
  
  # Calcula autocorrelação
  n <- length(series)
  x <- series - mean(series)
  acf_k <- sum(x[1:(n - lag)] * x[(lag + 1):n]) / sum(x^2)
  
  return(acf_k)
}




resc_Y<-function(Y,mu,Sigma,centred=TRUE){
  
  n<-nrow(Y);m<-ncol(Y)
  
  mu_Y <- matrix(1,n,1)%x%(cbind(apply(Y,2,mean)))
  
  muY_mat<-matrix(mu_Y,n,m,byrow=T)
  
  if(centred==T){
    
    Y_c<-Y - muY_mat
    
  }else{
    
    Y_c<-Y;mu<-rep(0,m)
    
  }
  
  
  covY<-cov(Y,use="complete.obs")
  
  Y_z<-t(solve(chol(covY)))%*%t(Y_c)
  
  mu_mat<-matrix(matrix(1,n,1)%x%mu,n,m,byrow=T)
  
  Y_s<-t(t(chol(Sigma))%*%Y_z) + mu_mat
  
  return(Y_s)
  
}


## Matrizes de covari?ncia

# Unstructured

make_corr_unstructured <- function(m, target_rho=0.4, mix=0.85, jitter=0.07, eps=1e-6){
  A <- matrix(rnorm(m*m),m,m)
  S <- crossprod(A)
  R <- cov2cor(S)
  off <- R - diag(m)
  mu <- mean(off[upper.tri(off)])
  sc <- target_rho/mu
  off <- off*sc
  R2 <- diag(m)+mix*off
  eig <- eigen(R2,sym=TRUE)
  lam <- pmax(eig$values,eps)
  Rpd <- eig$vectors%*%diag(lam)%*%t(eig$vectors)
  D <- sqrt(diag(Rpd))
  Rpd <- Rpd / outer(D,D)
  diag(Rpd)<-1
  (Rpd+t(Rpd))/2
}





# HU

HU.matrix<-function(sigma2,phi){
  
  m<-length(sigma2)
  
  U<-matrix(phi,m,m)
  
  diag(U)<-rep(1,m)
  
  
  mat_aux1<-sqrt(matrix(sigma2[matrix(rep(1:m,m),m,m,byrow=T)],m,m))
  mat_aux2<-sqrt(matrix(sigma2[matrix(rep(1:m,m),m,m,byrow=F)],m,m))
  
  U<-U*mat_aux1*mat_aux2
  
  return(U)
  
  
}



toeplitz.matrix<-function(sigma2,rho){
  
  m<-length(sigma2)
  
  R<-toeplitz(c(1,rho))
  
  mat_aux1<-sqrt(matrix(sigma2[matrix(rep(1:m,m),m,m,byrow=T)],m,m))
  mat_aux2<-sqrt(matrix(sigma2[matrix(rep(1:m,m),m,m,byrow=F)],m,m))
  
  return(R*mat_aux1*mat_aux2)
  
}




ARMAH11.matrix<-function(sigma2,phi,gama){
  
  
  vT<-length(sigma2)
  
  R<-gama*cov2cor(gama.arma(toeplitz(c(0,0:(vT-2))),1,phi,0))
  
  
  mat_aux1<-sqrt(matrix(sigma2[matrix(rep(1:vT,vT),vT,vT,byrow=T)],vT,vT))
  
  Sigma<-R*mat_aux1*t(mat_aux1)
  
  diag(Sigma)<-sigma2
  
  
  return(Sigma)
  
}



ARH1.matrix<-function(sigma2,rho){
  
  m<-length(sigma2)
  
  if(rho==0){
    
    return(diag(sigma2))
  
  }else{  
  Sigma<-gama.arma(toeplitz(0:(m-1)),1,rho,0)
  
  R<-(1/Sigma[1,1])*Sigma
  
  mat_aux1<-sqrt(matrix(sigma2[matrix(rep(1:m,m),m,m,byrow=T)],m,m))
  mat_aux2<-sqrt(matrix(sigma2[matrix(rep(1:m,m),m,m,byrow=F)],m,m))
  
  Sigma<-R*mat_aux1*mat_aux2
  
  diag(Sigma)<-sigma2
  
  return(Sigma)
  
  }
  
}

Prodacf<-function(rho,h){
  
  acf<-c()
  
  acf[1]<-1
  
  for(i in 2:h){
    
    acf[i]<-rho[i-1]*acf[i-1]
    
  }
  
  return(acf)
  
}



Prod.matrix<-function(vT,rho){
  
#  m<-length(sigma2)
  
  R<-toeplitz(Prodacf(rho,vT))
  
#  mat_aux1<-sqrt(matrix(sigma2[matrix(rep(1:m,m),m,m,byrow=T)],m,m))
#  mat_aux2<-sqrt(matrix(sigma2[matrix(rep(1:m,m),m,m,byrow=F)],m,m))
  
#  Sigma<-R*mat_aux1*mat_aux2
  
  
  return(R)
  
}


AD_matrix <- function(rho) {
  # rho = vetor de comprimento L-1
  L<-length(rho) + 1L
  
  # Cria matriz identidade
  R <- diag(1, L)
  
  # Preenche triângulo superior e copia para inferior
  for (i in 1:(L-1)) {
    for (j in (i+1):L) {
      prod_ij <- prod(rho[i:(j-1)])  # produto dos rhos do intervalo
      R[i, j] <- prod_ij
      R[j, i] <- prod_ij
    }
  }
  
  return(R)
}




## Uniform matrix at [0,1] interval

uniform.matrix<-function(n,I){
  
  U<-matrix(runif(n*I),n,I)
  
  return(U)
  
}


MAP<-function(x){
  dd<-density(x)
  return(dd$x[which.max(dd$y)])
}


## autocovariance function of an ARMA process

gama.arma<-function(k,sigma2,phi,theta){
  
  gama.0<- sigma2*((1 + 2*theta*phi + theta^2) / (1-phi^2))
  
  gama.1<-sigma2*((1+theta*phi)*(phi+theta)/(1-phi^2))
  
  gama.k<-(phi^(k-1))*gama.1
  
  return(gama.0*(k==0) + gama.1*(k==1) + gama.k*(k>=2))
  
}

#' Reconstrói a matriz de correlação a partir de um summary de MCMC
#'
#' @param summary_df Data frame resultante da leitura de um arquivo .rds com parâmetros nomeados como "C[i,j]"
#' @param value_col Nome da coluna a ser usada (default = "mean")
#' @return Matriz de correlação estimada
reconstruir_correlacao <- function(summary_df, value_col = "mean") {
  # Identifica linhas com o padrão "C[i,j]"
  idx_cor <- grep("^C\\[", rownames(summary_df))
  cor_rows <- rownames(summary_df)[idx_cor]
  
  # Extrai o maior índice para definir a dimensão da matriz
  max_idx <- max(as.numeric(unlist(regmatches(cor_rows, gregexpr("[0-9]+", cor_rows)))))
  T <- max_idx
  
  # Inicializa matriz
  R_hat <- matrix(NA, nrow = T, ncol = T)
  
  # Preenche matriz com os valores extraídos
  for (row in cor_rows) {
    pos <- as.numeric(unlist(regmatches(row, gregexpr("[0-9]+", row))))
    i <- pos[1]
    j <- pos[2]
    R_hat[i, j] <- summary_df[row, value_col]
  }
  
  # Garante simetria
  R_hat[lower.tri(R_hat)] <- t(R_hat)[lower.tri(R_hat)]
  
  # Atribui nomes
  rownames(R_hat) <- colnames(R_hat) <- paste0("V", 1:T)
  
  return(R_hat)
}


# ===========================================================================
# Train / Test split for multivariate longitudinal probit data
# (subject-level split, preserves within-subject dependence)
# ===========================================================================

make_train_test_split <- function(Y, X, test_frac = 0.2, seed = 1234) {
  
  set.seed(seed)
  
  # Number of subjects
  N <- dim(Y)[1]
  
  # Indices of subjects
  id_all <- seq_len(N)
  
  # Sample test subjects
  n_test <- ceiling(test_frac * N)
  id_test <- sample(id_all, size = n_test, replace = FALSE)
  id_train <- setdiff(id_all, id_test)
  
  # Split response
  Y_train <- Y[id_train, , drop = FALSE]
  Y_test  <- Y[id_test,  , drop = FALSE]
  
  # Split design array (N x vT x p)
  X_train <- X[id_train, , , drop = FALSE]
  X_test  <- X[id_test,  , , drop = FALSE]
  
  list(
    id_train = id_train,
    id_test  = id_test,
    Y_train  = Y_train,
    Y_test   = Y_test,
    X_train  = X_train,
    X_test   = X_test
  )
}

# Initial values by model
make_inits <- function(model_name, p, vT) {
  if (model_name %in% c("Toeplitz", "AD")) {
    return(function() list(beta = rep(0.1, p), rho_vec = rep(0.1, vT - 1)))
  }
  if (model_name == "AR1") {
    return(function() list(beta = rep(0.1, p), rho_ar1 = 0.1))
  }
  if (model_name == "ARMA11") {
    return(function() list(beta = rep(0.1, p), phi = 0.1, theta = 0.1))
  }
  if (model_name %in% c("Unstructured", "Independent")) {
    return(function() list(beta = rep(0.1, p)))
  }
  stop("Unknown model in make_inits().")
}

# Parameters to monitor (with optional log_lik)
make_pars <- function(model_name, log_lik = TRUE) {
  
  # Base parameter sets by model
  pars <- switch(model_name,
                 "Toeplitz"     = c("beta", "rho_vec"),
                 "AD"           = c("beta", "rho_vec"),
                 "AR1"          = c("beta", "rho_ar1"),
                 "ARMA11"       = c("beta", "phi", "theta"),
                 "Unstructured" = c("beta", "C"),
                 "Independent"  = c("beta"),
                 stop("Unknown model in make_pars().")
  )
  
  # Optionally include log_lik
  if (log_lik) {
    pars <- c(pars, "log_lik")
  }
  
  pars
}



# Compute model-implied marginal probabilities under probit:
# p = Phi(mu), where mu = X * beta (marginal variance = 1)
probit_probs <- function(X, beta) {
  # X: N x vT x p array
  N  <- dim(X)[1]
  vT <- dim(X)[2]
  mu <- matrix(NA_real_, N, vT)
  for (n in 1:N) mu[n, ] <- X[n, , ] %*% beta
  pnorm(mu)
}

# Yen's Q3 (correlation of Pearson residuals across persons)
q3_yen <- function(Y, P, eps = 1e-10) {
  # Y, P: N x vT matrices
  Y <- as.matrix(Y)
  P <- pmin(pmax(P, eps), 1 - eps)
  R <- (Y - P) / sqrt(P * (1 - P))
  Q3 <- stats::cor(R)
  diag(Q3) <- NA_real_
  Q3
}

# Optional: a global discrepancy from the Q3 matrix (choose one)
q3_global <- function(Q3, type = c("max_abs", "frobenius", "mean_abs")) {
  type <- match.arg(type)
  off <- Q3[upper.tri(Q3)]
  if (type == "max_abs") return(max(abs(off), na.rm = TRUE))
  if (type == "mean_abs") return(mean(abs(off), na.rm = TRUE))
  if (type == "frobenius") return(sqrt(sum(off^2, na.rm = TRUE)))
}

# --- Correlation builders (structured) ---------------------------------
Toep_matrix_R <- function(L, rho) {
  nr <- length(rho)
  R <- diag(1, L)
  for (i in 1:L) for (j in 1:L) {
    if (i != j) {
      d <- abs(i - j)
      R[i, j] <- if (nr == 1) rho[1] else rho[d]
    }
  }
  R
}

AR1_matrix_R <- function(L, rho) {
  outer(1:L, 1:L, function(i, j) rho^abs(i - j))
}

AD_matrix_R <- function(L, rho_vec) {
  R <- diag(1, L)
  for (i in 1:(L - 1)) for (j in (i + 1):L) {
    R[i, j] <- prod(rho_vec[i:(j - 1)])
    R[j, i] <- R[i, j]
  }
  R
}

ARMA11_matrix_R <- function(L, phi, theta) {
  denom <- 1 + theta^2 + 2 * phi * theta
  cst <- ((phi + theta) * (1 + phi * theta)) / denom
  R <- diag(1, L)
  for (i in 1:(L - 1)) for (j in (i + 1):L) {
    h <- j - i
    R[i, j] <- (phi^(h - 1)) * cst
    R[j, i] <- R[i, j]
  }
  R
}

# --- MVN simulator for test set ----------------------------------------
simulate_Yrep_test <- function(X_test, beta, R) {
  # returns Y_rep (N_test x vT)
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Package 'mvtnorm' is required: install.packages('mvtnorm').")
  }
  N_test <- dim(X_test)[1]
  vT     <- dim(X_test)[2]
  
  Yrep <- matrix(0L, N_test, vT)
  for (n in 1:N_test) {
    mu_n <- as.vector(X_test[n, , ] %*% beta)
    zrep <- mvtnorm::rmvnorm(1, mean = mu_n, sigma = R)
    Yrep[n, ] <- as.integer(zrep > 0)
  }
  Yrep
}


# ===========================================================================
# Posterior Holdout Predictive Checking (HPC) using Yen's Q3 as discrepancy
# ---------------------------------------------------------------------------
# This function implements the *adjusted / realized* Bayesian p-value:
#   p_HPC = E[ I{ d(theta, y_rep) > d(theta, y_new) } | y_train ]
# where:
#   theta ~ p(theta | y_train),
#   y_rep ~ p(y | theta)  (replicated holdout data),
#   y_new = y_test        (observed holdout data),
# and the *same* posterior draw theta is used in both terms.
# ---------------------------------------------------------------------------
# Inputs:
#   fit        : stanfit object fitted on training data
#   model_name : one of {"Toeplitz","AR1","AD","ARMA11","Independent","Unstructured"}
#   X_test     : array [N_test, vT, p]
#   Y_test     : matrix/data.frame [N_test, vT] with 0/1
#   vT         : number of time points
#   S          : number of posterior draws to use
#   seed       : RNG seed
#   cor_type   : required for structured models (1 Toeplitz, 2 AR1, 3 AD, 4 ARMA11)
#   eps        : small constant to avoid division by zero in Pearson residuals
#   global_type: global Q3 discrepancy ("max_abs","mean_abs","frobenius")
# Output:
#   list with observed Q3, replicated Q3, new-data Q3 (draw-specific),
#   pairwise adjusted p-values, global adjusted p-value, and draws used.
# ===========================================================================

hpc_q3_adjusted <- function(
    fit, model_name,
    X_test, Y_test, vT,
    S = 500, seed = 123,
    cor_type = NULL,
    eps = 1e-10,
    global_type = c("max_abs","mean_abs","frobenius")
) {
  
  global_type <- match.arg(global_type)
  set.seed(seed)
  
  # Basic checks
  Y_test <- as.matrix(Y_test)
  stopifnot(dim(X_test)[2] == vT, ncol(Y_test) == vT, nrow(Y_test) == dim(X_test)[1])
  
  # --- Extract posterior draws --------------------------------------------
  post <- rstan::extract(fit, permuted = TRUE)
  
  if (is.null(post$beta)) stop("Posterior draw 'beta' not found in fit.")
  beta_draws <- post$beta
  n_draws <- dim(beta_draws)[1]
  
  idx <- sample(seq_len(n_draws), size = min(S, n_draws), replace = FALSE)
  S_used <- length(idx)
  
  # --- Storage ------------------------------------------------------------
  Q3_rep_array <- array(NA_real_, dim = c(vT, vT, S_used))
  Q3_new_array <- array(NA_real_, dim = c(vT, vT, S_used))
  
  glob_rep <- numeric(S_used)
  glob_new <- numeric(S_used)
  
  # --- Main loop: adjusted/realized HPC -----------------------------------
  for (s in seq_along(idx)) {
    it <- idx[s]
    beta_s <- beta_draws[it, ]
    
    # Build correlation matrix R for this draw
    if (model_name %in% c("Toeplitz", "AR1", "AD", "ARMA11")) {
      
      if (is.null(cor_type)) {
        stop("For structured models, provide cor_type (1 Toeplitz, 2 AR1, 3 AD, 4 ARMA11).")
      }
      
      R <- switch(as.character(cor_type),
                  "1" = {
                    if (is.null(post$rho_vec)) stop("Posterior draw 'rho_vec' not found.")
                    Toep_matrix_R(vT, post$rho_vec[it, ])
                  },
                  "2" = {
                    if (is.null(post$rho_ar1)) stop("Posterior draw 'rho_ar1' not found.")
                    AR1_matrix_R(vT, post$rho_ar1[it])
                  },
                  "3" = {
                    if (is.null(post$rho_vec)) stop("Posterior draw 'rho_vec' not found.")
                    AD_matrix_R(vT, post$rho_vec[it, ])
                  },
                  "4" = {
                    if (is.null(post$phi) || is.null(post$theta)) stop("Posterior draws 'phi'/'theta' not found.")
                    ARMA11_matrix_R(vT, post$phi[it], post$theta[it])
                  },
                  stop("Invalid cor_type.")
      )
      
    } else if (model_name == "Unstructured") {
      
      R <- if (!is.null(post$Lcorr)) {
        L <- post$Lcorr[it, , ]
        L %*% t(L)
      } else if (!is.null(post$C)) {
        post$C[it, , ]
      } else {
        stop("Neither 'Lcorr' nor 'C' found in posterior draws.")
      }
      
    } else if (model_name == "Independent") {
      
      R <- diag(1, vT)
      
    } else {
      stop("Unknown model_name.")
    }
    
    # Draw-specific marginal probabilities (same theta for rep and new)
    P_s <- probit_probs(X_test, beta_s)
    
    # Replicated holdout data under theta_s
    Yrep <- simulate_Yrep_test(X_test, beta_s, R)
    
    # Q3 discrepancies under the SAME P_s
    Q3_rep <- q3_yen(Yrep,   P_s, eps = eps)
    Q3_new <- q3_yen(Y_test, P_s, eps = eps)
    
    Q3_rep_array[, , s] <- Q3_rep
    Q3_new_array[, , s] <- Q3_new
    
    # Global summaries
    glob_rep[s] <- q3_global(Q3_rep, type = global_type)
    glob_new[s] <- q3_global(Q3_new, type = global_type)
  }
  
  # --- Adjusted (realized) Bayesian p-values ------------------------------
  ppp_pair <- matrix(NA_real_, vT, vT)
  for (i in 1:vT) for (j in 1:vT) if (i != j) {
    # Adjusted/realized p-value: same theta draw in both terms
    ppp_pair[i, j] <- mean(Q3_rep_array[i, j, ] > Q3_new_array[i, j, ], na.rm = TRUE)
  }
  diag(ppp_pair) <- NA_real_
  
  ppp_global <- mean(glob_rep > glob_new, na.rm = TRUE)
  
  # For reporting: posterior mean of the "new-data" Q3 under draw-specific theta
  Q3_new_mean <- apply(Q3_new_array, c(1, 2), mean, na.rm = TRUE)
  diag(Q3_new_mean) <- NA_real_
  
  list(
    model_name  = model_name,
    cor_type    = cor_type,
    draws_used  = idx,
    global_type = global_type,
    
    Q3_rep = Q3_rep_array,
    Q3_new = Q3_new_array,
    Q3_new_mean = Q3_new_mean,
    
    ppp_pair   = ppp_pair,
    glob_rep   = glob_rep,
    glob_new   = glob_new,
    ppp_global = ppp_global
  )
}


Q3_expected <- function(Q3_rep_array, fun = c("median", "mean")) {
  fun <- match.arg(fun)
  if (fun == "median") {
    Q3_exp <- apply(Q3_rep_array, c(1, 2), median, na.rm = TRUE)
  } else {
    Q3_exp <- apply(Q3_rep_array, c(1, 2), mean, na.rm = TRUE)
  }
  diag(Q3_exp) <- NA_real_
  Q3_exp
}



Q3_discrepancy <- function(Q3_new_mean, Q3_rep_array, exp_fun = "median") {
  Q3_exp <- Q3_expected(Q3_rep_array, fun = exp_fun)
  D <- Q3_new_mean - Q3_exp
  diag(D) <- NA_real_
  D
}




plot_Q3_discrepancy_heatmap <- function(D,
                                        title = "",
                                        fill_limits = NULL,
                                        show_labels = TRUE,
                                        digits = 2) {
  
  vT <- nrow(D)
  
  df <- expand.grid(t1 = 1:vT, t2 = 1:vT)
  df$D <- as.vector(D)
  df$absD <- abs(df$D)
  
  # Labels with sign (+/-) and rounding
  if (show_labels) {
    df$lab <- ifelse(is.na(df$D), "",
                     formatC(df$D, format = "f", digits = digits))
  } else {
    df$lab <- ""
  }
  
  ggplot(df, aes(x = factor(t1), y = factor(t2), fill = absD)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = lab), size = 3) +
    scale_fill_gradient(
      low = "white",
      high = "black",
      limits = fill_limits,
      na.value = "white",
      name = "|ΔQ3|"
    ) +
    coord_fixed() +
    labs(title = title, x = "Time", y = "Time") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}



plot_Q3_discrepancy_independent <- function(
    hpc_independent,
    hpc_model,
    model_label = "Model",
    common_scale = TRUE,
    exp_fun = "median",
    show_labels = TRUE
) {
  
  D_ind <- Q3_discrepancy(
    Q3_new_mean = hpc_independent$Q3_new_mean,
    Q3_rep_array = hpc_independent$Q3_rep,
    exp_fun = exp_fun
  )
  
  D_mod <- Q3_discrepancy(
    Q3_new_mean = hpc_model$Q3_new_mean,
    Q3_rep_array = hpc_model$Q3_rep,
    exp_fun = exp_fun
  )
  
  lim <- NULL
  if (common_scale) {
    lim <- range(c(abs(D_ind), abs(D_mod)), na.rm = TRUE)
  }
  
  p1 <- plot_Q3_discrepancy_heatmap(
    D_ind,
    title = "Independent model",
    fill_limits = lim,
    show_labels = show_labels
  )
  
  p2 <- plot_Q3_discrepancy_heatmap(
    D_mod,
    title = model_label,
    fill_limits = lim,
    show_labels = show_labels
  )
  
  list(
    D_independent = D_ind,
    D_model = D_mod,
    plot_independent = p1,
    plot_model = p2
  )
}



plot_Q3_discrepancy_panel <- function(
    hpc_independent,
    hpc_model,
    model_label = "Model",
    exp_fun = "median",
    show_labels = TRUE
) {
  
  # Signed discrepancies
  D_ind <- Q3_discrepancy(
    Q3_new_mean = hpc_independent$Q3_new_mean,
    Q3_rep_array = hpc_independent$Q3_rep,
    exp_fun = exp_fun
  )
  
  D_mod <- Q3_discrepancy(
    Q3_new_mean = hpc_model$Q3_new_mean,
    Q3_rep_array = hpc_model$Q3_rep,
    exp_fun = exp_fun
  )
  
  # Common scale (VERY important for interpretation)
  lim <- range(c(abs(D_ind), abs(D_mod)), na.rm = TRUE)
  
  p_ind <- plot_Q3_discrepancy_heatmap(
    D_ind,
    title = "Independent",
    fill_limits = lim,
    show_labels = show_labels
  )
  
  p_mod <- plot_Q3_discrepancy_heatmap(
    D_mod,
    title = model_label,
    fill_limits = lim,
    show_labels = show_labels
  )
  
  # Combine with ONE shared legend
  (p_ind | p_mod) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")
}






