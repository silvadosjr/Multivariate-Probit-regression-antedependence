library(rstan)
library(rstansim)
library(doParallel)
library(foreach)

options(mc.cores = parallel::detectCores())

PC <- 'Usuário'
setwd(paste0('C:/Users/', PC, '/OneDrive/Documentos/Artigos/Multivariate Probit regression using antedependence/'))

library(here)

# Caminhos organizados
programas_dir  <- here("Programas")
pathDataSets   <- here("LatinSquareScenarios", "Toeplitz")
pathResults    <- here("ResultadosSimRobus")
source(file.path(programas_dir, "Aux_Functions.R"))

Scenario <- rbind(
  c(4,50,.9,.7), c(8,50,.4,.1), c(12,50,.7,.4),
  c(4,100,.7,.4), c(8,100,.9,.7), c(12,100,.4,.1),
  c(4,200,.4,.1), c(8,200,.7,.4), c(12,200,.9,.7)
)

R <- 20

# Registro do cluster paralelo
nCores <- parallel::detectCores() - 1  # deixe 1 núcleo livre
cl <- makeCluster(nCores)
registerDoParallel(cl)

for (j in 1:nrow(Scenario)) {
  
  vT        <- Scenario[j, 1]
  n         <- Scenario[j, 2]
  rho_start <- Scenario[j, 3]
  rho_end   <- Scenario[j, 4]
  
  t_age <- 1:vT - mean(1:vT)
  X_s   <- cbind(1, t_age, 1, t_age)
  X_ns  <- cbind(1, t_age, 0, 0)
  treat_counts <- rep(n / 2, 2)
  p <- 4
  
  X_array <- array(NA, dim = c(length(t_age), p, n))
  X_array[, , 1:treat_counts[1]] <- replicate(treat_counts[1], X_ns, simplify = "array")
  X_array[, , (treat_counts[1]+1):n] <- replicate(treat_counts[2], X_s, simplify = "array")
  X_array <- aperm(X_array, c(3, 1, 2))
  
  beta_verd <- c(1, .5, .8, .6)
  rho_verd <- round(seq(rho_start, to = rho_end, length.out = vT - 1), 2)
  
  corr_label <- if (rho_start == 0.4) {
    "Weak"
  } else if (rho_start == 0.7) {
    "Moderate"
  } else {
    "Strong"
  }
  
  # Paralelização das R réplicas
  foreach(r = 1:R, .packages = c("rstan", "here")) %dopar% {
    
    rds_path <- file.path(pathDataSets, paste0("SimMultProbit_Toeplitz_", vT, "_", n, "_", corr_label, "_", r, ".rds"))
    
    if (!file.exists(rds_path)) {
      message("Arquivo ausente: ", rds_path)
    } else {
      sim_Toeplitz <- readRDS(rds_path)
      
      data_list <- list(vT = vT, p = p, N = n, Y = sim_Toeplitz$Y, X = X_array,
                        sigma_beta = 1, sigma_rho = 1)
      
      params_Lcorr <- c('beta', 'C')
      params_AR1 <- c('beta', 'rho')
      
      ini_L <- function() list(beta = rep(0.1, p))
      ini_A <- function() list(beta = rep(0.1, p), rho = 0.1)
      
      burnInSteps <- 1000
      thinSteps   <- 15
      numSavedSteps <- 1000
      nChains <- 1
      nIter <- ceiling(burnInSteps + (numSavedSteps * thinSteps) / nChains)
      
      # Modelo Lcorr
      samp_Lcorr <- stan(
        data = data_list,
        file = file.path(programas_dir, "FitMultProbit_Lcorr.stan"),
        init = ini_L,
        chains = nChains,
        pars = params_Lcorr,
        iter = nIter,
        warmup = burnInSteps,
        thin = thinSteps,
        control = list(adapt_delta = 0.8, max_treedepth = 10),
        save_dso = TRUE,
        refresh = 0
      )
      
      # Modelo AR1
      samp_AR1 <- stan(
        data = data_list,
        file = file.path(programas_dir, "FitMultProbit_AR1.stan"),
        init = ini_A,
        chains = nChains,
        pars = params_AR1,
        iter = nIter,
        warmup = burnInSteps,
        thin = thinSteps,
        control = list(adapt_delta = 0.8, max_treedepth = 10),
        save_dso = TRUE,
        refresh = 0
      )
      
      # Salvar resultados
      saveRDS(summary(samp_Lcorr)$summary,
              file.path(pathResults, paste0("summaryMult_Lcorr_R", r, "_S", j, ".rds")))
      saveRDS(summary(samp_AR1)$summary,
              file.path(pathResults, paste0("summaryMult_AR1_R", r, "_S", j, ".rds")))
    }
  }
  
  message("Scenario ", j, " completed.")
}

stopCluster(cl)
