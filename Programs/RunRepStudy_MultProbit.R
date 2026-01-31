library(rstan)
library(rstansim)

options(mc.cores = parallel::detectCores())

PC <- 'Usuário'

setwd(paste0('C:/Users/', PC, '/OneDrive/Documentos/Artigos/Multivariate Probit regression using antedependence/'))

library(here)

# Escolha da estrutura de correlação
corr_type <- "Toeplitz"  # ou "AR1"

# Controle: quais modelos ajustar?
# Pode ser, por exemplo, c("Lcorr"), c("AR1","Ind"), c("Lcorr","AR1","Ind"), etc.
fit_models <- c("Ind")

# Caminhos organizados com `here()`
programas_dir  <- here("Programas")
pathDataSets   <- here("LatinSquareScenarios", corr_type)
pathResults    <- here("ResultadosSimRobus")

# Funções auxiliares
source(file.path(programas_dir, "Aux_Functions.R"))

## Correlation ranges
#[.10, .40] Weak; [.40, .70] moderate; [.70, .9] strong

Scenario <- rbind(
  c(4,  50, .9, .7),
  c(8,  50, .4, .1),
  c(12, 50, .7, .4),
  c(4,  100, .7, .4),
  c(8,  100, .9, .7),
  c(12, 100, .4, .1),
  c(4,  200, .4, .1),
  c(8,  200, .7, .4),
  c(12, 200, .9, .7)
)

R <- 20

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
  X_array[, , (treat_counts[1] + 1):n] <- replicate(treat_counts[2], X_s, simplify = "array")
  X_array <- aperm(X_array, c(3, 1, 2))
  
  beta_verd <- c(1, .5, .8, .6)
  
  rho_verd <- if (corr_type == "AR1") {
    rho_start
  } else if (corr_type == "Toeplitz") {
    round(seq(rho_start, to = rho_end, length.out = vT - 1), 2)
  } else {
    stop("corr_type must be either 'AR1' or 'Toeplitz'")
  }
  
  corr_label <- if (rho_start == 0.4) {
    "Weak"
  } else if (rho_start == 0.7) {
    "Moderate"
  } else {
    "Strong"
  }
  
  # Controle de parâmetros / MCMC
  params_Lcorr <- c('beta', 'C')
  params_AR1   <- c('beta', 'rho')
  params_Ind   <- c('beta')
  
  nChains      <- 1
  burnInSteps  <- 1000
  thinSteps    <- 15
  numSavedSteps <- 1000  # across all chains
  nIter <- ceiling(burnInSteps + (numSavedSteps * thinSteps) / nChains)
  
  # Função de início genérica: inclui beta e rho; Stan ignora parâmetros extras.
  ini <- function() {
    list(
      beta = rep(.1, p),
      rho  = rep(.1, max(vT - 1, 1)),  # para AR1 / Toeplitz se necessário
      C    = diag(vT)                  # se o modelo Lcorr usar algo tipo C (Cholesky/matriz)
    )
  }
  
  for (r in 1:R) {
    
    # Carregar dados simulados
    rds_path <- file.path(
      pathDataSets,
      paste0("SimMultProbit_", corr_type, "_", vT, "_", n, "_", corr_label, "_", r, ".rds")
    )
    sim_Toeplitz <- readRDS(rds_path)
    
    ##================================= Fitting models ========================================##
    
    data_list <- list(
      vT         = vT,
      p          = p,
      N          = n,
      Y          = sim_Toeplitz$Y,
      X          = X_array,
      sigma_beta = 1,
      sigma_rho  = 1
    )
    
    ## ----------------- Modelo Lcorr ----------------- ##
    if ("Lcorr" %in% fit_models) {
      
      begin <- Sys.time()
      print(begin)
      
      samp_Lcorr <- stan(
        data    = data_list,
        file    = file.path(programas_dir, "FitMultProbit_Lcorr.stan"),
        init    = ini,
        chains  = nChains,
        pars    = params_Lcorr,
        iter    = nIter,
        warmup  = burnInSteps,
        thin    = thinSteps,
        control = list(adapt_delta = 0.8, max_treedepth = 10),
        save_dso = TRUE
      )
      
      end <- Sys.time()
      print(end - begin)
      
      est_Lcorr <- summary(samp_Lcorr)
      Lcorr_file <- file.path(
        pathResults,
        paste0("summaryMult_Lcorr_R", r, "_S", j, ".rds")
      )
      saveRDS(est_Lcorr$summary, Lcorr_file)
    }
    
    ## ----------------- Modelo AR1 ----------------- ##
    if ("AR1" %in% fit_models) {
      
      begin <- Sys.time()
      print(begin)
      
      samp_AR1 <- stan(
        data    = data_list,
        file    = file.path(programas_dir, "FitMultProbit_AR1.stan"),
        init    = ini,
        chains  = nChains,
        pars    = params_AR1,
        iter    = nIter,
        warmup  = burnInSteps,
        thin    = thinSteps,
        control = list(adapt_delta = 0.8, max_treedepth = 10),
        save_dso = TRUE
      )
      
      end <- Sys.time()
      print(end - begin)
      
      est_AR1 <- summary(samp_AR1)
      AR1_file <- file.path(
        pathResults,
        paste0("summaryMult_AR1_R", r, "_S", j, ".rds")
      )
      saveRDS(est_AR1$summary, AR1_file)
    }
    
    ## ----------------- Modelo Independente ----------------- ##
    if ("Ind" %in% fit_models) {
      
      begin <- Sys.time()
      print(begin)
      
      # Supondo que o modelo independente só tenha beta (sem parâmetros de correlação)
      samp_Ind <- stan(
        data    = data_list,
        file    = file.path(programas_dir, "FitProbit_Ind.stan"),
        init    = ini,
        chains  = nChains,
        pars    = params_Ind,
        iter    = nIter,
        warmup  = burnInSteps,
        thin    = thinSteps,
        control = list(adapt_delta = 0.8, max_treedepth = 10),
        save_dso = TRUE
      )
      
      end <- Sys.time()
      print(end - begin)
      
      est_Ind <- summary(samp_Ind)
      Ind_file <- file.path(
        pathResults,
        paste0("summaryMult_Ind_R", r, "_S", j, ".rds")
      )
      saveRDS(est_Ind$summary, Ind_file)
    }
    
    message("Replication ", r, " completed.")
  }
  
  message("Scenario ", j, " completed.")
}
