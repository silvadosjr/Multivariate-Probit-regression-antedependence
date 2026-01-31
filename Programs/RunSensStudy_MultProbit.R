library(rstan)
library(rstansim)

setwd('C:/Users/Usuário/OneDrive/Documentos/Artigos/Multivariate Probit regression using antedependence/')

library(here)

options(mc.cores = parallel::detectCores())

# Escolha da estrutura de correlação
corr_type <- "AR1"  # ou "Toeplitz"

## Correlation ranges
#[.40, .70] Moderate; [.70, .9] Strong or Null

# Caminhos organizados com `here()`
programas_dir  <- here("Programas")
pathDataSets   <- here("LatinSquareScenarios", corr_type)
pathResults    <- here("SensStudy")

# Funções auxiliares
source(file.path(programas_dir, "Aux_Functions.R"))

# Priors
PriorSettings <- data.frame(
  Prior = paste0("Prior", 1:4),
  sigma_beta = c(1, 7, 10, 1),
  sigma_rho = c(1, 1, 10, .5)
)


# Cenários de simulação
Scenario <- rbind(
  c(4,50,.9,.7), c(8,50,0,0), c(12,50,.7,.4),
  c(4,100,.7,.4), c(8,100,.9,.7), c(12,100,0,0),
  c(4,200,0,0), c(8,200,.7,.4), c(12,200,.9,.7)
)

R <- 20  # número de réplicas

for (k in c(3)) {
  for (j in c(2,6,7)) {
    
    vT         <- Scenario[j, 1]
    n          <- Scenario[j, 2]
    rho_start  <- Scenario[j, 3]
    rho_end    <- Scenario[j, 4]
    
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
    
    rho_verd <- if (corr_type == "AR1") {
      rho_start
    } else if (corr_type == "Toeplitz") {
      round(seq(rho_start, to = rho_end, length.out = vT - 1), 2)
    } else {
      stop("corr_type must be either 'AR1' or 'Toeplitz'")
    }
    
    corr_label <- if (rho_start == 0) {
      "Null"
    } else if (rho_start == 0.7) {
      "Moderate"
    } else {
      "Strong"
    }
    
    for (r in seq_len(R)) {
      
      # Carregar dados simulados
      rds_path <- file.path(pathDataSets, paste0("SimMultProbit_", corr_type, "_", vT, "_", n, "_", corr_label, "_", r, ".rds"))
      sim_AR1  <- readRDS(rds_path)
      
      # Preparar lista de dados
      data_list <- list(
        vT = vT, p = p, N = n,
        Y = sim_AR1$Y, X = X_array,
        sigma_beta = PriorSettings$sigma_beta[k],
        sigma_rho  = PriorSettings$sigma_rho[k]
      )
      
      ini <- function() list(beta = rep(0.1, p), rho = 0.1)
      params_AR1 <- c("beta", "rho")
      
      nChains <- 1
      burnInSteps <- 1000
      thinSteps <- 15
      numSavedSteps <- 1000
      nIter <- ceiling(burnInSteps + (numSavedSteps * thinSteps) / nChains)
      
      begin <- Sys.time(); print(begin)
      
      samp_AR1 <- stan(
        data = data_list,
        file = file.path(programas_dir, "FitMultProbit_AR1.stan"),
        init = ini,
        chains = nChains,
        pars = params_AR1,
        iter = nIter,
        warmup = burnInSteps,
        thin = thinSteps,
        control = list(adapt_delta = 0.8, max_treedepth = 10),
        save_dso = TRUE
      )
      
      end <- Sys.time()
      print(end - begin)
      
      est_AR1 <- summary(samp_AR1)
      
      out_file <- file.path(pathResults, paste0("summaryMult_AR1_R", r, "_S", j, "_P", k, ".rds"))
      saveRDS(est_AR1$summary, out_file)
      
      message("Replication ", r, " completed.")
    }
    
    message("Scenario ", j, " completed.")
  }
  
  message("Prior ", k, " completed.")
}




