library(rstan)
library(rstansim)

setwd('C:/Users/Usuário/OneDrive/Documentos/Artigos/Multivariate Probit regression using antedependence/')

library(here)

options(mc.cores = parallel::detectCores())

# Define o tipo de correlação: "AR1" ou "Toeplitz"
corr_type <- "AR1"  

## Correlation ranges
#[.10, .40] Weak; [.40, .70] Moderate; [.70, .9] Strong or Null

# Define caminho para os scripts e resultados
programas_dir <- here("Programas")
pathSim <- here("LatinSquareScenarios", corr_type)

# Carrega funções auxiliares
source(file.path(programas_dir, "Aux_Functions.R"))

# Define cenários
Scenario <- rbind(
  c(4,50,.9,.7), c(8,50,0,0), c(12,50,.7,.4),
  c(4,100,.7,.4), c(8,100,.9,.7), c(12,100,0,0),
  c(4,200,0,0), c(8,200,.7,.4), c(12,200,.9,.7)
)

R <- 20

for (j in c(2,6,7)) {
  
  # Extração de parâmetros do cenário
  vT <- Scenario[j, 1]
  n <- Scenario[j, 2]
  rho_start <- Scenario[j, 3]
  rho_end <- Scenario[j, 4]
  
  # Construção da matriz de planejamento
  t_age <- 1:vT - mean(1:vT)
  X_s <- cbind(1, t_age, 1, t_age)
  X_ns <- cbind(1, t_age, 0, 0)
  
  treat_counts <- c(n / 2, n / 2)
  p <- 4
  
  X_array <- array(NA, dim = c(length(t_age), p, sum(treat_counts)))
  X_array[, , 1:treat_counts[1]] <- replicate(treat_counts[1], X_ns, simplify = "array")
  X_array[, , (treat_counts[1] + 1):sum(treat_counts)] <- replicate(treat_counts[2], X_s, simplify = "array")
  X_array <- aperm(X_array, c(3, 1, 2))
  
  # Parâmetros verdadeiros
  beta_verd <- c(1, .5, .8, .6)
  
  if (corr_type == "AR1") {
    rho_verd <- rho_start
  } else if (corr_type == "Toeplitz") {
    rho_verd <- round(seq(rho_start, to = rho_end, length.out = vT - 1), 2)
  } else {
    stop("corr_type must be either 'AR1' or 'Toeplitz'")
  }
  
  # Classificação da correlação
  corr_label <- if (rho_start == 0) {
    "Null"
  } else if (rho_start == 0.7) {
    "Moderate"
  } else {
    "Strong"
  }
  
  # Dados e parâmetros
  data_simAux <- list(vT = vT, p = length(beta_verd), n = n, X = X_array)
  param_sim1 <- list(beta = beta_verd, rho = rho_verd)
  
  # Simulação
  begin <- Sys.time(); print(begin)
  
  sim_out <- simulate_data(
    file =paste0(programas_dir,ifelse(corr_type == "AR1", "/simMultProbit_AR1.stan", "/simMultProbit_HT.stan")),
    data_name = paste0("SimMultProbit_", corr_type, "_", vT, "_", n, "_", corr_label),
    input_data = data_simAux,
    param_values = param_sim1,
    vars = c("Y"),
    nsim = R,
    path = pathSim
  )
  
  end <- Sys.time()
  print(end - begin)
}
