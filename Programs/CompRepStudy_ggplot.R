library(rstan)
library(rstansim)
library(ggplot2)
library(dplyr)
library(tidyr)


options(mc.cores = parallel::detectCores())

PC <- 'Usuário'
setwd(paste0('C:/Users/', PC, '/OneDrive/Documentos/Artigos/Multivariate Probit regression using antedependence/'))


library(here)

source(here("Programas", "Aux_Functions.R"))

corr_type <- "Toeplitz"
pathResults <- here("ResultadosSimRobus")

# Verdadeiros betas
beta_verd <- c(1, .5, .8, .6)
param_names <- paste0("beta[", 1:4, "]")

Scenario <- rbind(
  c(4, 50, .9, .7), c(8, 50, .4, .1), c(12, 50, .7, .4),
  c(4, 100, .7, .4), c(8, 100, .9, .7), c(12, 100, .4, .1),
  c(4, 200, .4, .1), c(8, 200, .7, .4), c(12, 200, .9, .7)
)

R <- 20
res_list <- list()
beta_list <- list()

for (j in 1:nrow(Scenario)) {
  vT         <- Scenario[j, 1]
  n          <- Scenario[j, 2]
  rho_start  <- Scenario[j, 3]
  rho_end    <- Scenario[j, 4]
  
  rho_verd <- round(seq(rho_start, to = rho_end, length.out = vT - 1), 2)
  C_verd <- toeplitz.matrix(rep(1, vT), rho_verd)
  
  corr_label <- if (rho_start == 0.4) {
    "Weak"
  } else if (rho_start == 0.7) {
    "Moderate"
  } else {
    "Strong"
  }
  
  nome_cenario <- paste0("T", vT, "_n", n, "_", corr_label)
  
  for (r in 1:R) {
    
    file_Lcorr <- file.path(pathResults, paste0("summaryMult_Lcorr_R", r, "_S", j, ".rds"))
    file_AR1   <- file.path(pathResults, paste0("summaryMult_AR1_R", r, "_S", j, ".rds"))
    
    if (file.exists(file_Lcorr) && file.exists(file_AR1)) {
      
      sum_Lcorr <- readRDS(file_Lcorr)
      sum_AR1   <- readRDS(file_AR1)
      
      ## MÉTRICAS lp__, distância de Stein
      C_Lcorr <- reconstruir_correlacao(sum_Lcorr)
      C_AR1   <- ARH1.matrix(rep(1, vT), sum_AR1["rho", "mean"])
      
      stein_L <- sum(diag(C_Lcorr %*% solve(C_verd))) - log(det(C_Lcorr %*% solve(C_verd))) - vT
      stein_A <- sum(diag(C_AR1 %*% solve(C_verd))) - log(det(C_AR1 %*% solve(C_verd))) - vT
      
      res_list[[length(res_list) + 1]] <- data.frame(
        Rep = r, Scenario = nome_cenario, Model = "Lcorr",
        ESS = sum_Lcorr["lp__", "n_eff"],
        Rhat = sum_Lcorr["lp__", "Rhat"],
        Stein = stein_L
      )
      res_list[[length(res_list) + 1]] <- data.frame(
        Rep = r, Scenario = nome_cenario, Model = "AR1",
        ESS = sum_AR1["lp__", "n_eff"],
        Rhat = sum_AR1["lp__", "Rhat"],
        Stein = stein_A
      )
      
      ## MÉTRICAS dos BETAS
      for (i in 1:4) {
        est <- sum_Lcorr[param_names[i], "mean"]
        sd_L <- sum_Lcorr[param_names[i], "sd"]
        beta_list[[length(beta_list) + 1]] <- data.frame(
          Rep = r,
          Scenario = nome_cenario,
          Model = "Lcorr",
          Parameter = param_names[i],
          Estimate = est,
          ABias = abs(est - beta_verd[i]),
          ViesR=(est - beta_verd[i])/beta_verd[i],
          SD=sd_L,
          RMSE = sqrt((est - beta_verd[i])^2)
        )
        est <- sum_AR1[param_names[i], "mean"]
        sd_A <- sum_AR1[param_names[i], "sd"]
        beta_list[[length(beta_list) + 1]] <- data.frame(
          Rep = r,
          Scenario = nome_cenario,
          Model = "AR1",
          Parameter = param_names[i],
          Estimate = est,
          ABias = abs(est - beta_verd[i]),
          ViesR=(est - beta_verd[i])/beta_verd[i],
          SD=sd_A,
          RMSE = sqrt((est - beta_verd[i])^2)
        )
      }
    }
  }
}

# Combina os data.frames
df_metrics <- do.call(rbind, res_list)
df_betas <- do.call(rbind, beta_list)

# Salva dados para uso posterior
saveRDS(df_metrics, here("ResultadosSimRobus", paste0("metrics_", corr_type, ".rds")))
saveRDS(df_betas, here("ResultadosSimRobus", paste0("betas_", corr_type, ".rds")))

# Gráficos (.eps)
df_long <- pivot_longer(df_metrics, cols = c("ESS", "Rhat", "Stein"),
                        names_to = "Metric", values_to = "Value")

df_long$Scenario <- factor(df_long$Scenario, levels = unique(df_long$Scenario))

cores_modelo <- c("Lcorr" = "gray40", "AR1" = "gray70")
tema_base <- theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 📈 Boxplot geral com facet
ggplot(df_long, aes(x = Scenario, y = Value, fill = Model)) +
  geom_boxplot(position = position_dodge(0.8), outlier.size = 0.8) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 1) +
  labs(x = "Cenário", y = "Valor da Métrica", fill = "Modelo") +
  tema_base +
  scale_fill_manual(values = cores_modelo)

# 📊 ESS
ggplot(filter(df_long, Metric == "ESS"),
       aes(x = Scenario, y = Value, fill = Model)) +
  geom_boxplot(position = position_dodge(0.8), outlier.size = 0.8) +
  scale_fill_manual(values = cores_modelo) +
  labs(title = "Tamanho Efetivo da Amostra (lp__)",
       x = "Cenário", y = "ESS", fill = "Modelo") +
  tema_base

# 📊 Rhat
ggplot(filter(df_long, Metric == "Rhat"),
       aes(x = Scenario, y = Value, fill = Model)) +
  geom_boxplot(position = position_dodge(0.8), outlier.size = 0.8) +
  geom_hline(yintercept = 1.1, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_manual(values = cores_modelo) +
  labs(title = "Fator de Redução Potencial (R̂) do lp__",
       x = "Cenário", y = "R̂", fill = "Modelo") +
  tema_base

# 📊 Distância de Stein
ggplot(filter(df_long, Metric == "Stein"),
       aes(x = Scenario, y = Value, fill = Model)) +
  geom_boxplot(position = position_dodge(0.8), outlier.size = 0.8) +
  scale_fill_manual(values = cores_modelo) +
  labs(title = "Distância de Stein entre Matrizes de Correlação",
       x = "Cenário", y = "Distância de Stein", fill = "Modelo") +
  tema_base


## Salvar figuras - gráficos

metrics <- unique(df_long$Metric)
for (m in metrics) {
  p <- ggplot(filter(df_long, Metric == m), aes(x = Scenario, y = Value, fill = Model)) +
    geom_boxplot(position = position_dodge(0.8), outlier.size = 0.8) +
    scale_fill_manual(values = cores_modelo) +
    labs(title = m, x = "Cenário", y = m, fill = "Modelo") +
    tema_base
  if (m == "Rhat") p <- p + geom_hline(yintercept = 1.1, linetype = "dashed", color = "black")
  ggsave(here("Figuras", paste0("Boxplot_", m, "_", corr_type, ".eps")),
         plot = p, width = 9, height = 6, device = "eps")
}


library(knitr)     # para kable()
library(xtable)
# Carregar resultados dos betas
df_betas <- readRDS(here("ResultadosSimRobus", "betas_Toeplitz.rds"))

# Verdadeiros valores (para garantir consistência)
beta_verd <- c("beta[1]" = 1.0, "beta[2]" = 0.5, "beta[3]" = 0.8, "beta[4]" = 0.6)

# Tabela no formato desejado
tabela_latex <- df_betas %>%
  group_by(Parameter, Model) %>%
  summarise(
    `True Value` = unique(beta_verd[Parameter]),
    Estimate = mean(Estimate),
    ABias = mean(ABias),
    RMSE = mean(RMSE),
    .groups = "drop"
  ) %>%
  arrange(Parameter, Model)
# Exibir com xtable para LaTeX
xt <- xtable(tabela_latex, digits = c(0, 0, 0, 1, 3, 3, 3),
             caption = "Summary of average estimates, absolute bias (ABias), and RMSE by parameter and model.",
             label = "tab:summary_betas")

print(xt, include.rownames = FALSE, booktabs = TRUE, caption.placement = "top")




## Gráficos betas

# Carregar resultados
df_betas <- readRDS(here("ResultadosSimRobus", "betas_Toeplitz.rds"))

# Adiciona coluna com o viés bruto
df_betas <- df_betas %>%
  mutate(Bias = Estimate - case_when(
    Parameter == "beta[1]" ~ 1,
    Parameter == "beta[2]" ~ 0.5,
    Parameter == "beta[3]" ~ 0.8,
    Parameter == "beta[4]" ~ 0.6
  ))

# 🎨 Tema e cores
tema_base <- theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
cores_modelo <- c("Lcorr" = "gray40", "AR1" = "gray70")

# 📊 Boxplot para Viés (bruto) dos betas
ggplot(df_betas, aes(x = Scenario, y = Bias, fill = Model)) +
  geom_boxplot(position = position_dodge(0.8), outlier.size = 0.7) +
  scale_fill_manual(values = cores_modelo) +
  labs(title = "",
       x = "Cenário", y = "Viés", fill = "Modelo") +
  tema_base

# 📊 Boxplot para RMSE dos betas
ggplot(df_betas, aes(x = Scenario, y = RMSE, fill = Model)) +
  geom_boxplot(position = position_dodge(0.8), outlier.size = 0.7) +
  scale_fill_manual(values = cores_modelo) +
  labs(title = "",
       x = "Cenário", y = "RMSE", fill = "Modelo") +
  tema_base



# Boxplot do SD a posteriori
ggplot(df_betas, aes(x = Scenario, y = SD, fill = Model)) +
  geom_boxplot(position = position_dodge(0.8), outlier.size = 0.7) +
  scale_fill_manual(values = c("Lcorr" = "gray40", "AR1" = "gray70")) +
  labs(title = "Desvio Padrão a Posteriori dos Coeficientes Beta",
       x = "Cenário", y = "SD", fill = "Modelo") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))







