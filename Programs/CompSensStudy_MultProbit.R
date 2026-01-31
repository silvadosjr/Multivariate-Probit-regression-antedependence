library(rstan)
library(rstansim)
library(ggplot2)
library(dplyr)
library(tidyr)

setwd('C:/Users/Usuário/OneDrive/Documentos/Artigos/Multivariate Probit regression using antedependence/')

library(here)

options(mc.cores = parallel::detectCores())

# Escolha da estrutura de correlação
corr_type <- "AR1"  # ou "Toeplitz"

## Correlation ranges
#[.10, .40] Weak; [.40, .70] Moderate; [.70, .9] Strong

# Caminhos organizados com `here()`
programas_dir  <- here("Programas")
pathResults    <- here("SensStudy")

# Funções auxiliares
source(file.path(programas_dir, "Aux_Functions.R"))


R <- 20  # número de réplicas


# Inicializa lista para armazenar dados de todas as métricas
res_list <- list()

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

for (k in c(2,3,4)){
  for (j in seq_len(nrow(Scenario))) {
    
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
    
    nome_cenario <- paste0("T", vT, "_n", n, "_", corr_label)
    
    for (r in seq_len(R)) {
      
      nome_arq <- paste0("R", r, "_S", j, "_P", k)
      
      file_AR1 <- file.path(pathResults, paste0("summaryMult_AR1_",nome_arq, ".rds"))
      
      # Leitura dos resumos
      if (file.exists(file_AR1)) {
        sum_AR1   <- readRDS(file_AR1)
        
        beta_AR1 <- sum_AR1[paste0('beta[',1:p,']'), 'mean']
        
        rho_AR1 <- sum_AR1['rho', 'mean']
        
         
        # Adiciona à lista os valores para cada métrica
        res_list[[length(res_list) + 1]] <- data.frame(
          Rep = r, Scenario = nome_cenario, Parameter = paste0('beta',1:p),
          Prior=paste0('Prior',k),
          ESS = sum_AR1[paste0('beta[',1:p,']'), "n_eff"],
          Rhat = sum_AR1[paste0('beta[',1:p,']'), "Rhat"],
          Bias = beta_AR1-beta_verd
        )
        res_list[[length(res_list) + 1]] <- data.frame(
          Rep = r, Scenario = nome_cenario, Parameter = 'rho',
          Prior=paste0('Prior',k),
          ESS = sum_AR1['rho', "n_eff"],
          Rhat = sum_AR1['rho', "Rhat"],
          Bias = rho_AR1-rho_verd
        )
      }
    
    }
}
    
}    

# Junta todos os resultados em um único dataframe
df_wide <- do.call(rbind, res_list)

df_wide$Prior[df_wide$Prior=='Prior4']<-'Prior1'

# Transforma em formato longo para ggplot
df_long <- pivot_longer(df_wide, cols = c("ESS", "Rhat", "Bias"),
                        names_to = "Metric", values_to = "Value")    
    
    
    
# # 📊 Boxplot para ESS (n_eff)
# ggplot(filter(df_long, Metric == "ESS"),
#        aes(x = Scenario, y = Value, fill = Parameter)) +
#   geom_boxplot(position = position_dodge(0.8), outlier.size = 0.8) +
#   labs(title = "Effective Sample Size",
#        x = "Scenario", y = "ESS", fill = "Parameter") +
#   theme_minimal(base_size = 13) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
#     
# 
# # 📊 Boxplot para PRSF (Rhat)
# ggplot(filter(df_long, Metric == "Rhat"),
#        aes(x = Scenario, y = Value, fill = Parameter)) +
#   geom_boxplot(position = position_dodge(0.8), outlier.size = 0.8) +
#   labs(title = "Potential Reduction Scale Factor",
#        x = "Scenario", y = "Rhat", fill = "Parameter") +
#   theme_minimal(base_size = 13) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# 
# # 📊 Boxplot para Abias 
# ggplot(filter(df_long, Metric == "Abias"),
#        aes(x = Scenario, y = Value, fill = Parameter)) +
#   geom_boxplot(position = position_dodge(0.8), outlier.size = 0.8) +
#   labs(title = "Absolute bias",
#        x = "Scenario", y = "Abias", fill = "Parameter") +
#   theme_minimal(base_size = 13) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))




## Resultado agregado

df_long2 <- df_long %>%
  mutate(ParameterGroup = ifelse(grepl("beta", Parameter), "beta", "rho"))

## Para beta

ggplot(filter(df_long2, Metric == "ESS", ParameterGroup == "beta"),
       aes(x = Scenario, y = Value, fill = Prior)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_grey(start = 0.3, end = 0.8) +
  labs(title = "Effective Sample Size (ESS)",
       x = "Scenario", y = "ESS", fill = "Prior") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(filter(df_long2, Metric == "Rhat", ParameterGroup == "beta"),
       aes(x = Scenario, y = Value, fill = Prior)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_hline(yintercept = 1.1, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_grey(start = 0.3, end = 0.8) +
  labs(title = "Potential Scale Reduction Factor (Rhat)",
       x = "Scenario", y = "Rhat") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(filter(df_long2, Metric == "Bias", ParameterGroup == "beta"),
       aes(x = Scenario, y = Value, fill = Prior)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_grey(start = 0.3, end = 0.8) +
  labs(title = "Bias",
       x = "Scenario", y = "Bias") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## Para rho


ggplot(filter(df_long2, Metric == "ESS", ParameterGroup == "rho"),
       aes(x = Scenario, y = Value, fill = Prior)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_grey(start = 0.3, end = 0.8) +
  labs(title = "Effective Sample Size (ESS)",
       x = "Scenario", y = "ESS") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(filter(df_long2, Metric == "Rhat", ParameterGroup == "rho"),
       aes(x = Scenario, y = Value, fill = Prior)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_hline(yintercept = 1.1, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_grey(start = 0.3, end = 0.8) +
  labs(title = "Potential Scale Reduction Factor (Rhat)",
       x = "Scenario", y = "Rhat") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




ggplot(filter(df_long2, Metric == "Bias", ParameterGroup == "rho"),
       aes(x = Scenario, y = Value, fill = Prior)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_grey(start = 0.3, end = 0.8) +
  labs(title = "Bias",
       x = "Scenario", y = "Bias") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



tab_aux<-df_long %>% filter(Scenario=='T4_n50_Strong',Metric=='Bias',Prior=='Prior4')

tab_aux[tab_aux$Value>2.5,]



## Tabela com resumos por parâmetro e priori


# Dicionário com valores verdadeiros de rho para cada cenário
rho_truth <- c(
  "T4_n50_Strong"   = 0.9,
  "T8_n50_Null"     = 0,
  "T12_n50_Moderate"= 0.7,
  "T4_n100_Moderate"= 0.7,
  "T8_n100_Strong"  = 0.9,
  "T12_n100_Null"   = 0,
  "T4_n200_Null"    = 0,
  "T8_n200_Moderate"= 0.7,
  "T12_n200_Strong" = 0.9
)

# Atribui valor verdadeiro (True) para cada parâmetro
df_corrigido <- df_wide %>%
  mutate(True = case_when(
    Parameter == "beta1" ~ 1,
    Parameter == "beta2" ~ 0.5,
    Parameter == "beta3" ~ 0.8,
    Parameter == "beta4" ~ 0.6,
    Parameter == "rho"   ~ rho_truth[Scenario],
    TRUE ~ NA_real_
  ),
  Estimativa = Bias + True)

# Calcula ABias e RMSE com base nas estimativas e verdadeiros valores
df_summary <- df_corrigido %>%
  group_by(Prior, Parameter, True) %>%
  summarise(
    Estimate=round(mean(Estimativa),3),
    ABias = round(mean(abs(Estimativa - True), na.rm = TRUE),3),
    RMSE  = round(sqrt(mean((Estimativa - True)^2, na.rm = TRUE)),3),
    .groups = 'drop'
  ) %>%
  arrange(Parameter, Prior)




## Salvando figuras

# Cria o diretório se não existir
dir.create("Figuras", showWarnings = FALSE)


## === Para beta === ##

# ESS - beta
g1 <- ggplot(filter(df_long2, Metric == "ESS", ParameterGroup == "beta"),
             aes(x = Scenario, y = Value, fill = Prior)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_grey(start = 0.3, end = 0.8) +
  labs(title = "Effective Sample Size (ESS)",
       x = "Scenario", y = "ESS", fill = "Prior") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Figuras/ESS_beta.eps", g1, device = "eps", width = 8, height = 5)

# Rhat - beta
g2 <- ggplot(filter(df_long2, Metric == "Rhat", ParameterGroup == "beta"),
             aes(x = Scenario, y = Value, fill = Prior)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_hline(yintercept = 1.1, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_grey(start = 0.3, end = 0.8) +
  labs(title = "Potential Scale Reduction Factor (Rhat)",
       x = "Scenario", y = "Rhat") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Figuras/Rhat_beta.eps", g2, device = "eps", width = 8, height = 5)

# Bias - beta
g3 <- ggplot(filter(df_long2, Metric == "Bias", ParameterGroup == "beta"),
             aes(x = Scenario, y = Value, fill = Prior)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_grey(start = 0.3, end = 0.8) +
  labs(title = "Bias",
       x = "Scenario", y = "Bias") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Figuras/Bias_beta.eps", g3, device = "eps", width = 8, height = 5)

## === Para rho === ##

# ESS - rho
g4 <- ggplot(filter(df_long2, Metric == "ESS", ParameterGroup == "rho"),
             aes(x = Scenario, y = Value, fill = Prior)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_grey(start = 0.3, end = 0.8) +
  labs(title = "Effective Sample Size (ESS)",
       x = "Scenario", y = "ESS") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Figuras/ESS_rho.eps", g4, device = "eps", width = 8, height = 5)

# Rhat - rho
g5 <- ggplot(filter(df_long2, Metric == "Rhat", ParameterGroup == "rho"),
             aes(x = Scenario, y = Value, fill = Prior)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_hline(yintercept = 1.1, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_grey(start = 0.3, end = 0.8) +
  labs(title = "Potential Scale Reduction Factor (Rhat)",
       x = "Scenario", y = "Rhat") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Figuras/Rhat_rho.eps", g5, device = "eps", width = 8, height = 5)

# Bias - rho
g6 <- ggplot(filter(df_long2, Metric == "Bias", ParameterGroup == "rho"),
             aes(x = Scenario, y = Value, fill = Prior)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_grey(start = 0.3, end = 0.8) +
  labs(title = "Bias",
       x = "Scenario", y = "Bias") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Figuras/Bias_rho.eps", g6, device = "eps", width = 8, height = 5)













