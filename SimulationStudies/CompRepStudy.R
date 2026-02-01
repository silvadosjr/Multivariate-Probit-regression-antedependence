###############################################################################
# Multivariate Probit Regression with Structured Dependence
# Robustness Study (Post-processing + Figures + LaTeX Table)
#
# Purpose
# -------
# This script post-processes the posterior summaries produced in the robustness
# simulation study, where data are generated under a given correlation structure
# (here: Toeplitz) and then fitted with competing models (e.g., AR1, Lcorr,
# Independence).
#
# For each scenario and replication, it:
#   1) Reads saved Stan summary tables (.rds) for each fitted model.
#   2) Reconstructs the fitted correlation matrix (C) when applicable.
#   3) Computes a Stein-type loss between the fitted correlation matrix and
#      the true correlation matrix used to generate the data.
#   4) Extracts diagnostics (ESS, Rhat) from lp__ and from each beta parameter.
#   5) Computes beta accuracy metrics (Bias, ABias, RMSE) replication-by-replication.
#   6) Saves consolidated data frames for (i) global metrics and (ii) betas.
#   7) Produces and saves boxplots (ESS, Rhat, Stein) and beta-specific boxplots
#      (ESS, Rhat, Bias, RMSE, SD).
#   8) Builds a LaTeX table (xtable) summarizing mean Estimate/ABias/RMSE by
#      parameter and fitted model.
#
# Inputs
# ------
# - Posterior summary files from fitted models, stored as:
#     SimulationStudies/ResultsSimRobus/summaryMult_<Model>_R<r>_S<j>.rds
#   where Model ∈ {Lcorr, AR1, Ind}, r ∈ {1..R}, j ∈ {1..nrow(Scenario)}.
#
# - Auxiliary functions from:
#     Programs/Aux_Functions.R
#   Expected helpers include:
#     - toeplitz.matrix(...)
#     - ARH1.matrix(...)
#     - reconstruir_correlacao(...)   (reconstruct correlation matrix from Lcorr fit)
#
# Outputs
# -------
# - Two consolidated objects saved as .rds:
#     ResultadosSimRobus/metrics_Toeplitz.rds   (df_metrics)
#     ResultadosSimRobus/betas_Toeplitz.rds     (df_betas)
#
# - Figures saved as EPS into Figuras/:
#     Boxplot_ESS_Toeplitz.eps
#     Boxplot_Rhat_Toeplitz.eps
#     Boxplot_Stein_Toeplitz.eps
#     Boxplot_ESS_beta_Toeplitz.eps
#     Boxplot_Rhat_beta_Toeplitz.eps
#     Boxplot_Bias_Toeplitz.eps
#     Boxplot_RMSE_Toeplitz.eps
#     Boxplot_SD_Toeplitz.eps
#
# - A LaTeX table printed to console (xtable) summarizing beta recovery.
#
# Notes / Assumptions
# -------------------
# - Default settings use 1 chain. With a single chain, Rhat may be NA or not
#   meaningful. If you rely on Rhat, consider running >= 2 chains.
#
# - The "Stein" criterion used here is:
#       Stein(C_hat, C_true) = tr(C_hat * C_true^{-1}) - log det(C_hat * C_true^{-1}) - vT
#   This is a common loss for covariance/correlation matrix comparison.
#
###############################################################################

# Packages -----------------------------------------------------------------
library(rstan)
library(rstansim)
library(ggplot2)
library(dplyr)
library(tidyr)

options(mc.cores = parallel::detectCores())

setwd('~/GitHub/Multivariate-Probit-regression-antedependence/')

library(here)
source(here("Programs", "Aux_Functions.R"))

# ------------------------------------------------------------
# Control: which models will be processed?
# ------------------------------------------------------------
#process_models <- c("AR1", "Ind")
# process_models <- c("Lcorr","Ind")
# process_models <- c("AR1")
 process_models <- c("AR1",'Lcorr')

corr_type <- "Toeplitz"
pathResults   <- here("SimulationStudies", "ResultsSimRobus")

# True betas
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

# ------------------------------------------------------------
#  Main LOOP 
# ------------------------------------------------------------

for (j in 1:nrow(Scenario)) {
  
  vT        <- Scenario[j, 1]
  n         <- Scenario[j, 2]
  rho_start <- Scenario[j, 3]
  rho_end   <- Scenario[j, 4]
  
  rho_verd <- round(seq(rho_start, to = rho_end, length.out = vT - 1), 2)
  C_verd <- toeplitz.matrix(rep(1, vT), rho_verd)
  
  corr_label <- if (rho_start == 0.4) "Weak" else if (rho_start == 0.7) "Moderate" else "Strong"
  nome_cenario <- paste0("T", vT, "_n", n, "_", corr_label)
  
  for (r in 1:R) {
    
    # paths
    file_Lcorr <- file.path(pathResults, paste0("summaryMult_Lcorr_R", r, "_S", j, ".rds"))
    file_AR1   <- file.path(pathResults, paste0("summaryMult_AR1_R", r, "_S", j, ".rds"))
    file_Ind   <- file.path(pathResults, paste0("summaryMult_Ind_R", r, "_S", j, ".rds"))
    
    # ---------------------------
    # Lcorr
    # ---------------------------
    if ("Lcorr" %in% process_models && file.exists(file_Lcorr)) {
      
      sum_Lcorr <- readRDS(file_Lcorr)
      C_Lcorr <- reconstruir_correlacao(sum_Lcorr)
      
      stein_L <- sum(diag(C_Lcorr %*% solve(C_verd))) -
        log(det(C_Lcorr %*% solve(C_verd))) - vT
      
      res_list[[length(res_list)+1]] <- data.frame(
        Rep = r, Scenario = nome_cenario, Model = "Unstructured",
        ESS = sum_Lcorr["lp__", "n_eff"],
        Rhat = sum_Lcorr["lp__", "Rhat"],
        Stein = stein_L
      )
      
      for (i in 1:4) {
        est <- sum_Lcorr[param_names[i], "mean"]
        sdv <- sum_Lcorr[param_names[i], "sd"]
        ess <- sum_Lcorr[param_names[i], 'n_eff']
        Rhat <- sum_Lcorr[param_names[i], 'Rhat']
        
        beta_list[[length(beta_list)+1]] <- data.frame(
          Rep = r, Scenario = nome_cenario, Model = "Unstructured",
          Parameter = param_names[i], Estimate = est,
          ABias = abs(est - beta_verd[i]),
          Bias = est - beta_verd[i],
          SD = sdv,
          RMSE = sqrt((est - beta_verd[i])^2),
          ESS = ess,
          Rhat = Rhat
        )
      }
    }
    
    # ---------------------------
    # AR1
    # ---------------------------
    if ("AR1" %in% process_models && file.exists(file_AR1)) {
      
      sum_AR1 <- readRDS(file_AR1)
      C_AR1 <- ARH1.matrix(rep(1, vT), sum_AR1["rho", "mean"])
      
      stein_A <- sum(diag(C_AR1 %*% solve(C_verd))) -
        log(det(C_AR1 %*% solve(C_verd))) - vT
      
      res_list[[length(res_list)+1]] <- data.frame(
        Rep = r, Scenario = nome_cenario, Model = "AR1",
        ESS = sum_AR1["lp__", "n_eff"],
        Rhat = sum_AR1["lp__", "Rhat"],
        Stein = stein_A
      )
      
      for (i in 1:4) {
        est <- sum_AR1[param_names[i], "mean"]
        sdv <- sum_AR1[param_names[i], "sd"]
        ess <- sum_AR1[param_names[i], 'n_eff']
        Rhat <- sum_AR1[param_names[i], 'Rhat']
        
        beta_list[[length(beta_list)+1]] <- data.frame(
          Rep = r, Scenario = nome_cenario, Model = "AR1",
          Parameter = param_names[i], Estimate = est,
          ABias = abs(est - beta_verd[i]),
          Bias = est - beta_verd[i],
          SD = sdv,
          RMSE = sqrt((est - beta_verd[i])^2),
          ESS = ess,
          Rhat = Rhat
        )
      }
    }
    
    # ---------------------------
    # Independente
    # ---------------------------
    if ("Ind" %in% process_models && file.exists(file_Ind)) {
      
      sum_Ind <- readRDS(file_Ind)
      
      # C independente = Identity matrix
      C_Ind <- diag(vT)
      
      stein_I <- sum(diag(C_Ind %*% solve(C_verd))) -
        log(det(C_Ind %*% solve(C_verd))) - vT
      
      res_list[[length(res_list)+1]] <- data.frame(
        Rep = r, Scenario = nome_cenario, Model = "Ind",
        ESS = sum_Ind["lp__", "n_eff"],
        Rhat = sum_Ind["lp__", "Rhat"],
        Stein = stein_I
      )
      
      for (i in 1:4) {
        est <- sum_Ind[param_names[i], "mean"]
        sdv <- sum_Ind[param_names[i], "sd"]
        ess <- sum_Ind[param_names[i], 'n_eff']
        Rhat <- sum_Ind[param_names[i], 'Rhat']
        
        beta_list[[length(beta_list)+1]] <- data.frame(
          Rep = r, Scenario = nome_cenario, Model = "Ind",
          Parameter = param_names[i], Estimate = est,
          ABias = abs(est - beta_verd[i]),
          Bias = est - beta_verd[i],
          SD = sdv,
          RMSE = sqrt((est - beta_verd[i])^2),
          ESS = ess,
          Rhat = Rhat
        )
      }
    }
    
  }
}

# ------------------------------------------------------------
# Consolidation
# ------------------------------------------------------------

df_metrics <- do.call(rbind, res_list)
df_betas <- do.call(rbind, beta_list)

saveRDS(df_metrics, here("ResultadosSimRobus", paste0("metrics_", corr_type, ".rds")))
saveRDS(df_betas, here("ResultadosSimRobus", paste0("betas_", corr_type, ".rds")))

# ------------------------------------------------------------
# Plots
# ------------------------------------------------------------

df_long <- pivot_longer(df_metrics, cols = c("ESS", "Rhat", "Stein"),
                        names_to = "Metric", values_to = "Value")

df_long$Scenario <- factor(df_long$Scenario, levels = unique(df_long$Scenario))

cores_modelo <- c("Unstructured" = "gray40", "AR1" = "gray70", "Ind" = "white")
tema_base <- theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ESS
ess_plot<-ggplot(filter(df_long, Metric == "ESS"),
       aes(x = Scenario, y = Value, fill = Model)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  scale_fill_manual(values = cores_modelo) +
  labs(title = '',
       x = "Scenario", y = "ESS", fill = "Model") +
  tema_base
ggsave(here("Figuras", paste0("Boxplot_ESS", "_", corr_type, ".eps")),
       plot = ess_plot, width = 9, height = 6, device = "eps")

# Rhat
Rhat_plot<-ggplot(filter(df_long, Metric == "Rhat"),
       aes(x = Scenario, y = Value, fill = Model)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_hline(yintercept = 1.1, linetype = "dashed") +
  scale_fill_manual(values = cores_modelo) +
  labs(title = '',
       x = "Scenario", y = "Rhat", fill = "Model") +
  tema_base
ggsave(here("Figuras", paste0("Boxplot_Rhat", "_", corr_type, ".eps")),
       plot = Rhat_plot, width = 9, height = 6, device = "eps")


# Stein
Stein_plot<-ggplot(filter(df_long, Metric == "Stein"),
       aes(x = Scenario, y = Value, fill = Model)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values = cores_modelo) +
  labs(title = '',
       x = "Scenario", y = "Stein", fill = "Model") +
  tema_base
ggsave(here("Figuras", paste0("Boxplot_Stein", "_", corr_type, ".eps")),
       plot = Stein_plot, width = 9, height = 6, device = "eps")


# ------------------------------------------------------------
# Boxplots of Betas
# ------------------------------------------------------------


ESS_beta_plot<-ggplot(df_betas, aes(x = Scenario, y = ESS, fill = Model)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  scale_fill_manual(values = cores_modelo) +
  labs(x = "Scenario", y = "ESS") +
  tema_base
ggsave(here("Figuras", paste0("Boxplot_ESS_beta", "_", corr_type, ".eps")),
       plot = ESS_beta_plot, width = 9, height = 6, device = "eps")

Rhat_beta_plot<-ggplot(df_betas, aes(x = Scenario, y = Rhat, fill = Model)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_hline(yintercept = 1.1, linetype = "dashed") +
  scale_fill_manual(values = cores_modelo) +
  labs(x = "Scenario", y = "Rhat") +
  tema_base
ggsave(here("Figuras", paste0("Boxplot_Rhat_beta", "_", corr_type, ".eps")),
       plot = Rhat_beta_plot, width = 9, height = 6, device = "eps")



bias_plot<-ggplot(df_betas, aes(x = Scenario, y = Bias, fill = Model)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = cores_modelo) +
  labs(x = "Scenario", y = "Bias") +
  tema_base
ggsave(here("Figuras", paste0("Boxplot_Bias", "_", corr_type, ".eps")),
       plot = bias_plot, width = 9, height = 6, device = "eps")



RMSE_plot<-ggplot(df_betas, aes(x = Scenario, y = RMSE, fill = Model)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values = cores_modelo) +
  labs(x = "Scenario", y = "RMSE") +
  tema_base
ggsave(here("Figuras", paste0("Boxplot_RMSE", "_", corr_type, ".eps")),
       plot = RMSE_plot, width = 9, height = 6, device = "eps")


SD_plot<-ggplot(df_betas, aes(x = Scenario, y = SD, fill = Model)) +
  geom_boxplot(position = position_dodge(0.8)) +
  scale_fill_manual(values = cores_modelo) +
  labs(x = "Scenario", y = "SD") +
  tema_base
ggsave(here("Figuras", paste0("Boxplot_SD", "_", corr_type, ".eps")),
       plot = SD_plot, width = 9, height = 6, device = "eps")




library(knitr)     # for kable()
library(xtable)
# Loading beta's result
df_betas <- readRDS(here("ResultadosSimRobus", "betas_Toeplitz.rds"))

# Actual values (to ensure consistency)
beta_verd <- c("beta[1]" = 1.0, "beta[2]" = 0.5, "beta[3]" = 0.8, "beta[4]" = 0.6)

# Summary Table
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
# Present with xtable for LaTeX
xt <- xtable(tabela_latex, digits = c(0, 0, 0, 1, 3, 3, 3),
             caption = "Summary of average estimates, absolute bias (ABias), and RMSE by parameter and model.",
             label = "tab:summary_betas")

print(xt, include.rownames = FALSE, booktabs = TRUE, caption.placement = "top")


