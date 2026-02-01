###############################################################################
# Multivariate Probit Regression with Structured Dependence
# Sensitivity Study: Post-processing, Diagnostics, and Figures
#
# Purpose
# -------
# This script aggregates posterior summaries produced by Stan fits in the
# sensitivity study (different priors) and computes diagnostic and accuracy
# metrics across scenarios and replications. In particular, it:
#
#   1) Reads per-fit summary tables saved as .rds files (one per replication).
#   2) Extracts posterior means, effective sample size (n_eff), and Rhat for:
#        - Regression coefficients: beta[1], ..., beta[p]
#        - Correlation parameter:  rho   (AR1 case)
#   3) Computes estimation error relative to known "true" values used in the
#      simulation (Bias, and later ABias and RMSE).
#   4) Produces boxplots (ESS, Rhat, Bias) by scenario and prior, separately
#      for beta and rho.
#   5) Saves figures to disk (EPS format) under the folder "Figures/".
#
# Inputs
# ------
# - Summary files produced previously by the fitting script, expected at:
#     SimulationStudies/SensStudy/summaryMult_AR1_R<r>_S<j>_P<k>.rds
#   where:
#     r = replication index (1..R)
#     j = scenario index (1..nrow(Scenario))
#     k = prior index (as used in the fitting loop)
#
# - The scenario grid (Scenario) and true parameters (beta_verd, rho_truth)
#   are defined inside this script so that errors can be computed.
#
# Outputs
# -------
# - Data frames in memory:
#     df_wide    : long-by-parameter table with ESS, Rhat, Bias
#     df_long    : pivoted long format for plotting
#     df_summary : aggregated table with posterior mean Estimate, ABias, RMSE
#
# - Figures (EPS) saved to:
#     Figuras/ESS_beta.eps
#     Figuras/Rhat_beta.eps
#     Figuras/Bias_beta.eps
#     Figuras/ESS_rho.eps
#     Figuras/Rhat_rho.eps
#     Figuras/Bias_rho.eps
#
# Notes / Assumptions
# -------------------
# - This script is currently configured for the AR(1) model output where
#   there is a single correlation parameter named 'rho' in the Stan summary.
# - For Toeplitz, the correlation parameter would typically be a vector
#   (e.g., rho[1:(vT-1)]), so the extraction and plotting blocks would need
#   to be adapted accordingly.
# - The horizontal lines in ESS (100) and Rhat (1.1) plots are common
#   diagnostic thresholds used to flag potential MCMC issues.
#
###############################################################################

# Packages -----------------------------------------------------------------

library(rstan)
library(rstansim)
library(ggplot2)
library(dplyr)
library(tidyr)

# Project root -------------------------------------------------------------

setwd('~/GitHub/Multivariate-Probit-regression-antedependence/')
library(here)

# Use all available CPU cores for Stan
options(mc.cores = parallel::detectCores())

# Correlation structure ----------------------------------------------------
corr_type <- "AR1"  # or "Toeplitz"

# Correlation ranges (for labeling scenarios)
# [.40, .70] Moderate; [.70, .90] Strong and Null

# Paths --------------------------------------------------------------------

programas_dir <- here("Programs")

# Folder containing posterior summaries from the sensitivity study fits
pathResults <- here("SimulationStudies", "SensStudy")

# Load auxiliary functions (kept for workflow consistency) -----------------
source(file.path(programas_dir, "Aux_Functions.R"))

# Simulation design parameters --------------------------------------------

R <- 20  # number of replications

# Container for results across all scenarios/priors/replications
res_list <- list()

# Prior configurations used during model fitting ---------------------------

PriorSettings <- data.frame(
  Prior = paste0("Prior", 1:4),
  sigma_beta = c(1, 7, 10, 1),
  sigma_rho  = c(1, 1, 10, .5)
)

# Scenario grid ------------------------------------------------------------
# Each row: (vT, n, rho_start, rho_end)
Scenario <- rbind(
  c(4,  50, .9, .7), c(8,  50, 0,  0), c(12,  50, .7, .4),
  c(4, 100, .7, .4), c(8, 100, .9, .7), c(12, 100, 0,  0),
  c(4, 200, 0,  0),  c(8, 200, .7, .4), c(12, 200, .9, .7)
)

# Main loops ---------------------------------------------------------------
# k : prior index (here using Prior2, Prior3, Prior4 from the stored files)
# j : scenario index
# r : replication index
for (k in c(2, 3, 4)) {
  
  for (j in seq_len(nrow(Scenario))) {
    
    # Extract scenario configuration
    vT        <- Scenario[j, 1]
    n         <- Scenario[j, 2]
    rho_start <- Scenario[j, 3]
    rho_end   <- Scenario[j, 4]
    
    # Build the design array X (same construction as in the fitting script)
    t_age <- 1:vT - mean(1:vT)
    X_s   <- cbind(1, t_age, 1, t_age)
    X_ns  <- cbind(1, t_age, 0, 0)
    
    treat_counts <- rep(n / 2, 2)
    p <- 4
    
    X_array <- array(NA, dim = c(length(t_age), p, n))
    X_array[, , 1:treat_counts[1]] <- replicate(treat_counts[1], X_ns, simplify = "array")
    X_array[, , (treat_counts[1] + 1):n] <- replicate(treat_counts[2], X_s, simplify = "array")
    X_array <- aperm(X_array, c(3, 1, 2))
    
    # True regression parameters (used to compute Bias / ABias / RMSE)
    beta_verd <- c(1, .5, .8, .6)
    
    # True correlation parameter(s) used in simulation (AR1 scalar)
    rho_verd <- if (corr_type == "AR1") {
      rho_start
    } else if (corr_type == "Toeplitz") {
      round(seq(rho_start, to = rho_end, length.out = vT - 1), 2)
    } else {
      stop("corr_type must be either 'AR1' or 'Toeplitz'")
    }
    
    # Scenario label used for grouping on plots and in file naming
    corr_label <- if (rho_start == 0) {
      "Null"
    } else if (rho_start == 0.7) {
      "Moderate"
    } else {
      "Strong"
    }
    
    nome_cenario <- paste0("T", vT, "_n", n, "_", corr_label)
    
    # Replication loop -----------------------------------------------------
    for (r in seq_len(R)) {
      
      # File naming pattern used by the fitting script
      nome_arq <- paste0("R", r, "_S", j, "_P", k)
      
      file_AR1 <- file.path(pathResults, paste0("summaryMult_AR1_", nome_arq, ".rds"))
      
      # Read posterior summary if it exists
      if (file.exists(file_AR1)) {
        
        sum_AR1 <- readRDS(file_AR1)
        
        # Extract posterior means
        beta_AR1 <- sum_AR1[paste0("beta[", 1:p, "]"), "mean"]
        rho_AR1  <- sum_AR1["rho", "mean"]
        
        # Store per-parameter metrics for beta
        res_list[[length(res_list) + 1]] <- data.frame(
          Rep      = r,
          Scenario = nome_cenario,
          Parameter = paste0("beta", 1:p),
          Prior    = paste0("Prior", k),
          ESS      = sum_AR1[paste0("beta[", 1:p, "]"), "n_eff"],
          Rhat     = sum_AR1[paste0("beta[", 1:p, "]"), "Rhat"],
          Bias     = beta_AR1 - beta_verd
        )
        
        # Store per-parameter metrics for rho
        res_list[[length(res_list) + 1]] <- data.frame(
          Rep      = r,
          Scenario = nome_cenario,
          Parameter = "rho",
          Prior    = paste0("Prior", k),
          ESS      = sum_AR1["rho", "n_eff"],
          Rhat     = sum_AR1["rho", "Rhat"],
          Bias     = rho_AR1 - rho_verd
        )
      }
    }
  }
}

# Combine results ----------------------------------------------------------

df_wide <- do.call(rbind, res_list)

# NOTE: This line relabels Prior4 as Prior1.
# This is only appropriate if, by design, Prior4 corresponds to the same
# hyperparameters as "Prior1" in your manuscript presentation.
df_wide$Prior[df_wide$Prior == "Prior4"] <- "Prior1"

# Reshape for ggplot: ESS/Rhat/Bias to long format -------------------------

df_long <- pivot_longer(
  df_wide,
  cols = c("ESS", "Rhat", "Bias"),
  names_to = "Metric",
  values_to = "Value"
)

# Aggregate beta parameters into a group label for plotting ----------------

df_long2 <- df_long %>%
  mutate(ParameterGroup = ifelse(grepl("beta", Parameter), "beta", "rho"))

# Plots: beta --------------------------------------------------------------

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

# Plots: rho ---------------------------------------------------------------

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

# Optional: quick check for extreme bias values ---------------------------
# Example shown for one scenario and one prior label.
tab_aux <- df_long %>%
  filter(Scenario == "T4_n50_Strong", Metric == "Bias", Prior == "Prior4")

tab_aux[tab_aux$Value > 2.5, ]

# Summary table: Estimate / ABias / RMSE ----------------------------------
# Map true rho values by scenario label (AR1 case)
rho_truth <- c(
  "T4_n50_Strong"     = 0.9,
  "T8_n50_Null"       = 0,
  "T12_n50_Moderate"  = 0.7,
  "T4_n100_Moderate"  = 0.7,
  "T8_n100_Strong"    = 0.9,
  "T12_n100_Null"     = 0,
  "T4_n200_Null"      = 0,
  "T8_n200_Moderate"  = 0.7,
  "T12_n200_Strong"   = 0.9
)

# Attach true values and compute posterior mean estimates from Bias
df_corrigido <- df_wide %>%
  mutate(
    True = case_when(
      Parameter == "beta1" ~ 1,
      Parameter == "beta2" ~ 0.5,
      Parameter == "beta3" ~ 0.8,
      Parameter == "beta4" ~ 0.6,
      Parameter == "rho"   ~ rho_truth[Scenario],
      TRUE ~ NA_real_
    ),
    Estimativa = Bias + True
  )

# Aggregate accuracy metrics by parameter and prior
df_summary <- df_corrigido %>%
  group_by(Prior, Parameter, True) %>%
  summarise(
    Estimate = round(mean(Estimativa), 3),
    ABias    = round(mean(abs(Estimativa - True), na.rm = TRUE), 3),
    RMSE     = round(sqrt(mean((Estimativa - True)^2, na.rm = TRUE)), 3),
    .groups  = "drop"
  ) %>%
  arrange(Parameter, Prior)

# Save figures -------------------------------------------------------------
# Create folder if needed
dir.create("Figures", showWarnings = FALSE)

## === beta figures === ##

g1 <- ggplot(filter(df_long2, Metric == "ESS", ParameterGroup == "beta"),
             aes(x = Scenario, y = Value, fill = Prior)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_grey(start = 0.3, end = 0.8) +
  labs(title = "Effective Sample Size (ESS)",
       x = "Scenario", y = "ESS", fill = "Prior") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Figures/ESS_beta.eps", g1, device = "eps", width = 8, height = 5)

g2 <- ggplot(filter(df_long2, Metric == "Rhat", ParameterGroup == "beta"),
             aes(x = Scenario, y = Value, fill = Prior)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_hline(yintercept = 1.1, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_grey(start = 0.3, end = 0.8) +
  labs(title = "Potential Scale Reduction Factor (Rhat)",
       x = "Scenario", y = "Rhat") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Figures/Rhat_beta.eps", g2, device = "eps", width = 8, height = 5)

g3 <- ggplot(filter(df_long2, Metric == "Bias", ParameterGroup == "beta"),
             aes(x = Scenario, y = Value, fill = Prior)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_grey(start = 0.3, end = 0.8) +
  labs(title = "Bias",
       x = "Scenario", y = "Bias") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Figures/Bias_beta.eps", g3, device = "eps", width = 8, height = 5)

## === rho figures === ##

g4 <- ggplot(filter(df_long2, Metric == "ESS", ParameterGroup == "rho"),
             aes(x = Scenario, y = Value, fill = Prior)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_grey(start = 0.3, end = 0.8) +
  labs(title = "Effective Sample Size (ESS)",
       x = "Scenario", y = "ESS") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Figures/ESS_rho.eps", g4, device = "eps", width = 8, height = 5)

g5 <- ggplot(filter(df_long2, Metric == "Rhat", ParameterGroup == "rho"),
             aes(x = Scenario, y = Value, fill = Prior)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_hline(yintercept = 1.1, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_grey(start = 0.3, end = 0.8) +
  labs(title = "Potential Scale Reduction Factor (Rhat)",
       x = "Scenario", y = "Rhat") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Figures/Rhat_rho.eps", g5, device = "eps", width = 8, height = 5)

g6 <- ggplot(filter(df_long2, Metric == "Bias", ParameterGroup == "rho"),
             aes(x = Scenario, y = Value, fill = Prior)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_grey(start = 0.3, end = 0.8) +
  labs(title = "Bias",
       x = "Scenario", y = "Bias") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Figures/Bias_rho.eps", g6, device = "eps", width = 8, height = 5)
