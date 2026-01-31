library(rstan)
library(rstansim)

options(mc.cores = parallel::detectCores())


PC<-'Usuário'

setwd(paste0('C:/Users/',PC,'/OneDrive/Documentos/Artigos/Multivariate Probit regression using antedependence/Programas/'))


source('Aux_Functions.R')

## Correlation ranges

#[.10, .40] Weak; [.40, .70] moderate; [.70, .9] strong

Scenario<-rbind(c(4,50,.9,.7),c(8,50,.4,.1),c(12,50,.7,.4),c(4,100,.7,.4),c(8,100,.9,.7),c(12,100,.4,.1),c(4,200,.4,.1),c(8,200,.7,.4),c(12,200,.9,.7))

R<-50


ESS_logPost_Lcorr<-matrix(NA,R,1)

ESS_logPost_AR1<-matrix(NA,R,1)


PSRF_logPost_Lcorr<-matrix(NA,R,1)

PSRF_logPost_AR1<-matrix(NA,R,1)

LStein_Lcorr<-matrix(NA,R,1)

LStein_AR1<-matrix(NA,R,1)





for(j in 1:nrow(Scenario)){
  
  
  ## Bulding design matrix
  
  vT<-Scenario[j,1]
  
  t_age<-1:vT
  
  X_s<-cbind(1,t_age,1,t_age)
  
  X_ns<-cbind(1,t_age,0,0)
  
  # Number of subjects
  
  n<-Scenario[j,2]
  
  # Number of subjects per group
  
  treat_counts<-c(n/2,n/2)
  
  
  p<-4
  
  # Criar um array tridimensional com a estrutura correta
  # Dimensões: (tempo x covariáveis x indivíduos)
  X_array <- array(NA, dim = c(length(t_age), p, sum(treat_counts)))
  
  # Preencher o array com as matrizes de planejamento para cada grupo
  X_array[, , 1:treat_counts[1]] <- replicate(treat_counts[1], X_ns, simplify = "array")
  X_array[, , (treat_counts[1]+1):sum(treat_counts)] <- replicate(treat_counts[2], X_s, simplify = "array")
  
  
  X_array<-aperm(X_array, c(3, 1, 2))
  
  
  ## Actual parameters
  
  # set.seed(258)
  
  beta_verd<-c(1,.5,.8,.6)
  
  rho_verd<-round(seq(Scenario[j,3],to=Scenario[j,4],length.out= vT-1),2)
  
  mu_verd_ns<-X_ns%*%beta_verd
  mu_verd_s<-X_s%*%beta_verd
  
  
  for (r in 1:R) {
    
    setwd("C:/Users/Usuário/OneDrive/Documentos/Artigos/Multivariate Probit regression using antedependence/ResultadosSim/")
    
    
    
    ESS_logPost_Lcorr[r,]<-readRDS(paste0("summaryMult_Lcorr_R",r,'_S',j,'.rds'))["lp__", "n_eff"]
    
    ESS_logPost_AR1[r,]<-readRDS(paste0("summaryMult_AR1_R",r,'_S',j,'.rds'))["lp__", "n_eff"]
    
    PSRF_logPost_Lcorr[r,]<-readRDS(paste0("summaryMult_Lcorr_R",r,'_S',j,'.rds'))["lp__", "Rhat"]
    
    PSRF_logPost_AR1[r,]<-readRDS(paste0("summaryMult_AR1_R",r,'_S',j,'.rds'))["lp__", "Rhat"]
    
    C_verd<-toeplitz.matrix(rep(1,vT),rho_verd)
    
    C_Lcorr<-reconstruir_correlacao(readRDS(paste0("summaryMult_Lcorr_R",r,'_S',j,'.rds')))
    
    rho_est<-readRDS(paste0("summaryMult_AR1_R",r,'_S',j,'.rds'))["rho", "mean"]
    
    C_AR1<-ARH1.matrix(rep(1,vT),rho_est)
    
    LStein_Lcorr[r,]<-sum(diag(C_Lcorr%*%solve(C_verd)))-log(det(C_Lcorr%*%solve(C_verd)))-vT
    LStein_AR1[r,]<-sum(diag(C_AR1%*%solve(C_verd)))-log(det(C_AR1%*%solve(C_verd)))-vT
    
  }
  
} 








