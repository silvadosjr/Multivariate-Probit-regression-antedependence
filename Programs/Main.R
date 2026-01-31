library(rstan)
library(rstansim)
library(shinystan)
library(devtools)
library(tidyr)
library(dplyr)


#devtools::install_github("ewan-keith/rstansim")


options(mc.cores = parallel::detectCores())


PC<-'Usuário'

setwd(paste0('C:/Users/',PC,'/OneDrive/Documentos/Artigos/Multivariate Probit regression using antedependence/Programas/'))

pathSim<-paste0('C:/Users/',PC,'/OneDrive/Documentos/Artigos/Multivariate Probit regression using antedependence/simData/')


source('Aux_Functions.R')

## ================= Six cities data ==========================================##

six_cities <- read.csv(paste0('C:/Users/',PC,'/OneDrive/Documentos/Artigos/Multivariate Probit regression using antedependence/Dados/six_cities_expanded.csv'))


## Bulding design matrix

vT<-4

#t_age<-7:10-9

t_age<-1:vT-mean(1:vT)

X_s<-cbind(1,t_age,1,t_age)

X_ns<-cbind(1,t_age,0,0)


# Number of subjects per group
#treat_counts <- table(six_cities$maternal_smoking)

treat_counts<-c(25,25)

n<-sum(treat_counts)

p<-4

# Criar um array tridimensional com a estrutura correta
# Dimensões: (tempo x covariáveis x indivíduos)
X_array <- array(NA, dim = c(length(t_age), p, sum(treat_counts)))

# Preencher o array com as matrizes de planejamento para cada grupo
X_array[, , 1:treat_counts[1]] <- replicate(treat_counts[1], X_ns, simplify = "array")
X_array[, , (treat_counts[1]+1):sum(treat_counts)] <- replicate(treat_counts[2], X_s, simplify = "array")


X_array<-aperm(X_array, c(3, 1, 2))


## Actual parameters

#set.seed(258)

beta_verd<-c(1,.5,.8,.6)

rho_verd<-round(seq(.9, to = .7, length.out = vT - 1), 2)

#R_true<-ARH1.matrix(rep(1,vT),rho_verd)

#R_true<-toeplitz.matrix(rep(1,vT),rho_verd)

#R_true<-AD_matrix(vT,rho_verd)

R_true<-make_corr_unstructured(vT,target_rho = rho_verd[1])

#R_true<-diag(vT)

mu_verd_ns<-X_ns%*%beta_verd
mu_verd_s<-X_s%*%beta_verd


## ===================== Simulated data ======================================##

data_simAux=list(vT=vT,p=length(beta_verd),n=n,X=X_array) 

param_sim1=list(beta=beta_verd,R=R_true)


begin<-Sys.time()
print(begin)

sim_AR1<-simulate_data(file = "simMultProbit.stan",data_name = "simMultProbit_Uns",
                              input_data = data_simAux,param_values = param_sim1,
                              vars =c('Y','R'),nsim = 1,
                              path = pathSim)

end<-Sys.time()
print(end-begin)


sim_AR1<-readRDS(paste0(pathSim,'/simMultProbit_AR1_1.rds'))

sim_Toep<-readRDS(paste0(pathSim,'/simMultProbit_Toep_1.rds'))

sim_AD<-readRDS(paste0(pathSim,'/simMultProbit_AD_1.rds'))

sim_Ind<-readRDS(paste0(pathSim,'/simMultProbit_Ind_1.rds'))

sim_Uns<-readRDS(paste0(pathSim,'/simMultProbit_Uns_1.rds'))

Y<-sim_Toep$Y


##================================= Fitting models ========================================##

#Y<-six_cities[,-5]

# Escolha da estrutura de correlação
# 1 = Toeplitz
# 2 = AR(1)
# 3 = AD(1)
# 4 = ARMA(1,1)

data_list=list(vT=vT,p=p,N=n,Y=Y,X=X_array,sigma_beta=10,sigma_rho=1,cor_type=3) 


params <- c('beta','rho_vec')
nChains = 1
burnInSteps = 1000
thinSteps = 10
numSavedSteps = 1000  #across all chains
nIter = ceiling(burnInSteps + (numSavedSteps * thinSteps)/nChains)
nIter

#nIter=11000


ini<-function() {list(beta=rep(.1,p),rho_vec=rep(.1,vT-1))}

begin<-Sys.time()
print(begin)


samp1 <- stan(data = data_list, file = "FitMultProbit_Structured.stan",init = ini,
              chains = nChains, pars = params, iter = nIter,
              warmup = burnInSteps, thin = thinSteps,control = list(adapt_delta = 0.9,max_treedepth=15),
              save_dso = T)


end<-Sys.time()
print(end-begin)



#saveRDS(samp1,paste0(pathSim,'ResultsFit_Lcorr_Real.rds'))



samp1<-readRDS(paste0(pathSim,'ResultsFit_AR1_Real.rds'))


launch_shinystan(samp1)


fit1 <- rstan::extract(samp1,permuted=T,inc_warmup=F)



# # Calcular autocorrelação de lag 5 para beta[1], considerando todas as cadeias
# autocorr_lag(samp1, "rho", lag = 15)
# 
# # Ou para a cadeia 1 apenas
# autocorr_lag(samp1, "beta[1]", lag = 5, chain = 3)


beta_amos<-fit1$beta

rho_amos<-fit1$rho_vec

#C_amos<-fit1$C

#(C_est=apply(C_amos, c(2,3), mean))


beta_est<-apply(beta_amos,2,mean)

rho_est<-apply(rho_amos,2,mean)

#rho_est<-mean(rho_amos)

C_est<-AD_matrix(vT,rho_est)

#C_est<-toeplitz.matrix(rep(1,vT),rho_est)

#C_verd<-ARH1.matrix(rep(1,vT),rho_verd)

R_true<-sim_Toep$R

LStein<-sum(diag(C_est%*%solve(R_true)))-log(det(C_est%*%solve(R_true)))-vT


require(boa)

intcred_beta<-apply(beta_amos,2,function(x){quantile(x,c(0.025,0.975))})
hpd_beta<-apply(beta_amos,2,function(x){boa.hpd(x,0.05)})


intcred_rho<-quantile(rho_amos,c(0.025,0.975))
hpd_rho<-boa.hpd(rho_amos,0.05)


IC<-cbind(intcred_beta,intcred_rho)



## Summary of posterior estimates


true_value<-c(beta_verd,rho_verd)

media<-c(apply(beta_amos,2,mean),mean(rho_amos))

mediana<-c(apply(beta_amos,2,median),median(rho_amos))

moda<-c(apply(beta_amos,2,MAP),MAP(rho_amos))

resu<-data.frame(true_value,media,mediana,moda,t(IC))

colnames(resu)<-c('True','EAP','MeAP','MoAP','2.5%','97.5%')
row.names(resu)<-c(paste0('beta',0:(p-1)),'rho')








