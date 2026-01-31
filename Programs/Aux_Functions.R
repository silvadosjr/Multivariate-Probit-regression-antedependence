
## parâmetros globais

# Y: matriz de dados (n x T) em que "n" é o total de observações por tempo e "T" é o número de tempos (serve para dados balanceados).

# sigma2: vetor de variâncias marginais

# mu: vetor médias marginais

# Sigma: matriz de covariância



## Para rescalonar os dados simulados com vetor de médias "mu" e matriz de covariância "Sigma". Retorna a matriz de dados "Y", reescalonada.

#Y<-Y_C[,-1];mu<-mu_verd1;Sigma<-Sigma_verd

autocorr_lag <- function(stanfit, param, lag = 1, chain = NULL) {
  # Extrai amostras brutas (não permutadas)
  samples <- rstan::extract(stanfit, pars = param, permuted = FALSE)  # array: iter x chains x parameters
  
  # Seleciona dimensão correta: supõe parâmetro escalar
  if (is.null(chain)) {
    # Agrega todas as cadeias
    series <- as.vector(samples[, , 1])
  } else {
    if (!chain %in% 1:dim(samples)[2]) stop("Número de cadeia inválido.")
    series <- samples[, chain, 1]
  }
  
  # Verifica se lag é válido
  if (lag >= length(series)) stop("Lag maior ou igual ao tamanho da cadeia.")
  
  # Calcula autocorrelação
  n <- length(series)
  x <- series - mean(series)
  acf_k <- sum(x[1:(n - lag)] * x[(lag + 1):n]) / sum(x^2)
  
  return(acf_k)
}




resc_Y<-function(Y,mu,Sigma,centred=TRUE){
  
  n<-nrow(Y);m<-ncol(Y)
  
  mu_Y <- matrix(1,n,1)%x%(cbind(apply(Y,2,mean)))
  
  muY_mat<-matrix(mu_Y,n,m,byrow=T)
  
  if(centred==T){
    
    Y_c<-Y - muY_mat
    
  }else{
    
    Y_c<-Y;mu<-rep(0,m)
    
  }
  
  
  covY<-cov(Y,use="complete.obs")
  
  Y_z<-t(solve(chol(covY)))%*%t(Y_c)
  
  mu_mat<-matrix(matrix(1,n,1)%x%mu,n,m,byrow=T)
  
  Y_s<-t(t(chol(Sigma))%*%Y_z) + mu_mat
  
  return(Y_s)
  
}


## Matrizes de covari?ncia

# Unstructured

make_corr_unstructured <- function(m, target_rho=0.4, mix=0.85, jitter=0.07, eps=1e-6){
  A <- matrix(rnorm(m*m),m,m)
  S <- crossprod(A)
  R <- cov2cor(S)
  off <- R - diag(m)
  mu <- mean(off[upper.tri(off)])
  sc <- target_rho/mu
  off <- off*sc
  R2 <- diag(m)+mix*off
  eig <- eigen(R2,sym=TRUE)
  lam <- pmax(eig$values,eps)
  Rpd <- eig$vectors%*%diag(lam)%*%t(eig$vectors)
  D <- sqrt(diag(Rpd))
  Rpd <- Rpd / outer(D,D)
  diag(Rpd)<-1
  (Rpd+t(Rpd))/2
}





# HU

HU.matrix<-function(sigma2,phi){
  
  m<-length(sigma2)
  
  U<-matrix(phi,m,m)
  
  diag(U)<-rep(1,m)
  
  
  mat_aux1<-sqrt(matrix(sigma2[matrix(rep(1:m,m),m,m,byrow=T)],m,m))
  mat_aux2<-sqrt(matrix(sigma2[matrix(rep(1:m,m),m,m,byrow=F)],m,m))
  
  U<-U*mat_aux1*mat_aux2
  
  return(U)
  
  
}



toeplitz.matrix<-function(sigma2,rho){
  
  m<-length(sigma2)
  
  R<-toeplitz(c(1,rho))
  
  mat_aux1<-sqrt(matrix(sigma2[matrix(rep(1:m,m),m,m,byrow=T)],m,m))
  mat_aux2<-sqrt(matrix(sigma2[matrix(rep(1:m,m),m,m,byrow=F)],m,m))
  
  return(R*mat_aux1*mat_aux2)
  
}




ARMAH11.matrix<-function(sigma2,phi,gama){
  
  
  vT<-length(sigma2)
  
  R<-gama*cov2cor(gama.arma(toeplitz(c(0,0:(vT-2))),1,phi,0))
  
  
  mat_aux1<-sqrt(matrix(sigma2[matrix(rep(1:vT,vT),vT,vT,byrow=T)],vT,vT))
  
  Sigma<-R*mat_aux1*t(mat_aux1)
  
  diag(Sigma)<-sigma2
  
  
  return(Sigma)
  
}



ARH1.matrix<-function(sigma2,rho){
  
  m<-length(sigma2)
  
  if(rho==0){
    
    return(diag(sigma2))
  
  }else{  
  Sigma<-gama.arma(toeplitz(0:(m-1)),1,rho,0)
  
  R<-(1/Sigma[1,1])*Sigma
  
  mat_aux1<-sqrt(matrix(sigma2[matrix(rep(1:m,m),m,m,byrow=T)],m,m))
  mat_aux2<-sqrt(matrix(sigma2[matrix(rep(1:m,m),m,m,byrow=F)],m,m))
  
  Sigma<-R*mat_aux1*mat_aux2
  
  diag(Sigma)<-sigma2
  
  return(Sigma)
  
  }
  
}

Prodacf<-function(rho,h){
  
  acf<-c()
  
  acf[1]<-1
  
  for(i in 2:h){
    
    acf[i]<-rho[i-1]*acf[i-1]
    
  }
  
  return(acf)
  
}



Prod.matrix<-function(vT,rho){
  
#  m<-length(sigma2)
  
  R<-toeplitz(Prodacf(rho,vT))
  
#  mat_aux1<-sqrt(matrix(sigma2[matrix(rep(1:m,m),m,m,byrow=T)],m,m))
#  mat_aux2<-sqrt(matrix(sigma2[matrix(rep(1:m,m),m,m,byrow=F)],m,m))
  
#  Sigma<-R*mat_aux1*mat_aux2
  
  
  return(R)
  
}


AD_matrix <- function(L, rho) {
  # L = dimensão
  # rho = vetor de comprimento L-1
  
  if (length(rho) != L - 1)
    stop("O vetor rho deve ter comprimento L - 1.")
  
  # Cria matriz identidade
  R <- diag(1, L)
  
  # Preenche triângulo superior e copia para inferior
  for (i in 1:(L-1)) {
    for (j in (i+1):L) {
      prod_ij <- prod(rho[i:(j-1)])  # produto dos rhos do intervalo
      R[i, j] <- prod_ij
      R[j, i] <- prod_ij
    }
  }
  
  return(R)
}




## Uniform matrix at [0,1] interval

uniform.matrix<-function(n,I){
  
  U<-matrix(runif(n*I),n,I)
  
  return(U)
  
}


MAP<-function(x){
  dd<-density(x)
  return(dd$x[which.max(dd$y)])
}


## autocovariance function of an ARMA process

gama.arma<-function(k,sigma2,phi,theta){
  
  gama.0<- sigma2*((1 + 2*theta*phi + theta^2) / (1-phi^2))
  
  gama.1<-sigma2*((1+theta*phi)*(phi+theta)/(1-phi^2))
  
  gama.k<-(phi^(k-1))*gama.1
  
  return(gama.0*(k==0) + gama.1*(k==1) + gama.k*(k>=2))
  
}

#' Reconstrói a matriz de correlação a partir de um summary de MCMC
#'
#' @param summary_df Data frame resultante da leitura de um arquivo .rds com parâmetros nomeados como "C[i,j]"
#' @param value_col Nome da coluna a ser usada (default = "mean")
#' @return Matriz de correlação estimada
reconstruir_correlacao <- function(summary_df, value_col = "mean") {
  # Identifica linhas com o padrão "C[i,j]"
  idx_cor <- grep("^C\\[", rownames(summary_df))
  cor_rows <- rownames(summary_df)[idx_cor]
  
  # Extrai o maior índice para definir a dimensão da matriz
  max_idx <- max(as.numeric(unlist(regmatches(cor_rows, gregexpr("[0-9]+", cor_rows)))))
  T <- max_idx
  
  # Inicializa matriz
  R_hat <- matrix(NA, nrow = T, ncol = T)
  
  # Preenche matriz com os valores extraídos
  for (row in cor_rows) {
    pos <- as.numeric(unlist(regmatches(row, gregexpr("[0-9]+", row))))
    i <- pos[1]
    j <- pos[2]
    R_hat[i, j] <- summary_df[row, value_col]
  }
  
  # Garante simetria
  R_hat[lower.tri(R_hat)] <- t(R_hat)[lower.tri(R_hat)]
  
  # Atribui nomes
  rownames(R_hat) <- colnames(R_hat) <- paste0("V", 1:T)
  
  return(R_hat)
}

