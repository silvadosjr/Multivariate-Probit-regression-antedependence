functions{
  
    int sum2d(array[,] int a){
    int s = 0;
    for (i in 1:size(a))
      s += sum(a[i]);
    return s;
  }
  
  // Gera matriz de covariância com estrutura ARH(1)
  matrix ARH1_matrix(vector sigma2, real rho, int nrows) {
    matrix[nrows, nrows] R;
    matrix[nrows, nrows] sigma_mat;

    R = diag_matrix(rep_vector(1.0, nrows));
    for (t in 2:nrows) {
      for (k in 1:(t - 1)) {
        real val = pow(rho, t - k);
        R[t, k] = val;
        R[k, t] = val;
      }
    }

    // Constrói matriz de variâncias
    sigma_mat = rep_matrix(sqrt(sigma2)', nrows);
    return R .* sigma_mat .* sigma_mat';  // Matriz de covariância
  }
}
data {
  int<lower=1> vT;              // Número de tempos
  int<lower=1> p;               // Número de covariáveis
  int<lower=1> N;               // Número de indivíduos
  array[N] matrix[vT, p] X;           // Array com matrizes de delineamento
  array[N, vT] int<lower=0, upper=1> Y;
  real<lower=0> sigma_beta; 
  real<lower=0> sigma_rho;
}
transformed data {
  
  // For latent variable representation
  int<lower=0> N_pos;
  array[sum2d(Y)] int<lower=1, upper=N> n_pos;
  array[size(n_pos)] int<lower=1, upper=vT> d_pos;
  int<lower=0> N_neg;
  array[(N * vT) - size(n_pos)] int<lower=1, upper=N> n_neg;
  array[size(n_neg)] int<lower=1, upper=vT> d_neg;
  

  N_pos = size(n_pos);
  N_neg = size(n_neg);
  {
    int i;
    int j;
    i = 1;
    j = 1;
    for (n in 1:N) {
      for (d in 1:vT) {
        if (Y[n,d] == 1) {
          n_pos[i] = n;
          d_pos[i] = d;
          i += 1;
        } else {
          n_neg[j] = n;
          d_neg[j] = d;
          j += 1;
        }
      }
    }
  }
}
parameters {
  vector[p] beta;                      // Coeficientes de regressão
  real<lower= -.98, upper=.98> rho_ar1;         // Parâmetro AR(1)
  vector<lower=0>[N_pos] Z_pos;
  vector<upper=0>[N_neg] Z_neg;
}

transformed parameters {
  matrix[N,vT] Z;
  matrix[vT, vT] Sigma;               // Matriz de covariância ARH(1)
  matrix[vT, vT] L;                   // Cholesky de Sigma
  matrix[vT, vT] D_var;               // Matriz diagonal com variâncias
  matrix[vT, vT] L_std;               // L com diagonal padronizada
  matrix[vT, vT] T_mat;               // Inversa de L_std
  matrix[vT, vT] Phi_mat;             // Coef. atoregressivos generalizados 

  Sigma = ARH1_matrix(rep_vector(1.0, vT), rho_ar1, vT);
  L = cholesky_decompose(Sigma);
  vector[vT] d = diagonal(L);
  L_std = L * inverse(diag_matrix(d));       // L com diagonal unitária
  D_var = diag_matrix(d .* d);               // Variâncias (D^2)
  T_mat = inverse(L_std);
  Phi_mat = diag_matrix(rep_vector(1.0, vT)) - T_mat;
  
  for (n in 1:N_pos)
    Z[n_pos[n], d_pos[n]] = Z_pos[n];
  for (n in 1:N_neg)
    Z[n_neg[n], d_neg[n]] = Z_neg[n];
}
model{
 
 // Prior
 
 beta~normal(0,sigma_beta);
 
 rho_ar1~normal(0,sigma_rho);
 
 // Likelihood
 
  for (j in 1:N) {
    vector[vT] eta;
    eta[1] = dot_product(X[j, 1, ], beta);
    target+=normal_lpdf(Z[j,1]|eta[1],sqrt(D_var[1,1]));

    for (t in 2:vT) {
      eta[t] = dot_product(X[j, t, ], beta);
      for (k in 1:(t - 1)) {
        eta[t] += Phi_mat[t, k] * (Z[j, k] - dot_product(X[j, k, ], beta));
      }
      target+=normal_lpdf(Z[j,t]|eta[t],sqrt(D_var[t,t]));
    }
  }
}


