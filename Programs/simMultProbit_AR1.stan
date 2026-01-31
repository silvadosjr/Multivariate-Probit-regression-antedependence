functions{
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
  int<lower=1> n;               // Número de indivíduos
  matrix[vT, p] X[n];           // Array com matrizes de delineamento
}

parameters {
  vector[p] beta;                      // Coeficientes de regressão
  real<lower=-.98, upper=.98> rho;         // Parâmetro AR(1)
}

transformed parameters {
  matrix[vT, vT] Sigma;               // Matriz de covariância ARH(1)
  matrix[vT, vT] L;                   // Cholesky de Sigma
  matrix[vT, vT] D_var;               // Matriz diagonal com variâncias
  matrix[vT, vT] L_std;               // L com diagonal padronizada
  matrix[vT, vT] T_mat;               // Inversa de L_std
  matrix[vT, vT] Phi_mat;             // Coef. atoregressivos generalizados 

  Sigma = ARH1_matrix(rep_vector(1.0, vT), rho, vT);
  L = cholesky_decompose(Sigma);
  vector[vT] d = diagonal(L);
  L_std = L * inverse(diag_matrix(d));       // L com diagonal unitária
  D_var = diag_matrix(d .* d);               // Variâncias (D^2)
  T_mat = inverse(L_std);
  Phi_mat = diag_matrix(rep_vector(1.0, vT)) - T_mat;
}

generated quantities {
  matrix[n, vT] Z;                         // Variáveis latentes
  int<lower=0, upper=1> Y[n, vT];          // Respostas observadas

  for (j in 1:n) {
    vector[vT] eta;
    eta[1] = dot_product(X[j, 1, ], beta);
    Z[j, 1] = normal_rng(eta[1], sqrt(D_var[1, 1]));
    Y[j, 1] = Z[j, 1] > 0;

    for (t in 2:vT) {
      eta[t] = dot_product(X[j, t, ], beta);
      for (k in 1:(t - 1)) {
        eta[t] += Phi_mat[t, k] * (Z[j, k] - dot_product(X[j, k, ], beta));
      }
      Z[j, t] = normal_rng(eta[t], sqrt(D_var[t, t]));
      Y[j, t] = Z[j, t] > 0;
    }
  }
}

