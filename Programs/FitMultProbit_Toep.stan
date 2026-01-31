functions{
  
    int sum2d(array[,] int a){
    int s = 0;
    for (i in 1:size(a))
      s += sum(a[i]);
    return s;
  }
  
   // Gera matriz de correlação com estrutura Toeplitz
  matrix Toep_matrix(int L, vector rhovec){
    matrix[L,L] R;
    int nr = rows(rhovec);
    for (i in 1:L) {
      for (j in 1:L) {
        if (i == j) {
          R[i,j] = 1;
        } else {
          int d = abs(i - j);
          if (nr == 1) {
            // HU: uma única correlação para todo par
            R[i,j] = rhovec[1];
          } else {
            // HT: correlação depende do lag d (1..L-1)
            R[i,j] = rhovec[d];
          }
        }
      }
    }
    return R;
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
  vector<lower=-.95, upper=.95>[vT-1] rho;         // correlation vector
  vector<lower=0>[N_pos] Z_pos;
  vector<upper=0>[N_neg] Z_neg;
}

transformed parameters {
  matrix[N,vT] Z;
  matrix[vT, vT] R;               // Matriz de correlação Toeplitz
  matrix[vT, vT] L;                   // Cholesky de Sigma
  matrix[vT, vT] D_var;               // Matriz diagonal com variâncias
  matrix[vT, vT] L_std;               // L com diagonal padronizada
  matrix[vT, vT] T_mat;               // Inversa de L_std
  matrix[vT, vT] Phi_mat;             // Coef. atoregressivos generalizados 
  
  
  R = Toep_matrix(vT, rho);
  L = cholesky_decompose(R);
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
 
 rho~normal(0,sigma_rho);
 
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



