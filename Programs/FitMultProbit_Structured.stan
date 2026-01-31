functions {

  //===========================================================
  // SUM OF 2D ARRAY
  //===========================================================
  int sum2d(array[,] int a){
    int s = 0;
    for (i in 1:size(a))
      s += sum(a[i]);
    return s;
  }

  //===========================================================
  // 1. TOEPLITZ: HU (1 parâmetro) ou HT (vT−1 parâmetros)
  // rho: length 1  => HU
  // rho: length vT−1 => HT
  //===========================================================
  matrix Toep_matrix(int L, vector rho) {
    matrix[L,L] R;
    int nr = rows(rho);

    for (i in 1:L)
      for (j in 1:L) {
        if (i == j)
          R[i,j] = 1;
        else {
          int d = abs(i-j);
          R[i,j] = (nr == 1 ? rho[1] : rho[d]);
        }
      }
    return R;
  }

  //===========================================================
  // 2. AR(1): Corr(i,j) = rho^|i-j|
  //===========================================================
  matrix AR1_matrix(int L, real rho) {
    matrix[L,L] R;
    for (i in 1:L)
      for (j in 1:L)
        R[i,j] = pow(rho, abs(i-j));
    return R;
  }

  //===========================================================
  // 3. ANTEDEPENDÊNCIA AD(1)
  // Corr(i,j) = Π_{r=min(i,j)}^{max(i,j)-1} rho[r]
  //===========================================================
  matrix AD_matrix(int L, vector rho) {
    matrix[L,L] R;

    for (i in 1:L)
      R[i,i] = 1;

    for (i in 1:(L-1))
      for (j in (i+1):L) {
        real prod = 1;
        for (r in i:(j-1))
          prod *= rho[r];
        R[i,j] = prod;
        R[j,i] = prod;
      }

    return R;
  }

  //===========================================================
  // 4. ARMA(1,1) correlation matrix
  // Corr(h) = 
  //   1, if h=0
  //   phi^(h-1) * ((phi + theta)(1 + phi*theta)) /
  //                 (1 + theta^2 + 2*phi*theta)
  //===========================================================
  matrix ARMA11_matrix(int L, real phi, real theta) {
    matrix[L,L] R;

    real denom = 1 + theta^2 + 2*phi*theta;
    real c = ((phi + theta) * (1 + phi*theta)) / denom;

    for (i in 1:L)
      R[i,i] = 1;

    for (i in 1:(L-1))
      for (j in (i+1):L) {
        int h = j - i;
        real cor = pow(phi, h-1) * c;
        R[i,j] = cor;
        R[j,i] = cor;
      }

    return R;
  }

}
data {
  int<lower=1> N;
  int<lower=1> vT;
  int<lower=1> p;

  array[N] matrix[vT,p] X;
  array[N,vT] int<lower=0,upper=1> Y;

  // Escolha da estrutura de correlação
  // 1 = Toeplitz
  // 2 = AR(1)
  // 3 = AD(1)
  // 4 = ARMA(1,1)
  int<lower=1,upper=4> cor_type;

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
  vector[p] beta;

  // Toeplitz ou AD: precisaremos de vT-1 rhos
  vector<lower=-0.98, upper=0.98>[vT-1] rho_vec;

  // AR(1): um único rho
  real<lower=-0.98, upper=0.98> rho_ar1;

  // ARMA(1,1)
  real<lower=-0.98, upper=0.98> phi;
  real<lower=-0.98, upper=0.98> theta;

  // Latentes Z
  vector<lower=0>[N_pos]  Z_pos;
  vector<upper=0>[N_neg]  Z_neg;
}

transformed parameters {
  matrix[N,vT] Z;
  matrix[vT, vT] R;               // Matriz de correlação Toeplitz
  matrix[vT, vT] L;                   // Cholesky de Sigma
  matrix[vT, vT] D_var;               // Matriz diagonal com variâncias
  matrix[vT, vT] L_std;               // L com diagonal padronizada
  matrix[vT, vT] T_mat;               // Inversa de L_std
  matrix[vT, vT] Phi_mat;             // Coef. atoregressivos generalizados 
  
  
  // Seleção da estrutura de correlação
  if (cor_type == 1) {
    R = Toep_matrix(vT, rho_vec);
  }
  else if (cor_type == 2) {
    R = AR1_matrix(vT, rho_ar1);
  }
  else if (cor_type == 3) {
    R = AD_matrix(vT, rho_vec);
  }
  else if (cor_type == 4) {
    R = ARMA11_matrix(vT, phi, theta);
  }
  
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
 
 //-------- Priors --------//
  beta      ~ normal(0, sigma_beta);

  rho_vec   ~ normal(0, sigma_rho);
  rho_ar1   ~ normal(0, sigma_rho);
  phi       ~ normal(0, sigma_rho);
  theta     ~ normal(0, sigma_rho);
 
 //-------- Likelihood --------//
 
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



