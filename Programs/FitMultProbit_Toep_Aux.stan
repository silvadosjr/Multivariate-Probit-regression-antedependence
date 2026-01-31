functions{
  // Soma todos os elementos de um array 2D
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
  int<lower=1> vT;                 // Número de tempos
  int<lower=1> p;                  // Número de covariáveis
  int<lower=1> N;                  // Número de indivíduos
  array[N] matrix[vT, p] X;        // Array com matrizes de delineamento
  array[N, vT] int<lower=0, upper=1> Y;
  real<lower=0> sigma_beta;        // DP do prior de beta
  real<lower=0> sigma_rho;         // DP do prior de rho
}

transformed data {
  // Índices da representação truncada dos latentes
  int<lower=0> N_pos;
  array[sum2d(Y)] int<lower=1, upper=N>   n_pos;
  array[size(n_pos)] int<lower=1, upper=vT> d_pos;
  int<lower=0> N_neg;
  array[(N * vT) - size(n_pos)] int<lower=1, upper=N>   n_neg;
  array[size(n_neg)] int<lower=1, upper=vT> d_neg;
  
  N_pos = size(n_pos);
  N_neg = size(n_neg);
  {
    int i = 1;
    int j = 1;
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
  vector[p] beta;                               // Coeficientes de regressão
  vector<lower=-0.98, upper=0.98>[vT-1] rho;   // Correlações Toeplitz (lags 1..vT-1)
  
  // Latentes truncados
  vector<lower=0>[N_pos]  Z_pos;
  vector<upper=0>[N_neg]  Z_neg;
}

transformed parameters {
  matrix[N, vT] Z;               // Latentes completos
  matrix[vT, vT] R;              // Matriz de correlação Toeplitz
  matrix[vT, vT] L;              // Cholesky de R
  vector[vT] d;                  // Diagonal de L
  matrix[vT, vT] L_std;          // L com diagonal unitária
  matrix[vT, vT] T_mat;          // Inversa de L_std
  matrix[vT, vT] Phi_mat;        // Coeficientes de antedependência
  vector[vT] sd_eps;             // DP das inovações (variâncias condicionais^0.5)
  
  // Reconstrói Z a partir dos vetores truncados
  for (n in 1:N_pos)
    Z[n_pos[n], d_pos[n]] = Z_pos[n];
  for (n in 1:N_neg)
    Z[n_neg[n], d_neg[n]] = Z_neg[n];
  
  // Matriz de correlação e decomposição de Cholesky
  R = Toep_matrix(vT, rho);
  L = cholesky_decompose(R);
  d = diagonal(L);
  
  // L com diagonal = 1 (normalização)
  L_std = L * diag_matrix(1 ./ d);
  
  // T_mat = L_std^{-1} usando solver triangular (evita inverse())
  T_mat = mdivide_left_tri_low(diag_matrix(rep_vector(1.0, vT)), L_std);
  
  // Phi = I - T (coeficientes autoregressivos generalizados)
  Phi_mat = diag_matrix(rep_vector(1.0, vT)) - T_mat;
  
  // DP das inovações (condicionais); variância condicional = d[t]^2
  sd_eps = d;
}

model {
  // Priors
  beta ~ normal(0, sigma_beta);
  rho  ~ normal(0, sigma_rho);
  
  // Likelihood via representação de antedependência
  for (n in 1:N) {
    vector[vT] eta;   // médias condicionais
    
    // t = 1
    eta[1] = dot_product(X[n, 1, ], beta);
    Z[n, 1] ~ normal(eta[1], sd_eps[1]);
    
    // t >= 2
    for (t in 2:vT) {
      real mu_t;
      mu_t = dot_product(X[n, t, ], beta);
      // soma dos termos de correlação (antedependência)
      for (k in 1:(t - 1)) {
        mu_t += Phi_mat[t, k] * (Z[n, k] - dot_product(X[n, k, ], beta));
      }
      eta[t] = mu_t;
      Z[n, t] ~ normal(eta[t], sd_eps[t]);
    }
  }
}
