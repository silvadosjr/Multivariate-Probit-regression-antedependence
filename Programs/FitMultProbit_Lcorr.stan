functions{
  
    int sum2d(int[,] a){
    int s = 0;
    for (i in 1:size(a))
      s += sum(a[i]);
    return s;
  }
  
  real jeffreys_prior_nu(real nu, int K) {
    real acc = 0;
    for (k in 1:(K - 1)) {
      real term1 = trigamma(nu + 0.5 * (K - k - 1));
      real term2 = 2 * trigamma(2 * nu + K - k - 1);
      acc += term1 - term2;
    }
    return 0.5 * log(2 * acc);  // log of sqrt(...) = 0.5 * log(...)
  }
}
data {
  int<lower=1> vT;              // Número de tempos
  int<lower=1> p;               // Número de covariáveis
  int<lower=1> N;               // Número de indivíduos
  matrix[vT, p] X[N];           // Array com matrizes de delineamento
  int<lower=0,upper=1> Y[N,vT];
  real<lower=0> sigma_beta; 
}
transformed data {
  
  // For latent variable representation
  int<lower=0> N_pos;
  int<lower=1,upper=N> n_pos[sum2d(Y)];
  int<lower=1,upper=vT> d_pos[size(n_pos)];
  int<lower=0> N_neg;
  int<lower=1,upper=N> n_neg[(N*vT) - size(n_pos)];
  int<lower=1,upper=vT> d_neg[size(n_neg)];
  

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
 cholesky_factor_corr[vT] Lcorr;
// real<lower=0> nu; // for correlation matrix prior
  vector<lower=0>[N_pos] Z_pos;
  vector<upper=0>[N_neg] Z_neg;
}

transformed parameters {
  matrix[N,vT] Z;
  for (n in 1:N_pos)
    Z[n_pos[n], d_pos[n]] = Z_pos[n];
  for (n in 1:N_neg)
    Z[n_neg[n], d_neg[n]] = Z_neg[n];
}
model{
 
 // Prior
 
 beta~normal(0,sigma_beta);
 
// target += jeffreys_prior_nu(nu, vT);
 Lcorr~lkj_corr_cholesky(4);
 
 // Likelihood
 
  for (n in 1:N) {
    vector[vT] eta;
    eta = X[n, , ]*beta;
    target+=multi_normal_cholesky_lpdf(Z[n,]|eta,Lcorr);
  }
}
generated quantities{
 corr_matrix[vT] C;
 C=multiply_lower_tri_self_transpose(Lcorr);
}


