functions{
  
    int sum2d(array[,] int a){
    int s = 0;
    for (i in 1:size(a))
      s += sum(a[i]);
    return s;
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
 
 // Likelihood
 
  for (j in 1:N) {
    vector[vT] eta;
    for (t in 1:vT) {
    eta[t] = dot_product(X[j, t, ], beta);
    target+=normal_lpdf(Z[j,t]|eta[t],1);
    }
  }
}


