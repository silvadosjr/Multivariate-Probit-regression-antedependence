functions{

  //===========================================================
  // Utility: sum all entries of a 2D integer array
  // Used to count the total number of Y==1 cells (positives).
  //===========================================================
  int sum2d(array[,] int a){
    int s = 0;
    for (i in 1:size(a))
      s += sum(a[i]);
    return s;
  }
}

data {
  // Number of repeated measures / time points per subject
  int<lower=1> vT;

  // Number of regression coefficients (covariates)
  int<lower=1> p;

  // Number of subjects
  int<lower=1> N;

  // Covariate array: for each subject j, X[j] is a vT x p matrix
  array[N] matrix[vT, p] X;

  // Binary outcomes Y[j,t] in {0,1}
  array[N, vT] int<lower=0, upper=1> Y;

  // Prior scales (hyperparameters)
  real<lower=0> sigma_beta;
}

transformed data {

  //===========================================================
  // Latent-variable bookkeeping for probit likelihood
  //
  // Probit augmentation:
  //   Y[j,t] = 1  <=>  Z[j,t] > 0
  //   Y[j,t] = 0  <=>  Z[j,t] <= 0
  //
  // We store indices of positive and negative cells to enforce
  // truncation via parameter bounds:
  //   Z_pos >= 0 and Z_neg <= 0
  //===========================================================

  int<lower=0> N_pos;
  array[sum2d(Y)] int<lower=1, upper=N> n_pos;
  array[size(n_pos)] int<lower=1, upper=vT> d_pos;

  int<lower=0> N_neg;
  array[(N * vT) - size(n_pos)] int<lower=1, upper=N> n_neg;
  array[size(n_neg)] int<lower=1, upper=vT> d_neg;

  // Total counts of positives/negatives
  N_pos = size(n_pos);
  N_neg = size(n_neg);

  // Fill index arrays for Y==1 (positives) and Y==0 (negatives)
  {
    int i;
    int j;
    i = 1;
    j = 1;
    for (n in 1:N) {
      for (d in 1:vT) {
        if (Y[n, d] == 1) {
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
  // Regression coefficients
  vector[p] beta;

  // Latent variables split by sign constraints implied by Y
  vector<lower=0>[N_pos] Z_pos;   // entries where Y==1 (force Z >= 0)
  vector<upper=0>[N_neg] Z_neg;   // entries where Y==0 (force Z <= 0)
}

transformed parameters {
  // Full latent matrix Z[j,t]
  matrix[N, vT] Z;

  // Place constrained latent values into their original (subject,time) positions
  for (n in 1:N_pos)
    Z[n_pos[n], d_pos[n]] = Z_pos[n];

  for (n in 1:N_neg)
    Z[n_neg[n], d_neg[n]] = Z_neg[n];
}

model{

  //===========================================================
  // Priors
  //===========================================================
  beta ~ normal(0, sigma_beta);

  //===========================================================
  // Likelihood: independent probit with latent augmentation
  //
  // For each subject j and time t:
  //   Z[j,t] ~ Normal( X[j,t,]*beta , 1 )
  // with truncation enforced by Z_pos/Z_neg bounds.
  //
  // This corresponds to:
  //   P(Y[j,t]=1 | X) = Phi( X[j,t,]*beta )
  // assuming independence across t given covariates.
  //===========================================================
  for (j in 1:N) {
    vector[vT] eta;

    for (t in 1:vT) {
      eta[t] = dot_product(X[j, t, ], beta);
      target += normal_lpdf(Z[j, t] | eta[t], 1);
    }
  }
}
generated quantities{
  
vector[N] log_lik;

  for (j in 1:N) {
    vector[vT] eta;
    real ll = 0;

    for (t in 1:vT) {
      eta[t] = dot_product(X[j, t, ], beta);
      ll += normal_lpdf(Z[j, t] | eta[t], 1);
    }
    
  log_lik[j]=ll;  
    
  }  
}

