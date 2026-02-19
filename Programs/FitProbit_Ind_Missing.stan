data {
  // Dimensions
  int<lower=1> vT;                 // number of time points
  int<lower=1> p;                  // number of covariates
  int<lower=1> N;                  // number of subjects

  // Design matrices: X[n] is vT x p
  array[N] matrix[vT, p] X;

  // Responses: 0/1 observed, -1 missing
  array[N, vT] int Y;

  // Prior scale for regression coefficients
  real<lower=0> sigma_beta;

  // LKJ shape (fixed here; you can move to data if you want)
  // real<lower=0> lkj_shape;
}

transformed data {
  // Count observed positives/negatives and missing cells
  int N_pos = 0;
  int N_neg = 0;
  int N_mis = 0;

  for (n in 1:N) {
    for (t in 1:vT) {
      if (Y[n, t] == 1) N_pos += 1;
      else if (Y[n, t] == 0) N_neg += 1;
      else if (Y[n, t] == -1) N_mis += 1;
      else reject("Y must be 0, 1, or -1.");
    }
  }

  // Index arrays
  array[N_pos] int<lower=1, upper=N>  n_pos;
  array[N_pos] int<lower=1, upper=vT> t_pos;

  array[N_neg] int<lower=1, upper=N>  n_neg;
  array[N_neg] int<lower=1, upper=vT> t_neg;

  array[N_mis] int<lower=1, upper=N>  n_mis;
  array[N_mis] int<lower=1, upper=vT> t_mis;

  // Fill index arrays
  {
    int ip = 1;
    int ineg = 1;
    int im = 1;

    for (n in 1:N) {
      for (t in 1:vT) {
        if (Y[n, t] == 1) {
          n_pos[ip] = n; t_pos[ip] = t; ip += 1;
        } else if (Y[n, t] == 0) {
          n_neg[ineg] = n; t_neg[ineg] = t; ineg += 1;
        } else {
          n_mis[im] = n; t_mis[im] = t; im += 1;
        }
      }
    }
  }
}

parameters {
  vector[p] beta;

  // Unstructured correlation for latent Z across time
  cholesky_factor_corr[vT] Lcorr;

  // Latent variables, split by sign constraints implied by Y
  vector<lower=0>[N_pos] Z_pos;     // Y==1  => Z >= 0
  vector<upper=0>[N_neg] Z_neg;     // Y==0  => Z <= 0
  vector[N_mis]          Z_mis;     // missing => unconstrained
}

transformed parameters {
  matrix[N, vT] Z;

  // Initialize (avoid uninitialized entries)
  for (n in 1:N)
    for (t in 1:vT)
      Z[n, t] = 0;

  // Place constrained/unconstrained latents back into Z
  for (i in 1:N_pos) Z[n_pos[i], t_pos[i]] = Z_pos[i];
  for (i in 1:N_neg) Z[n_neg[i], t_neg[i]] = Z_neg[i];
  for (i in 1:N_mis) Z[n_mis[i], t_mis[i]] = Z_mis[i];
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

