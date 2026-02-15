functions{

  //===========================================================
  // Utility: sum all entries of a 2D int array
  // Used to count the total number of Y==1 cells (positives).
  //===========================================================
  int sum2d(int[,] a){
    int s = 0;
    for (i in 1:size(a))
      s += sum(a[i]);
    return s;
  }

  //===========================================================
  // Optional: Jeffreys prior for an LKJ "shape" parameter nu
  //
  // If you want to estimate a hyperparameter nu controlling the amount
  // of shrinkage toward the identity in an LKJ prior, a Jeffreys-type
  // prior can be used. This function returns log p(nu) up to a constant,
  // based on trigamma terms (polygamma of order 1).
  //
  // NOTE: In the current model, nu is commented out, so this function is
  // not used unless you uncomment nu in parameters and add it to the model.
  //===========================================================
  real jeffreys_prior_nu(real nu, int K) {
    real acc = 0;
    for (k in 1:(K - 1)) {
      real term1 = trigamma(nu + 0.5 * (K - k - 1));
      real term2 = 2 * trigamma(2 * nu + K - k - 1);
      acc += term1 - term2;
    }
    // Returns log(sqrt(2*acc)) = 0.5*log(2*acc)
    return 0.5 * log(2 * acc);
  }
}

data {
  // Number of repeated measures / time points
  int<lower=1> vT;

  // Number of regression coefficients (covariates)
  int<lower=1> p;

  // Number of subjects
  int<lower=1> N;

  // Design matrices: X[n] is vT x p for subject n
  matrix[vT, p] X[N];

  // Binary responses
  int<lower=0,upper=1> Y[N, vT];

  // Prior scale for beta
  real<lower=0> sigma_beta;
}

transformed data {

  //===========================================================
  // Latent-variable bookkeeping for probit likelihood
  //
  // Probit augmentation:
  //   Y[n,t] = 1  <=>  Z[n,t] > 0
  //   Y[n,t] = 0  <=>  Z[n,t] <= 0
  //
  // We store indices of positive and negative entries so we can represent
  // Z with two parameter vectors that are constrained by sign:
  //   Z_pos >= 0 and Z_neg <= 0
  //===========================================================

  int<lower=0> N_pos;
  int<lower=1,upper=N>  n_pos[sum2d(Y)];
  int<lower=1,upper=vT> d_pos[size(n_pos)];

  int<lower=0> N_neg;
  int<lower=1,upper=N>  n_neg[(N*vT) - size(n_pos)];
  int<lower=1,upper=vT> d_neg[size(n_neg)];

  N_pos = size(n_pos);
  N_neg = size(n_neg);

  // Fill (subject,time) index vectors for Y==1 and Y==0
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

  // Cholesky factor of the correlation matrix for Z[n,]
  // This implies an unstructured (dense) correlation across time.
  cholesky_factor_corr[vT] Lcorr;

  // Optional: LKJ shape parameter (commented out here)
  // real<lower=0> nu;

  // Latent Z values split by sign constraints
  vector<lower=0>[N_pos] Z_pos;   // entries where Y==1
  vector<upper=0>[N_neg] Z_neg;   // entries where Y==0
}

transformed parameters {

  // Full latent matrix Z[n,t]
  matrix[N, vT] Z;

  // Reconstruct Z from constrained vectors
  for (n in 1:N_pos)
    Z[n_pos[n], d_pos[n]] = Z_pos[n];

  for (n in 1:N_neg)
    Z[n_neg[n], d_neg[n]] = Z_neg[n];
}

model{

  //===========================================================
  // Priors
  //===========================================================

  // Regression prior
  beta ~ normal(0, sigma_beta);

  // LKJ prior for the correlation matrix via its Cholesky factor.
  // Larger shape values concentrate mass near the identity matrix.
  //
  // Here shape = 4 is a fairly informative prior favoring weak correlations.
  // If you want to learn the shape parameter nu, uncomment nu and use:
  //   target += jeffreys_prior_nu(nu, vT);
  //   Lcorr ~ lkj_corr_cholesky(nu);
  Lcorr ~ lkj_corr_cholesky(4);

  //===========================================================
  // Likelihood (latent probit with unstructured correlation)
  //
  // For each subject n:
  //   Z[n,] ~ MVN( eta[n,],  C )
  // where C is a correlation matrix (unit variances) and its Cholesky
  // factor is Lcorr.
  //
  // NOTE: multi_normal_cholesky_lpdf expects the Cholesky factor of the
  // covariance matrix. Since C is a correlation matrix, its Cholesky factor
  // is directly Lcorr (unit marginal variances).
  //
  // Truncation implied by Y is enforced by bounds on Z_pos/Z_neg.
  //===========================================================
  for (n in 1:N) {
    vector[vT] eta;

    // Linear predictor across time for subject n
    eta = X[n, , ] * beta;

    // Multivariate normal contribution with correlated errors
    target += multi_normal_cholesky_lpdf(Z[n,] | eta, Lcorr);
  }
}

generated quantities{
  // Recover the full correlation matrix C from its Cholesky factor
  corr_matrix[vT] C;
  C = multiply_lower_tri_self_transpose(Lcorr);
}

