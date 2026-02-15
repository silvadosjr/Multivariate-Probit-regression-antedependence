functions {

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

  //===========================================================
  // 1) Toeplitz correlation matrix (HU or HT)
  //
  // If rho has length 1:    compound symmetry / HU
  //   Corr(i,j) = rho[1] for i!=j
  //
  // If rho has length L-1:  Toeplitz-by-lag / HT
  //   Corr(i,j) = rho[|i-j|]
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
  // 2) AR(1) correlation matrix
  // Corr(i,j) = rho^|i-j|
  //===========================================================
  matrix AR1_matrix(int L, real rho) {
    matrix[L,L] R;
    for (i in 1:L)
      for (j in 1:L)
        R[i,j] = pow(rho, abs(i-j));
    return R;
  }

  //===========================================================
  // 3) Antedependence AD(1) correlation matrix
  //
  // Corr(i,j) = Π_{r=i}^{j-1} rho[r] for i<j
  // (and symmetric for i>j).
  //
  // Here rho is a vector of length L-1 with "local" dependence
  // parameters between consecutive occasions.
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
  // 4) ARMA(1,1) correlation matrix
  //
  // For a stationary ARMA(1,1), the lag-h autocorrelation is:
  //   Corr(0) = 1
  //   Corr(h) = phi^(h-1) * c, for h>=1
  // where
  //   c = ((phi + theta)(1 + phi*theta)) / (1 + theta^2 + 2*phi*theta)
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
  // Number of subjects
  int<lower=1> N;

  // Number of repeated measures / time points per subject
  int<lower=1> vT;

  // Number of regression coefficients
  int<lower=1> p;

  // Covariate array: for each subject j, X[j] is a vT x p matrix
  array[N] matrix[vT,p] X;

  // Binary outcomes: Y[j,t] in {0,1}
  array[N,vT] int<lower=0,upper=1> Y;

  // Correlation structure selector:
  // 1 = Toeplitz (HU/HT depending on rho length)
  // 2 = AR(1)
  // 3 = AD(1)
  // 4 = ARMA(1,1)
  int<lower=1,upper=4> cor_type;

  // Prior scales (hyperparameters)
  real<lower=0> sigma_beta;
  real<lower=0> sigma_rho;
}

transformed data {

  //===========================================================
  // Latent-variable bookkeeping for probit likelihood
  //
  // Under probit, we introduce Z[j,t] such that:
  //   Y[j,t] = 1  <=>  Z[j,t] > 0
  //   Y[j,t] = 0  <=>  Z[j,t] <= 0
  //
  // We therefore store indices of positive and negative cells
  // to impose truncation constraints via parameter bounds:
  //   Z_pos >= 0 and Z_neg <= 0
  //===========================================================

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

  // Regression coefficients
  vector[p] beta;

  // Generic vector for correlation parameters that require (vT-1) values:
  // - Toeplitz-by-lag (HT): rho_vec[lag]
  // - Antedependence AD(1): rho_vec[t] between t and t+1
  vector<lower=-0.98, upper=0.98>[vT-1] rho_vec;

  // Scalar AR(1) parameter
  real<lower=-0.98, upper=0.98> rho_ar1;
  

  // ARMA(1,1) parameters
  real<lower=-0.98, upper=0.98> phi;
  real<lower=-0.98, upper=0.98> theta;

  // Latent Z values, split into positive and negative parts
  // to enforce truncation consistent with Y.
  vector<lower=0>[N_pos]  Z_pos;
  vector<upper=0>[N_neg]  Z_neg;
}

transformed parameters {

  // Full latent matrix Z[j,t]
  matrix[N,vT] Z;

  // Selected correlation matrix and derived quantities
  matrix[vT, vT] R;       // correlation matrix implied by cor_type
  matrix[vT, vT] L;       // Cholesky factor of R (lower-triangular)
  matrix[vT, vT] D_var;   // diagonal matrix of conditional variances
  matrix[vT, vT] L_std;   // standardized Cholesky with unit diagonal
  matrix[vT, vT] T_mat;   // inverse of L_std
  matrix[vT, vT] Phi_mat; // generalized autoregressive coefficients

  //===========================================================
  // Build the correlation matrix R according to cor_type
  //===========================================================
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

  //===========================================================
  // Convert R into an antedependence-style conditional form
  //
  // Let L be the Cholesky of R: R = L * L'
  // Let d = diag(L).
  //
  // Define L_std = L * diag(d)^{-1}, so diag(L_std)=1.
  //
  // Then Z can be written through a triangular system with:
  //   T_mat = inv(L_std)
  // and
  //   Phi_mat = I - T_mat
  //
  // The diagonal D_var = diag(d^2) provides the conditional variances.
  //===========================================================
  L = cholesky_decompose(R);
 
  vector[vT] d = diagonal(L);
  L_std = L * inverse(diag_matrix(d));   // unit diagonal
  D_var = diag_matrix(d .* d);           // conditional variances
  T_mat = inverse(L_std);
  Phi_mat = diag_matrix(rep_vector(1.0, vT)) - T_mat;

  //===========================================================
  // Populate the full latent matrix Z from Z_pos and Z_neg
  //===========================================================
  for (n in 1:N_pos)
    Z[n_pos[n], d_pos[n]] = Z_pos[n];

  for (n in 1:N_neg)
    Z[n_neg[n], d_neg[n]] = Z_neg[n];
}

model{

  //===========================================================
  // Priors
  //===========================================================
  beta    ~ normal(0, sigma_beta);

  // Correlation parameters:
  // NOTE: even though only one structure is used at a time (cor_type),
  // all parameters receive priors. Since unused parameters do not appear
  // in the likelihood, their posterior will equal their prior.
  if(cor_type==1 || cor_type==3){
  rho_vec ~ normal(0, sigma_rho);
  }
  else if(cor_type==2){
  rho_ar1 ~ normal(0, sigma_rho);
  }
  else if(cor_type==4){
  phi     ~ normal(0, sigma_rho);
  theta   ~ normal(0, sigma_rho);
  }

  //===========================================================
  // Likelihood (latent-variable probit with conditional factorization)
  //
  // For each subject j, the latent vector Z[j,1:vT] follows a multivariate
  // normal with mean X*beta and correlation R. Instead of using a dense MVN,
  // we exploit the Cholesky-based conditional representation:
  //
  //   Z[j,1] ~ N(eta[1], D_var[1,1])
  //
  // For t>=2:
  //   Z[j,t] | Z[j,1:(t-1)] ~ N( eta[t], D_var[t,t] )
  // where eta[t] includes "innovations" based on Phi_mat and previous
  // residuals (Z - X*beta).
  //
  // The truncation implied by Y is enforced by the bounds on Z_pos/Z_neg.
  //===========================================================
  for (j in 1:N) {

    vector[vT] eta;

    // t = 1
    eta[1] = dot_product(X[j, 1, ], beta);
    target += normal_lpdf(Z[j, 1] | eta[1], sqrt(D_var[1, 1]));
    

    // t = 2..vT
    for (t in 2:vT) {

      // Start with the marginal mean at time t
      eta[t] = dot_product(X[j, t, ], beta);

      // Add dependence adjustment from previous occasions
      // based on generalized AR coefficients Phi_mat
      for (k in 1:(t - 1)) {
        eta[t] += Phi_mat[t, k] * (Z[j, k] - dot_product(X[j, k, ], beta));
      }

      // Conditional density contribution
      target += normal_lpdf(Z[j, t] | eta[t], sqrt(D_var[t, t]));
      
      
    }
  }
}

