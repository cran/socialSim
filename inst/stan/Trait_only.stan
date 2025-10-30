data {
  // Phenotype dataset
  int<lower=0> n_obs;        // number of observations
  int<lower=0> n_ind;        // number of individuals

  array[n_obs] int<lower=0> individual;  // Individual ID
  array[n_obs] int<lower=0> opponent;  // Individual ID
  vector[n_obs] xj;                       // Opponent trait (covariate)
  vector[n_obs] z;                        // Phenotypic observations
}

transformed data {
  real mean_z = mean(z);
  vector[n_obs] z_c;
  for (i in 1:n_obs) {
    z_c[i] = z[i] - mean_z;
  }
}

parameters {
  // Fixed effects
  real B_0s;          // intercept
  real psi;           // population-level slope

  // Random effects
  vector[n_ind] zI;             // latent intercepts
  real<lower=0> sigma_I;        // SD for intercepts
  real<lower=0> sigma_e;        // residual SD
}

transformed parameters {
  vector[n_ind] I_alpha;
  I_alpha = zI * sigma_I;       // scaled intercepts
}

model {
  // Priors
  B_0s ~ normal(0, 1);
  psi  ~ normal(0, 1);
  zI ~ normal(0, 1);
  sigma_I ~ normal(0, 1);
  sigma_e ~ normal(0, 1);

  // Linear predictor
  vector[n_obs] e_z;
  e_z = B_0s + I_alpha[individual] + psi * xj;

  // Likelihood
  z_c ~ normal(e_z, sigma_e);
}

generated quantities {
  real B_0 = B_0s + mean_z;
  real<lower=0> Sigma2_intercept = square(sigma_I);
  real<lower=0> Sigma2_phi;
  real<lower=0> Sigma2_e = square(sigma_e);
  real<lower=0> Sigma2_x;

  vector[n_ind] xmean;
  vector[n_ind] phi;
  real cov_int_phi;

  // Compute xmean and phi per individual
  for (i in 1:n_ind) {
    real xsum = 0;
    int n = 0;
    for (j in 1:n_obs) {
      if (opponent[j] == i) {
        xsum += xj[j];
        n += 1;
          }
        }
    xmean[i] = xsum / n;
    phi[i] = psi * xmean[i]; }
  // Among-individual variance of x
  Sigma2_x = variance(xmean);
  // Covariance between phi and intercepts
  {
  real mean_phi = mean(phi);
  real mean_alpha = mean(I_alpha[individual]);
  real cov = 0;
  for (i in 1:n_ind) {
    cov += (phi[i] - mean_phi) * (I_alpha[i] - mean_alpha);
    }
    cov_int_phi = cov / (n_ind - 1);
    }
    Sigma2_phi = variance(phi);

}

