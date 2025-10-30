data {
//Phenotype dataset

  // Number of clusters (an integer)
   int<lower=0> n_obs; // number of observations
   int<lower=0> n_ind; //number of individuals

  array[n_obs] int<lower=0> individual;  // Individual ID
  array[n_obs] int<lower=0> opponent;    // Opponent ID
  vector[n_obs] xj;                  // Observed covariate
  vector[n_obs] z;                   // Phenotypic observations

}
transformed data {
  real mean_z = mean(z);
  array[n_obs] real z_c;
  for (i in 1:n_obs) {
    z_c[i]  = z[i]  - mean_z;
}
}
 parameters {
// Define parameters to estimate
  // Parameters for model of phenotype
     real B_0s; //intercept
     real psi; //slope


   // Random effects
   matrix[3,n_ind]         zI; //(intercepts and slopes, res_impact for each individual)
   vector<lower=0>[3]      sigma_I; // sd  intercepts, slopes, res_impact
   cholesky_factor_corr[3] L;  // factor to estimate covariance int-slopes
   real<lower=0> sigma_e;
 }

 transformed parameters{
    matrix[n_ind,3] I; //  Unscaled blups intercept and slope and res_impact for each individual
    I = (diag_pre_multiply(sigma_I, L) * zI)'; // get the unscaled value

 }

model {
// Create vector of predicted values
  B_0s ~ normal(0, 1);
  psi ~ normal(0, 1);


 // Random effects distribution
    to_vector(zI) ~ normal(0, 1);
    sigma_I ~ normal(0, 1);
    L ~ lkj_corr_cholesky(2);

    sigma_e ~ normal(0,1);

 // Likelihood function
    vector[n_ind] I_alpha = col(I,1); //intercepts
    vector[n_ind] I_psi = col(I,2); // responsiveness
    vector[n_ind] I_epsilon = col(I,3); // residual impact

    vector[n_obs] e_z;

     e_z  = B_0s  +
            I_alpha[individual] +
            (psi + I_psi[individual]) .* xj +
            I_epsilon[opponent];
     z_c ~ normal(e_z, sigma_e);
 }

generated quantities {
  real B_0 = B_0s + mean_z;

  real<lower=0> Sigma2_intercept;
  real<lower=0> Sigma2_psi;
  real<lower=0> Sigma2_epsilon;
  real<lower=0> Sigma2_phi;
  real<lower=0> Sigma2_e;
  real<lower=0> Sigma2_x;

  real cov_1; // psi-int
  real cov_2; // int-resimpact
  real cov_3; // psi-resimpact

  matrix[3, 3] Omega_I;

  vector[n_ind] xmean;
  vector[n_ind] phi;
  real cov_int_phi;

  // Variances
  Sigma2_intercept = square(sigma_I[1]);
  Sigma2_psi       = square(sigma_I[2]);
  Sigma2_epsilon   = square(sigma_I[3]);
  Sigma2_e         = square(sigma_e);

  // Correlation matrix
  Omega_I = L * L';
  cov_1 = Omega_I[1,2] * sqrt(Sigma2_psi * Sigma2_intercept);
  cov_2 = Omega_I[1,3] * sqrt(Sigma2_epsilon * Sigma2_intercept);
  cov_3 = Omega_I[2,3] * sqrt(Sigma2_psi * Sigma2_epsilon);

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
    phi[i] = I[i, 3] + psi * xmean[i];
  }

  // Among-individual variance of x
  Sigma2_x = variance(xmean);

  // Covariance between phi and intercepts
  {
    real mean_phi = mean(phi);
    real mean_alpha = mean(col(I, 1));
    real cov = 0;
    for (i in 1:n_ind) {
      cov += (phi[i] - mean_phi) * (I[i, 1] - mean_alpha);
    }
    cov_int_phi = cov / (n_ind - 1);
  }

  Sigma2_phi = variance(phi);

 }
