data {
//Phenotype dataset

  // Number of clusters (an integer)
   int<lower=0> n_obs; // number of observations for phenotypes or aantal rows
   int<lower=0> n_ind; //number of individuals

  array[n_obs] int<lower=0> individual;  // Individual ID
  array[n_obs] int<lower=0> opponent;    // Opponent ID
  array[n_obs] real z;                   // Phenotypic observations
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

   // Random effects
   matrix[2,n_ind]         zI; //(intercepts and slopes, res_impact for each individual)
   vector<lower=0>[2]      sigma_I; // sd  intercepts, slopes, res_impact
   cholesky_factor_corr[2] L;  // factor to estimate covariance int-slopes
   real<lower=0> sigma_e;
 }

 transformed parameters{
    matrix[2,n_ind] I;     //  Unscaled blups intercept and slope and res_impact for each individual
    array[n_obs] real e_z; // predicted values for covariate
    I  = diag_pre_multiply(sigma_I, L) * zI; // get the unscaled value


   for (i in 1:n_obs) {
     e_z[i]  = B_0s  +  I[1, individual[i]] + I[2, opponent[i]];
   }
 }

model {
// Create vector of predicted values
  to_vector([B_0s]) ~ normal(0, 1);
   sigma_e ~ normal(0, 1);

 // Random effects distribution
    to_vector(zI) ~ normal(0, 1);
    to_vector(sigma_I) ~ normal(0, 1);
    L ~ lkj_corr_cholesky(2);

 // Likelihood function
    for (i in 1:n_obs)
     z_c[i]~normal(e_z[i], sigma_e);

 }

generated quantities{
real B_0 = B_0s + mean_z;
real<lower=0> Sigma2_intercept;
real<lower=0> Sigma2_phi;
real<lower=0> Sigma2_e;

real cor_2; //   int-resimpact
real cov_2; //   int-resimpact
real cov_int_phi;


matrix[2, 2]  Omega_I;

Sigma2_intercept=sigma_I[1]^2;
Sigma2_phi=sigma_I[2]^2;
Sigma2_e=sigma_e^2;

Omega_I = L * L';
cor_2 = Omega_I[1,2];
cov_2 = Omega_I[1,2]*sqrt(Sigma2_phi*Sigma2_intercept);
cov_int_phi = cov_2;
 }
