
data {
  //Phenotype dataset

  // Number of clusters (an integer)
  int<lower=0> n_obs; // number of observations
  int<lower=0> n_ind; //number of individuals

  array[n_obs] int<lower=0> individual;  // Individual ID
  array[n_obs] int<lower=0> opponent;    // Opponent ID
  array[n_obs] real xj;                  // Observed covariate
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
  real B_x; //intercept covariate model
  real B_0s; //intercept
  real psi; //slope


  // Random effects
  matrix[3,n_ind]         zI; //(intercepts and slopes, res_impact for each individual)
  vector<lower=0>[3]      sigma_I; // sd  intercepts, slopes, res_impact
  cholesky_factor_corr[3] L;  // factor to estimate covariance int-slopes
  real<lower=0> sigma_e;
  real<lower=0> sigma_ex;
}

transformed parameters{
  matrix[n_ind,3] I; //  Unscaled blups intercept and slope and res_impact for each individual
  I = (diag_pre_multiply(sigma_I, L) * zI)'; // get the unscaled value

 }

model {
// Create vector of predicted values
  B_0s ~ normal(0, 1);
  B_x ~ normal(0, 1);
  psi ~ normal(0, 1);

 // Random effects distribution
    to_vector(zI) ~ normal(0, 1);
    sigma_I ~ normal(0, 1);
    L ~ lkj_corr_cholesky(2);

    sigma_e ~ normal(0,1);
    sigma_ex ~ normal(0,1);

 // Likelihood function
    vector[n_ind] I_alpha = col(I,1); //intercepts
    vector[n_ind] I_epsilon = col(I,2); // residual impact
    vector[n_ind] I_x = col(I,3); // residual impact

    vector[n_obs] e_x;
    vector[n_obs] e_z;

     e_x = B_x + I_x[opponent];
     xj ~ normal(e_x, sigma_ex);

     e_z  = B_0s  +
            I_alpha[individual] +
            psi.*I_x[opponent] +
            I_epsilon[opponent];
     z_c ~ normal(e_z, sigma_e);
 }

generated quantities{
real B_0 = B_0s + mean_z;

real<lower=0> Sigma2_intercept;
real<lower=0> Sigma2_epsilon;
real<lower=0> Sigma2_x;
real<lower=0> Sigma2_phi;
real<lower=0> Sigma2_e;
real<lower=0> Sigma2_ex;

real cov_2; //   int-resimpact
real cov_4; //   x-int
real cov_6; //   x-resimpact
real cov_int_phi;

matrix[3, 3]  Omega_I;

Sigma2_intercept=sigma_I[1]^2;
Sigma2_epsilon=sigma_I[2]^2;
Sigma2_x=sigma_I[3]^2;

Sigma2_e=sigma_e^2;
Sigma2_ex=sigma_ex^2;

Omega_I = L * L';
  cov_2 = Omega_I[1,2]*sqrt(Sigma2_epsilon*Sigma2_intercept);
  cov_4 = Omega_I[1,3]*sqrt(Sigma2_x*Sigma2_intercept);
  cov_6 = Omega_I[2,3]*sqrt(Sigma2_x*Sigma2_epsilon);
  cov_int_phi = cov_2 + psi*cov_4;
  Sigma2_phi=sigma_I[2]^2 + psi^2*sigma_I[3]^2;

}

