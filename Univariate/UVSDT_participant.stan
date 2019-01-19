//UVSDT model for Signal/Noise SDT for 6pt rating scale with participant-random-effects
//Combination of Trippas et al. (2018), Selker et al. (preprint), Pratte & Rouder (2011)
//no item effects

functions {
  vector predictData(real mu, real sigma, vector cr) {
    vector[6] theta;
    theta[1] = normal_cdf(0, mu, sigma);
    theta[2] = normal_cdf(cr[1], mu, sigma) - normal_cdf(0, mu, sigma);
    theta[3] = normal_cdf(cr[2], mu, sigma) - normal_cdf(cr[1], mu, sigma);
    theta[4] = normal_cdf(cr[3], mu, sigma) - normal_cdf(cr[2], mu, sigma);
    theta[5] = normal_cdf(1, mu, sigma) - normal_cdf(cr[3], mu, sigma);
    theta[6] = 1 - normal_cdf(1, mu, sigma);
    return theta;
  }
}

data { 
  int<lower=1> numberParticipants;
  //matrix with participants in rows, frequency of ratings 1-6 in columns
  int<lower=0> testData_signalitems[numberParticipants, 6]; 
  int<lower=0> testData_noiseitems[numberParticipants, 6];
}

transformed data {
  int<lower=1> k;
  int<lower=1> numberParameters;
  numberParameters = 4; //key parameters (mu_s,sigma_s,mu_n,sigma_n); excludes criteria
  k = 6; //rating bins
}

parameters {
  // Initialize parameter space for hierarchical model, top-down
  // base choice: anchor SDT model by criteria and allow mu_s and mu_n to vary
  // alternative: fix mu_n = 0, sigma_n = 1 as typical, allow criteria to vary but then define criteria relative to mu_n 
  // alternative: fix mu_n = 0, https://osf.io/v3b76/ only estimate spread and shift of criteria
  
  
  //Prep hyperparameters
  real grand_mu[2]; //-> for mu_s, mu_n
  real<lower=0> grand_sigma[2]; // -> for sigma_s, sigma_n
  
  // Assumption: no parent distribution for criteria
  // Participant effects: as variance-covariance matrix across parameters
  // participant effect includes correlations of participant effects across parameters
  // Alternative: drop correlations of effect between parameters and introduce separate participant effects for signal and noise items (https://doi.org/10.1016/j.jmp.2010.08.007)
  
  matrix[numberParameters,numberParticipants] alpha_t_subj; // prep for multiplication
  cholesky_factor_corr[numberParameters] Lower_Omega_subj; //for prior chol.corr correlation
  vector<lower=0>[numberParameters] sigma_alpha_subj; //for participant effect (w/o corr)
  
  // Prep criteria prior parameters
  ordered[3] mu_crits; //-> for participant-specific mu_crits for crit = 2,3,4, crit1 = 0, crit5 = 1
  real<lower=0> sigma_crits[2]; // -> for participant-specific sigma_crits for crit2, crit4 with crit3 reasonably unvariable
  
  // Prep mid-criteria on unit scale, crit = 2,3,4
  ordered[3] crit_un[numberParticipants];

} 

transformed parameters {
  
  // Apply simplex algoritm to predicted data
  simplex[k] theta_signalitems[numberParticipants];
  simplex[k] theta_noiseitems[numberParticipants];

  // Prep matrix for crit(2,3,4) for all participants
  ordered[3] crit[numberParticipants];
  
  // prep alpha
  matrix[numberParticipants,numberParameters] alpha_subj; 
  // participant effect includes correlations of participan effect across parameters
  alpha_subj  = (diag_pre_multiply(sigma_alpha_subj, Lower_Omega_subj) * alpha_t_subj)'; 

  
  for (j in 1:numberParticipants) {
    
    //Transform unit scale criteria
    crit[j, 1] = Phi(crit_un[j,1]);
    crit[j, 2] = Phi(crit_un[j,2]);
    crit[j, 3] = Phi(crit_un[j,3]);

    // add individual displacement to group level parameters:
    theta_signalitems[j] = predictData(
        grand_mu[1] + alpha_subj[j,1], 
        exp(log(grand_sigma[1]) + alpha_subj[j,2]), 
        crit[j]
    );

    theta_noiseitems[j] = predictData(
        grand_mu[2] + alpha_subj[j,3] , 
        exp(log(grand_sigma[2]) + alpha_subj[j,4]), 
        crit[j]
    );
  }
}

model {
  // Set the priors for hyperparameters
  // grand mu and sigma
  to_vector(grand_mu) ~ cauchy(0.5, 4); 
  to_vector(grand_sigma) ~  cauchy(1, 10);
  
  // participant effects
  Lower_Omega_subj ~ lkj_corr_cholesky(1); 
  to_vector(alpha_t_subj) ~ normal(0, 1); 
  sigma_alpha_subj ~ cauchy(0, 4);

  // criteria
  mu_crits ~ normal(0, 1); //(2,3,4) # normals transform into flat uniform in uniform space
  sigma_crits ~ cauchy(0, 4); //(2,4)

  for (m in 1:3) { //crit2,3,4, crit1 = 0, crit5 = 1
    if (m < 2) crit_un[,m] ~ normal(mu_crits[m], sigma_crits[m]);
    if (m == 2) crit_un[,m] ~ normal(mu_crits[m], 0.1); 
    if (m > 2) crit_un[,m] ~ normal(mu_crits[m], sigma_crits[m-1]);
  }
 
  
  // Observed counts (binned freq rather than trial-by-trial basis)
   for (j in 1:numberParticipants) {
    testData_signalitems[j,1:k] ~ multinomial(theta_signalitems[j]); 
    testData_noiseitems[j,1:k] ~ multinomial(theta_noiseitems[j]);
  }
}

generated quantities {
  corr_matrix[numberParameters] Omega_subj;
  real dprime;
  real bias;
  real sigma_ratio;
  // ^ if these are init'ed lower down, they throw parsing errors
  
  // correlations of participant effects between parameters: LL'
  Omega_subj = multiply_lower_tri_self_transpose(Lower_Omega_subj);
  
  // generic output
  dprime = grand_mu[1] - grand_mu[2];
  bias = (grand_mu[1] + grand_mu[2])/2;
  sigma_ratio = grand_sigma[1]/grand_sigma[2];
  

}
