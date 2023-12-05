functions {
  #include generalized_gamma_functions.stan
}
data {
  int n;
  int n_cohort;
  array[n] int cohort;
  array[n] int time;
  array[n] int censored;
  
  // Hyperparameters for priors.
  real<lower=0> mu_df;
  real<lower=0> mu_loc;
  real<lower=0> mu_scale;
  
  real<lower=0> sigma_df;
  real<lower=0> sigma_loc;
  real<lower=0> sigma_scale;
  
  real<lower=0> k_df;
  real<lower=0> k_loc;
  real<lower=0> k_scale;
  
  // Data for predicting the survival curves
  int nt;
  array[nt] real pred_time;
  
}
parameters {
  vector<lower=0>[n_cohort] mu;
  vector<lower=0>[n_cohort] k;
  vector<lower=0>[n_cohort] sigma;
}
model {
  // TODO: Vectorize pdfs, cdfs, ccdfs
  mu ~ student_t(mu_df, mu_loc, mu_scale);
  sigma ~ student_t(sigma_df, sigma_loc, sigma_scale);
  k ~ student_t(k_df, k_loc, k_scale);
  
  for (i in 1:n){
    
      if(censored[i]==1){
          // If the time is the largest it can be (i.e. is truncation time)
          // Increment the target by the log probability that the survival time is above truncation_time
          target += lawless_generalized_gamma_lccdf(time[i]| k[cohort[i]], mu[cohort[i]], sigma[cohort[i]]);
        
      } else {
        
          target += lawless_generalized_gamma_lpdf(time[i] |k[cohort[i]], mu[cohort[i]], sigma[cohort[i]]);
          
        }
      }
        
} 
generated quantities {
  array[nt, n_cohort] real survival_function;
  
  for(i_time in 1:nt) {
    for(j_cohort in 1:n_cohort) {
      
      survival_function[i_time, j_cohort] = 1 - lawless_generalized_gamma_cdf(pred_time[i_time] | k[j_cohort], mu[j_cohort], sigma[j_cohort]);
      
    }
  }
  
}

