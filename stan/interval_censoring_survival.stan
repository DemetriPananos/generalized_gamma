functions {
  #include generalized_gamma_functions.stan
}
data {
  int n;
  int n_cohort;
  array[n] int cohort;
  array[n] int time_lower;
  array[n] int time_upper;
  array[n] int wt;
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
          // If the lower time is the largest it can be (and there is no upper time as a consequence)
          // Increment the target by the log probability that the survival time is above truncation_time
          target += wt[i] * lawless_generalized_gamma_lccdf(time_lower[i]| k[cohort[i]], mu[cohort[i]], sigma[cohort[i]]);
        
      } else {
        if(time_lower[i]==0) {
          // time_lower == 0 is a special case since.  At time=0, the cdf is 0.  Can't take the log of this
          // so the log_diff_exp won't work.  Simply increment by the log cdf.
          target += wt[i] * lawless_generalized_gamma_lcdf(time_upper[i]| k[cohort[i]], mu[cohort[i]], sigma[cohort[i]]);
          
        } else{
          
          target += wt[i] * log_diff_exp(
            lawless_generalized_gamma_lcdf(time_upper[i] |k[cohort[i]], mu[cohort[i]], sigma[cohort[i]]),
            lawless_generalized_gamma_lcdf(time_lower[i] | k[cohort[i]], mu[cohort[i]], sigma[cohort[i]])
            );
            
        }
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

