functions {
  #include generalized_gamma_functions.stan
}
data {
  int n;
  int n_grp;
  array[n] int grp;
  array[n] int time_lower;
  array[n] int time_upper;
  array[n] int wt;
  array[n] int censored;
  
  real df;
  int nt;
  array[nt] real pred_time;
  
}
parameters {
  vector<lower=0>[n_grp] mu;
  vector<lower=0>[n_grp] k;
  vector<lower=0>[n_grp] sigma;
}
model {
  mu ~ student_t(df, 0, 1);
  sigma ~ cauchy(0, 1);
  k ~  student_t(df, 0, 1);
  
  for (i in 1:n){
    
      if(censored[i]==1){
        target += wt[i] * lawless_generalized_gamma_lccdf(time_upper[i]| k[grp[i]], mu[grp[i]], sigma[grp[i]]);
      } else {
        if(time_lower[i]==0) {
          target += wt[i] * lawless_generalized_gamma_lcdf(time_upper[i]| k[grp[i]], mu[grp[i]], sigma[grp[i]]);
        } else{
          target += wt[i] * log_diff_exp(
            lawless_generalized_gamma_lcdf(time_upper[i] |k[grp[i]], mu[grp[i]], sigma[grp[i]]),
            lawless_generalized_gamma_lcdf(time_lower[i] | k[grp[i]], mu[grp[i]], sigma[grp[i]])
            );
        }
      }
        
    } 
}
generated quantities {
  array[nt, n_grp] real survival_function;
  
  for(i in 1:nt) {
    for(j in 1:n_grp) {
      
      survival_function[i, j] = 1 - lawless_generalized_gamma_cdf(pred_time[i] | k[j], mu[j], sigma[j]);
      
    }
  }
  
}

