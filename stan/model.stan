functions {
  #include generalized_gamma_functions.stan
}
data {
  int n;
  array[n] int time_lower;
  array[n] int time_upper;
  array[n] int wt;
  
  int nt;
  array[nt] real time;
  real df;
  
}

parameters {
  real mu;
  real<lower=0> k;
  real<lower=0> sigma;
}
model {
  mu ~ student_t(df, 0, 1);
  sigma ~ cauchy(0, 1);
  k ~  student_t(df, 0, 1);
  
  for (i in 1:n){
      
      // The if time_lower == 0 then log(cdf(0)) is undefined since 
      // cdf(0) = 0 for survival time distributions
      // Instead, if time == 0, then use cdf time_upper.
      // This will get the probability of seeing a survival between 0 and time_upper
      if(time_lower[i]==0) {
        
        target += wt[i] * lawless_generalized_gamma_lcdf(time_upper[i]| k, mu, sigma);
        
      } else{
        
        // The probability of dying between time_lower and time_upper
        // is given by the difference in cdfs.  Since Stan works on the log scale,
        // Use log_diff_exp to take the difference of the log_cdfs in exp space
        // and return to log space.
        
        target += wt[i] * log_diff_exp(
          lawless_generalized_gamma_lcdf(time_upper[i] | k, mu, sigma),
          lawless_generalized_gamma_lcdf(time_lower[i] | k, mu, sigma)
          );
          
      }
        
    } 
}
generated quantities {
  array[nt] real cumulative_hazard;
  array[nt] real survival_function;
  
  for(i in 1:nt) {
    survival_function[i] = 1 - lawless_generalized_gamma_cdf(time[i] | k, mu, sigma);
    cumulative_hazard[i] = -log(survival_function[i]);
  }
  
}
