functions {
  #include generalized_gamma_functions.stan
}
data {
  int n;
  array[n] int time_lower;
  array[n] int time_upper;
  array[n] int wt;
  array[n] int censored;
  
  int nt;
  array[nt] real pred_time;
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
    
      if(censored[i]==1){
        target += wt[i] * lawless_generalized_gamma_lccdf(time_upper[i]| k, mu, sigma);
      } else {
        if(time_lower[i]==0) {
          target += wt[i] * lawless_generalized_gamma_lcdf(time_upper[i]| k, mu, sigma);
        } else{
          target += wt[i] * log_diff_exp(
            lawless_generalized_gamma_lcdf(time_upper[i] | k, mu, sigma),
            lawless_generalized_gamma_lcdf(time_lower[i] | k, mu, sigma)
            );
        }
      }
        
    } 
}
generated quantities {
  array[nt] real pdf;
  array[nt] real cumulative_hazard;
  array[nt] real survival_function;
  
  for(i in 1:nt) {
    pdf[i] = lawless_generalized_gamma_pdf(pred_time[i] , k, mu, sigma);
    survival_function[i] = 1 - lawless_generalized_gamma_cdf(pred_time[i] | k, mu, sigma);
    cumulative_hazard[i] = -log(survival_function[i]);
  }
  
}
