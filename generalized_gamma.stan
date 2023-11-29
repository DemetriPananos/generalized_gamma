functions {
  #include generalized_gamma_functions.stan
}
data {
  int n;
  int nn;
  array[n] real x;
  real Q;
  real mu;
  real sigma;
}
parameters {}
model {}
generated quantities {
  array[n] real log_pdf;
  array[n] real log_cdf;
  array[nn] real gg_draws;
  
  for(i in 1:n) {
    log_pdf[i] = generalized_gamma_lpdf(x[i]|mu, sigma, Q);
    log_cdf[i] = generalized_gamma_lcdf(x[i]|mu, sigma, Q);
  }
  
  for(i in 1:nn) {
    gg_draws[i] = generalized_gamma_rng(mu, sigma, Q);
  }
}
