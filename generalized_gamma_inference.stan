functions {
  #include generalized_gamma_functions.stan
}
data {
  int n;
  array[n] real x;
}
parameters {
  real mu;
  real<lower=0> sigma;
  real Q;
}
model {
  mu ~ student_t(30, 0, 1);
  sigma ~ gamma(3, 3);
  Q ~ student_t(30, 0, 1);
  
  for(i in 1:n){
    target += generalized_gamma_lpdf(x[i]|mu, sigma, Q);
  }
}
generated quantities {
  array[n] real x_ppc;
  for(i in 1:n) {
    x_ppc[i] = generalized_gamma_rng(mu, sigma, Q);
  }
}
