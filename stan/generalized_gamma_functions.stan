real generalized_gamma_pdf(real x, real alpha, real beta, real nu) {
  real d;
  d = nu/(beta*tgamma(alpha/nu)) *
    (x/beta)^(alpha-1) * exp(-(x/beta)^nu);
  return d;
}

real generalized_gamma_lpdf(real x, real alpha, real beta, real nu) {
  real d;
  d = log(nu) - log(beta) - lgamma(alpha/nu) +
    (alpha-1)*(log(x) - log(beta)) - (x/beta)^nu;
  return d;
}

real generalized_gamma_cdf(real x, real alpha, real beta, real nu) {
  real d;
  d = gamma_p(alpha/nu, (x/beta)^nu);
  return d;
}

real generalized_gamma_lccdf(real x, real alpha, real beta, real nu) {
  return log_diff_exp(0, generalized_gamma_lcdf(x| alpha, beta, nu));
}


real generalized_gamma_lcdf(real x, real alpha, real beta, real nu) {
  real d;
  d = log(generalized_gamma_cdf(x| alpha, beta, nu));
  return d;
}


real lawless_generalized_gamma_pdf(real x, real k, real mu, real sigma) {
  real w;
  real d;
  w = (log(x)-mu)/sigma;
  d = (k^(k-.5))/sigma/tgamma(k) * exp(
    sqrt(k)*w - k*exp(k^(-.5)*w)
  )/x;
  return d;
}


real lawless_generalized_gamma_lpdf(real x, real k, real mu, real sigma) {
  real y;
  real w;
  real d;
  y = log(x);
  w = (y-mu)/sigma;
  d = (k-.5)*log(k) - log(sigma) - lgamma(k) +
    (sqrt(k)*w - k*exp(1/sqrt(k)*w)) - y;
  return d;
}



real lawless_generalized_gamma_cdf(real x, real k, real mu, real sigma) {
  real w;
  real d;
  w = (log(x) - mu)/sigma;
  d = gamma_p(k, k*exp(1/sqrt(k)*w));
  return d;
}


real lawless_generalized_gamma_lccdf(real x, real k, real mu, real sigma) {
 
 return log_diff_exp(0, lawless_generalized_gamma_lcdf(x| k, mu, sigma));
}

real lawless_generalized_gamma_lcdf(real x, real k, real mu, real sigma) {
  real w;
  real d;
  w = (log(x) - mu)/sigma;
  d = log(gamma_p(k, k*exp(1/sqrt(k)*w)));
  return d;
}