// Implementation stolen from flexsurv
// https://github.com/chjackson/flexsurv/blob/3b76de9481f96f4107124a3e3a427a6766e6b60e/src/gengamma.h

// TODO: Proper warnings for negagive params
real generalized_gamma_lpdf(real x, real mu, real sigma, real Q) {
  
  real lpdf;
  
  if (Q!=0) {
  real y = log(x);
  real w = (y-mu) / sigma;
  real abs_q = abs(Q);
  real qw = Q*w;
  real qi = 1.0 / (Q*Q);
  lpdf = -log(sigma * x) + 
         log(abs_q) * (1 - 2*qi) + 
         qi * (qw - exp(qw)) - 
         lgamma(qi);
  } else {
   lpdf = lognormal_lpdf(x | mu, sigma);
  }
  
  return lpdf;
  
}

// TODO: Proper warnings for negagive params
real generalized_gamma_lcdf(real x, real mu, real sigma, real Q) {
  
  real y = log(x);
  real w = (y-mu) / sigma;
  real lcdf = 0;
  
  if (Q != 0) {
  real qq = 1.0 / (Q*Q);
  real expnu = exp(Q*w)*qq;
  
    if (Q > 0) {
      lcdf = gamma_lcdf(expnu | qq, 1);
    } else {
      lcdf = gamma_lccdf(expnu | qq, 1);
    }
    
  } else {
    lcdf = lognormal_lcdf(x | mu, sigma);
  } 
  
  
  return lcdf;
  
}


// Implementation from flexsurv
// https://github.com/chjackson/flexsurv/blob/3b76de9481f96f4107124a3e3a427a6766e6b60e/R/GenGamma.R
real generalized_gamma_rng(real mu, real sigma, real Q) {
  
  real y;
  real w;
  real qq = Q*Q;
  if (Q!=0) {
    real qi = 1.0/qq;
    real gamma_var = gamma_rng(qi, 1.0);
    w = log(qq * gamma_var) / Q;
    y = exp(mu + w*sigma);
  } else {
    y = lognormal_rng(mu, sigma);
  }
  
  return y;
  
  
}