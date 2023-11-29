library(tidyverse)
library(cmdstanr)
library(flexsurv)
library(posterior)


n <- 50
mu <- log(20)
sigma <- 0.1
Q <- 2
nn <- 10000
x <- if_else(Q<0, 50, 25) * ppoints(n)

stan_data <- list(n=n, x=x, mu=mu, sigma=sigma,Q=Q, nn=nn)
model <- cmdstan_model('generalized_gamma.stan')
fit <- model$sample(stan_data, 
                    fixed_param = T, 
                    parallel_chains = 1, 
                    chains = 1,
                    iter_warmup = 0,
                    iter_sampling = 1)


fit$draws('log_pdf')


# Note that Stan oeprates on the log scale
target_log_pdf <- dgengamma(x, mu=mu, sigma=sigma, Q=Q, log=TRUE)
target_log_cdf <- pgengamma(x, mu=mu, sigma=sigma, Q=Q, log=TRUE)

stan_log_pdf <- fit$draws('log_pdf') %>% 
                as_draws_matrix()

stan_log_cdf <- fit$draws('log_cdf') %>% 
  as_draws_matrix()

gg_draws <- fit$draws('gg_draws')

par(mfrow = c(3, 2))

plot(x, target_log_pdf, type='l', main='Log Density')
points(x, stan_log_pdf, col='red')


plot(x, target_log_cdf, type='l', main='Log Cumulative Density')
points(x, stan_log_cdf, col='red')



plot(x, exp(target_log_pdf), type='l', main='Density')
points(x, exp(stan_log_pdf), col='red')


plot(x, exp(target_log_cdf), type='l', main='Cumulative Density')
points(x, exp(stan_log_cdf) , col='red')



z <- seq(0, max(gg_draws), length.out = 1000)
plot(z, dgengamma(z, mu=mu, sigma=sigma, Q=Q, log=FALSE), type='l')
hist(gg_draws, probability = T, add=T)
lines(z, dgengamma(z, mu=mu, sigma=sigma, Q=Q, log=FALSE), type='l')


plot(ecdf(gg_draws))
lines(z, pgengamma(z, mu=mu, sigma=sigma, Q=Q, log=FALSE),col='red')



