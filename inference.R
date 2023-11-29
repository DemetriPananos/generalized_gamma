library(tidyverse)
library(cmdstanr)
library(flexsurv)
library(posterior)
library(tidybayes)
# Pull From The Prior -----------------------------------------------------

mu <- rt(1, 3.5)
sigma <- rgamma(1,  3, 3)
Q <- rt(1, 3.5)

n <- 100
x <- rgengamma(n, mu, sigma, Q)

stan_data <- list(x=x, n=n)
model <- cmdstan_model('generalized_gamma_inference.stan')
fit <- model$sample(stan_data, 
                    parallel_chains = parallel::detectCores(),
                    chains = 12)

fit$cmdstan_diagnose()


draws <- fit %>% 
         spread_draws(mu, sigma, Q)

par(mfrow = c(2, 2))
hist(draws$mu, main=expression(mu))
abline(v=mu, col='red')

hist(draws$sigma, main=expression(sigma))
abline(v=sigma, col='red')


hist(draws$Q, main=expression(Q))
abline(v=Q, col='red')


# Check for bias in estimates ---------------------------------------------

sim_fit_check <- function(i, model){
  mu <- rt(1, 30)
  sigma <- rgamma(1, 3, 3)
  Q <- rt(1, 30)
  
  truth <- tibble(
    .variable = c('Q','sigma','mu'),
    truth = c(Q, sigma, mu)
  )
  
  n <- 100
  x <- rgengamma(n, mu, sigma, Q)
  
  stan_data <- list(x=x, n=n)
  fit <- model$sample(stan_data, 
                      parallel_chains = parallel::detectCores(),
                      chains = 12)
  
  diags = fit$sampler_diagnostics()
  is_divergent = sum(diags[, , 'divergent__'])>0
  
  draws <- fit %>% 
           gather_draws(mu, sigma, Q) %>% 
           mean_qi %>% 
           inner_join(truth) %>% 
           mutate(is_divergent)
  
  dir_loc = glue::glue('simulation_results/results_{i}')
  dir.create(dir_loc, showWarnings = FALSE)
  file_loc = file.path(dir_loc, 'results.rds')
  write_rds(draws, file=file_loc)
}

results <- map(1:50, ~sim_fit_check(i=.x, model=model), .id = 'simulation')




# Crawl results and plot --------------------------------------------------

files <- list.files('simulation_results/', recursive = T)
d <- map_dfr(files, ~read_rds(file.path('simulation_results/', .x)))


d %>% 
  ggplot(aes(truth, .value, ymin=.lower, ymax=.upper, color=is_divergent)) + 
  geom_pointrange() + 
  geom_abline() + 
  facet_wrap(~.variable, scales = 'free', nrow = 2)


