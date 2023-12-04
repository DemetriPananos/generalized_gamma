library(tidyverse)
library(tidybayes)
library(cmdstanr)
library(bayesplot)

n <- 100
raw_time <-rweibull(n,10, 20)

hist(raw_time)

time_lower <- floor(raw_time)
time_upper <- ceiling(raw_time)


d <- tibble(time_lower, time_upper) %>% 
  group_by(time_lower, time_upper) %>% 
  summarise(wt = n())

d

model <- cmdstan_model(stan_file = 'stan/model.stan')

stan_data <- compose_data(d)

time <- seq(0,max(time_upper), 0.1)
nt <- length(time)

stan_data$df <- 35
stan_data$nt <- nt
stan_data$time <- time

fit <- model$sample(stan_data, chains = 4, parallel_chains = 4, adapt_delta = 0.99)

fit$cmdstan_diagnose()

np <- nuts_params(fit)

mcmc_pairs(fit$draws(),
           pars = c('mu','sigma','k'),
           transformations = log,
           np=np)


mcmc_parcoord(fit$draws(),
           pars = c('mu','sigma','k'),
           transformations = log,
           np=np)

fit %>% 
  spread_draws(survival_function[i]) %>% 
  mutate(
    time = time[i]
  ) %>% 
  ggplot(aes(time, survival_function)) + 
  stat_lineribbon()+ 
  stat_function(
    fun=pweibull,
    args=list(shape=10, scale=20, lower.tail=F),
    color='red',
    size=1
  ) +
  scale_fill_brewer(palette = 'Blues')+ 
  theme_minimal()
