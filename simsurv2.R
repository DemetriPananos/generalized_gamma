library(tidyverse)
library(tidybayes)
library(cmdstanr)
library(bayesplot)

time <- c(1:7, 7)
censored <- c(rep(0, 7), 1)
wt <- c(131, 126 ,90 ,60 ,42 ,34 ,26, 491)

d <- tibble(
  time_lower = time-1,
  time_upper = time,
  wt = wt,
  censored=censored
)


model <- cmdstan_model(stan_file = 'stan/model_2.stan')

stan_data <- compose_data(d)

pred_time <- seq(0.1,12, 0.1)
nt <- length(pred_time)

stan_data$df <- 30
stan_data$nt <- nt
stan_data$pred_time <- pred_time

fit <- model$sample(stan_data, chains = 10, parallel_chains = 10, adapt_delta=0.99)

fit$cmdstan_diagnose()

np <- nuts_params(fit)

mcmc_pairs(fit$draws(),
           pars = c('mu','sigma','k'),
           transformations = log,
           np=np)


fit %>% 
  spread_draws(survival_function[i]) %>% 
  mutate(time = pred_time[i]) %>% 
  ggplot(aes(time, survival_function)) + 
  stat_lineribbon() + 
  geom_point(
    data= tibble(time=1:7, active = c(863 ,743 ,653 ,593 ,551 ,517 ,491)),
    aes(time, active/1000),
    inherit.aes = F,
    color='black',
    fill='white',
    shape=21,
    size = 3
      ) + 
  theme_classic() + 
  scale_fill_brewer(palette = 'Blues')


fit %>% 
  spread_draws(cumulative_hazard[i]) %>% 
  mutate(time=pred_time[i]) %>% 
  ggplot(aes(time, cumulative_hazard)) + 
  stat_lineribbon() + 
  theme_minimal() +
  scale_fill_brewer(palette = 'Blues')


