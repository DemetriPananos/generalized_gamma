library(cmdstanr)
library(tidyverse)
library(tidybayes)
library(survival)

model_data <- veteran %>% 
              select(trt, time, status) %>% 
              rename(censored=status,
                     cohort = trt) %>% 
              mutate(
                censored = as.integer(abs(censored-1)),
                cohort = factor(cohort)
              )

stan_data <- compose_data(model_data)

stan_data$mu_df <- 30
stan_data$mu_loc <- 0
stan_data$mu_scale <- 1


# Encourage an exponential fit.
stan_data$sigma_df <- 30
stan_data$sigma_loc <- 1
stan_data$sigma_scale <- 1


stan_data$k_df <- 30
stan_data$k_loc <- 1
stan_data$k_scale <- 1


stan_data$nt <- 100
stan_data$pred_time <- 250 * ppoints(stan_data$nt)

model <- cmdstanr::cmdstan_model('stan/right_censoring_survival.stan')

fit <- model$sample(stan_data, parallel_chains = 4, adapt_delta = 0.99)

np <- nuts_params(fit)

mcmc_parcoord(fit$draws(),
              regex_pars = c('mu','k','sigma'), 
              np=np, 
              transformations = )


pred <- fit %>% 
  recover_types(model_data) %>% 
  spread_draws(survival_function[i, cohort], ndraws = 500) %>% 
  mutate(time = stan_data$pred_time[i]) %>% 
  group_by(time, cohort) %>% 
  mean_qi %>% 
  mutate(cohort = str_c('cohort_',cohort)) %>% 
  select(time, cohort, survival_function) %>% 
  spread(cohort, survival_function)


plot(
  pred$time, 
  pred$cohort_1, 
  type='l',
  col='red'
)

lines(
  pred$time, 
  pred$cohort_2, 
  col='blue'
)

# plot KM over the survival curves
km_1 <- survfit(Surv(time, status) ~ strata(trt), data=filter(veteran, trt==1))
km_2 <- survfit(Surv(time, status) ~ strata(trt), data=filter(veteran, trt==2))


lines(km_1$time, km_1$surv, col='red', lty=2)
lines(km_2$time, km_2$surv, col='blue', lty=2)


