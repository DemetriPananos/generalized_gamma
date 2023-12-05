library(tidyverse)
library(cmdstanr)
library(tidybayes)
library(bayesplot)


# Data from 
# “How to project customer retention.” 
# Journal of Interactive Marketing 21.1 (2007): 76-90.
otime <- 0:12
regular_surv <- c(1000, 631, 468, 382, 326, 289, 262, 241, 223, 207, 194, 183, 173)
highend_surv <- c(1000, 869, 743, 653, 593, 551, 517, 491, 468, 445, 427, 409, 394)

# Construct a dataframe containing how many users churned in each period.
time <- c(1:11, 11)
regular_wt <- regular_surv[1:11] - regular_surv[2:12]
regular_wt <- c(regular_wt, 1000 - sum(regular_wt))
highend_wt <- highend_surv[1:11] - highend_surv[2:12]
highend_wt <- c(highend_wt, 1000 - sum(highend_wt))
censored <- c(rep(0, length(regular_wt)-1), 1)

model_data <- tibble(time_upper=time, censored, regular=regular_wt, highend=highend_wt) %>% 
  pivot_longer(regular:highend, names_to = 'cohort', values_to = 'wt') %>% 
  mutate(
    cohort=factor(cohort),
    time_lower=time_upper-1
  )


model <- cmdstan_model('stan/interval_censoring_survival.stan')

stan_data <- compose_data(model_data)

stan_data$mu_df <- 30
stan_data$mu_loc <- 0
stan_data$mu_scale <- 1

stan_data$sigma_df <- 1
stan_data$sigma_loc <- 0
stan_data$sigma_scale <- 1


stan_data$k_df <- 10
stan_data$k_loc <- 0
stan_data$k_scale <- 1


stan_data$nt <- 100
stan_data$pred_time <- 12 * ppoints(stan_data$nt)

fit <- model$sample(stan_data, parallel_chains = 4, chains = 4)

fit$cmdstan_diagnose()


pred <- fit %>% 
        recover_types(model_data) %>% 
        spread_draws(survival_function[i, cohort]) %>% 
        mutate(time=stan_data$pred_time[i]) %>% 
        mean_qi

regular_pred <- filter(pred, cohort=='regular')
highend_pred <- filter(pred, cohort=='highend')

plot(otime, regular_surv, col='red', pch=19, ylim = c(0, 1000))
points(otime, highend_surv, col='blue', pch=19)

lines(regular_pred$time, 1000 * regular_pred$survival_function, col='red')
lines(regular_pred$time, 1000 * regular_pred$survival_function.lower, col='red', lty=2)
lines(regular_pred$time, 1000 * regular_pred$survival_function.upper, col='red', lty=2)


lines(highend_pred$time, 1000 * highend_pred$survival_function, col='blue')
lines(highend_pred$time, 1000 * highend_pred$survival_function.lower, col='blue', lty=2)
lines(highend_pred$time, 1000 * highend_pred$survival_function.upper, col='blue', lty=2)
