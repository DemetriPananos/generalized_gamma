library(tidyverse)
library(cmdstanr)
library(tidybayes)
library(bayesplot)
time <- 1:12
regular <- c(1000, 631, 468, 382, 326, 289, 262, 241, 223, 207, 194, 183, 173)
highend <- c(1000, 869, 743, 653, 593, 551, 517, 491, 468, 445, 427, 409, 394)

d <- tibble(
  time_lower=c(time-1, max(time-1)),
  time_upper=c(time, max(time)),
  censored = c(rep(0, length(time)), 1),
  regular=c(abs(na.omit(diff(regular))), min(regular) ),
  highend=c(abs(na.omit(diff(highend))), min(highend))
)  %>% 
  pivot_longer(regular:highend, names_to = 'grp', values_to = 'wt') %>% 
  mutate(i = seq_along(grp),
         grp = factor(grp)) %>% 
  group_by(grp) %>% 
  arrange(grp, time_lower, censored) %>% 
  mutate(
    y = cumsum(wt)
  )

d

model <- cmdstan_model('stan/model_3.stan')

stan_data <- compose_data(d)
stan_data$df <- 15 
stan_data$nt <- 250
stan_data$pred_time <- 12*ppoints(stan_data$nt)


fit <- model$sample(stan_data, adapt_delta = 0.99, parallel_chains = 4)

np <- nuts_params(fit)

fit %>% 
  recover_types(d) %>% 
  spread_draws(survival_function[i, grp], ndraws = 100) %>% 
  mutate(
    time = stan_data$pred_time[i]
  ) %>% 
  ggplot(aes(time, survival_function, group=grp)) + 
  stat_lineribbon(point_interval = mean_qi)


fit %>% 
  recover_types(d) %>% 
  spread_draws(cumulative_incidence[i, grp], ndraws = 100) %>% 
  mutate(
    time = stan_data$pred_time[i]
  ) %>% 
  ggplot(aes(time, cumulative_incidence, group=grp)) + 
  stat_lineribbon(point_interval = mean_qi) +
geom_point(
    data=filter(d, censored<1), 
    aes(time_upper, y/1000),
    inherit.aes = F, 
    color='red',
    size=2
  ) + 
  theme(aspect.ratio = 1) + 
  scale_fill_brewer(palette = 'Blues')
