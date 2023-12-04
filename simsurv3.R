library(tidyverse)
library(tidybayes)
library(cmdstanr)
library(bayesplot)

n <- 60
raw_time <-rgamma(n, 3, 1)

hist(raw_time)

time_lower <- floor(raw_time)
time_upper <- ceiling(raw_time)


d <- tibble(time_lower, time_upper) %>% 
  group_by(time_lower, time_upper) %>% 
  summarise(wt = n()) %>% 
  ungroup %>% 
  mutate(
    time_mid = (time_lower + time_upper)/2,
    pct_surviving = 1 - cumsum(wt)/sum(wt)
  )


train_test_split <- function(d, cut_time) {
  
  
  total_units <- sum(d$wt)
  
  dtrain <- filter(d, time_upper <= cut_time) %>% 
            mutate(censored=0)
  dtest <- filter(d, time_upper > cut_time)
  
  censored_info <- tibble(
    time_lower = cut_time - 1,
    time_upper = cut_time,
    censored=1,
    wt = total_units - sum(dtrain$wt)
  )
  
  stan_data <- bind_rows(dtrain, censored_info) %>% 
               select(time_lower, time_upper, wt, censored) %>% 
               compose_data()
  
  pred_time <- seq(0.1,max(dtest$time_upper), 0.1)
  nt <- length(pred_time)
  
  stan_data$df <- 30
  stan_data$nt <- nt
  stan_data$pred_time <- pred_time
  
  list(train=dtrain, test = dtest, stan_data=stan_data)
}


data <- train_test_split(d, 4)
model <- cmdstan_model(stan_file = 'stan/model_2.stan')

fit <- model$sample(data$stan_data, parallel_chains = parallel::detectCores(), adapt_delta=0.99)


fit %>% 
  spread_draws(survival_function[i]) %>% 
  mutate(
    time = data$stan_data$pred_time[i]
  ) %>% 
  ggplot(aes(time, survival_function)) + 
  stat_lineribbon() + 
  scale_fill_brewer(palette = 'Blues') + 
  geom_point(
    data=data$train, 
    aes(time_upper, pct_surviving),
    color='red', 
    inherit.aes = F
  ) + 
  geom_point(
    data=data$test, 
    aes(time_upper, pct_surviving),
    color='red', 
    inherit.aes = F,
    shape = 21
  ) + 
  stat_function(
    fun=pgamma,
    args = list(shape=3, rate=1, lower.tail=F),
    color='red'
  )

