library(cmdstanr)
library(tidyverse)
library(tidybayes)
library(posterior)

X = matrix(rep(0,5), nrow = 1)


model_data = list(n_subjects = 1, 
                  X = X)


model_data$t_pred = c(
  runif(1, 24, 36),
  runif(1, 154, 168)
)
model_data$n_pred = length(model_data$t_pred)
model_data$dose_times = seq(0, 168, 12)
model_data$doses = rep(5.0, length(model_data$dose_times))
model_data$n_doses = length(model_data$dose_times)


model = cmdstan_model('models/generate_repeated_dose_patients.stan')


fit = model$sample(model_data, fixed_param = T, iter_sampling = 1)

tframe = tibble(time = model_data$t_pred) %>% 
         mutate(j = seq_along(time))

fit$draws('C_obs') %>% 
  as_draws_df() %>% 
  spread_draws(C_obs[i,j]) %>% 
  left_join(tframe)->f

f %>% 
  ggplot(aes(time/24, C_obs, color = factor(i)))+
  geom_line(aes(group = interaction(i,.draw)))


model2 = cmdstan_model('Untitled.stan')



model2_data = list(
  X = as.numeric(X),
  n = nrow(f),
  C_obs = f$C_obs,
  times = f$time,
  n_doses = model_data$n_doses,
  doses = model_data$doses,
  dose_times = model_data$dose_times,
  tpred = seq(0, 168, 0.05),
  nt = length(seq(0, 168, 0.05))
)

tframe = tibble(t = model2_data$tpred) %>% 
         mutate(i = seq_along(t))
fit2 = model2$sample(model2_data, chains=4, parallel_chains = 4, adapt_delta = 0.99)


fit2$draws('ypred') %>% 
  as_draws_df() %>% 
  spread_draws(ypred[i], n = 200) %>% 
  left_join(tframe) %>% 
  ggplot(aes(t, ypred))+
  geom_line(aes(group = .draw), alpha = 0.1)+
  geom_point(
    data = tibble(t = model2_data$times, 
                  ypred = model2_data$C_obs)
  )
