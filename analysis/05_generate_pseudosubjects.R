library(cmdstanr)


X = matrix(rep(0,10), nrow = 2)
X[2, 5] = 1

model_data = list(n_subjects = 2, 
                  X = X)


model_data$t_pred = seq(1, 168, 1)
model_data$n_pred = length(model_data$t_pred)
model_data$dose_times = seq(0, 168, 12)
model_data$doses = rep(2.5, length(model_data$dose_times))
model_data$n_doses = length(model_data$dose_times)


model = cmdstan_model('models/generate_repeated_dose_patients.stan')


fit = model$sample(model_data, fixed_param = T, iter_sampling = 1)

tframe = tibble(time = model_data$t_pred) %>% 
         mutate(j = seq_along(time))

fit$draws('C') %>% 
  as_draws_df() %>% 
  spread_draws(C[i,j]) %>% 
  left_join(tframe) %>% 
  ggplot(aes(time/24, C, color = factor(i)))+
  geom_line(aes(group = interaction(i,.draw)))
