

model <- cmdstan_model(stan_file = 'models/big_gamma.stan')
model_data<-combined_data<-readRDS('fit_data/combined_Data.RDS')

fit<-model$sample(data = model_data, parallel_chains = 4, adapt_delta = .99)
draws = as_draws_df(fit$draws())

d = tibble(subjectids = model_data$subjectids,
           time = model_data$time,
           yobs = model_data$yobs) %>% 
    mutate(i = seq_along(subjectids))



draws %>% 
  spread_draws(Cppc[i]) %>% 
  mean_qi() %>% 
  bind_cols(d) %>% 
  filter(subjectids<=36) %>% 
  ggplot(aes(time,Cppc))+
  geom_line()+
  geom_ribbon(aes(ymin = .lower, ymax = .upper))+
  geom_point(aes(time, yobs), inherit.aes = F, color = 'red')+
  facet_wrap(~subjectids)



draws %>% 
  spread_draws(Cppc[i], n = 10) %>% 
  left_join(d) %>% 
  filter(subjectids<=36) %>% 
  ggplot(aes(time,Cppc))+
  geom_line(aes(group = .draw))+
  facet_wrap(~subjectids)
