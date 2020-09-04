library(cmdstanr)
library(tidyverse)
library(tidybayes)  
library(posterior)
library(bayesplot)
library(patchwork)
library(lme4)


TEST_SIZE = 100
CHAINS = 4
ADAPT_DELTA = .99



raw_data = read_csv('data/combined_data.csv') %>% 
  drop_na() %>% 
  filter(study == 'new') %>% 
  mutate(sex = if_else(sex=='male',1, 0),
         yobs = yobs/1000,
         subjectids = factor(subjectids),
         i = seq_along(subjectids)) %>% 
  mutate_at(vars(age, weight, creatinine), ~as.vector(scale(.)))



model_data = compose_data(raw_data)


model = cmdstan_model(stan_file = 'model_dev_7.stan')

fit = model$sample(model_data)

linmod = lm(log(yobs) ~ age + sex + weight + creatinine + D, data = raw_data)
summary(linmod)

rmse = Metrics::rmse(log(raw_data$yobs), predict(linmod))
f = raw_data %>% modelr::add_predictions(linmod) %>% mutate(pred = exp(pred)*exp(rmse^2/2))


f = fit$draws('C') %>% 
  as_draws_df() %>% 
  spread_draws(C[i]) %>% 
  mean_qi %>% 
  bind_cols(f) 


Metrics::rmse(f$yobs, f$C)

Metrics::rmse(f$yobs, f$pred)
