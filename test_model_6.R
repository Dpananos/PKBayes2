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
       mutate(sex = if_else(sex=='male',1, 0),
              yobs = yobs/1000,
              subjectids = factor(subjectids),
              i = seq_along(subjectids)) 

subjects = distinct(raw_data, subjectids, .keep_all = T) %>% mutate(i = seq_along(subjectids))

old = raw_data %>% 
  filter(study=='old')

old_means = old %>% select(age, weight, creatinine) %>% summarise_all(mean) %>% as.numeric
old_sds = old %>% select(age, weight, creatinine) %>% summarise_all(mean) %>% as.numeric

raw_data[, c('age','weight','creatinine')] = scale(raw_data[, c('age','weight','creatinine')], old_means, old_sds)


ii = sample(289:nrow(raw_data), size = TEST_SIZE)

test =raw_data[ii,] %>% mutate(i = seq_along(i))
train =raw_data[-ii,]


model_data = compose_data(train)
test_model_data = test %>% 
  rename(
    test_time = time,
    test_sex = sex,
    test_weight = weight,
    test_creatinine = creatinine,
    test_age = age,
    test_D = D
  ) %>% 
  compose_data(.n_name = n_prefix('n_test'))


model = cmdstan_model(stan_file = 'model_dev_6.stan')
fit = model$sample(c(model_data, test_model_data), chains = CHAINS, adapt_delta = ADAPT_DELTA)
draws = fit$draws() %>% as_draws_df()


results = draws %>% 
  spread_draws(COOS[i]) %>% 
  mean_qi %>% 
  left_join(test)  


Metrics::mae(results$yobs, results$COOS)

#With linear model
d=filter(train, study=='new')
linmod = lm(log(yobs) ~ age + sex + weight + creatinine + D, data = d)
rmse = Metrics::rmse(log(d$yobs), predict(linmod))
preds = exp(predict(linmod, newdata = d))*exp(rmse^2/2)


plot(test$yobs, exp(predict(linmod, newdata = test)))


test$pred_bayes = results$COOS
test$pred_linmod =  exp(predict(linmod, newdata = test))

Metrics::mae(test$yobs, test$pred_linmod)


test %>% 
  ggplot()+
  geom_point(aes(yobs, pred_bayes), color = 'red')+
  geom_point(aes(yobs, pred_linmod), color = 'black')+
  geom_abline()


test %>% 
  ggplot(aes(pred_bayes, pred_linmod))+
  geom_point()+
  geom_abline()
