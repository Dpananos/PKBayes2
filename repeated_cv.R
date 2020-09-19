library(cmdstanr)
library(posterior)
library(tidybayes)
library(tidyverse)
library(patchwork)
library(table1)
library(lme4)
source("analysis/02_prepare_data.R")

do_it <-function(){
    
  CHAINS = 4
  
  data<-make_split(100)
  model = cmdstan_model(stan_file = 'models/big_model.stan')
  fit <- model$sample(data$model_data, chains = CHAINS, parallel_chains = 4)
  
  
  
  draws = as_draws_df(fit$draws('CPRED'))
  
  predictions<-draws %>% 
    spread_draws(CPRED[i]) %>% 
    mean_qi
  predictions$yobs<-combined_data$test_yobs
  
  a = Metrics::rmse(predictions$CPRED, predictions$yobs)
  
  linmod_data <- data$train_data %>% filter(from_old==1)
  linmod_test_data <- data$test_data
  colnames(linmod_test_data) = str_replace(colnames(linmod_test_data), 'test_','')
  
  linmod <- lm(log(yobs) ~ time + sex + age + weight + creatinine + D , data = linmod_data)
  
  lm_predictions<-predict(linmod)
  sigma = Metrics::rmse( log(linmod_data$yobs), lm_predictions)
  
  lm_predictions<-exp(predict(linmod, newdata = linmod_test_data))*exp(sigma^2/2)
  
  b = Metrics::rmse(lm_predictions, linmod_test_data$yobs)
  
  tibble(bayes_error = a, linmod_error = b)
  
}



results<-rerun(1000, do_it()) %>% 
          map_dfr(~ .x)
  