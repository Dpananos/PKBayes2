library(cmdstanr)
library(tidyverse)
library(tidybayes)
library(tidymodels)
library(posterior)
library(Metrics)

concentration_data <- read_csv("data/combined_data.csv") %>% 
  mutate(
    sex = if_else(sex == 'male', 0, 1),
    yobs = yobs/1000, #Convert from ng/ml to mg/L
    study = if_else(study=='old',0,1)
  ) %>% 
  rename(from_new = study)


new_subjects<-filter(concentration_data, from_new==1)
old_subjects<-filter(concentration_data, from_new==0)




bayes_fit_and_predict<-function(splits){
  
  #Compute how to center the predictors
  train<-old_subjects %>% 
         bind_rows(analysis(splits)) %>% 
         mutate(subjectids = factor(subjectids))
  
  unique_subjects<-distinct(train, subjectids, .keep_all = T)
  
  covar_mean<- unique_subjects %>% 
                select(age, creatinine, weight) %>% 
                summarise_all(mean) %>% 
                as.numeric
  
  
  covar_sd<- unique_subjects %>% 
    select(age, creatinine, weight) %>% 
    summarise_all(sd) %>% 
    as.numeric
  
  train[,c('age','creatinine','weight')] = scale(train[,c('age','creatinine','weight')], center = covar_mean, scale = covar_sd)
  
  training_data_list <- compose_data(train)
  
  test<-assessment(splits) %>% 
        mutate(subjectids = factor(subjectids))
  
  test[,c('age','creatinine','weight')] = scale(test[,c('age','creatinine','weight')], center = covar_mean, scale = covar_sd)
  colnames(test) <- str_c('test_', colnames(test))
  
  testing_data_list <- compose_data(test)
  testing_data_list['test_n'] <- testing_data_list$n
  testing_data_list <- within(testing_data_list, rm(n))
  
  
  combined_data<-c(training_data_list, testing_data_list)
  model<- cmdstan_model(stan_file = 'models/big_model.stan')
  fit <- model$sample(combined_data, chains = 4, parallel_chains = 4) 
  
  draws = as_draws_df(fit$draws())
  
  
  predictions<-draws %>% 
    spread_draws(CPRED[i]) %>% 
    mean_qi
  predictions$yobs<-combined_data$test_yobs

  error = rmse(predictions$CPRED, predictions$yobs)
  
  return(error)
}



