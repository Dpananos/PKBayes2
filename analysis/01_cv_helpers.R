library(cmdstanr)
library(tidyverse)
library(tidybayes)
library(tidymodels)
library(posterior)
library(Metrics)

# This script creates some tools I use to measure out of sample performance.
# My strategy is to use repeated CV from the {rsample} library
# The functions below resample the new subejcts into test and train, then append the old subjects 
# To the training data.  We need the old subjects in their entirety so we can estimate the within subject variance.

# load in the data from step 00.
# Will have a fixed effect for study so convert it to a numeric.
# Rescale concentrations to be mg/L since priors are on L scale and Drug amount is on mg
concentration_data <- read_csv("data/combined_data.csv") %>% 
  mutate(
    sex = if_else(sex == 'male', 0, 1),
    yobs = yobs/1000,
    study = if_else(study=='old',0,1)
  ) %>% 
  rename(from_new = study)


# Splice out new and old subjects.  Will do vfold_cv on new and append old subjects to the training set.
# Need the old subjects since they will be used to estimate variance of random effects.
new_subjects<-filter(concentration_data, from_new==1)
old_subjects<-filter(concentration_data, from_new==0)


bayes_fit_and_predict<-function(splits, filename){
  
  # splits is an rsample split
  # extract training and just append the old subjects to it.
  # Convert the subject ids to a factor for use in compose_data eventually
  # Converting it to a factor before causes headaches because not all levels are present in the factor.
  train<-analysis(splits) %>% 
         bind_rows(old_subjects) %>% 
         mutate(subjectids = factor(subjectids))
  
  # Since old subjects appear multiple times, we cant let them affect the covariate centering.
  # Drop all repeated people.  Keep one observation.
  unique_subjects<-distinct(train, subjectids, .keep_all = T)
  
  
  # Compute mean and sds of predictors where neccesary.
  covar_mean<- unique_subjects %>% 
                select(age, creatinine, weight) %>% 
                summarise_all(mean) %>% 
                as.numeric
  
  
  covar_sd<- unique_subjects %>% 
            select(age, creatinine, weight) %>% 
            summarise_all(sd) %>% 
            as.numeric
  
  # Now, rescale the training data.
  train[,c('age','creatinine','weight')] = scale(train[,c('age','creatinine','weight')], center = covar_mean, scale = covar_sd)
  
  # Turn the data into a list to pass to cmdstanr
  training_data_list <- compose_data(train)
  
  # Now, we do the same for the test set.
  test<-assessment(splits) %>% 
        mutate(subjectids = factor(subjectids))
  
  test[,c('age','creatinine','weight')] = scale(test[,c('age','creatinine','weight')], center = covar_mean, scale = covar_sd)
  colnames(test) <- str_c('test_', colnames(test))
  
  testing_data_list <- compose_data(test)
  testing_data_list['test_n'] <- testing_data_list$n
  testing_data_list <- within(testing_data_list, rm(n))
  
  # Combine the training and test set to pass to Stan.  Fir the model, extract predictions, compute OOS error.
  combined_data<-c(training_data_list, testing_data_list)
  model<- cmdstan_model(stan_file = 'models/model_for_repeated_cv.stan')
  fit <- model$sample(combined_data, chains = 4, parallel_chains = 4) 
  
  # Only extract CPRED.  Can speed things up a little.
  draws = as_draws_df(fit$draws('CPRED'))
  
  predictions<-draws %>% 
              spread_draws(CPRED[i]) %>% 
              mean_qi
  
  predictions$yobs<-combined_data$test_yobs

  error = rmse(predictions$CPRED, predictions$yobs)
  
  # Stan often causes R to crash for some reason.
  # In case R crashes (fuck me, right?) I'll just write the output to a text file
  fl = paste0('data/bayes_error/',filename,'.txt')
  con = file(fl)
  writeLines(as.character(error), con)
  close(con)
}

linmod_fit_and_predict<-function(splits){
  
  training = analysis(splits)
  val = assessment(splits)
  
  linmod <- lm(log(yobs) ~ time + sex + age + weight + creatinine + D , data = training)
  
  lm_predictions<-predict(linmod)
  sigma = Metrics::rmse(log(training$yobs), lm_predictions)
  lm_predictions<-exp(predict(linmod, newdata = val))*exp(sigma^2/2)
  
  error = rmse(lm_predictions, val$yobs)
  
  return(error)
}
