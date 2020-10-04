library(tidyverse)
library(cmdstanr)
library(tidybayes)
library(posterior)
library(rstan)

# This script fits the full model.

#Load in the data and turn categorical variables into indicators. 
#Rescale outcome from ng/ml to mg/L
concentration_data <- read_csv("data/combined_data.csv") %>% 
  mutate(
    sex = if_else(sex == 'male', 0, 1),
    yobs = yobs/1000,
    study = if_else(study=='old',0,1)
  ) %>% 
  rename(from_new = study)

# Standardize continuous variables for each patient.
# Model accepts covariates for patients seperately from observed outcome data
# Compose two different dataframes (one for covs, another for obs data) seperately
subject_data = concentration_data %>%
  distinct(subjectids, .keep_all = T) %>% 
  select(sex, weight, creatinine, age, D, from_new) %>% 
  compose_data()

# Here is the obs data
obs_data = concentration_data %>% 
  select(yobs, time, subjectids) %>% 
  compose_data()

subject_data$n_subjectids = subject_data$n
subject_data = within(subject_data, rm(n))

model_data = c(subject_data, obs_data)

model_data$t_pred = seq(1, 168, 1)
model_data$n_pred = length(model_data$t_pred)
model_data$dose_times = seq(0, 168, 12)
model_data$doses = rep(2.5, length(model_data$dose_times))
model_data$n_doses = length(model_data$dose_times)



model = cmdstan_model(stan_file = 'models/model_for_all_data.stan')

fit = model$sample(model_data, parallel_chains = 4)

library(fitdistrplus)


fit$draws('beta_a') %>% 
  as_draws_df() %>% 
  gather_draws(beta_a[i]) %>% 
  ungroup %>% 
  group_by(i) %>% 
  nest %>% 
  pull(data) %>% 
  map(., ~fitdist(.x$.value, 'norm'))
