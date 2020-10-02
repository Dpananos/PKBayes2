library(tidyverse)
library(tidybayes)
library(cmdstanr)
library(posterior)
concentration_data <- read_csv("data/combined_data.csv") %>% 
  mutate(
    sex = if_else(sex == 'male', 0, 1),
    yobs = yobs/1000,
    study = if_else(study=='old',0,1)
  ) %>% 
  rename(from_new = study)


subject_data = concentration_data %>%
               distinct(subjectids, .keep_all = T) %>% 
               select(sex, weight, creatinine, age, D, from_new) %>% 
               mutate_at(vars(age, weight, creatinine), ~as.vector(scale(.))) %>% 
               compose_data()

obs_data = concentration_data %>% 
           select(yobs, time, subjectids) %>% 
           compose_data()

subject_data$n_subjectids = subject_data$n
subject_data = within(subject_data, rm(n))

model_data = c(subject_data, obs_data)

model_data$t_pred = seq(0.5,48, 2)
model_data$n_pred = length(model_data$t_pred)
model_data$dose_times = seq(0, 48, 12)
model_data$doses = rep(1, length(model_data$dose_times))
model_data$n_doses = length(model_data$dose_times)



model = cmdstan_model(stan_file = 'models/repeated_dosing.stan')

fit = model$sample(model_data, parallel_chains = 4)

