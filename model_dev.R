library(cmdstanr)
library(tidyverse)
library(tidybayes)  
library(posterior)
library(patchwork)
library(bayesplot)

num_cores <- parallel::detectCores()
# ---- Helper function ----

plot_random_effects<-function(draws, covar, random_effect){
  
  x<- enquo(covar)
  y<- enquo(random_effect)
  
  draws %>% 
    spread_draws((!!y)[i]) %>% 
    mean_qi %>% 
    left_join(subjects) %>% 
    ggplot(aes(!!x, !!y, color = sex, ymin = .lower, ymax = .upper))+
    geom_pointrange()+
    geom_smooth(method = 'lm')
}

plot_rfx_grid<-function(draws, random_effect, title){
  
  y <- enquo(random_effect)
  
  a <- plot_random_effects(draws, weight, !!y)
  b <- plot_random_effects(draws, age , !!y)
  c <- plot_random_effects(draws, creatinine , !!y)
  
  (a + b)/(c + plot_spacer()) + plot_annotation(title = title)
  
}

# ---- Rerun old model ----

combined_data<-read_csv('data/combined_data.csv') %>% 
               filter(study=='old') %>% 
               mutate(yobs=yobs/1000, i = seq_along(subjectids))

subjects<- distinct(combined_data, subjectids, .keep_all = T) %>% 
            mutate(i=subjectids)

model_1_data <- combined_data %>% 
              select(subjectids, time, yobs, D) %>% 
              compose_data(.n_name = n_prefix("N"))

model_1_data$n_subjects <- n_distinct(combined_data$subjectids)

model <- cmdstan_model(stan_file = 'model_dev_1.stan')

fit <- model$sample(data=model_1_data, 
                    parallel_chains=num_cores, 
                    chains=num_cores)

draws <- as_draws_df(fit$draws())

plot_rfx_grid(draws, z_cl, 'Clearance')
  
# ---- Rerun old model, but include some covariate effects on 

model_2_data <- combined_data %>% 
  select(subjectids, time, yobs, D) %>% 
  compose_data(.n_name = n_prefix("N"))

model_2_data$n_subjects <- n_distinct(combined_data$subjectids)
model_2_data$sex<- if_else(subjects$sex=='male',1,0)
model_2_data$weight<- as.numeric(scale(subjects$weight))
model_2_data$creatinine<- as.numeric(scale(subjects$creatinine))
model_2_data$D<- subjects$D

model <- cmdstan_model(stan_file = 'model_dev_2.stan')

fit <- model$sample(data=model_2_data, parallel_chains=num_cores, chains=num_cores)

draws <- as_draws_df(fit$draws())

plot_rfx_grid(draws, z_cl, 'Clearance')
plot_rfx_grid(draws, z_t, 'tmax')
plot_rfx_grid(draws, z_alpha, 'alpha')

mcmc_hist(fit$draws(), regex_pars = 'beta')

# ---- Summarize priors ----

draws %>% 
  gather_draws(sigma) %>% 
  pull(.value) %>% 
  fitdistrplus::fitdist(., 'gamma') %>% plot


# ---- New Data ----

d<-read_csv('data/combined_data.csv')

old = filter(d, study=='old')

old_subjects<- distinct(old, subjectids, .keep_all = T)


new = filter(d, study=='new') %>% 
  mutate(subjectids =  factor(subjectids - min(subjectids)+1),
         weight = (weight - mean(old_subjects$weight))/sd(old_subjects$weight),
         age = (age - mean(old_subjects$age))/sd(old_subjects$age),
         creatinine  = (creatinine  - mean(old_subjects$creatinine ))/sd(old_subjects$creatinine ),
         sex = if_else(sex=='male', 1, 0),
         yobs = yobs*0.001
  ) %>% 
  select(-study) 

model_data = compose_data(new, .n_name=n_prefix('N'))

model_data
model = cmdstan_model(stan_file='model_dev_3.stan')


fit = model$sample(model_data, chains = 4, parallel_chains = num_cores)

fit$draws('C') %>% 
  as_draws_df() %>% 
  spread_draws(C[i]) %>% 
  mean_qi %>% 
  mutate(i = factor(i)) %>% 
  left_join(new, by = c('i'='subjectids')) %>% 
  ggplot(aes(yobs, C))+
  geom_pointrange(aes(ymin = .lower, ymax = .upper))+
  geom_abline()+
  labs(x='observed', y='predicted')

# ---- New Model, only rfx on clearance ----

model = cmdstan_model(stan_file='model_dev_4.stan')
fit = model$sample(model_data, chains = 4, parallel_chains = num_cores)

fit$draws('C') %>% 
  as_draws_df() %>% 
  spread_draws(C[i]) %>% 
  mean_qi %>% 
  mutate(i = factor(i)) %>% 
  left_join(new, by = c('i'='subjectids')) %>% 
  ggplot(aes(yobs, C))+
  geom_pointrange(aes(ymin = .lower, ymax = .upper))+
  geom_abline()+
  labs(x='observed', y='predicted')

draws = fit$draws() %>% as_draws_df()


plot_rfx_grid(draws, t, '')
