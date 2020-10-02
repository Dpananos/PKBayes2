library(tidyverse)
library(tidymodels)
source('analysis/02_prepare_data.R')
`%notin%` <- Negate(`%in%`)

#Get the folds I've already done up to this point
flnames<-list.files('data/bayes_error/')

#Remove the extension
last_computed<-sub('\\.txt', '', flnames )

# Create folds. Seed ensures reproducibility
set.seed(0)
folds<-vfold_cv(new_subjects, v = 10, repeats = 100 ) %>% 
       mutate(filename = str_c(id, id2, sep = "_")) %>% 
       rownames_to_column('.id')


# If R crashed, just start from the last point
folds %>% 
  filter(filename %notin% last_computed) %>% 
  mutate(bayes_rmse = map2(splits, filename, bayes_fit_and_predict))


results<-list.files('data/bayes_error/', full.names = T)

bayes_error = map_dfr(results , ~{read_file(.x) %>%
                    str_remove(., '\n') %>% 
                    as.numeric() %>% 
                    tibble(error = .)}, .id = '.id') 



folds %>% 
  inner_join(bayes_error) %>% 
  mutate(linmod_error = map_dbl(splits, linmod_fit_and_predict)) %>% 
  group_by(id) %>% 
  summarise(bayes = mean(error), linmod = mean(linmod_error)) %>% 
  mutate(delta = linmod - bayes) %>% 
  summarise(mean(1000*delta))
  
