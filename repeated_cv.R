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
       mutate(filename = str_c(id, id2, sep = "_"))


# If R crashed, just start from the last point
folds %>% 
  filter(filename %notin% last_computed) %>% 
  mutate(bayes_rmse = map2(splits, filename, bayes_fit_and_predict))
