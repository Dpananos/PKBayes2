library(tidyverse)
library(here)


old<-read_csv('data/rommels_apixiban_data.csv')  %>% 
  transmute(
    time = Time,
    subject = as.factor(Subject),
    yobs = Concentration,
    sex = str_to_lower(Sex), 
    age = Age, 
    weight = Weight,
    creatinine = Creatinine,
    D = 2.5,
    study = 'old'
  )


new<-read_csv('data/apixaban_2020_data.csv') %>% 
  transmute(
    subject = `New Subject ID`,
    age = enrollment_age,
    sex = str_to_lower(sex),
    weight = weight_kg,
    time = hrs_post_dose,
    yobs = concentration_ng_per_ml,
    creatinine = creatinine_umol_per_L,
    D = dose_mg,
    study='new')


df <-  bind_rows(old, new)
  


model_data <- df %>% 
  distinct(subject) %>% 
  mutate(subjectids=seq_along(subject)) %>% 
  right_join(df) %>% 
  select(-subject)


write_csv(model_data, 'data/combined_data.csv')
