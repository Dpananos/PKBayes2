library(tidyverse)
library(here)

# This script combines rommel's old data with the new data from the clinic.
# I use the result of this script in future model fits (e.g validating OOS model prediction capabilities)
# I only use the sex, age, weight, creatinine variables from each data set.  I also include the dose (to be passed to models)
# as well as an indicator for which study the observation came from.  I use a fixed effect for study in the models.


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
