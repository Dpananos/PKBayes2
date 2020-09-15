library(tidyverse)
library(tidybayes)
library(magrittr)

set.seed(19920908)

# How many patients from the new data do I want to hold out?
# There are ~200, so hold out half
TEST_SIZE = 100


# Load in the combined data.  
# Recode sex variable and rescale the concentrations from ng/ml to mg/L so that priors are on the same scale
concentration_data <- read_csv("data/combined_data.csv") %>% 
  mutate(
    sex = if_else(sex == 'male', 0, 1),
    yobs = yobs/1000 #Convert from ng/ml to mg/L
  )

# A dataframe for unique subjects. 
# Patients from the old data appear multiple times.
# Want the scaling not to be affected by multiple observations, so drop all duplicates and save to new df.
unique_subjects <- distinct(concentration_data, subjectids, .keep_all = T)


# Compute means and sds to pass to base::scale
covar_means <- unique_subjects %>% 
                select(age, creatinine, weight) %>% 
                summarise_all(mean) %>% 
                as.numeric


covar_sd <- unique_subjects %>% 
            select(age, creatinine, weight) %>% 
            summarise_all(sd) %>% 
            as.numeric

# Scale the covariates.  There may be a more sane way to do this, e.g. scaling  based on old data.
concentration_data[, c('age','creatinine','weight')] <- scale(concentration_data[, c('age','creatinine','weight')], center = covar_means, scale = covar_sd)



# Reserve some subjects from the new data for test.
# Sample only the new data for hold out.
test_subjects<-unique_subjects %>% 
                filter(study=='new') %>% 
                sample_n(TEST_SIZE) %>% 
                mutate(i = seq_along(subjectids))

write_csv(test_subjects, 'fit_data/test_subjects.csv')


train_subjects<- anti_join(unique_subjects, test_subjects, on = 'subjectids') %>% 
                 mutate(i = seq_along(subjectids))

write_csv(train_subjects, 'fit_data/train_subjects.csv')

# Scale the covariates in the hold out set so we can make predictions while we fit.
scaled_test_subjects<-test_subjects
scaled_test_subjects[, c('age','creatinine','weight')] <- scale(scaled_test_subjects[, c('age','creatinine','weight')], center = covar_means, scale = covar_sd)

# Rename all columns test_$x.  This is what they are named in the model
colnames(scaled_test_subjects) <- str_c('test_', colnames(scaled_test_subjects))

pred_data<-compose_data(scaled_test_subjects)
names(pred_data)[12]<-'test_n'

train_subjects<-anti_join(concentration_data, test_subjects, by = 'subjectids') %>% 
  mutate(subjectids = factor(subjectids))


# Concatenate train and validation into a single set to be passed to stan model.
model_data<-compose_data(train_subjects)
combined_data = c(model_data, pred_data)
saveRDS(combined_data, 'fit_data/combined_Data.RDS')