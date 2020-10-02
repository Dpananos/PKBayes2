import cmdstanpy
import pandas as pd
import arviz as az
import numpy as np
from sklearn.preprocessing import StandardScaler

df = pd.read_csv('data/combined_data.csv')

def prepare_data(df):

    # Recode sex to be a binary voariable.  1==male    
    df['sex'] = df.sex.apply(lambda x: 1 if x=='male' else 0)

    # Include a fixed effect for study
    df['from_new'] = df.study.apply(lambda x: 1 if x=='new' else 0)

    # Recode subjectids to be unique numeric identifiers
    # Python is 0 indexed and Stan is 1 index, so add 1
    df['subjectids'] = pd.Categorical(df.subjectids).codes + 1

    conc_data = df.loc[:, ['subjectids','yobs','time']]
    subj_data = df.drop_duplicates('subjectids').loc[:, ['age','sex', 'weight','creatinine','from_new','D','subjectids']]
    subj_data.loc[:, ['age','creatinine','weight']] = StandardScaler().fit_transform(subj_data.loc[:, ['age','creatinine','weight']].values)


    model_data = {
        'n': conc_data.shape[0],
        'yobs': (conc_data.yobs/1000).tolist(),
        'time': conc_data.time.tolist(),
        'subjectids': (conc_data.subjectids).tolist(),
        'D': subj_data.D.tolist(),
        
        'n_subjectids': subj_data.shape[0],
        'age': subj_data.age.tolist(),
        'sex': subj_data.sex.tolist(),
        'weight': subj_data.weight.tolist(),
        'creatinine': subj_data.creatinine.tolist(),
        'from_new': subj_data.from_new.tolist(),
        
        't_pred' : np.arange(0.5, 75, 0.25).tolist(),
        'n_pred' : np.arange(0.5, 75, 0.25).size,
        'dose_times' : [0, 12, 24, 36, 48, 60, 72],
        'doses' : [1, 1, 1, 1, 2, 2, 2],
        'n_doses': 7
    }

    return model_data

model_data = prepare_data(df)

mod = cmdstanpy.CmdStanModel(stan_file='models/repeated_dosing.stan')

fit = mod.sample(model_data, show_progress = True, parallel_chains = 4)