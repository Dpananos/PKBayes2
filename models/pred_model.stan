data{
  int n; //Total number of observations
  vector[n] time; //time at which subjects were observed?  Length N
  vector[n] sex;
  vector[n] weight;
  vector[n] creatinine;
  vector[n] age;
  vector[n] D;
  
  
}
parameters{
  
  real mu_cl;
  real mu_tmax;
  real phi;
  real mu_alpha;
  
  real beta_cl_sex;
  real beta_cl_weight;
  real beta_cl_creatinine;
  real beta_cl_age;
  
}

generated quantities{
    vector<lower=0>[n] Clp = exp(mu_cl + beta_cl_sex*sex + beta_cl_weight*weight + beta_cl_creatinine*creatinine + beta_cl_age*age);
    real<lower=0> tp = exp(mu_tmax);
    real<lower=0> alphap = inv_logit(mu_alpha);
    
    real<lower=0> kap = log(alphap)/(tp * (alphap-1));
    real<lower=0> kep = alphap * kap;
    vector<lower=0>[n] delayed_timep = time - 0.5*phi;
    
    vector<lower=0>[n] Cp = (0.5*D ./ Clp) * (kep * kap) / (kep - kap) .* (exp(-kap * delayed_timep) -exp(-kep * delayed_timep));
}