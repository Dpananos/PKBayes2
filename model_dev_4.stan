functions {
  real conc_curve(real D, real t, real Cl, real ka, real ke){
    real C = (0.5*D / Cl) * (ke * ka) / (ke - ka) * (exp(-ka * t) - exp(-ke * t));
    
    return C;
  }
}
data{
  int N; //Total number of observations
  int N_test;
  int subjectids[N]; //Subject idendification number as an integer.  Mine go from 1 - 36
  
  int N_subjectids; //How many unique subjects do I have?
  
  int N_test_subjectids_test;
  int subjectids_test[N_test];
  
  vector[N] time; //time at which subjects were observed?  Length N
  real yobs[N]; //Observed concentraitons
  
  vector[N_subjectids] sex;
  vector[N_subjectids] weight;
  vector[N_subjectids] age;
  vector[N_subjectids] creatinine;
  vector[N_subjectids] D;
  
  vector[N_test] test_time;
  vector[N_test_subjectids_test] sex_test;
  vector[N_test_subjectids_test] weight_test;
  vector[N_test_subjectids_test] age_test;
  vector[N_test_subjectids_test] creatinine_test;
  vector[N_test_subjectids_test] D_test;
}
parameters{
  
  real<lower=0>  mu_cl;
  real<lower=0> s_cl;
  vector[N_subjectids] z_cl;
  
  real<lower=0> mu_tmax;
  
  real<lower=0, upper=1> phi;
  real<lower=0, upper=1> kappa;
  real<lower=0, upper=1> delays;
  
  real<lower=0> sigma;
  
  real mu_alpha;
  
  
  real beta_cl_sex;
  real beta_cl_weight;
  real beta_cl_creatinine;
  
  real beta_t_weight;
  
  
}
transformed parameters{
  vector<lower=0>[N_subjectids] Cl = exp(mu_cl + z_cl*s_cl + beta_cl_sex*sex + beta_cl_weight*weight + beta_cl_creatinine*creatinine);
  vector<lower=0>[N_subjectids] t = exp(mu_tmax + beta_t_weight*weight);
  real alpha = inv_logit(mu_alpha);
  vector<lower=0>[N_subjectids] ka = log(alpha) ./(t * (alpha-1));
  vector<lower=0>[N_subjectids] ke = alpha * log(alpha)./(t * (alpha-1));
  vector<lower=0>[N] delayed_time = time - 0.5*delays;
  
  vector<lower=0>[N] C = (0.5*D[subjectids] ./ Cl[subjectids]) .* (ke[subjectids] .* ka[subjectids]) ./ (ke[subjectids] - ka[subjectids]) .* (exp(-ka[subjectids] .* delayed_time) -exp(-ke[subjectids] .* delayed_time));
}
model{
  //updated model priors
  mu_tmax ~ normal(0.87, 0.04);
  
  mu_cl ~ normal(0.47,0.04);
  s_cl ~ gamma(70.7, 346.66);
  z_cl ~ normal(0,1);
  
  
  mu_alpha ~ normal(-1.47,0.12);
  
  phi ~ beta(72.4, 57.8);
  kappa ~ beta(25,26);
  delays ~ beta(phi/kappa, (1-phi)/kappa);
  
  beta_cl_sex ~ normal(0.44, 0.07);
  beta_cl_weight ~ normal(0.18, 0.04);
  beta_cl_creatinine ~ normal(0.02, 0.03);
  beta_t_weight ~ normal(0.03, 0.04);
  
  sigma ~ gamma(282.3, 2241.62);
  yobs ~ lognormal(log(C), sigma);
}
generated quantities{
  vector[N] posterior_predict;
  
  vector[N_test_subjectids_test] Cl_test = exp(mu_cl + beta_cl_sex*sex_test + beta_cl_weight*weight_test + beta_cl_creatinine*creatinine_test);
  vector<lower=0>[N_test_subjectids_test] t_test = exp(mu_tmax + beta_t_weight*weight_test);
  vector<lower=0>[N_test_subjectids_test] ka_test = log(alpha)./ (t_test * (alpha-1));
  vector<lower=0>[N_test_subjectids_test] ke_test = alpha * log(alpha)./ (t_test * (alpha-1));
  vector<lower=0>[N_test] delayed_test_time = test_time - 0.5*phi;
  
  vector[N_test] C_test = (0.5*D_test[subjectids_test] ./ Cl_test[subjectids_test]) .* (ke_test[subjectids_test] .* ka_test[subjectids_test]) ./ (ke_test[subjectids_test] - ka_test[subjectids_test]) .* (exp(-ka_test[subjectids_test] .* delayed_test_time) -exp(-ke_test[subjectids_test] .* delayed_test_time));
  

  
  for (i in 1:N){
    posterior_predict[i] = lognormal_rng(log(C[i]),sigma);
  }
  

}
