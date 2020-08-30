
data{
  int N; //Total number of observations
  int subjectids[N]; //Subject idendification number as an integer.  Mine go from 1 - 36
  int n_subjects; //How many unique subjects do I have?
  vector[N] time; //time at which subjects were observed?  Length N
  real yobs[N]; //Observed concentraitons
  
  vector[n_subjects] sex;
  vector[n_subjects] weight;
  vector[n_subjects] creatinine;
  vector[n_subjects] D;
}
parameters{
  
  real<lower=0>  mu_cl;
  real<lower=0> s_cl;
  vector[n_subjects] z_cl;
  
  real<lower=0> mu_tmax;
  real<lower=0> s_t;
  vector[n_subjects] z_t;
  
  
  
  real<lower=0, upper=1> phi;
  real<lower=0, upper=1> kappa;
  vector<lower=0, upper=1>[n_subjects] delays;
  
  real<lower=0> sigma;
  
  real mu_alpha;
  real<lower=0> s_alpha;
  vector[n_subjects] z_alpha;
  
  real beta_cl_sex;
  real beta_cl_weight;
  real beta_cl_creatinine;
  
  real beta_t_weight;
  
}
transformed parameters{
  vector<lower=0>[n_subjects] Cl = exp(mu_cl + z_cl*s_cl + beta_cl_sex*sex + beta_cl_weight*weight + beta_cl_creatinine*creatinine);
  vector<lower=0>[n_subjects] t = exp(mu_tmax + z_t*s_t + beta_t_weight*weight);
  vector<lower=0, upper=1>[n_subjects] alpha = inv_logit(mu_alpha + s_alpha*z_alpha);
  vector<lower=0>[n_subjects] ka = log(alpha)./(t .* (alpha-1));
  vector<lower=0>[n_subjects] ke = alpha .* log(alpha)./(t .* (alpha-1));
  vector<lower=0>[N] delayed_time = time - 0.5*delays[subjectids];
  
  vector<lower=0>[N] C = (0.5*D[subjectids] ./ Cl[subjectids]) .* (ke[subjectids] .* ka[subjectids]) ./ (ke[subjectids] - ka[subjectids]) .* (exp(-ka[subjectids] .* delayed_time) -exp(-ke[subjectids] .* delayed_time));
}
model{
  mu_tmax ~ normal(log(3.3), 0.25);
  s_t ~ gamma(10, 100);
  z_t ~ normal(0,1);
  
  mu_cl ~ normal(log(3.3),0.15);
  s_cl ~ gamma(15,100);
  z_cl ~ normal(0,1);
  
  
  mu_alpha ~ normal(0,1);
  s_alpha ~ gamma(10, 100);
  z_alpha ~ normal(0,1);
  
  phi ~ beta(20,20);
  kappa ~ beta(20,20);
  delays ~ beta(phi/kappa, (1-phi)/kappa);
  
  beta_cl_sex ~ student_t(3,0,2.5);
  beta_cl_weight ~ student_t(3,0,2.5);
  beta_cl_creatinine ~ student_t(3,0,2.5);
  beta_t_weight ~ student_t(3,0,2.5);
  
  sigma ~ lognormal(log(0.1), 0.2);
  yobs ~ lognormal(log(C), sigma);
}
