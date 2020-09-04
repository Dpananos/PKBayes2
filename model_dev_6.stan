functions{
    real conc_curve(real D, real t, real Cl, real ka, real ke){
      real C = (0.5*D / Cl) * (ke * ka) / (ke - ka) * (exp(-ka * t) - exp(-ke * t));
      
      return C;
  }
}
data{
  int n; //Total number of observations
  int subjectids[n]; //Subject idendification number as an integer.  Mine go from 1 - 36
  int n_subjectids; //How many unique subjects do I have?
  vector[n] time; //time at which subjects were observed?  Length N
  real yobs[n]; //Observed concentraitons
  
  vector[n] sex;
  vector[n] weight;
  vector[n] creatinine;
  vector[n] age;
  vector[n] D;

  int n_test;
  vector[n_test] test_time;
  vector[n_test] test_sex;
  vector[n_test] test_weight;
  vector[n_test] test_creatinine;
  vector[n_test] test_age;
  vector[n_test] test_D;

}
parameters{
  
  real<lower=0>  mu_cl;
  real<lower=0> s_cl;                                                   
  vector[n_subjectids] z_cl;
  
  real<lower=0> mu_tmax;
  real<lower=0> s_t;
  vector[n_subjectids] z_t;

  real<lower=0, upper=1> phi;
  real<lower=0, upper=1> kappa;
  vector<lower=0, upper=1>[n_subjectids] delays;
  
  real<lower=0> sigma;
  
  real mu_alpha;
  real<lower=0> s_alpha;
  vector[n_subjectids] z_alpha;
  
  real beta_cl_sex;
  real beta_cl_weight;
  real beta_cl_creatinine;
  real beta_cl_age;
  
  real beta_t_weight;
  
}
transformed parameters{
  vector<lower=0>[n] Cl = exp(mu_cl + z_cl[subjectids]*s_cl + beta_cl_sex*sex + beta_cl_weight*weight + beta_cl_creatinine*creatinine + beta_cl_age*age);
  vector<lower=0>[n] t = exp(mu_tmax + z_t[subjectids]*s_t);
  vector<lower=0, upper=1>[n] alpha = inv_logit(mu_alpha + z_alpha[subjectids]*s_alpha);
  vector<lower=0>[n]ka = log(alpha)./(t .* (alpha-1));
  vector<lower=0>[n] ke = alpha .* log(alpha)./(t .* (alpha-1));
  vector<lower=0>[n] delayed_time = time - 0.5*delays[subjectids];
   
  vector<lower=0>[n] C = (0.5*D ./ Cl) .* (ke .* ka) ./ (ke - ka) .* (exp(-ka .* delayed_time) -exp(-ke .* delayed_time));
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
generated quantities{

  vector<lower=0>[n_test] cl_test = exp(mu_cl + beta_cl_sex*test_sex + beta_cl_weight*test_weight + beta_cl_creatinine*test_creatinine + beta_cl_age*test_age );
  real<lower=0> t_test = exp(mu_tmax);
  real<lower=0, upper=1> alpha_test = inv_logit(mu_alpha);
  real<lower=0>ka_test = log(alpha_test) /(t_test * (alpha_test-1));
  real<lower=0> ke_test = alpha_test * t_test;
  vector<lower=0>[n_test] delayed_time_test = test_time - 0.5*phi;

  vector[n_test] COOS;

  for (i in 1:n_test){
    COOS[i] = conc_curve(test_D[i], delayed_time_test[i], cl_test[i], ka_test, ke_test);
  }

}