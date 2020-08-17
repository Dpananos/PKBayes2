data{
  int N; //Total number of observations
  int subjectids[N]; //Subject idendification number as an integer.  Mine go from 1 - 36
  int n_subjects; //How many unique subjects do I have?
  vector[N] time; //time at which subjects were observed?  Length N
  real yobs[N]; //Observed concentraitons
}
parameters{
  
  real<lower=0>  mu_cl;
  real<lower=0> s_cl;
  vector[n_subjects] z_cl;
  
  real<lower=0> mu_tmax;
  real<lower=0> s_t;
  vector[n_subjects] z_t;
  
  vector<lower=0, upper=1>[n_subjects] alpha ;
  
  real<lower=0, upper=1> phi;
  real<lower=0, upper=1> kappa;
  vector<lower=0, upper=1>[n_subjects] delays;
  
  real<lower=0> sigma;
  
  real<lower=0, upper=1> alpha_a;
  real<lower=0, upper=1> alpha_b;
  
}
transformed parameters{
  vector<lower=0>[n_subjects] Cl = exp(mu_cl + z_cl*s_cl);
  vector<lower=0>[n_subjects] t = exp(mu_tmax + z_t*s_t);
  vector<lower=0>[n_subjects] ka = log(alpha)./(t .* (alpha-1));
  vector<lower=0>[n_subjects] ke = alpha .* ka;
  vector<lower=0>[N] delayed_time = time - 0.5*delays[subjectids];
  
  vector<lower=0>[N] C = (2.5 ./ Cl[subjectids]) .* (ke[subjectids] .* ka[subjectids]) ./ (ke[subjectids] - ka[subjectids]) .* (exp(-ka[subjectids] .* delayed_time) -exp(-ke[subjectids] .* delayed_time));
}
model{
  mu_tmax ~ normal(log(3.3), 0.25);
  s_t ~ gamma(10, 100);
  z_t ~ normal(0,1);
  
  mu_cl ~ normal(log(3.3),0.15);
  s_cl ~ gamma(15,100);
  z_cl ~ normal(0,1);

  alpha_a ~ beta(1,1);
  alpha_b ~ beta(2,2);
  alpha ~ beta(alpha_a/alpha_b, (1-alpha_a)/alpha_b);
  
  phi ~ beta(20,20);
  kappa ~ beta(20,20);
  delays ~ beta(phi/kappa, (1-phi)/kappa);
  
  sigma ~ lognormal(log(0.1), 0.2);
  yobs ~ lognormal(log(C), sigma);
}
