data{
  int N; //Total number of observations
  int n;
  int subjectids[N]; //Subject idendification number as an integer.  
  int n_subjects; //Number of unique subjects
  vector[N] times; //Times at which subjects were observed?  Length N
  real yobs[N]; //Observed concentraitons
  vector[N] D; //Size of the dose in mg
  
  int pcl;
  int pt;
  matrix[n,pcl] Xcl; //Covariates.  No intercept, mean centered.
  matrix[n,pt] Xt;
}
parameters{
  //Mu is on the log scale
  //Because x ~ lognormal -> log(x) ~ normal
  //I model on the log scale and then take exp to give me a lognormal rv
  real<lower=0>  mu_cl;
  real<lower=0> s_cl;
  vector[pcl] beta_cl;
  vector[n_subjects] z_cl;
  
  real<lower=0> mu_tmax;
  real<lower=0> s_t;
  vector[pt] beta_tmax;
  vector[n_subjects] z_t;
  
  real<lower=0, upper=1> mua;
  real<lower=0> mub;
  vector<lower=0, upper=1>[n_subjects] alpha ;
  
  real<lower=0, upper=1> phi;
  real<lower=0, upper=1> kappa;
  vector<lower=0, upper=1>[n_subjects] delays;
  
  real<lower=0> sigma;
}
transformed parameters{
  vector[n_subjects] Cl = exp(mu_cl + z_cl*s_cl + Xcl*beta_cl);
  vector[n_subjects] t = exp(mu_tmax + z_t*s_t + Xt*beta_tmax);
  vector[n_subjects] ka = log(alpha)./(t .* (alpha-1));
  vector[n_subjects] ke = alpha .* ka;
  //Since first non-zero observation is t = 0.5, we know delay is no larger than 0.5
  //since delay is beta, if we multiply by 0.5 then the largest the delay can be is 0.5
  vector[N] delayed_times = times - 0.5*delays[subjectids];
  
  //Vectorized concentration function
  vector[N] C = (0.5*D ./ Cl[subjectids]) .* (ke[subjectids] .* ka[subjectids]) ./ (ke[subjectids] - ka[subjectids]) .* (exp(-ka[subjectids] .* delayed_times) -exp(-ke[subjectids] .* delayed_times));
}
model{
  mu_tmax ~ normal(0.85, 0.06);
  s_t ~ gamma(77.7, 308.8);
  z_t ~ normal(0,1);
  
  mu_cl ~ normal(1.17, 0.05);
  s_cl ~ gamma(67.2, 312.22);
  z_cl ~ normal(0,1);

  mua ~ beta(47,200);
  mub ~ gamma(10,100);
  alpha ~ beta_proportion(mua,mub);
  phi ~ beta(76.3, 60.23);
  kappa ~ beta(24.24, 26.31);
  delays ~ beta(phi/kappa, (1-phi)/kappa);
  sigma ~ gamma(314.46, 1943.91);
  
  beta_tmax ~ normal(0, 0.07);
  
  beta_cl[1] ~ normal(0.31, 0.07);
  beta_cl[2] ~ normal(.02, 0.04);
  beta_cl[3] ~ normal(0.17, 0.04);
  beta_cl[4] ~ normal(0.02, 0.04);

  
  yobs ~ lognormal(log(C), sigma);
}
