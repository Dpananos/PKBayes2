functions {
  real[] one_comp_lin_elim_abs(real t,
                               real[] y,
                               real[] theta,
                               real[] x_r,
                               int[] x_i) {
    real dydt[1];
    real ka = theta[1]; // Dosing rate in 1/day
    real ke = theta[2];   // Elimination rate in 1/day
    real D = x_r[1];
    real Cl = x_r[2];
    real elim = ke * y[1];
    
    dydt[1] = 0.5*D*ke*ka*exp(-ka*t)/Cl  - elim;
    
    for( i in 1:8){
      
      dydt[1] = dydt[1] + (0.5*D*ke*ka*exp(-ka*(t-12*i))/Cl)*step(t-12*i);
      
    }

    

    return dydt;
  }
}
data{
  int N; //Total number of observations
  int subjectids[N]; //Subject idendification number as an integer.  Mine go from 1 - 36
  int N_subjectids; //How many unique subjects do I have?
  vector[N] time; //time at which subjects were observed?  Length N
  real yobs[N]; //Observed concentraitons
  
  vector[N_subjectids] sex;
  vector[N_subjectids] weight;
  vector[N_subjectids] age;
  vector[N_subjectids] creatinine;
  vector[N_subjectids] D;
}
parameters{
  
  real<lower=0>  mu_cl;
  real<lower=0> s_cl;
  vector[N_subjectids] z_cl;
  
  real<lower=0> mu_tmax;
  real<lower=0> s_t;
  vector[N_subjectids] z_t;
  
  
  
  real<lower=0, upper=1> phi;
  real<lower=0, upper=1> kappa;
  vector<lower=0, upper=1>[N_subjectids] delays;
  
  real<lower=0> sigma;
  
  real mu_alpha;
  real<lower=0> s_alpha;
  vector[N_subjectids] z_alpha;
  
  
  real beta_cl_sex;
  real beta_cl_weight;
  real beta_cl_creatinine;
  
  real beta_t_weight;
  
  
}
transformed parameters{
  vector<lower=0>[N_subjectids] Cl = exp(mu_cl + z_cl*s_cl + beta_cl_sex*sex + beta_cl_weight*weight + beta_cl_creatinine*creatinine);
  vector<lower=0>[N_subjectids] t = exp(mu_tmax + z_t*s_t + beta_t_weight*weight);
  vector<lower=0, upper=1>[N_subjectids] alpha = inv_logit(mu_alpha + s_alpha*z_alpha);
  vector<lower=0>[N_subjectids] ka = log(alpha)./(t .* (alpha-1));
  vector<lower=0>[N_subjectids] ke = alpha .* log(alpha)./(t .* (alpha-1));
  vector<lower=0>[N] delayed_time = time - 0.5*delays[subjectids];
  
  vector<lower=0>[N] C = (0.5*D[subjectids] ./ Cl[subjectids]) .* (ke[subjectids] .* ka[subjectids]) ./ (ke[subjectids] - ka[subjectids]) .* (exp(-ka[subjectids] .* delayed_time) -exp(-ke[subjectids] .* delayed_time));
}
model{
  //updated model priors
  mu_tmax ~ normal(0.87, 0.04);
  s_t ~ gamma(81, 363);
  z_t ~ normal(0,1);
  
  mu_cl ~ normal(0.47,0.04);
  s_cl ~ gamma(70.7, 346.66);
  z_cl ~ normal(0,1);
  
  
  mu_alpha ~ normal(-1.47,0.12);
  s_alpha ~ gamma(9.14, 89.27);
  z_alpha ~ normal(0,1);
  
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
  for (i in 1:N){
    posterior_predict[i] = lognormal_rng(log(C[i]),sigma);
  }
}
