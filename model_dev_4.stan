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
  
  vector<lower=0, upper=1>[N_subjectids] alpha ;
  
  real<lower=0, upper=1> phi;
  real<lower=0, upper=1> kappa;
  vector<lower=0, upper=1>[N_subjectids] delays;
  
  real<lower=0> sigma;
  
  real<lower=0, upper=1> alpha_a;
  real<lower=0, upper=1> alpha_b;
  
  real beta_sex;
  real beta_weight;
  real beta_age;
  real beta_creatinine ;
  
  vector<lower=0>[N_subjectids] y0;
  
}
transformed parameters{
  vector<lower=0>[N_subjectids] Cl = exp(mu_cl + z_cl*s_cl + beta_sex*sex + beta_weight*weight + beta_age*age + beta_creatinine *creatinine );
  vector<lower=0>[N_subjectids] t = exp(mu_tmax + z_t*s_t);
  vector<lower=0>[N_subjectids] ka = log(alpha)./(t .* (alpha-1));
  vector<lower=0>[N_subjectids] ke = alpha .* ka;
  vector<lower=0>[N] delayed_time = time - 0.5*delays[subjectids];
  
  vector<lower=0>[N] C = (0.5*D[subjectids] ./ Cl[subjectids]) .* (ke[subjectids] .* ka[subjectids]) ./ (ke[subjectids] - ka[subjectids]) .* (exp(-ka[subjectids] .* delayed_time) -exp(-ke[subjectids] .* delayed_time)) + D[subjectids].*y0[subjectids];
}
model{
  mu_tmax ~ normal(0.87, 0.05);
  mu_cl ~ normal(0.48,0.05);
  beta_sex ~ normal(0.43, 0.08);
  beta_weight ~ normal(0.19, 0.04);
  beta_age ~ normal(0,1);
  beta_creatinine  ~ normal(0,1);
  
  z_t ~ normal(0,1);
  z_cl ~ normal(0,1);
  
  
  s_t ~ gamma(76.44, 300);
  s_cl ~ gamma(74, 365.23);
  sigma ~ lognormal(308.75, 1918);
  
  alpha_a ~ beta(38, 140);
  alpha_b ~ beta(3.5, 41.7);
  alpha ~ beta(alpha_a/alpha_b, (1-alpha_a)/alpha_b);
  
  phi ~ beta(73.57, 59.25);
  kappa ~ beta(24.94, 26.61);
  delays ~ beta(phi/kappa, (1-phi)/kappa);
  
  //This is an awful approximation.  Fix this.
  y0 ~ lognormal(-4.6,0.25);
  yobs ~ lognormal(log(C), sigma);
}
