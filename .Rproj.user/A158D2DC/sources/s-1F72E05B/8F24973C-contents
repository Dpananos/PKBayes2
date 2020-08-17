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
  int n_subjects; //How many unique subjects do I have?
  vector[N] time; //time at which subjects were observed?  Length N
  real yobs[N]; //Observed concentraitons
  
  vector[n_subjects] sex;
  vector[n_subjects] weight;
  vector[n_subjects] D;
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
  
  real beta_sex;
  real beta_weight;
  
}
transformed parameters{
  vector<lower=0>[n_subjects] Cl = exp(mu_cl + z_cl*s_cl + beta_sex*sex + beta_weight*weight);
  vector<lower=0>[n_subjects] t = exp(mu_tmax + z_t*s_t);
  vector<lower=0>[n_subjects] ka = log(alpha)./(t .* (alpha-1));
  vector<lower=0>[n_subjects] ke = alpha .* ka;
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
  
  alpha_a ~ beta(1,1);
  alpha_b ~ beta(5,5);
  alpha ~ beta(alpha_a/alpha_b, (1-alpha_a)/alpha_b);
  
  phi ~ beta(20,20);
  kappa ~ beta(20,20);
  delays ~ beta(phi/kappa, (1-phi)/kappa);
  
  beta_sex ~ student_t(3,0,2.5);
  beta_weight ~ student_t(3,0,2.5);
  
  sigma ~ lognormal(log(0.1), 0.2);
  yobs ~ lognormal(log(C), sigma);
}
generated quantities{
  real popcl = exp(mu_cl);
  real popa = alpha_a;
  real poptmax = exp(mu_tmax);
  real popka = log(popa)/(poptmax*(popa-1));
  real popke = popa*popka;
  real Crepeat[8,1];
  
  {
    real theta[2] = {popka, popke};
    real x_r[2] = {1, popcl};
    int x_i[0];
    Crepeat = integrate_ode_bdf(one_comp_lin_elim_abs, {0.0}, 0.0, {12.0,24.0,36.0,48.0,60.0,72.0,84.0,96.0}, theta, x_r, x_i);
  }
  
  
  
}
