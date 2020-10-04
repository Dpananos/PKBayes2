functions{
  vector pk_ode(real t, vector y, real Cl, real ke, real ka,  int n_doses, vector doses, vector dose_times){
    
    vector[1] dydt;
    dydt[1] = -ke*y[1];
    
    for(i in 1:n_doses){
      if(t>dose_times[i]){
        dydt[1] = dydt[1] + 0.5*doses[i] * ke * ka * exp( -ka * (t-dose_times[i]) ) / Cl;
      }
    }
    
    return dydt;
  }
  
  real concentration(real t, real Cl, real ke, real ka){
    real y = (0.5 / Cl) * (ke * ka) / (ke - ka) * (exp(-ka * t) - exp(-ke * t));
    return y;
  }
}

data{
  int n; //Total number of observations
  int subjectids[n]; //Subject idendification number as an integer.  Mine go from 1 - 36
  int n_subjectids; //How many unique subjects do I have?
  vector[n] time; //time at which subjects were observed?  Length N
  real yobs[n]; //Observed concentraitons
  
  //Covars
  vector[n_subjectids] sex;
  vector[n_subjectids] weight;
  vector[n_subjectids] creatinine;
  vector[n_subjectids] age;
  vector[n_subjectids] D;
  vector[n_subjectids] from_new;
  
  //for gen quants
  int n_doses;
  int n_pred;
  real t_pred[n_pred];
  vector[n_doses] doses;
  vector[n_doses] dose_times;
}
transformed data{
  vector[n_subjectids] scaled_weight = (weight - mean(weight))/sd(weight);
  vector[n_subjectids] scaled_age = (age - mean(age))/sd(age);
  vector[n_subjectids] scaled_creatinine = (creatinine - mean(creatinine))/sd(creatinine);
  matrix[n_subjectids, 5] X = [sex', scaled_weight', scaled_creatinine', scaled_age', from_new']';
  real t0 = 0;
  vector[1] y0 = [0]';
}
parameters{
  
  real  mu_cl;
  real<lower=0> s_cl;                                                   
  vector[n_subjectids] z_cl;
  
  real mu_tmax;
  real<lower=0> s_t;
  vector[n_subjectids] z_t;
  
  real<lower=0, upper=1> phi;
  real<lower=0, upper=1> kappa;
  vector<lower=0, upper=1>[n_subjectids] delays;
  
  real<lower=0> sigma;
  
  real mu_alpha;
  real<lower=0> s_alpha;
  vector[n_subjectids] z_alpha;
  
  vector[5] beta_cl;
  vector[5] beta_t;
  vector[5] beta_a;
}
transformed parameters{
  vector<lower=0>[n_subjectids] Cl = exp(mu_cl + z_cl*s_cl + X*beta_cl);
  vector<lower=0>[n_subjectids] tmax = exp(mu_tmax + z_t*s_t + X*beta_t);
  vector<lower=0, upper=1>[n_subjectids] alpha = inv_logit(mu_alpha + z_alpha*s_alpha + X*beta_a);
  vector<lower=0>[n_subjectids]ka = log(alpha)./(tmax .* (alpha-1));
  vector<lower=0>[n_subjectids] ke = alpha .* log(alpha)./(tmax .* (alpha-1));
  vector<lower=0>[n] delayed_time = time - 0.5*delays[subjectids];
  
  vector<lower=0>[n] C = (0.5*D[subjectids] ./ Cl[subjectids]) .* (ke[subjectids] .* ka[subjectids]) ./ (ke[subjectids] - ka[subjectids]) .* (exp(-ka[subjectids] .* delayed_time) -exp(-ke[subjectids] .* delayed_time));
}
model{
  //See Byon et. al 2019
  mu_tmax ~ normal(log(3.3), 0.1);
  s_t ~ gamma(5, 100);
  z_t ~ normal(0,1);
  
  mu_cl ~ normal(log(3.3),0.15);
  s_cl ~ gamma(15,100);
  z_cl ~ normal(0,1);
  
  
  mu_alpha ~ normal(-0.25,0.5);
  s_alpha ~ gamma(10, 100);
  z_alpha ~ normal(0,1);
  
  
  phi ~ beta(20,20);
  kappa ~ beta(20,20);
  delays ~ beta(phi/kappa, (1-phi)/kappa);
  
  beta_cl ~ normal(0,0.25);
  beta_t ~ normal(0, 0.25);
  beta_a ~ normal(0, 0.25);
  
  
  sigma ~ lognormal(log(0.1), 0.2);
  yobs ~ lognormal(log(C), sigma);
}
// generated quantities{
//   
//   matrix[n_subjectids, n_pred] C_repeat;
//   matrix[n_subjectids, n_pred] C_repeat_ppc;
//   
//   for(i in 1:n_subjectids){
//     C_repeat[i,] = to_row_vector(ode_rk45(pk_ode, y0, t0, t_pred, Cl[i], ke[i], ka[i],  n_doses, doses,  dose_times )[,1]);
//     C_repeat_ppc[i,] = to_row_vector(lognormal_rng( log(C_repeat[i,]), sigma ));  
//   }
//   
//   
//   
//   
// }
// 
