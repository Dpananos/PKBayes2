functions{
  // This is a helper function to compute the concentration curve given PK parameters
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
  
  //Covars
  vector[n] sex;
  vector[n] weight;
  vector[n] creatinine;
  vector[n] age;
  vector[n] D;
  vector[n] from_new;
  
  
    //Covars
  int test_n;
  vector[test_n] test_time;
  vector[test_n] test_sex;
  vector[test_n] test_weight;
  vector[test_n] test_creatinine;
  vector[test_n] test_age;
  vector[test_n] test_D;
  vector[test_n] test_from_new;
  
  
}
transformed data{
  matrix[n, 5] X = [sex', weight', creatinine', age', from_new']';
  matrix[test_n, 5] test_X = [test_sex', test_weight', test_creatinine', test_age', test_from_new']';
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
  vector<lower=0>[n] Cl = exp(mu_cl + z_cl[subjectids]*s_cl + X*beta_cl);
  vector<lower=0>[n] t = exp(mu_tmax + z_t[subjectids]*s_t + X*beta_t);
  vector<lower=0, upper=1>[n] alpha = inv_logit(mu_alpha + z_alpha[subjectids]*s_alpha + X*beta_a);
  vector<lower=0>[n]ka = log(alpha)./(t .* (alpha-1));
  vector<lower=0>[n] ke = alpha .* log(alpha)./(t .* (alpha-1));
  vector<lower=0>[n] delayed_time = time - 0.5*delays[subjectids];
  
  vector<lower=0>[n] C = (0.5*D ./ Cl) .* (ke .* ka) ./ (ke - ka) .* (exp(-ka .* delayed_time) -exp(-ke .* delayed_time));
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
generated quantities{
 real z1;
 real z2;
 real z3;
 
 vector[test_n] CL ;
 vector[test_n] T ;
 vector[test_n] A;
 vector[test_n] KA ;
 vector[test_n] KE ;
 
 vector[test_n] TIME = test_time - 0.5*beta_rng(phi/kappa, (1-phi)/kappa);
 
 vector[test_n] CPRED;
 
 
 for (i in 1:test_n){
   z1 = normal_rng(0,1);
   z2 = normal_rng(0,1);
   z3 = normal_rng(0,1);
   
   CL[i] = exp(mu_cl + s_cl*z1 + test_X[i]*beta_cl);
   T[i] = exp(mu_tmax + s_t*z2 + test_X[i]*beta_t);
   A[i] = inv_logit(mu_alpha + s_alpha*z3 + test_X[i]*beta_a );
   KA[i] = log(A[i]) / (T[i] * (A[i]-1));
   KE[i] = A[i] * KA[i];
   CPRED[i] = conc_curve(test_D[i], TIME[i], CL[i], KA[i], KE[i]);
 }
 
}

