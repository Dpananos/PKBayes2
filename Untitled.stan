functions{
  vector heaviside(vector t){
    
    vector[size(t)] y;
    for(i in 1:size(t)){
        y[i] = t[i]<0 ? 0 : 1;
    }
    return y;
  }
  

  vector conc(real D, vector t, real Cl, real ka, real ke){
    
    return heaviside(t) .* (exp(-ka*t) - exp(-ke*t)) * (0.5 * D * ke * ka ) / (Cl *(ke - ka));
  }
  
}
data{
  int n;
  vector[n] C_obs;
  vector[n] times;
  row_vector[5] X;
  
  
  
  int n_doses;
  vector[n_doses] doses;
  vector[n_doses] dose_times;
  
  int nt;
  vector[nt] tpred;
}

parameters{
  real mu_cl;
  vector[5] beta_cl;
  real<lower=0> s_cl;
  real z_cl;
  
  real mu_t;
  vector[5] beta_t;
  real <lower=0> s_t;
  real z_t;
  
  
  real mu_alpha;
  real<lower=0> s_alpha;
  vector[5] beta_alpha;
  real z_a;
  
  real<lower=0> sigma;
}
transformed parameters{
    real<lower=0, upper=1> alpha = inv_logit(mu_alpha + X*beta_alpha + z_a*s_alpha);
    real<lower=0> tmax = exp(mu_t + X*beta_t + z_t*s_t);
    real<lower=0> cl = exp(mu_cl + X*beta_cl + z_cl*s_cl);
    real<lower=0> ka = log(alpha)/(tmax * (alpha-1));
    real<lower=0, upper = ka> ke = alpha*ka;
    vector<lower=0>[n] C = rep_vector(0.0, n);
    for(i in 1:n_doses){
      C += conc(doses[i], times - dose_times[i], cl, ka, ke);
    }
}
model{
    mu_cl ~ normal(0.51, 0.08);
    s_cl ~ gamma(240.72, 628.91);
    beta_cl[1] ~ normal(-0.23, 0.07);
    beta_cl[2] ~ normal(0.09, 0.03);
    beta_cl[3] ~ normal(-0.18, 0.07);
    beta_cl[4] ~ normal(-0.11, 0.04);
    beta_cl[5] ~ normal(-0.89, 0.1);
    z_cl ~ std_normal();
    
    
    beta_t[1] ~ normal(-0.06, 0.07);
    beta_t[2] ~ normal(0.09, 0.04);
    beta_t[3] ~ normal(0.11, 0.07);
    beta_t[4] ~ normal(0.04, 0.04);
    beta_t[5] ~ normal(0.07, 0.1);
    s_t ~ gamma(64.73, 301.12);
    z_t ~ std_normal();
    
    
    mu_alpha ~ normal(-1.18, 0.20);
    s_alpha  ~ gamma(0.52, 90.33);
    beta_alpha[1] ~ normal(-0.16, 0.003);
    beta_alpha[2] ~ normal(0.31, 0.17);
    beta_alpha[3] ~ normal(0.11, 0.07);
    beta_alpha[4] ~ normal(0, 0.09);
    beta_alpha[5] ~ normal(-0.24, 0.24);
    z_a ~ std_normal();
    
    sigma ~  gamma(333.88, 1920.92);
    C_obs ~ lognormal(log(C), sigma);
    
}
generated quantities{
  vector[nt] ypred = rep_vector(0.0, nt);
  
   for(i in 1:n_doses){
    ypred += conc(doses[i], tpred - dose_times[i], cl, ka, ke);
  }
}