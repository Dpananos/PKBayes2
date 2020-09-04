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

}
parameters{
  
  real<lower=0>  mu_cl;
  
  real<lower=0> mu_tmax;

  real mu_ke;

  
  real<lower=0> sigma;
  
  
  real beta_sex;
  real beta_weight;
  real beta_creatinine;
  real beta_age;
  
}
transformed parameters{
  real<lower=0> Cl = exp(mu_cl);
  real<lower=0> t = exp(mu_tmax);
  real ka = 0.42;
  vector<upper=ka>[n] ke = ka*inv_logit(mu_ke + beta_sex*sex + beta_weight*weight + beta_creatinine*creatinine + beta_age*age);
   
  vector<lower=0>[n] C = (0.5 * D / Cl) .* (ke * ka) ./ (ke - ka) .* (exp(-ka * time) - exp(-ke .* time));
}
model{
  mu_tmax ~ normal(0,2);

  
  mu_cl ~ normal(0,2);

  mu_ke ~ std_normal();

  
  beta_sex ~ student_t(3,0,2.5);
  beta_weight ~ student_t(3,0,2.5);
  beta_creatinine ~ student_t(3,0,2.5);
  beta_age ~ student_t(3, 0, 2.5);
  
  sigma ~ cauchy(0,1);
  yobs ~ lognormal(log(C), sigma);
}
