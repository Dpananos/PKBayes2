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
}

data{
  int n_subjects;
  matrix[n_subjects, 5] X;
  
  int n_doses;
  int n_pred;
  real t_pred[n_pred];
  vector[n_doses] doses;
  vector[n_doses] dose_times;

  
}
transformed data{
  real t0 = 0;
  vector[1] y0 = [0]';
}
generated quantities{
  
  real mu_cl = normal_rng(0.51, 0.08);
  real mu_tmax = normal_rng(1.07, 0.07);
  real mu_alpha = normal_rng(-1.18, 0.20);
  
  real s_cl = gamma_rng(240.72, 628.91);
  real s_tmax = gamma_rng(64.73, 301.12);
  real s_alpha = gamma_rng(0.52, 90.33);
  
  real b1_cl = normal_rng(-0.23, 0.07);
  real b2_cl = normal_rng(0.09, 0.03);
  real b3_cl = normal_rng(-0.18, 0.07);
  real b4_cl = normal_rng(-0.11, 0.04);
  real b5_cl = normal_rng(-0.89, 0.1);
  vector[5] beta_cl = to_vector([b1_cl, b2_cl, b3_cl, b4_cl, b5_cl]);
  
  real b1_t = normal_rng(-0.06, 0.07);
  real b2_t = normal_rng(0.09, 0.04);
  real b3_t = normal_rng(0.11, 0.07);
  real b4_t = normal_rng(0.04, 0.04);
  real b5_t = normal_rng(0.07, 0.1);
  vector[5] beta_t = to_vector([b1_t, b2_t, b3_t, b4_t, b5_t]);
    
  real b1_alpha = normal_rng(-0.16, 0.003);
  real b2_alpha = normal_rng(0.31, 0.17);
  real b3_alpha = normal_rng(0.11, 0.07);
  real b4_alpha = normal_rng(0, 0.09);
  real b5_alpha = normal_rng(-0.24, 0.24);
  vector[5] beta_alpha = to_vector([b1_alpha, b2_alpha, b3_alpha, b4_alpha, b5_alpha]);

  vector[n_subjects] cl = to_vector(exp(normal_rng(mu_cl + X*beta_cl, s_cl)));
  vector[n_subjects] tmax = to_vector(exp(normal_rng(mu_tmax + X*beta_t, s_tmax)));
  vector[n_subjects] alpha =  to_vector(inv_logit(normal_rng(mu_alpha + X*beta_alpha, s_alpha)));
  vector[n_subjects] ka = log(alpha)./(tmax .* (alpha-1));
  vector[n_subjects] ke = alpha .* ka;
  
  real sigma = gamma_rng(333.88, 1920.92);
  
  matrix[n_subjects, n_pred] C;
  matrix[n_subjects, n_pred] C_obs;
  
  for(i in 1:n_subjects){
      C[i,] = to_row_vector(ode_rk45(pk_ode, y0, t0, t_pred, cl[i], ke[i], ka[i],  n_doses, doses,  dose_times )[,1]);
      C_obs[i,] = to_row_vector(lognormal_rng( log(C[i,]) , sigma ));  
      
  }
}