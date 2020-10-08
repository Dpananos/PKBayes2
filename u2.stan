functions{
  
  vector heaviside(vector t){
    return (t + fabs(t)) ./ (2 * t);
  }
  


  vector conc(real D, vector t, real Cl, real ka, real ke){
    
    return heaviside(t) .* (exp(-ka*t) - exp(-ke*t)) * (0.5 * D * ke * ka ) / (Cl *(ke - ka));
  }
}
data{
  int n;
  vector[n] t;
  vector[3] dose_times;
}
generated quantities{
  vector[n] C = rep_vector(0.0, n);
  
  for(i in 1:3){
    C+=conc(1, t-dose_times[i], 3.0, 1.1, 0.4);
  }
}