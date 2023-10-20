data {
  int<lower=0> N;  // nber of cells
  int<lower=0> A;  // nber of age groups
  int<lower=0> Y_obs;  // nber of observed years
  real<lower=0> age[A];  // ages
  int<lower=0> G; // nber of subpopulation (races or states)
  int<lower=0> n[N];  // women 
  int<lower=0> y[N]; // women childless
}

parameters {
  
  // U
  matrix<lower=0.7, upper=1>[Y_obs, G] U;
  real<lower=0> sigma_U;
  
  // L
  matrix<lower=0, upper=0.05>[Y_obs, G] L;
  real<lower=0> sigma_L;
  
  // k
  matrix<lower=0, upper=1>[Y_obs, G] k;
  real<lower=0> sigma_k;
  
  // p
  matrix<lower=1>[Y_obs, G] c;
  real<lower=0> sigma_c;
  
}

transformed parameters {
  
  matrix<lower=0, upper=1>[A, Y_obs] p[G];


 for (g in 1:G) {
    for (t in 1:Y_obs) {
      for (x in 1:A) {
        
        // Janoscheck
      p[g][x , t]  = U[t, g] - (U[t, g] - L[t, g]) * exp(-k[t, g] * pow(age[x], c[t, g]));
      }
    }
  }
}

model {
  
  {
    matrix[A*Y_obs, G] p_ord;
    
    for (g in 1:G) {
      p_ord[, g] = to_vector(p[g]);
    }
    y ~ binomial(n, to_vector(p_ord)); // to_vector(): col major order
  }
  
  
  // Priors
  
  // x0
  for (g in 1:G) {
      // RW(2) with common variance parameter
      k[3:Y_obs, g] ~ normal(2 * k[2:(Y_obs-1), g] - k[1:(Y_obs-2), g], sigma_k);
      // RW(2) with common variance parameter
      U[3:Y_obs, g] ~ normal(2 * U[2:(Y_obs-1), g] - U[1:(Y_obs-2), g], sigma_U);
      // RW(2) with common variance parameter
      L[3:Y_obs, g] ~ normal(2 * L[2:(Y_obs-1), g] - L[1:(Y_obs-2), g], sigma_L);
      // RW(2) with common variance parameter
      c[3:Y_obs, g] ~ normal(2 * c[2:(Y_obs-1), g] - c[1:(Y_obs-2), g], sigma_c);
      
    }
  sigma_k ~ normal(0,0.1);
  sigma_U ~ normal(0,0.1);
  sigma_L ~ normal(0,0.1);
  sigma_c ~ normal(0, 1);
}
