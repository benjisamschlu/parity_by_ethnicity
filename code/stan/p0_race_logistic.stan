data {
  int<lower=0> N;  // nber of cells
  int<lower=0> A;  // nber of age groups
  int<lower=0> Y_obs;  // nber of observed years
  vector<lower=0>[A] age;  // ages
  int<lower=0> R; // nber of races
  int<lower=0> n[N];  // women 
  int<lower=0> y[N]; // women childless
}

parameters {
  matrix<lower=0, upper=1>[Y_obs, R] L;
  matrix<lower=0, upper=1>[Y_obs, R] k;
  matrix<lower=0, upper=45>[Y_obs, R] x0;
  vector<lower=0, upper=45>[Y_obs] MU_x0;
  vector[R] eta;
  
  real<lower=0> sigma_x0;
  real<lower=0> sigma_MU_x0;
}

transformed parameters {
  
  matrix<lower=0, upper=1>[A, Y_obs] p[R];
  
  
  
  for (r in 1:R) {
    for (t in 1:Y_obs) {
      
      p[r][ , t]  = L[t, r] + (1 - L[t, r]) ./ (1 + exp(k[t, r] * (age - x0[t, r])));  // non-zero asymptote logistic fct (4 pars)

    }
  }
}

model {
  
  {
    matrix[A*Y_obs, R] p_ord;
    
    for (r in 1:R) {
      p_ord[, r] = to_vector(p[r]);
    }
    y ~ binomial(n, to_vector(p_ord)); // to_vector(): col major order
  }
  
  
  // Priors
  to_vector(L) ~ normal(0,1);
  to_vector(k) ~ normal(0,1);
  
  for (r in 1:R) {
      
      x0[, r] ~ normal(MU_x0 + eta[r], sigma_x0);
    }
  
  // RW2 on MU_x0
  MU_x0[3:Y_obs] ~ normal(2 * MU_x0[2:(Y_obs-1)] - MU_x0[1:(Y_obs-2)], sigma_MU_x0);
  
  // Race specific effect around RW2
  eta ~ normal(0,1); 
  
  sigma_x0 ~ normal(0,1);
  sigma_MU_x0 ~ normal(0,1);

}
