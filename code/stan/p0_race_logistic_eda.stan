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
  matrix<lower=0, upper=1>[Y_obs, R] U;
  matrix<lower=0, upper=1>[Y_obs, R] L;
  matrix<lower=0, upper=45>[Y_obs, R] x0;
  matrix<lower=0>[Y_obs, R] k;
  
  real<lower=0> sigma_U;
  real<lower=0> sigma_L;
  real<lower=0> sigma_x0;
  real<lower=0> sigma_k;
}

transformed parameters {
  matrix[A, Y_obs] p[R]; 
  for (t in 1:Y_obs) {
    for (r in 1:R) {
      
      p[r][ , t]  = L[t, r] + (U[t, r]-L[t, r])./(1+exp(k[t, r]*(age-x0[t, r])));  // non-zero asymptote logistic fct (4 pars)

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
  to_vector(U) ~ normal(0,1);
  to_vector(L) ~ normal(0,1);
  to_vector(x0) ~ normal(0,1);
  to_vector(k) ~ normal(0,1);
  
}
