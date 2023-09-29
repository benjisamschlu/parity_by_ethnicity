data {
  int<lower=0> N;  // nber of cells
  int<lower=0> A;  // nber of age groups
  int<lower=0> Y_obs;  // nber of observed years
  vector<lower=0>[A] age;  // ages
  int<lower=0> G; // nber of groups (races or states)
  int<lower=0> n[N];  // women 
  int<lower=0> y[N]; // women childless
}

parameters {
  matrix<lower=0, upper=1>[Y_obs, G] L;
  matrix<lower=0, upper=45>[Y_obs, G] x0;
  matrix<lower=0, upper=1>[Y_obs, G] k;
  
}

transformed parameters {
  matrix[A, Y_obs] p[G]; 
  for (t in 1:Y_obs) {
    for (g in 1:G) {
      
      p[g][ , t]  = L[t, g] + (1 - L[t, g]) ./ (1+exp(k[t, g]*(age-x0[t, g])));  // non-zero asymptote logistic fct (4 pars)

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
  to_vector(L) ~ normal(0.2,0.1);
  to_vector(x0) ~ normal(30,5);
  to_vector(k) ~ normal(0.4,0.1);
  
}
