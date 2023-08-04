data {
  int<lower=0> N;  // nber of cells
  int<lower=0> A;  // nber of age groups
  int<lower=0> Y_obs;  // nber of observed years
  int<lower=0> Y;  // nber of years to estimate
  vector<lower=0>[A] age;  // ages
  int<lower=0> t_i[Y_obs]; // identify observed years: allows to estimate missing years
  int<lower=0> n[N];  // women 
  int<lower=0> y[N]; // women childless
}

parameters {
  vector<lower=0, upper=1>[Y] U;
  vector<lower=0, upper=1>[Y] L;
  vector<lower=15, upper=45>[Y] x0;
  vector<lower=0>[Y] k;
  
  real<lower=0> sigma_U;
  real<lower=0> sigma_L;
  real<lower=0> sigma_x0;
  real<lower=0> sigma_k;
}

transformed parameters {
  matrix[A, Y] p; 
  for (t in 1:Y) {
    
    p[ , t]  = L[t] + (U[t]-L[t])./(1+exp(k[t]*(age-x0[t])));  // non-zero asymptote logistic fct (4 pars)
  }
}

model {
  y ~ binomial(n, to_vector(p[,t_i])); // to_vector(): col major order
  
  // Priors
  U[3:Y] ~ normal(2 * U[2:(Y-1)] - U[1:(Y-2)], sigma_U);
  L[3:Y] ~ normal(2 * L[2:(Y-1)] - L[1:(Y-2)], sigma_L);
  x0[3:Y] ~ normal(2 * x0[2:(Y-1)] - x0[1:(Y-2)], sigma_x0);
  k[3:Y] ~ normal(2 * k[2:(Y-1)] - k[1:(Y-2)], sigma_k);
  
  [sigma_U, sigma_L] ~ normal(0,0.1);
  sigma_k ~ normal(0,1);
  sigma_x0 ~ normal(0, 10);
}
