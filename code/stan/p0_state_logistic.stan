data {
  int<lower=0> N;  // nber of cells
  int<lower=0> A;  // nber of age groups
  int<lower=0> Y_obs;  // nber of observed years
  vector<lower=0>[A] age;  // ages
  int<lower=0> G; // nber of subpopulation (races or states)
  int<lower=0> n[N];  // women 
  int<lower=0> y[N]; // women childless
}

parameters {
  
  // x0
  matrix<lower=0, upper=45>[Y_obs, G] x0;
  vector<lower=0, upper=45>[Y_obs] MU_x0;
  real<lower=0> sigma_x0;
  real<lower=0> sigma_MU_x0;
  vector[G] eta;
  
  // k
  real<lower=0, upper = 1> mu_k;
  real<lower=0> sigma_k;
  vector[G] eps_k;

  // L
  real<lower=0, upper = 1> mu_L;
  real<lower=0> sigma_L;
  vector[G] eps_L;
  
}

transformed parameters {
  
  vector[G] L;
  vector[G] k;
  
  matrix<lower=0, upper=1>[A, Y_obs] p[G];
  
  k = mu_k + sigma_k * eps_k;
  
  L = mu_L + sigma_L * eps_L;
  
  for (g in 1:G) {
    for (t in 1:Y_obs) {
      
      p[g][ , t]  = L[g] + (1 - L[g]) ./ (1 + exp(k[g] * (age - x0[t, g])));  // non-zero asymptote logistic fct (4 pars)

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
      
      x0[, g] ~ normal(MU_x0 + eta[g], sigma_x0);
    }
  sigma_x0 ~ normal(0,1);
  // Group specific effect around RW2
  eta ~ normal(0,1); 
  // RW2 on MU_x0
  MU_x0[1:2] ~ normal(30,5);
  MU_x0[3:Y_obs] ~ normal(2 * MU_x0[2:(Y_obs-1)] - MU_x0[1:(Y_obs-2)], sigma_MU_x0);
  sigma_MU_x0 ~ normal(0,1);

  // k  
  mu_k ~ normal(0.35, 0.025);
  sigma_k ~ normal(0,0.05);
  to_vector(eps_k) ~ normal(0,1);
  
  // L  
  mu_L ~ normal(0.2, 0.025);
  sigma_L ~ normal(0,0.05);
  to_vector(eps_L) ~ normal(0,1);
}
