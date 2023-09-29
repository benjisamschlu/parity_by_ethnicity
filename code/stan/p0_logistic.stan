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
  real<lower=0> sigma_x0;
  
  // k
  matrix<lower=0, upper=1>[Y_obs, G] k;
  real<lower=0> sigma_k;

  // L
  matrix<lower=0, upper=1>[Y_obs, G] U;
  real<lower=0> sigma_U;
  
}

transformed parameters {
  
  matrix<lower=0, upper=1>[A, Y_obs] theta[G];
  matrix<lower=0, upper=1>[A, Y_obs] p[G]; 


 for (g in 1:G) {
    for (t in 1:Y_obs) {
      
      // logistic fct (3 pars)
      theta[g][ , t]  = U[t, g] ./ (1 + exp(-k[t, g] * (age - x0[t, g])));
      
      // complement of logistic fct
      p[g][, t] = 1 - theta[g][, t];
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
      x0[3:Y_obs, g] ~ normal(2 * x0[2:(Y_obs-1), g] - x0[1:(Y_obs-2), g], sigma_x0);
    }
  sigma_x0 ~ normal(0,1);

  // k  
  for (g in 1:G) {
      // RW(1) with common variance parameter
      k[2:Y_obs, g] ~ normal(k[1:(Y_obs-1), g], sigma_k);
    }
  sigma_k ~ normal(0,0.1);
  
  // L  
  for (g in 1:G) {
      // RW(1) with common variance parameter
      U[2:Y_obs, g] ~ normal(U[1:(Y_obs-1), g], sigma_U);
    }
  sigma_U ~ normal(0,0.1);

}
