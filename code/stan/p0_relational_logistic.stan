data {
  int<lower=0> N;  // nber of cells
  int<lower=0> A;  // nber of age groups
  int<lower=0> Y_obs;  // nber of observed years
  vector<lower=0>[A] age;  // ages
  int<lower=0> G; // nber of subpopulation (races or states)
  int<lower=0> K; // nber of knots
  matrix[A, K] B; // linear splines basis
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

  // U
  matrix<lower=0, upper=1>[Y_obs, G] U;
  real<lower=0> sigma_U;
  
  // L
  matrix<lower=0, upper=1>[Y_obs, G] L;
  real<lower=0> sigma_L;
  
  // intercept
  matrix[Y_obs, G] a;
  real<lower=0> sigma_a;
  
  // slope
  matrix[Y_obs, G] b;
  real<lower=0> sigma_b;
  
}

transformed parameters {
  
  matrix<lower=0, upper=1>[A, Y_obs] theta[G];
  matrix[A, Y_obs] omega[G]; 
  matrix<lower=0, upper=1>[A, Y_obs] p[G]; 


 for (g in 1:G) {
    for (t in 1:Y_obs) {
      
      // logistic fct (4 pars)
      theta[g][ , t]  = L[t, g] + ((U[t, g] - L[t,g]) ./ (1 + exp(k[t, g] * (age - x0[t, g]))));
      
      // linear transformation
      omega[g][, t] = a[t, g] + b[t, g] * logit(theta[g][, t]);
    }
  }
  // convert omega to [0,1] scale
  p = inv_logit(omega);
  
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
      // RW(2) with common variance parameter
      k[3:Y_obs, g] ~ normal(2 * k[2:(Y_obs-1), g] - k[1:(Y_obs-2), g], sigma_k);
      // RW(2) with common variance parameter
      U[3:Y_obs, g] ~ normal(2 * U[2:(Y_obs-1), g] - U[1:(Y_obs-2), g], sigma_U);
      // RW(2) with common variance parameter
      L[3:Y_obs, g] ~ normal(2 * L[2:(Y_obs-1), g] - L[1:(Y_obs-2), g], sigma_L);
      // RW(2) with common variance parameter
      a[3:Y_obs, g] ~ normal(2 * a[2:(Y_obs-1), g] - a[1:(Y_obs-2), g], sigma_a);
      // RW(2) with common variance parameter
      b[3:Y_obs, g] ~ normal(2 * b[2:(Y_obs-1), g] - b[1:(Y_obs-2), g], sigma_b);
      
    }
  sigma_x0 ~ normal(0,1);
  sigma_k ~ normal(0,0.1);
  sigma_U ~ normal(0,0.1);
  sigma_L ~ normal(0,0.1);
  sigma_a ~ normal(0, 0.1);
  sigma_b ~ normal(0, 0.1);
  
 
}
