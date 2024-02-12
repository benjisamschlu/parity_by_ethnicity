data {
  int<lower=0> N;  // nber of cells
  int<lower=0> T_all;  // years to be estimated
  int<lower=0> A_all;  // nber of ages: 15-50
  real<lower=0> age[A_all];  // ages
  int<lower=0> G; // nber of subpopulation (races or states)
  int<lower=0> i[N]; // indexing variable
  int<lower=0> y[N]; // # mother
  int<lower=0> n[N];  // # women 
  vector[N] se; // SE of mother prop. on logit scale
}

parameters {
  
  // U
  matrix<lower=0.7, upper=1>[T_all, G] U;
  real<lower=0> sigma_U;
  
  // L
  matrix<lower=0, upper=0.05>[T_all, G] L;
  real<lower=0> sigma_L;
  
  // k
  matrix<lower=0, upper=1>[T_all, G] k;
  real<lower=0> sigma_k;
  
  // p
  matrix<lower=1>[T_all, G] c;
  real<lower=0> sigma_c;
  
  // data error
  vector[N] eps;
  
}

transformed parameters {
  
  matrix<lower=0, upper=1>[A_all, T_all] mu[G];
  matrix<lower=0, upper=1> [A_all * T_all, G] mu_ord;
  vector[N] phi;


 for (g in 1:G) {
    for (t in 1:T_all) {
      for (x in 1:A_all) {
        
        // Janoscheck sigmoidal function
      mu[g][x , t]  = U[t, g] - (U[t, g] - L[t, g]) * exp(-k[t, g] * pow(age[x], c[t, g]));
      }
    }
    // Reorder mu to be able to get a vector
    mu_ord[, g] = to_vector(mu[g]); // to_vector(): col major order
  }
  
  phi = logit(to_vector(mu_ord)[i]) + eps .* se;
}

model {
  
  // Likelihood 
  y ~ binomial(n, inv_logit(phi)); 
  
  
  // Priors
  
  // se depend on data point (combination of survey, group, year, and age)
  eps ~ normal(0, 1);
  
  for (g in 1:G) {
      // RW(2) with common variance parameter
      k[3:T_all, g] ~ normal(2 * k[2:(T_all-1), g] - k[1:(T_all-2), g], sigma_k);
      // RW(2) with common variance parameter
      U[3:T_all, g] ~ normal(2 * U[2:(T_all-1), g] - U[1:(T_all-2), g], sigma_U);
      // RW(2) with common variance parameter
      L[3:T_all, g] ~ normal(2 * L[2:(T_all-1), g] - L[1:(T_all-2), g], sigma_L);
      // RW(2) with common variance parameter
      c[3:T_all, g] ~ normal(2 * c[2:(T_all-1), g] - c[1:(T_all-2), g], sigma_c);
      
    }
  sigma_k ~ normal(0,0.1);
  sigma_U ~ normal(0,0.1);
  sigma_L ~ normal(0,0.1);
  sigma_c ~ normal(0, 1);
}
