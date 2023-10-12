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
  
  // k
  matrix<lower=0, upper=1>[Y_obs, G] k;

  // U
  matrix<lower=0, upper=1>[Y_obs, G] U;
  
}

transformed parameters {
  
  matrix<lower=0, upper=1>[A, Y_obs] theta[G];
  matrix<lower=0, upper=1>[A, Y_obs] p[G]; 


 for (g in 1:G) {
    for (t in 1:Y_obs) {
      
      // logistic fct (3 pars)
      theta[g][ , t]  = U[t, g] ./ (1 + exp(-k[t, g] * (age - x0[t, g])));
      
      // add splines on real scale
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
  
  
  // EDA Priors
  
  // x0
  to_vector(x0) ~ normal(0, 1);
  to_vector(k) ~ normal(0, 1);
  to_vector(U) ~ normal(0, 1);

}
