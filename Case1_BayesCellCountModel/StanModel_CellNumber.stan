data {
  array[30] int C_bl;
  array[30] int C_w4;
  array[28] int C_end;
  array[30] int Treatment;
}

parameters {
  // GRand means for three timepoints
  real<lower=0> mu0_bl; 
  real<lower=0> mu0_w4; 
  real<lower=0> mu0_end;
  
  // Coefficients for drugs at three timepoints
  vector[5] bT_bl; 
  vector[5] bT_w4;
  vector[5] bT_end;
  
  // Hyperparameters for bT
  real mu_T_bl;
  real mu_T_w4;
  real mu_T_end;
  real<lower=0> sigma_T_bl;
  real<lower=0> sigma_T_w4;
  real<lower=0> sigma_T_end;
  
  // Overdispersion parameter of negative binomial distribution
  real<lower=0> thi_bl; 
  real<lower=0> thi_w4; 
  real<lower=0> thi_end; 
}

transformed parameters{
  // Set vehicle group as 0
  vector[6] bT_adj_bl;
  vector[6] bT_adj_w4;
  vector[6] bT_adj_end;
  
  bT_adj_bl = append_row(bT_bl, 0.0);
  bT_adj_w4 = append_row(bT_w4, 0.0);
  bT_adj_end = append_row(bT_end, 0.0);
}

model {
  //Computed mean for negative binomial distribution
  vector[30] mu_bl;
  vector[30] mu_w4;
  vector[28] mu_end;
  
  //Priors for contrl (vehicle) mean. Share the same priors
  mu0_bl ~ uniform(0, 17.7);
  mu0_w4 ~ uniform(0, 17.7);
  mu0_end ~ uniform(0, 17.7);
  
  //Hyperpriors for treatment coefficients
  mu_T_bl ~ normal(0, 7); 
  mu_T_w4 ~ normal(0, 7);
  mu_T_end ~ normal(0, 7);
  sigma_T_bl ~ normal(0, 2); 
  sigma_T_w4 ~ normal(0, 2);
  sigma_T_end ~ normal(0, 2);
  
  //Priors for treatment coefficients
  bT_bl ~ normal(mu_T_bl, sigma_T_bl); 
  bT_w4 ~ normal(mu_T_w4, sigma_T_w4);
  bT_end ~ normal(mu_T_end, sigma_T_end);
  
  //Priors for overdispersion
  thi_bl ~ normal(0, 1);
  thi_w4 ~ normal(0, 1);
  thi_end ~ normal(0, 1);
  
  //Compute means. One timepoint (loop) at a time.
  //(1) Baseline
  for(i in 1:30){
    mu_bl[i] = mu0_bl + bT_adj_bl[Treatment[i]];
    mu_bl[i] = exp(mu_bl[i]);
  }
  //(2) Week4
  for(i in 1:30){
    mu_w4[i] = mu0_w4 + bT_adj_w4[Treatment[i]];
    mu_w4[i] = exp(mu_w4[i]);
  }
  //(2) Endpoint
  for(i in 1:28){
    mu_end[i] = mu0_end + bT_adj_end[Treatment[i]];
    mu_end[i] = exp(mu_end[i]);
  }
  
  //Likelihood function
  C_bl ~ neg_binomial_2(mu_bl, thi_bl);
  C_w4 ~ neg_binomial_2(mu_w4, thi_w4);
  C_end ~ neg_binomial_2(mu_end, thi_end);
}

generated quantities{
  real n_inj; //Number injected
  vector[6] n_bl; //Posterior cell number in one femur at baseline
  vector[6] n_w4; //Posterior cell number in one femur at week4
  vector[6] n_end; //Posterior cell number in whole bone marrow at endpoint
  
  n_inj = poisson_rng(1e6); //Injected 1e6 CD3 depleted cells
  
  for(i in 1:6){
    n_bl[i] = neg_binomial_2_rng(exp(mu0_bl + bT_adj_bl[i]), thi_bl);
  }
  for(i in 1:6){
    n_w4[i] = neg_binomial_2_rng(exp(mu0_w4 + bT_adj_w4[i]), thi_w4);
  }
  for(i in 1:6){
    n_end[i] = neg_binomial_2_rng(exp(mu0_end + bT_adj_end[i]), thi_end);
  }
}
