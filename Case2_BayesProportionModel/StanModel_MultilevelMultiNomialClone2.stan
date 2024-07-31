//Multinomial regression modelled as a series of poisson likelihood
data{
  //Number of rows in d_CloneWide
  int<lower=0> N;
  
  //Number of Treatments passed into model
  int<lower=0> N_Treatment;
  
  // Number of cells in clones
  array[N] int Clone0;
  array[N] int Clone1;
  array[N] int Clone2;
  array[N] int Clone3;
  array[N] int Clone4;
  array[N] int Clone5;
  array[N] int Clone6;
  array[N] int Clone7;
  // array[N] int Clone8;
  // array[N] int Clone9;
  
  //Treatment
  array[N] int Treatment;
  
  //Sample
  // array[N] int Sample;
}

parameters{
  vector[N_Treatment] b_Clone0;
  vector[N_Treatment] b_Clone1;
  vector[N_Treatment] b_Clone2;
  vector[N_Treatment] b_Clone3;
  vector[N_Treatment] b_Clone4;
  vector[N_Treatment] b_Clone5;
  vector[N_Treatment] b_Clone6;
  vector[N_Treatment] b_Clone7;
  // vector[N_Treatment] b_Clone8;
  // vector[N_Treatment] b_Clone9;
  
  //hyperpriors
  vector[N_Treatment] mu0_Clone0;
  vector[N_Treatment] mu0_Clone1;
  vector[N_Treatment] mu0_Clone2;
  vector[N_Treatment] mu0_Clone3;
  vector[N_Treatment] mu0_Clone4;
  vector[N_Treatment] mu0_Clone5;
  vector[N_Treatment] mu0_Clone6;
  vector[N_Treatment] mu0_Clone7;
  // vector[N_Treatment] mu0_Clone8;
  // vector[N_Treatment] mu0_Clone9;
  
  //hyperpriors
  vector<lower=0>[N_Treatment] sigma_Clone0;
  vector<lower=0>[N_Treatment] sigma_Clone1;
  vector<lower=0>[N_Treatment] sigma_Clone2;
  vector<lower=0>[N_Treatment] sigma_Clone3;
  vector<lower=0>[N_Treatment] sigma_Clone4;
  vector<lower=0>[N_Treatment] sigma_Clone5;
  vector<lower=0>[N_Treatment] sigma_Clone6;
  vector<lower=0>[N_Treatment] sigma_Clone7;
  // vector<lower=0>[N_Treatment] sigma_Clone8;
  // vector<lower=0>[N_Treatment] sigma_Clone9;
}

model{
  //Vector storing lambda
  vector[N] mu_Clone0;
  vector[N] mu_Clone1;
  vector[N] mu_Clone2;
  vector[N] mu_Clone3;
  vector[N] mu_Clone4;
  vector[N] mu_Clone5;
  vector[N] mu_Clone6;
  vector[N] mu_Clone7;
  // vector[N] mu_Clone8;
  // vector[N] mu_Clone9;
  
  mu0_Clone0 ~ normal(0, 10);
  mu0_Clone1 ~ normal(0, 10);
  mu0_Clone2 ~ normal(0, 10);
  mu0_Clone3 ~ normal(0, 10);
  mu0_Clone4 ~ normal(0, 10);
  mu0_Clone5 ~ normal(0, 10);
  mu0_Clone6 ~ normal(0, 10);
  mu0_Clone7 ~ normal(0, 10);
  // mu0_Clone8 ~ normal(0, 10);
  // mu0_Clone9 ~ normal(0, 10);
  
  sigma_Clone0 ~ gamma(3, 1);
  sigma_Clone1 ~ gamma(3, 1);
  sigma_Clone2 ~ gamma(3, 1);
  sigma_Clone3 ~ gamma(3, 1);
  sigma_Clone4 ~ gamma(3, 1);
  sigma_Clone5 ~ gamma(3, 1);
  sigma_Clone6 ~ gamma(3, 1);
  sigma_Clone7 ~ gamma(3, 1);
  // sigma_Clone8 ~ gamma(3, 1);
  // sigma_Clone9 ~ gamma(3, 1);
  
  //Priors in the log scale
  b_Clone0 ~ normal(mu0_Clone0, sigma_Clone0);
  b_Clone1 ~ normal(mu0_Clone1, sigma_Clone1);
  b_Clone2 ~ normal(mu0_Clone2, sigma_Clone2);
  b_Clone3 ~ normal(mu0_Clone3, sigma_Clone3);
  b_Clone4 ~ normal(mu0_Clone4, sigma_Clone4);
  b_Clone5 ~ normal(mu0_Clone5, sigma_Clone5);
  b_Clone6 ~ normal(mu0_Clone6, sigma_Clone6);
  b_Clone7 ~ normal(mu0_Clone7, sigma_Clone7);
  // b_Clone8 ~ normal(mu0_Clone8, sigma_Clone8);
  // b_Clone9 ~ normal(mu0_Clone9, sigma_Clone9);
  
  for(i in 1:N){
    mu_Clone0[i] = b_Clone0[Treatment[i]];
    mu_Clone1[i] = b_Clone1[Treatment[i]];
    mu_Clone2[i] = b_Clone2[Treatment[i]];
    mu_Clone3[i] = b_Clone3[Treatment[i]];
    mu_Clone4[i] = b_Clone4[Treatment[i]];
    mu_Clone5[i] = b_Clone5[Treatment[i]];
    mu_Clone6[i] = b_Clone6[Treatment[i]];
    mu_Clone7[i] = b_Clone7[Treatment[i]];
    // mu_Clone8[i] = b_Clone8[Treatment[i]];
    // mu_Clone9[i] = b_Clone9[Treatment[i]];

  }
  
  Clone0 ~ poisson_log(mu_Clone0);
  Clone1 ~ poisson_log(mu_Clone1);
  Clone2 ~ poisson_log(mu_Clone2);
  Clone3 ~ poisson_log(mu_Clone3);
  Clone4 ~ poisson_log(mu_Clone4);
  Clone5 ~ poisson_log(mu_Clone5);
  Clone6 ~ poisson_log(mu_Clone6);
  Clone7 ~ poisson_log(mu_Clone7);
  // Clone8 ~ poisson_log(mu_Clone8);
  // Clone9 ~ poisson_log(mu_Clone9);
}

generated quantities{
  //Calculate posterior samples. Number of cells in each of the clone.
  vector[N_Treatment] n_Clone0; 
  vector[N_Treatment] n_Clone1; 
  vector[N_Treatment] n_Clone2; 
  vector[N_Treatment] n_Clone3; 
  vector[N_Treatment] n_Clone4; 
  vector[N_Treatment] n_Clone5; 
  vector[N_Treatment] n_Clone6; 
  vector[N_Treatment] n_Clone7; 
  // vector[N_Treatment] n_Clone8; 
  // vector[N_Treatment] n_Clone9; 
  
  for(i in 1:N_Treatment){
    n_Clone0[i] = poisson_log_rng(b_Clone0[i]);
    n_Clone1[i] = poisson_log_rng(b_Clone1[i]);
    n_Clone2[i] = poisson_log_rng(b_Clone2[i]);
    n_Clone3[i] = poisson_log_rng(b_Clone3[i]);
    n_Clone4[i] = poisson_log_rng(b_Clone4[i]);
    n_Clone5[i] = poisson_log_rng(b_Clone5[i]);
    n_Clone6[i] = poisson_log_rng(b_Clone6[i]);
    n_Clone7[i] = poisson_log_rng(b_Clone7[i]);
    // n_Clone8[i] = poisson_log_rng(b_Clone8[i]);
    // n_Clone9[i] = poisson_log_rng(b_Clone9[i]);
  }
}
