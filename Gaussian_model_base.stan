 // Gaussian model fit - basic model for all aggergated data
 // Jess K Hopf
 // June 2023
 
 // Associated file: Bayes_analysis_vx.Rmd
 
 // Set up the model in stan language so that it can be complied and parsed to Rstan
 
 
 data {
  int<lower = 0> n;  // number of data points (reserves)
  vector[n] Y;       // values for each reserve (observed data)
  real theta;        // prior 
  real tau;          // prior 
  }
  
  parameters {
  real mu;                // mean 
  real<lower = 0> sigma;  //std dev (prob redundant as exp dist is >0)
  }
  
  model {
  Y ~ normal(mu,sigma^2);      // likelihood model
  mu ~ normal(theta, tau^2);   // prior for mu (from Lester data)
  sigma ~ exponential(1); //weakly informative prior. sigma is >0. // For exp dist: expected mean of stddev = 1/l, l = 1/stddev  /tau
  }


