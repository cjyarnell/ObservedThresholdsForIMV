// for doing the seemingly-unrelated regressions for ObsCrit project
// based on https://mc-stan.org/docs/2_29/stan-users-guide/multivariate-outcomes.html

data{
  int<lower=0> N; // number of patients
  int<lower=1> D; // number of thresholds
  int<lower=0> K; // number of covariates
  int<lower=2> M; // number of categories (in this case 3)
  array[D] matrix[N,K] x; // predictors
  // outcome y, 1 = never met, 2 = met not followed, 3 = followed
  array[N,D] int<lower=1,upper=M> y; 
}

transformed data {
    int m = M-1;
}

parameters {
  array[m] vector[K] beta_raw; // covariates are same for every threshold
  array[m] row_vector[D] alpha; // intercepts vary
  array[m] cholesky_factor_corr[D] L_Omega; // covariance
  array[m] matrix[N,D] z_raw; //latent
}

transformed parameters {
  array[m] matrix[N,D] z;
  array[m] vector[K] beta = beta_raw;
  
  // zero out coefficients unique to meeting / following
   beta[1,17:22] = rep_vector(0, 6);
   beta[2,15:16] = rep_vector(0, 2);
  
    {
      for (i in 1:m){
        matrix[N,D] L_Oz = z_raw[i] * (L_Omega[i]') ;
        for (d in 1:D){
        //  
          z[i,,d] = rep_vector(alpha[i,d],N) + col(L_Oz,d) + x[d] * beta[i];    
        }
      }
    }
}

model {

// priors
  for (i in 1:m){
    L_Omega[i] ~ lkj_corr_cholesky(4);
    to_vector(z_raw[i]) ~ std_normal();
    beta_raw[i] ~ std_normal();
    alpha[i] ~ std_normal();
  }
  
  
  {
    vector[M] temp;
    for (n in 1:N){
      for (d in 1:D){
          // do this as a for loop if more categorical outcomes
          // we do it this way because we want met/not followed as baseline
          temp[1] = 0;
          temp[2] = z[1,n,d];
          temp[3] = z[1,n,d] + z[2,n,d];
          //
          target += categorical_logit_lupmf(y[n,d] | temp);
      }
    }
  }
}

generated quantities {
}