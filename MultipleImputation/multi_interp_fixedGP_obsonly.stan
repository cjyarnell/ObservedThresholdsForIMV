//
// this stan program interpolates using the Gaussian process for one patient
// parameters found from MLE or the posterior mean of MCMC

data {
  int<lower=1> N_obs;
  int<lower=1> D; //  output dimensions (continuous first, then binary)
  
  // observation times
  array[N_obs] real x_obs;
  array[N_obs] int<lower=1,upper=D> d_obs;
  matrix[N_obs, D] d_obs_matrix;

  // continuous outcomes
  int<lower=0> N_c; // number of continuous observations
  int<lower=0> D_c; // number of continuous outcome dimensions 
  array[N_c] real y_c; // outcome values
 
  // binary outcomes 
  int<lower=0> N_b; // number of continuous observations
  int<lower=0> D_b; // number of binary outcome dimensions
  array[N_b] int<lower=0,upper=1> y_b; // outcome values
  
  // predictions
  int<lower=1> N_pred; // total hourly predictions
  array[N_pred] real x_pred; // prediction times
  array[N_pred] int<lower=1,upper=D> d_pred; // dimension of each prediction
  
  // Gaussian process
  real<lower=0> rho;
  vector[D] alpha;
  vector[D_c] sigma;
  vector[D_b] offset;
  matrix[D,D] L_C;
}

transformed data{
  
  real delta = 1e-6;

  int N = N_obs + N_pred;
  array[N] real x;
  
  for (n in 1:N_obs){
    x[n] = x_obs[n];
  }
  for (n in 1:N_pred){
    x[N_obs+n] = x_pred[n];
  }

  // precalculated covariance matrix
  matrix[N,N] L_K;
  
  L_K =gp_exp_quad_cov(x, 1.0, rho);
  
    for (n in 1:N) {
      L_K[n, n] = L_K[n, n] + delta;
    }
  
  L_K = cholesky_decompose(L_K);
  
  matrix[D,D] alpha_L_C = diag_pre_multiply(alpha, L_C)';
  
  // for use in model block
  vector[D] d_ones = rep_vector(1.0, D);
  
  // extended offset vector
  
  vector[D] offset2 = append_row(rep_vector(0.0, D_c), offset);
  
}

parameters {
  
  matrix[N,D] eta;

}

transformed parameters{

    matrix[N,D]   f = L_K * eta * alpha_L_C;
    vector[N_c] f_c = block(f,1,1,N_c,D) 
                   .* block(d_obs_matrix,1,1,N_c,D) 
                    * d_ones;
    vector[N_b] f_b = block(f,N_c+1,1,N_b,D) 
                   .* block(d_obs_matrix, N_c+1,1,N_b,D) 
                    * d_ones;
}

model {
  
  to_vector(eta) ~ std_normal();
  
  // use d_ones and .* to isolate observed latent variables
  target += normal_lupdf(y_c | f_c , sigma[d_obs[1:N_c]]) + 
            bernoulli_logit_lpmf(y_b | f_b + offset2[d_obs[(N_c+1):N_obs]]);
  
}

generated quantities{
  
  array[N_pred] real y_pred;
  {
      for (n in 1:N_pred){
        if(d_pred[n] < 5){
            y_pred[n] = normal_rng(f[N_obs + n,d_pred[n]], sigma[d_pred[n]]);
        }
        if (d_pred[n] == 6){
            y_pred[n] = normal_rng(f[N_obs + n,d_pred[n]], sigma[d_pred[n]]);
        } 
        if (d_pred[n] == 10){
            y_pred[n] = bernoulli_logit_rng(f[N_obs + n, d_pred[n]] + offset2[d_pred[n]]);
        } 
     }
  }
}
