// multioutput GP with Hilbert space approximation, 
// multiple patients multiple covariates

functions {
  
  real partial_sum_joint_lpdf(real[] slice_k,
                        int start,
                        int end,
                        real[] y_c,
                        int[] y_b,
                        int[] d_obs,
                        matrix d_obs_matrix,
                        matrix[] PHI,
                        vector diagSPD,
                        int M,
                        matrix[] eta,
                        matrix L_C,
                        vector alpha,
                        vector sigma,
                        vector mu_hat,
                        int[] N_k,
                        int[] N_kc,
                        int[] N_kb,
                        int D,
                        vector d_ones
                        ) {
   // first row of first patient relative to whole dataset
    int startn = sum(N_k[1:start])- N_k[start] + 1;
    int startc = sum(N_kc[1:start])- N_kc[start] + 1;
    int startb = sum(N_kb[1:start])- N_kb[start] + 1;
    
    int length = sum(N_k[start:end]);
    int lengthc = sum(N_kc[start:end]);
    int lengthb = sum(N_kb[start:end]);

    matrix[D,end-start+1] mu;
    vector[D] mu_temp;
    vector[length] f_temp;
    array[lengthc] int d_obs_temp;
    array[lengthc] int fc_ind;
    array[lengthb] int fb_ind;
    matrix[M,D] SPD_eta;    
    int pos = 0;
    int posc = 0;
    int posb = 0;
 
    for (k in start:end){
  
    // mean    
      mu_temp = mu_hat;
      
    // HSGP
      SPD_eta = diag_pre_multiply(diagSPD, eta[k]);

      f_temp[(pos + 1):(pos + N_k[k])] = mu_temp[segment(d_obs, pos + startn, N_k[k])]
                + (((PHI[k,1:N_k[k],] * SPD_eta * diag_pre_multiply(alpha, L_C)')
                .* block(d_obs_matrix, pos + startn, 1, N_k[k], D)) * d_ones);
                
      d_obs_temp[(posc + 1):(posc + N_kc[k])] = segment(d_obs, pos + startn, N_kc[k]);
      
      for (nc in 1:N_kc[k]){fc_ind[posc + nc] = nc + pos;}
      for (nb in 1:N_kb[k]){fb_ind[posb + nb] = nb + N_kc[k] + pos;}

      pos  = pos + N_k[k];
      posc = posc + N_kc[k];
      posb = posb + N_kb[k];
    }
    // joint likelihood
      return normal_lupdf(segment(y_c, startc, lengthc) | f_temp[fc_ind], sigma[d_obs_temp]) 
           + bernoulli_logit_lupmf(segment(y_b, startb, lengthb) | f_temp[fb_ind]);
 } 

// from https://github.com/gabriuma/basis_functions_approach_to_GP/blob/master/uni_dimensional/simulated_dataset/stancode_BF_1dim.stan
 real lambda(real L, int m) {
		real lam;
		lam = ((m*pi())/(2*L))^2;
				
		return lam;
	}
	real spd(real alpha, real rho, real w) {
		real S;
		S = (alpha^2) * sqrt(2*pi()) * rho * exp(-0.5*(rho^2)*(w^2));
				
		return S;
	}
	vector phi(real L, int m, vector x) {
		vector[rows(x)] fi;
		fi = 1/sqrt(L) * sin(m*pi()/(2*L) * (x+L));
				
		return fi;
	}
}
data {
  int<lower=1> N_obs;
  int<lower=1> K; // number of patients - train set first, test set second
  int<lower=1> D; //  output dimensions (continuous first, then binary)
  
  // number of observations per patient
  // organize data sequentially by patient
  array[K] int<lower=1> N_k;
  array[K] int<lower=0> N_kc; // continuous obs
  array[K] int<lower=0> N_kb; // binary obs

  // observation times
  array[N_obs] real x_obs;
  array[N_obs] int<lower=1,upper=D> d_obs;
  matrix[N_obs, D] d_obs_matrix;
  array[N_obs] int<lower=1,upper=K> k_obs;
  
  // continuous outcomes
  int<lower=0> N_c; // number of continuous observations
  int<lower=0> D_c; // number of continuous outcome dimensions 
  array[N_c] real y_c; // outcome values
 
  // binary outcomes 
  int<lower=0> N_b; // number of continuous observations
  int<lower=0> D_b; // number of binary outcome dimensions
  array[N_b] int<lower=0,upper=1> y_b; // outcome values
 
  // for Hilbert space basis function approximation
  int M; // number of basis functions
  real L; // boundary condition factor
  
}

transformed data {
  
// make dummy arrays for reduce sum
  array[K] real id_real = rep_array(1.0, K);
  
  vector[D] d_ones = rep_vector(1.0, D);
  
  real delta = 1e-5;
// make an array of matrices big enough to accommodate each patient
  array[K] int<lower=0> phi_size;

// build a counter for observations by patient
  array[N_obs] int obs_by_k;
  {
    int posk = 0;
    for (k in 1:K){
      for (n in 1:N_k[k]){
        obs_by_k[posk+n] = n;
      }
      posk = posk + N_k[k];
    }
  }

// build time-based covariance matrices  

   for (k in 1:K){
    phi_size[k] = N_k[k];
    }

  array[K] matrix[max(phi_size), M] PHI;

  // build the approximate matrices patient by patient
  {
    int pos = 1;
    for (k in 1:K){
      for (m in 1:M){
        PHI[k,1:(phi_size[k]),m] = 
          phi(L, m, segment(to_vector(x_obs), pos, N_k[k]));
      }
      pos = pos + N_k[k];
    }
  }
// define offsets and multipliers
  // rho
  real logrho_offset = -2;
  real logrho_sd     = 0.1;
  
  // alpha
  real logalpha1_offset = -0.5;
  real logalpha1_sd = 0.5;
  real logalpha2_offset = -0.5;
  real logalpha2_sd = 0.5;
  real logalpha3_offset = -0.5;
  real logalpha3_sd = 0.5;
  real logalpha4_offset = -0.5;
  real logalpha4_sd = 0.5;
  real logalpha5_offset = -0.5;
  real logalpha5_sd = 0.5;
  real logalpha6_offset = -0.5;
  real logalpha6_sd = 0.5;
  real logalpha7_offset = -0.5;
  real logalpha7_sd = 0.5;
  real logalpha8_offset = -0.5;
  real logalpha8_sd = 0.5;
  real logalpha9_offset = -0.5;
  real logalpha9_sd = 0.5;
  real logalpha10_offset = -0.5;
  real logalpha10_sd = 0.5;
  real logalpha11_offset = 0.7;
  real logalpha11_sd = 0.5;
  real logalpha12_offset = 0.64;
  real logalpha12_sd = 0.5;
  real logalpha13_offset = 1.2;
  real logalpha13_sd = 0.5;
  real logalpha14_offset = 0.42;
  real logalpha14_sd = 0.5;
  
  // noise
  real logsigma1_offset = -0.34;
  real logsigma1_sd = 0.5;
  real logsigma2_offset = -0.35;
  real logsigma2_sd = 0.5;
  real logsigma3_offset = -35;
  real logsigma3_sd = 0.5;
  real logsigma4_offset = -0.34;
  real logsigma4_sd = 0.5;
  real logsigma5_offset = -0.35;
  real logsigma5_sd = 0.5;
  real logsigma6_offset = -0.35;
  real logsigma6_sd = 0.5;
  real logsigma7_offset = -0.35;
  real logsigma7_sd = 0.5;
  real logsigma8_offset = -0.35;
  real logsigma8_sd = 0.5;
  real logsigma9_offset =  -0.35;
  real logsigma9_sd = 0.5;
  real logsigma10_offset =  -0.35;
  real logsigma10_sd = 0.5;


  // offset
  real offset1_offset = -1.4;
  real offset1_sd = 0.5;
  real offset2_offset = -3.5;
  real offset2_sd = 0.5;
  real offset3_offset = -2.45;
  real offset3_sd = 0.5;
  real offset4_offset = -3;
  real offset4_sd = 0.5;

  // continuous global means (centered already)
  vector[D_c] zeros = rep_vector(0.0, D_c);
  
  // prior sd
  
  real prior_alpha_contsd = 0.5;
  real prior_alpha_binsd = 1;
}

parameters {

  // HSGP 
  real<offset=logrho_offset,multiplier = logrho_sd> logrho;
  
  real<offset=logalpha1_offset,multiplier=logalpha1_sd> logalpha1;
  real<offset=logalpha2_offset,multiplier=logalpha2_sd> logalpha2;
  real<offset=logalpha3_offset,multiplier=logalpha3_sd> logalpha3;
  real<offset=logalpha4_offset,multiplier=logalpha4_sd> logalpha4;
  real<offset=logalpha5_offset,multiplier=logalpha5_sd> logalpha5;
  real<offset=logalpha6_offset,multiplier=logalpha6_sd> logalpha6;
  real<offset=logalpha7_offset,multiplier=logalpha7_sd> logalpha7;
  real<offset=logalpha8_offset,multiplier=logalpha8_sd> logalpha8;
  real<offset=logalpha9_offset,multiplier=logalpha9_sd> logalpha9;
  real<offset=logalpha10_offset,multiplier=logalpha10_sd> logalpha10;
  real<offset=logalpha11_offset,multiplier=logalpha11_sd> logalpha11;
  real<offset=logalpha12_offset,multiplier=logalpha12_sd> logalpha12;
  real<offset=logalpha13_offset,multiplier=logalpha13_sd> logalpha13;
  real<offset=logalpha14_offset,multiplier=logalpha14_sd> logalpha14;


  array[K] matrix[M,D] eta;
  
  // noise
  real<offset=logsigma1_offset, multiplier=logsigma1_sd> logsigma1;
  real<offset=logsigma2_offset, multiplier=logsigma2_sd> logsigma2;
  real<offset=logsigma3_offset, multiplier=logsigma3_sd> logsigma3;
  real<offset=logsigma4_offset, multiplier=logsigma4_sd> logsigma4;
  real<offset=logsigma5_offset, multiplier=logsigma5_sd> logsigma5;
  real<offset=logsigma6_offset, multiplier=logsigma6_sd> logsigma6;
  real<offset=logsigma7_offset, multiplier=logsigma7_sd> logsigma7;
  real<offset=logsigma8_offset, multiplier=logsigma8_sd> logsigma8;
  real<offset=logsigma9_offset, multiplier=logsigma9_sd> logsigma9;
  real<offset=logsigma10_offset, multiplier=logsigma10_sd> logsigma10;

  // offset
  real<offset=offset1_offset, multiplier=offset1_sd> offset1;
  real<offset=offset2_offset, multiplier=offset2_sd> offset2;
  real<offset=offset3_offset, multiplier=offset3_sd> offset3;
  real<offset=offset4_offset, multiplier=offset4_sd> offset4;

  // correlation matrices
  cholesky_factor_corr[D] L_C; // correlation matrix across outputs

}

transformed parameters{
  real<lower=0> rho =exp(logrho);
  vector<lower=0>[D] alpha = exp([ logalpha1
                                 , logalpha2
                                 , logalpha3
                                 , logalpha4
                                 , logalpha5
                                 , logalpha6
                                 , logalpha7
                                 , logalpha8
                                 , logalpha9
                                 , logalpha10
                                 , logalpha11
                                 , logalpha12
                                 , logalpha13
                                 , logalpha14
                                 ]');

  // noise and offset 
  vector<lower=0>[D_c] sigma = exp([
      logsigma1
    , logsigma2
    , logsigma3
    , logsigma4
    , logsigma5
    , logsigma6
    , logsigma7
    , logsigma8
    , logsigma9
    , logsigma10
  ]');
  
// noise and offset 
  vector[D_b] offset = [
      offset1
    , offset2
    , offset3
    , offset4
  ]';
  
   vector[M] diagSPD;
    for (m in 1:M){
      diagSPD[m] = sqrt(spd(1.0, rho, sqrt(lambda(L, m))));
    }
    
    vector[D] mu_hat = append_row(zeros, offset);
    
}

model {

// PRIORS

  // HSGP
  logrho ~ normal(logrho_offset,logrho_sd); 
  
  logalpha1  ~ normal(logalpha1_offset , prior_alpha_contsd);
  logalpha2  ~ normal(logalpha2_offset , prior_alpha_contsd);
  logalpha3  ~ normal(logalpha3_offset , prior_alpha_contsd);
  logalpha4  ~ normal(logalpha4_offset , prior_alpha_contsd);
  logalpha5  ~ normal(logalpha5_offset , prior_alpha_contsd);
  logalpha6  ~ normal(logalpha6_offset , prior_alpha_contsd);
  logalpha7  ~ normal(logalpha7_offset , prior_alpha_contsd);
  logalpha8  ~ normal(logalpha8_offset , prior_alpha_contsd);
  logalpha9  ~ normal(logalpha9_offset , prior_alpha_contsd);
  logalpha10  ~ normal(logalpha10_offset , prior_alpha_contsd);
  logalpha11  ~ normal(logalpha11_offset , prior_alpha_binsd);
  logalpha12  ~ normal(logalpha12_offset , prior_alpha_binsd);
  logalpha13  ~ normal(logalpha13_offset , prior_alpha_binsd);
  logalpha14  ~ normal(logalpha14_offset , prior_alpha_binsd);

  for(k in 1:K){
    to_vector(eta[k]) ~ std_normal();
  }

  // noise in continuous variables
  logsigma1  ~ normal(logsigma1_offset , logsigma1_sd);
  logsigma2  ~ normal(logsigma2_offset , logsigma2_sd);
  logsigma3  ~ normal(logsigma3_offset , logsigma3_sd);
  logsigma4  ~ normal(logsigma4_offset , logsigma4_sd);
  logsigma5  ~ normal(logsigma5_offset , logsigma5_sd);
  logsigma6  ~ normal(logsigma6_offset , logsigma6_sd);
  logsigma7  ~ normal(logsigma7_offset , logsigma7_sd);
  logsigma8  ~ normal(logsigma8_offset , logsigma8_sd);
  logsigma9  ~ normal(logsigma9_offset , logsigma9_sd);
  logsigma10  ~ normal(logsigma10_offset , logsigma10_sd);

  // offset in binary variables
  offset1  ~ normal(offset1_offset , offset1_sd);
  offset2  ~ normal(offset2_offset , offset2_sd);
  offset3  ~ normal(offset3_offset , offset3_sd);
  offset4  ~ normal(offset4_offset , offset4_sd);


  // correlation across variables
  L_C ~ lkj_corr_cholesky(3);

  
// MODEL       
  

    target += reduce_sum(partial_sum_joint_lpdf,
                        id_real, 1,
                        y_c, y_b,
                        d_obs, d_obs_matrix,
                        PHI, diagSPD, M, eta, L_C, alpha, sigma,
                        mu_hat,
                        N_k, N_kc, N_kb, D, d_ones);
}

generated quantities{
  
}

