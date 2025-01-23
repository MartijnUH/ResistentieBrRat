// Direct modeling of proportions using a Dirichlet-Multinomial framework
functions {
  real partial_sum(array[] int site, int start, int end, int T, int K,
                   array[,] int y, array[,] int index, array[,] vector prop) {
    real lp = 0;
    for (i in start:end) {
      for (t in 1:T) {
        if(sum(y[index[i, t]]) > 0 ){
          lp += multinomial_lpmf(y[index[i, t]] | to_vector(prop[i, t, ]));
        }
      }
    }
    return lp;
  }
  
  vector lambda_nD(array[] real L, array[] int m, int D) {
    vector[D] lam;
    for(i in 1:D){
      lam[i] = ((m[i]*pi())/(2*L[i]))^2; }
      
      return lam;
  }
  
  real spd_nD_onelscale(real alpha, real rho, vector w, int D) {
    real S;
    S = alpha^2 * sqrt(2*pi())^D * rho^D * 
        exp(-0.5*rho^2 * (w' * w));
    
    return S;
  }

  vector phi_nD(array[] real L, array[] int m, matrix x) {
    int c = cols(x);
    int r = rows(x);

    matrix[r, c] fi;
    vector[r] fi1;
    for (i in 1:c) {
      fi[, i] = 1 / sqrt(L[i]) * sin(m[i] * pi() * (x[, i] + L[i]) / (2 * L[i]));
    }
    fi1 = fi[, 1];
    for (i in 2:c) {
      fi1 = fi1 .* fi[, i];
    }
    return fi1;
  }
  
  vector HSGP_2D(int N, matrix phi, real rho_space, real alpha_space, 
                 vector eta_space, array[] real L, array[,] int indices, int M_nD) {
    vector[N] f;
    vector[M_nD] diagSPD;

    for (m in 1:M_nD) {
      diagSPD[m] =  sqrt(spd_nD_onelscale(alpha_space, rho_space, sqrt(lambda_nD(L, indices[m,], 2)), 2));
    }
    f = phi * (diagSPD .* eta_space);
    return f;
  }
  
  real generalized_inverse_gaussian_lpdf(real x, int p, real a, real b) {
    return p * 0.5 * log(a / b)
      - log(2 * modified_bessel_second_kind(p, sqrt(a * b)))
      + (p - 1) * log(x)
      - (a * x + b / x) * 0.5;
  }
}

data {
  int<lower=1> N;                           // Number of presence points
  int<lower=1> R;                           // Number of grid cells
  int<lower=1> T;                           // Number of time points
  int<lower=1> K;                           // Number of categories

  array[N, K] int<lower=0> y;               // Observed counts
  simplex[K] mu;                            // Observed proportions in t==1
  array[R, T] int<lower=1, upper=N> index;  // pointer for value corresponding to i,t

  // Space
  int n_bf;                                 // Number of basis functions
  array[n_bf, 2] int xy_design;             // Design matrix for basis functions
  matrix[R, 2] xy_std;                      // Standardized coordinates
  array[2] real xy_lim;                     // Spatial domain limits

  // Time
  row_vector[T] year;
}

transformed data {
  array[R] int site = rep_array(0, R);      // Dummy for site index
  array[K] real logit_mu = to_array_1d(logit(mu)); // E on logit scale
  matrix[R, n_bf] phi;                      // Reduced dimensional space
  int N0 = 0;                               // N with all zeroes
  to_array_1d
  for (i in 1:R) {
    for (t in 1:T) {
      
      if (sum(y[index[i, t]]) == 0) {
        (y[index[i, t]] + 1e-6) / sum(y[index[i, t]] + 1e-6)
        N0 += 1;
      }
    }
  }
  for (n in 1:n_bf) {
    phi[, n] = phi_nD(xy_lim, xy_design[n, ], xy_std);
  }
}

parameters {
  array[K] real bs0;                          // Slope for time effect per mutation
  array[K] real<lower=0> alpha_u;             // Marginal variance HSGP per mutation
  array[K] real<lower=0> rho_u;               // Lengthscale HSGP per mutation
  array[K] vector[n_bf] eta_u;                // Loadings HSGP per mutation
  array[2] real<lower=0> alpha_b;             // Marginal variance HSGP per mutation
  array[2] real<lower=0> rho_b;               // Lengthscale HSGP per mutation
  array[2] vector[n_bf] eta_b;                // Loadings HSGP per mutation
}

transformed parameters {
  array[K] vector[R] U;                       // Spatial random effects
  array[K] vector[R] Bs;                      // Spatial random slopes
  array[K] vector[R] B;                       // Slope coefficients
  array[R, T] simplex[K] prop;                // Proportions of each mutation (softmax)

  for (k in 1:K) {
    B[k] = rep_vector(bs0[k], R);
    U[k] = HSGP_2D(R, phi, rho_u[k], alpha_u[k], eta_u[k], xy_lim, xy_design, n_bf);
    if (k == 1 || k == 2) {
      Bs[k] = HSGP_2D(R, phi, rho_b[k], alpha_b[k], eta_b[k], xy_lim, xy_design, n_bf);
      B[k] += Bs[k];
    }
  }

  for (i in 1:R) {
    for (t in 1:T) {
      vector[K] logits;
      for (k in 1:K)
        logits[k] = logit_mu[k] + (B[k][i] * year[t]) + U[k][i];
      prop[i, t, ] = softmax((0.999 * logits) + 1e-6);
    }
  }
}

model {
  int grainsize = 1;

  // Priors
  to_vector(bs0) ~ student_t(3, 0, 5);
  to_vector(alpha_u) ~ student_t(3, 0, 5);
  to_vector(alpha_b) ~ student_t(3, 0, 5);

  for (k in 1:K) {
    rho_u[k] ~ generalized_inverse_gaussian(3, 11, 0.1);
    eta_u[k] ~ normal(0, 1);
    sum(U[k]) ~ normal(0, 0.001 * R);
  
    if (k == 1 || k == 2) {
      rho_b[k] ~ generalized_inverse_gaussian(3, 15, 0.05);
      eta_b[k] ~ normal(0, 1);
      sum(B[k]) ~ normal(0, 0.001 * R);
    }
  }

  // Partial likelihood
  target += reduce_sum(partial_sum, site, grainsize, T, K, y, index, prop);
}

generated quantities { 
  array[N, K] int y_new;                     // Simulated counts
  real fit = 0;                              // Overall fit statistic for observed data
  real fit_new = 0;                          // Overall fit statistic for simulated data
  array[N-N0] real log_lik;                  // Log-likelihood for each observation
  
  {
    int counter = 0;
    array[R, T, K] real E = rep_array(0, R, T, K); 
    array[R, T, K] real chisq_obs = rep_array(0, R, T, K);  
    array[R, T, K] real chisq_sim = rep_array(0, R, T, K);  
    
    for (i in 1:R) {
      for (t in 1:T) {
        if (sum(y[index[i, t]]) > 0) {
          counter += 1;
          // Simulate new counts based on proportions
          y_new[counter, ] = multinomial_rng(prop[i, t, ], sum(y[index[i, t]]));

          // Compute expected counts based on proportions
          for (k in 1:K) {
            E[i, t, k] = prop[i, t, k] * sum(y[index[i, t]]);
          }
          
          // Compute chi-squared discrepancy for observed data
          for (k in 1:K) {
            chisq_obs[i, t, k] = square(y[index[i, t]][k] - E[i, t, k]) / (E[i, t, k] + 0.5);
            chisq_sim[i, t, k] = square(y_new[counter, k] - E[i, t, k]) / (E[i, t, k] + 0.5);
            
            // Add to the overall fit statistics
            fit += chisq_obs[i, t, k];
            fit_new += chisq_sim[i, t, k];
          }
          log_lik[counter] = multinomial_lpmf(y[index[i, t]] | to_vector(prop[i, t, ]));
        }
      }
    }
  }
}

