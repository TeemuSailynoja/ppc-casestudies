data {
  int<lower=2> K;                   // Number of classes
  int<lower=0> N;                   // Total number of observations
  array[N] int<lower=1, upper=K> y; // Target classes
  vector[N] x;                      // Observed predictor values
  real<lower=0> sigma;              // User supplied standard deviation
  int correct_sigma;                // How to handle sigma, see below.
}

transformed data{
  vector[N] x_st;
  real Sigma;

  // Standardize data
  x_st = (x - mean(x)) / sd(x);

  // Maybe remember to scale sigma accordingly
  if (correct_sigma == 1) {
    Sigma = sigma / sd(x);
  } else {
    Sigma = sigma;
  }
}

parameters {
  // Inferred means.
  ordered[K] c;
  simplex[K] p_c;
}

model {
  // Prior
  c ~ normal(0,1);
  p_c ~ dirichlet(rep_vector(1,K));

  // Likelihood
  for (n in 1:N) {
    target += normal_lpdf(x_st[n] | c[y[n]], Sigma);
  }
}

generated quantities {
  // Posterior predictive sample
  vector[N] yrep;
  // For each observation, posterior predictive mass of classes
  array[N] vector[K] ppm;

  for (n in 1:N) {
    for (k in 1:K) {
      ppm[n, k] = normal_lpdf(x_st[n] | c[k], Sigma);
    }
    ppm[n, ] = softmax(ppm[n, ]);
    yrep[n] = categorical_rng(ppm[n, ]);
  }
}
