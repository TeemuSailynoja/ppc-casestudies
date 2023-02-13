
  data {
    int N; // number of observations
    int D; // number of features
    int N_classes; // number of classes
    matrix [N, D] X; // observation data
    array[N] int <lower = 1, upper = N_classes> y; // target values {1,..., N_classes}
  }
  
  transformed data {
    matrix[D + 1, N] X_stn;
    X_stn[D + 1, ] = rep_row_vector(1, N);
    for (d in 1:D) {
      X_stn[d,] = to_row_vector((X[, d] - mean(X[, d])) / sd(X[, d]));
    }
  }
  
  parameters {
    matrix[N_classes, D + 1] W;
  }
  
  transformed parameters {
    matrix[N_classes, N] Beta;
    for (c in 1:N_classes) {
      Beta[c, ] =  W[c, ] * X_stn;
    }
  }
  
  model {
    for (d in 1:(D + 1)) {
      for (c in 1:N_classes) {
        target += normal_lpdf(W[c, d] | 0, 1);
      }
    }
    for (n in 1:N) {
      target += categorical_logit_lpmf(y[n] | Beta[,n]);
    }
  }
  
  generated quantities {
    vector[N] yrep;
    for (n in 1:N) {
      yrep[n] = categorical_logit_rng(Beta[,n]);
    }
    matrix[N,N_classes] lpd;
    for (n in 1:N) {
      for (c in 1:N_classes) {
        lpd[n, c] = categorical_logit_lpmf(c | Beta[,n]);
      }
    }
  }
  
