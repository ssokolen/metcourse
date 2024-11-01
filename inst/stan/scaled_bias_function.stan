// 
data {
  // Data points
  int<lower = 1> n;

  // Time points
  int<lower = 1> nt;

  // Time points to scale
  int scal[nt];
  int<lower = 1> sample[n];

  // Compounds
  int<lower = 1> nc;
  int<lower = 1> compound[n];

  // Basis
  int nb; 
  matrix[n, nb] B;

  // Response
  real y[n];
}

parameters {
  // B-spline parameters
  vector[nb] alpha;

  // Bias estimates
  real bias[nt];

  // Overall standard deviation of estimated response
  real<lower = 0> sigma[nc]; 
}

transformed parameters {
  // Scaling parameter
  vector[nt] S;

  // Estimated response
  vector[n] y_hat;

  y_hat = B*alpha;
  
  for (i in 1:nt) {
    if (scal[i] == 1) {
      S[i] = 1 + bias[i];
    }
    else {
      S[i] = 1;
    }
  }

  for (i in 1:n) {
    y_hat[i] = S[sample[i]]*y_hat[i];
  }
}

model {
  bias ~ normal(0.0, 0.1);
  for (i in 1:n) {
    y[i] ~ normal(y_hat[i], sigma[compound[i]]);
  }
}

