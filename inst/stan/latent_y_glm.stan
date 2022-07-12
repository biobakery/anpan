data {
  int<lower=0> N;
  int<lower=1> d;
  matrix[N,d] X;
  matrix[N,N] dist_mat;
  real<lower=0> y_scale;
  real<lower=0> sigma_factor;
}

transformed data{
  matrix[N,N] A = (-1/2.0) * (dist_mat.^2);
  matrix[N,N] C = identity_matrix(N) - ((1.0/N) * (ones_vector(N) * ones_row_vector(N)));
  matrix[N,N] G = C * A * C;
  real sst = trace(G);
}

parameters {
  vector[d] beta;
  real<lower=0> sigma;
  vector[N] y;
}

model {
  // real gdiff = sum(G - (y * y'));

  vector[N*N] gdiff = to_vector(G - (y*y'));

  y ~ normal(X * beta, sigma);
  gdiff ~ normal(0, .01);
  // I can't think of the proper way to enforce G = y*y'
  // Just set the elements of the difference to be very close to 0

  //priors

  beta ~ normal(0, (1.0 / 2) * y_scale);
  // sigma ~ exponential(1/y_scale);
  sigma ~ inv_gamma(sigma_factor + 1, sigma_factor * y_scale);
  y ~ normal(0, y_scale);

}

generated quantities{
  real ssh = ((X * beta)' * (X*beta));
  real ssr = sst - ssh;
  real F = (ssh / d) / (ssr / (N - d)); // There's no -1 because no intercept
  real R2 = 1 - ssr / sst;
  real adj_R2 = 1 - (1-R2)*(N)/(N-d);

}
