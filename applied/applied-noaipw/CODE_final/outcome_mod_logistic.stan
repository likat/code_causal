data {
  int<lower=0> N;
  int<lower=0> num_knot;
  int<lower=0> numpred;
  array[N] int<lower=0> y;  // 1 if positive, 0 if negative
  // vector<lower=0>[N] wts;
  matrix[N,numpred] pred_mat;
  matrix[N,num_knot] spline_mat;
}

parameters {
  real<lower=0.0001> sigma_spline;
  vector[numpred] beta;
  vector[num_knot] gamma;
}

transformed parameters {
  vector<lower=0, upper=1>[N] p;
  p = inv_logit(pred_mat*beta + spline_mat*gamma);
}

model {
  y ~ bernoulli(p);
  gamma ~ normal(0, sigma_spline);
}
