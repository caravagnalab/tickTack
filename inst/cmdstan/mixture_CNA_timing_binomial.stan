
data {
  int N;
  array[N] int NV;
  array[N] int DP;

  array[2] real peaks;
  // real beta_dispersion;
  // real purity;
  // int multiplicities[n_groups];
}

parameters {
  simplex[2] omega;
}

model {
  vector[2] contributions;
  // priors
  omega ~ dirichlet(rep_vector(2.0, 2));

  // likelihood
  for (i in 1:N) {
    for (k in 1:2) {
      contributions[k] = log(omega[k]) + binomial_lpmf(NV[i] | DP[i], peaks[k]);
    }
    target += log_sum_exp(contributions);
  }
}
