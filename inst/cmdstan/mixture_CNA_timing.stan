
data {
  int N;
  array[N] int NV;
  array[N] int DP;
  int n_groups;
  // real purity;
  // int multiplicities[n_groups];
}

parameters {
  ordered[n_groups] peaks;
  simplex[n_groups] omega;
}

model {
  vector[n_groups] contributions;
  // priors
  omega ~ dirichlet(rep_vector(2.0, n_groups));
  peaks ~ beta(.5,.5);
  // peaks = 1 + 2;

  // likelihood
  for (i in 1:N) {
    for (k in 1:n_groups) {
      contributions[k] = log(omega[k]) + binomial_lpmf(NV[i] | DP[i], peaks[k]);
    }
    target += log_sum_exp(contributions);
  }
}
