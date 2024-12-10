## ----eval=FALSE---------------------------------------------------------------
#  # Example pseudo-code for filtering
#  probs <- c(alpha / 2, 1 - alpha / 2)
#  for (p in peaks) {
#    quantiles <- qbinom(probs, DP, p) # or qbb for beta-binomial
#    if (NV >= quantiles[1] && NV <= quantiles[2]) {
#      # Mutation is accepted
#    }
#  }

## ----eval=FALSE---------------------------------------------------------------
#  data {
#    int N;
#    array[N] int NV;
#    array[N] int DP;
#  
#    array[2] real peaks;
#  }

## ----eval=FALSE---------------------------------------------------------------
#  parameters {
#    simplex[2] omega;
#  }

## ----eval=FALSE---------------------------------------------------------------
#  model {
#    vector[2] contributions;
#  
#    omega ~ dirichlet(rep_vector(2.0, 2));
#  
#    for (i in 1:N) {
#      for (k in 1:2) {
#        contributions[k] = log(omega[k]) + binomial_lpmf(NV[i] | DP[i], peaks[k]);
#      }
#      target += log_sum_exp(contributions);
#    }
#  }

