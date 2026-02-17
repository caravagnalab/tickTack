# fit_variational Function

This function performs the inference using the ADVI algorithm. Repeat
inference if it fails and repeat to avoid local minima, taking the best
run.

## Usage

``` r
fit_variational_h(
  input_data,
  local_executable = FALSE,
  purity,
  max_attempts = 2,
  initialization = NULL,
  INIT = TRUE,
  initial_iter = 100,
  grad_samples = 200,
  elbo_samples = 200,
  tolerance = 0.01,
  tmp_file_path = NULL,
  cmd_version_old = FALSE,
  eta = NULL,
  adapt_engaged = FALSE,
  adapt_iter = NULL,
  algorithm = NULL
)
```

## Arguments

- input_data:

  list: List of 7: \$S: int, \$N: int, \$karyotype: num (0 or 1),
  \$seg_assignment: num, \$peaks:List of N of num (1:2), \$NV: num,
  \$DP: num

- local_executable:

  FALSE

- purity:

  sample purity

- max_attempts:

  num: max number of repeated inference for ADVI

- initialization:

  list: List of 4: \$w: num (1:S, 1:3), \$tau: num (1:K), \$phi: num
  (1:K), \$kappa: num

- INIT:

  logical: boolean variable to set the initialization phase to TRUE or
  FALSE

- initial_iter:

  description

- grad_samples:

  description

- elbo_samples:

  description

- tolerance:

  num: tolerance in the ELBO optimization procedure

- tmp_file_path:

  output_dir getOption("cmdstanr_output_dir") path of the directory
  where to save the temporary files during the inference

- cmd_version_old:

  version of cmdstanr for the draws parameter invariational method

- eta:

  NULL

- adapt_engaged:

  FALSE

- adapt_iter:

  NULL

- algorithm:

  NULL

## Value

best_fit
