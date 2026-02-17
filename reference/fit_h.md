# Get the fit to obtain the clocks posteriors

Obtain the approximate posterior of the clocks for each model fit with
the number of components up to k_max.

## Usage

``` r
fit_h(
  x,
  local_executable = FALSE,
  max_attempts = 2,
  INIT = TRUE,
  tolerance = 1e-04,
  possible_k = c("2:1", "2:2", "2:0"),
  alpha = 0.05,
  min_mutations_number = 4,
  n_components = 0,
  initial_iter = 200,
  grad_samples = 10,
  elbo_samples = 200,
  tmp_file_path = NULL,
  cmd_version_old = FALSE,
  eta = NULL,
  adapt_engaged = FALSE,
  adapt_iter = NULL,
  algorithm = NULL
)
```

## Arguments

- x:

  list: A CNAqc object.

- local_executable:

  FALSE

- max_attempts:

  num: max number of repeated inference for ADVI

- INIT:

  logical: boolean variable to set the initialization phase to TRUE or
  FALSE

- tolerance:

  num: tolerance in the ELBO optimization procedure

- possible_k:

  chr: "2:1" "2:2" "2:0"

- alpha:

  num: (type double) confidence interval level to choose the data that
  fall in the expected binomial intervals

- min_mutations_number:

  num: (type double) minimum number of accepted mutations for a segment
  to be included in the inference

- n_components:

  number of components specified from user

- initial_iter:

  description

- grad_samples:

  description

- elbo_samples:

  description

- tmp_file_path:

  output_dir getOption("cmdstanr_output_dir") path of the directory
  where to save the temporary files with the info on the elbo
  evaluations during the VI inference

- cmd_version_old:

  version of cmdstanr for the draws parameter in variational method

- eta:

  NULL

- adapt_engaged:

  FALSE

- adapt_iter:

  NULL

- algorithm:

  NULL

## Value

results_and_data = list(data = input_data_list, results = results,
output_files_list = output_files_list)
