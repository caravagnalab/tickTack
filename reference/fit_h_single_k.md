# Get the fit to obtain the clocks posteriors

Obtain the approximate posterior of the clocks for each model fit with
the number of components up to k_max.

## Usage

``` r
fit_h_single_k(
  x,
  single_k,
  max_attempts = 2,
  INIT = TRUE,
  tolerance = 1e-04,
  possible_k = c("2:1", "2:2", "2:0"),
  alpha = 0.05,
  min_mutations_number = 2,
  n_components = 0,
  initial_iter = 200,
  grad_samples = 10,
  elbo_samples = 200
)
```

## Arguments

- x:

  list: A CNAqc object.

- single_k:

  perform a single inference with only K components

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

## Value

results_and_data = list(data = input_data_list, results = results,
output_files_list = output_files_list)
