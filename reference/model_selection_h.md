# model_selection Function

Perform model selection among the models fit with varying number of
mixture components (number of clocks).

## Usage

``` r
model_selection_h(results, n_components = 0)
```

## Arguments

- results:

  list of 4: \$data, \$draws_and_summary, \$log_lik_matrix_list and
  \$elbo_iterations

- n_components:

  number of components specified from user

## Value

result_model_selection: list \$best_fit, \$best_K,
\$model_selection_tibble, \$entropy_list)
