# model_selection Function

Perform model selection among the models fit with varying number of
mixture components (number of clocks).

## Usage

``` r
get_k_max_k_means(input_data, purity, alpha = 0.05)
```

## Arguments

- input_data:

  input data

- purity:

  sample purity

- alpha:

  description

## Value

optimal_k: suggested number from c-means
