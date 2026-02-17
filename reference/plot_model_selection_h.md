# Plot the behavior of the model selection scores vs the number of clusters

Plot the BIC, Log Likelihood, ICL scores obtatined after the fit for
each number of cluster.

## Usage

``` r
plot_model_selection_h(model_selection_tibble, best_K)
```

## Arguments

- model_selection_tibble:

  a tibble with 3 scores and k_max values, one for each inference

- best_K:

  the integer corresponding to the best number of components

## Value

model_selection_plot
