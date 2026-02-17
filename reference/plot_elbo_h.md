# Plot elbo Function

Plot the ELBO behaviour for each fit fixing the number of cluster.

## Usage

``` r
plot_elbo_h(elbo_iteration, K = 2)
```

## Arguments

- elbo_iteration:

  data.frame': \#iterations until convergence of elbo obs. of 2
  variables: iteration: int, elbo : num one of the element list
  elbo_iteration obtained from fit_h "results"

- K:

  number of components in the fit

## Value

p
