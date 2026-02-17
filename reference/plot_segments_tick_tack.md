# Plot segmentation results from tickTack analysis

Plot segmentation results from tickTack analysis

## Usage

``` r
plot_segments_tick_tack(x, colour_by = "clock_mean", K = 1)
```

## Arguments

- x:

  An object containing segmentation results.

- colour_by:

  A character string specifying the variable to color segments by
  (default: 'clock_mean').

- K:

  Numeric, number of clusters (default: 1).

## Value

A ggplot object visualizing tickTack segmentation.
