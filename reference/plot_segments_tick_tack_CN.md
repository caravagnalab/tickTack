# Plot copy number segments from CNAqc

Plot copy number segments from CNAqc

## Usage

``` r
plot_segments_tick_tack_CN(
  cnaqc_x,
  K = K,
  max_alleles = 6,
  chromosomes = paste0("chr", c(1:22))
)
```

## Arguments

- cnaqc_x:

  A CNAqc object containing CNA data.

- K:

  Numeric, number of clusters.

- max_alleles:

  Numeric, maximum number of alleles to plot (default: 6).

- chromosomes:

  A character vector specifying which chromosomes to plot.

## Value

A ggplot object showing CNA segments.
