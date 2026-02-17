# Plot timing of clonal peaks in cancer genome sequencing data

This function generates a plot showing the timing of clonal peaks in
cancer genome sequencing data.

## Usage

``` r
plot_timing(fit_results, segments, colour_by = "karyotype", ref = "GRCh38")
```

## Arguments

- fit_results:

  A list containing the results of the `fit_timing` function,
  specifically `summarized_results`.

- segments:

  A data frame containing segment information with columns `chr`,
  `from`, `to`, `Major`, and `minor`.

- colour_by:

  A character string specifying the variable to color the plot by.
  Default is "karyotype".

- ref:

  Reference genome desired. Either 'GRCh38' or 'hg19'

## Value

A ggplot object showing the timing of clonal peaks.
