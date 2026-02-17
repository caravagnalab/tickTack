# Plot

Plot

## Usage

``` r
plot_cnaqc(
  x,
  chromosomes = paste0("chr", c(1:22)),
  add_mobster = FALSE,
  max_Y_height = 6,
  cn = "absolute",
  highlight = x$most_prevalent_karyotype,
  highlight_QC = FALSE
)
```

## Arguments

- x:

  A CNAqc object containing mutation data.

- chromosomes:

  A character vector of chromosome names to filter for (default: 'chr1'
  to 'chr22', 'X', 'Y').

- add_mobster:

  TRUE or FALSE

- max_Y_height:

  heiht plot

- cn:

  type of cn

- highlight:

  which karyotype to highlight

- highlight_QC:

  FALSE

## Value

Plot.
