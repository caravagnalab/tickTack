# Extract driver mutations from mutation data

Extract driver mutations from mutation data

## Usage

``` r
get_drivers(x, chromosomes = paste0("chr", c(1:22, "X", "Y")), which = "VAF")
```

## Arguments

- x:

  A CNAqc object containing mutation data.

- chromosomes:

  A character vector of chromosome names to filter for (default: 'chr1'
  to 'chr22', 'X', 'Y').

- which:

  A character string specifying whether to extract 'VAF' or 'CCF'
  values.

## Value

A data frame of driver mutations filtered by chromosome and data type.
