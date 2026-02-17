# Extract mutations.

Getter to obtain mutation calls from an object.

## Usage

``` r
Mutations(x, cna = c("clonal", "subclonal"), type = c("SNV", "indel"))
```

## Arguments

- x:

  A CNAqc object.

- cna:

  `"clonal"` for clonal CNAs, `"subclonal"` for subclonal CNAs.

- type:

  `"SNV"` for single-nucleotide variants, `"indel"` for
  insertion-deletions.

## Value

A tibble with the data.
