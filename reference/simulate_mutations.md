# simulate_mutations Function

This function allows you to obtain sample mutations to use for CN Timing
inference.

## Usage

``` r
simulate_mutations(
  karyotype,
  time_interval,
  tau,
  l,
  mu,
  w,
  segment_id = "segment_id"
)
```

## Arguments

- karyotype:

  karyotype

- time_interval:

  time

- tau:

  tau generated

- l:

  length of the segment considered

- mu:

  mutation rate

- w:

  cell division rate

- segment_id:

  segment id

## Value

df with all the mutations associated with a simulated CN event
