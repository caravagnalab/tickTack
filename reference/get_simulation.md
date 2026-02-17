# get_simulation Function

This function allows you to obtain simulated data for copy number
events.

## Usage

``` r
get_simulation(
  taus,
  karyotypes,
  purity = 0.9,
  time_interval = 20,
  l = 1e+07,
  mu = 1e-04,
  w = 0.01,
  coverage = 100
)
```

## Arguments

- taus:

  = a vector of doubles between 0-1

- karyotypes:

  = a vector of the same length of taus of strings of the type "2:1",
  "2:0", "2:2"

- purity:

  sample purity

- time_interval:

  time interval not in 0-1, real number

- l:

  length of the segment considered

- mu:

  mutation rate

- w:

  cell division rate

- coverage:

  average number of reads that align to a reference base

## Value

all the mutations associated with the simulated CN events
