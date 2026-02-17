# get_taus_karyo Function

This function allows you to obtain simulated data for copy number events
(tau and karyo input for get_simulation or simulate_and_fit) inputting
the number of events and the proportions of taus and karyos.

## Usage

``` r
get_simulation_tickTack(
  number_clocks,
  number_events,
  purity,
  coverage,
  epsilon,
  seed = 123,
  vector_karyo = c("2:0", "2:1", "2:2"),
  time_interval = 7,
  w = 0.01,
  mu = 1e-04,
  l = 2e+07
)
```

## Arguments

- number_clocks:

  unique tau to be present in the simulation

- number_events:

  number of events to be simulated

- purity:

  sample purity

- coverage:

  average number of reads that align to a reference base

- epsilon:

  minimum distance between taus in terms of pseudo time

- seed:

  simulation seed

- vector_karyo:

  possible karyotypes to be included in the inference

- time_interval:

  time interval

- w:

  cell division rate

- mu:

  mutation rate

- l:

  length of the segments

## Value

list of simulated data and x object on which to perform the inference
