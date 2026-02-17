# get_taus_karyo gerstrung

This function allows you to obtain simulated data for copy number events
(tau and karyo input for get_simulation or simulate_and_fit) inputting
the number of events and the proportions of taus and karyos.

## Usage

``` r
get_taus_karyo_gerstrung(
  number_events,
  vector_tau,
  vector_karyo,
  weigths_tau,
  weights_karyo,
  chromosomes = c(1:22)
)
```

## Arguments

- number_events:

  number of events to be simulated

- vector_tau:

  unique tau to be present in the simulation

- vector_karyo:

  unique karyotype to be present in the simulation

- weigths_tau:

  proportion of tau to generate for each unique value

- weights_karyo:

  proportion of karyo to generate for each unique karyo

- chromosomes:

  possible chromosomes interested by the simulated CN events

## Value

list of simulated taus and karyo together with respective chromosomes
