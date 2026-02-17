# plot_timing_h Function

This function obtains the list of input data to be used in the stan
model.

## Usage

``` r
plot_timing_h(
  results,
  K,
  colour_by = "karyotype",
  split_contiguous_segments = TRUE
)
```

## Arguments

- results:

  list(data = input_data_list, results = results, output_files_list =
  output_files_list)

- K:

  mun: number of clocks

- colour_by:

  chr: default = "karyotype"

- split_contiguous_segments:

  option to plot segments' setalarion lines

## Value

p : plot of inference results with credibility intervals in the
chromosome absolute positions
