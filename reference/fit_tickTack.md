# Fit the tickTack Model for Timing CNA Events

Fits the tickTack variational inference model to a CNAqc object to
estimate the timing of copy number alteration (CNA) events. The function
smooths nearby segments, filters mutations and segments based on quality
criteria, and runs variational inference across a range of cluster
numbers `K`, selecting the best model via BIC.

## Usage

``` r
fit_tickTack(
  x,
  tolerance = 1e-04,
  possible_k = c("2:1", "2:2", "2:0"),
  alpha = 0.05,
  min_mutations_number = 10,
  max_distance_smooth = 5e+06,
  min_segment_length = 1e+06,
  n_components = NULL
)
```

## Arguments

- x:

  A CNAqc object containing somatic mutations and copy number segments,
  along with purity information. The object must have mutations
  accessible via `Mutations(x)` and CNAs via `CNA(x)`.

- tolerance:

  Numeric. Relative tolerance for the variational inference convergence
  criterion. Default is `1e-4`.

- possible_k:

  Character vector. Karyotypes to consider during filtering. Each entry
  should be a string of the form `"major:minor"` (e.g., `"2:1"`,
  `"2:2"`, `"2:0"`). Default is `c("2:1", "2:2", "2:0")`.

- alpha:

  Numeric. Significance level used when filtering mutations within
  expected VAF peaks. Default is `0.05`.

- min_mutations_number:

  Integer. Minimum number of mutations required for a segment to be
  included in the inference. Default is `10`.

- max_distance_smooth:

  Numeric. Maximum genomic distance (in base pairs) between adjacent
  segments for them to be merged during smoothing. Default is `5e6`.

- min_segment_length:

  Numeric. Minimum segment length (in base pairs) for a CNA segment to
  be retained. Default is `1e6`.

- n_components:

  Integer or `NULL`. Maximum number of timing clusters \\K\\ to
  evaluate. If `NULL` (default), the range is automatically set to
  `1:floor(sqrt(S))`, where `S` is the number of accepted segments.

## Value

The input CNAqc object `x`, augmented with a `results_timing` field.
This field is a named list containing:

- `data`:

  Filtered input data, including accepted CNA segments and the prepared
  Stan input list.

- `draws_and_summary`:

  Named list (by `K`) of variational inference draws and a per-segment
  summary tibble with columns `clock_mean`, `clock_median`, `clock_low`
  (5th percentile), and `clock_high` (95th percentile) of the assigned
  timing clock \\\tau\\.

- `log_lik_matrix_list`:

  Named list (by `K`) of log-likelihood draw matrices used for BIC
  computation.

- `elbo_iterations`:

  Named list (by `K`) of data frames tracking ELBO values across
  variational inference iterations.

- `bic_values`:

  Numeric vector of BIC values for each `K` in `range_k`.

- `best_K`:

  Integer. The value of `K` that minimises the BIC.

## Details

The function proceeds in the following steps:

1.  Nearby CNA segments are merged using
    [`smooth_segments`](https://rdrr.io/pkg/CNAqc/man/smooth_segments.html).

2.  Mutations and segments are filtered by karyotype, minimum mutation
    count, VAF peak membership, and minimum segment length.

3.  For each `K` in `range_k`, a tickTack Stan model is fit via ADVI
    (automatic differentiation variational inference).

4.  Each segment is hard-assigned to the timing clock \\\tau_k\\ with
    the highest posterior probability.

5.  Model selection is performed by computing BIC as \\-2 \cdot
    \bar{\ell} + (2K - 1) \cdot \log(N)\\, where \\\bar{\ell}\\ is the
    mean total log-likelihood across draws and \\N\\ is the total number
    of mutations.

6.  The `K` minimising BIC is stored as `best_K`.

If inference fails to converge for a given `K`, a warning is emitted and
results for that `K` are omitted.

## See also

[`smooth_segments`](https://rdrr.io/pkg/CNAqc/man/smooth_segments.html),
[`prepare_input_data`](https://caravagnalab.github.io/tickTack/reference/prepare_input_data.md),
`get_model`

## Examples

``` r
if (FALSE) { # \dontrun{
# x is a CNAqc object with purity, mutations, and CNA segments
x_timed <- fit_tickTack(
  x,
  tolerance          = 1e-4,
  possible_k         = c("2:1", "2:2", "2:0"),
  alpha              = 0.05,
  min_mutations_number = 10,
  max_distance_smooth  = 5e6,
  min_segment_length   = 1e6,
  n_components       = NULL
)

# Inspect the best number of clocks
x_timed$results_timing$best_K

# Per-segment timing summary for the best K
best <- x_timed$results_timing$best_K
x_timed$results_timing$draws_and_summary[[as.character(best)]]$summarized_results
} # }
```
