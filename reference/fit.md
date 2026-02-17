# Fit timing of clonal peaks in cancer genome sequencing data

This function fits the timing of clonal peaks in cancer genome
sequencing data using either a beta-binomial or binomial model.

## Usage

``` r
fit(
  segments,
  mutations,
  purity,
  possible_k = c("2:1", "2:2", "2:0"),
  alpha = 0.05,
  min_mutations_number = 2,
  beta_binomial = FALSE,
  beta_binomial_disp = 0.01,
  tmp_file_path = NULL
)
```

## Arguments

- segments:

  A data frame containing segment information with columns `chr`,
  `from`, `to`, `Major`, and `minor`.

- mutations:

  A data frame containing mutation information with columns `chr`,
  `from`, `to`, `DP`, and `NV`.

- purity:

  A numeric value representing the tumor purity.

- possible_k:

  A character vector of possible karyotypes in the format "Major:minor".
  Default is c("2:1", "2:2", "2:0").

- alpha:

  A numeric value for the significance level. Default is 0.05.

- min_mutations_number:

  An integer specifying the minimum number of mutations required for
  analysis. Default is 2.

- beta_binomial:

  A logical value indicating whether to use the beta-binomial model.
  Default is FALSE.

- beta_binomial_disp:

  A numeric value for the beta-binomial dispersion parameter. Default is
  0.01.

- tmp_file_path:

  path of the directory where to save the temporary files during the
  inference

## Value

A list containing two tibbles: `inference_results` and
`summarized_results`. Returns NULL if no results are obtained.
