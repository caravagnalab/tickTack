# prepare_input_data Function

This function obtains the list of input data to be used in the stan
model.

## Usage

``` r
prepare_input_data(
  mutations,
  segments,
  purity,
  possible_k = c("2:1", "2:2", "2:0"),
  alpha = 0.05,
  min_mutations_number = 2
)
```

## Arguments

- mutations:

  list: output from CNAqc Mutations(x), where x is a CNAqc object

- segments:

  list: output from CNAqc CNA(x), where x is a CNAqc object

- purity:

  num: (type = double) sample purity between 0 and 1

- possible_k:

  chr: "2:1" "2:2" "2:0"

- alpha:

  num: (type double) confidence interval level to choose the data that
  fall in the expected binomial intervals

- min_mutations_number:

  num: (type double) minimum number of accepted mutations for a segment
  to be included in the inference

## Value

accepted_data: list \$input_data:list: List of 7: \$S: int, \$N: int,
\$karyotype: num (0 or 1), \$seg_assignment: num, \$peaks:List of N of
num (1:2), \$NV: num, \$DP: num \$accepted_cna: tibble (S × 5)(S3:
tbl_df/tbl/data.frame): \$segment_original_indx: int, \$segment_name:
chr, \$segment_id: num, \$karyotype: chr, \$chr: chr
