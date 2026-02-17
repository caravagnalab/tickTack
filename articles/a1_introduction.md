# 1. Introduction

``` r
library(tickTack)
require(dplyr)
#> Loading required package: dplyr
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
```

The `tickTack` package provides tools for timing the occurrence of copy
number alterations (CNA). This vignette introduces the structure of the
input data required for `tickTack` functions, using the example dataset
[`tickTack::pcawg_example`](https://caravagnalab.github.io/tickTack/reference/pcawg_example.md).

## Input Data Structure

The input data consists of three components stored as named elements in
a list: `mutations`, `cna`, and `metadata`. Below is a description of
each component.

### 1. Mutations

The `mutations` component is a tibble containing information about
somatic mutations. Each row represents a mutation, with the following
columns:

- **`chr`**: Chromosome where the mutation occurs.
- **`from`** and **`to`**: Start and end positions of the mutation on
  the chromosome.
- **`ref`** and **`alt`**: Reference and alternate alleles.
- **`DP`**: Depth of sequencing coverage at the mutation site.
- **`NV`**: Number of reads supporting the variant.
- **`VAF`**: Variant allele frequency, calculated as `NV / DP`.
- **`sample`**: Unique identifier for the sample.

For example, the first few rows of `mutations` look like this:

``` r
tickTack::pcawg_example$mutations %>% head()
#> # A tibble: 6 × 9
#>   chr      from      to ref   alt      DP    NV   VAF sample                    
#>   <chr>   <dbl>   <dbl> <chr> <chr> <dbl> <dbl> <dbl> <chr>                     
#> 1 chr1  1015594 1015594 C     C        99    16 0.162 00db1b95-8ca3-4cc4-bb46-6…
#> 2 chr1  1866371 1866371 C     C       250    67 0.268 00db1b95-8ca3-4cc4-bb46-6…
#> 3 chr1  1921712 1921712 C     C        62    20 0.323 00db1b95-8ca3-4cc4-bb46-6…
#> 4 chr1  2049858 2049858 G     G       118    15 0.127 00db1b95-8ca3-4cc4-bb46-6…
#> 5 chr1  2357842 2357842 C     C        84    12 0.143 00db1b95-8ca3-4cc4-bb46-6…
#> 6 chr1  2771915 2771915 G     G        90    27 0.3   00db1b95-8ca3-4cc4-bb46-6…
```

### 2. Copy Number Alterations (CNA)

The `cna` component is a tibble that describes regions of the genome
with alterations in copy number. Each row represents a genomic segment,
with the following columns:

- **`chr`**: Chromosome of the segment.
- **`from`** and **`to`**: Start and end positions of the segment.
- **`Major`** and **`minor`**: Major and minor allele copy numbers.
- **`CCF`** Cancer cell fraction for the segment.
- **`total_cn`**: Total copy number (sum of **`Major`** and
  **`minor`**).

Here is the preview of the `cna` data:

``` r
tickTack::pcawg_example$cna %>% head()
#> # A tibble: 6 × 7
#>   chr       from        to Major minor   CCF total_cn
#>   <chr>    <dbl>     <dbl> <dbl> <dbl> <dbl>    <dbl>
#> 1 chr1     10001    790008     2     2 1            4
#> 2 chr1    790009  13212499     2     2 1            4
#> 3 chr1  13212500  33458785     2     2 1            4
#> 4 chr1  33458786  33564126     2     2 0.194        4
#> 5 chr1  33564127  56834601     2     2 1            4
#> 6 chr1  56834602 121499999     2     2 1            4
```

### 3. Metadata

The `metadata` component is a tibble containing sample-level
information, with the following columns:

- **`sample`**: Unique identifier for the sample.
- **`purity`**: Tumor purity, representing the proportion of cancer
  cells in the sample.
- **`ploidy`**: Average ploidy of the sample.
- **`purity_conf_mad`**: Confidence interval for the purity estimate.
- **`wgd_status`**: Whole genome doubling status (e.g., **`wgd`** or
  **`no wgd`**).
- **`wgd_uncertain`**: Logical indicating uncertainty in the
  **`wgd_status`**.

An example of the `metadata is shown below`:

``` r
tickTack::pcawg_example$metadata
#> # A tibble: 1 × 6
#>   sample                  purity ploidy purity_conf_mad wgd_status wgd_uncertain
#>   <chr>                    <dbl>  <dbl>           <dbl> <chr>      <lgl>        
#> 1 00db1b95-8ca3-4cc4-bb4…  0.406   4.08           0.009 wgd        FALSE
```

## Types of Copy Number Alterations (CNAs)

The `cna` component can include different types of copy number segments,
categorized as follows:

### Clonal Simple CNAs

These are straightforward alterations with specific copy number
configurations:

- **`1:0`**: Loss of heterozygosity (LOH).
- **`2:0`**: Copy neutral LOH.
- **`1:1`**: Diploid heterozygous (assumed to be the normal reference).
- **`2:1`**: Trisomy.
- **`2:2`**: Tetraploidy.

### Clonal Complex and Subclonal CNAs

Clonal complex CNAs include any clonal CNA that is not considered
simple. On the other hand, subclonal CNAs involve a mixture of two
subclones, where each subclone is defined by a simple CNA.

### Timing Eligibility

Only **clonal CNAs** with specific configurations can be timed:

- **`2:1`**: Trisomy.
- **`2:0`**: Copy neutral LOH.
- **`2:2`**: Tetraploidy.

> Note : Only **clonal CNAs** are accepted by `tickTack`. Subclonal CNAs
> or other complex configurations are not supported for timing analysis.

## Reference genome coordinates

tickTack uses a genome coordinates reference system to convert relative
relative to absolute coordinates, a step required to plot segments
across the whole genome. For instance, if a mutation maps to position
$100$ of chromosome `chr2`, its absolute coordinate is $100 + L$ where
$L$ is the length of `chr1`. The reference system adopted by tickTack
needs therefore to report the length of each chromosome, plus the
information regarding the boundary of each centromere.

> Note: mapping of mutations onto segments is independent of the
> reference genome, and it will work as far as both mutation and CNA
> segments are mapped to the same reference.

tickTack supports two coordinates reference genomes:

- `hg19` or `GRCh37`;
- `hg38` or `GRCh38` (default),

for which two dataframes are stored inside the package.

``` r
tickTack::chr_coordinates_hg19
#> # A tibble: 24 × 6
#>    chr      length       from         to centromerStart centromerEnd
#>    <chr>     <int>      <dbl>      <dbl>          <dbl>        <dbl>
#>  1 chr1  249250621          0  249250621      121535434    124535434
#>  2 chr2  243199373  249250621  492449994      341576792    344576792
#>  3 chr3  198022430  492449994  690472424      582954848    585954848
#>  4 chr4  191154276  690472424  881626700      740132541    743132541
#>  5 chr5  180915260  881626700 1062541960      928032341    931032341
#>  6 chr6  171115067 1062541960 1233657027     1121372126   1124372126
#>  7 chr7  159138663 1233657027 1392795690     1291711358   1294711358
#>  8 chr8  146364022 1392795690 1539159712     1436634577   1439634577
#>  9 chr9  141213431 1539159712 1680373143     1586527391   1589527391
#> 10 chr10 135534747 1680373143 1815907890     1719628078   1722628078
#> # ℹ 14 more rows

tickTack::chr_coordinates_GRCh38
#> # A tibble: 24 × 6
#>    chr      length       from         to centromerStart centromerEnd
#>    <chr>     <dbl>      <dbl>      <dbl>          <dbl>        <dbl>
#>  1 chr1  248956422          0  248956422      122026459    124849229
#>  2 chr2  242193529  248956422  491149951      341144567    341144567
#>  3 chr3  198295559  491149951  689445510      581922409    582703370
#>  4 chr4  190214555  689445510  879660065      739157571    739157571
#>  5 chr5  181538259  879660065 1061198324      926145965    929381368
#>  6 chr6  170805979 1061198324 1232004303     1119752212   1119752212
#>  7 chr7  159345973 1232004303 1391350276     1290173956   1293382091
#>  8 chr8  145138636 1391350276 1536488912     1435384020   1435384020
#>  9 chr9  138394717 1536488912 1674883629     1579878547   1579878547
#> 10 chr10 133797422 1674883629 1808681051     1714570311   1716429449
#> # ℹ 14 more rows
```
