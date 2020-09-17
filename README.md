
<!-- README.md is generated from README.Rmd. Please edit that file -->

# crc.partialid

<!-- badges: start -->

<!-- badges: end -->

The goal of crc.partialid is to implement the partial identification
analysis approach for capture-recapture (CRC) experiments, introduced by
\[Sun, Jinghao, Luk Van Baelen, Els Plettinckx, and Forrest W. Crawford.
“Partial identification and dependence-robust confidence intervals for
capture-recapture surveys.” arXiv preprint arXiv:2008.00127 (2020)\].

## Installation

Currently, you can install the development version of crc.partialid from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Jinghao-Sun/crc.partialid")
```

## Example

This is a basic example which shows you how to estimate the population
size with crc.partialid.

### Import Data

We use 3 capture samples of People Who Inject Drugs (PWID) in Brussels,
Belgium, to estimate the total population size of PWID in Brusells. See
\[Els Plettinckx, Forrest W. Crawford, J ́erˆome Antoine, Gremeaux Lies,
and Luk Van Baelen. How many people injected drugs over the last 12
months in Belgium? Estimates based on a capture-recapture and multiplier
method. Working Paper, 2020.\]

The format of the data needs to be a data.frame or tibble that contains
aggregated capture histories. For a CRC study with k samples, the
data.frame contains one row per observed capture history followed by its
frequency. the data.frame has 2^k - 1 rows and k+1 columns, where the
first k columns can only have 0 or 1 (1 indicates captured in this
sample), and the last column is the frequency.

``` r
library(crc.partialid)
print(mem_list)
#>   Fieldwork_study Low_threshold Crisis_intervention   N
#> 1               1             0                   0  89
#> 2               1             1                   0  24
#> 3               0             1                   0 103
#> 4               1             1                   1  27
#> 5               1             0                   1  29
#> 6               0             1                   1  13
#> 7               0             0                   1  21
```

### Interval estimates for PWID population size

  - Suppose from empirical knowledge, one assumes that the pairwise
    dependence (odds ratio) between samples (1, 2), (1, 3) and (2, 3)
    are all between 1 and 5. Then we can encode the above constraints
    with the following conditions:

<!-- end list -->

``` r
r_cond = c(1, 1, 2)
t_cond = c(2, 3, 3)
eta_cond = c(1, 1, 1)
xi_cond = c(5, 5, 5)
```

We can use ci.tib or ci.pl to compute 95% confidence intervals under
such constraints. For example,

``` r
ci.tib(mem_list, restriction_type = "pairwise", alpha = 0.05, r_cond = r_cond, 
       t_cond = t_cond, eta_cond = eta_cond, xi_cond = xi_cond, B = 1000,
       lb = 0.5, ub = 10, tol = 1e-2, search_step = 30, tasks = 2, verbose = TRUE)

ci.pl(mem_list, restriction_type = "pairwise", alpha = 0.05, r_cond = r_cond, 
       t_cond = t_cond, eta_cond = eta_cond, xi_cond = xi_cond, lb = 0.5,
      ub = 30, tol = 1e-3)
```

  - Sometimes, one can only have knowledge about the dependence between
    some but not all pairs, i.e. odds ratio for sample (1, 2) is between
    1 and 10 and odds ratio for samples (1, 3) is between 0.5 and 6:

<!-- end list -->

``` r
r_cond = c(1, 1)
t_cond = c(2, 3)
eta_cond = c(1, 0.5)
xi_cond = c(10, 6)
```

  - When there isn’t much empirical or domain knowledge, restrictions on
    the highest order interaction term in a loglinear model can be
    useful to check the robustness, which is forced to be 0 in the usual
    CRC models. For example, setting \(gamma = 0.1\), we thus restrict
    \(\lambda_c\) between -0.1 and 0.1.

<!-- end list -->

``` r
ci.tib(mem_list, restriction_type = "highest", gamma = 0.1, alpha = 0.05, B = 500,lb = 0.5, ub = 60, tol = 1e-3, search_step = 30, tasks = 2, verbose = TRUE)

ci.pl(mem_list, restriction_type = "highest", gamma = 0.1, lb = 0.05, ub = 60, tol = 1e-3)
```

## Remarks

  - Varying the values of parameters “lb” and “ub” to avoid being
    trapped in local minima.
  - Test Inversion Bootstrap Confidence Interval (ci.tib) needs a large
    amount of resamples (B needs to be large, like more than 2000) to
    get stable and precise estimates. Thus, we allow parallel computing
    to draw bootstrap resamples by specifying the parameter “tasks”.
