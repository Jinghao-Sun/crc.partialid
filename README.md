
<!-- README.md is generated from README.Rmd. Please edit that file -->

# crc.partialid

<!-- badges: start -->

<!-- badges: end -->

The goal of R package **crc.partialid** is to implement the partial
identification analysis approach for capture-recapture (CRC)
experiments, introduced by [Sun, Jinghao, Luk Van Baelen, Els
Plettinckx, and Forrest W. Crawford. “Partial identification and
dependence-robust confidence intervals for capture-recapture surveys.”
arXiv preprint arXiv:2008.00127
(2020)](https://arxiv.org/abs/2008.00127).

## Installation

Currently, you can install the development version of **crc.partialid**
from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Jinghao-Sun/crc.partialid")
```

## Example 1

This is a basic example which shows you how to estimate the population
size with **crc.partialid**.

### Import Data

We use 3 capture samples of People Who Inject Drugs (PWID) in Brussels,
Belgium, to estimate the total population size of PWID in Brusells. See
\[Els Plettinckx, Forrest W. Crawford, Jérôme Antoine, Gremeaux Lies,
and Luk Van Baelen. How many people injected drugs over the last 12
months in Belgium? Estimates based on a capture-recapture and multiplier
method. Working Paper, 2020.\]

The format of the data needs to be a data.frame or tibble that contains
aggregated capture histories. For a CRC study with k samples, the
data.frame contains one row per observed capture history followed by its
frequency. The data.frame has 2^k - 1 rows and k+1 columns, where the
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

We can use *ci.tib* or *ci.pl* to compute 95% (i.e., 1 - alpha)
confidence intervals under such constraints. (For the difference of
these two methods, see the [methodology
paper](https://arxiv.org/abs/2008.00127).) For example,

``` r
CITIB = ci.tib(mem_list, restriction_type = "pairwise", alpha = 0.05, r_cond = r_cond, 
       t_cond = t_cond, eta_cond = eta_cond, xi_cond = xi_cond, B = 500,
       lb = 0.5, ub = 10, tol = 1, search_step = 30, tasks = -1, verbose = TRUE)
#> [1] 548.2308
#> [1] 153
#> [1] 3060
#> [1] "start left"
#> [1] 153.0000 548.2308
#> [1] 350.6154
#> [1] 350.6154 548.2308
#> [1] 449.4231
#> [1] 350.6154 449.4231
#> [1] 400.0192
#> [1] 400.0192 449.4231
#> [1] 424.7212
#> [1] 424.7212 449.4231
#> [1] 437.0721
#> [1] 424.7212 437.0721
#> [1] 430.8966
#> [1] 424.7212 430.8966
#> [1] 427.8089
#> [1] 427.8089 430.8966
#> [1] 429.3528
#> [1] 429.3528 430.8966
#> [1] 430.1247
#> [1] "start right"
#> [1]  548.2308 3060.0000
#> [1] 1804.115
#> [1]  548.2308 1804.1154
#> [1] 1176.173
#> [1]  548.2308 1176.1731
#> [1] 862.2019
#> [1] 548.2308 862.2019
#> [1] 705.2163
#> [1] 705.2163 862.2019
#> [1] 783.7091
#> [1] 783.7091 862.2019
#> [1] 822.9555
#> [1] 783.7091 822.9555
#> [1] 803.3323
#> [1] 783.7091 803.3323
#> [1] 793.5207
#> [1] 783.7091 793.5207
#> [1] 788.6149
#> [1] 783.7091 788.6149
#> [1] 786.162
#> [1] 783.7091 786.1620
#> [1] 784.9356
#> [1] 783.7091 784.9356
#> [1] 784.3224
#> CI 429.3528 784.3224

CIPL = ci.pl(mem_list, restriction_type = "pairwise", alpha = 0.05, r_cond = r_cond, 
       t_cond = t_cond, eta_cond = eta_cond, xi_cond = xi_cond, lb = 0.5,
      ub = 10, tol = 1, verbose = TRUE, xtol_rel = 1e-3, maxeval = 300)
#> 
#> Call:
#> nloptr::nloptr(x0 = m0_init, eval_f = eval_l, lb = lb, ub = ub, 
#>     eval_g_eq = NULL, opts = opts)
#> 
#> 
#> Minimization using NLopt version 2.4.2 
#> 
#> NLopt solver status: 4 ( NLOPT_XTOL_REACHED: Optimization stopped because 
#> xtol_rel or xtol_abs (above) was reached. )
#> 
#> Number of Iterations....: 22 
#> Termination conditions:  xtol_rel: 1e-04 maxeval: 100 
#> Number of inequality constraints:  0 
#> Number of equality constraints:    0 
#> Optimal value of objective function:  -931.053268255555 
#> Optimal value of controls: 551.2998
#> 
#> 
#> [1]  352.1499 -919.8478
#> [1]  451.7249 -929.0307
#> [1]  501.5123 -930.5622
#> [1]  476.6186 -929.8793
#> [1]  464.1717 -929.5581
#> [1]  457.9483 -929.3073
#> [1]  454.8366 -929.1717
#> [1]  453.2807 -929.1034
#> [1]  454.0587 -929.1378
#> [1] 1805.6499 -907.5869
#> [1] 1178.4749 -921.2863
#> [1]  864.8873 -927.3271
#> [1]  708.0936 -929.7870
#> [1]  786.4905 -928.6553
#> [1]  747.292 -929.241
#> [1]  766.8912 -928.9486
#> [1]  757.0916 -929.1117
#> [1]  752.1918 -929.0279
#> [1]  749.7419 -929.1969
#> [1]  750.9669 -929.1946
#> [1]  751.5793 -929.1887
#> [1] 453.2807 752.1918
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

CITIB = ci.tib(mem_list, restriction_type = "pairwise", alpha = 0.05, r_cond = r_cond, 
       t_cond = t_cond, eta_cond = eta_cond, xi_cond = xi_cond, B = 500,
       lb = 0.5, ub = 10, tol = 1, search_step = 30, tasks = -1, verbose = TRUE)

CIPL = ci.pl(mem_list, restriction_type = "pairwise", alpha = 0.05, r_cond = r_cond, 
       t_cond = t_cond, eta_cond = eta_cond, xi_cond = xi_cond, lb = 0.5,
      ub = 10, tol = 1, verbose = TRUE, xtol_rel = 1e-3, maxeval = 300)
```

  - When there isn’t much empirical or domain knowledge, restrictions on
    the highest order interaction term in a loglinear model can be
    useful to check the robustness, which is forced to be 0 in the usual
    CRC models. For example, setting \(gamma = 0.1\), we thus restrict
    \(\lambda_c\) between -0.1 and 0.1.

<!-- end list -->

``` r
CITIB = ci.tib(mem_list, restriction_type = "highest", gamma = 0.1, alpha = 0.05, B = 500,lb = 0.5, ub = 60, tol = 1e-3, search_step = 30, tasks = 2, verbose = TRUE)

CIPL = ci.pl(mem_list, restriction_type = "highest", gamma = 0.1, lb = 0.05, ub = 60, tol = 1e-3, xtol_rel = 1e-3, maxeval = 300, verbose = TRUE)
```

## Example 2

## Remarks

  - Varying the values of parameters “lb” and “ub” to avoid being
    trapped in local minima.
  - Varying optimization parameters “tol”, “xtol\_rel”, and “maxeval” to
    improve precison and stability of estimation.
  - Test Inversion Bootstrap Confidence Interval (ci.tib) needs a large
    amount of resamples (B needs to be large, like more than 2000) to
    get stable and precise estimates. Thus, we allow parallel computing
    to draw bootstrap resamples by specifying the parameter “tasks”.

## TODO

  - Test more instances for k = 2, 3, 4 and with different number of
    pairwise restrictions.
