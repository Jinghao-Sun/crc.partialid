
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

## Example: People Who Inject Drugs (PWID) in Brussels, Belgium

We use 3 capture samples of People Who Inject Drugs (PWID) in Brussels,
Belgium, to estimate the total population size of PWID in Brusells. See
\[Els Plettinckx, Forrest W. Crawford, Jérôme Antoine, Gremeaux Lies,
and Luk Van Baelen. How many people injected drugs over the last 12
months in Belgium? Estimates based on a capture-recapture and multiplier
method. Working Paper, 2020.\]

### Import Data

The format of the data needs to be a data.frame or tibble that contains
aggregated capture histories. For a CRC study with k samples, the
data.frame contains one row per observed capture history followed by its
frequency. The data.frame has 2^k - 1 rows and k+1 columns, where the
first k columns can only have 0 or 1 (1 indicates captured in this
sample), and the last column is the frequency.

``` r
library(crc.partialid)
print(PWID)
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
CITIB = ci.tib(PWID, restriction_type = "pairwise", alpha = 0.05, r_cond = r_cond, 
       t_cond = t_cond, eta_cond = eta_cond, xi_cond = xi_cond, B = 500,
       lb = 0.5, ub = 10, tol = 1, search_step = 30, tasks = -1, verbose = TRUE)
#> [1] 548.2308
#> [1] 153
#> [1] 3060
#> [1] "Finding left endpoint of confidence interval using bisection..."
#> [1] 350.6154
#> [1] 449.4231
#> [1] 400.0192
#> [1] 424.7212
#> [1] 437.0721
#> [1] 430.8966
#> [1] 433.9844
#> [1] 435.5282
#> [1] 434.7563
#> [1] "Finding right endpoint of confidence interval using bisection..."
#> [1] 1804.115
#> [1] 1176.173
#> [1] 862.2019
#> [1] 705.2163
#> [1] 783.7091
#> [1] 822.9555
#> [1] 803.3323
#> [1] 793.5207
#> [1] 788.6149
#> [1] 786.162
#> [1] 784.9356
#> [1] 784.3224
#> 1 - 0.05  confidence interval: 433.9844 784.3224

CIPL = ci.pl(PWID, restriction_type = "pairwise", alpha = 0.05, r_cond = r_cond, 
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
#> Number of Iterations....: 24 
#> Termination conditions:  xtol_rel: 1e-04 maxeval: 100 
#> Number of inequality constraints:  0 
#> Number of equality constraints:    0 
#> Optimal value of objective function:  -931.004106554509 
#> Optimal value of controls: 573.0189
#> 
#> 
#> [1] "Finding left endpoint of confidence interval using bisection..."
#> Current Value, Current, Objective:  363.0094 -921.4101 
#> Current Value, Current, Objective:  468.0141 -929.6612 
#> Current Value, Current, Objective:  415.5118 -926.807 
#> Current Value, Current, Objective:  441.763 -928.5132 
#> Current Value, Current, Objective:  454.8886 -929.1556 
#> Current Value, Current, Objective:  448.3258 -928.8484 
#> Current Value, Current, Objective:  451.6072 -929.0097 
#> Current Value, Current, Objective:  453.2479 -929.0642 
#> Current Value, Current, Objective:  454.0682 -929.1261 
#> [1] "Finding right endpoint of confidence interval using bisection..."
#> Current Value, Current, Objective:  1816.509 -906.6307 
#> Current Value, Current, Objective:  1194.764 -919.0036 
#> Current Value, Current, Objective:  883.8915 -926.762 
#> Current Value, Current, Objective:  728.4552 -929.3221 
#> Current Value, Current, Objective:  806.1733 -928.2319 
#> Current Value, Current, Objective:  767.3143 -928.5034 
#> Current Value, Current, Objective:  747.8847 -929.1337 
#> Current Value, Current, Objective:  757.5995 -928.728 
#> Current Value, Current, Objective:  752.7421 -928.8461 
#> Current Value, Current, Objective:  750.3134 -928.9272 
#> Current Value, Current, Objective:  749.0991 -928.9087 
#> Current Value, Current, Objective:  748.4919 -928.9231 
#> 1- 0.05  confidence interval: 453.2479 748.4919
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

CITIB = ci.tib(PWID, restriction_type = "pairwise", alpha = 0.05, r_cond = r_cond, 
       t_cond = t_cond, eta_cond = eta_cond, xi_cond = xi_cond, B = 500,
       lb = 0.5, ub = 10, tol = 1, search_step = 30, tasks = -1, verbose = TRUE)

CIPL = ci.pl(PWID, restriction_type = "pairwise", alpha = 0.05, r_cond = r_cond, 
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
CITIB = ci.tib(PWID, restriction_type = "highest", gamma = 0.1, alpha = 0.05, B = 500,lb = 0.5, ub = 60, tol = 1e-3, search_step = 30, tasks = 2, verbose = TRUE)

CIPL = ci.pl(PWID, restriction_type = "highest", gamma = 0.1, lb = 0.05, ub = 60, tol = 1e-3, xtol_rel = 1e-3, maxeval = 300, verbose = TRUE)
```

## Example: People with a Polish nationality living in the Netherlands in 2009.

Two registers of Statistics Netherlands, the GBA and the HKS on people
with a Polish nationality living in the Netherlands in 2009. See [Van
der Heijden, P.G.M., M.J.L.F. Cruy, and G. van Gils. 2011. Aantallen
Geregistreerde en Nietgeregistreerde Burgers uit MOE-landen die in
Nederland Verblijven. Rapportage Schattingen 2008 en 2009. The Number of
Registered and Non-registered Citizens from MOE-countries Residing in
the Netherlands. Reporting Estimations 2008 and 2009. The Hague:
Ministry of Social Affairs and
Employment.](https://www.rijksoverheid.nl/documenten/rapporten/2013/01/14/aantallen-geregistreerde-en-niet-geregistreerde-burgers-uit-moe-landen-die-in-nederland-verblijven)

``` r
print(Polish)
#>   GBA HKS     N
#> 1   1   0 39488
#> 2   0   1  1445
#> 3   1   1   374
```

Considering the Odds Ratio between 2 samples are between 0.5 and 2. Then
the parameter gamma is log(2) in function *ci.tib* and *ci.pl*.

``` r
# CITIB = ci.tib(Polish, restriction_type = "highest", gamma = 0.1, alpha = 0.05, 
#               B = 500,lb = 0.5, ub = 60, tol = 1e-3, search_step = 30, 
#               tasks = 2, verbose = TRUE)

CIPL = ci.pl(Polish, restriction_type = "highest", gamma = log(2), lb = 0.05, 
             ub = 60, tol = 1e-3, xtol_rel = 1e-3, maxeval = 300, verbose = TRUE)
#> [1] "Infeasible!"
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
#> Number of Iterations....: 12 
#> Termination conditions:  xtol_rel: 1e-04 maxeval: 100 
#> Number of inequality constraints:  0 
#> Number of equality constraints:    0 
#> Optimal value of objective function:  -389353.499227705 
#> Optimal value of controls: 127981.7
#> 
#> 
#> [1] "Finding left endpoint of confidence interval using bisection..."
#> Current Value, Current, Objective:  65023.52 -389084 
#> Current Value, Current, Objective:  96502.6 -389336.5 
#> Current Value, Current, Objective:  112242.1 -389352.6 
#> Current Value, Current, Objective:  104372.4 -389347.9 
#> Current Value, Current, Objective:  108307.3 -389350.6 
#> Current Value, Current, Objective:  110274.7 -389351.9 
#> Current Value, Current, Objective:  109291 -389351.3 
#> Current Value, Current, Objective:  109782.8 -389351.6 
#> Current Value, Current, Objective:  109536.9 -389351.6 
#> Current Value, Current, Objective:  109659.9 -389351.6 
#> Current Value, Current, Objective:  109598.4 -389351.6 
#> Current Value, Current, Objective:  109567.7 -389351.6 
#> Current Value, Current, Objective:  109552.3 -389351.6 
#> Current Value, Current, Objective:  109560 -389351.6 
#> Current Value, Current, Objective:  109556.1 -389351.6 
#> Current Value, Current, Objective:  109558 -389351.6 
#> Current Value, Current, Objective:  109557.1 -389351.6 
#> Current Value, Current, Objective:  109556.6 -389351.6 
#> Current Value, Current, Objective:  109556.4 -389351.6 
#> Current Value, Current, Objective:  109556.5 -389351.6 
#> Current Value, Current, Objective:  109556.5 -389351.6 
#> Current Value, Current, Objective:  109556.6 -389351.6 
#> Current Value, Current, Objective:  109556.6 -389351.6 
#> Current Value, Current, Objective:  109556.6 -389351.6 
#> Current Value, Current, Objective:  109556.5 -389351.6 
#> Current Value, Current, Objective:  109556.6 -389351.6 
#> Current Value, Current, Objective:  109556.6 -389351.6 
#> [1] "Finding right endpoint of confidence interval using bisection..."
#> Current Value, Current, Objective:  1303201 -388261.8 
#> Current Value, Current, Objective:  715591.3 -388653.6 
#> Current Value, Current, Objective:  421786.5 -389346.4 
#> Current Value, Current, Objective:  274884.1 -389353.5 
#> Current Value, Current, Objective:  348335.3 -389353.4 
#> Current Value, Current, Objective:  385060.9 -389351.4 
#> Current Value, Current, Objective:  366698.1 -389352.9 
#> Current Value, Current, Objective:  375879.5 -389352.2 
#> Current Value, Current, Objective:  380470.2 -389351.8 
#> Current Value, Current, Objective:  382765.5 -389351.6 
#> Current Value, Current, Objective:  383913.2 -389351.5 
#> Current Value, Current, Objective:  383339.4 -389351.6 
#> Current Value, Current, Objective:  383052.4 -389351.6 
#> Current Value, Current, Objective:  383195.9 -389351.6 
#> Current Value, Current, Objective:  383267.6 -389351.6 
#> Current Value, Current, Objective:  383303.5 -389351.6 
#> Current Value, Current, Objective:  383321.4 -389351.6 
#> Current Value, Current, Objective:  383312.5 -389351.6 
#> Current Value, Current, Objective:  383317 -389351.6 
#> Current Value, Current, Objective:  383319.2 -389351.6 
#> Current Value, Current, Objective:  383318.1 -389351.6 
#> Current Value, Current, Objective:  383318.6 -389351.6 
#> Current Value, Current, Objective:  383318.4 -389351.6 
#> Current Value, Current, Objective:  383318.5 -389351.6 
#> Current Value, Current, Objective:  383318.4 -389351.6 
#> Current Value, Current, Objective:  383318.5 -389351.6 
#> Current Value, Current, Objective:  383318.5 -389351.6 
#> Current Value, Current, Objective:  383318.5 -389351.6 
#> Current Value, Current, Objective:  383318.5 -389351.6 
#> Current Value, Current, Objective:  383318.5 -389351.6 
#> Current Value, Current, Objective:  383318.5 -389351.6 
#> Current Value, Current, Objective:  383318.5 -389351.6 
#> 1- 0.05  confidence interval: 109556.6 383318.5
```

## Remarks

  - Varying the values of parameters “lb” and “ub” to avoid being
    trapped in local minima.
  - Varying optimization parameters “tol”, “xtol\_rel”, and “maxeval” to
    improve precison and stability of estimation.
  - Test Inversion Bootstrap Confidence Interval (ci.tib) needs a large
    amount of resamples (B needs to be large, like more than 2000) to
    get stable and precise estimates. Thus, we allow parallel computing
    to draw bootstrap resamples by specifying the parameter “tasks”.
