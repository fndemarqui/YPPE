
<!-- README.md is generated from README.Rmd. Please edit that file -->

# YPPE

<!-- badges: start -->
<!-- badges: end -->

The R package YPPE provides semiparametric modeling of lifetime data
with crossing survival curves via Yang and Prentice model with piecewise
exponential baseline distribution. Details about the model can be found
in Demarqui and Mayrink (2019) \<doi.org/10.1214/20-BJPS471\>. Model
fitting carried out via likelihood-based and Bayesian approaches. The
package also provides point and interval estimation for the crossing
survival times.

## Installation

You can install the released version of YPPE from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("YPPE")
```

And the development version from [GitHub](https://github.com/) with:

``` r
install.packages("remotes")
remotes::install_github("fndemarqui/YPPE")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(YPPE)
#> Loading required package: survival

data(gastric)

# MLE approach:
mle <- yppe(Surv(time, status)~trt, data=gastric, approach = "mle", init = 0)
summary(mle)
#> Call:
#> yppe(formula = Surv(time, status) ~ trt, data = gastric, approach = "mle", 
#>     init = 0)
#> 
#> Short-term coefficients:
#>     Estimate StdErr z.value  p.value   
#> trt   1.8373 0.6483   2.834 0.004597 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Long-term coefficients:
#>     Estimate   StdErr z.value   p.value    
#> trt -1.01753  0.30036 -3.3877 0.0007049 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> --- 
#> loglik = 47.65404   AIC = 62.69192

# Bayesian approach:
bayes <- yppe(Surv(time, status)~trt, data=gastric, approach = "bayes", verbose = FALSE)
#> 
#> SAMPLING FOR MODEL 'yppe' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.001536 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 15.36 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 27.8873 seconds (Warm-up)
#> Chain 1:                25.3193 seconds (Sampling)
#> Chain 1:                53.2066 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'yppe' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.001448 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 14.48 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 27.8573 seconds (Warm-up)
#> Chain 2:                24.9122 seconds (Sampling)
#> Chain 2:                52.7695 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'yppe' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 0.001441 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 14.41 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 26.8641 seconds (Warm-up)
#> Chain 3:                25.2976 seconds (Sampling)
#> Chain 3:                52.1617 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'yppe' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 0.001508 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 15.08 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 27.2158 seconds (Warm-up)
#> Chain 4:                25.2288 seconds (Sampling)
#> Chain 4:                52.4446 seconds (Total)
#> Chain 4:
summary(bayes)
#> Call:
#> yppe(formula = Surv(time, status) ~ trt, data = gastric, approach = "bayes", 
#>     verbose = FALSE)
#> 
#> Short-term coefficients:
#>      mean se_mean    sd 2.5%  25%   50%   75% 97.5%    n_eff  Rhat
#> trt 1.843   0.014 0.612  0.7 1.42 1.827 2.226 3.144 2019.152 1.001
#> 
#> Long-term coefficients:
#>       mean se_mean   sd   2.5%   25%    50%    75%  97.5%  n_eff Rhat
#> trt -0.938   0.005 0.32 -1.526 -1.16 -0.955 -0.733 -0.277 3545.9    1
#> 
#> --- 
#> Inference for Stan model: yppe.
#> 4 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=4000.
```
