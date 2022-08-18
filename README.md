
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
mle <- yppe(Surv(time, status)~trt, data=gastric, 
            approach = "mle", init = 0, n_int = 10)
summary(mle)
#> Call:
#> yppe(formula = Surv(time, status) ~ trt, data = gastric, n_int = 10, 
#>     approach = "mle", init = 0)
#> 
#> Short-term coefficients:
#>     Estimate  StdErr z.value  p.value   
#> trt  1.77113 0.61843  2.8639 0.004185 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Long-term coefficients:
#>     Estimate   StdErr z.value   p.value    
#> trt -0.98230  0.29576 -3.3213 0.0008961 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> --- 
#> loglik = 6.751216   AIC = 10.49757

# Bayesian approach:
bayes <- yppe(Surv(time, status)~trt, data=gastric, 
              approach = "bayes", n_int = 10, 
              refresh = FALSE)

summary(bayes)
#> Call:
#> yppe(formula = Surv(time, status) ~ trt, data = gastric, n_int = 10, 
#>     approach = "bayes", refresh = FALSE)
#> 
#> Short-term coefficients:
#>      mean se_mean    sd  2.5%  25%   50%   75% 97.5%    n_eff  Rhat
#> trt 1.817   0.014 0.628 0.681 1.38 1.776 2.215 3.093 1918.249 1.001
#> 
#> Long-term coefficients:
#>       mean se_mean    sd   2.5%    25%    50%    75%  97.5%    n_eff Rhat
#> trt -0.965   0.006 0.305 -1.538 -1.177 -0.965 -0.766 -0.336 2491.661    1
#> 
#> --- 
#> Inference for Stan model: yppe.
#> 4 chains, each with iter=2000; warmup=1000; thin=1; 
#> post-warmup draws per chain=1000, total post-warmup draws=4000.
```
