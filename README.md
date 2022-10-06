
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dirinla

<!-- badges: start -->

[![R-CMD-check](https://github.com/inlabru-org/dirinla/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/inlabru-org/dirinla/actions/workflows/R-CMD-check.yaml)
[![dirinla status
badge](https://inlabru-org.r-universe.dev/badges/dirinla)](https://inlabru-org.r-universe.dev)
<!-- badges: end -->

The goal of **dirinla** is to analyze compositional data with a
Dirichlet regression using the integrated nested Laplace approximation
via the [R-INLA package](https://www.r-inla.org/).

Package documentation can be found at
<https://inlabru-org.github.io/dirinla/>

## Installation

It is not yet in CRAN, but you can install the latest bugfix release of
dirinla from [github](https://github.com/inlabru-org/dirinla) with:

``` r
# install.packages("remotes")
remotes::install_github("inlabru-org/dirinla", ref = "master")
```

The latest development version can be installed via
[inlabru-org.r-universe.dev](https://inlabru-org.r-universe.dev/ui#builds):

``` r
# Enable universe(s) by inlabru-org
options(repos = c(
  inlabruorg = 'https://inlabru-org.r-universe.dev',
  INLA = 'https://inla.r-inla-download.org/R/testing',
  CRAN = 'https://cloud.r-project.org'))

# Install the package
install.packages('dirinla')
```

or directly from github:

``` r
# install.packages("remotes")
remotes::install_github("inlabru-org/dirinla", ref = "devel")
```

## Example

This is a basic example which shows you how to solve a common problem:

### Loading libraries

``` r
library(dirinla)
library(INLA)
library(DirichletReg)
```

### Simulating from a Dirichlet likelihood

``` r
set.seed(1000)
N <- 50 #number of data
V <- as.data.frame(matrix(runif((4) * N, 0, 1), ncol = 4)) #Covariates
names(V) <- paste0('v', 1:4)

formula <- y ~ 1 + v1 | 1 + v2 | 1 + v3 | 1 + v4
(names_cat <- formula_list(formula))
#> $`category 1`
#> [1] "intercept" "v1"       
#> 
#> $`category 2`
#> [1] "intercept" "v2"       
#> 
#> $`category 3`
#> [1] "intercept" "v3"       
#> 
#> $`category 4`
#> [1] "intercept" "v4"

x <- c(-1.5, 1, -3, 1.5,
       2, -3 , -1, 5)

mus <- exp(x) / sum(exp(x))
C <- length(names_cat)
A_construct <-
  data_stack_dirich(y = as.vector(rep(NA, N * C)),
                    covariates = names_cat,
                    data       = V,
                    d          = C,
                    n          = N)

A_construct[1:8, ]
#> 8 x 8 sparse Matrix of class "dgCMatrix"
#>                                                     
#> [1,] 1 0.3278787 . .         . .         . .        
#> [2,] . .         1 0.7267993 . .         . .        
#> [3,] . .         . .         1 0.5993679 . .        
#> [4,] . .         . .         . .         1 0.3224284
#> [5,] 1 0.7588465 . .         . .         . .        
#> [6,] . .         1 0.6820559 . .         . .        
#> [7,] . .         . .         1 0.4516818 . .        
#> [8,] . .         . .         . .         1 0.5613199

eta <- A_construct %*% x
alpha <- exp(eta)
alpha <- matrix(alpha,
                ncol  = C,
                byrow = TRUE)
y_o <- rdirichlet(N, alpha)
colnames(y_o) <- paste0("y", 1:C)
head(y_o)
#>                y1           y2          y3        y4
#> [1,] 1.139109e-04 2.413110e-01 0.345644051 0.4129310
#> [2,] 7.342592e-02 3.633687e-03 0.193128806 0.7298116
#> [3,] 2.079875e-02 4.007038e-05 0.323599848 0.6555613
#> [4,] 1.361278e-05 3.894290e-02 0.172308411 0.7887351
#> [5,] 2.339991e-02 4.035880e-03 0.055003377 0.9175608
#> [6,] 8.910796e-01 6.743970e-10 0.000663061 0.1082573
```

### Fitting a simple model

``` r
y <- y_o
model.inla <- dirinlareg(
  formula  = y ~ 1 + v1 | 1 + v2 | 1 + v3 | 1 + v4,
  y        = y,
  data.cov = V,
  prec     = 0.0001,
  verbose  = TRUE)
#> 
#>  
#>  ---------------------- Looking for the mode ----------------- 
#>  
#>  Iter = 1, |grad| = 824.56, log.post = -457.73, |x_new - x_old| = 13.92583, |f_new - f_old| = 368.13749
#> Iter = 2, |grad| = 136.74, log.post = -577.55, |x_new - x_old| = 3.65634, |f_new - f_old| = 119.82746
#> Iter = 3, |grad| = 120.98, log.post = -671.77, |x_new - x_old| = 3.31293, |f_new - f_old| = 94.22001
#> Iter = 4, |grad| = 96.8, log.post = -739.43, |x_new - x_old| = 2.81732, |f_new - f_old| = 67.65798
#> Iter = 5, |grad| = 70.76, log.post = -785.42, |x_new - x_old| = 2.1453, |f_new - f_old| = 45.99234
#> Iter = 6, |grad| = 48.13, log.post = -809.29, |x_new - x_old| = 1.3377, |f_new - f_old| = 23.86892
#> Iter = 7, |grad| = 24.92, log.post = -814.96, |x_new - x_old| = 0.55663, |f_new - f_old| = 5.67249
#> Iter = 8, |grad| = 5.97, log.post = -815.24, |x_new - x_old| = 0.10262, |f_new - f_old| = 0.27719
#> Iter = 9, |grad| = 0.3, log.post = -815.24, |x_new - x_old| = 0.00994, |f_new - f_old| = 0.00062
#> Iter = 10, |grad| = 0, log.post = -815.24, |x_new - x_old| = 0.00011, |f_new - f_old| = 0
#> Iter = 11, |grad| = 0, log.post = -815.24, |x_new - x_old| = 0, |f_new - f_old| = 0
#> 
#> Great news! The mode has been properly located!
#>  
#> Real Hessian has been used 10 times 
#> Expected Hessian has been used 40 times 
#> 
#>  ----------------------    INLA call    ----------------- 
#> INLA-Iter = 1, fixed.effects = 8, hyperparameters = 0 ---> PASS
#> 
#>  ---------------------- Obtaining linear predictor -----------------

summary(model.inla)
#> 
#> Call: 
#>  dirinlareg(formula = y ~ 1 + v1 | 1 + v2 | 1 + v3 | 1 + v4, y = y, 
#>     data.cov = V, prec = 1e-04, verbose = TRUE)
#> 
#>  
#> ---- FIXED EFFECTS ---- 
#> ======================================================================= 
#> y1
#> ----------------------------------------------------------------------- 
#>              mean     sd 0.025quant 0.5quant 0.975quant    mode
#> intercept -1.5293 0.2922    -2.1020  -1.5293    -0.9565 -1.5293
#> v1         0.9977 0.5146    -0.0108   0.9977     2.0063  0.9977
#> ======================================================================= 
#> y2
#> ----------------------------------------------------------------------- 
#>              mean     sd 0.025quant 0.5quant 0.975quant    mode
#> intercept -2.8538 0.2783    -3.3993  -2.8538     -2.308 -2.8538
#> v2         0.7187 0.4870    -0.2357   0.7187      1.673  0.7187
#> ======================================================================= 
#> y3
#> ----------------------------------------------------------------------- 
#>             mean     sd 0.025quant 0.5quant 0.975quant   mode
#> intercept  1.902 0.2448      1.422    1.902      2.382  1.902
#> v3        -3.045 0.3738     -3.777   -3.045     -2.312 -3.045
#> ======================================================================= 
#> y4
#> ----------------------------------------------------------------------- 
#>              mean     sd 0.025quant 0.5quant 0.975quant    mode
#> intercept -0.7033 0.3145     -1.320  -0.7033   -0.08698 -0.7033
#> v4         4.5247 0.4313      3.679   4.5247    5.37012  4.5247
#> ======================================================================= 
#> 
#> ---- HYPERPARAMETERS ---- 
#> 
#>  No hyperparameters in the model 
#> ======================================================================= 
#> DIC = 1555.1536 , WAIC = 1039.5914 , LCPO = 779.5894 
#> Number of observations: 50
#> Number of Categories: 4
```

### Predicting for v1 = 0.25, v2 = 0.5, v3 = 0.5, v4 = 0.1

``` r
model.prediction <-
  predict(model.inla,
                  data.pred = data.frame(v1 = 0.25,
                                         v2 = 0.5,
                                         v3 = 0.5,
                                         v4 = 0.1))
#> 
#>  
#>  ---------------------- Predicting ----------------- 
#>  
#> 
model.prediction$summary_predictive_means
#> $y1
#>            Min.    1st Qu.    Median     Mean   3rd Qu.     Max.
#> [1,] 0.02366363 0.08237956 0.1045064 0.110114 0.1321715 0.306817
#> 
#> $y2
#>             Min.    1st Qu.     Median       Mean    3rd Qu.      Max.
#> [1,] 0.005978446 0.02318953 0.03082928 0.03341511 0.04057716 0.1574298
#> 
#> $y3
#>           Min.   1st Qu.    Median      Mean   3rd Qu.    Max.
#> [1,] 0.2386448 0.4920715 0.5580837 0.5553907 0.6213676 0.86177
#> 
#> $y4
#>            Min.   1st Qu.    Median      Mean   3rd Qu.     Max.
#> [1,] 0.08145196 0.2404386 0.2950362 0.3010802 0.3541525 0.626968
```
