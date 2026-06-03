# Finding the mode of the full posterior distribution

`predict.dirinlaregmodel` computes the posterior predictive distribution
for some given values of the covariates

## Usage

``` r
# S3 method for class 'dirinlaregmodel'
predict(object, data.pred.cov, ...)
```

## Arguments

- object:

  dirinlaregmodel object.

- data.pred.cov:

  Data.frame with the covariate values for the variables to predict.

- ...:

  Other arguments.

## Value

model dirinlaregmodel object

## Author

Joaquín Martínez-Minaya <jomarminaya@gmail.com>

## Examples

``` r
if (dirinla_safe_inla() &&
    requireNamespace("DirichletReg", quietly = TRUE)) {
### In this example, we show how to fit a model using the dirinla package ###
### --- 1. Loading the libraries --- ####
library(INLA)
library(DirichletReg)


### --- 2. Simulating from a Dirichlet likelihood --- ####
set.seed(1000)
N <- 50 #number of data
V <- as.data.frame(matrix(runif((4) * N, 0, 1), ncol = 4)) #Covariates
names(V) <- paste0('v', 1:4)

formula <- y ~ 1 + v1 | 1 + v2 | 1 + v3 | 1 + v4
(names_cat <- formula_list(formula))

x <- c(-1.5, 1, -3, 1.5,
       2, -3 , -1, 5)

mus <- exp(x) / sum(exp(x))
C <- length(names_cat)
data_stack_construct <-
  data_stack_dirich(y = as.vector(rep(NA, N * C)),
                    covariates = names_cat,
                    data       = V,
                    d          = C,
                    n          = N)

A_construct <- data_stack_construct
A_construct[1:8, ]

eta <- A_construct %*% x
alpha <- exp(eta)
alpha <- matrix(alpha,
                ncol  = C,
                byrow = TRUE)
y_o <- rdirichlet(N, alpha)
colnames(y_o) <- paste0("y", 1:C)
head(y_o)


### --- 3. Fitting the model --- ####
y <- y_o
model.inla <- dirinlareg(
  formula  = y ~ 1 + v1 | 1 + v2 | 1 + v3 | 1 + v4,
  y        = y,
  data.cov = V,
  prec     = 0.0001,
  verbose  = FALSE)


summary(model.inla)
### --- 4. Predicting for v1 = 0.25, v2 = 0.5, v3 = 0.5, v4 = 0.1 --- ####
model.prediction <- predict(model.inla,
                data.pred.cov= data.frame(v1 = 0.25,
                                       v2 = 0.5,
                                       v3 = 0.5,
                                       v4 = 0.1))
model.prediction$summary_predictive_means
}
#> Loading required package: Matrix
#> 
#> Loading required package: Formula
#> 
#>  
#>  ---------------------- Looking for the mode ----------------- 
#>  
#>  
#>  ----------------------    INLA call    ----------------- 
#> 
#>  ---------------------- Obtaining linear predictor ----------------- 
#> 
#> Call: 
#>  dirinlareg(formula = y ~ 1 + v1 | 1 + v2 | 1 + v3 | 1 + v4, y = y, 
#>     data.cov = V, prec = 1e-04, verbose = FALSE)
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
#> DIC = 1555.1541 , WAIC = 1039.5909 , LCPO = 779.5896 
#> Number of observations: 50
#> Number of Categories: 4
#> 
#>  
#>  ---------------------- Predicting ----------------- 
#>  
#>  
#> $y1
#>            Min.    1st Qu.    Median      Mean 3rd Qu.      Max.
#> [1,] 0.02366365 0.08238598 0.1044688 0.1100701 0.13194 0.3068168
#> 
#> $y2
#>             Min.    1st Qu.    Median       Mean   3rd Qu.      Max.
#> [1,] 0.005978457 0.02321143 0.0308166 0.03341568 0.0405673 0.1574295
#> 
#> $y3
#>          Min.   1st Qu.   Median      Mean   3rd Qu.      Max.
#> [1,] 0.238645 0.4918204 0.558138 0.5552991 0.6211048 0.8617697
#> 
#> $y4
#>            Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#> [1,] 0.08145221 0.2405209 0.2952006 0.3012151 0.3543549 0.6269677
#> 
```
