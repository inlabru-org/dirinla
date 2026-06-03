# Fitting a Dirichlet regression

`dirinlareg` Main function to do a Dirichlet Regression

## Usage

``` r
dirinlareg(
  formula,
  y,
  data.cov,
  share = NULL,
  x0 = NULL,
  tol0 = 1e-05,
  tol1 = 0.1,
  k0 = 20,
  k1 = 5,
  a = 0.5,
  strategy = "ls-quasi-newton",
  prec = prec,
  verbose = FALSE,
  cores = 1,
  sim = 1000,
  prediction = FALSE,
  data.pred.cov = NULL,
  ...
)
```

## Arguments

- formula:

  object of class formula indicating the response variable and the
  covariates of the Dirichlet regression

- y:

  matrix containing the response variable \\R^{nxd}\\, being n number of
  individuals and d the number of categories

- data.cov:

  data.frame with the covarites, only the covariates!

- share:

  parameters to be fitted jointly.

- x0:

  initial optimization value

- tol0:

  tolerance

- tol1:

  tolerance for the gradient such that \|grad\| \< tol1 \* max(1, \|f\|)

- k0:

  number of iterations

- k1:

  number of iterations including the calling to inla

- a:

  step length in the optimization algorithm

- strategy:

  strategy to use to optimize

- prec:

  precision for the prior of the fixed effects

- verbose:

  if TRUE all the computing process is shown. Default is FALSE

- cores:

  Number of cores for parallel computation. The package parallel is
  used.

- sim:

  Simulations to call inla.posterior.sample and extract linear
  predictor, alphas and mus. The bigger it is, better is the
  approximation, but more computational time.

- prediction:

  if TRUE we will predict with the new values of the covariates given in
  data.pred.cov.

- data.pred.cov:

  data.frame with the values for the covariates where we want to
  predict.

- ...:

  arguments for the inla command

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
y_o <- DirichletReg::rdirichlet(N, alpha)
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
}
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
```
