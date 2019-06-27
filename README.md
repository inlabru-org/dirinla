
dirinla
=======

The goal of **dirinla** is to analyze compositional data with a Dirichlet regression using the integrated nested Laplace approximation via the [R-INLA package](http://www.r-inla.org).

Installation
------------

It is not yet in CRAN, but you can install the latest bugfix release of dirinla from [bitbucket](https://bitbucket.org/) with:

``` r
# install.packages("remotes")
remotes::install_bitbucket("joaquin-martinez-minaya/dirinla", ref="master")
```

You can install the development version of dirinla from [bitbucket](https://bitbucket.org/) with:

``` r
# install.packages("remotes")
remotes::install_bitbucket("joaquin-martinez-minaya/dirinla", ref="devel")
```

Example
-------

This is a basic example which shows you how to solve a common problem:

### Loading libraries

``` r
library(dirinla)
library(INLA)
library(DirichletReg)
library(ggplot2)
library(gridExtra)
```

### Simulating from a Dirichlet likelihood

``` r
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

A_construct <- data_stack_construct$A
A_construct[1:8, ]

eta <- A_construct %*% x
alpha <- exp(eta)
alpha <- matrix(alpha,
                ncol  = C,
                byrow = TRUE)
y_o <- rdirichlet(N, alpha)
colnames(y_o) <- paste0("y", 1:C)
head(y_o)
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

summary(model.inla)
```

### Predicting for v1 = 0.25, v2 = 0.5, v3 = 0.5, v4 = 0.1

``` r
model.prediction <-
  predict_dirinla(model.inla,
                  data.pred = data.frame(v1 = 0.25,
                                         v2 = 0.5,
                                         v3 = 0.5,
                                         v4 = 0.1))
```
