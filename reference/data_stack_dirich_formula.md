# Preparing the data

`data_stack_dirich_formula` prepares the data using inla.stack from the
package INLA.

## Usage

``` r
data_stack_dirich_formula(y, covariates, share = NULL, data, d, n)
```

## Arguments

- y:

  Response variable in a matrix format.

- covariates:

  String with the name of covariates.

- share:

  Covariates to share in all the cateogries. Not implemented yet.

- data:

  Data.frame which contains all the covariates.

- d:

  Number of categories.

- n:

  Number of locations.

## Value

List with two objects

- Object of class inla.stack

- Object with class formula

## Author

Joaquín Martínez-Minaya <jomarminaya@gmail.com>

## Examples

``` r
n <- 100
d <- 4

V <- matrix(rnorm(4*n, 0, 1), ncol=4)
V <- as.data.frame(V)
names(V) <- c('v1', 'v2', 'v3', 'v4' )

covariates <- names(V)


formula <- y ~ 1 + v1 + v2 | 1 + v1 | 1 + v1

names_cat <- formula_list(formula)

data_stack_construct <- data_stack_dirich(y          = as.vector(rep(NA, n*d)),
                                          covariates = names_cat,
                                          share      = NULL,
                                          data       = V,
                                          d          = d,
                                          n          = n )
```
