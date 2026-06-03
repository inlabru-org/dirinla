# Formula in to list

`formula_list` reads the formula and generates a list with the name of
the covariates used in each category

## Usage

``` r
formula_list(form, y = NULL)
```

## Arguments

- form:

  Object of class formula.

- y:

  Matrix containing the response variable \\R^{nxd}\\, being n number of
  individuals and d the number of categories.

## Value

A list with the names of the variables used in each category.

## Author

Joaquín Martínez-Minaya <jomarminaya@gmail.com>

## Examples

``` r
formula <- y ~ 1 + v1 + v2 | -1 + v1 | 0 + v2
formula_list(formula)
#> $`category 1`
#> [1] "intercept" "v1"        "v2"       
#> 
#> $`category 2`
#> [1] "v1"
#> 
#> $`category 3`
#> [1] "v2"
#> 
```
