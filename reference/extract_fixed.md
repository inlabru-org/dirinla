# Extracting marginals fixed of an inla object

`extract_fixed` is a function to extract summary and marginals
distribution corresponding to the fixed effects

## Usage

``` r
extract_fixed(inla_model, names_cat)
```

## Arguments

- inla_model:

  Object of inla class.

- names_cat:

  List generated with extract_formula.

## Value

summary_fixed Summary of fixed effects for each category.

marginals_fixed Marginals for each parameter estimated.

## Author

Joaquín Martínez-Minaya <jomarminaya@gmail.com>
