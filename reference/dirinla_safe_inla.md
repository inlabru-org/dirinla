# Load INLA safely for examples and tests

Loads the INLA package with
[`requireNamespace("INLA", quietly = TRUE)`](https://rdrr.io/r/base/ns-load.html),
and optionally checks and sets the multicore `num.threads` INLA option.

## Usage

``` r
dirinla_safe_inla(multicore = NULL, quietly = FALSE)
```

## Arguments

- multicore:

  logical; if `TRUE`, multiple cores are allowed, and the INLA
  `num.threads` option is not checked or altered. If `FALSE`, forces
  `num.threads="1:1"`. Default: NULL, checks if running in testthat or
  non-interactively, in which case sets `multicore=FALSE`, otherwise
  `TRUE`.

- quietly:

  logical; if `TRUE`, prints diagnostic messages. Default: FALSE.

## Value

logical; `TRUE` if INLA was loaded safely, otherwise FALSE

## Details

Code adapted from `inlabru::bru_safe_inla`

## Author

Finn Lindgren

## Examples

``` r
if (FALSE) { # \dontrun{
if (dirinla_safe_inla()) {
  # Run inla dependent calculations
}
} # }
```
