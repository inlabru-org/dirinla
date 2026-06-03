# Block diagonal matrix creation

Fast version of
[`Matrix::.bdiag()`](https://rdrr.io/pkg/Matrix/man/bdiag.html) – for
the case of *many* (k x k) matrices: Copyright (C) 2016 Martin Maechler,
ETH Zurich

## Usage

``` r
bdiag_m(lmat)
```

## Arguments

- lmat:

  `list(<mat1>, <mat2>, ....., <mat_N>)` where each `mat_j` is a `k x k`
  'matrix'

## Value

a sparse (N*k x N*k) matrix of class
[Matrix::dgCMatrix](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html).
