# Computing expected Hessian in eta

`H0_matrix_eta_x` computes the expected Hessian in eta of -loglikelihood

## Usage

``` r
H0_matrix_eta_x(eta, d, cores)
```

## Arguments

- eta:

  Linear predictor resulting of the product \\A x\\.

- d:

  Dimension.

- cores:

  Number of cores for parallel computation. The package parallel is
  used.

## Value

Expected Hessian in eta.

## Author

Joaquín Martínez-Minaya <jomarminaya@gmail.com>
