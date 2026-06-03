# Computing additional diagonal part for the real Hessian H = H0 + diag

`H_matrix_eta_diag` computes the expected Hessian in eta of
-loglikelihood

## Usage

``` r
H_matrix_eta_diag(eta, d, y)
```

## Arguments

- eta:

  eta vector to compute the expected Hessian.

- d:

  Dimension

- y:

  Data corresponding to the i-individual

## Value

Elements of the diagonal such as H = H0 + diag

## Author

Joaquín Martínez-Minaya <jomarminaya@gmail.com>
