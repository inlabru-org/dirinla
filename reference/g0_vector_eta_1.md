# Computing gradient vector in eta

`g0_vector_eta` computes the gradient of -loglikelihood

## Usage

``` r
g0_vector_eta_1(A = A, x, y)
```

## Arguments

- A:

  Matrix which links eta with the latent field, i.e., eta = A x.

- x:

  Vector with the elements of the latent field, i.e., eta = A x.

- y:

  Vector with the response variable.

## Value

A numeric vector with the gradient in eta.

## Author

Joaquín Martínez-Minaya <jomarminaya@gmail.com>
