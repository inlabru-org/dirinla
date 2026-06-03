# Dirichlet log posterior function

`dirichlet_log_pos_x` returns the -log posterior Dirichlet distribution
asumming multivariate normal prior with precision matrix Qx for elements
of the latent field.

## Usage

``` r
dirichlet_log_pos_x(A = A, x, Qx = Qx, y)
```

## Arguments

- A:

  A matrix which links eta with the latent field, i.e., eta = A x.

- x:

  Vector with the elements of the latent field, i.e., eta = A x.

- Qx:

  Precision matrix for the priors of the latent field.

- y:

  Vector with the response variable.

## Value

A real value showing the -log posterior density is returned

## Author

Joaquín Martínez-Minaya <jomarminaya@gmail.com>
