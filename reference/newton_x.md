# Newton-Raphson algorithm

`newton_x` computes optimization algorithms to find the mode of the
posterior. Line search strategy with Armijo conditions is implemented

## Usage

``` r
newton_x(A, x_hat, gk, Hk, a, Qx, strategy, y, d = d)
```

## Arguments

- A:

  Matrix which links eta with the latent field, i.e., eta = A x

- x_hat:

  Vector with the elements of the latent field, i.e., eta_hat = A x_hat

- gk:

  Gradient in eta.

- Hk:

  Hessian in eta.

- a:

  Step length.

- Qx:

  Precision matrix for the prior of the latent field.

- strategy:

  Strategy to use to optimize. Now, line search strategy with
  quasi-newton algorithm is the only one avaliable.

- y:

  Vector with the response variable

- d:

  Number of categories.

## Value

g0 : Gradient in x_hat_new. A numeric vector with the gradient in
x_hat_new.

x_hat_new: New value of x after apply one iteration.

## Author

Joaquín Martínez-Minaya <jomarminaya@gmail.com>
