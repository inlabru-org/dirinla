# Finding the mode of the full posterior distribution

`look_for_mode_x` computes optimization algorithms to find the mode of
the posterior

## Usage

``` r
look_for_mode_x(
  A = A,
  x0,
  tol0,
  tol1,
  k0,
  a = 0.5,
  y,
  d,
  n,
  strategy = "ls-quasi-newton",
  Qx,
  verbose,
  cores
)
```

## Arguments

- A:

  Matrix which links latent field with linear predictor.

- x0:

  Initial optimization value.

- tol0:

  Tolerance for \|x_new - x_old\| and \|f_new - f_old\|.

- tol1:

  Tolerance for the gradient such that \|grad\| \< tol1 \* max(1, \|f\|)

- k0:

  Number of iterations.

- a:

  Step length in the algorithm.

- y:

  Response variable. Number of columns correspond to the number of
  categories.

- d:

  Number of categories.

- n:

  Number of individuals.

- strategy:

  Strategy to use to optimize.

- Qx:

  Prior precision matrix for the fixed effects.

- verbose:

  By default is FALSE. If TRUE, the computation process is shown in the
  scream.

- cores:

  Number of cores for parallel computation. The package parallel is
  used.

## Value

x_hat Matrix with the x of the iterations.

Hk Hessian in eta. This Hessian is a combination of the real Hessian
(when it is positive definite) and the expected Hessian (when the real
Hessian is not positive definite).

gk Gradient in eta.

Lk Cholesky decomposition matrix.

eta Linear predictor.

z New pseudo observation conditioned to eta.

## Author

Joaquín Martínez-Minaya <jomarminaya@gmail.com>
