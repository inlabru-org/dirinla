# Extracting posterior distributions from the linear predictor

`extract_linear_predictor` extracts the posterior distribution from the
linear predictor

## Usage

``` r
extract_linear_predictor(
  inla_model,
  n,
  d,
  Lk_eta,
  names_cat = names_cat,
  sim,
  verbose,
  cores
)
```

## Arguments

- inla_model:

  An object of class inla.

- n:

  Number of observations.

- d:

  Number of categories.

- Lk_eta:

  Cholesky decomposition of the Hessian matrix.

- names_cat:

  List generated with extract_formula.

- sim:

  simulations for the function inla.posterior.sample

- verbose:

  if TRUE all the computing process is shown. Default is FALSE

- cores:

  number of cores to be used in the computations

## Value

summary_linear_predictor List containing a summary of the marginal
posterior distributions of the linear predictor.

marginals_linear_predictor List containing simulations of marginal
posterior distributions of the linear predictor.

summary_alphas List containing a summary of the marginal posterior
distributions of the alphas.

marginals_alphas List containing simulations of the marginal posterior
distributions of the alphas.

summary_precision List containing a summary of the marginal posterior
distributions of the precision.

marginals_precision List containing simulations of the marginal
posterior distributions of the precision.

summary_means List containing a summary of the marginal posterior
distributions of the means.

marginals_means List containing the simulations of the marginal
posterior distributions of the means.

## Author

Joaquín Martínez-Minaya <jomarminaya@gmail.com>
