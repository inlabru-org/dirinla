# Defining a new class

`dirinlaregmodel` is a new object class

## Usage

``` r
dirinlaregmodel(
  call = NULL,
  formula = NULL,
  summary_fixed = NULL,
  marginals_fixed = NULL,
  summary_random = NULL,
  marginals_random = NULL,
  summary_hyperpar = NULL,
  marginals_hyperpar = NULL,
  summary_linear_predictor = NULL,
  marginals_linear_predictor = NULL,
  summary_alphas = NULL,
  marginals_alphas = NULL,
  summary_precision = NULL,
  marginals_precision = NULL,
  summary_means = NULL,
  marginals_means = NULL,
  summary_predictive_alphas = NULL,
  marginals_predictive_alphas = NULL,
  summary_predictive_means = NULL,
  marginals_predictive_means = NULL,
  summary_predictive_precision = NULL,
  marginals_predictive_precision = NULL,
  dic = NULL,
  waic = NULL,
  cpo = NULL,
  nobs = NULL,
  ncat = NULL,
  y = NULL,
  data.cov = NULL
)
```

## Arguments

- call:

  The call of the function dirinlareg.

- formula:

  Formula introduced by the user.

- summary_fixed:

  List containing a summary of the marginal posterior distributions of
  the fixed effects.

- marginals_fixed:

  List containing the marginal posterior distributions of the fixed
  effects.

- summary_random:

  List containing a summary of the marginal posterior distributions of
  the random effects.

- marginals_random:

  List containing the marginal posterior distributions of the random
  effects.

- summary_hyperpar:

  List containing a summary of the marginal posterior distributions of
  the hyperparameters.

- marginals_hyperpar:

  List containing the marginal posterior distributions of the
  hyperparameters.

- summary_linear_predictor:

  List containing a summary of the marginal posterior distributions of
  the linear predictor.

- marginals_linear_predictor:

  List containing the marginal posterior distributions of the linear
  predictor.

- summary_alphas:

  List containing a summary of the marginal posterior distributions of
  the alphas.

- marginals_alphas:

  List containing the marginal posterior distributions of the alphas.

- summary_precision:

  List containing a summary of the marginal posterior distributions of
  the precision.

- marginals_precision:

  List containing the marginal posterior distributions of the precision.

- summary_means:

  List containing a summary of the marginal posterior distributions of
  the means.

- marginals_means:

  List containing the marginal posterior distributions of the means.

- summary_predictive_alphas:

  List containing a summary of the marginal posterior predictive
  distribution of the alphas.

- marginals_predictive_alphas:

  List containing the marginal posterior predictive distribution of the
  alphas.

- summary_predictive_means:

  List containing a summary fo the marginal posterior predictive
  distribution of the means.

- marginals_predictive_means:

  List containing the marginal posterior predictive distribution of the
  means.

- summary_predictive_precision:

  List containing a summary of the marginal posterior predictive
  distribution of the precision.

- marginals_predictive_precision:

  List containing the marginal posterior predictive distribution of the
  precision.

- dic:

  List containing the inla output for dic.

- waic:

  List containing the inla output for waic.

- cpo:

  List containing the inla output for cpo.

- nobs:

  Number of observations.

- ncat:

  Number of categories.

- y:

  matrix containing the response variable \\R^{nxd}\\, being n number of
  individuals and d the number of categories

- data.cov:

  data.frame with the covarites, only the covariates!

## Value

object of list and dirinlaregmodel class.
