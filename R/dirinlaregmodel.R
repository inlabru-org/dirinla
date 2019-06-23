#' Defining a new class
#'
#' `dirinlaregmodel` is a new object class
#'
#' @param call: The call of the function dirinlareg.
#' @param summary_fixed: List containing a summary of the marginal posterior distributions of the fixed effects.
#' @param marginals_fixed: List containing the marginal posterior distributions of the fixed effects.
#' @param summary_linear_predictor: List containing a summary of the marginal posterior distributions of the linear predictor.
#' @param marginals_linear_predictor: List containing the marginal posterior distributions of the linear predictor.
#' @param summary_alphas: List containing a summary of the marginal posterior distributions of the alphas.
#' @param marginals_alphas: List containing the marginal posterior distributions of the alphas.
#' @param summary_precision: List containing a summary of the marginal posterior distributions of the precision.
#' @param marginals_precision: List containing the marginal posterior distributions of the precision.
#' @param summary_means: List containing a summary of the marginal posterior distributions of the means.
#' @param marginals_means: List containing the marginal posterior distributions of the means.
#' @param dic: List containing the inla output for dic.
#' @param waic: List containing the inla output for waic.
#' @param cpo: List containing the inla output for cpo.
#' @param nobs: Number of observations.
#' @param ncat: Number of categories.
#'
#' @return object of list and dirinlaregmodel class.

dirinlaregmodel <- function(call = NULL, summary_fixed = NULL, marginals_fixed = NULL, summary_linear_predictor = NULL,
    marginals_linear_predictor = NULL, summary_alphas = NULL, marginals_alphas = NULL, summary_precision = NULL,
    marginals_precision = NULL, summary_means = NULL, marginals_means = NULL,
    summary_predictive_alphas      = NULL,
    marginals_predictive_alphas    = NULL,
    summary_predictive_means       = NULL,
    marginals_predictive_means     = NULL,
    summary_predictive_precision   = NULL,
    marginals_predictive_precision = NULL,
    dic = NULL, waic = NULL, cpo = NULL,
    nobs = NULL, ncat = NULL) {
    res <- list(call                           = call,
                summary_fixed                  = summary_fixed,
                marginals_fixed                = marginals_fixed,
                summary_linear_predictor       = summary_linear_predictor,
                marginals_linear_predictor     = marginals_linear_predictor,
                summary_alphas                 = summary_alphas,
                marginals_alphas               = marginals_alphas,
                summary_precision              = summary_precision,
                marginals_precision            = marginals_precision,
                summary_means                  = summary_means,
                marginals_means                = marginals_means,
                summary_predictive_alphas      = summary_predictive_alphas,
                marginals_predictive_alphas    = marginals_predictive_alphas,
                summary_predictive_means       = summary_predictive_means,
                marginals_predictive_means     = marginals_predictive_means,
                summary_predictive_precision   = summary_predictive_precision,
                marginals_predictive_precision = marginals_predictive_precision,
                dic                            = dic,
                waic                           = waic,
                cpo                            = cpo,
                nobs                           = nobs,
                ncat                           = ncat)

    ## Set the name for the class
    class(res) <- append(class(res), "dirinlaregmodel")
    return(res)
}





