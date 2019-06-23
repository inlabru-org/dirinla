predict_dirinla <- function(model, data.pred)
{
  if(!any(names(data.pred)=="intercept"))
  {
    data.pred <- cbind(intercept=1, data.pred)
  }
  sim <- model$marginals_fixed %>% purrr::map(function(x){
    sapply(x, inla.rmarginal, n=10000)})
  sim <- sim %>% purrr::map(function(x)as.matrix(t(x)))
  data.pred.cov <- lapply(model$marginals_fixed, function(x){
    as.matrix(dplyr::select(data.pred, names(x)))})
  predictive_etas <- Map('%*%', data.pred.cov, sim)
  predictive_alphas <- predictive_etas %>% purrr::map(function(x)exp(x))

  predictive_precision <- Reduce("+", predictive_alphas)
  predictive_means <- predictive_alphas %>% purrr::map(function(x)x/predictive_precision)

  #Outputs

  summary_predictive_alphas <- predictive_alphas %>% purrr::map(function(x) t(apply(x,1,summary)))
  marginals_predictive_alphas <- predictive_alphas %>% purrr::map(function(x){
    apply(x, 1, function(y){
      dens <- density(y, n=1000)
      data.frame(x=dens$x, y=dens$y)
    })
  })

  summary_predictive_means <- predictive_means %>% purrr::map(function(x) t(apply(x,1,summary)))
  marginals_predictive_means <- predictive_means %>% purrr::map(function(x){
    apply(x, 1, function(y){
      dens <- density(y, n=1000)
      data.frame(x=dens$x, y=dens$y)
    })
  })

  summary_predictive_precision <- t(apply(predictive_precision, 1, summary))
  marginals_predictive_precision <- apply(predictive_precision, 1, function(y){
    dens <- density(y, n=1000)
    data.frame(x = dens$x, y = dens$y)
  })



  #Checking Reduce("+", predictive_means)
  structure(list(call                           = model$call,
                 summary_fixed                  = model$summary_fixed,
                 marginals_fixed                = model$marginals_fixed,
                 summary_linear_predictor       = model$summary_linear_predictor,
                 marginals_linear_predictor     = model$marginals_linear_predictor,
                 summary_alphas                 = model$summary_alphas,
                 marginals_alphas               = model$marginals_alphas,
                 summary_precision              = model$summary_precision,
                 marginals_precision            = model$marginals_precision,
                 summary_means                  = model$summary_means,
                 marginals_means                = model$marginals_means,
                 summary_predictive_alphas      = summary_predictive_alphas,
                 marginals_predictive_alphas    = marginals_predictive_alphas,
                 summary_predictive_means       = summary_predictive_means,
                 marginals_predictive_means     = marginals_predictive_means,
                 summary_predictive_precision   = summary_predictive_precision,
                 marginals_predictive_precision = marginals_predictive_precision,
                 dic                            = model$dic,
                 waic                           = model$waic,
                 cpo                            = model$cpo,
                 nobs                           = model$n,
                 ncat                           = model$d
  ), class = "dirinlaregmodel")
}

