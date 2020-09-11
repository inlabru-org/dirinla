#' Finding the mode of the full posterior distribution
#'
#' `predict.dirinlaregmodel` computes the posterior predictive distribution for some given values of the covariates
#'
#' @param model dirinlaregmodel object.
#' @param data.pred.cov Data.frame with the covariate values for the variables to predict.
#'
#' @return model dirinlaregmodel object
#'
#' @examples
#' ### In this example, we show how to fit a model using the dirinla package ###
#' ### --- 1. Loading the libraries --- ####
#' library(dirinla)
#' library(INLA)
#' library(DirichletReg)
#'
#'
#' ### --- 2. Simulating from a Dirichlet likelihood --- ####
#' set.seed(1000)
#' N <- 50 #number of data
#' V <- as.data.frame(matrix(runif((4) * N, 0, 1), ncol = 4)) #Covariates
#' names(V) <- paste0('v', 1:4)
#'
#' formula <- y ~ 1 + v1 | 1 + v2 | 1 + v3 | 1 + v4
#' (names_cat <- formula_list(formula))
#'
#' x <- c(-1.5, 1, -3, 1.5,
#'        2, -3 , -1, 5)
#'
#' mus <- exp(x) / sum(exp(x))
#' C <- length(names_cat)
#' data_stack_construct <-
#'   data_stack_dirich(y = as.vector(rep(NA, N * C)),
#'                     covariates = names_cat,
#'                     data       = V,
#'                     d          = C,
#'                     n          = N)
#'
#' A_construct <- data_stack_construct$A
#' A_construct[1:8, ]
#'
#' eta <- A_construct %*% x
#' alpha <- exp(eta)
#' alpha <- matrix(alpha,
#'                 ncol  = C,
#'                 byrow = TRUE)
#' y_o <- rdirichlet(N, alpha)
#' colnames(y_o) <- paste0("y", 1:C)
#' head(y_o)
#'
#'
#' ### --- 3. Fitting the model --- ####
#' y <- y_o
#' model.inla <- dirinlareg(
#'   formula  = y ~ 1 + v1 | 1 + v2 | 1 + v3 | 1 + v4,
#'   y        = y,
#'   data.cov = V,
#'   prec     = 0.0001,
#'   verbose  = TRUE)
#'
#'
#' summary(model.inla)
#' ### --- 4. Predicting for v1 = 0.25, v2 = 0.5, v3 = 0.5, v4 = 0.1 --- ####
#' model.prediction <- predict.dirinla(model.inla,
#'                 data.pred.cov= data.frame(v1 = 0.25,
#'                                        v2 = 0.5,
#'                                        v3 = 0.5,
#'                                        v4 = 0.1))
#' model.prediction$summary_predictive_means
#'
#' @export
#' @author Joaquín Martínez-Minaya <\email{jomarminaya@@gmail.com}>
predict.dirinlaregmodel <- function(model, data.pred.cov)
{
  if(!any(names(data.pred.cov)=="intercept"))
  {
    data.pred.cov<- cbind(intercept=1, data.pred.cov)
  }

  if(any(sapply(model$marginals_fixed, function(x){names(x)}) %in% colnames(data.pred.cov) == FALSE)){
    stop("Names of the variables in the data.frame does not match with the variables in the formula.")
  }

  cat(paste0("\n \n ----------------------", " Predicting ", "----------------- \n \n "))
  sim <- model$marginals_fixed %>% purrr::map(function(x){
    sapply(x, inla.rmarginal, n=10000)})
  sim <- sim %>% purrr::map(function(x)as.matrix(t(x)))
  data.pred.cov <- lapply(model$marginals_fixed, function(x){
    as.matrix(dplyr::select(data.pred.cov, names(x)))})
  predictive_etas <- Map('%*%', data.pred.cov, sim)
  predictive_alphas <- predictive_etas %>% purrr::map(function(x)exp(x))

  predictive_precision <- Reduce("+", predictive_alphas)
  predictive_means <- predictive_alphas %>% purrr::map(function(x)x/predictive_precision)

  #Outputs

  summary_predictive_alphas <- predictive_alphas %>% purrr::map(function(x) t(apply(x,1,summary)))
  # marginals_predictive_alphas <- predictive_alphas %>% purrr::map(function(x){
  #   apply(x, 1, function(y){
  #     dens <- density(y, n=1000)
  #     data.frame(x=dens$x, y=dens$y)
  #   })
  # })

  summary_predictive_means <- predictive_means %>% purrr::map(function(x) t(apply(x,1,summary)))
  # marginals_predictive_means <- predictive_means %>% purrr::map(function(x){
  #   apply(x, 1, function(y){
  #     dens <- density(y, n=1000)
  #     data.frame(x=dens$x, y=dens$y)
  #   })
  # })

  summary_predictive_precision <- t(apply(predictive_precision, 1, summary))
  # marginals_predictive_precision <- apply(predictive_precision, 1, function(y){
  #   dens <- density(y, n=1000)
  #   data.frame(x = dens$x, y = dens$y)
  # })



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
                 marginals_predictive_alphas    = predictive_alphas,
                 summary_predictive_means       = summary_predictive_means,
                 marginals_predictive_means     = predictive_means,
                 summary_predictive_precision   = summary_predictive_precision,
                 marginals_predictive_precision = predictive_precision,
                 dic                            = model$dic,
                 waic                           = model$waic,
                 cpo                            = model$cpo,
                 nobs                           = model$n,
                 ncat                           = model$d
  ), class = "dirinlaregmodel")
}

