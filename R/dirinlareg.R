#' Fitting a Dirichlet regression
#'
#' `dirinlareg` Main function to do a Dirichlet Regression
#'
#' @param formula object of class formula indicating the response variable and the covariates of the Dirichlet regression
#' @param y matrix containing the response variable R^{nxd}, being n number of individuals
#' and d the number of categories
#' @param data.cov data.frame with the covarites, only the covariates!
#' @param share parameters to be fitted jointly.
#' @param x0 initial optimization value
#' @param tol0 tolerance
#' @param tol1 tolerance for the gradient such that |grad| < tol1 * max(1, |f|)
#' @param k0 number of iterations
#' @param a step length in the optimization algorithm
#' @param strategy strategy to use to optimize
#' @param prec precision for the prior of the fixed effects
#' @param verbose if TRUE all the computing process is shown. Default is FALSE
#' @param cores Number of cores for parallel computation. The package parallel is used.
#' @param sim Simulations to call inla.posterior.sample and extract linear predictor, alphas and mus.
#' The bigger it is, better is the approximation, but more computational time.
#' @param prediction if TRUE we will predict with the new values of the covariates given in data.pred.cov.
#' @param data.pred.cov data.frame with the values for the covariates where we want to predict.
#' @param ... arguments for the inla command
#'
#' @return model dirinlaregmodel object
#'
#' @examples
#' #' ### In this example, we show how to fit a model using the dirinla package ###
#' ### --- 1. Loading the libraries --- ####
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
#'   verbose  = FALSE)
#'
#'
#' summary(model.inla)
#'
#' @export
#' @import stringr
#' @import INLA
#' @import samplingDataCRT
#' @importFrom purrr map
#' @author Joaquín Martínez-Minaya <\email{joaquin.martinez-minaya@@uv.es}>
dirinlareg <- function (formula,
                        y,
                        data.cov,
                        share    = NULL,
                        x0       = NULL,
                        tol0     = 1e-8,
                        tol1     = 0.01,
                        k0       = 100,
                        a        = 0.5,
                        strategy = "ls-quasi-newton",
                        prec     = prec,
                        verbose  = FALSE,
                        cores    = 1,
                        sim      = 1000,
                        prediction  = FALSE,
                        data.pred.cov = NULL,
                        ...)
{


  ### --- Some checkings --- ####
  this.call <- match.call()

  if(dim(y)[1] !=dim(data.cov)[1])
  {
    stop("The dimension of the response is different from the covariates")
  }


  d <- dim(y)[2] #dimension
  n <- dim(y)[1] #amount data

  ### --- Checking if there is some initial value
  # if(is.null(x0)){
  #   m <- 3
  #   x0 <- runif(d*m)
  # }

  ### --- 1. Reading the formula --- ####
  names_cat <- formula_list(formula, y = y)


  ### --- 2. Checking if fixed or random --- ####
  # Looking for the covariates in the data.frame for each category
  #data.cov.cat <- lapply(names_cat, function(x){dplyr::select(data.cov, x)} )

  ### --- 2. Store the data in a inla.stack --- ####
  data_stack <- data_stack_dirich(y          = as.vector(t(y)),
                                  covariates = names_cat,
                                  share      = share,
                                  data       = data.cov,
                                  d          = d,
                                  n          = n )

  A <- inla.stack.A(data_stack)


  ### --- 3. First step : Looking for the mode --- ####
  # m is the number of parameters to approximate
  m <- dim(A)[2]

  #Checking if there is initial condition
  if(is.null(x0)){
    #x0 <- rep(0, m)
    x0 <- rep(0, dim(A)[2])
  }

  # Creating precision matrix for the prior
  #Qx <- Matrix(diag(prec, m))
  #Check this prior. As we are giving priors for any realizations of the Gaussian field
  Qx <- Matrix(diag(prec, dim(A)[2]))
  cat(paste0("\n \n ----------------------", " Looking for the first mode ", "----------------- \n \n "))

# tol0 <- 1
# tol1 <- 1

  x_hat1 <- look_for_mode_x(A        = A,
                            x0       = as.vector(x0),
                            tol0     = tol0,
                            tol1     = tol1,
                            k0       = k0,
                            a        = a,
                            y        = y,
                            d        = d,
                            n        = n,
                            strategy = strategy,
                            Qx       = Qx,
                            verbose  = verbose,
                            cores    = cores)
####################################################################################
####################################################################################
  ### --- 4. Second step : Include it in INLA.  --- ####
  # La forma de introducir los datos en INLA la vamos a cambiar, ya que es un poco lioso
  # a la hora de introducir efectos aleatorios multiplicar por la matriz A
####################################################################################
####################################################################################
    ### --- The Hessian in the mode --- ###
  Hk_eta <- x_hat1$Hk

  ### --- Cholesky decomposition --- ###
  Lk_eta <- x_hat1$Lk

  ### --- The gradient in the mode --- ###
  gk_eta <- x_hat1$gk

  ### --- eta --- ###
  eta_eta <- x_hat1$eta

  ### --- The new variables conditioned to eta --- ###
  z_eta <- x_hat1$z


  ##############################################################################





  #############################################################################

  #
  data_stack_2 <- data_stack_dirich_formula(y          = as.vector(z_eta),
                                            covariates = names_cat,
                                            share      = share,
                                            data       = data.cov,
                                            d          = d,
                                            n          = n )


  ### ------- 4.1. Using the A matrix --- ####
  ### Create the formula
  formula.inla <- "y ~ -1 + "


  ### Names to introduce in INLA
  names_inla <- names(data_stack_2$effects$data)

  #Including fixed effects
  pos_fixed <- names_inla %>% stringr::str_starts("cat")
  names_inla_fixed <- names_inla[pos_fixed]
  formula.inla.pred <- character()
  if(length(names_inla_fixed) >=1){
    formula.inla.pred <- c(formula.inla.pred, paste0("f(",
                                                     names_inla_fixed,
                                                     ", model = 'linear')"))

    formula.inla.pred <- stringr::str_c(formula.inla.pred, collapse=" + ")
  }

  ############################################################################3
  ####### Revisar: incluyendo por defecto un efecto aleatorio compartido ######
  #Including random effects
  # names_inla_random <- names_inla[!pos_fixed]
  # if(length(names_inla_random) >=1)
  # {
  #   terms_random <- sapply(names_inla_random, function(x){
  #     unlist(names_cat) %>% grep(pattern = x, .) %>% unlist(names_cat)[.] -> res
  #   res[1] %>% as.character()
  #   })
  #
  #   #Check for the effects which has the index
  #   terms_random %>% paste(., collapse = "+") %>%
  #     paste(formula.inla.pred, ., sep = "+") -> formula.inla.pred
  #
  #
  # }
  #formula.inla.pred <- paste(formula.inla.pred, "f(iid1, model = 'iid')", sep = "+")
  formula.inla <- as.formula(paste(formula.inla, formula.inla.pred, collapse = " " ))



  cat(paste0("\n ----------------------", "    INLA call    ", "----------------- \n"))

  ### This is for the linear predictor using lincomb
  #Inverse of Lk_eta
  # Lk_eta_inv <- solve(t(Lk_eta))
  # all3 <- sapply(1:dim(Lk_eta_inv)[2], function(x){
  #   a <- inla.make.lincomb(APredictor = Lk_eta_inv[x,])
  #   names(a) <- c(paste0("lc", x))
  #   a})

  #Faltaría multiplicar pos las realizaciones del efecto aleatorio
  mod0 <- inla(formula.inla,
             family            = "gaussian",
             data              = inla.stack.data(data_stack_2),
             control.predictor = list(A = t(Lk_eta), compute = TRUE),
             control.compute   = list(config = TRUE, #Compute marginals
                                      dic    = TRUE,
                                      waic   = TRUE,
                                      cpo    = TRUE),
             #control.inla      = list(strategy = "gaussian"),
             control.family    = list(hyper =
                                        list(prec =
                                               list(initial = log(1),
                                                    fixed   = TRUE))),
             verbose           = verbose,
             num.threads = cores)

  #Checking
  # ## --- New modes
  #x_hat2 <- sapply(mod0$summary.random, function(x){x$mean})
  #rownames(x_hat2) <- paste0("category", 1:d)
  #colnames(x_hat2)[1:(m)] <- names_cat[[1]]

  summary(mod0)


  ### --- 5 Extracting posterior distributions --- ####
  ### --- 5.1. Extracting fixed effects --- ####
  fixed_effects <- extract_fixed(inla_model = mod0,
                                 names_cat  = names_cat)

  cat(paste0("\n ----------------------", " Obtaining linear predictor ", "----------------- \n"))

  ### --- 5.2. Extracting linear predictor --- ####
  linear_predictor <- extract_linear_predictor(inla_model = mod0,
                                               n          = n,
                                               d          = d,
                                               Lk_eta     = Lk_eta,
                                               names_cat  = names_cat,
                                               sim        = sim,
                                               verbose    = verbose,
                                               cores      = cores)

  ### --- Prediction --- ####
  if(prediction == TRUE)
  {
    model.inla <- structure(list(call                           = this.call,
                   summary_fixed                  = fixed_effects$summary_fixed,
                   marginals_fixed                = fixed_effects$marginals_fixed,
                   summary_linear_predictor       = linear_predictor$summary_linear_predictor,
                   marginals_linear_predictor     = linear_predictor$marginals_linear_predictor,
                   summary_alphas                 = linear_predictor$summary_alphas,
                   marginals_alphas               = linear_predictor$marginals_alphas,
                   summary_precision              = linear_predictor$summary_precision,
                   marginals_precision            = linear_predictor$marginals_precision,
                   summary_means                  = linear_predictor$summary_means,
                   marginals_means                = linear_predictor$marginals_means,
                   summary_predictive_alphas      = NULL,
                   marginals_predictive_alphas    = NULL,
                   summary_predictive_means       = NULL,
                   marginals_predictive_means     = NULL,
                   summary_predictive_precision   = NULL,
                   marginals_predictive_precision = NULL,
                   dic                            = mod0$dic,
                   waic                           = mod0$waic,
                   cpo                            = mod0$cpo,
                   nobs                           = n,
                   ncat                           = d,
                   y                              = y,
                   data.cov                       = data.cov
    ), class = "dirinlaregmodel")

    model.prediction <-
      predict.dirinlaregmodel(model.inla,
              data.pred.cov = data.pred.cov)
    summary_predictive_alphas       <- model.prediction$summary_predictive_alphas
    marginals_predictive_alphas    <- model.prediction$marginals_predictive_alphas
    summary_predictive_means       <- model.prediction$summary_predictive_means
    marginals_predictive_means     <- model.prediction$marginals_predictive_means
    summary_predictive_precision   <- model.prediction$summary_predictive_precision
    marginals_predictive_precision <- model.prediction$marginals_predictive_precision
  }else{
    summary_predictive_alphas       <- NULL
    marginals_predictive_alphas    <- NULL
    summary_predictive_means       <- NULL
    marginals_predictive_means     <- NULL
    summary_predictive_precision   <- NULL
    marginals_predictive_precision <- NULL
  }


  structure(list(call                           = this.call,
                 summary_fixed                  = fixed_effects$summary_fixed,
                 marginals_fixed                = fixed_effects$marginals_fixed,
                 summary_linear_predictor       = linear_predictor$summary_linear_predictor,
                 marginals_linear_predictor     = linear_predictor$marginals_linear_predictor,
                 summary_alphas                 = linear_predictor$summary_alphas,
                 marginals_alphas               = linear_predictor$marginals_alphas,
                 summary_precision              = linear_predictor$summary_precision,
                 marginals_precision            = linear_predictor$marginals_precision,
                 summary_means                  = linear_predictor$summary_means,
                 marginals_means                = linear_predictor$marginals_means,
                 summary_predictive_alphas      = summary_predictive_alphas,
                 marginals_predictive_alphas    = marginals_predictive_alphas,
                 summary_predictive_means       = summary_predictive_means,
                 marginals_predictive_means     = marginals_predictive_means,
                 summary_predictive_precision   = summary_predictive_precision,
                 marginals_predictive_precision = marginals_predictive_precision,
                 dic                            = mod0$dic,
                 waic                           = mod0$waic,
                 cpo                            = mod0$cpo,
                 nobs                           = n,
                 ncat                           = d,
                 y                              = y,
                 data.cov                       = data.cov
                ), class = "dirinlaregmodel")
}




