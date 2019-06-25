#' Fitting a Dirichlet regression
#'
#' `dirinlareg` Main function to do a Dirichlet Regression
#'
#' @param y: matrix containing the response variable R^{nxd}, being n number of individuals
#' and d the number of categories
#' @param data.cov: data.frame with the covarites, only the covariates!
#' @param share: parameters to be fitted jointly.
#' @param x0: initial optimization value
#' @param tol0: tolerance
#' @param tol1: tolerance for the gradient such that |grad| < tol1 * max(1, |f|)
#' @param k0: number of iterations
#' @param a: step length in the optimization algorithm
#' @param strategy: strategy to use to optimize
#' @param Qx: matrix of precision for the prior Gaussian field
#' @param ...:
#'
#' @return model: inla object
#' @return mean : posteriors means of the parameters corresponding to the latent variables
#'
#' @examples
#' ### In this example, we show how to fit a model using the dirinla package ###
#' ### --- 1. Loading the libraries --- ####
#' library(dirinla)
#' library(INLA)
#' library(DirichletReg)
#' library(ggplot2)
#' library(gridExtra)
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
#'
#'
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
  names_cat <- formula_list(formula)

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
    x0 <- rep(0, m)
  }

  # Creating precision matrix for the prior
  Qx <- Matrix(diag(prec, m))

  if(verbose==TRUE)
  {
    cat(paste0("\n \n ----------------------", " Looking for the mode ", "----------------- \n \n "))
  }


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
                            verbose  = verbose)


  ### --- 4. Second step : Include it in INLA.  --- ####

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



  data_stack_2 <- data_stack_dirich(y          = as.vector(z_eta),
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

  formula.inla.pred <- paste0("f(",
                              names_inla,
                              ", model = 'iid', hyper = list(theta = list(initial = log(",
                              prec,
                              "), fixed = TRUE)))")
  formula.inla.pred <- stringr::str_c(formula.inla.pred, collapse=" + ")
  formula.inla <- as.formula(stringr::str_c(formula.inla, formula.inla.pred, collapse = " " ))



  if(verbose==TRUE)
  {
    cat(paste0("\n ----------------------", " Fitting it using INLA ", "----------------- \n"))
  }

  mod0<-inla(formula.inla,
             family            = "gaussian",
             data              = inla.stack.data(data_stack_2),
             control.predictor = list(A = t(Lk_eta) %*% inla.stack.A(data_stack_2)),
             control.compute   = list(config = TRUE, #Compute marginals
                                      dic    = TRUE,
                                      waic   = TRUE,
                                      cpo    = TRUE),
             control.inla      = list(strategy = "laplace"),
             control.family    = list(hyper =
                                        list(prec =
                                               list(initial = log(1),
                                                    fixed   = TRUE))),
             verbose           = verbose)

  #Checking
  # ## --- New modes
  #x_hat2 <- sapply(mod0$summary.random, function(x){x$mean})
  #rownames(x_hat2) <- paste0("category", 1:d)
  #colnames(x_hat2)[1:(m)] <- names_cat[[1]]


  ### --- 5 Extracting posterior distributions --- ####
  ### --- 5.1. Extracting fixed effects --- ####
  fixed_effects <- extract_fixed(inla_model = mod0,
                                 names_cat  = names_cat)

  ### --- 5.2. Extracting linear predictor --- ####
  linear_predictor <- extract_linear_predictor(inla_model = mod0,
                                               n          = n,
                                               d          = d,
                                               Lk_eta     = Lk_eta)


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
                 ncat                           = d
                ), class = "dirinlaregmodel")
}




