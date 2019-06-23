#' Draw samples from the prior distribution
#'
#' `draw_data` Function to draw samples for a multivariate
#' normal prior distribution.
#'
#' @param n: the number of simulations.
#' @param x: parameters of the latent field.
#' @param d: the dimension of the response variable
#' @param V: matrix with covariates corresponding with the parameters
#'
#' @return data.frame with the simulations
#'
#' @example
#' ## Generating 100 samples from a Dirichlet distribution with two covariates in
#' each category and the intercept.
#' y_ic ~ Dirichlet(alpha_i1, alpha_i2, alpha_i3)
#'        log(alpha_i1) = x[1] * v1 + x[1+3] * v2 + x[1+2*3]
#'        log(alpha_i2) = x[2] * v1 + x[2+3] * v2 + x[2+2*3]
#'        log(alpha_i3) = x[3] * v1 + x[3+3] * v2 + x[3+2*3]
#'
#' Parameters
#' x <- c(-2, -2, -2,
#'        + 2, + 2, + 2,
#'        1 ,1 ,1)
#'
#' Simulations
#' n <- 100
#'
#' Generate covariates
#' V <- as.data.frame(matrix(rnorm(2*n, 0, 1), ncol=2))
#' names(V) <- paste0('v', 1:2)
#'
#' #Dimension
#' d <- 3
#'
#' y_o <- draw_data(x = x,
#'                  d = d,
#'                  V = V)
#' boxplot(y_o)
#'
#' @author Joaquín Martínez-Minaya <\email{joaquin.martinez-minaya@@uv.es}>



draw_data <- function(x, d, V) {
    print("The number of simulations will be the number of data of the covariates")
    n <- dim(V)[1]
    x <- Matrix(x, ncol = 1)
    
    ### --- Data stack to construct the data
    data_stack_construct <- data_stack_dirich(y = rep(NA, d * n), covariates = names(V), data = V, d = d, n = n)
    A_construct <- data_stack_construct$A
    
    # Ordering the data with covariates --- ####
    
    eta <- A_construct %*% x
    alpha <- exp(eta)
    alpha <- matrix(alpha, ncol = d, byrow = TRUE)
    y_o <- rdirichlet(n, alpha)
    y_o
    
}

