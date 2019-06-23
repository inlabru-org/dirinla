#' Draw samples from the prior distribution
#'
#' `draw_prior` draws samples from a multivariate normal prior distribution.
#'
#' @param n: Number of samples.
#' @param Q: Precision matrix of the multivariate normal distribution.
#' @param ...: Arguments of the mvrnorm function.
#'
#' @return Data.frame with the simulations.
#' @example
#' draw_prior(10, Q = Matrix(diag(rep(0.1, 3))))
#' @author Joaquín Martínez-Minaya <\email{joaquin.martinez-minaya@@uv.es}>

draw_prior <- function(n, Q) {
    Sigma <- sqrt(solve(Q))
    len <- dim(Q)[1]
    mvrnorm(n = n, mu = rep(0, len), Sigma = Sigma)
}



