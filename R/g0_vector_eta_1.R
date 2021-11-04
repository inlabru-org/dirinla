#' Computing gradient vector in eta
#'
#' `g0_vector_eta` computes the gradient of -loglikelihood
#'
#' @param A Matrix which links eta with the latent field, i.e., eta = A x.
#' @param x Vector with the elements of the latent field, i.e., eta = A x.
#' @param y Vector with the response variable.
#'
#' @return A numeric vector with the gradient in eta.
#'
#' @author Joaquín Martínez-Minaya <\email{joaquin.martinez-minaya@@uv.es}>

g0_vector_eta_1 <- function(A = A, x, y) {
    eta <- Matrix::Matrix(as.numeric(A %*% x), ncol = dim(y)[2], byrow = TRUE)
    #g0 <- numeric()
    g0 <- t(apply(eta, 1, function(eta)exp(eta) * (digamma(exp(eta)) - digamma(sum(exp(eta))))))
    #sum_exp <- sum(exp(eta))
    #g0 <- exp(eta) * (digamma(exp(eta)) - digamma(sum_exp))
    #g0 <- matrix(g0, nrow = length(g0))
    g0_eta <- -as.numeric(t(-g0 + exp(eta) * log(y)))
    g0_eta
}
##########################################
