#' Computing additional diagonal part for the real Hessian H = H0 + diag
#'
#' `H_matrix_eta_diag` computes the expected Hessian in eta of -loglikelihood
#'
#' @param eta eta vector to compute the expected Hessian.
#' @param d Dimension
#' @param y Data corresponding to the i-individual
#'
#' @return Elements of the diagonal  such as H = H0 + diag
#' @importFrom Rfast rowsums
#' @author Joaquín Martínez-Minaya <\email{jomarminaya@@gmail.com}>
H_matrix_eta_diag <- function(eta, d, y) {

  sum_exp <- Rfast::rowsums(exp(eta))
  -(-exp(eta) * (digamma(exp(eta)) - digamma(sum_exp))  +
      exp(eta) * log(y))
}
