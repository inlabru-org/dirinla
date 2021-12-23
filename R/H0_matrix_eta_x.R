#' Computing expected Hessian in eta
#'
#' `H0_matrix_eta_x` computes the expected Hessian in eta of -loglikelihood
#'
#' @param eta Linear predictor resulting of the product Ax.
#' @param d Dimension.
#' @param cores Number of cores for parallel computation. The package parallel is used.
#'
#' @return Expected Hessian in eta.
#' @import parallel
#' @author Joaquín Martínez-Minaya <\email{jomarminaya@@gmail.com}>
H0_matrix_eta_x <- function(eta, d, cores) {
  #eta <- Matrix::Matrix(as.numeric(A %*% x), ncol = d, byrow = TRUE)
  a1 <- apply(eta, 1, H0_matrix_eta1, d)
  parallel::mclapply(a1, "[[", 1, mc.cores = cores)
}


#' Computing expected Hessian for a vector eta
#'
#' `H0_matrix_eta_1` computes the expected Hessian in eta of -loglikelihood
#'
#' @param eta eta vector to compute the expected Hessian.
#' @param d Dimension
#'
#' @return Expected Hessian in eta.
#' @author Joaquín Martínez-Minaya <\email{jomarminaya@@gmail.com}>
H0_matrix_eta1 <- function(eta, d) {
  sum_exp <- sum(exp(eta))

  ### --- The rest of the elements --- ###
  H0 <- exp(eta %*% matrix(rep(1,d), ncol = d) +
              matrix(rep(1,d), ncol = 1) %*% t(eta)) *
    trigamma(sum_exp)

  ### --- Elements of the diagonal --- ###
  diag(H0) <- -exp(2 * eta) * (trigamma(exp(eta)) - trigamma(sum_exp))


  # t(A) %*% H0 %*% A
  list(-H0)
}

