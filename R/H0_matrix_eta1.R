#' Computing expected Hessian in eta
#'
#' `H0_matrix_x` computes the expected Hessian in eta of -loglikelihood
#'
#' @param A Matrix which links eta with the latent field, i.e., eta = A x.
#' @param x Vector with the elements of the latent field, i.e., eta = A x.
#' @param d Dimension
#' @param cores Number of cores for parallel computation. The package parallel is used.
#'
#' @return Expected Hessian in eta.
#' @import parallel
#' @author Joaquín Martínez-Minaya <\email{jomarminaya@@gmail.com}>

H0_matrix_eta_x <- function(A, x, d, cores) {
  eta <- Matrix::Matrix(as.numeric(A %*% x), ncol = d, byrow = TRUE)
  a1 <- apply(eta, 1, H0_matrix_eta1, d)
  parallel::mclapply(a1, "[[", 1, mc.cores = cores)


  # eta_list <- split(eta, rep(1:nrow(eta), each = ncol(eta)))
  # a <- parallel::mclapply(eta_list, H0_matrix_eta1, d, mc.cores = 4)

}



H0_matrix_eta1 <- function(eta, d) {
  #eta <- as.numeric(A %*% x)
  sum_exp <- sum(exp(eta))

  ### --- Elements of the diagonal --- ###
  H0 <- -diag(exp(2 * eta) * (trigamma(exp(eta)) - trigamma(sum_exp)))

  ### --- The rest of the elements --- ###
  for (i in 1:length(eta)) {
    for (j in 1:length(eta)) {
      if (i != j) {
        H0[i, j] <- +exp(sum(eta[i], eta[j])) * trigamma(sum_exp)
      }
    }
  }

  # t(A) %*% H0 %*% A
  list(-H0)
}

