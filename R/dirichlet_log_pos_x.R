#' Dirichlet log posterior function
#'
#' `dirichlet_log_pos_x` returns the -log posterior Dirichlet distribution asumming
#' multivariate normal prior with precision matrix Qx for elements of the latent field.
#'
#'
#' @param A: A matrix which links eta with the latent field, i.e., eta = A x.
#' @param x: Vector with the elements of the latent field, i.e., eta = A x.
#' @param Qx: Precision matrix for the priors of the latent field.
#' @param y: Vector with the response variable.
#'
#' @return A real value showing the -log posterior density is returned
#'
#'
#' @author Joaquín Martínez-Minaya <\email{joaquin.martinez-minaya@@uv.es}>
dirichlet_log_pos_x <- function(A = A, x, Qx = Qx, y) {
    eta_hat <- A %*% x
    eta <- matrix(eta_hat,
                  ncol  = dim(y)[2],
                  nrow  = dim(y)[1],
                  byrow = TRUE)
    as.numeric(sum(apply(eta, 1, function(z)
      {
      log_beta_mult_eta(exp(z))}) -
          apply( (log(y) * (exp(eta) - 1) ), 1, sum)) +
          t(x) %*% Qx %*% x)
}

