#' Computing the real Hessian in eta
#'
#'
#' `H_matrix_eta` computes the real Hessian in eta of -loglikelihood
#'
#' @param A Matrix which links eta with the latent field, i.e., eta = A x.
#' @param x Vector with the elements of the latent field, i.e., eta = A x.
#' @param y Vector with dimension 1 x d.
#'
#' @return Expected Hessian in eta.
#'
#' @author Joaquín Martínez-Minaya <\email{joaquin.martinez-minaya@@uv.es}>
H_matrix_eta <- function(A = A, x, y) {
    eta <- as.numeric(A %*% x)
    sum_exp <- sum(exp(eta))


    ### --- Elements of the diagonal --- ###
    H <- diag(-exp(eta) * (digamma(exp(eta)) - digamma(sum_exp)) - exp(2 * eta) * (trigamma(exp(eta)) - trigamma(sum_exp)) +
        exp(eta) * log(y))

    ### --- The rest of the elements --- ###
    for (i in 1:length(eta)) {
        for (j in 1:length(eta)) {
            if (i != j) {
                H[i, j] <- + exp(sum(eta[i], eta[j])) * trigamma(sum_exp)
            }
        }
    }
    -H
}
