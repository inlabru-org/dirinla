#' Computing expected Hessian in eta
#'
#' `H0_matrix_eta` computes the expected Hessian in eta of -loglikelihood
#'
#' @param A: Matrix which links eta with the latent field, i.e., eta = A x.
#' @param x: Vector with the elements of the latent field, i.e., eta = A x.
#'
#' @return Expected Hessian in eta.
#'
#' @author Joaquín Martínez-Minaya <\email{joaquin.martinez-minaya@@uv.es}>
H0_matrix_eta <- function(A = A, x = x) {
    eta <- as.numeric(A %*% x)
    sum_exp <- sum(exp(eta))
    
    ### --- Elemenths of the diagonal --- ###
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
    -H0
}

