#' Newton-Raphson algorithm
#'
#' `newton_x` computes optimization algorithms to find the mode of the
#' posterior. Line search strategy with Armijo conditions is implemented
#'
#' @param A: Matrix which links eta with the latent field, i.e., eta = A x
#' @param x_hat: Vector with the elements of the latent field, i.e., eta_hat = A x_hat
#' @param gk: Gradient in eta.
#' @param Hk: Hessian in eta.
#' @param a: Step length.
#' @param Qx: Precision matrix for the prior of the latent field.
#' @param strategy: Strategy to use to optimize. Now, line search strategy with quasi-newton algorithm is the only one avaliable.
#' @param y: Vector with the response variable
#' @param d: Number of categories.
#'
#' @return g0 : Gradient in x_hat_new. A numeric vector with the gradient in x_hat_new.
#' @return x_hat_new: New value of x after apply one iteration.
#'
#' @author Joaquín Martínez-Minaya <\email{joaquin.martinez-minaya@@uv.es}>
newton_x <- function(A, x_hat, gk, Hk, a, Qx, strategy, y, d = d) {
    x_hat <- as.numeric(x_hat)
    
    ### --- Cholesky decomposition --- ###
    Lk <- t(chol(Hk))
    
    Hkx <- t(A) %*% Hk %*% A
    gkx <- t(A) %*% gk
    
    
    ### --- Quasi-newton without line search
    if (strategy == "quasi-newton") {
        # Qx_z <- Hkx + Qx x_hat_new <- as.numeric(t(matrix(x_hat, ncol = 1) - a * ( solve(Qx_z, gkx + (Qx %*%
        # x_hat))))) gk_x_new <- gkx + (Qx %*% x_hat_new)
        
        ### --- Line search direction quasi-newton and Armijo conditions --- ###
    } else if (strategy == "ls-quasi-newton") {
        f0 <- dirichlet_log_pos_x(A = A, x_hat, Qx, y)
        g0 <- gkx + (Qx %*% x_hat)
        Qx_z <- Hkx + Qx
        # p <- -as.vector(solve(as.matrix(Qx_z), as.matrix(g0)))
        p <- -as.vector(solve(Qx_z, g0))
        
        aa <- 1
        found <- FALSE
        
        j <- 1
        
        
        while (!found & (j < 20)) {
            x_hat_new <- x_hat + aa * p
            f <- suppressWarnings(dirichlet_log_pos_x(A = A, x_hat_new, Qx, y))
            
            ## Checking if f is finite. If not, we reduce aa
            if (is.finite(f)) {
                found <- as.logical(f < f0 - aa * 1e-06 * abs(sum(p * g0)/sqrt(sum(p * p))))
                aa <- aa * a
                j <- j + 1
            } else {
                aa <- aa/2
            }
        }
        gk_x_new <- gkx + (Qx %*% x_hat_new)
    }
    
    output <- list(gk_x_new, x_hat_new)
    names(output) <- c("gk", "x_hat_new")
    output
}
