#' Finding the mode of the full posterior distribution
#'
#' `look_for_mode_x` computes optimization algorithms to find the mode of the posterior
#'
#' @param A Matrix which links latent field with linear predictor.
#' @param x0 Initial optimization value.
#' @param tol0 Tolerance for |x_new - x_old| and |f_new - f_old|.
#' @param tol1 Tolerance for the gradient such that |grad| < tol1 * max(1, |f|)
#' @param k0 Number of iterations.
#' @param a Step length in the algorithm.
#' @param y Response variable. Number of columns correspond to the number of categories.
#' @param d Number of categories.
#' @param n Number of individuals.
#' @param strategy Strategy to use to optimize.
#' @param Qx Prior precision matrix for the fixed effects.
#' @param verbose By default is FALSE. If TRUE, the computation process is shown in the scream.
#' @param cores Number of cores for parallel computation. The package parallel is used.

#' @return x_hat Matrix with the x of the iterations.
#' @return Hk Hessian in eta. This Hessian is a combination of the real Hessian (when
#' it is positive definite) and the expected Hessian (when the real Hessian is not positive
#' definite).
#' @return gk Gradient in eta.
#' @return Lk Cholesky decomposition matrix.
#' @return eta Linear predictor.
#' @return z New pseudo observation conditioned to eta.
#'
#'
#' @author Joaquín Martínez-Minaya <\email{joaquin.martinez-minaya@@uv.es}>


look_for_mode_x <- function(A = A,
                            x0,
                            tol0,
                            tol1,
                            k0,
                            a = 0.5,
                            y,
                            d,
                            n,
                            strategy = "ls-quasi-newton",
                            Qx,
                            verbose,
                            cores) {
  x_hat <- matrix(ncol = length(x0), nrow = 1)
  x_hat[1, ] <- x0
  k <- 1
  less <- FALSE
  ### --- Looking for using the expected Hessian --- ###
  while ((!less == TRUE) && (k<k0)) {
    ### --- Call the function to define new variables in each iteration --- ###
    Hk_list <- list()
    gk <- numeric()
    gk <- g0_vector_eta_1(
        A = A,
        x = x_hat[k, ],
        y = y
      )

    Hk_list <- H0_matrix_eta_x(eta = Matrix::Matrix(as.numeric(A %*% x_hat[k, ]),
                                                    ncol = d,
                                                    byrow = TRUE),
                               d,
                               cores = cores)

    if (any(!is.finite(gk))) {
      gk[which(!is.finite(gk))] <- exp(100)
      print("You should apply the transformation DR_data")
    }
    ### Expected Hessian ####
    ### --- Diagonal matrix dim=dn x dn with the expected hessian in the diagonal
    #Hk <- Matrix(blockMatrixDiagonal(Hk_list))
    #Hk <- Matrix::bdiag_m(Hk_list)
    #Hk <- blockMatrixDiagonal(Hk_list) Too much time
    #Hk <- Matrix::bdiag(Hk_list)
    Hk <- bdiag_m(Hk_list)

    ### --- Iterative method --- ###
    x_g0_hat_new <- newton_x(
      A     = A,
      x_hat = x_hat[k, ],
      gk,
      Hk,
      a,
      Qx,
      strategy = strategy,
      y = y
    )
    x_hat_new <- x_g0_hat_new$x_hat_new
    gk_new <- x_g0_hat_new$gk
    x_hat <- rbind(x_hat, x_hat_new)

    ### --- Loglikelihoods
    f_new <- dirichlet_log_pos_x(
      A = A,
      x = x_hat_new,
      Qx = Qx,
      y = y
    )

    f_old <- dirichlet_log_pos_x(
      A = A,
      x = x_hat[k, ],
      Qx = Qx,
      y = y
    )


    ### Checking condition
    less <- (as.logical(abs(sum(x_hat_new - x_hat[k, ])) < tol0) && #Condition in x
      as.logical(abs(f_new - f_old) < tol0)) && #Condition in f
      as.logical(abs(sum(gk)) < tol1 * max(1, abs(f_new)))

    if(verbose == TRUE)
    {
      cat(paste0(
        "Iter = ", k,
        ", |grad| = ", round(abs(sum(gk)), 2),
        ", log.post = ", round(f_new, 2),
        ", |x_new - x_old| = ", round(abs(sum(x_hat_new - x_hat[k, ])), 5),
        ", |f_new - f_old| = ", round(abs(f_new - f_old), 5),
        "\n"
        #", x = ", paste(round(x_hat_new, 2), collapse = " ")
      ))
    }

    # print(eigen(Hk)$values)
    k <- k + 1
  }


  #cat("Computing Hessian in the mode \n")
  ### --- Compute the new variables using the real Hessian when it is positive definite --- ###
  ### --- if not, we use the expected Hessian                                           --- ###
  count_exp <- 0 #Count the times that expected hessian is used

  Hk_list <- H0_list <- H_list <- list()
  gk <- numeric()
  gk <- g0_vector_eta_1(
    A = A,
    x = x_hat[k, ],
    y = y
  )

  eta = matrix(as.numeric(A %*% x_hat[k, ]),
                       ncol = d,
                       byrow = TRUE)

  H0_list <- H0_matrix_eta_x(eta = eta,
                             d,
                             cores = cores)

  H_diag <- H_matrix_eta_diag(eta = eta, d = d, y = y)


  lapply(1:length(H0_list), function(x){
    H_list[[x]] <- H0_list[[x]] + diag(H_diag[x,])
    real <- TRUE
    if(any(eigen(H_list[[x]])$values < 0))
    {
      H_list[[x]] <- H0_list[[x]]
      real <- FALSE
    }
     list(block = H_list[[x]], real = real)
    }) -> H_list


  Hk_list <- lapply(H_list, "[[", 1)
  count_exp <- table(unlist(lapply(H_list, "[[", 2)))[1]


  # for (i in 1:n)
  # {
  #
  #   ## Auxiliar matrix to determina if the real Hessian is positive definite
  #   Hk_list[[i]] <- H_matrix_eta(
  #     A = A[(d * (i - 1) + 1):(i * d), ],
  #     x_hat[k, ],
  #     y = y[i, ]
  #   )
  #   if (any(eigen(Hk_list[[i]])$values < 0)) {
  #     Hk_list[[i]] <- H0_matrix_eta(
  #       A = A[(d * (i - 1) + 1):(i * d), ],
  #       x_hat[k, ]
  #     )
  #     #cat("Expected hessian is used \n")
  #     count_exp <- count_exp + 1
  #   }
  # }

  if(verbose==TRUE)
  {
    if(less == TRUE)
    {
      cat("\nGreat news! The mode has been properly located!")
    }else{
      cat("\nBad news! You should increase the number of iterations (k0), or increase tol0 and tol1.")
    }

    cat(paste0("\n \nReal Hessian has been used ", n - count_exp , " times \n"))
    cat(paste0("Expected Hessian has been used ", count_exp, " times \n"))

  }
  ### --- Diagonal matrix dim=dn x dn with the real hessian in the diagonal
  #Hk <- Matrix(blockMatrixDiagonal(Hk_list))
  Hk <- bdiag_m(Hk_list)

  ### --- Compute the new variables z --- ###
  Lk <- chol(Hk)
  Lk <- t(Lk)
  gk <- gk

  ### --- Copute z --- ###
  eta <- A %*% x_hat[k, ]
  z <- t(Lk) %*% (eta) - solve(Lk) %*% gk

  ### --- Returning --- ###
  list(
    x_hat = as.vector(x_hat[k, ]),
    Hk    = Hk,
    gk    = gk,
    Lk    = Lk,
    eta   = as.vector(eta),
    z     = z
  )
}

