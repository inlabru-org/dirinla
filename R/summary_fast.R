#' Summary using the packages Rfast and Rfast2 of a matrix by rows
#'
#' `summary_fast` summarise a matrix by rows
#'
#' @param A matrix to be summarised
#'
#' @return A matrix whose columns are "mean", "min", "q0.025", "q0.25", "q0.5", "q0.75", "q0.975", "max"
#'
#' @importFrom Rfast rowmeans rowMinsMaxs rowMedians
#' @importFrom Rfast2 rowQuantile
#'
#' @example
#' A <- matrix(rnorm(10000), ncol = 1000)
#' summary_fast(A)
#'
#' @author Joaquín Martínez-Minaya <\email{jomarminaya@@gmail.com}>
summary_fast <- function(A){
    result <- t(
        rbind(
            Rfast::rowmeans(A),
            Rfast::rowMinsMaxs(A),
            Rfast::rowMedians(A),
            t(Rfast2::rowQuantile(A,
                                probs = c(0.025, 0.25, 0.75, 0.975)))))
    colnames(result) <- c("mean", "min", "max", "q0.5", "q0.025", "q0.25", "q0.75", "q0.975")
    result <- result[,c("mean", "min", "q0.025", "q0.25", "q0.5", "q0.75", "q0.975", "max")]
    result
}

