#' Summary using the packages Rfast and Rfast2 of a matrix by rows
#'
#' `summary_fast` summarise a matrix by rows
#'
#' @param A matrix to be summarised
#'
#' @return A matrix whose columns are "mean", "sd", "0.025quant", "0.5quant", "0.975quant"
#'
#' @importFrom Rfast rowmeans rowMinsMaxs rowMedians
#' @importFrom Rfast2 rowQuantile
#'
#' @examples
#' A <- matrix(rnorm(10000), ncol = 1000)
#' summary_fast(A)
#' @export
#' @author Joaquín Martínez-Minaya <\email{jomarminaya@@gmail.com}>
summary_fast <- function(A){
    result <- t(
        rbind(
            Rfast::rowmeans(A),
            Rfast::rowVars(A, std = TRUE),
            #Rfast::rowMinsMaxs(A),
            #Rfast::rowMedians(A),
            t(Rfast2::rowQuantile(A,
                                probs = c(0.025, 0.5, 0.975)))))
    colnames(result) <- c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")
    result
}


