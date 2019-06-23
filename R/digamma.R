#' Computing the function digamma
#'
#' `digamma_def` is the function digamma appropiate for really small values
#'
#' @param x: Argument to applied the function digamma.
#' @param ...: Rest of arguments used in the case of digamma functions.
#'
#' @return Result of applying digamma function
#'
#' @author Joaquín Martínez-Minaya <\email{joaquin.martinez-minaya@@uv.es}>
digamma_red <- function(x, ...) {
    pos1 <- which(x <= 0.01)
    pos2 <- which(x > 0.01)
    result <- numeric()
    result[pos1] <- digamma(x[pos1] + 1, ...) - 1/x[pos1]
    result[pos2] <- digamma(x[pos2], ...)
    result
}

#' Computing the function trigamma
#'
#' `trigamma_red` is the function trigamma appropiate for really small values
#'
#' @param x: Argument to applied the function trigamma.
#' @param ...: Rest of arguments used in the case of digamma functions.
#'
#' @return Result of applying trigamma function.
#'
#' @author Joaquín Martínez-Minaya <\email{joaquin.martinez-minaya@@uv.es}>

trigamma_red <- function(x, ...) {
    # We use the equatlity trigamma(x+1) = trigamma(x) - 1/x^2
    pos1 <- which(x <= 0.01)
    pos2 <- which(x > 0.01)
    result <- numeric()
    result[pos1] <- trigamma(x[pos1] + 1, ...) + 1/x[pos1]^2
    result[pos2] <- trigamma(x[pos2], ...)
    result
}
