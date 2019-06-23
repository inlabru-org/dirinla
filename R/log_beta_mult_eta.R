#' Calculating the log beta function in eta
#'
#' `beta_mult_eta` computes the log beta function in eta
#'
#' @param x: Vector of elements.
#' @param y: Vector with dimension 1 x d.
#'
#' @return Numeric value.
#'
#' @author Joaquín Martínez-Minaya <\email{joaquin.martinez-minaya@@uv.es}>


log_beta_mult_eta <- function(x) {
    num <- sum(sapply(x, lgamma))
    den <- lgamma(sum(x))
    log_beta_mult <- num - den
    log_beta_mult
}
