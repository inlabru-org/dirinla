#' Load INLA safely for examples and tests
#'
#' Loads the INLA package with `requireNamespace("INLA", quietly = TRUE)`, and
#' optionally checks and sets the multicore `num.threads` INLA option.
#'
#' @param multicore logical; if `TRUE`, multiple cores are allowed, and the
#' INLA `num.threads` option is not checked or altered.
#' If `FALSE`, forces `num.threads="1:1"`. Default: NULL, checks
#' if running in testthat or non-interactively, in which case sets
#' `multicore=FALSE`, otherwise `TRUE`.
#' @param quietly logical; if `TRUE`, prints diagnostic messages. Default: FALSE.
#' @return logical; `TRUE` if INLA was loaded safely, otherwise FALSE
#' @details Code adapted from `inlabru::bru_safe_inla`
#' @author Finn Lindgren
#' @keywords internal
#' @export
#' @examples
#' \dontrun{
#' if (dirinla_safe_inla()) {
#'   # Run inla dependent calculations
#' }
#' }
#'
dirinla_safe_inla <- function(multicore = NULL,
                      quietly = FALSE) {
  if (requireNamespace("INLA", quietly = TRUE)) {
    if (is.null(multicore)) {
      multicore <-
        !identical(Sys.getenv("TESTTHAT"), "true") ||
        interactive()
    }
    if (!multicore) {
      n.t <- tryCatch(
        INLA::inla.getOption("num.threads"),
        error = function(e) {
          e
        }
      )
      if (inherits(n.t, "simpleError")) {
        if (!quietly) {
          message("inla.getOption() failed. INLA not installed correctly.")
        }
        return(FALSE)
      }
      if (!quietly) {
        message(paste0("Current num.threads is '", n.t, "'."))
      }
      if (!identical(n.t, "1:1")) {
        if (!quietly) {
          message(paste0(
            "Setting INLA option num.threads to '1:1'.",
            " Previous value '", n.t, "'."
          ))
        }
        INLA::inla.setOption(num.threads = "1:1")
      } else {
        if (!quietly) {
          message("No num.threads change needed.")
        }
      }
    }
    TRUE
  } else {
    if (!quietly) {
      message("INLA not loaded safely.")
    }
    FALSE
  }
}

inla_install_info <- function(fun) {
  paste0("To run '", fun, "', please install the INLA package:\n",
         "install.packages('INLA', repos = 'https://inla.r-inla-download.org/R/testing')")
}

