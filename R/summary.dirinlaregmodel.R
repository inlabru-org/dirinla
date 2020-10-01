#' Summary of dirinlaregmodel objects
#'
#' `summary.dirinlaregmodel` is a function which gives a summary of a dirinlaregmodel object
#'
#' @param object Object of class dirinlaregmodel.
#' @param ... Other arguments.
#' @return Print summary.
#' @method summary dirinlaregmodel
#' @export
#' @author Joaquín Martínez-Minaya <\email{jomarminaya@@gmail.com}>

summary.dirinlaregmodel <- function(object, ...) {
    cat("\n")
    ## Call
    cat("Call: \n ")
    print(object$call)

    cat("\n \n")
    ## Summary of the fixed effects
    cat("---- FIXED EFFECTS ---- \n")

    for (i in 1:length(object$summary_fixed)) {
        cat("======================================================================= \n")
        cat(paste0(names(object$summary_fixed)[i], "\n"))
        #cat(paste0("Category "), i, "\n")
        cat("----------------------------------------------------------------------- \n")
        print(object$summary_fixed[[i]], digits = 4)
    }
    cat("======================================================================= \n")
    cat("\n")

    ### Model selection
    cat(paste0("DIC = ", round(object$dic$dic, 4)), ", ")
    cat(paste0("WAIC = ", round(object$waic$waic, 4)), ", ")
    cat(paste0("LCPO = ", round(-sum(log(object$cpo$cpo)), 4), " \n"))

    ### Observations and categories
    cat(paste0("Number of observations: ", object$nobs, "\n"))
    cat(paste0("Number of Categories: ", object$ncat, "\n"))
}





