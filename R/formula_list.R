#' Formula in to list
#'
#' `formula_list` reads the formula and generates a list with the name of the covariates used in each category
#'
#' @param form Object of class formula.
#' @param y Matrix containing the response variable R^{nxd}, being n number of individuals
#' and d the number of categories.
#'
#' @return A list with the names of the variables used in each category.
#'
#' @examples
#' formula <- y ~ 1 + v1 + v2 | -1 + v1 | 0 + v2
#' formula_list(formula)
#'
#' @export
#' @author Joaquín Martínez-Minaya <\email{jomarminaya@@gmail.com}>

formula_list <- function(form,y = NULL) {


    ### --- 1. Reading the formula --- ####
    oformula <- form

    response_variable <- oformula[[2L]]
    oformula_str <- format(oformula[[3L]]) %>% paste0(collapse="")


    # Split using the symbol | diferenciating between categories
    oformula_str <- strsplit(oformula_str, split = "\\|")[[1]]

    oformula_str <- gsub(" ", "", oformula_str, fixed = TRUE)

    # Clasified variables by category
    names_cat <- sapply(oformula_str, strsplit, split = "\\+")
    if(is.null(y)){
        names(names_cat) <- paste0("category ", 1:length(names_cat))
    }else{
        names(names_cat) <- colnames(y)
    }
    # Checking if there are intercepts
    lapply(names_cat, function(x) {
        if (any(grepl("^1$", x) == TRUE)) {
            x[grep("^1$", x)] <- "intercept"
        } else if (any(grepl("^-1$", x) == TRUE)){
            x <- x[-grep("^-1$", x)]
        }else if(any(grepl("^0$", x) == TRUE)){
            x <- x[-grep("^0$", x)]
        }else{
            x <- c("intercept", x)
        }
        x
    })
}

