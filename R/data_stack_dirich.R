#' Preparing the data
#'
#' `data_stack_dirich` prepares the data using inla.stack from the package INLA.
#'
#' @param y Response variable in a matrix format.
#' @param covariates String with the name of covariates.
#' @param share Covariates to share in all the cateogries. Not implemented yet.
#' @param data Data.frame which contains all the covariates.
#' @param d Number of categories.
#' @param n Number of locations.
#'
#' @return Object of class inla.stack
#'
#' @examples
#' n <- 100
#' d <- 4
#'
#' V <- matrix(rnorm(4*n, 0, 1), ncol=4)
#' V <- as.data.frame(V)
#' names(V) <- c('v1', 'v2', 'v3', 'v4' )
#'
#' covariates <- names(V)
#'
#'
#' formula <- y ~ 1 + v1 + v2 | 1 + v1 | 1 + v1
#'
#' names_cat <- formula_list(formula)
#'
#' data_stack_construct <- data_stack_dirich(y          = as.vector(rep(NA, n*d)),
#'                                           covariates = names_cat,
#'                                           share      = NULL,
#'                                           data       = V,
#'                                           d          = d,
#'                                           n          = n )
#'
#' @export
#' @import dplyr
#' @import Matrix
#' @author Joaquín Martínez-Minaya <\email{joaquin.martinez-minaya@@uv.es}>
data_stack_dirich <- function(y, covariates, share = NULL, data, d, n) {
    data <- cbind(intercept = rep(1, dim(data)[1]), data)

    A <- list()
    effects <- list()

    notcommon <- covariates
    covariates
    ### Fixed effect
    covariates %>% lapply(., function(x){
        logic1 <- x %>% str_starts("f\\(") %>% !.
        x[logic1]
    }) -> notcommon




    ### Fixed effects
    data_cov_notcommon <- lapply(notcommon, dplyr::select, .data = data)

    A.names <- names(A)
    for (j in 1:length(data_cov_notcommon)) {
        if (length(notcommon[[j]]) != 0) {
            pos <- rep(0, d)
            pos[j] <- 1
            for (k in 1:dim(data_cov_notcommon[[j]])[2]) {
                A <- c(A, Matrix::Matrix(kronecker(
                  data_cov_notcommon[[j]][, k], pos)))
                new.effect <- list(1)
                names(new.effect) <- paste0(
                  "cat", j, "_",
                  names(data_cov_notcommon[[j]][k]))
                effects <- c(effects, new.effect)
                A.names <- c(A.names,
                             paste0("cat", j, "-",
                                    names(data_cov_notcommon[[j]][k])))
            }
        }
    }
    names(A) <- A.names

    ### Random effects
    covariates %>% lapply(., function(x){
        logic1 <- x %>% str_starts("f\\(") %>% x[.]
    }) -> random_eff

    if(any(random_eff %>% sapply(., length) >=1))
    {
        ### iid
        B <- diag(1, dim(data)[1])

        #Not sharing
        Biid <- list()
        for (j in 1:length(data_cov_notcommon)) {
            pos <- rep(0, d)
            pos[j] <- 1
            Biid[[j]] <- Matrix::Matrix(kronecker(B, pos))
        }
        effectsiid <- lapply(1:d, function(x){1:dim(data)[1]})
        names(effectsiid) <- paste0("id", 1:d)

        #sharing
        Biid <- Matrix::Matrix(kronecker(B, rep(1, d)))
        effectsiid <- list(id1 = 1:dim(data)[1])
    }




    inla.stack(data = list(y = y), A = c(A, Biid), effects = c(effects, effectsiid))
}
