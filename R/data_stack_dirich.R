#' Preparing the data
#'
#' `data_stack_dirich` prepares the data using inla.stack from the package INLA.
#'
#' @param y Response variable in a matrix format.
#' @param covariates String with the name of covariates.
#' @param share Covariates to share in all the cateogries. TODO
#' @param data Data.frame which contains all the covariates.
#' @param d Number of categories.
#' @param n Number of locations.
#'
#' @return Matrix A such as eta = A %*%x
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
#' @author Joaquín Martínez-Minaya <\email{jomarminaya@@gmail.com}>

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

    #############################################################################
    #### If we want to include more steps for random effects here is the part ###
    #############################################################################
    #TODO: tHIS is constructed just for one shared random effect
          # - Include more random effects
          # - Include not shared random effects

    ### Random effects
    covariates %>% lapply(., function(x){
        logic1 <- x %>% str_starts("f\\(") %>% x[.]
    }) -> random_eff

    if(any(random_eff %>% sapply(., length) >=1))
    {
        ### Extract arguments from formula. We use INLA
        random_eff_args <- lapply(random_eff, function(x){
            x %>% paste0("INLA::", .) %>% parse(text = .) %>% eval(.)
        })

        ### Check if they are going to share components
        #All terms are equal
        cond1 <- unlist(lapply(random_eff_args, function(x) x$term)) %>%
            sapply(., function(x) x == .[1]) %>%
            all(.)

        #All models are equal
        cond2 <- unlist(lapply(random_eff_args, function(x) x$model)) %>%
            sapply(., function(x) x == .[1]) %>%
            all(.)

        if(cond1 && cond2){
            cat("Shared random effect")
            # different elements for different index value (non-balanced effect)
            data %>% dplyr::select(random_eff_args[[1]]$term) %>%
                dplyr::pull() %>%
                table(.) -> index
            lapply(1:length(index), function(x){
                vec <- rep(0, length(data[[random_eff_args[[1]]$term]]))
                vec[data[[random_eff_args[[1]]$term]] == x] <- 1
                Matrix::Matrix(vec)
            }) %>%
            do.call(cbind, .) -> Biid

            #Categories
            Biid <- Matrix::Matrix(kronecker(Biid, rep(1, d)))
         }
        A <- c(A, Biid)
    }

    A <- do.call(cbind, A)

    A
}
