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
#' @import dplyr purrr
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
            list1 <- lapply(x, function(x1){
                form1 <- paste0("INLA::", x1) %>% parse(text = .) %>% eval(.)
                form1
                })
            names(list1) <- purrr::map(list1, "term")
            list1
        })


        ### Check if same index is used in different categories.
        ## All the names for index
        index_random_names <- purrr::map(random_eff_args, names) %>% unlist(.) %>% unique(.)

        ## Checking where are common effects
        index_mat <- purrr::map(random_eff_args, names) %>%
            lapply(., function(x){(index_random_names %in% x) %>% as.numeric()}) %>%
            do.call(rbind, .) %>% Matrix(.)
        colnames(index_mat) <- index_random_names

        ## For common effect we have to check it they have the same arguments!
        #Esto falta por hacer

        ## Constructing the matrix for the effects (withouth having in mind categories)
        index_random_names %>% lapply(., function(x){
            dat <- data.frame(x = 1:n,
                              y = data[,x])
            A_random <- Matrix(data = 0, nrow = n, ncol = length(data[,x] %>% table(.)),
                               dimnames = list(as.character(1:n), names(table(data[, x]))))
            A_random[as.matrix(dat)] <- 1
            A_random
        }) -> A_random
        names(A_random) <- index_random_names


        ### Mixing categories with random effects, index_mat and A_random
        index_random_names %>% lapply(., function(x){
            kronecker(A_random[[x]], index_mat[,x] )
        }) -> A_random
        names(A_random) <- index_random_names

        #Joining_matrix
        A <- c(A, A_random)
    }

    A <- do.call(cbind, A)

    A
}
