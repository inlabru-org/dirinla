#' Preparing the data
#'
#' `data_stack_dirich_formula` prepares the data using inla.stack from the package INLA.
#'
#' @param y Response variable in a matrix format.
#' @param covariates String with the name of covariates.
#' @param share Covariates to share in all the cateogries. Not implemented yet.
#' @param data Data.frame which contains all the covariates.
#' @param d Number of categories.
#' @param n Number of locations.
#'
#' @return List with two objects
#' - Object of class inla.stack
#' - Object with class formula
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
data_stack_dirich_formula <- function(y, covariates, share = NULL, data, d, n) {
  data <- cbind(intercept = rep(1, dim(data)[1]), data)


  covariatesall <- covariates

  ### Fixed effect
  covariates %>% lapply(., function(x){
    logic1 <- x %>% str_starts("f\\(") %>% !.
    x[logic1]
  }) -> covariatesall

  if(length(unlist(covariatesall)) >=1){
    ### Prepare covariates
    1:length(covariatesall)  %>%
      lapply(., function(x){
        data_x <- data %>% dplyr::select(covariatesall[[x]])
        categories <- paste0("cat", x,  "_", covariatesall[[x]])
        index <- rep(NA, d)
        index[x] <- 1
        data_x <- kronecker(as.matrix(data_x), index)
        colnames(data_x) <- categories
        data_x
      }) %>%
      do.call(cbind.data.frame, .) -> effects

    ### Matrix A
    A <- 1

    ### Terms for the formula
    names_inla_fixed <- names(effects)
    formula.inla.pred <- character()

    formula.inla.pred <- c(formula.inla.pred, paste0("f(",
                                                     names_inla_fixed,
                                                     ", model = 'linear')"))

    formula.inla.pred <- stringr::str_c(formula.inla.pred, collapse=" + ")
  }




#
#
#
#
#
#   ### Fixed effects
#   data_cov_covariatesall <- lapply(covariatesall, dplyr::select, .data = data)
#
#   A.names <- names(A)
#   for (j in 1:length(data_cov_covariatesall)) {
#     if (length(covariatesall[[j]]) != 0) {
#       pos <- rep(0, d)
#       pos[j] <- 1
#       for (k in 1:dim(data_cov_covariatesall[[j]])[2]) {
#         A <- c(A, Matrix::Matrix(kronecker(
#           data_cov_covariatesall[[j]][, k], pos)))
#         new.effect <- list(1)
#         names(new.effect) <- paste0(
#           "cat", j, "_",
#           names(data_cov_covariatesall[[j]][k]))
#         effects <- c(effects, new.effect)
#         A.names <- c(A.names,
#                      paste0("cat", j, "-",
#                             names(data_cov_covariatesall[[j]][k])))
#       }
#     }
#   }
#   names(A) <- A.names
#
  #############################################################################
  #### If we want to include more steps for random effects here is the part ###
  #############################################################################


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
      #sharing
      effectsiid <- kronecker(data[[random_eff_args[[1]]$term]], rep(1, d)) %>% data.frame(.)
      colnames(effectsiid) <- random_eff_args[[1]]$term
      formula.inla.pred <- paste(formula.inla.pred, random_eff[[1]], sep = "+")
    }
  }else{
    effectsiid <- rep(NA, n*d)
  }

  formula.inla <- covariatesall %>% names() %>% str_remove(., "1") %>% .[1] %>% paste0(., " ~ -1 + ")
  formula.inla <- as.formula(paste(formula.inla, formula.inla.pred, collapse = " " ))
  effects <- cbind(effects, effectsiid)
  if(all(is.na(effects[,dim(effects)[2]])))
  {
    effects <- effects[,-dim(effects)[2]]

  }

  list(inla.stack(data = list(y = y), A = c(A), effects = effects),
       formula.inla = formula.inla)
}
