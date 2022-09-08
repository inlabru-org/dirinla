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
  if (!safe_inla()) {
    stop(inla_install_info("data_stack_dirich_formula"))
  }

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
    index_mat[index_mat == 0] <- NA

    ## For common effect we have to check it they have the same arguments!
    #Esto falta por hacer


    ### Mixing categories with random effects, index_mat and A_random
    index_random_names %>% lapply(., function(x){
      kronecker(data[,x], index_mat[,x] )
    }) -> effects_random
    names(effects_random) <- index_random_names


    #Formula
    formula.inla.pred2 <- random_eff %>% unlist(.) %>% unique(.) %>% paste(., collapse = "+")
    formula.inla.pred <- paste(formula.inla.pred, formula.inla.pred2, sep = "+")

  }else{
    effects_random <- rep(NA, n*d)
  }

  ### Mixing two formulas
  #formula.inla <- covariatesall %>% names() %>% str_remove(., "1") %>% .[1] %>% paste0(., " ~ -1 + ")
  formula.inla <- "y" %>% paste0(., " ~ -1 + ")

  formula.inla <- as.formula(paste(formula.inla, formula.inla.pred, collapse = " " ))
  effects <- cbind(effects, effects_random)
  if(all(is.na(effects[,dim(effects)[2]])))
  {
    effects <- effects[,-dim(effects)[2]]
  }

  A <- 1
  #Giving back the inla.stack
  list(stk = INLA::inla.stack(data = list(y = y), A = c(A), effects = effects),
       formula.inla = formula.inla)
}
