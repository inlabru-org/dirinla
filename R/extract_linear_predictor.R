#' Extracting posterior distributions from the linear predictor
#'
#' `extract_linear_predictor` extracts the posterior distribution from the linear predictor
#'
#' @param inla_model An object of class inla.
#' @param n Number of observations.
#' @param d Number of categories.
#' @param Lk_eta Cholesky decomposition of the Hessian matrix.
#' @param names_cat List generated with extract_formula.
#' @param sim simulations for the function inla.posterior.sample
#' @param verbose if TRUE all the computing process is shown. Default is FALSE
#' @param cores number of cores to be used in the computations
#'
#'
#' @return summary_linear_predictor List containing a summary of the marginal posterior distributions of the linear predictor.
#' @return marginals_linear_predictor List containing simulations of marginal posterior distributions of the linear predictor.
#' @return summary_alphas List containing a summary of the marginal posterior distributions of the alphas.
#' @return marginals_alphas List containing simulations of the marginal posterior distributions of the alphas.
#' @return summary_precision List containing a summary of the marginal posterior distributions of the precision.
#' @return marginals_precision List containing simulations of the marginal posterior distributions of the precision.
#' @return summary_means List containing a summary of the marginal posterior distributions of the means.
#' @return marginals_means List containing the simulations of the marginal posterior distributions of the means.
#'
#' @importFrom stats density as.formula sd
#' @import dplyr
#' @export
#' @author Joaquín Martínez-Minaya <\email{jomarminaya@@gmail.com}>
extract_linear_predictor <- function(inla_model, n, d, Lk_eta, names_cat = names_cat,
                                     sim, verbose, cores) {
    ### --- 1. Simulating in order to get posterior distributions of the linear predictor --- ####
    names_list <- list(1:(n*d))
    names(names_list) <- c("APredictor")
    p_mod <- INLA::inla.posterior.sample(sim, inla_model,
                                   selection = names_list,
                                   verbose = verbose,
                                   num.threads = cores)


    #a <- proc.time()
    p_mod <- Matrix(sapply(p_mod, "[[", "latent"))
    #p_mod <- Matrix(sapply(p_mod, function(x) x$latent[1:(n * d)]))

    #b <- proc.time() - a

    # a <- proc.time()
    # p_mod1 <- inla.posterior.sample.eval(function(...){APredictor},
    #                                      p_mod, return.matrix = TRUE)
    # b <- proc.time() - a

    # a <- proc.time()
     #p_mod2 <- Matrix(sapply(p_mod, function(x) x$latent[1:(n * d)]))  #L^t eta
    # b2 <- proc.time() - a

    p_predictor <- Matrix::solve(t(Lk_eta), p_mod, sparse = TRUE)
    rm(p_mod)
    #p_predictor2 <- Matrix::solve(t(Lk_eta), p_mod, sparse = FALSE)

    p_predictor <- as.matrix(p_predictor)

    # Computing predictor
    p_predictor <- lapply(1:d, function(j) {
        p_predictor[seq(j, dim(p_predictor)[1], by = d), ]
    })

    names(p_predictor) <- names(names_cat)

    ### --- 2. Linear predictor: summary and marginals ---  ###
    #### ----- 2.1. Summary of the linear predictor --- ###
    ####
    summary_linear_predictor <- p_predictor %>% purrr::map(function(A) {
        summary_fast(A)
    })

    ### --- 3. alphas: summary and marginals --- ##
    ## ----- 3.1. Simulation of alphas --- ####
    a <- proc.time()
    p_alphas <- p_predictor %>% purrr::map(exp)

    ### ----- 3.2. Summary of alphas --- ####
    summary_alphas <- p_alphas %>% purrr::map(function(A) {
        summary_fast(A)
    })

    ### ----- 3.4. Computing simulation of precision (sum of alphas) --- ####
    alpha0 <- Reduce("+", p_alphas)

    ### ----- 3.5. Summary of precision (alpha0)
    summary_precision <- summary_fast(alpha0)

    ### --- 4. Means: posterior means --- #### --- 4.1. Simulations --- ####
    p_means <- p_alphas %>% purrr::map(function(x) x/alpha0)

    ### ----- 4.2. Summary --- ####
    summary_means <- p_means %>% purrr::map(function(A) {
        summary_fast(A)
    })

    # Reduce('+', summary_means) #Summ up to one




    ### --- 5. Returning parameters --- ####
    list(summary_linear_predictor   = summary_linear_predictor,
        marginals_linear_predictor = p_predictor,
        summary_alphas              = summary_alphas,
        marginals_alphas            = p_alphas,
        summary_precision           = summary_precision,
        marginals_precision         = alpha0,
        summary_means               = summary_means,
        marginals_means             = p_means)

}
