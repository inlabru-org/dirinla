#' Extracting posterior distributions from the linear predictor
#'
#' `extract_linear_predictor` extracts the posterior distribution from the linear predictor
#'
#' @param inla_model An object of class inla.
#' @param n Number of observations.
#' @param d Number of categories.
#' @param Lk_eta Cholesky decomposition of the Hessian matrix.
#'
#' @return summary_linear_predictor List containing a summary of the marginal posterior distributions of the linear predictor.
#' @return marginals_linear_predictor List containing the marginal posterior distributions of the linear predictor.
#' @return summary_alphas List containing a summary of the marginal posterior distributions of the alphas.
#' @return marginals_alphas List containing the marginal posterior distributions of the alphas.
#' @return summary_precision List containing a summary of the marginal posterior distributions of the precision.
#' @return marginals_precision List containing the marginal posterior distributions of the precision.
#' @return summary_means List containing a summary of the marginal posterior distributions of the means.
#' @return marginals_means List containing the marginal posterior distributions of the means.
#' @return simulations of the linear predictor.
#'
#' @importFrom stats density as.formula sd
#' @author Joaquín Martínez-Minaya <\email{joaquin.martinez-minaya@@uv.es}>
#'
extract_linear_predictor <- function(inla_model, n, d, Lk_eta) {
    ### --- 1. Simulating in order to get posterior distributions of the linear predictor --- ####
    p_mod <- inla.posterior.sample(10000, inla_model)
    p_mod <- Matrix(sapply(p_mod, function(x) x$latent[1:(n * d)]))  #L^t eta
    p_predictor <- solve(t(Lk_eta), p_mod)

    # Computing predictor
    p_predictor <- lapply(1:d, function(j) {
        p_predictor[seq(j, dim(p_predictor)[1], by = d), ]
    })
    names(p_predictor) <- paste0("Category ", 1:d)

    ### --- 2. Linear predictor: summary and marginals --- #### ----- 2.1. Summary of the linear predictor ---
    ### ####
    summary_linear_predictor <- p_predictor %>% map(function(A) {
        as.data.frame(t(apply(A, 1, function(x) summary(x))))
    })

    ### ----- 2.2. Marginals of the linear predictor --- ####
    marginals_linear_predictor <- p_predictor %>% map(function(cat) {
        apply(cat, 1, function(sim) {
            dens <- density(sim, n = 1000)
            data.frame(x = dens$x, y = dens$y)
        })
    })


    ### --- 3. alphas: summary and marginals --- #### ----- 3.1. Simulation of alphas --- ####
    p_alphas <- p_predictor %>% map(exp)

    ### ----- 3.2. Summary of alphas --- ####
    summary_alphas <- p_alphas %>% map(function(A) {
        as.data.frame(t(apply(A, 1, function(x) summary(x))))
    })

    ### ----- 3.3. Marginals of alphas --- ####
    marginals_alphas <- p_alphas %>% map(function(cat) {
        apply(cat, 1, function(sim) {
            dens <- density(sim, n = 75)
            data.frame(x = dens$x, y = dens$y)
        })
    })

    ### ----- 3.4. Computing simulation of precision (sum of alphas) --- ####
    alpha0 <- Reduce("+", p_alphas)

    ### ----- 3.5. Summary of precision (alpha0)
    summary_precision <- as.data.frame(t(apply(alpha0, 1, summary)))

    ### ----- 3.6. Marginals of precision
    marginals_precision <- apply(alpha0, 1, function(sim) {
        dens <- density(sim, n = 75)
        data.frame(x = dens$x, y = dens$y)
    })

    ### --- 4. Means: posterior means --- #### --- 4.1. Simulations --- ####
    p_means <- p_alphas %>% map(function(x) x/alpha0)

    ### ----- 4.2. Summary --- ####
    summary_means <- p_means %>% map(function(A) {
        as.data.frame(t(apply(A, 1, function(x) summary(x))))
    })

    # Reduce('+', summary_means) #Summ up to one


    ### ----- 4.3. Marginals of means --- ####
    marginals_means <- p_means %>% map(function(cat) {
        apply(cat, 1, function(sim) {
            dens <- density(sim, n = 75)
            data.frame(x = dens$x, y = dens$y)
        })
    })




    ### --- 5. Returning parameters --- ####
    list(summary_linear_predictor = summary_linear_predictor, marginals_linear_predictor = marginals_linear_predictor,
        summary_alphas = summary_alphas, marginals_alphas = marginals_alphas, summary_precision = summary_precision,
        marginals_precision = marginals_precision, summary_means = summary_means, marginals_means = marginals_means,
        p_predictor = p_predictor)

}
