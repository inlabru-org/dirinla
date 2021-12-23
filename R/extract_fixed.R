#' Extracting marginals fixed of an inla object
#'
#' `extract_fixed` is a function to extract summary and marginals distribution corresponding to the fixed effects
#'
#' @param inla_model Object of inla class.
#' @param names_cat List generated with extract_formula.
#'
#' @return summary_fixed Summary of fixed effects for each category.
#' @return marginals_fixed Marginals for each parameter estimated.
#'
#' @export
#' @importFrom dplyr select
#' @importFrom plyr ldply
#' @author Joaquín Martínez-Minaya <\email{joaquin.martinez-minaya@@uv.es}>
extract_fixed <- function(inla_model, names_cat) {
    summary_fixed <- list()  #List to store the summary of fixed effects for the different categories
    #names_inla <- names(inla_model$summary.random)  #
    names_inla <- rownames(inla_model$summary.fixed)
    marginals_fixed <- list()

    for (i in 1:length(names_cat)) {
        # Auxiliar variables
        names_cat_ind <- names_cat[[i]]
        names_cat_ind <- names_inla %>% stringr::str_detect(., paste0("cat", i, "_")) %>% names_inla[.]
        names_cov <- NULL
        summary_fixed_i <- data.frame(matrix(ncol =7 ))
        summary_fixed_i <- summary_fixed_i[-1,]
        marginals_fixed_i <- list()

            ### Summary
            summary_fixed_i <- rbind(summary_fixed_i, inla_model$summary.fixed[names_cat_ind,])

            ### Marginals
            marginals_fixed_i <- c(marginals_fixed_i, inla_model$marginals.fixed[names_cat_ind])

            ### Names cov
            names_cov <- names_cat_ind %>% stringr::str_remove(., paste0("cat", i, "_"))



        #     # Look for common variables (not common parameters)
        #     if (length(grep(paste0("^", names_cat_ind[j], "$"), names_inla)) != 0) {
        #         pos <- grep(paste0("^", names_cat_ind[j], "$"), names_inla)
        #         ### Summary
        #         summary_fixed_i <- c(summary_fixed_i, list(inla_model$summary.random[[pos]][i, ]))
        #
        #         ### Marginals
        #         marginals_fixed_i <- c(marginals_fixed_i, list(inla_model$marginals.random[[pos]][[i]]))
        #
        #         ### Names cov
        #         names_cov <- c(names_cov, names_cat_ind[j])
        #
        #     } else {
        #         # Look for not common covariates
        #         pos <- grep(paste0("^", "cat", i, "_", names_cat_ind[j], "$"), names_inla)
        #
        #         ### Summary
        #         summary_fixed_i <- c(summary_fixed_i, inla_model$summary.random[pos])
        #
        #         ### Marginals
        #         marginals_fixed_i <- c(marginals_fixed_i, inla_model$marginals.random[[pos]])
        #
        #         ### Names cov
        #         names_cov <- c(names_cov, names_cat_ind[j])
        #
        #     }
        # }

        ### summary fixed Transform the list in data.frame
        #summary_fixed_i <- plyr::ldply(summary_fixed_i, data.frame)
        # Removing .id , ID, kld summary_fixed_i <- summary_fixed_i[,- c(1,2,9)]
        #summary_fixed_i <- dplyr::select(summary_fixed_i, mean, sd, X0.025quant, X0.5quant, X0.975quant, mode)

        summary_fixed_i <- summary_fixed_i %>%
            dplyr::select(.data$mean, .data$sd,
                          .data$`0.025quant`, .data$`0.5quant`, .data$`0.975quant`,
                          .data$mode)

        # Give the names to the covariates
        rownames(summary_fixed_i) <- names_cov

        colnames(summary_fixed_i) <- c("mean", "sd", "0.025quant", "0.5quant", "0.975quant", "mode")

        # Store all the fixed effects corresponding to the i category in the summary_fixed list()
        summary_fixed <- c(summary_fixed, list(summary_fixed_i))


        ### marginals fixed
        names(marginals_fixed_i) <- names_cov
        marginals_fixed <- c(marginals_fixed, list(marginals_fixed_i))
    }

    names(summary_fixed) <- names(marginals_fixed) <- names(names_cat)

    #names(summary_fixed) <- names(marginals_fixed) <- paste0("category ", 1:length(names_cat))

    # Return
    list(summary_fixed = summary_fixed, marginals_fixed = marginals_fixed)

}


