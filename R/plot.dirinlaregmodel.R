#' plot of dirinlaregmodel objects
#'
#' `summary.dirinlaregmodel` Method which plots a dirinlaregmodel object
#'
#' @param object Object of class dirinlareg.
#' @return Plotting the posterior of the fixed effects.
#' @export
#'
#' @import ggplot2
#' @importFrom ggtern ggtern
#' @importFrom gridExtra grid.arrange
#' @author Joaquín Martínez-Minaya <\email{jomarminaya@@gmail.com}>
plot.dirinlaregmodel <- function(object) {
  if(dim(object$y)[2]> 3)
  {
    warning("Dimension is greater than 3 -> Ternary diagram has not been plotted.")
  }else{
    nombres <- names(object$summary_means)
    datos <- as.data.frame(sapply(object$summary_alphas, function(x){x[,"mean"]}))

    #Simulating from response variable
    alpha <- as.matrix(datos)
    y_resp <- as.data.frame(rdirichlet(dim(datos)[1], alpha))
    colnames(y_resp) <- colnames(datos)
    a <- ggtern::ggtern(data = y_resp,
           aes_string( x = nombres[1],
                       y = nombres[2],
                       z = nombres[3])) +
      ggtern::stat_density_tern(geom='polygon',
                        n = 200,
                        aes(fill=..level..,
                            alpha = ..level..),
                        base = "identity") +
      theme_rgbw() +
      guides(color = "none", fill = "none", alpha = "none") +
      geom_point(data = as.data.frame(object$y),
                 aes_string(x = nombres[1],
                            y = nombres[2],
                            z = nombres[3]),
                 size = 0.2) +
      ggtitle("Fitted Density vs Original data") +
      scale_fill_gradient(low='blue',high='red')
    print(a)
  }

  for(x in 1:length(object$marginals_fixed))
  {
    p1 <- list()
    for(i in 1:length(object$marginals_fixed[[x]]))
    {
      dens <- as.data.frame(object$marginals_fixed[[x]][[i]])
      p1[[i]] <- ggplot2::ggplot(dens,
                                 aes(x = x,
                                     y = y)) +
        ggplot2::geom_line(size = 0.6, col = "red4") +
        #xlim(c(min(dens$x[dens$group=="R-JAGS"]), max(dens$x[dens$group=="R-JAGS"]))) +
        ggplot2::theme_light() + #Show axes
        ggplot2::xlab(names(object$marginals_fixed[[x]])[i]) + #xlab
        ggplot2::ylab("f()")
    }
    args <- c(p1, list(ncol = 2, top = nombres[x]))
  do.call(gridExtra::grid.arrange,
          args)
        }
}
