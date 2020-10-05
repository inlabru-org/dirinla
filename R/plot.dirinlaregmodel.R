#' plot of dirinlaregmodel xs
#'
#' `plot.dirinlaregmodel` Method which plots a dirinlaregmodel x
#'
#' @param x Object of class dirinlaregmodel.
#' @param ... Other arguments.
#' @return Plotting the posterior of the fixed effects.
#' @export
#'
#' @import ggplot2
#' @importFrom ggtern ggtern
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices devAskNewPage
#' @author Joaquín Martínez-Minaya <\email{jomarminaya@@gmail.com}>
plot.dirinlaregmodel <- function(x, ...) {
  nombres <- names(x$summary_means)
  if(dim(x$y)[2]> 3)
  {
    warning("Dimension is greater than 3 -> Ternary diagram has not been plotted.")
  }else{
    datos <- as.data.frame(sapply(x$summary_alphas, function(x){x[,"mean"]}))

    #Simulating from response variable
    alpha <- as.matrix(datos)
    y_resp <- as.data.frame(DirichletReg::rdirichlet(dim(datos)[1], alpha))
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
      ggtern::theme_rgbw() +
      guides(color = "none", fill = "none", alpha = "none") +
      geom_point(data = as.data.frame(x$y),
                 aes_string(x = nombres[1],
                            y = nombres[2],
                            z = nombres[3]),
                 size = 0.2) +
      ggtitle("Fitted Density vs Original data") +
      scale_fill_gradient(low='blue',high='red')
    print(a)
  }
  devAskNewPage(ask=TRUE)

  for(j in 1:length(x$marginals_fixed))
  {
    p1 <- list()
    for(i in 1:length(x$marginals_fixed[[j]]))
    {
      dens <- as.data.frame(x$marginals_fixed[[j]][[i]])
      p1[[i]] <- ggplot2::ggplot(dens,
                                 aes(x = x,
                                     y = dens$y)) +
        ggplot2::geom_line(size = 0.6, col = "red4") +
        #xlim(c(min(dens$x[dens$group=="R-JAGS"]), max(dens$x[dens$group=="R-JAGS"]))) +
        ggplot2::theme_light() + #Show axes
        ggplot2::xlab(names(x$marginals_fixed[[j]])[i]) + #xlab
        ggplot2::ylab("f()")
    }
    args <- c(p1, list(ncol = 2, top = nombres[j]))
  do.call(gridExtra::grid.arrange,
          args)
  devAskNewPage(ask=TRUE)
        }
}
