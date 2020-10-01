### In this example, we show how to fit a model using the dirinla package ###
### --- 1. Loading the libraries --- ####
library(dirinla)
library(INLA)
library(DirichletReg)
library(ggplot2)
library(gridExtra)


### --- 2. Simulating from a Dirichlet likelihood --- ####
set.seed(1000)
N <- 100 #number of data
V <- as.data.frame(matrix(runif((4) * N, 0, 1), ncol = 4)) #Covariates
names(V) <- paste0('v', 1:4)

formula <- y ~ 1 + v1 | 1 + v2 | 1 + v3
(names_cat <- formula_list(formula))

intercepts <-
x <- c(-1.5, 1, #Cat 1
       -2, 2.3, #Cat 2
       0 , -1.9) #Cat 3

mus <- exp(x) / sum(exp(x))
C <- length(names_cat)
data_stack_construct <-
  data_stack_dirich(y = as.vector(rep(NA, N * C)),
                    covariates = names_cat,
                    data       = V,
                    d          = C,
                    n          = N)

A_construct <- data_stack_construct$A
A_construct[1:8, ]

eta <- A_construct %*% x
alpha <- exp(eta)
alpha <- matrix(alpha,
                ncol  = C,
                byrow = TRUE)
y_o <- rdirichlet(N, alpha)
colnames(y_o) <- paste0("y", 1:C)
head(y_o)


### --- 3. Fitting the model --- ####
y <- y_o
model.inla <- dirinlareg(
  formula  = y ~ 1 + v1 | 1 + v2 | 1 + v3,
  y        = y,
  data.cov = V,
  prec     = 0.0001,
  verbose  = FALSE)


summary(model.inla)
plot(model.inla)



model.inla$summary_fixed

### --- 4. Plotting marginal posterior distributions of the parameters --- ####
### ----- 4.1. Marginal posterior distributions of the intercepts --- ####
p1 <- list()
beta0 <- expression(paste("p(", beta[0], "|", "y)"))

for (i in 1:length(model.inla$marginals_fixed))
{

  #Data combining jags (1) and inla (2)
  dens <- as.data.frame(model.inla$marginals_fixed[[i]][[1]])

  ### Intercept
  p1[[i]] <- ggplot(dens,
                    aes(x = x,
                        y = y
                    )) +
    geom_line(size = 0.6) +
    xlim(quantile(dens$x, probs=c(0.05, 0.95))) +
    theme_bw() + #Show axes
    xlab(expression(beta[0])) + #xlab
    ylab(beta0) #ylab


  #Real value
  p1[[i]] <- p1[[i]] + geom_vline(xintercept = x[seq(1,6, by=2)][i], col = "red")


  p1[[i]] <- p1[[i]] + ggtitle(colnames(y)[i]) +
    theme(
      plot.title = element_text(color = "black",
                                size  = 12,
                                face  = "bold.italic",
                                hjust = 0.5))
}


grid.arrange(p1[[1]], p1[[2]], p1[[3]], ncol = 3)

### ----- 4.2. Marginal posterior distributions of the slopes --- ####
p2 <- list()
beta1 <- expression(paste("p(", beta[1], "|", "y)"))

for (i in 1:length(model.inla$marginals_fixed))
{

  #Data combining jags (1) and inla (2)
  dens <- as.data.frame(model.inla$marginals_fixed[[i]][[2]])

  ### Intercept
  p2[[i]] <- ggplot(dens,
                    aes(x = x,
                        y = y
                    )) +
    geom_line(size = 0.6) +
    xlim(quantile(dens$x, probs=c(0.05, 0.95))) +
    theme_bw() + #Show axes
    xlab(expression(beta[1])) + #xlab
    ylab(beta1) #ylab


  #Real value
  p2[[i]] <- p2[[i]] + geom_vline(xintercept = x[seq(2,6, by=2)][i], col = "red")


  p2[[i]] <- p2[[i]] +
    ggtitle(" ") +
    theme(
      plot.title = element_text(color = "black",
                                size  = 15,
                                face  = "bold.italic",
                                hjust = 0.5))
}


grid.arrange(p2[[1]], p2[[2]], p2[[3]], ncol = 3)


pdf("intercept_slopes_example.pdf", width=8, height=4)
grid.arrange(p1[[1]], p1[[2]], p1[[3]],
             p2[[1]], p2[[2]], p2[[3]],
             ncol = 3)
dev.off()


### --- 5. Posterior predictive density for $y$ in the simplex
set.seed(4)
nombres <- names(model.inla$summary_means)
datos <- as.data.frame(sapply(model.inla$summary_alphas, function(x){x[,"mean"]}))

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
  ggtern::theme_rgbw() +
  guides(color = "none", fill = "none", alpha = "none") +
  geom_point(data = as.data.frame(model.inla$y),
             aes_string(x = nombres[1],
                        y = nombres[2],
                        z = nombres[3]),
             size = 0.2) +
  #ggtitle("Fitted Density vs Original data") +
  scale_fill_gradient(low='blue',high='red')
pdf("example_tern.pdf", width=5, height=5)
 print(a)
dev.off()

### --- 5. Predicting for v1 = 0.25, v2 = 0.5, v3 = 0.5, v4 = 0.1 --- ####
model.prediction <-
  predict(model.inla,
                  data.pred.cov = data.frame(v1 = 0.2 ,
                                         v2 = 0.5,
                                         v3 = -0.1))
model.prediction$summary_predictive_means

### --- 6. We can also predict directly --- ####
model.inla <- dirinlareg(
  formula  = y ~ 1 + v1 | 1 + v2 | 1 + v3,
  y        = y,
  data.cov = V,
  prec     = 0.0001,
  verbose  = TRUE,
  prediction = TRUE,
  data.pred.cov = data.frame(v1 = 0.2 ,
                             v2 = 0.5,
                             v3 = -0.1))


model.prediction$summary_predictive_means

