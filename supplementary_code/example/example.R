### In this example, we show how to fit a model using the dirinla package ###
### --- 1. Loading the libraries --- ####
library(dirinla)
library(INLA)
library(DirichletReg)
library(ggplot2)
library(gridExtra)


### --- 2. Simulating from a Dirichlet likelihood --- ####
set.seed(1000)
N <- 50 #number of data
V <- as.data.frame(matrix(runif((4) * N, 0, 1), ncol = 4)) #Covariates
names(V) <- paste0('v', 1:4)

formula <- y ~ 1 + v1 | 1 + v2 | 1 + v3 | 1 + v4
(names_cat <- formula_list(formula))

x <- c(-1.5, 1, -3, 1.5,
       2, -3 , -1, 5)

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
  formula  = y ~ 1 + v1 | 1 + v2 | 1 + v3 | 1 + v4,
  y        = y,
  data.cov = V,
  prec     = 0.0001,
  verbose  = TRUE)


summary(model.inla)

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
  p1[[i]] <- p1[[i]] + geom_vline(xintercept = x[seq(1,4, by=1)][i], col = "red")


  p1[[i]] <- p1[[i]] + ggtitle(paste0("Category ", i)) +
    theme(
      plot.title = element_text(color = "black",
                                size  = 12,
                                face  = "bold.italic",
                                hjust = 0.5))
}


grid.arrange(p1[[1]], p1[[2]], p1[[3]], p1[[4]], ncol = 4)

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
  p2[[i]] <- p2[[i]] + geom_vline(xintercept = x[seq(5,8, by=1)][i], col = "red")


  p2[[i]] <- p2[[i]] +
    ggtitle(" ") +
    theme(
      plot.title = element_text(color = "black",
                                size  = 15,
                                face  = "bold.italic",
                                hjust = 0.5))
}


grid.arrange(p2[[1]], p2[[2]], p2[[3]], p2[[4]], ncol = 4)


pdf("intercept_slopes_example.pdf", width=10, height=4)
grid.arrange(p1[[1]], p1[[2]], p1[[3]], p1[[4]],
             p2[[1]], p2[[2]], p2[[3]], p2[[4]],
             ncol = 4)
dev.off()


### --- 5. Predicting for v1 = 0.25, v2 = 0.5, v3 = 0.5, v4 = 0.1 --- ####
model.prediction <-
  predict_dirinla(model.inla,
                  data.pred = data.frame(v1 = 0.25,
                                         v2 = 0.5,
                                         v3 = 0.5,
                                         v4 = 0.1))
model.prediction$summary_predictive_means
