### In this script simulations with N=50, 100, 5000, 1000, 10000 are conducted #
# in order to check    ###
# how the package dirinla works. We compare with r-jags with enough amount of #
# simulations to guarantee convergence of the method, and with long r-jags,   #
# which has a large amount of simulations. The model that we try to fit is :  #
# --- Y \sim Dirich(exp(eta_1), exp(eta_2), ... , exp(eta_4))                 #
# --- --- eta_1 = X1 beta1 + f(iid1, model = "iid"),               #
# --- --- eta_2 = X2 beta2 + f(iid2, model = "iid"),               #
# --- --- eta_3 = X3 beta3 + f(iid3, model = "iid"),               #
# --- --- eta_4 = X4 beta4 + f(iid4, model = "iid"),               #
# ----------------------------------------------------------------------------#

### --- 1. Libraries ---- #####
### Needed
library(dirinla)
library(INLA)

#Option
library(DirichletReg) #Dirichlet
library(ggplot2)

#Checking with mcmc
library(rjags)
library(R2jags)

#Generating latex table
library(xtable)

library(dplyr)


### --- 2. Simulation data --- ####
cat("n = ", n, " - levels_factor = ", levels_factor," -----> Simulating data \n")

set.seed(100)
if(is.na(levels_factor)){
  levels_factor <- n
}
cat_elem <- n/levels_factor
cat(paste0(n, "-", levels_factor, "\n"))
#Covariates
V <- as.data.frame(matrix(runif((10)*n, -1, 1), ncol=10))
names(V) <- paste0('v', 1:(10))
tau0 <- 1

#Same covariate
# V <- as.data.frame(cbind(V[,1], V[,1], V[,1], V[,1]))
# #V <- as.data.frame(matrix(rnorm((10)*n, 0, 1), ncol=10))
# names(V) <- paste0('v', 1:(4))

### 4 random effects
iid1 <- iid2  <- rep(1:levels_factor, rep(n/levels_factor, levels_factor))
#Desorder index 3
# pos <- sample(1:length(iid3))
# iid3 <- iid3[pos]

V <- cbind(V, iid1, iid2)

# Formula that we want to fit
formula <- y ~ -1 + v1 + f(iid1, model = 'iid') |
  -1 + v2 + f(iid1, model = 'iid') |
  -1 + v3 + f(iid2, model = 'iid') |
  -1 + v4 + f(iid2, model = 'iid')
names_cat <- formula_list(formula)

x <- c(-1.5, 2,
       1, -3)

#random effect
prec_w <- c(4, 9)
(sd_w <- 1/sqrt(prec_w))

w1 <- rnorm(levels_factor, sd = sqrt(1/prec_w[1])) %>% rep(., rep(n/levels_factor, levels_factor))
w2 <- w1
w3 <- rnorm(levels_factor, sd = sqrt(1/prec_w[2])) %>% rep(., rep(n/levels_factor, levels_factor))
w4 <- w3


#w3 <- w3[pos]
x <- c(x, c(unique(w1),
            unique(w3)))


d <- length(names_cat)
A_construct <- data_stack_dirich(y          = as.vector(rep(NA, n*d)),
                                 covariates = names_cat,
                                 share      = NULL,
                                 data       = V,
                                 d          = d,
                                 n          = n )

# Ordering the data with covariates --- ###
eta <- A_construct %*% x
alpha <- exp(eta)
alpha <- matrix(alpha,
                ncol  = d,
                byrow = TRUE)
y_o <- rdirichlet(n, alpha)
colnames(y_o) <- paste0("y", 1:d)


y <- y_o


### --- 3. Fitting the model --- ####
formula  = y ~
  -1 + v1 + f(iid1, model = 'iid', hyper=list(theta=(list(prior="pc.prec", param=c(10, 0.01))))) |
  -1 + v2 + f(iid1, model = 'iid', hyper=list(theta=(list(prior="pc.prec", param=c(10, 0.01))))) |
  -1 + v3 + f(iid2, model = 'iid', hyper=list(theta=(list(prior="pc.prec", param=c(10, 0.01))))) |
  -1 + v4 + f(iid2, model = 'iid', hyper=list(theta=(list(prior="pc.prec", param=c(10, 0.01)))))
# formula  = y ~
#   -1 + v1 + f(iid1, model = 'iid', hyper=list(theta=list(prior="loggamma",param=c(1,0.1)))) |
#   -1 + v2 + f(iid1, model = 'iid', hyper=list(theta=list(prior="loggamma",param=c(1,0.1)))) |
#   -1 + v3 + f(iid2, model = 'iid', hyper=list(theta=list(prior="loggamma",param=c(1,0.1)))) |
#   -1 + v4 + f(iid2, model = 'iid', hyper=list(theta=list(prior="loggamma",param=c(1,0.1))))
model.inla <- dirinlareg(formula = formula,
                         y        = y,
                         data.cov = V,
                         prec     = 0.01,
                         verbose  = FALSE,
                         control.inla = list(strategy = "laplace", int.strategy = "grid")
)













### --- 4. Summary of the model --- ####
summary(model.inla)
plot(model.inla)
