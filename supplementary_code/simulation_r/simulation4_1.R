### In this script simulations with N=50, 100, 5000, 1000, 10000 are conducted #
# in order to check    ###
# how the package dirinla works. We compare with r-jags with enough amount of #
# simulations to guarantee convergence of the method, and with long r-jags,   #
# which has a large amount of simulations. The model that we try to fit is :  #
# --- Y \sim Dirich(exp(eta_1), exp(eta_2), ... , exp(eta_4))                 #
# --- --- eta_1 = beta_{01} + X1 beta1,                                       #
# --- --- eta_2 = beta_{02} + X2 beta2,                                       #
# --- --- eta_3 = beta_{03} + X3 beta3,                                       #
# --- --- eta_4 = beta_{04} + X4 beta4,                                       #
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


### --- 2. Simulating using a random effect --- ####
n <- 100
cat("n = ", n, " -----> Simulating data \n")
set.seed(1000)

#Covariates
V <- as.data.frame(matrix(runif((10)*n, 0, 1), ncol=10))
names(V) <- paste0('v', 1:(10))
id1 <- 1:dim(V)[1]

# Formula that we want to fit
formula <- y ~ 1 + v1 + f(id1, model = 'iid') | 1 + v2 + f(id1, model = 'iid') | 1 + v3 + f(id1, model = 'iid') | 1 + v4 + f(id1, model = 'iid')
names_cat <- formula_list(formula)

# Parameters to fit
# x <- c(-1.5, 1, -3, 1.5,
#        2, -3 , -1, 5)


x <- c(-1.5, 2,
       1, -3,
       -3, -1,
       1.5, 5)

#random effect
prec_w <- 10
w <- rnorm(n, sd = sqrt(1/prec_w))

x <- c(x, w)

d <- length(names_cat)
data_stack_construct <- data_stack_dirich(y          = as.vector(rep(NA, n*d)),
                                          covariates = names_cat,
                                          share      = NULL,
                                          data       = V,
                                          d          = d,
                                          n          = n )

# Ordering the data with covariates --- ###
A_construct <- data_stack_construct$A
eta <- A_construct %*% x
alpha <- exp(eta)
alpha <- matrix(alpha,
                ncol  = d,
                byrow = TRUE)
y_o <- rdirichlet(n, alpha)
colnames(y_o) <- paste0("y", 1:d)


y <- y_o

### ----- 3.2. Fitting the model with INLA --- ####
cat(paste0("n = ", n, " -----> Fitting using INLA \n"))
t <- proc.time() # Measure the time
model.inla <- dirinlareg( formula  = y ~ 1 + v1 + f(id1, model = 'iid') | 1 + v2 + f(id1, model = 'iid') | 1 + v3 + f(id1, model = 'iid') | 1 + v4 + f(id1, model = 'iid'),
                          y        = y,
                          data.cov = V,
                          prec     = 0.0001,
                          verbose  = TRUE)

t_inla <- proc.time()-t    # Stop the time
t_inla <- t_inla[3]
summary(model.inla)


