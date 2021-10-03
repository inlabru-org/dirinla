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

library(dplyr)

### --- 2. Simulating using a random effect --- ####
n <- 500
cat("n = ", n, " -----> Simulating data \n")
set.seed(100)
cat_elem <- 5

#Covariates
V <- as.data.frame(matrix(runif((10)*n, 0, 1), ncol=10))
names(V) <- paste0('v', 1:(10))
#iid1 <- 1:n

iid1 <- rep(1:(n/cat_elem), rep(cat_elem, n/cat_elem))
V <- cbind(V, iid1)
# Formula that we want to fit
formula <- y ~ 1 + v1 + f(iid1, model = 'iid') | 1 + v2 + f(iid1, model = 'iid') | 1 + v3 + f(iid1, model = 'iid') | 1 + v4 + f(iid1, model = 'iid')
names_cat <- formula_list(formula)

# formula <- y ~ 1 + v1  | 1 + v2  | 1 + v3  | 1 + v4
# names_cat <- formula_list(formula)


# Parameters to fit
# x <- c(-1.5, 1, -3, 1.5,
#        2, -3 , -1, 5)


x <- c(-1.5, 2,
       1, -3,
       -3, -1,
       1.5, 5)

#random effect
prec_w <- 10
w <- rnorm(n/cat_elem, sd = sqrt(1/prec_w)) %>% rep(., rep(cat_elem, n/cat_elem))
#w <- rnorm(n, sd = sqrt(1/prec_w))
unique(w)

x <- c(x, unique(w))

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
#y <- DirichletReg::DR_data(y)
summary(y)



### --- 3. Fitting the model with INLA --- ####
cat(paste0("n = ", n, " -----> Fitting using INLA \n"))
t <- proc.time() # Measure the time
model.inla <- dirinlareg( formula  = y ~ 1 + v1 + f(iid1, model = 'iid') | 1 + v2 + f(iid1, model = 'iid') | 1 + v3 + f(iid1, model = 'iid') | 1 + v4 + f(iid1, model = 'iid'),
                          y        = y,
                          data.cov = V,
                          prec     = 0.01,
                          verbose  = TRUE)


t_inla <- proc.time()-t    # Stop the time
t_inla <- t_inla[3]
summary(model.inla)
model.inla$summary_hyperpar


### --- 4. Fitting the model using JAGS --- ####
ni <- 2000
nt <- 5
nb <- 200
nc <- 3

## Data set
data_jags <- list(y = y,
                  N = dim(y)[1],
                  d = d,
                  V = V)

## Initial values
inits <- function(){list(beta0 = rnorm(d, 0, 1),
                         beta1 = rnorm(d, 0, 1),
                         tau1  = runif(1, 0.1, 20))}

## Parameters of interest
parameters <- c('beta0', 'beta1', 'tau1')

cat("
    model {
    #model
    for (i in 1:N){
    y[i, ] ~ ddirch(alpha[i,])

    log(alpha[i,1]) <- beta0[1] + beta1[1]*V[i,1] + w[i]
    log(alpha[i,2]) <- beta0[2] + beta1[2]*V[i,2] + w[i]
    log(alpha[i,3]) <- beta0[3] + beta1[3]*V[i,3] + w[i]
    log(alpha[i,4]) <- beta0[4] + beta1[4]*V[i,4] + w[i]
    }

    for(j in 1:N){
        w[j] ~ dnorm(0, tau1)
    }

    #priors
    tau1 ~ dunif(0.1,100) ## convert to precision


    for (c in 1:d)
    {
    beta0[c]  ~ dnorm(0, 0.01)
    beta1[c] ~ dnorm(0, 0.01)
    }
    }", file="model_cov.jags" )



## Call jags
t <- proc.time() # Measure the time
model.jags <- jags(data_jags,
                   inits,
                   parameters,
                   "model_cov.jags",
                   n.chains          = nc,
                   n.thin            = nt,
                   n.iter            = ni,
                   n.burnin          = nb,
                   working.directory = getwd()) #
t_jags <- proc.time()-t    # Stop the time
t_jags <- t_jags[3]
print(model.jags)

hist(model.jags$BUGSoutput$sims.matrix[, c("tau1")])

