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
cat("n = ", n, " -----> Simulating data \n")
set.seed(100)
n <- 100
levels_factor <- n
cat_elem <- n/levels_factor

#Covariates
V <- as.data.frame(matrix(runif((10)*n, -1, 1), ncol=10))
#V <- as.data.frame(matrix(rnorm((10)*n, 0, 1), ncol=10))
names(V) <- paste0('v', 1:(10))

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
prec_w <- c(4, 3)
(sd_w <- 1/sqrt(prec_w))

w1 <- rnorm(levels_factor, sd = sqrt(1/prec_w[1])) %>% rep(., rep(n/levels_factor, levels_factor))
w2 <- w1
w3 <- rnorm(levels_factor, sd = sqrt(1/prec_w[2])) %>% rep(., rep(n/levels_factor, levels_factor))
w4 <- w2


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





###
model.inla.2 <- dirinlareg( formula  = y ~ -1 + v1 + f(iid1, model = 'iid', hyper=list(theta=(list(prior="pc.prec", param=c(10, 0.01))))) |
                              -1 + v2 + f(iid1, model = 'iid', hyper=list(theta=(list(prior="pc.prec", param=c(10, 0.01))))) |
                              -1 + v3 + f(iid2, model = 'iid', hyper=list(theta=(list(prior="pc.prec", param=c(10, 0.01))))) |
                              -1 + v4 + f(iid2, model = 'iid', hyper=list(theta=(list(prior="pc.prec", param=c(10, 0.01))))),
                            y        = y,
                            data.cov = V,
                            prec     = 0.01,
                            verbose  = TRUE)

t_inla_2 <- proc.time()-t    # Stop the time
t_inla_2 <- t_inla_2[3]
summary(model.inla.2)

model.inla.2$summary_hyperpar









### --- 3. Comparing posterior distributions. jags vs INLA --- ####
### ----- 3.1. Fitting the model with jags with pc-prior for sd --- ####
cat(paste0("n = ", n, " -----> Fitting using SHORT JAGS \n"))

x_pc <- seq(0.01,15, length.out = 1000)
y_pc <- inla.pc.dprec(x_pc, u = 10, alpha = 0.01, log = FALSE)
prec_pc <- list(x_pc = x_pc, y_pc = y_pc)


## MCMC configuration
ni <- 5000
nt <- 5
nb <- 1000
nc <- 3

## Data set
table(V$iid1) %>%length() -> niv_length1
table(V$iid2) %>%length() -> niv_length2
data_jags <- list(y = y,
                  N = dim(y)[1],
                  d = d,
                  V = V,
                  niv1 = V$iid1,
                  niv2 = V$iid2,
                  niv_length1 = niv_length1,
                  niv_length2 = niv_length2,
                  x_pc = prec_pc$x_pc,
                  y_pc = prec_pc$y_pc)

## Initial values
inits <- function(){list(beta1 = rnorm(d, 0, 1))}

## Parameters of interest
parameters <- c('beta1', 'tau1', 'tau2')

cat("
    model {
    #model
    for (i in 1:N){
    y[i, ] ~ ddirch(alpha[i,])

    log(alpha[i,1]) <- beta1[1]*V[i,1] + w1[niv1[i]]
    log(alpha[i,2]) <- beta1[2]*V[i,2] + w1[niv1[i]]
    log(alpha[i,3]) <- beta1[3]*V[i,3] + w2[niv2[i]]
    log(alpha[i,4]) <- beta1[4]*V[i,4] + w2[niv2[i]]
    }

    for(j in 1:niv_length1){
        w1[j] ~ dnorm(0, tau1)
    }

     for(j in 1:niv_length2){
        w2[j] ~ dnorm(0, tau2)
    }

    #priors
    tau1 <- x_pc[index1]
    index1 ~ dcat(y_pc)

    #priors
    tau2 <- x_pc[index2]
    index2 ~ dcat(y_pc)

    for (c in 1:d)
    {
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



hist(model.jags$BUGSoutput$sims.list$tau1, breaks = 100, freq = FALSE)
lines(model.inla.2$marginals_hyperpar$`Precision for iid1`)

hist(model.jags$BUGSoutput$sims.list$tau2, breaks = 100, freq = FALSE)
lines(model.inla.2$marginals_hyperpar$`Precision for iid2`)


