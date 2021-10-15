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
n <- 50
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

V <- cbind(V, iid1, iid2, iid3, iid4)

# Formula that we want to fit
formula <- y ~ -1 + v1 + f(iid1, model = 'iid') |
  -1 + v2 + f(iid1, model = 'iid') |
  -1 + v3 + f(iid2, model = 'iid') |
  -1 + v4 + f(iid2, model = 'iid')
names_cat <- formula_list(formula)

x <- c(-1.5, 2,
       1, -3)

#random effect
prec_w <- c(1, 20)
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
### ----- 3.1. Fitting the model with jags --- ####
cat(paste0("n = ", n, " -----> Fitting using SHORT JAGS \n"))

if(file.exists(paste0("model_jags_", n,".RDS"))){
  model.jags <- readRDS(paste0("model_jags_", n,".RDS"))
  if(n < 1000)
  {
    simulation <- readRDS("simulation4_50-500.RDS")
    t_jags <- simulation[[paste0("n",n)]]$times[1]
  }else{
    simulation <- readRDS("simulation4_1000-10000.RDS")
    t_jags <- simulation[[paste0("n",n)]]$times[1]
  }

}else{
  ## MCMC configuration
  ni <- 5000
  nt <- 5
  nb <- 1000
  nc <- 3

  ## Data set
  table(V$iid2) %>%length() -> niv_length
  data_jags <- list(y = y,
                    N = dim(y)[1],
                    d = d,
                    V = V,
                    niv = V$iid2,
                    niv_length = niv_length)

  ## Initial values
  inits <- function(){list(beta1 = rnorm(d, 0, 1),
                           sd1  = runif(1, 0.1, 10))}

  ## Parameters of interest
  parameters <- c('beta1', 'sd1')

  cat("
    model {
    #model
    for (i in 1:N){
    y[i, ] ~ ddirch(alpha[i,])

    log(alpha[i,1]) <- beta1[1]*V[i,1] + w[niv[i]]
    log(alpha[i,2]) <- beta1[2]*V[i,2] + w[niv[i]]
    log(alpha[i,3]) <- beta1[3]*V[i,3] + w[niv[i]]
    log(alpha[i,4]) <- beta1[4]*V[i,4] + w[niv[i]]
    }

    for(j in 1:niv_length){
        w[j] ~ dnorm(0, tau1)
    }

    #priors
    tau1 <- 1/(sd1*sd1)
    sd1 ~ dunif(0, 10)


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
}





### ----- 3.2. Fitting the model with INLA --- ####
cat(paste0("n = ", n, " -----> Fitting using INLA \n"))
if(file.exists(paste0("model_inla_", n,".RDS"))){
  model.inla <- readRDS(paste0("model_inla_", n,".RDS"))
  if(n< 1000)
  {
    simulation <- readRDS("simulation4_50-500.RDS")
    t_inla <- simulation[[paste0("n",n)]]$times[2]
  }else{
    simulation <- readRDS("simulation4_1000-10000.RDS")
    t_inla <- simulation[[paste0("n",n)]]$times[2]
  }
}else{
  t <- proc.time() # Measure the time
  model.inla <- dirinlareg( formula  = y ~ -1 + v1 + f(iid2, model = 'iid', hyper=list(theta=list(prior="loggamma",param=c(1,0.1)))) |
                              -1 + v2 + f(iid2, model = 'iid') |
                              -1 + v3 + f(iid2, model = 'iid') |
                              -1 + v4 + f(iid2, model = 'iid'),
                            y        = y,
                            data.cov = V,
                            prec     = 0.01,
                            verbose  = TRUE)

  t_inla <- proc.time()-t    # Stop the time
  t_inla <- t_inla[3]
  summary(model.inla)
}

### ----- 3.2. Fitting the model with INLA pc prior--- ####
cat(paste0("n = ", n, " -----> Fitting using INLA \n"))
if(file.exists(paste0("model_inla_", n,".RDS"))){
  model.inla.2 <- readRDS(paste0("model_inla_pc_", n,".RDS"))
  if(n< 1000)
  {
    simulation <- readRDS("simulation4_50-500.RDS")
    t_inla <- simulation[[paste0("n",n)]]$times[2]
  }else{
    simulation <- readRDS("simulation4_1000-10000.RDS")
    t_inla <- simulation[[paste0("n",n)]]$times[2]
  }
}else{
  t <- proc.time() # Measure the time
  model.inla.2 <- dirinlareg( formula  = y ~ -1 + v1 + f(iid1, model = 'iid', hyper=list(theta=(list(prior="pc.prec", param=c(10, 0.01))))) |
                                -1 + v2 + f(iid2, model = 'iid', hyper=list(theta=(list(prior="pc.prec", param=c(10, 0.01))))) |
                                -1 + v3 + f(iid3, model = 'iid', hyper=list(theta=(list(prior="pc.prec", param=c(10, 0.01))))) |
                                -1 + v4 + f(iid4, model = 'iid', hyper=list(theta=(list(prior="pc.prec", param=c(10, 0.01))))),
                              y        = y,
                              data.cov = V,
                              prec     = 0.01,
                              verbose  = TRUE)

  t_inla_2 <- proc.time()-t    # Stop the time
  t_inla_2 <- t_inla_2[3]
  summary(model.inla.2)
}
