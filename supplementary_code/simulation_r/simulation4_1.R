### In this script simulations with N=50, 100, 5000, 1000, 10000 are conducted #
# in order to check    ###
# how the package dirinla works. We compare with r-jags with enough amount of #
# simulations to guarantee convergence of the method, and with long r-jags,   #
# which has a large amount of simulations. The model that we try to fit is :  #
# --- Y \sim Dirich(exp(eta_1), exp(eta_2), ... , exp(eta_4))                 #
# --- --- eta_1 = beta_{01} + X1 beta1 + f(iid, model = "iid"),               #
# --- --- eta_2 = beta_{02} + X2 beta2 + f(iid, model = "iid"),               #
# --- --- eta_3 = beta_{03} + X3 beta3 + f(iid, model = "iid"),               #
# --- --- eta_4 = beta_{04} + X4 beta4 + f(iid, model = "iid"),               #
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
n <- 50
cat("n = ", n, " -----> Simulating data \n")
set.seed(100)
cat_elem <- 25

#Covariates
V <- as.data.frame(matrix(runif((10)*n, 0, 1), ncol=10))
names(V) <- paste0('v', 1:(10))
#iid2 <- 1:n

iid2 <- rep(1:(n/cat_elem), rep(cat_elem, n/cat_elem))
V <- cbind(V, iid2)
#Desordenamos para el índice
#V <- V[sample(1:dim(V)[1]),]
# V$iid2
# Formula that we want to fit
formula <- y ~ -1 + v1 + f(iid2, model = 'iid') | -1 + v2 + f(iid2, model = 'iid') | -1 + v3 + f(iid2, model = 'iid') | -1 + v4 + f(iid2, model = 'iid')
names_cat <- formula_list(formula)

# formula <- y ~ -1 + v1 + f(iid2, model = 'iid') | -1 + v2 + f(iid2, model = 'iid') | -1 + v3 + f(iid2, model = 'iid') | -1 + v4 + f(iid2, model = 'iid')
# names_cat <- formula_list(formula)

# formula <- y ~ 1 + v1  | 1 + v2  | 1 + v3  | 1 + v4
# names_cat <- formula_list(formula)


# Parameters to fit
# x <- c(-1.5, 1, -3, 1.5,
#        2, -3 , -1, 5)


x <- c(-1.5, 2,
       1, -3)
      # -3, -1,
      # 1.5, 5)


# x <- c(-1.5, 2,
#        1, -3)

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
model.inla <- dirinlareg( formula  = y ~ -1 + v1 + f(iid2, model = 'iid', hyper = list(theta=list(prior="loggamma",param=c(1,0.01)))) |
                             # hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.01))))
                              -1 + v2 + f(iid2, model = 'iid') |
                              -1 + v3 + f(iid2, model = 'iid') |
                              -1 + v4 + f(iid2, model = 'iid'),
                          y        = y,
                          data.cov = V,
                          prec     = 0.01,
                          verbose  = TRUE)

# model.inla <- dirinlareg( formula  = y ~ 1 + v1  | 1 + v2  | 1 + v3  | 1 + v4,
#                           y        = y,
#                           data.cov = V,
#                           prec     = 0.01,
#                           verbose  = TRUE)


t_inla <- proc.time()-t    # Stop the time
t_inla <- t_inla[3]
summary(model.inla)
model.inla$summary_hyperpar
plot(model.inla$marginals_hyperpar$`Precision for iid2`, col = "blue", xlim = c(0,100), type = "l")
abline(v= prec_w)

### --- 4. Fitting the model using JAGS --- ####
ni <- 3000
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
                         tau1  = runif(1, 0.1, 20))}

## Parameters of interest
parameters <- c('beta1', 'tau1')

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
    tau1 ~ dgamma(1, 0.01) ## convert to precision


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



### --- 5. Sensitivity analysis --- ####

### ----- 5.1. Fitting 4 models --- ####
model.inla <- dirinlareg( formula  = y ~ -1 + v1 + f(iid2, model = 'iid', hyper = list(theta=list(prior="loggamma",param=c(1,1)))) |
                              -1 + v2 + f(iid2, model = 'iid') |
                              -1 + v3 + f(iid2, model = 'iid') |
                              -1 + v4 + f(iid2, model = 'iid'),
                          y        = y,
                          data.cov = V,
                          prec     = 0.01,
                          verbose  = TRUE)

model.inla2 <- dirinlareg( formula  = y ~ -1 + v1 + f(iid2, model = 'iid', hyper = list(theta=list(prior="loggamma",param=c(1,0.1)))) |
                              -1 + v2 + f(iid2, model = 'iid') |
                              -1 + v3 + f(iid2, model = 'iid') |
                              -1 + v4 + f(iid2, model = 'iid'),
                          y        = y,
                          data.cov = V,
                          prec     = 0.01,
                          verbose  = TRUE)


model.inla3 <- dirinlareg( formula  = y ~ -1 + v1 + f(iid2, model = 'iid', hyper = list(theta=list(prior="loggamma",param=c(1,0.01)))) |
                              -1 + v2 + f(iid2, model = 'iid') |
                              -1 + v3 + f(iid2, model = 'iid') |
                              -1 + v4 + f(iid2, model = 'iid'),
                          y        = y,
                          data.cov = V,
                          prec     = 0.01,
                          verbose  = TRUE)


model.inla4 <- dirinlareg( formula  = y ~ -1 + v1 + f(iid2, model = 'iid', hyper = list(theta=list(prior="loggamma",param=c(1,0.001)))) |
                               -1 + v2 + f(iid2, model = 'iid') |
                               -1 + v3 + f(iid2, model = 'iid') |
                               -1 + v4 + f(iid2, model = 'iid'),
                           y        = y,
                           data.cov = V,
                           prec     = 0.01,
                           verbose  = TRUE)

### ----- 5.2. Priors used --- ####
x <- seq(0, 100, 0.01)
prior1 <- dgamma(x, 1, 1)
prior2 <- dgamma(x, 1, 0.1)
prior3 <- dgamma(x, 1, 0.01)
prior4 <- dgamma(x, 1, 0.001)
plot(x, prior1)


### ----- 5.3. Ploting priors against posteriors --- ####
pdf("priors_vs_posteriors_tau_slope_100.pdf", width = 10, height = 8)
par(mfrow = c(2,2))
plot(x, prior1, type = "l", col = "red", xlab = expression(tau), main= "Gamma(1,1)")
lines(model.inla$marginals_hyperpar$`Precision for iid2`, col = "blue")
abline(v = prec_w)
legend(60, 0.8, legend=c("Prior", "Posterior", "Real Value"),
       col=c("red", "blue", "black"), lty=1, cex=0.8)


plot(x, prior2, type = "l", col = "red", xlab = expression(tau),
     #ylim = c(0,0.1),
     ylim = c(0,0.2),
     main= "Gamma(1,0.1)")
lines(model.inla2$marginals_hyperpar$`Precision for iid2`, col = "blue")
abline(v = prec_w)

plot(x, prior3, type = "l",
     #ylim = c(0, 0.01),
     ylim = c(0,0.2),
     col = "red", xlab = expression(tau),  main= "Gamma(1,0.01)")
lines(model.inla3$marginals_hyperpar$`Precision for iid2`, col = "blue")
abline(v = prec_w)

plot(x, prior4, type = "l",
     #ylim = c(0, 0.002),
     ylim = c(0,0.2),
     col = "red", xlab = expression(tau),  main= "Gamma(1,0.001)")
lines(model.inla4$marginals_hyperpar$`Precision for iid2`, col = "blue")
abline(v = prec_w)

dev.off()

### --- 6. Some plots --- ####
t_jags
t_inla
par(mfrow=c(1,1))
plot(model.inla$marginals_hyperpar$`Precision for iid2`, col = "blue", xlim = c(0,40), type = "l")
hist(model.jags$BUGSoutput$sims.matrix[, c("tau1")], breaks = 500, freq = FALSE, xlim = c(0,40), add = TRUE)
hist(model.jags$BUGSoutput$sims.matrix[, c("tau1")], breaks = 500, freq = FALSE, xlim = c(0,40))

abline(v = 10, lwd = 3, col = "blue")



hist(model.jags$BUGSoutput$sims.matrix[,c("beta1[1]")],  freq = FALSE, breaks = 100)
lines(model.inla$marginals_fixed$y1$v1, col = "red", xlim = c(0,40), type = "l")

hist(model.jags$BUGSoutput$sims.matrix[,c("beta0[3]")],  freq = FALSE, breaks = 30)
lines(model.inla$marginals_fixed$y3$intercept, col = "red", xlim = c(0,40), type = "l")

model.inla3 <- model.inla
plot(model.inla2$marginals_hyperpar$`Precision for iid2`, col = "red", xlim = c(0,40), type = "l")
lines(model.inla3$marginals_hyperpar$`Precision for iid2`, col = "blue", xlim = c(0,40), type = "l")

