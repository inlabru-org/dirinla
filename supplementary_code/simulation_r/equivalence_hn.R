### Equivalence between half normal from INLA and half normal from jags.

library(INLA)
library(rjags)
library(R2jags)
library(dplyr)
### Plotting Half normal
## INLA implementation in terms of tau
# HN.prior <<- "expression:
#   tau0 = 0.01;
#   sigma = exp(-theta/2);
#   log_dens = log(2) - 0.5 * log(2 * pi) + 0.5 * log(tau0);
#   log_dens = log_dens - 0.5 * tau0 * sigma^2;
#   log_dens = log_dens - log(2) - theta / 2;
#   return(log_dens);
# "
#inla in terms of sd
hn_inla <- function(sigma)
{
  tau0 = 1
  log_dens = log(2) - 0.5 * log(2 * pi) + 0.5 * log(tau0)
  log_dens = log_dens - 0.5 * tau0 * sigma^2
  return(exp(log_dens))
}

sigma <- seq(-10, 10, 0.01)
plot(sigma, hn_inla(sigma))
lines(sigma, dnorm(sigma, sd = 1))

## rjags
## MCMC configuration
nt <- 5
nc <- 3
ni <- 100000
nb <- 5000

## Data set. Fake dataset (1 observation)
y <- 0.0001

data_jags <- list(y = y,
                  N = 1)

## Initial values
inits <- function(){list(sigma = runif(1, 0, 0.5))}

#Data
data_jags <- list(y = y)

parameters <- c('sigma')

cat("
    model {
    #model
    y ~ dnorm(sigma, 0.000000000000001)

    sigma ~ dnorm(0,1)T(0,)
    }", file="model_cov2.jags" )





## Call jags
t <- proc.time() # Measure the time
model.jags.2 <- jags(data_jags,
                     inits,
                     parameters,
                     "model_cov2.jags",
                     n.chains          = nc,
                     n.thin            = nt,
                     n.iter            = ni,
                     n.burnin          = nb,
                     working.directory = getwd()) #
t_jags_2<-proc.time()-t    # Stop the time
t_jags_2 <- t_jags_2[3]



#Comparing both
model.jags.2$BUGSoutput$sims.array[,,2] %>% hist(., freq = FALSE, breaks = 100)
lines(sigma, hn_inla(sigma))

