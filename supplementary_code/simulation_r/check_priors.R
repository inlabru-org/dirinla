library(R2jags)


y <- rnorm(50, sd = 0.1)
hist(y)

#Model1
formula  = y ~ 1

model.inla1 <- inla(formula = formula,
                    data   = data.frame(y),
                   control.family = list(hyper = list(
                     theta=(list(prior="pc.prec", param=c(500, 0.01))))))

#Model 2
### ----- 3.2. Fitting the model with INLA half gaussian --- ####
HN.prior <- "expression:
  tau0 = 0.01;
  sigma = exp(-theta/2);
  log_dens = log(2) - 0.5 * log(2 * pi) + 0.5 * log(tau0);
  log_dens = log_dens - 0.5 * tau0 * sigma^2;
  log_dens = log_dens - log(2) - theta / 2;
  return(log_dens);
"

half.normal  <- list(prec = list(prior = HN.prior))

model.inla2 <- inla(formula = formula,
                   data   = data.frame(y),
                   control.family = list(hyper = half.normal),
                   control.inla = list(strategy = "laplace", int.strategy = "grid"))



#### Comparing
prec1 <- model.inla1$marginals.hyperpar$`Precision for the Gaussian observations`
prec2 <- model.inla2$marginals.hyperpar$`Precision for the Gaussian observations`

sigma1 <- inla.tmarginal(function(x)sqrt(1/x), prec1)
sigma2 <- inla.tmarginal(function(x)sqrt(1/x), prec2)


plot(prec1, type = "l")
lines(prec2, type = "l", col = "red")

plot(sigma1, type = "l")
lines(sigma2, type = "l", col = "red")



####
ni <- 200000
nb <- 2000
nt <- 5

nc <- 3
# #
#     ni <- 1000
#     nb <- 100



data_jags <- list(y = y,
                  N = length(y))

## Initial values
inits <- function(){list(mu1 = rnorm(1, 0, 1),
                         sigma1 = runif(1, 0,1))}

## Parameters of interest
parameters <- c('mu1', 'sigma1')

cat("
    model {
    #model
    for (i in 1:N){
    y[i] ~ dnorm(mu, tau1)
    }


    mu ~ dnorm(0, 0.01)
    tau1 <- 1/(sigma1*sigma1)
   sigma1 ~ dnorm(0,0.01)T(0,)

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



model.jags$BUGSoutput$sims.list$sigma1

plot(sigma1, type = "l")
lines(sigma2, type = "l", col = "red")
model.jags$BUGSoutput$sims.list$sigma1 %>% as.matrix(.) %>%   hist(., breaks = 100, add = TRUE, freq = FALSE)








