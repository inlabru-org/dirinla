### In this script simulations with N=100 are conducted in order to check    ###
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

### --- 2. Simulation data --- ####
set.seed(1000)
n <- 100

#Covariates
V <- as.data.frame(matrix(runif((10)*n, 0, 1), ncol=10))
names(V) <- paste0('v', 1:(10))

# Formula that we want to fit
formula <- y ~ 1 + v1 | 1 + v2 | 1 + v3 | 1 + v4

names_cat <- formula_list(formula)

# Parameters to fit
x <- c(-1.5, 1, -3, 1.5,
       2, -3 , -1, 5)
#x <- c(-2.4, 4, -3.1, 4)

mus <- exp(x)/sum(exp(x))
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

### --- 3. Comparing posterior distributions. jags vs INLA --- ####
### ----- 3.1. Fitting the model with jags --- ####
## MCMC configuration
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
                         beta1 = rnorm(d, 0, 1))}

## Parameters of interest
parameters <- c('beta0', 'beta1')

cat("
    model {
    #model
    for (i in 1:N){
    y[i, ] ~ ddirch(alpha[i,])

    log(alpha[i,1]) <- beta0[1] + beta1[1]*V[i,1]
    log(alpha[i,2]) <- beta0[2] + beta1[2]*V[i,2]
    log(alpha[i,3]) <- beta0[3] + beta1[3]*V[i,3]
    log(alpha[i,4]) <- beta0[4] + beta1[4]*V[i,4]
    }

    #priors
    for (c in 1:d)
    {
    beta0[c]  ~ dnorm(0, 0.0001)
    beta1[c] ~ dnorm(0, 0.0001)
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
t_jags<-proc.time()-t    # Stop the time
print(model.jags)





### ----- 3.2. Fitting the model with INLA --- ####
t <- proc.time() # Measure the time
model.inla <- dirinlareg( formula  = y ~ 1 + v1 | 1 + v2 | 1 + v3 | 1 + v4  ,
                          y        = y,
                          data.cov = V,
                          prec     = 0.0001,
                          verbose  = TRUE)

t_inla <- proc.time()-t    # Stop the time
summary(model.inla)

### ----- 3.3. Fitting the model with long jags --- ####
## MCMC configuration
ni <- 1000000
nt <- 5
nb <- 100000
nc <- 3

## Data set
data_jags <- list(y = y,
                  N = dim(y)[1],
                  d = d,
                  V = V)

## Initial values
inits <- function(){list(beta0 = rnorm(d, 0, 1))}

## Parameters of interest
parameters <- c('beta0', 'beta1')

cat("
    model {
    #model
    for (i in 1:N){
    y[i, ] ~ ddirch(alpha[i,])

    log(alpha[i,1]) <- beta0[1] + beta1[1]*V[i,1]
    log(alpha[i,2]) <- beta0[2] + beta1[2]*V[i,2]
    log(alpha[i,3]) <- beta0[3] + beta1[3]*V[i,3]
    log(alpha[i,4]) <- beta0[4] + beta1[4]*V[i,4]
    }

    #priors
    for (c in 1:d)
    {
    beta0[c]  ~ dnorm(0, 0.0001)
    beta1[c] ~ dnorm(0, 0.0001)    }
    }", file="model_cov.jags" )



## Call jags
t <- proc.time() # Measure the time
model.jags.2 <- jags(data_jags,
                     inits,
                     parameters,
                     "model_cov.jags",
                     n.chains          = nc,
                     n.thin            = nt,
                     n.iter            = ni,
                     n.burnin          = nb,
                     working.directory = getwd()) #
t_jags_2<-proc.time()-t    # Stop the time

### ----- 3.4. Comparing methodologies --- ####
### ------- 3.4.1. intercepts --- ####
### ------- 3.4.1.1. Using plots --- ####
pdf("example_simulation2_intercepts_100.pdf", width = 20, height = 5)
beta0_not <- expression(beta[0])
p_beta0_not <- expression(paste("p(", beta[0], "|", "y)"))
par(mfrow=c(1,4))
for (i in 1:4)
{
  plot(density(model.jags$BUGSoutput$sims.list$beta0[,i]),
       col = "orange",
       lwd  = 2,
       type = "l",
       ylim= c(0, max(model.inla$marginals_fixed[[i]][[1]][,2],
                      density(model.jags$BUGSoutput$sims.list$beta0[,i])$y)),
       xlab = beta0_not,
       ylab = p_beta0_not,
       main = beta0_not)

  lines(model.inla$marginals_fixed[[i]][[1]],
        col="red",
        lwd=2)
  lines(density(model.jags.2$BUGSoutput$sims.list$beta0[,i]),
        col = "blue")
  abline(v = x[i],
         col = "black",
         lwd = 2)

  legend("topright", legend=c("R-jags", "dirinla", "R-long-jags"),
         col = c("orange", "red", "blue"),
         lty = 1,
         lwd = 2)
}
dev.off()


### ------- 3.4.1.2. Using ggplot --- ####
## Intercept
p1 <- list()
beta0 <- expression(paste("p(", beta[0], "|", "y)"))

for (i in 1:length(model.inla$marginals_fixed))
{
  #jags1
  dens <- density(model.jags$BUGSoutput$sims.matrix[,i])
  dens <- as.data.frame(cbind(dens$x, dens$y))
  colnames(dens) <- c("x", "y")

  #jags2
  dens2 <- density(model.jags.2$BUGSoutput$sims.matrix[,i])
  dens2 <- as.data.frame(cbind(dens2$x, dens2$y))
  colnames(dens2) <- c("x", "y")

  #Data combining jags (1) and inla (2)
  dens <- rbind(cbind(dens, group = 1),
                cbind(as.data.frame(model.inla$marginals_fixed[[i]][[1]]), group = 2),
                cbind(dens2, group = 3))
  dens$group <- factor(dens$group,
                       labels = c("R-JAGS", "dirinla", "long R-JAGS"))

  ### Intercept
  p1[[i]] <- ggplot(dens,
                    aes(x = x,
                        y = y,
                        group = group
                        #colour = factor(group)
                    ))  +
    geom_line(size = 0.6, aes(linetype = group ,
                              color    = group)) +
    xlim(c(min(dens$x[dens$group=="R-JAGS"]), max(dens$x[dens$group=="R-JAGS"]))) +
    theme_bw() + #Show axes
    xlab(expression(beta[0])) + #xlab
    ylab(beta0) #ylab


  #Frequentist approach
  p1[[i]] <- p1[[i]] + geom_vline(xintercept = x[seq(1,4, by=1)][i])



  ### --- legend --- ###
  p1[[i]]<- p1[[i]] + theme(legend.position   = c(0.2, 0.8),
                            legend.title      = element_blank(),
                            legend.background = element_rect(colour = "gray"),
                            legend.key        = element_rect(colour = "white", fill="white"),
                            legend.key.size   = unit(0.5, "cm")) +
    theme(legend.text = element_text(size = 9)) +
    # scale_fill_manual(labels=c("R-JAGS", "dirinla", "long R-JAGS"),
    #                   values = c("darkgreen", "red4", "blue4" )) +
    scale_colour_manual (
      values= c("darkgreen", "red4", "blue4")) +
    scale_linetype_manual(labels=c("R-JAGS", "dirinla", "long R-JAGS"),
                          values=c("dotted", "twodash", "solid"))

  if(i!=1)
  {
    p1[[i]] <- p1[[i]] + theme(legend.position="none")
  }

  p1[[i]] <- p1[[i]] + ggtitle(paste0("Category ", i)) +
    theme(
      plot.title = element_text(color = "black",
                                size  = 15,
                                face  = "bold.italic",
                                hjust = 0.5))
}


pdf("example_simulation2_intercepts_100.pdf", width = 18, height = 4)
grid.arrange(p1[[1]], p1[[2]], p1[[3]], p1[[4]], ncol = 4)
dev.off()



### ------- 3.4.2. slopes --- ####
### ------- 3.4.2.1. Using plots --- ####
pdf("example_simulation2_slope_100.pdf", width = 10, height = 3)
beta1_not <- expression(beta[1])
p_beta1_not <- expression(paste("p(", beta[1], "|", "y)"))
par(mfrow=c(1,4))
for (i in 1:4)
{
  plot(density(model.jags$BUGSoutput$sims.list$beta1[,i]),
       col = "orange",
       lwd  = 2,
       type = "l",
       ylim= c(0, max(model.inla$marginals_fixed[[i]][[2]][,2],
                      density(model.jags$BUGSoutput$sims.list$beta1[,i])$y)),
       # xlim = c(min(model.inla$marginals_fixed[[i]][[2]][,1],
       #              density(model.jags$BUGSoutput$sims.list$beta1[,i])$x),
       #          max(model.inla$marginals_fixed[[i]][[2]][,1],
       #              density(model.jags$BUGSoutput$sims.list$beta1[,i])$x)),
       xlab = beta1_not,
       ylab = p_beta1_not,
       main = beta1_not)

  lines(model.inla$marginals_fixed[[i]][[2]],
        col="red",
        lwd=2)
  lines(density(model.jags.2$BUGSoutput$sims.list$beta1[,i]),
        col = "blue")
  abline(v = x[i + 4],
         col = "black",
         lwd = 2)

  legend("topright", legend=c("R-jags", "dirinla", "R-long-jags"),
         col = c("orange", "red", "blue"),
         lty = 1,
         lwd = 2)
}
dev.off()


### ------- 3.4.2.2. Using ggplot --- ####
## Intercept
p2 <- list()
beta1 <- expression(paste("p(", beta[1], "|", "y)"))

for (i in 1:length(model.inla$marginals_fixed))
{
  #jags1
  dens <- density(model.jags$BUGSoutput$sims.matrix[,i + d])
  dens <- as.data.frame(cbind(dens$x, dens$y))
  colnames(dens) <- c("x", "y")

  #jags2
  dens2 <- density(model.jags.2$BUGSoutput$sims.matrix[,i + d])
  dens2 <- as.data.frame(cbind(dens2$x, dens2$y))
  colnames(dens2) <- c("x", "y")

  #Data combining jags (1) and inla (2)
  dens <- rbind(cbind(dens, group = 1),
                cbind(as.data.frame(model.inla$marginals_fixed[[i]][[2]]), group = 2),
                cbind(dens2, group = 3))
  dens$group <- factor(dens$group,
                       labels = c("R-JAGS", "dirinla", "long R-JAGS"))

  ### Intercept
  p2[[i]] <- ggplot(dens,
                    aes(x = x,
                        y = y,
                        group = group
                        #colour = factor(group)
                    ))  +
    geom_line(size = 0.6, aes(linetype = group ,
                              color    = group)) +
    xlim(c(min(dens$x[dens$group=="R-JAGS"]), max(dens$x[dens$group=="R-JAGS"]))) +
    theme_bw() + #Show axes
    xlab(expression(beta[1])) + #xlab
    ylab(beta1) #ylab


  #Frequentist approach
  p2[[i]] <- p2[[i]] + geom_vline(xintercept = x[seq(5,8, by=1)][i])



  ### --- legend --- ###
  p2[[i]]<- p2[[i]] + theme(legend.position   = c(0.2, 0.8),
                            legend.title      = element_blank(),
                            legend.background = element_rect(colour = "gray"),
                            legend.key        = element_rect(colour = "white", fill="white"),
                            legend.key.size   = unit(0.5, "cm")) +
    theme(legend.text = element_text(size = 9)) +
    # scale_fill_manual(labels=c("R-JAGS", "dirinla", "long R-JAGS"),
    #                   values = c("darkgreen", "red4", "blue4" )) +
    scale_colour_manual (
      values= c("darkgreen", "red4", "blue4")) +
    scale_linetype_manual(labels=c("R-JAGS", "dirinla", "long R-JAGS"),
                          values=c("dotted", "twodash", "solid"))

  if(i!=5)
  {
    p2[[i]] <- p2[[i]] + theme(legend.position="none")
  }

  p2[[i]] <- p2[[i]] + ggtitle("") +
    theme(
      plot.title = element_text(color = "black",
                                size  = 15,
                                face  = "bold.italic",
                                hjust = 0.5))
}


pdf("example_simulation2_slopes_100.pdf", width = 18, height = 4)
grid.arrange(p2[[1]], p2[[2]], p2[[3]], p2[[4]], ncol = 4)
dev.off()



pdf("examples_simualtion2_slopes_intercepts_100.pdf", width = 15, height = 6)
grid.arrange(p1[[1]], p1[[2]], p1[[3]], p1[[4]],
             p2[[1]], p2[[2]], p2[[3]], p2[[4]], ncol = 4)
dev.off()
