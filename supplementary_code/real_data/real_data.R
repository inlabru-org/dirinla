### In this script an example extracted from Aitchison 1986 is analyzed:   ###
# Galcialtills. We compare dirinla package  with r-jags with enough amount of #
# simulations to guarantee convergence of the method, and with long r-jags,   #
# which has a large amount of simulations. The model that we try to fit is :  #
# --- Y \sim Dirich(exp(eta_1), exp(eta_2), ... , exp(eta_4))                 #
# --- --- eta_1 = beta_{01} + Pcount beta1,                                   #
# --- --- eta_2 = beta_{02} + Pcount beta2,                                   #
# --- --- eta_3 = beta_{03} + Pcount beta3,                                   #
# --- --- eta_4 = beta_{04} + Pcount beta4,                                   #
# ----------------------------------------------------------------------------#
setwd("~/GIT1/dirinla/supplementary_code/real_data")
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

### --- 2. Reading the data --- ####
Glc <- GlacialTills
### --- 3. Transforming the data --- ####
Glc$Y <- DR_data(Glc[,1:4]/100, trafo = TRUE)[,1:4] #package DirichletReg
Glc$Pcount <- Glc$Pcount/100


### --- 4. Comparing posterior distributions. Jags vs INLA --- ####
### ----- 4.1. Fitting the model with jags --- ####
## MCMC configuration
ni <- 2000
nt <- 5
nb <- 200
nc <- 3

## Data set
data_jags <- list(y = Glc$Y,
                  N = dim(Glc$Y)[1],
                  d = dim(Glc$Y)[2],
                  V = Glc$Pcount)

## Initial values
inits <- function(){list(beta0 = rnorm(4, 0, 1),
                         beta1 = rnorm(4, 0, 1))}

## Parameters of interest
parameters <- c('beta0', 'beta1')

cat("
    model {
    #model
    for (i in 1:N){
    y[i, ] ~ ddirch(alpha[i,])

    log(alpha[i,1]) <- beta0[1] + beta1[1]*V[i]
    log(alpha[i,2]) <- beta0[2] + beta1[2]*V[i]
    log(alpha[i,3]) <- beta0[3] + beta1[3]*V[i]
    log(alpha[i,4]) <- beta0[4] + beta1[4]*V[i]
    }

    #priors
    for (c in 1:d)
    {
    beta0[c]  ~ dnorm(0, 0.01)
    beta1[c] ~ dnorm(0, 0.01)
    }
    }", file="model_real.jags" )



## Call jags
t <- proc.time() # Measure the time
model.jags <- jags(data_jags,
                   inits,
                   parameters,
                   "model_real.jags",
                   n.chains          = nc,
                   n.thin            = nt,
                   n.iter            = ni,
                   n.burnin          = nb,
                   working.directory = getwd()) #
t_jags<-proc.time()-t    # Stop the time
print(model.jags)




### ----- 4.2. Fitting the model with INLA --- ####
t <- proc.time() # Measure the time
model.inla <- dirinlareg( formula  = y ~ 1 + Pcount | 1 + Pcount | 1 + Pcount | 1 + Pcount  ,
                          y        = Glc$Y,
                          data.cov = Glc,
                          prec     = 0.0001,
                          verbose  = TRUE)

t_inla <- proc.time()-t    # Stop the time
summary(model.inla)

### ----- 4.3. Fitting the model with long jags --- ####
## MCMC configuration
ni <- 100000
nt <- 5
nb <- 10000
nc <- 3
# ni <- 100
# nb <- 10

## Data set
data_jags <- list(y = Glc$Y,
                  N = dim(Glc$Y)[1],
                  d = dim(Glc$Y)[2],
                  V = Glc$Pcount)

## Initial values
inits <- function(){list(beta0 = rnorm(4, 0, 1),
                         beta1 = rnorm(4, 0, 1))}

## Parameters of interest
parameters <- c('beta0', 'beta1')

cat("
    model {
    #model
    for (i in 1:N){
    y[i, ] ~ ddirch(alpha[i,])

    log(alpha[i,1]) <- beta0[1] + beta1[1]*V[i]
    log(alpha[i,2]) <- beta0[2] + beta1[2]*V[i]
    log(alpha[i,3]) <- beta0[3] + beta1[3]*V[i]
    log(alpha[i,4]) <- beta0[4] + beta1[4]*V[i]
    }

    #priors
    for (c in 1:d)
    {
    beta0[c]  ~ dnorm(0, 0.01)
    beta1[c] ~ dnorm(0, 0.01)
    }
    }", file="model_real.jags" )



## Call jags
t <- proc.time() # Measure the time
model.jags.2 <- jags(data_jags,
                     inits,
                     parameters,
                     "model_real.jags",
                     n.chains          = nc,
                     n.thin            = nt,
                     n.iter            = ni,
                     n.burnin          = nb,
                     working.directory = getwd()) #
t_jags_2<-proc.time()-t    # Stop the time

saveRDS(file = paste0("model_jags_real", ".RDS"), model.jags)
saveRDS(file = paste0("model_jags_long_real", ".RDS"), model.jags.2)
saveRDS(file = paste0("model_inla_real.RDS"), model.inla)

### ----- 4.1. Computational times --- ####
times <- c(t_jags[3], t_inla[3], t_jags_2[3])

### ----- 4.2. (E(INLA) - E(JAGS2))/SD(JAGS2) and variance ratios --- ####
ratio1_beta0 <- ratio2_beta0 <- ratio1_beta1 <- ratio2_beta1 <- numeric()
for (i in 1:4)
{
  mean_jags_2 <- mean(model.jags.2$BUGSoutput$sims.list$beta0[,i])
  sd_jags_2 <- sd(model.jags.2$BUGSoutput$sims.list$beta0[,i])
  # mean_jags_1 <- mean(model.jags$BUGSoutput$sims.list$beta0[,1])
  # sd_jags_1 <- sd(model.jags$BUGSoutput$sims.list$beta0[,1])
  mean_inla <- model.inla$summary_fixed[[i]]$mean[1]

  ratio1_beta0 <- c(ratio1_beta0, c(mean_inla - mean_jags_2)/sd_jags_2)
  ratio2_beta0 <- c(ratio2_beta0, sd(inla.rmarginal(10000, model.inla$marginals_fixed[[i]]$intercept))^2/sd_jags_2^2)
}

for (i in 1:4)
{
  mean_jags_2 <- mean(model.jags.2$BUGSoutput$sims.list$beta1[,i])
  sd_jags_2 <- sd(model.jags.2$BUGSoutput$sims.list$beta1[,i])
  # mean_jags_1 <- mean(model.jags$BUGSoutput$sims.list$beta0[,1])
  # sd_jags_1 <- sd(model.jags$BUGSoutput$sims.list$beta0[,1])
  mean_inla <- model.inla$summary_fixed[[i]]$mean[2]

  ratio1_beta1 <- c(ratio1_beta1, c(mean_inla - mean_jags_2)/sd_jags_2)
  ratio2_beta1 <- c(ratio2_beta1, sd(inla.rmarginal(10000, model.inla$marginals_fixed[[i]][[2]]))^2/(sd_jags_2^2))
}

### ----- 4.3. Mean and sd of the posterior distributions --- ####
### Intercepts
result_beta0 <- numeric()
for(i in 1:4)
{
  result_beta0 <- rbind(result_beta0,
                        t(matrix(c(model.jags$BUGSoutput$summary[paste0("beta0[", i,"]"), c("mean", "sd")],
                                   model.inla$summary_fixed[[i]][1,c("mean", "sd")],
                                   model.jags.2$BUGSoutput$summary[paste0("beta0[", i,"]"), c("mean", "sd")]))))
}
rownames(result_beta0) <- paste0("beta0", 1:4)
colnames(result_beta0) <- c(paste0("JAGS", c("_mean", "_sd")),
                            paste0("INLA", c("_mean", "_sd")),
                            paste0("LONG_JAGS", c("_mean", "_sd")))


### Beta1
result_beta1 <- numeric()
for(i in 1:4)
{
  result_beta1 <- rbind(result_beta1,
                        t(matrix(c(model.jags$BUGSoutput$summary[paste0("beta1[", i,"]"), c("mean", "sd")],
                                   model.inla$summary_fixed[[i]][2,c("mean", "sd")],
                                   model.jags.2$BUGSoutput$summary[paste0("beta1[", i,"]"), c("mean", "sd")]))))
}
rownames(result_beta1) <- paste0("beta1", 1:4)
colnames(result_beta1) <- c(paste0("JAGS", c("_mean", "_sd")),
                            paste0("INLA", c("_mean", "_sd")),
                            paste0("LONG_JAGS", c("_mean", "_sd")))

total <- list(times = times,
     intercepts = result_beta0,
     slopes     = result_beta1,
     ratio1_intercepts = ratio1_beta0,
     ratio1_slopes = ratio1_beta1,
     ratio2_intercepts = ratio2_beta0,
     ratio2_slopes = ratio2_beta1)

### --- 5. Plotting ---
### ----- 5.1. intercepts --- ####
## Intercept
p1 <- list()
beta0 <- expression(paste("p(", beta[0], "|", "y)"))

for (i in 1:length(model.inla$marginals_fixed))
{
  #jags1
  dens <- density(model.jags$BUGSoutput$sims.matrix[,i], adjust = 2)
  dens <- as.data.frame(cbind(dens$x, dens$y))
  colnames(dens) <- c("x", "y")

  #jags2
  dens2 <- density(model.jags.2$BUGSoutput$sims.matrix[,i], adjust = 2)
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
 # p1[[i]] <- p1[[i]] + geom_vline(xintercept = x[seq(1,8, by=2)][i])



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
                          values=c("dotted", "twodash",  "solid"))

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

#
# pdf("example_simulation2_intercepts_50.pdf", width = 18, height = 4)
#   gridExtra::grid.arrange(p1[[1]], p1[[2]], p1[[3]], p1[[4]], ncol = 4)
# dev.off()



### ----- 5.2. slopes --- ####
p2 <- list()
beta1 <- expression(paste("p(", beta[1], "|", "y)"))

for (i in 1:length(model.inla$marginals_fixed))
{
  #jags1
  dens <- density(model.jags$BUGSoutput$sims.matrix[,i + d], adjust = 2)
  dens <- as.data.frame(cbind(dens$x, dens$y))
  colnames(dens) <- c("x", "y")

  #jags2
  dens2 <- density(model.jags.2$BUGSoutput$sims.matrix[,i + d], adjust = 2)
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
  # p2[[i]] <- p2[[i]] + geom_vline(xintercept = x[seq(2,8, by=2)][i])



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
                          values=c("dotted", "twodash",  "solid"))

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


# pdf("example_simulation2_slopes_50.pdf", width = 18, height = 4)
#   gridExtra::grid.arrange(p2[[1]], p2[[2]], p2[[3]], p2[[4]], ncol = 4)
# dev.off()



pdf(paste0("examples_simualtion2_slopes_intercepts_", n ,".pdf"), width = 15, height = 6)
gridExtra::grid.arrange(p1[[1]], p1[[2]], p1[[3]], p1[[4]],
                        p2[[1]], p2[[2]], p2[[3]], p2[[4]], ncol = 4)
dev.off()
list(times = times,
     intercepts = result_beta0,
     slopes     = result_beta1,
     ratio1_intercepts = ratio1_beta0,
     ratio1_slopes = ratio1_beta1,
     ratio2_intercepts = ratio2_beta0,
     ratio2_slopes = ratio2_beta1)


# pdf("example_real_slopes_50.pdf", width = 18, height = 4)
# grid.arrange(p2[[1]], p2[[2]], p2[[3]], p2[[4]], ncol = 4)
# dev.off()



pdf("examples_real_slopes_intercepts.pdf", width = 15, height = 6)
gridExtra::grid.arrange(p1[[1]], p1[[2]], p1[[3]], p1[[4]],
             p2[[1]], p2[[2]], p2[[3]], p2[[4]], ncol = 4)
dev.off()



saveRDS(total, file = "real_data.RDS")

### --- 4. Extracting tables for the paper --- ####
results <- readRDS(file = "real_data.RDS")

result_ratio1 <- t(matrix(c(results$ratio1_intercepts, results$ratio1_slopes)))
colnames(result_ratio1) <- c(paste0("beta0", 1:4), paste0("beta1", 1:4))

result_ratio2 <- t(matrix(sqrt(c(results$ratio2_intercepts, results$ratio2_slopes))))
colnames(result_ratio2) <- c(paste0("beta0", 1:4), paste0("beta1", 1:4))


#Latex
xtable(result_ratio1, digits = 4)
xtable(result_ratio2, digits = 4)
