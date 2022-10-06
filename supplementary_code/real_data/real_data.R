### In this script an example extracted from https://zenodo.org/record/2552025#.Yta8D3ZByUl   ###
# We compare dirinla package  with r-jags with enough amount of #
# simulations to guarantee convergence of the method, and with long r-jags,   #
# which has a large amount of simulations. The model that we try to fit is :  #
# --- Y \sim Dirich(exp(eta_1), exp(eta_2), ... , exp(eta_4))                 #
# --- --- eta_1 = beta_{01} + bio1 beta1[1] + bio12 beta2[1],                                   #
# --- --- eta_2 = beta_{02} + bio1 beta1[2] + bio12 beta2[2],                                   #
# --- --- eta_3 = beta_{03} + bio1 beta1[3] + bio12 beta2[3],                                   #
# --- --- eta_4 = beta_{04} + bio1 beta1[4] + bio12 beta2[4],                                   #
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

library(dplyr)

### --- 2. Reading the data --- ####
#data <- read.csv("supplementary_code/real_data_2/arabidopsis/ath_accessions.csv")
data <- read.csv("ath_accessions.csv")

data[,-c(1:7)] <- scale(data[,-c(1:7)])
Glc <- data
### --- 3. Transforming the data --- ####
Glc$Y <- as.matrix(data[, paste0("gc", 1:4)]) #package DirichletReg

### --- 4. Comparing posterior distributions. Jags vs INLA --- ####
### ----- 4.1. Fitting the model with jags --- ####
## MCMC configuration
ni <- 20000
nt <- 5
nb <- 2000
nc <- 3

## Data set
data_jags <- list(y = Glc$Y,
                  N = dim(Glc$Y)[1],
                  d = dim(Glc$Y)[2],
                  V = data[,c("bio1", "bio12")])

## Initial values
inits <- function(){list(beta0 = rnorm(4, 0, 1),
                         beta1 = rnorm(4, 0, 1),
                         beta2 = rnorm(4, 0, 1))}

## Parameters of interest
parameters <- c('beta0', 'beta1', ' beta2')

cat("
    model {
    #model
    for (i in 1:N){
    y[i, ] ~ ddirch(alpha[i,])

    log(alpha[i,1]) <- beta0[1] + beta1[1]*V[i,1] + beta2[1]*V[i,2]
    log(alpha[i,2]) <- beta0[2] + beta1[2]*V[i,1] + beta2[2]*V[i,2]
    log(alpha[i,3]) <- beta0[3] + beta1[3]*V[i,1] + beta2[3]*V[i,2]
    log(alpha[i,4]) <- beta0[4] + beta1[4]*V[i,1] + beta2[4]*V[i,2]
    }

    #priors
    for (c in 1:d)
    {
    beta0[c]  ~ dnorm(0, 0.01)
    beta1[c] ~ dnorm(0, 0.01)
    beta2[c] ~ dnorm(0, 0.01)
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

model.inla <- dirinlareg( formula  = y ~ 1 + bio1 + bio12 | 1 + bio1 + bio12  | 1 +bio1 + bio12 | 1 + bio1 + bio12 ,
                          y        = Glc$Y,
                          data.cov = Glc,
                          prec     = 0.01,
                          verbose  = TRUE)

t_inla <- proc.time()-t    # Stop the time
summary(model.inla)



y2 <- DR_data(Glc$Y)
mod_freq <- DirichletReg::DirichReg(y2 ~ 1 + bio1 + bio12 | 1 + bio1 + bio12  | 1 +bio1 + bio12 | 1 + bio1 + bio12 ,
data = Glc)
### ----- 4.3. Fitting the model with long jags --- ####
## MCMC configuration
ni <- 1000000
nt <- 10
nb <- 100000
nc <- 3



inits <- function(){list(beta0 = rnorm(4, 0, 1),
                         beta1 = rnorm(4, 0, 1),
                         beta2 = rnorm(4, 0, 1))}

## Parameters of interest
parameters <- c('beta0', 'beta1', ' beta2')

cat("
    model {
    #model
    for (i in 1:N){
    y[i, ] ~ ddirch(alpha[i,])

    log(alpha[i,1]) <- beta0[1] + beta1[1]*V[i,1] + beta2[1]*V[i,2]
    log(alpha[i,2]) <- beta0[2] + beta1[2]*V[i,1] + beta2[2]*V[i,2]
    log(alpha[i,3]) <- beta0[3] + beta1[3]*V[i,1] + beta2[3]*V[i,2]
    log(alpha[i,4]) <- beta0[4] + beta1[4]*V[i,1] + beta2[4]*V[i,2]
    }

    #priors
    for (c in 1:d)
    {
    beta0[c]  ~ dnorm(0, 0.01)
    beta1[c] ~ dnorm(0, 0.01)
    beta2[c] ~ dnorm(0, 0.01)
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
saveRDS(file = paste0("times_real.RDS"), times)

times <- readRDS("times_real.RDS")

### ----- 4.2. (E(INLA) - E(JAGS2))/SD(JAGS2) and variance ratios --- ####
ratio1_beta0 <- ratio2_beta0 <- ratio1_beta1 <- ratio2_beta1 <- ratio1_beta2 <- ratio2_beta2 <- numeric()
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


for (i in 1:4)
{
  mean_jags_2 <- mean(model.jags.2$BUGSoutput$sims.list$beta2[,i])
  sd_jags_2 <- sd(model.jags.2$BUGSoutput$sims.list$beta2[,i])
  # mean_jags_1 <- mean(model.jags$BUGSoutput$sims.list$beta0[,1])
  # sd_jags_1 <- sd(model.jags$BUGSoutput$sims.list$beta0[,1])
  mean_inla <- model.inla$summary_fixed[[i]]$mean[3]

  ratio1_beta2 <- c(ratio1_beta2, c(mean_inla - mean_jags_2)/sd_jags_2)
  ratio2_beta2 <- c(ratio2_beta2, sd(inla.rmarginal(10000, model.inla$marginals_fixed[[i]][[3]]))^2/(sd_jags_2^2))
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


### Beta2
result_beta_2 <- numeric()
for(i in 1:4)
{
  result_beta_2 <- rbind(result_beta_2,
                        t(matrix(c(model.jags$BUGSoutput$summary[paste0("beta1[", i,"]"), c("mean", "sd")],
                                   model.inla$summary_fixed[[i]][2,c("mean", "sd")],
                                   model.jags.2$BUGSoutput$summary[paste0("beta1[", i,"]"), c("mean", "sd")]))))
}
rownames(result_beta_2) <- paste0("beta1", 1:4)
colnames(result_beta_2) <- c(paste0("JAGS", c("_mean", "_sd")),
                            paste0("INLA", c("_mean", "_sd")),
                            paste0("LONG_JAGS", c("_mean", "_sd")))

total <- list(times = times,
     intercepts = result_beta0,
     slopes     = result_beta1,
     ratio1_intercepts = ratio1_beta0,
     ratio1_slopes = ratio1_beta1,
     ratio2_intercepts = ratio2_beta0,
     ratio2_slopes = ratio2_beta1,
     ratio1_slopes2 = ratio1_beta2,
     ratio2_slopes2 = ratio2_beta2
)

### --- 5. Plotting ---
model.inla <- readRDS("model_inla_real.RDS")
model.jags <- readRDS("model_jags_real.RDS")
model.jags.2 <- readRDS("model_jags_long_real.RDS")

### ----- 5.1. intercepts --- ####
## Intercept
p1 <- list()
beta0 <- expression(paste("p(", beta[0], "|", "y)"))

for (i in 1:length(model.inla$marginals_fixed))
{
  #jags1
  dens <- density(model.jags$BUGSoutput$sims.list$beta0[,i], adjust = 2)
  dens <- as.data.frame(cbind(dens$x, dens$y))
  colnames(dens) <- c("x", "y")

  #jags2
  dens2 <- density(model.jags.2$BUGSoutput$sims.list$beta0[,i], adjust = 2)

  dens2 <- as.data.frame(cbind(dens2$x, dens2$y))
  colnames(dens2) <- c("x", "y")

  #Data combining jags (1) and inla (2)
  dens <- rbind(cbind(dens, group = 1),
                cbind(as.data.frame(inla.smarginal(model.inla$marginals_fixed[[i]][[1]])), group = 2),
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
#  p1[[i]] <- p1[[i]] + geom_vline(xintercept = mod_freq$coefficients[seq(1,12, by = 3)][i])



  ### --- legend --- ###
  p1[[i]]<- p1[[i]] + theme(legend.position   = c(0.2, 0.8),
                            legend.title      = element_blank(),
                            legend.background = element_rect(colour = "gray"),
                            legend.key        = element_rect(colour = "white", fill="white"),
                            legend.key.size   = unit(0.5, "cm"),
                            legend.text       = element_text(size = 15),
                            axis.title        = element_text(size = 14))
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

  p1[[i]] <- p1[[i]] + ggtitle(colnames(Glc[,-c(1:3)])[i]) +
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
d <- 4

for (i in 1:length(model.inla$marginals_fixed))
{
  #jags1
  dens <- density(model.jags$BUGSoutput$sims.list$beta1[,i], adjust = 2)
  dens <- as.data.frame(cbind(dens$x, dens$y))
  colnames(dens) <- c("x", "y")

  #jags2
  dens2 <- density(model.jags.2$BUGSoutput$sims.list$beta1[,i], adjust = 2)
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
  #p2[[i]] <- p2[[i]] + geom_vline(xintercept = mod_freq$coefficients[seq(2,12, by = 3)][i])



  ### --- legend --- ###
  p2[[i]]<- p2[[i]] + theme(legend.position   = c(0.2, 0.8),
                            legend.title      = element_blank(),
                            legend.background = element_rect(colour = "gray"),
                            legend.key        = element_rect(colour = "white", fill="white"),
                            legend.key.size   = unit(0.5, "cm"),
                            legend.text       = element_text(size = 15),
                            axis.title        = element_text(size = 14)) +
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

p3 <- list()
beta1 <- expression(paste("p(", beta[2], "|", "y)"))
d <- 4

for (i in 1:length(model.inla$marginals_fixed))
{
  #jags1
  dens <- density(model.jags$BUGSoutput$sims.list$beta2[,i], adjust = 2)
  dens <- as.data.frame(cbind(dens$x, dens$y))
  colnames(dens) <- c("x", "y")

  #jags2
  dens2 <- density(model.jags.2$BUGSoutput$sims.list$beta2[,i], adjust = 2)
  dens2 <- as.data.frame(cbind(dens2$x, dens2$y))
  colnames(dens2) <- c("x", "y")

  #Data combining jags (1) and inla (2)
  dens <- rbind(cbind(dens, group = 1),
                cbind(as.data.frame(model.inla$marginals_fixed[[i]][[3]]), group = 2),
                cbind(dens2, group = 3))
  dens$group <- factor(dens$group,
                       labels = c("R-JAGS", "dirinla", "long R-JAGS"))

  ### Intercept
  p3[[i]] <- ggplot(dens,
                    aes(x = x,
                        y = y,
                        group = group
                        #colour = factor(group)
                    ))  +
    geom_line(size = 0.6, aes(linetype = group ,
                              color    = group)) +
    xlim(c(min(dens$x[dens$group=="R-JAGS"]), max(dens$x[dens$group=="R-JAGS"]))) +
    theme_bw() + #Show axes
    xlab(expression(beta[2])) + #xlab
    ylab(beta1) #ylab


  #Frequentist approach
 # p3[[i]] <- p3[[i]] + geom_vline(xintercept = mod_freq$coefficients[seq(3,12, by = 3)][i])



  ### --- legend --- ###
  p3[[i]]<- p3[[i]] + theme(legend.position   = c(0.2, 0.8),
                            legend.title      = element_blank(),
                            legend.background = element_rect(colour = "gray"),
                            legend.key        = element_rect(colour = "white", fill="white"),
                            legend.key.size   = unit(0.5, "cm"),
                            legend.text       = element_text(size = 15),
                            axis.title        = element_text(size = 14)) +
    # scale_fill_manual(labels=c("R-JAGS", "dirinla", "long R-JAGS"),
    #                   values = c("darkgreen", "red4", "blue4" )) +
    scale_colour_manual (
      values= c("darkgreen", "red4", "blue4")) +
    scale_linetype_manual(labels=c("R-JAGS", "dirinla", "long R-JAGS"),
                          values=c("dotted", "twodash",  "solid"))

  if(i!=5)
  {
    p3[[i]] <- p3[[i]] + theme(legend.position="none")
  }

  p3[[i]] <- p3[[i]] + ggtitle("") +
    theme(
      plot.title = element_text(color = "black",
                                size  = 15,
                                face  = "bold.italic",
                                hjust = 0.5))
}
# pdf("example_simulation2_slopes_50.pdf", width = 18, height = 4)
   gridExtra::grid.arrange(p2[[1]], p2[[2]], p2[[3]], p2[[4]], ncol = 4)
# dev.off()






# pdf("example_real_slopes_50.pdf", width = 18, height = 4)
# grid.arrange(p2[[1]], p2[[2]], p2[[3]], p2[[4]], ncol = 4)
# dev.off()



# pdf("examples_real_slopes_intercepts2.pdf", width = 15, height = 10)
# gridExtra::grid.arrange(p1[[1]], p1[[2]], p1[[3]], p1[[4]],
#              p2[[1]], p2[[2]], p2[[3]], p2[[4]],
#              p3[[1]], p3[[2]], p3[[3]], p3[[4]], ncol = 4)
# dev.off()

pl_combined <-
  ((p1[[1]] | p1[[2]] | p1[[3]] | p1[[4]]) /
   (p2[[1]] | p2[[2]] | p2[[3]] | p2[[4]]) /
   (p3[[1]] | p3[[2]] | p3[[3]] | p3[[4]])) +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "bottom",
                 legend.key.width = unit(1.5,"cm")) &
  guides(col = guide_legend(nrow=1,byrow=TRUE),
         linetype = guide_legend(override.aes = list(size = 1.1)))

pdf("examples_real_slopes_intercepts.pdf", width = 15, height = 10)
  pl_combined
dev.off()


saveRDS(total, file = "real_data.RDS")



### --- 5. Computing also ratios for shortJAGS --- ####
ratios_jags <- function()
{
  model.jags <- readRDS(paste0("model_jags_", "real",".RDS"))
  model.inla <- readRDS(paste0("model_inla_", "real",".RDS"))
  model.jags.2 <- readRDS(paste0("model_jags_long_", "real",".RDS"))

  ratio1_beta0_jags <- ratio2_beta0_jags <- ratio1_beta1_jags <- ratio2_beta1_jags <- ratio1_beta2_jags <- ratio2_beta2_jags <- numeric()
  for (i in 1:4)
  {
    mean_jags_2 <- mean(model.jags.2$BUGSoutput$sims.list$beta0[,i])
    sd_jags_2 <- sd(model.jags.2$BUGSoutput$sims.list$beta0[,i])
    mean_jags_1 <- mean(model.jags$BUGSoutput$sims.list$beta0[,i])
    sd_jags_1 <- sd(model.jags$BUGSoutput$sims.list$beta0[,i])

    ratio1_beta0_jags <- c(ratio1_beta0_jags, (mean_jags_1 - mean_jags_2)/sd_jags_2)
    ratio2_beta0_jags <- c(ratio2_beta0_jags, (sd_jags_1^2)/(sd_jags_2^2))
  }

  for (i in 1:4)
  {
    mean_jags_2 <- mean(model.jags.2$BUGSoutput$sims.list$beta1[,i])
    sd_jags_2 <- sd(model.jags.2$BUGSoutput$sims.list$beta1[,i])
    mean_jags_1 <- mean(model.jags$BUGSoutput$sims.list$beta1[,i])
    sd_jags_1 <- sd(model.jags$BUGSoutput$sims.list$beta1[,i])

    ratio1_beta1_jags <- c(ratio1_beta1_jags, (mean_jags_1 - mean_jags_2)/sd_jags_2)
    ratio2_beta1_jags <- c(ratio2_beta1_jags, (sd_jags_1^2)/(sd_jags_2^2))
  }

  for (i in 1:4)
  {
    mean_jags_2 <- mean(model.jags.2$BUGSoutput$sims.list$beta2[,i])
    sd_jags_2 <- sd(model.jags.2$BUGSoutput$sims.list$beta2[,i])
    mean_jags_1 <- mean(model.jags$BUGSoutput$sims.list$beta2[,i])
    sd_jags_1 <- sd(model.jags$BUGSoutput$sims.list$beta2[,i])

    ratio1_beta2_jags <- c(ratio1_beta2_jags, (mean_jags_1 - mean_jags_2)/sd_jags_2)
    ratio2_beta2_jags <- c(ratio2_beta2_jags, (sd_jags_1^2)/(sd_jags_2^2))
  }

  list(ratio1_beta0_jags = ratio1_beta0_jags,
       ratio2_beta0_jags = ratio2_beta0_jags,
       ratio1_beta1_jags = ratio1_beta1_jags,
       ratio2_beta1_jags = ratio2_beta1_jags,
       ratio1_beta2_jags = ratio1_beta2_jags,
       ratio2_beta2_jags = ratio2_beta2_jags)
}

res_ratios <- ratios_jags()

saveRDS(res_ratios, file = "simulation_real_ratios_jags.RDS")


### --- 6. Extracting tables for the paper --- ####
results <- readRDS(file = "real_data.RDS")
res_ratios <- readRDS(file = "simulation_real_ratios_jags.RDS")

#dirinla ratios
result_ratio1 <- t(matrix(c(results$ratio1_intercepts, results$ratio1_slopes, results$ratio1_slopes2)))
colnames(result_ratio1) <- c(paste0("beta0", 1:4), paste0("beta1", 1:4), paste0("beta2", 1:4))

result_ratio2 <- t(matrix(sqrt(c(results$ratio2_intercepts, results$ratio2_slopes, results$ratio2_slopes2))))
colnames(result_ratio2) <- c(paste0("beta0", 1:4), paste0("beta1", 1:4), paste0("beta2", 1:4))

#jags ratios
result_ratio1_jags <- matrix(c(res_ratios$ratio1_beta0_jags, res_ratios$ratio1_beta1_jags, res_ratios$ratio1_beta2_jags), ncol = 12)
colnames(result_ratio1_jags) <- c(paste0("beta0", 1:4), paste0("beta1", 1:4), paste0("beta2", 1:4))

result_ratio2_jags <- matrix(sqrt(c(res_ratios$ratio2_beta0_jags, res_ratios$ratio2_beta1_jags, res_ratios$ratio2_beta2_jags )), ncol = 12)
colnames(result_ratio2_jags) <- c(paste0("beta0", 1:4), paste0("beta1", 1:4), paste0("beta2", 1:4))


#Latex
library(xtable)
xtable(result_ratio1, digits = 4)
xtable(result_ratio1_jags, digits = 4)

xtable(result_ratio2, digits = 4)
xtable(result_ratio2_jags, digits = 4)

total <- rbind(result_ratio1, result_ratio1_jags, result_ratio2, result_ratio2_jags)
xtable(t(total), digits = 4)



