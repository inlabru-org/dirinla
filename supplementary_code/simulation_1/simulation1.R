### In this script simulations with N=50, 100, 5000, 1000, 10000 are conducted #
# in order to check    ###
# how the package dirinla works. We compare with r-jags with enough amount of #
# simulations to guarantee convergence of the method, and with long r-jags,   #
# which has a large amount of simulations. The model that we try to fit is :  #
# --- Y \sim Dirich(exp(eta_1), exp(eta_2), ... , exp(eta_4))                 #
# --- --- eta_1 = beta_{01},                                                  #
# --- --- eta_2 = beta_{02},                                                  #
# --- --- eta_3 = beta_{03},                                                  #
# --- --- eta_4 = beta_{04},                                                  #
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

library(xtable)
library(dplyr)
library(cowplot)


#setwd("supplementary_code/simulation_1")

### --- 2. Function for simulation --- ####
simulations_just_intercepts <- function(n)
{
  ### --- 2. Simulation data --- ####
  set.seed(1000)

  cat("n = ", n, " -----> Simulating data \n")
  #Covariates
  V <- as.data.frame(matrix(runif((10)*n, 0, 1), ncol=10))
  names(V) <- paste0('v', 1:(10))

  # Formula that we want to fit
  formula <- y ~ 1 | 1 | 1 | 1

  names_cat <- formula_list(formula, y= NULL)

  # Parameters to fit
  x <- c(-2.4, 1.2, -3.1, 1.3)
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
  A_construct <- data_stack_construct
  eta <- A_construct %*% x
  alpha <- exp(eta)
  alpha <- matrix(alpha,
                  ncol  = d,
                  byrow = TRUE)
  y_o <- rdirichlet(n, alpha)
  colnames(y_o) <- paste0("y", 1:d)


  y <- y_o

  ### --- 3. Fitting models: jags vs INLA --- ####
  ### ----- 3.1. Fitting the model with jags --- ####
  cat(paste0("n = ", n, " -----> Fitting using SHORT JAGS \n"))
  if(file.exists(paste0("model_jags_", n,".RDS"))){
    model.jags <- readRDS(paste0("model_jags_", n,".RDS"))
    if(n < 1000)
    {
      simulation <- readRDS("simulation1_50-500.RDS")
      t_jags <- simulation[[paste0("n",n)]]$times[1]
    }else{
      simulation <- readRDS("simulation1_1000-10000.RDS")
      t_jags <- simulation[[paste0("n",n)]]$times[1]
    }

  }else{
    ## MCMC configuration
    ni <- 2000
    nt <- 5
    nb <- 200
    nc <- 3

    ## Data set
    data_jags <- list(y = y,
                      N = dim(y)[1],
                      d = d)

    ## Initial values
    inits <- function(){list(beta0 = rnorm(d, 0, 1))}

    ## Parameters of interest
    parameters <- c('beta0', 'alpha[1,1]', 'alpha[1,2]', 'alpha[1,3]', 'alpha[1,4]', 'mu')

    cat("
    model {
    #model
    for (i in 1:N){
    y[i, ] ~ ddirch(alpha[i,])

    log(alpha[i,1]) <- beta0[1]
    log(alpha[i,2]) <- beta0[2]
    log(alpha[i,3]) <- beta0[3]
    log(alpha[i,4]) <- beta0[4]

    }

    alpha0 <- alpha[1,1] + alpha[1,2] + alpha[1,3] + alpha[1,4]
    mu[1] <- alpha[1,1]/alpha0
    mu[2] <- alpha[1,2]/alpha0
    mu[3] <- alpha[1,3]/alpha0
    mu[4] <- alpha[1,4]/alpha0

    #priors
    for (c in 1:d)
    {
    beta0[c]  ~ dnorm(0, 0.0001)
    }
    }", file="model.jags" )



    ## Call jags
    t <- proc.time() # Measure the time
    model.jags <- jags(data_jags,
                       inits,
                       parameters,
                       "model.jags",
                       n.chains          = nc,
                       n.thin            = nt,
                       n.iter            = ni,
                       n.burnin          = nb,
                       working.directory = getwd()) #
    t_jags<-proc.time()-t    # Stop the time
    t_jags <- t_jags[3]
    print(model.jags)
  }
    #Rhat
    # summary(model.jags$BUGSoutput$summary[,8])
  # summary(model.jags$BUGSoutput$summary[,9])




  ### ----- 3.2. Fitting the model with INLA --- ####
  cat(paste0("n = ", n, " -----> Fitting using INLA \n"))
  if(file.exists(paste0("model_inla_", n,".RDS"))){
    model.inla <- readRDS(paste0("model_inla_", n,".RDS"))
    if(n< 1000)
    {
      simulation <- readRDS("simulation1_50-500.RDS")
      t_inla <- simulation[[paste0("n",n)]]$times[2]
    }else{
      simulation <- readRDS("simulation1_1000-10000.RDS")
      t_inla <- simulation[[paste0("n",n)]]$times[2]
    }
  }else{
    t <- proc.time() # Measure the time
    model.inla <- dirinlareg( formula  = y ~ 1 | 1 | 1 | 1  ,
                              y        = y,
                              data.cov = V,
                              prec     = 0.0001,
                              verbose  = TRUE,
                              sim       = 4000)

    t_inla <- proc.time()-t    # Stop the time
    t_inla <- t_inla[3]
    summary(model.inla)
  }
  ### ----- 3.3. Fitting the model with long jags --- ####
  cat(paste0("n = ", n, " -----> Fitting using long JAGS \n"))
  if(file.exists(paste0("model_jags_long_", n,".RDS"))){
    model.jags.2 <- readRDS(paste0("model_jags_long_", n,".RDS"))
    if(n< 1000)
    {
      simulation <- readRDS("simulation1_50-500.RDS")
      t_jags_2 <- simulation[[paste0("n",n)]]$times[3]
    }else{
      simulation <- readRDS("simulation1_1000-10000.RDS")
      t_jags_2 <-simulation[[paste0("n",n)]]$times[3]
    }
  }else{
    ## MCMC configuration
    ni <- 1000000
    nt <- 5
    nb <- 100000
    nc <- 3


    ## Data set
    data_jags <- list(y = y,
                      N = dim(y)[1],
                      d = d)

    ## Initial values
    inits <- function(){list(beta0 = rnorm(d, 0, 1))}

    ## Parameters of interest
    parameters <- c('beta0', 'alpha[1,1]', 'alpha[1,2]', 'alpha[1,3]', 'alpha[1,4]', 'mu')

    cat("
    model {
    #model
    for (i in 1:N){
    y[i, ] ~ ddirch(alpha[i,])

    log(alpha[i,1]) <- beta0[1]
    log(alpha[i,2]) <- beta0[2]
    log(alpha[i,3]) <- beta0[3]
    log(alpha[i,4]) <- beta0[4]

    }
    alpha0 <- alpha[1,1] + alpha[1,2] + alpha[1,3] + alpha[1,4]
    mu[1] <- alpha[1,1]/alpha0
    mu[2] <- alpha[1,2]/alpha0
    mu[3] <- alpha[1,3]/alpha0
    mu[4] <- alpha[1,4]/alpha0
    #priors
    for (c in 1:d)
    {
    beta0[c]  ~ dnorm(0, 0.0001)
    }
    }", file="model.jags" )



    ## Call jags
    t <- proc.time() # Measure the time
    model.jags.2 <- jags(data_jags,
                         inits,
                         parameters,
                         "model.jags",
                         n.chains          = nc,
                         n.thin            = nt,
                         n.iter            = ni,
                         n.burnin          = nb,
                         working.directory = getwd()) #
    t_jags_2<-proc.time()-t    # Stop the time
    t_jags_2 <- t_jags_2[3]
  }
  ### ----- 3.4. Saving models --- ####
  saveRDS(file = paste0("model_jags_", n,".RDS"), model.jags)
  saveRDS(file = paste0("model_jags_long_", n, ".RDS"), model.jags.2)
  saveRDS(file = paste0("model_inla_", n, ".RDS"), model.inla)

  t_inla
  t_jags
  t_jags_2


  ### --- 4. Comparing methodologies --- ####
  cat(paste0("n = ", n, " -----> Comparing methodologies \n"))

  ### ----- 4.1. Computational times --- ####
  times <- c(t_jags, t_inla, t_jags_2)

  ### ----- 4.2. (E(INLA) - E(JAGS2))/SD(JAGS2) and variance ratios --- ####
  ratio1_beta0 <- ratio2_beta0 <- ratio1_mu <- ratio2_mu <-  numeric()
  for (i in 1:4)
  {
    mean_jags_2 <- mean(model.jags.2$BUGSoutput$sims.list$beta0[,i])
    sd_jags_2 <- sd(model.jags.2$BUGSoutput$sims.list$beta0[,i])
    # mean_jags_1 <- mean(model.jags$BUGSoutput$sims.list$beta0[,1])
    # sd_jags_1 <- sd(model.jags$BUGSoutput$sims.list$beta0[,1])
    mean_inla <- model.inla$summary_fixed[[i]]$mean

    ratio1_beta0 <- c(ratio1_beta0, c(mean_inla - mean_jags_2)/sd_jags_2)
    ratio2_beta0 <- c(ratio2_beta0, sd(inla.rmarginal(10000, model.inla$marginals_fixed[[i]]$intercept))^2/sd_jags_2^2)
  }

  for (i in 1:4)
  {
    mean_jags_2 <- mean(model.jags.2$BUGSoutput$sims.list$mu[,i])
    sd_jags_2 <- sd(model.jags.2$BUGSoutput$sims.list$mu[,i])
    # mean_jags_1 <- mean(model.jags$BUGSoutput$sims.list$beta0[,1])
    # sd_jags_1 <- sd(model.jags$BUGSoutput$sims.list$beta0[,1])
    mean_inla <- model.inla$summary_means[[i]][1, "mean"]

    ratio1_mu <- c(ratio1_mu, c(mean_inla - mean_jags_2)/sd_jags_2)
    ratio2_mu <- c(ratio2_mu, sd(model.inla$marginals_means[[i]])^2/sd_jags_2^2)
  }

  ### ----- 4.3. Mean and sd of the posterior distributions --- ####
  ### Intercepts
  result <- numeric()
  for(i in 1:4)
  {
    result <- rbind(result,
                    t(matrix(c(model.jags$BUGSoutput$summary[paste0("beta0[", i,"]"), c("mean", "sd")],
                           model.inla$summary_fixed[[i]][,c("mean", "sd")],
                           model.jags.2$BUGSoutput$summary[paste0("beta0[", i,"]"), c("mean", "sd")]))))
  }
  rownames(result) <- paste0("beta0", 1:4)
  colnames(result) <- c(paste0("JAGS", c("_mean", "_sd")),
                        paste0("INLA", c("_mean", "_sd")),
                        paste0("LONG_JAGS", c("_mean", "_sd")))


  ### Mus
  result_mus <- numeric()
  for(i in 1:4)
  {
    result_mus <- rbind(result_mus,
                    t(matrix(c(model.jags$BUGSoutput$summary[paste0("mu[", i,"]"), c("mean", "sd")],
                               model.inla$summary_means[[i]][1,c("mean", "sd")],
                               model.jags.2$BUGSoutput$summary[paste0("mu[", i,"]"), c("mean", "sd")]))))
  }
  rownames(result_mus) <- paste0("mu", 1:4)
  colnames(result_mus) <- c(paste0("JAGS", c("_mean", "_sd")),
                        paste0("INLA", c("_mean", "_sd")),
                        paste0("LONG_JAGS", c("_mean", "_sd")))
  ### --- 5. Plotting --- ####
  cat(paste0("n = ", n, " -----> Plotting \n"))

  ### ----- 5.1. intercepts --- ####
  ## Intercept
  p2 <- list()
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
                  cbind(as.data.frame(model.inla$marginals_fixed[[i]][[1]]), group = 2),
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
      xlab(expression(beta[0])) + #xlab
      ylab(beta0) #ylab


    #Frequentist approach
    p2[[i]] <- p2[[i]] + geom_vline(xintercept = x[seq(1,4, by=1)][i])



    ### --- legend --- ###
    p2[[i]]<- p2[[i]] + theme(legend.position   = c(0.2, 0.8),
                              legend.title      = element_blank(),
                              legend.background = element_rect(colour = "gray"),
                              legend.key        = element_rect(colour = "white", fill="white"),
                              legend.key.size   = unit(0.5, "cm")) +
      theme(legend.text = element_text(size = 9)) +
      # scale_fill_manual(labels=c("R-JAGS", "R-INLA", "long R-JAGS"),
      #                   values = c("darkgreen", "red4", "blue4" )) +
      scale_colour_manual (
        values= c("darkgreen", "red4", "blue4")) +
      scale_linetype_manual(labels=c("R-JAGS", "dirinla", "long R-JAGS"),
                            values=c("dotted", "twodash", "solid" ))

    if(i!=1)
    {
      p2[[i]] <- p2[[i]] + theme(legend.position="none")
    }

    p2[[i]] <- p2[[i]] + ggtitle(paste0("Category ", i)) +
      theme(
        plot.title = element_text(color = "black",
                                  size  = 15,
                                  face  = "bold.italic",
                                  hjust = 0.5))
  }


  # pdf("example_simulation1_beta0_50.pdf", width = 18, height = 4)
  # gridExtra::grid.arrange(p2[[1]], p2[[2]], p2[[3]], p2[[4]], ncol = 4)
  # dev.off()






  ### ----- 5.2. Means --- ####
  # ### ------- 5.2.1. Plot --- ####
  # pdf(paste0("example_simulation_means_", n, ".pdf"), width = 20, height = 5)
  # mu_not <- expression(mu)
  # p_mu_not <- expression(paste("p(", mu, "|", "y)"))
  # par(mfrow=c(1,4))
  # for (i in 1:4)
  # {
  #   # hist(model.inla$marginals_means[[i]],
  #   #      col="red",
  #   #      lwd=2, freq = FALSE, breaks = 50)
  #
  #   plot(density(model.inla$marginals_means[[i]], adjust = 2.5),
  #        col = "green",
  #        lwd = 4, type = "l")
  #
  #   lines(density(model.jags$BUGSoutput$sims.list$mu[,i], adjust = 2.5),
  #         col = "orange",
  #         lwd  = 4)
  #   # ylim= c(0, max(model.inla$marginals_means[[i]][[1]][,2],
  #   #                density(model.jags.2$BUGSoutput$sims.list$mu[,1,i])$y)),
  #   # xlab = mu_not,
  #   # ylab = p_mu_not,
  #   # main = mu_not
  #
  #
  #
  #
  #
  #   lines(density(model.jags.2$BUGSoutput$sims.list$mu[,i], adjust = 2.5),
  #         col = "blue")
  #
  #   # hist(model.jags.2$BUGSoutput$sims.list$mu[,i],
  #   #      col = "blue", breaks = 20, freq = FALSE)
  #   abline(v = mus[i],
  #          col = "black",
  #          lwd = 2)
  #   #
  #   #   legend("topright", legend=c("R-jags", "dirinla", "R-long-jags"),
  #   #          col = c("orange", "red", "blue"),
  #   #          lty = 1,
  #   #          lwd = 2)
  # }
  # dev.off()
  #
  #
  #
  #
  #

  ### ------- 5.2.2. ggplot --- ####
  ## Intercept
  p3 <- list()
  mu_not <- expression(paste("p(", mu, "|", "y)"))
  if(n>=500)
  {
    adjust_inla <- 5
  }else{
    adjust_inla <- 3
  }

  for (i in 1:length(model.inla$marginals_fixed))
  {
    #jags1
    dens <- density(model.jags$BUGSoutput$sims.list$mu[,i], adjust = 2)
    dens <- as.data.frame(cbind(dens$x, dens$y))
    colnames(dens) <- c("x", "y")

    #jags2
    dens2 <- density(model.jags.2$BUGSoutput$sims.list$mu[,i], adjust = 2)
    dens2 <- as.data.frame(cbind(dens2$x, dens2$y))
    colnames(dens2) <- c("x", "y")

    #inla
    #jags2
    dens3 <- density(model.inla$marginals_means[[i]], adjust = adjust_inla)
    dens3 <- as.data.frame(cbind(dens3$x, dens3$y))
    colnames(dens3) <- c("x", "y")

    #Data combining jags (1) and inla (2)
    dens <- rbind(cbind(dens, group = 1),
                  cbind(dens3, group = 2),
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
      xlab(expression(mu)) + #xlab
      ylab(mu_not) #ylab


    #Frequentist approach
    p3[[i]] <- p3[[i]] + geom_vline(xintercept = mus[seq(1,4, by=1)][i])



    ### --- legend --- ###
    p3[[i]]<- p3[[i]] + theme(legend.position   = c(0.2, 0.8),
                              legend.title      = element_blank(),
                              legend.background = element_rect(colour = "gray"),
                              legend.key        = element_rect(colour = "white", fill="white"),
                              legend.key.size   = unit(0.5, "cm")) +
      theme(legend.text = element_text(size = 9)) +
      # scale_fill_manual(labels=c("R-JAGS", "R-INLA", "long R-JAGS"),
      #                   values = c("darkgreen", "red4", "blue4" )) +
      scale_colour_manual (
        values= c("darkgreen", "red4", "blue4")) +
      scale_linetype_manual(labels=c("R-JAGS", "dirinla", "long R-JAGS"),
                            values=c("dotted",  "twodash", "solid"))

    if(i!=5)
    {
      p3[[i]] <- p3[[i]] + theme(legend.position="none")
    }

    # p3[[i]] <- p3[[i]] + ggtitle(paste0("Category ", i)) +
    #   theme(
    #     plot.title = element_text(color = "black",
    #                               size  = 15,
    #                               face  = "bold.italic",
    #                               hjust = 0.5))
  }


  # pdf(paste0("example_simulation1_mus_", n, ".pdf"), width = 18, height = 4)
  # gridExtra::grid.arrange(p3[[1]], p3[[2]], p3[[3]], p3[[4]], ncol = 4)
  # dev.off()

  pdf(paste0("examples_simulation1_beta0_mus_", n, ".pdf"), width = 15, height = 6)
  gridExtra::grid.arrange(p2[[1]], p2[[2]], p2[[3]], p2[[4]],
               p3[[1]], p3[[2]], p3[[3]], p3[[4]], ncol = 4)
  dev.off()


  list(times = times,
       intercepts       = result,
       mus              = result_mus,
       ratio1_intercept = ratio1_beta0,
       ratio2_intercept = ratio2_beta0,
       ratio1_mu        = ratio1_mu,
       ratio2_mu        = ratio2_mu)
}


### --- 3. Calling the function --- ####
n <- c(50, 100, 500)
# a <- parallel::mclapply(n, simulations_just_intercepts,
#                         mc.cores = 3)

a <- lapply(n, simulations_just_intercepts)

names(a) <- paste0("n", n)
saveRDS(b, file = "simulation1_50-500.RDS")
a <- readRDS(file = "simulation1_50-500.RDS")


n <- c(1000, 10000)
# a <- parallel::mclapply(n, simulations_just_intercepts,
#                         mc.cores = 2)
a <- lapply(n, simulations_just_intercepts)

names(a) <- paste0("n", n)
saveRDS(a, file = "simulation1_1000-10000.RDS")
a <- readRDS(file = "simulation1_1000-10000.RDS")

### --- 4. Computing also ratios for shortJAGS --- ####
ratios_jags <- function(n)
{
  print(n)
  model.jags <- readRDS(paste0("model_jags_", n,".RDS"))
  #model.inla <- readRDS(paste0("model_inla_", n,".RDS"))
  model.jags.2 <- readRDS(paste0("model_jags_long_", n,".RDS"))

  ratio1_beta0_jags <- ratio2_beta0_jags <- ratio1_mu_jags <- ratio2_mu_jags <- numeric()
  for (i in 1:4)
  {
    mean_jags_2 <- mean(model.jags.2$BUGSoutput$sims.list$beta0[,i])
    sd_jags_2 <- sd(model.jags.2$BUGSoutput$sims.list$beta0[,i])
    mean_jags_1 <- mean(model.jags$BUGSoutput$sims.list$beta0[,i])
    sd_jags_1 <- sd(model.jags$BUGSoutput$sims.list$beta0[,i])

    ratio1_beta0_jags <- c(ratio1_beta0_jags, c(mean_jags_1 - mean_jags_2)/sd_jags_2)
    ratio2_beta0_jags <- c(ratio2_beta0_jags, (sd_jags_1^2)/(sd_jags_2^2))
  }

  for (i in 1:4)
  {
    mean_jags_2 <- mean(model.jags.2$BUGSoutput$sims.list$mu[,i])
    sd_jags_2 <- sd(model.jags.2$BUGSoutput$sims.list$mu[,i])
    mean_jags_1 <- mean(model.jags$BUGSoutput$sims.list$mu[,i])
    sd_jags_1 <- sd(model.jags$BUGSoutput$sims.list$mu[,i])

    ratio1_mu_jags <- c(ratio1_mu_jags, c(mean_jags_1 - mean_jags_2)/sd_jags_2)
    ratio2_mu_jags <- c(ratio2_mu_jags, (sd_jags_1^2/sd_jags_2^2))
  }
  list(ratio1_beta0_jags = ratio1_beta0_jags,
       ratio2_beta0_jags = ratio2_beta0_jags,
       ratio1_mu_jags = ratio1_mu_jags,
       ratio2_mu_jags = ratio2_mu_jags)

}

n <- c(50, 100, 500, 1000, 10000)
res_ratios <- lapply(n, ratios_jags)
names(res_ratios) <- paste0("n", n)
res_ratios

saveRDS(res_ratios, file = "simulation1_ratios_jags.RDS")


### --- 5. Extracting tables for the paper --- ####
results1 <- readRDS(file = "simulation1_50-500.RDS")
results2 <- readRDS(file = "simulation1_1000-10000.RDS")
results <- c(results1, results2)

res_ratios <- readRDS("simulation1_ratios_jags.RDS")


#Computational times
result_time <- rbind(results$n50$times,
                     results$n100$times,
                     results$n500$times,
                     results$n1000$times,
                     results$n10000$times)
colnames(result_time) <- c("R-JAGS", "dirinla", "long R-JAGS")
rownames(result_time) <- paste0( c(50, 100, 500, 1000, 10000))
result_time

### ratios dirinla
result_ratio1 <- cbind(rbind(round(results$n50$ratio1_intercept, 4),
                             round(results$n100$ratio1_intercept, 4),
                             round(results$n500$ratio1_intercept, 4),
                             round(results$n1000$ratio1_intercept, 4),
                             round(results$n10000$ratio1_intercept, 4)),
                       rbind(round(results$n50$ratio1_mu, 4),
                             round(results$n100$ratio1_mu, 4),
                             round(results$n500$ratio1_mu, 4),
                             round(results$n1000$ratio1_mu, 4),
                             round(results$n10000$ratio1_mu, 4)))
colnames(result_ratio1) <- c(paste0("beta0", 1:4), paste0("mu", 1:4))
rownames(result_ratio1) <- paste0( c(50, 100, 500, 1000, 10000))

result_ratio2 <- cbind(rbind(round(sqrt(results$n50$ratio2_intercept), 4),
                             round(sqrt(results$n100$ratio2_intercept), 4),
                             round(sqrt(results$n500$ratio2_intercept), 4),
                             round(sqrt(results$n1000$ratio2_intercept), 4),
                             round(sqrt(results$n10000$ratio2_intercept), 4)),
                       rbind(round(sqrt(results$n50$ratio2_mu), 4),
                             round(sqrt(results$n100$ratio2_mu), 4),
                             round(sqrt(results$n500$ratio2_mu), 4),
                             round(sqrt(results$n1000$ratio2_mu), 4),
                             round(sqrt(results$n10000$ratio2_mu), 4)))
colnames(result_ratio2) <- c(paste0("beta0", 1:4), paste0("mu", 1:4))
rownames(result_ratio2) <- paste0( c(50, 100, 500, 1000, 10000))

### ratios short jags
result_ratio1_jags <- cbind(rbind(round(res_ratios$n50$ratio1_beta0_jags, 4),
                                  round(res_ratios$n100$ratio1_beta0_jags, 4),
                                  round(res_ratios$n500$ratio1_beta0_jags, 4),
                                  round(res_ratios$n1000$ratio1_beta0_jags, 4),
                                  round(res_ratios$n10000$ratio1_beta0_jags, 4)),
                            rbind(round(res_ratios$n50$ratio1_mu_jags, 4),
                                  round(res_ratios$n100$ratio1_mu_jags, 4),
                                  round(res_ratios$n500$ratio1_mu_jags, 4),
                                  round(res_ratios$n1000$ratio1_mu_jags, 4),
                                  round(res_ratios$n10000$ratio1_mu_jags, 4)))
colnames(result_ratio1_jags) <- c(paste0("beta0", 1:4), paste0("mu", 1:4))
rownames(result_ratio1_jags) <- paste0( c(50, 100, 500, 1000, 10000))

result_ratio2_jags <- cbind(rbind(round(sqrt(res_ratios$n50$ratio2_beta0_jags), 4),
                                  round(sqrt(res_ratios$n100$ratio2_beta0_jags), 4),
                                  round(sqrt(res_ratios$n500$ratio2_beta0_jags), 4),
                                  round(sqrt(res_ratios$n1000$ratio2_beta0_jags), 4),
                                  round(sqrt(res_ratios$n10000$ratio2_beta0_jags), 4)),
                            rbind(round(sqrt(res_ratios$n50$ratio2_mu_jags), 4),
                                  round(sqrt(res_ratios$n100$ratio2_mu_jags), 4),
                                  round(sqrt(res_ratios$n500$ratio2_mu_jags), 4),
                                  round(sqrt(res_ratios$n1000$ratio2_mu_jags), 4),
                                  round(sqrt(res_ratios$n10000$ratio2_mu_jags), 4)))
colnames(result_ratio2_jags) <- c(paste0("beta0", 1:4), paste0("mu", 1:4))
rownames(result_ratio2_jags) <- paste0( c(50, 100, 500, 1000, 10000))







#Latex
library(xtable)
xtable(result_time, digits = 4)
xtable(result_ratio1, digits = 4)
xtable(result_ratio1_jags, digits = 4)

xtable(result_ratio2, digits = 4)
xtable(result_ratio2_jags, digits = 4)


### --- 6. Plotting results --- ####


### ----- 6.1. function to plot ratios --- ####
plot_ratios <- function(result_ratio_tot = result_ratio1_tot,
                        param = c("beta01", "beta02", "beta03",
                                  "beta04", "method", "N"),
                        names_param = c('beta[0][1]', 'beta[0][2]', 'beta[0][3]', 'beta[0][4]'),
                        exp1        = expression(ratio[1]),
                        int_plot    = 0,
                        label1      = "a")
{
  result_ratio_tot %>%
    dplyr::select_at(., param) %>%
    tidyr::pivot_longer(1:(length(param)-2), names_to = "betas") -> result_ratio_tot_betas

  #result_ratio_tot_betas$value <- abs(result_ratio_tot_betas$value)
  result_ratio_tot_betas$method <- ordered(result_ratio_tot_betas$method, levels = c("R-JAGS", "dirinla"))
  #result_ratio_tot_betas$N <- factor(paste0("N = ", result_ratio_tot_betas$N), levels = c("N = 50", "N = 100", "N = 500", "N = 1000", "N = 10000"))
  #result_ratio_tot_betas$N <- factor(result_ratio_tot_betas$N)

  result_ratio_tot_betas$newbetas <- as.factor(result_ratio_tot_betas$betas)
  levels(result_ratio_tot_betas$newbetas) <- names_param
  #result_ratio_tot_betas$N <- as.numeric(result_ratio_tot_betas$N)

  result_ratio_tot_betas %>%
    ggplot(data = .) +
    geom_point(aes(x = N, y = value, col = method, shape = method), size = 3) +
    theme_bw() +
    # geom_line(aes(x = log(N), y = value, col = method, shape = method)) +
    theme(legend.position   = c(0.2, 0.92),
          legend.title      = element_blank(),
          legend.background = element_rect(colour = "gray"),
          legend.key        = element_rect(colour = "white", fill="white"),
          legend.key.size   = unit(0.5, "cm"),
          axis.text         = element_text(size=12)) +
    facet_wrap('newbetas', ncol = 4, labeller = label_parsed) +
    theme(strip.text=element_text(face='bold', size=12, color='black'),
          strip.background=element_rect(fill='white')) +
    theme(legend.text = element_text(size = 9)) +
    # scale_fill_manual(labels=c("R-JAGS", "R-INLA", "long R-JAGS"),
    #                   values = c("darkgreen", "red4", "blue4" )) +
    scale_shape_manual(values=c(16, 17)) +
    scale_colour_manual (
      values= c("darkgreen", "red4")) +
    geom_hline(yintercept = int_plot,  linetype="dashed") +
    scale_x_log10() +
    ylab(exp1) -> ratio_beta#xlab
  ratio_beta#xlab
}


### ----- 6.2. ratio 1 --- ####
result_ratio1 <- as.data.frame(result_ratio1)
result_ratio1_jags <- as.data.frame(result_ratio1_jags)

result_ratio1_tot <- rbind(
  cbind(N = rownames(result_ratio1_jags) %>% as.numeric(),
        result_ratio1_jags, method = "R-JAGS"),
  cbind(N = rownames(result_ratio1) %>% as.numeric(),
        result_ratio1, method = "dirinla")) %>%
  data.frame()

ratios1 <- plot_ratios(result_ratio_tot = result_ratio1_tot,
            param = c("beta01", "beta02", "beta03", "beta04", "mu1", "mu2", "mu3",
                      "mu4", "method", "N"),
            names_param = c('beta[0][1]', 'beta[0][2]', 'beta[0][3]', 'beta[0][4]',
                            'mu[1]', 'mu[2]', 'mu[3]', 'mu[4]'),
            exp1 = expression(ratio[1]))

pdf(paste0("simulation1_beta0_ratio1", ".pdf"), width = 15, height = 6)
# gridExtra::grid.arrange(ratio1_beta,
#                         ratio1_mu, ncol = 1)
ratios1
dev.off()

### ----- 6.2. ratio 2 --- ####
result_ratio2 <- as.data.frame(result_ratio2)
result_ratio2_jags <- as.data.frame(result_ratio2_jags)

result_ratio2_tot <- rbind(
  cbind(N = rownames(result_ratio2_jags) %>% as.numeric(),
        result_ratio2_jags, method = "R-JAGS"),
  cbind(N = rownames(result_ratio2) %>% as.numeric(),
        result_ratio2, method = "dirinla")) %>%
  data.frame()


ratios2 <- plot_ratios(result_ratio_tot = result_ratio2_tot,
                       param = c("beta01", "beta02", "beta03", "beta04", "mu1", "mu2", "mu3",
                                 "mu4", "method", "N"),
                       names_param = c('beta[0][1]', 'beta[0][2]', 'beta[0][3]', 'beta[0][4]',
                                       'mu[1]', 'mu[2]', 'mu[3]', 'mu[4]'),
                       exp1 = expression(ratio[2]),
                       int_plot = 1)

pdf(paste0("simulation1_beta0_ratio2", ".pdf"), width = 15, height = 6)
ratios2
dev.off()

### ----- 6.3. Putting both together --- ####
pdf(paste0("simulation1_ratios", ".pdf"), width = 15, height = 12)
plot_grid(ratios1, ratios2, labels=c('a','b'), ncol = 1)
dev.off()




### ----- 6.4. Times --- ####
result_time2 <- cbind(result_time, N = as.numeric(rownames(result_time))) %>% as.data.frame(.)
result_time2 %>%
  tidyr::pivot_longer(1:3, names_to = "method") -> result_time2
result_time2$method <- ordered(result_time2$method, levels = c("R-JAGS", "dirinla", "long R-JAGS"))
#result_time2$N <- factor(result_time2$N)


times2 <- result_time2 %>%
  ggplot(data = .) +
  geom_point(aes(x = N, y = value, col = method, shape = method), size = 3) +
  theme_bw() +
  theme(legend.position   = c(0.15, 0.80),
        legend.title      = element_blank(),
        legend.background = element_rect(colour = "gray"),
        legend.key        = element_rect(colour = "white", fill="white"),
        legend.key.size   = unit(0.5, "cm"),
        axis.text         = element_text(size=12)) +
  theme(legend.text = element_text(size = 9)) +
  # scale_fill_manual(labels=c("R-JAGS", "R-INLA", "long R-JAGS"),
  #                   values = c("darkgreen", "red4", "blue4" )) +
  scale_shape_manual(values=c(16, 17, 18)) +
  scale_colour_manual (
    values= c("darkgreen", "red4", "blue4")) +
  ylab("Time (sec)") +
  scale_x_log10() +
  scale_y_log10()

pdf(paste0("simulation1_times", ".pdf"), width = 6, height = 4)
times2
dev.off()
