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
setwd("~/GIT1/dirinla/supplementary_code/simulation_1")

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
  A_construct <- data_stack_construct$A
  eta <- A_construct %*% x
  alpha <- exp(eta)
  alpha <- matrix(alpha,
                  ncol  = d,
                  byrow = TRUE)
  y_o <- rdirichlet(n, alpha)
  colnames(y_o) <- paste0("Category", 1:d)


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
                              verbose  = FALSE,
                              sim       = 5000)

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
    # ni <- 10000
    # nb <- 100
    ni <- 1000000
    nt <- 5
    nb <- 100000
    nc <- 3
    # ni <- 10000
    # nb <- 100

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
#n <- c(50, 100, 500, 1000, 5000, 10000)
n <- c(50, 100, 500)

a <- parallel::mclapply(n, simulations_just_intercepts,
                        mc.cores = 3)
a <- lapply(n, simulations_just_intercepts)
names(a) <- paste0("n", n)
saveRDS(b, file = "simulation1_50-500.RDS")
a <- readRDS(file = "simulation1_50-500.RDS")
a$n50$times


### --- 4. Extracting tables for the paper --- ####
results <- readRDS(file = "simulation1_50-500.RDS")
results$n50$times
results$n50$intercept
results$n50$times
results$n50$times
results$n50$ratio1_intercept
results$n50$ratio1_mu

sqrt(results$n50$ratio2_intercept)
sqrt(results$n50$ratio2_mu)

#Computational times
result_time <- rbind(results$n50$times,
                     results$n100$times,
                     results$n500$times)
colnames(result_time) <- c("R-JAGS", "dirinla", "long R-JAGS")
rownames(result_time) <- paste0( c(50, 100, 500))
result_time

result_ratio1 <- cbind(rbind(round(results$n50$ratio1_intercept, 4),
                             round(results$n100$ratio1_intercept, 4),
                             round(results$n500$ratio1_intercept, 4)),
                       rbind(round(results$n50$ratio1_mu, 4),
                             round(results$n100$ratio1_mu, 4),
                             round(results$n500$ratio1_mu, 4)))
colnames(result_ratio1) <- c(paste0("beta0", 1:4), paste0("beta1", 1:4))
rownames(result_ratio1) <- paste0( c(50, 100, 500))

result_ratio2 <- cbind(rbind(round(sqrt(results$n50$ratio2_intercept), 4),
                             round(sqrt(results$n100$ratio2_intercept), 4),
                             round(sqrt(results$n500$ratio2_intercept), 4)),
                       rbind(round(sqrt(results$n50$ratio2_mu), 4),
                             round(sqrt(results$n100$ratio2_mu), 4),
                             round(sqrt(results$n500$ratio2_mu), 4)))
colnames(result_ratio2) <- c(paste0("beta0", 1:4), paste0("mu", 1:4))
rownames(result_ratio2) <- paste0( c(50, 100, 500))

result_ratio2

#Latex
xtable(result_time, digits = 4)
xtable(result_ratio1, digits = 4)
xtable(result_ratio2, digits = 4)

