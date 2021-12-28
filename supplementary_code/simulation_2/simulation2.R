### In this script simulations with N=50, 100, 5000, 1000, 10000 are conducted #
# in order to check    ###
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

#Generating latex table
library(xtable)

### --- 2. Function for simulation --- ####
simulations_with_slopes <- function(n)
{
  ### --- 2. Simulation data --- ####
  cat("n = ", n, " -----> Simulating data \n")
  set.seed(1000)

  #Covariates
  V <- as.data.frame(matrix(runif((10)*n, 0, 1), ncol=10))
  names(V) <- paste0('v', 1:(10))

  # Formula that we want to fit
  formula <- y ~ 1 + v1 | 1 + v2 | 1 + v3 | 1 + v4

  names_cat <- formula_list(formula)

  # Parameters to fit
  # x <- c(-1.5, 1, -3, 1.5,
  #        2, -3 , -1, 5)


  x <- c(-1.5, 2,
         1, -3,
         -3, -1,
         1.5, 5)


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



  ### --- 3. Comparing posterior distributions. jags vs INLA --- ####
  ### ----- 3.1. Fitting the model with jags --- ####
  cat(paste0("n = ", n, " -----> Fitting using SHORT JAGS \n"))

  if(file.exists(paste0("model_jags_", n,".RDS"))){
    model.jags <- readRDS(paste0("model_jags_", n,".RDS"))
    if(n < 1000)
    {
      simulation <- readRDS("simulation2_50-500.RDS")
      t_jags <- simulation[[paste0("n",n)]]$times[1]
    }else{
      simulation <- readRDS("simulation2_1000-10000.RDS")
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
      simulation <- readRDS("simulation2_50-500.RDS")
      t_inla <- simulation[[paste0("n",n)]]$times[2]
    }else{
      simulation <- readRDS("simulation2_1000-10000.RDS")
      t_inla <- simulation[[paste0("n",n)]]$times[2]
    }
  }else{
    t <- proc.time() # Measure the time
    model.inla <- dirinlareg( formula  = y ~ 1 + v1 | 1 + v2 | 1 + v3 | 1 + v4  ,
                              y        = y,
                              data.cov = V,
                              prec     = 0.0001,
                              verbose  = TRUE)

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
      simulation <- readRDS("simulation2_50-500.RDS")
      t_jags_2 <- simulation[[paste0("n",n)]]$times[3]
    }else{
      simulation <- readRDS("simulation2_1000-10000.RDS")
      t_jags_2 <-simulation[[paste0("n",n)]]$times[3]
    }
  }else{
    ## MCMC configuration
    #  ni <- 10000
    #  nb <- 1000
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
    p1[[i]] <- p1[[i]] + geom_vline(xintercept = x[seq(1,8, by=2)][i])



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
    p2[[i]] <- p2[[i]] + geom_vline(xintercept = x[seq(2,8, by=2)][i])



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
}


### --- 3. Calling the function --- ####
n <- c(50, 100, 500)
# a <- parallel::mclapply(n, simulations_with_slopes,
#                         mc.cores = 3)

a <- lapply(n, simulations_with_slopes)

names(a) <- paste0("n", n)
saveRDS(a, file = "simulation2_50-500.RDS")
a <- readRDS(file = "simulation2_50-500.RDS")

n <- c(1000, 10000)
# a <- parallel::mclapply(n, simulations_with_slopes,
#                         mc.cores = 2)
a <- lapply(n, simulations_with_slopes)

names(a) <- paste0("n", n)
saveRDS(a, file = "simulation2_1000-10000.RDS")
b <- readRDS(file = "simulation2_1000-10000.RDS")
results <- c(a,b)

### --- 4. Computing also ratios for shortJAGS --- ####
ratios_jags <- function(n)
{
  print(n)
  model.jags <- readRDS(paste0("model_jags_", n,".RDS"))
  model.inla <- readRDS(paste0("model_inla_", n,".RDS"))
  model.jags.2 <- readRDS(paste0("model_jags_long_", n,".RDS"))

  ratio1_beta0_jags <- ratio2_beta0_jags <- ratio1_beta1_jags <- ratio2_beta1_jags <- numeric()
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
    mean_jags_2 <- mean(model.jags.2$BUGSoutput$sims.list$beta1[,i])
    sd_jags_2 <- sd(model.jags.2$BUGSoutput$sims.list$beta1[,i])
    mean_jags_1 <- mean(model.jags$BUGSoutput$sims.list$beta1[,i])
    sd_jags_1 <- sd(model.jags$BUGSoutput$sims.list$beta1[,i])

    ratio1_beta1_jags <- c(ratio1_beta1_jags, c(mean_jags_1 - mean_jags_2)/sd_jags_2)
    ratio2_beta1_jags <- c(ratio2_beta1_jags, (sd_jags_1^2/sd_jags_2^2))
  }
  list(ratio1_beta0_jags = ratio1_beta0_jags,
       ratio2_beta0_jags = ratio2_beta0_jags,
       ratio1_beta1_jags = ratio1_beta1_jags,
       ratio2_beta1_jags = ratio2_beta1_jags)
}

n <- c(50, 100, 500, 1000, 10000)
res_ratios <- lapply(n, ratios_jags)
names(res_ratios) <- paste0("n", n)
res_ratios




### --- 5. Extracting tables for the paper --- ####
results1 <- readRDS(file = "simulation2_50-500.RDS")
results2 <- readRDS(file = "simulation2_1000-10000.RDS")
results <- c(results1, results2)
results$n50$times
results$n50$intercepts
results$n50$times
results$n50$times
results$n50$ratio1_intercepts
results$n50$ratio1_slopes

sqrt(results$n50$ratio2_intercepts)
sqrt(results$n50$ratio2_slopes)

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
result_ratio1 <- cbind(rbind(round(results$n50$ratio1_intercepts, 4),
                             round(results$n100$ratio1_intercepts, 4),
                             round(results$n500$ratio1_intercepts, 4),
                             round(results$n1000$ratio1_intercepts, 4),
                             round(results$n10000$ratio1_intercepts, 4)),
                       rbind(round(results$n50$ratio1_slopes, 4),
                             round(results$n100$ratio1_slopes, 4),
                             round(results$n500$ratio1_slopes, 4),
                             round(results$n1000$ratio1_slopes, 4),
                             round(results$n10000$ratio1_slopes, 4)))
colnames(result_ratio1) <- c(paste0("beta0", 1:4), paste0("beta1", 1:4))
rownames(result_ratio1) <- paste0( c(50, 100, 500, 1000, 10000))

result_ratio2 <- cbind(rbind(round(sqrt(results$n50$ratio2_intercepts), 4),
                             round(sqrt(results$n100$ratio2_intercepts), 4),
                             round(sqrt(results$n500$ratio2_intercepts), 4),
                             round(sqrt(results$n1000$ratio2_intercepts), 4),
                             round(sqrt(results$n10000$ratio2_intercepts), 4)),
                       rbind(round(sqrt(results$n50$ratio2_slopes), 4),
                             round(sqrt(results$n100$ratio2_slopes), 4),
                             round(sqrt(results$n500$ratio2_slopes), 4),
                             round(sqrt(results$n1000$ratio2_slopes), 4),
                             round(sqrt(results$n10000$ratio2_slopes), 4)))
colnames(result_ratio2) <- c(paste0("beta0", 1:4), paste0("beta1", 1:4))
rownames(result_ratio2) <- paste0( c(50, 100, 500, 1000, 10000))

### ratios short jags
result_ratio1_jags <- cbind(rbind(round(res_ratios$n50$ratio1_beta0_jags, 4),
                             round(res_ratios$n100$ratio1_beta0_jags, 4),
                             round(res_ratios$n500$ratio1_beta0_jags, 4),
                             round(res_ratios$n1000$ratio1_beta0_jags, 4),
                             round(res_ratios$n10000$ratio1_beta0_jags, 4)),
                       rbind(round(res_ratios$n50$ratio1_beta1_jags, 4),
                             round(res_ratios$n100$ratio1_beta1_jags, 4),
                             round(res_ratios$n500$ratio1_beta1_jags, 4),
                             round(res_ratios$n1000$ratio1_beta1_jags, 4),
                             round(res_ratios$n10000$ratio1_beta1_jags, 4)))
colnames(result_ratio1_jags) <- c(paste0("beta0", 1:4), paste0("beta1", 1:4))
rownames(result_ratio1_jags) <- paste0( c(50, 100, 500, 1000, 10000))

result_ratio2_jags <- cbind(rbind(round(sqrt(res_ratios$n50$ratio2_beta0_jags), 4),
                             round(sqrt(res_ratios$n100$ratio2_beta0_jags), 4),
                             round(sqrt(res_ratios$n500$ratio2_beta0_jags), 4),
                             round(sqrt(res_ratios$n1000$ratio2_beta0_jags), 4),
                             round(sqrt(res_ratios$n10000$ratio2_beta0_jags), 4)),
                       rbind(round(sqrt(res_ratios$n50$ratio2_beta1_jags), 4),
                             round(sqrt(res_ratios$n100$ratio2_beta1_jags), 4),
                             round(sqrt(res_ratios$n500$ratio2_beta1_jags), 4),
                             round(sqrt(res_ratios$n1000$ratio2_beta1_jags), 4),
                             round(sqrt(res_ratios$n10000$ratio2_beta1_jags), 4)))
colnames(result_ratio2_jags) <- c(paste0("beta0", 1:4), paste0("beta1", 1:4))
rownames(result_ratio2_jags) <- paste0( c(50, 100, 500, 1000, 10000))







#Latex
library(xtable)
xtable(result_time, digits = 4)
xtable(result_ratio1, digits = 4)
xtable(result_ratio1_jags, digits = 4)

xtable(result_ratio2, digits = 4)
xtable(result_ratio2_jags, digits = 4)



plot(result_time[,1])
lines(result_time[,2], col = "red")







