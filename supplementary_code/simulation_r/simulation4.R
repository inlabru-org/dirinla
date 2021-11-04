### In this script simulations with N=50, 100, 5000, 1000, 10000 are conducted #
# in order to check    ###
# how the package dirinla works. We compare with r-jags with enough amount of #
# simulations to guarantee convergence of the method, and with long r-jags,   #
# which has a large amount of simulations. The model that we try to fit is :  #
# --- Y \sim Dirich(exp(eta_1), exp(eta_2), ... , exp(eta_4))                 #
# --- --- eta_1 = X1 beta1 + f(iid1, model = "iid"),               #
# --- --- eta_2 = X2 beta2 + f(iid1, model = "iid"),               #
# --- --- eta_3 = X3 beta3 + f(iid2, model = "iid"),               #
# --- --- eta_4 = X4 beta4 + f(iid2, model = "iid"),               #
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
library(parallel)

### --- 2. Function for simulation --- ####
simulations_with_slopes_iid <- function(n, levels_factor = NA)
{
  ### --- 2. Simulation data --- ####
  cat("n = ", n, " - levels_factor = ", levels_factor," -----> Simulating data \n")

  set.seed(100)
  if(is.na(levels_factor)){
    levels_factor <- n
  }
  cat_elem <- n/levels_factor
  cat(paste0(n, "-", levels_factor, "\n"))
  #Covariates
  V <- as.data.frame(matrix(runif((10)*n, -1, 1), ncol=10))
  #V <- as.data.frame(matrix(rnorm((10)*n, 0, 1), ncol=10))
  names(V) <- paste0('v', 1:(10))

  ### 4 random effects
  iid1 <- iid2  <- rep(1:levels_factor, rep(n/levels_factor, levels_factor))
  #Desorder index 3
  # pos <- sample(1:length(iid3))
  # iid3 <- iid3[pos]

  V <- cbind(V, iid1, iid2)

  # Formula that we want to fit
  formula <- y ~ -1 + v1 + f(iid1, model = 'iid') |
    -1 + v2 + f(iid1, model = 'iid') |
    -1 + v3 + f(iid2, model = 'iid') |
    -1 + v4 + f(iid2, model = 'iid')
  names_cat <- formula_list(formula)

  x <- c(-1.5, 2,
         1, -3)

  #random effect
  prec_w <- c(4, 9)
  (sd_w <- 1/sqrt(prec_w))

  w1 <- rnorm(levels_factor, sd = sqrt(1/prec_w[1])) %>% rep(., rep(n/levels_factor, levels_factor))
  w2 <- w1
  w3 <- rnorm(levels_factor, sd = sqrt(1/prec_w[2])) %>% rep(., rep(n/levels_factor, levels_factor))
  w4 <- w3


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



  ### --- 3. Comparing posterior distributions. jags vs INLA --- ####
  ### ----- 3.1. Fitting the model with jags --- ####
  cat(paste0("n = ", n, " - levels_factor = ", levels_factor," -----> Fitting using SHORT JAGS \n"))

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
    ni <- 20000
    nt <- 5
    #nb <- 20
    nb <- 2000
    nc <- 3
#
    # ni <- 1000
    # nb <- 100

    x_pc <- seq(0.01, 20, length.out = 1000)
    y_pc <- inla.pc.dprec(x_pc, u = 10, alpha = 0.01, log = FALSE)
    prec_pc <- list(x_pc = x_pc, y_pc = y_pc)

    ## Data set
    table(V$iid1) %>%length() -> niv_length1
    table(V$iid2) %>%length() -> niv_length2
    data_jags <- list(y = y,
                      N = dim(y)[1],
                      d = d,
                      V = V,
                      niv1 = V$iid1,
                      niv2 = V$iid2,
                      niv_length1 = niv_length1,
                      niv_length2 = niv_length2,
                      x_pc = prec_pc$x_pc,
                      y_pc = prec_pc$y_pc)

    ## Initial values
    inits <- function(){list(beta1 = rnorm(d, 0, 1))}

    ## Parameters of interest
    parameters <- c('beta1', 'tau1', 'tau2')

    cat("
    model {
    #model
    for (i in 1:N){
    y[i, ] ~ ddirch(alpha[i,])

    log(alpha[i,1]) <- beta1[1]*V[i,1] + w1[niv1[i]]
    log(alpha[i,2]) <- beta1[2]*V[i,2] + w1[niv1[i]]
    log(alpha[i,3]) <- beta1[3]*V[i,3] + w2[niv2[i]]
    log(alpha[i,4]) <- beta1[4]*V[i,4] + w2[niv2[i]]
    }

    for(j in 1:niv_length1){
        w1[j] ~ dnorm(0, tau1)
    }

     for(j in 1:niv_length2){
        w2[j] ~ dnorm(0, tau2)
    }

    #priors
    tau1 <- x_pc[index1]
    index1 ~ dcat(y_pc)

    #priors
    tau2 <- x_pc[index2]
    index2 ~ dcat(y_pc)

    for (c in 1:d)
    {
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
  cat(paste0("n = ", n, " - levels_factor = ", levels_factor," -----> Fitting using INLA \n"))
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
    model.inla <- dirinlareg( formula  = y ~
                                 -1 + v1 + f(iid1, model = 'iid', hyper=list(theta=list(prior="loggamma",param=c(1,0.1)))) |
                                 -1 + v2 + f(iid1, model = 'iid', hyper=list(theta=list(prior="loggamma",param=c(1,0.1)))) |
                                 -1 + v3 + f(iid2, model = 'iid', hyper=list(theta=list(prior="loggamma",param=c(1,0.1)))) |
                                 -1 + v4 + f(iid2, model = 'iid', hyper=list(theta=list(prior="loggamma",param=c(1,0.1)))),
                y        = y,
                data.cov = V,
                prec     = 0.01,
                verbose  = FALSE,
                control.inla = list(strategy = "laplace", int.strategy = "grid"),
)

    t_inla <- proc.time()-t    # Stop the time
    t_inla <- t_inla[3]
    summary(model.inla)
  }

  ### ----- 3.2. Fitting the model with INLA half gaussian --- ####
  HN.prior <<- "expression:
  tau0 = 1;
  sigma = exp(-theta/2);
  log_dens = log(2) - 0.5 * log(2 * pi) + 0.5 * log(tau0);
  log_dens = log_dens - 0.5 * tau0 * sigma^2;
  log_dens = log_dens - log(2) - theta / 2;
  return(log_dens);
"



  half.normal  <<- list(prec = list(prior = HN.prior))
  cat(paste0("n = ", n, " - levels_factor = ", levels_factor," -----> Fitting using INLA \n"))
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
    # formula  = y ~ -1 + v1 + f(iid1, model = 'iid', hyper=list(theta=(list(prior="pc.prec", param=c(10, 0.01))))) |
    #   -1 + v2 + f(iid1, model = 'iid', hyper=list(theta=(list(prior="pc.prec", param=c(10, 0.01))))) |
    #   -1 + v3 + f(iid2, model = 'iid', hyper=list(theta=(list(prior="pc.prec", param=c(10, 0.01))))) |
    #   -1 + v4 + f(iid2, model = 'iid', hyper=list(theta=(list(prior="pc.prec", param=c(10, 0.01)))))
    model.inla.2 <- dirinlareg( formula  = y ~ -1 + v1 + f(iid1, model = 'iid', hyper = half.normal, constr=TRUE) |
                                  -1 + v2 + f(iid1, model = 'iid', hyper = half.normal, constr = TRUE) |
                                  -1 + v3 + f(iid2, model = 'iid', hyper = half.normal, constr = TRUE) |
                                  -1 + v4 + f(iid2, model = 'iid', hyper = half.normal, constr = TRUE),
                                y        = y,
                                data.cov = V,
                                prec     = 0.01,
                                verbose  = FALSE,
                                control.inla = list(strategy = "laplace", int.strategy = "grid"))

    t_inla_2 <- proc.time()-t    # Stop the time
    t_inla_2 <- t_inla_2[3]
    summary(model.inla.2)
  }


  ### ----- 3.3. Fitting the model with long jags --- ####
  cat(paste0("n = ", n, " - levels_factor = ", levels_factor," -----> Fitting using long JAGS \n"))

  if(file.exists(paste0("model_jags_long_", n,".RDS"))){
    model.jags.2 <- readRDS(paste0("model_jags_long_", n,".RDS"))
    if(n< 1000)
    {
      simulation <- readRDS("simulation4_50-500.RDS")
      t_jags_2 <- simulation[[paste0("n",n)]]$times[3]
    }else{
      simulation <- readRDS("simulation4_1000-10000.RDS")
      t_jags_2 <-simulation[[paste0("n",n)]]$times[3]
    }
  }else{
    ## MCMC configuration
     # ni <- 100000
     # nb <- 10000
    #ni <- 1000
    nt <- 5
    #nb <- 10
    nc <- 3


 # ni <- 1000
 # nb <- 100
    ## Data set
    table(V$iid1) %>%length() -> niv_length1
    table(V$iid2) %>%length() -> niv_length2
    data_jags <- list(y = y,
                      N = dim(y)[1],
                      d = d,
                      V = V,
                      niv1 = V$iid1,
                      niv2 = V$iid2,
                      niv_length1 = niv_length1,
                      niv_length2 = niv_length2
                      #x_pc = prec_pc$x_pc,
                      #y_pc = prec_pc$y_pc
                      )

    ## Initial values
    inits <- function(){list(beta1 = rnorm(d, 0, 1),
                             sigma1 = runif(1, 0,1),
                             sigma2 = runif(1, 0,1))}

    ## Parameters of interest
    parameters <- c('beta1', 'sigma1', 'sigma2')

    cat("
    model {
    #model
    for (i in 1:N){
    y[i, ] ~ ddirch(alpha[i,])

    log(alpha[i,1]) <- beta1[1]*V[i,1] + w1[niv1[i]]
    log(alpha[i,2]) <- beta1[2]*V[i,2] + w1[niv1[i]]
    log(alpha[i,3]) <- beta1[3]*V[i,3] + w2[niv2[i]]
    log(alpha[i,4]) <- beta1[4]*V[i,4] + w2[niv2[i]]
    }

    for(j in 1:niv_length1){
        w1[j] ~ dnorm(0, tau1)
    }

     for(j in 1:niv_length2){
        w2[j] ~ dnorm(0, tau2)
    }

    #priors
    # tau1 <- x_pc[index1]
    # index1 ~ dcat(y_pc)

    # #priors
    # tau2 <- x_pc[index2]
    # index2 ~ dcat(y_pc)

    tau1 <- 1/(sigma1*sigma1)
    tau2 <- 1/(sigma2*sigma2)

    sigma1 ~ dnorm(0,1)T(0,)
    sigma2 ~ dnorm(0,1)T(0,)


    for (c in 1:d)
    {
    beta1[c] ~ dnorm(0, 0.01)
    }
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
  saveRDS(file = paste0("model_jags_", n,"_", levels_factor, ".RDS"), model.jags)
  saveRDS(file = paste0("model_jags_long_", n,"_", levels_factor, ".RDS"), model.jags.2)
  saveRDS(file = paste0("model_inla_", n,"_", levels_factor, ".RDS"), model.inla)
  saveRDS(file = paste0("model_inla_2_", n,"_", levels_factor, ".RDS"), model.inla.2)



  ### --- 4. Comparing methodologies --- ####
  cat(paste0("n = ", n, " - levels_factor = ", levels_factor," -----> Comparing methodologies \n"))

  ### ----- 4.1. Computational times --- ####
  times <- c(t_jags, t_inla, t_jags_2, t_inla_2)

  ### ----- 4.2. (E(INLA) - E(JAGS2))/SD(JAGS2) and variance ratios --- ####
  ratio1_beta1 <- ratio2_beta1 <-  ratio1_beta1_pc <-  ratio2_beta1_pc <- numeric()
  # for (i in 1:4)
  # {
  #   mean_jags_2 <- mean(model.jags.2$BUGSoutput$sims.list$beta0[,i])
  #   sd_jags_2 <- sd(model.jags.2$BUGSoutput$sims.list$beta0[,i])
  #   mean_jags_1 <- mean(model.jags$BUGSoutput$sims.list$beta0[,i])
  #   sd_jags_1 <- sd(model.jags$BUGSoutput$sims.list$beta0[,i])
  #   mean_inla <- model.inla$summary_fixed[[i]]$mean[1]
  #
  #   ratio1_beta0 <- c(ratio1_beta0, c(mean_inla - mean_jags_2)/sd_jags_2)
  #   ratio2_beta0 <- c(ratio2_beta0, sd(inla.rmarginal(10000, model.inla$marginals_fixed[[i]]$intercept))^2/sd_jags_2^2)
  # }

  for (i in 1:4)
  {
    mean_jags_2 <- mean(model.jags.2$BUGSoutput$sims.list$beta1[,i])
    sd_jags_2 <- sd(model.jags.2$BUGSoutput$sims.list$beta1[,i])
    #mean_jags_1 <- mean(model.jags$BUGSoutput$sims.list$beta0[,i])
    #sd_jags_1 <- sd(model.jags$BUGSoutput$sims.list$beta0[,i])
    mean_inla <- model.inla$summary_fixed[[i]]$mean[1]
    mean_inla_2 <- model.inla.2$summary_fixed[[i]]$mean[1]

    ratio1_beta1 <- c(ratio1_beta1, c(mean_inla - mean_jags_2)/sd_jags_2)
    ratio2_beta1 <- c(ratio2_beta1, sd(inla.rmarginal(10000, model.inla$marginals_fixed[[i]][[1]]))^2/(sd_jags_2^2))
    ratio1_beta1_pc <- c(ratio1_beta1_pc, c(mean_inla_2 - mean_jags_2)/sd_jags_2)
    ratio2_beta1_pc <- c(ratio2_beta1_pc, sd(inla.rmarginal(10000, model.inla.2$marginals_fixed[[i]][[1]]))^2/(sd_jags_2^2))

  }
#
#   mean_jags_2_sigma <- c("tau1", "tau2") %>% lapply(., function(x) mean(1/sqrt(model.jags.2$BUGSoutput$sims.list[[c(x)]]))) %>% unlist(.)
#   sd_jags_2_sigma <- c("tau1", "tau2") %>% lapply(., function(x) sd(1/sqrt(model.jags.2$BUGSoutput$sims.list[[c(x)]]))) %>% unlist(.)
#
#   mean_jags_2_sigma_log <- c("tau1", "tau2") %>% lapply(., function(x) mean(log(1/sqrt(model.jags.2$BUGSoutput$sims.list[[c(x)]])))) %>% unlist(.)
#   sd_jags_2_sigma_log <- c("tau1", "tau2") %>% lapply(., function(x) sd(log(1/sqrt(model.jags.2$BUGSoutput$sims.list[[c(x)]])))) %>% unlist(.)
#
#
#
  mean_jags_2_sigma <- c("sigma1", "sigma2") %>% lapply(., function(x) mean(model.jags.2$BUGSoutput$sims.list[[c(x)]])) %>% unlist(.)
  sd_jags_2_sigma <- c("sigma1", "sigma2") %>% lapply(., function(x) sd(model.jags.2$BUGSoutput$sims.list[[c(x)]])) %>% unlist(.)

  mean_jags_2_sigma_log <- c("sigma1", "sigma2") %>% lapply(., function(x) mean(log(model.jags.2$BUGSoutput$sims.list[[c(x)]]))) %>% unlist(.)
  sd_jags_2_sigma_log <- c("sigma1", "sigma2") %>% lapply(., function(x) sd(log(model.jags.2$BUGSoutput$sims.list[[c(x)]]))) %>% unlist(.)

  #inla_sigma
  inla_sigma <- lapply(1:2, function(x) inla.tmarginal(function(x) 1/sqrt(x), inla.smarginal(model.inla$marginals_hyperpar[[x]], factor = 100)))
  inla_sigma_2 <- lapply(1:2, function(x) inla.tmarginal(function(x) 1/sqrt(x), inla.smarginal(model.inla.2$marginals_hyperpar[[x]], factor = 100)))
  names(inla_sigma) <- c("sigma1", "sigma2")
  names(inla_sigma_2) <- c("sigma1", "sigma2")

  #inla sigma log
  inla_sigma_log <- lapply(1:2, function(x) inla.tmarginal(function(x) log(1/sqrt(x)), inla.smarginal(model.inla$marginals_hyperpar[[x]], factor = 100)))
  inla_sigma_log_2 <- lapply(1:2, function(x) inla.tmarginal(function(x) log(1/sqrt(x)), inla.smarginal(model.inla.2$marginals_hyperpar[[x]], factor = 100)))
  names(inla_sigma_log) <- c("log(sigma1)", "log(sigma2)")
  names(inla_sigma_log_2) <- c("log(sigma1)", "log(sigma2)")



  # inla_sigma_sim <- lapply(1:2, function(x) 1/sqrt(inla.rmarginal(10000, model.inla$marginals_hyperpar[[x]])))
  # inla_sigma_sim_2 <- lapply(1:2, function(x) 1/sqrt(inla.rmarginal(10000, model.inla.2$marginals_hyperpar[[x]])))
  # names(inla_sigma_sim) <- c("sigma1", "sigma2")
  # names(inla_sigma_sim_2) <- c("sigma1", "sigma2")

  #Ratios sigma
  ratio1_sigma <- (lapply(inla_sigma, function(x) inla.zmarginal(x, silent = TRUE)[[1]]) %>% unlist() - mean_jags_2_sigma)/sd_jags_2_sigma
  ratio2_sigma <- (lapply(inla_sigma, function(x) inla.zmarginal(x, silent = TRUE)[[2]]) %>% unlist())^2/(sd_jags_2_sigma^2)


  ratio1_sigma_pc <- (lapply(inla_sigma_2, function(x) inla.zmarginal(x, silent = TRUE)[[1]]) %>% unlist() - mean_jags_2_sigma)/sd_jags_2_sigma
  ratio2_sigma_pc <- (lapply(inla_sigma_2, function(x) inla.zmarginal(x, silent = TRUE)[[2]]) %>% unlist())^2/(sd_jags_2_sigma^2)

  #Ratios logarithm
  ratio1_sigma_log <- (lapply(inla_sigma_log, function(x) inla.zmarginal(x, silent = TRUE)[[1]]) %>% unlist() - mean_jags_2_sigma_log)/sd_jags_2_sigma_log
  ratio2_sigma_log <- (lapply(inla_sigma_log, function(x) inla.zmarginal(x, silent = TRUE)[[2]]) %>% unlist())^2/(sd_jags_2_sigma_log^2)


  ratio1_sigma_log_pc <- (lapply(inla_sigma_log_2, function(x) inla.zmarginal(x, silent = TRUE)[[1]]) %>% unlist() - mean_jags_2_sigma_log)/sd_jags_2_sigma_log
  ratio2_sigma_log_pc <- (lapply(inla_sigma_log_2, function(x) inla.zmarginal(x, silent = TRUE)[[2]]) %>% unlist())^2/(sd_jags_2_sigma_log^2)






  ### ----- 4.3. Mean and sd of the posterior distributions --- ####
  ### Intercepts
  #result_beta0 <- numeric()
  # for(i in 1:4)
  # {
  #   result_beta0 <- rbind(result_beta0,
  #                         t(matrix(c(model.jags$BUGSoutput$summary[paste0("beta0[", i,"]"), c("mean", "sd")],
  #                                    model.inla$summary_fixed[[i]][1,c("mean", "sd")],
  #                                    model.jags.2$BUGSoutput$summary[paste0("beta0[", i,"]"), c("mean", "sd")]))))
  # }
  # rownames(result_beta0) <- paste0("beta0", 1:4)
  # colnames(result_beta0) <- c(paste0("JAGS", c("_mean", "_sigma")),
  #                             paste0("INLA", c("_mean", "_sigma")),
  #                             paste0("LONG_JAGS", c("_mean", "_sigma")))


  ### Beta1
  result_beta1 <- numeric()
  for(i in 1:4)
  {
    result_beta1 <- rbind(result_beta1,
                          t(matrix(c(model.jags$BUGSoutput$summary[paste0("beta1[", i,"]"), c("mean", "sd")],
                                     model.inla$summary_fixed[[i]][1,c("mean", "sd")],
                                     model.jags.2$BUGSoutput$summary[paste0("beta1[", i,"]"), c("mean", "sd")],
                                     model.inla.2$summary_fixed[[i]][1,c("mean", "sd")]))))
  }
  rownames(result_beta1) <- paste0("beta1", 1:4)
  colnames(result_beta1) <- c(paste0("JAGS", c("_mean", "_sigma")),
                              paste0("INLA", c("_mean", "_sigma")),
                              paste0("LONG_JAGS", c("_mean", "_sigma")),
                              paste0("INLA_PC", c("_mean", "_sigma")))
  ### --- 5. Plotting --- ####
  ### ----- 5.1. intercepts --- ####
  # ## Intercept
  # p1 <- list()
  # beta0 <- expression(paste("p(", beta[0], "|", "y)"))
  #
  # for (i in 1:length(model.inla$marginals_fixed))
  # {
  #   #jags1
  #   dens <- density(model.jags$BUGSoutput$sims.matrix[,i], adjust = 2)
  #   dens <- as.data.frame(cbind(dens$x, dens$y))
  #   colnames(dens) <- c("x", "y")
  #
  #   #jags2
  #   dens2 <- density(model.jags.2$BUGSoutput$sims.matrix[,i], adjust = 2)
  #   dens2 <- as.data.frame(cbind(dens2$x, dens2$y))
  #   colnames(dens2) <- c("x", "y")
  #
  #   #Data combining jags (1) and inla (2)
  #   dens <- rbind(cbind(dens, group = 1),
  #                 cbind(as.data.frame(model.inla$marginals_fixed[[i]][[1]]), group = 2),
  #                 cbind(dens2, group = 3))
  #   dens$group <- factor(dens$group,
  #                        labels = c("R-JAGS", "dirinla", "long R-JAGS"))
  #
  #   ### Intercept
  #   p1[[i]] <- ggplot(dens,
  #                     aes(x = x,
  #                         y = y,
  #                         group = group
  #                         #colour = factor(group)
  #                     ))  +
  #     geom_line(size = 0.6, aes(linetype = group ,
  #                               color    = group)) +
  #     xlim(c(min(dens$x[dens$group=="R-JAGS"]), max(dens$x[dens$group=="R-JAGS"]))) +
  #     theme_bw() + #Show axes
  #     xlab(expression(beta[0])) + #xlab
  #     ylab(beta0) #ylab
  #
  #
  #   #Frequentist approach
  #   p1[[i]] <- p1[[i]] + geom_vline(xintercept = x[seq(1,8, by=2)][i])
  #
  #
  #
  #   ### --- legend --- ###
  #   p1[[i]]<- p1[[i]] + theme(legend.position   = c(0.2, 0.8),
  #                             legend.title      = element_blank(),
  #                             legend.background = element_rect(colour = "gray"),
  #                             legend.key        = element_rect(colour = "white", fill="white"),
  #                             legend.key.size   = unit(0.5, "cm")) +
  #     theme(legend.text = element_text(size = 9)) +
  #     # scale_fill_manual(labels=c("R-JAGS", "dirinla", "long R-JAGS"),
  #     #                   values = c("darkgreen", "red4", "blue4" )) +
  #     scale_colour_manual (
  #       values= c("darkgreen", "red4", "blue4")) +
  #     scale_linetype_manual(labels=c("R-JAGS", "dirinla", "long R-JAGS"),
  #                           values=c("dotted", "twodash",  "solid"))
  #
  #   if(i!=1)
  #   {
  #     p1[[i]] <- p1[[i]] + theme(legend.position="none")
  #   }
  #
  #   p1[[i]] <- p1[[i]] + ggtitle(paste0("Category ", i)) +
  #     theme(
  #       plot.title = element_text(color = "black",
  #                                 size  = 15,
  #                                 face  = "bold.italic",
  #                                 hjust = 0.5))
  # }

  #
  # pdf("example_simulation4_intercepts_50.pdf", width = 18, height = 4)
  #   gridExtra::grid.arrange(p1[[1]], p1[[2]], p1[[3]], p1[[4]], ncol = 4)
  # dev.off()



  ### ----- 5.2. slopes --- ####
  p2 <- list()
  beta1 <- expression(paste("p(", beta[1], "|", "y)"))

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
                  cbind(dens2, group = 3),
                  cbind(as.data.frame(model.inla.2$marginals_fixed[[i]][[1]]), group = 4)
    )
    dens$group <- factor(dens$group,
                         labels = c("R-JAGS", "dirinla", "long R-JAGS", "dirinla pc"))

    ###
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
    p2[[i]] <- p2[[i]] + geom_vline(xintercept = x[seq(1,4, by=1)][i])



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
        values= c("darkgreen", "red4", "blue4", "orange2")) +
      scale_linetype_manual(labels=c("R-JAGS", "dirinla", "long R-JAGS", "dirinla pc"),
                            values=c("dotted", "twodash",  "solid", "longdash"))

    if(i!=1)
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



  #### Hyperparemeters

  #jags1
  #logSigma1 and logSigma2
  lapply(c("tau1", "tau2"), function(x){
    model.jags$BUGSoutput$sims.matrix[,c(x)] %>%
      .^(-1/2) %>%
      log(.) %>%
      density(., adjust = 2) %>%
      .[1:2] %>%
      do.call(data.frame, .) -> dens
    colnames(dens) <- c("x", "y")
    dens
  }) -> dens_log
  names(dens_log) <- c("log(sigma1)", "log(sigma2)")


  #Sigma1 and Sigma2
  lapply(c("tau1", "tau2"), function(x){
    model.jags$BUGSoutput$sims.matrix[,c(x)] %>%
      .^(-1/2) %>%
      log(.) %>%
      density(., adjust = 2) %>%
      .[1:2] %>%
      inla.tmarginal(fun = function(x)exp(x), .) %>%
      do.call(data.frame, .) -> dens
    colnames(dens) <- c("x", "y")
    dens
  }) -> dens
  names(dens) <- c("sigma1", "sigma2")
  #
  # dens <- density(model.jags$BUGSoutput$sims.matrix[,c("sigma1")], adjust = 2)
  # dens <- as.data.frame(cbind(dens$x, dens$y))



  #jags2
  #logSigma1 and logSigma2
  lapply(c("sigma1", "sigma2"), function(x){
    model.jags.2$BUGSoutput$sims.matrix[,c(x)] %>%
      #.^(-1/2) %>%
      log(.) %>%
      density(., adjust = 2) %>%
      .[1:2] %>%
      do.call(data.frame, .) -> dens
    colnames(dens) <- c("x", "y")
    dens
  }) -> dens_log2
  names(dens_log2) <- c("log(sigma1)", "log(sigma2)")



  lapply(c("sigma1", "sigma2"), function(x){
    model.jags.2$BUGSoutput$sims.matrix[,c(x)] %>%
      #.^(-1/2) %>%
      log(.)%>%
      density(., adjust = 2) %>%
      .[1:2] %>%
      inla.tmarginal(fun = function(x)exp(x), .) %>%
      do.call(data.frame, .) -> dens
    colnames(dens) <- c("x", "y")
    dens
  }) -> dens2
  names(dens2) <- c("sigma1", "sigma2")

# plot(dens2$sigma1, xlim = c(0.1, 2))
# plot(dens2$sigma2, xlim = c(0.1, 2))
# (dens2$sigma1$x < 0.05) %>% table(.)

# lapply(dens2, function(a){
#   a %>% dplyr::filter(x > 0.07)
# }) -> dens2

  dens <- c(dens, dens_log)
  dens2 <- c(dens2, dens_log2)
  inla_sigma <- c(inla_sigma, inla_sigma_log)
  inla_sigma_2 <- c(inla_sigma_2, inla_sigma_log_2)


  ### ----- 5.3. hyperpar --- ####
  sd_name <- list(expression(paste("p(", sigma_1, "|", "y)")),
                  expression(paste("p(", sigma_2, "|", "y)")),
                  expression(paste("p(", "log(sigma_1)", "|", "y)")),
                  expression(paste("p(", "log(sigma_2)", "|", "y)")))
  sd_name2 <- list(expression(sigma_1), expression(sigma_2),
                   expression("log(sigma_1)"),
                   expression("log(sigma_2)"))

  p3 <- list()
  for(i in 1:length(dens))
  {
    #Data combining jags (1) and inla (2)
    dens_sigma <- rbind(cbind(dens[[i]], group = 1),
                         cbind(as.data.frame(inla_sigma[[i]]), group = 2),
                         cbind(dens2[[i]], group = 3),
                         cbind(as.data.frame(inla_sigma_2[[i]]), group = 4))
    dens_sigma$group <- factor(dens_sigma$group,
                                labels = c("R-JAGS", "dirinla", "long R-JAGS", "dirinla pc"))


    p3[[i]] <- ggplot(dens_sigma,
                 aes(x = x,
                     y = y,
                     group = group
                     #colour = factor(group)
                 ))  +
      geom_line(size = 0.6, aes(linetype = group ,
                                color    = group)) +
      xlim(c(min(dens_sigma$x[dens_sigma$group=="R-JAGS"]), max(dens_sigma$x[dens_sigma$group=="R-JAGS"]))) +
      theme_bw() + #Show axes
      xlab(sd_name2[[i]]) + #xlab
      ylab(sd_name[[i]]) #ylab


    #Frequentist approach
    constant <- c(1/sqrt(prec_w), log(1/sqrt(prec_w)))
    p3[[i]] <- p3[[i]] + geom_vline(xintercept = constant[i])


    #xmax <- ifelse(n <= 100, 100, 50)
    xmax <- 2
    ### --- legend --- ###
    p3[[i]] <- p3[[i]] + theme(legend.position   = c(0.8, 0.8),
                     legend.title      = element_blank(),
                     legend.background = element_rect(colour = "gray"),
                     legend.key        = element_rect(colour = "white", fill="white"),
                     legend.key.size   = unit(0.5, "cm")) +
      theme(legend.text = element_text(size = 9)) +
      # scale_fill_manual(labels=c("R-JAGS", "dirinla", "long R-JAGS"),
      #                   values = c("darkgreen", "red4", "blue4" )) +
      scale_colour_manual (
        values= c("darkgreen", "red4", "blue4", "orange2")) +
      scale_linetype_manual(labels=c("R-JAGS", "dirinla", "long R-JAGS", "dirinla pc"),
                            values=c("dotted", "twodash",  "solid", "longdash"))


    if(i!=1)
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




  # pdf("example_simulation4_slopes_50.pdf", width = 18, height = 4)
  #   gridExtra::grid.arrange(p2[[1]], p2[[2]], p2[[3]], p2[[4]], ncol = 4)
  # dev.off()



  # pdf(paste0("examples_simulation4_slopes_", n ,"_", levels_factor, ".pdf"), width = 15, height = 4)
  # gridExtra::grid.arrange(
  #                         p2[[1]], p2[[2]], p2[[3]], p2[[4]], ncol = 4)
  # dev.off()
  #
  # pdf(paste0("examples_simulation4_sigma_", n ,"_", levels_factor,".pdf"), width = 15, height = 4)
  # gridExtra::grid.arrange(p3[[1]], p3[[2]], p3[[3]], p3[[4]], ncol = 4)
  # dev.off()

  png(paste0("examples_simulation4_slopes_", n ,"_", levels_factor, ".png"), width = 2000, height = 600, res = 150)
  gridExtra::grid.arrange(
    p2[[1]], p2[[2]], p2[[3]], p2[[4]], ncol = 4)
  dev.off()

  png(paste0("examples_simulation4_sigma_", n ,"_", levels_factor,".png"), width = 2000, height = 600, res = 120)
  gridExtra::grid.arrange(p3[[1]], p3[[2]], p3[[3]], p3[[4]], ncol = 4)
  dev.off()

  check_jags1 <- print(model.jags)
  check_jags2 <- print(model.jags.2)

  res_check_jags1 <- check_jags1$summary %>% data.frame(.) %>% select(Rhat, n.eff)
  res_check_jags2 <- check_jags2$summary %>% data.frame(.) %>% select(Rhat, n.eff)


  list(times = times,
       #intercepts = result_beta0,
       slopes     = result_beta1,
       #ratio1_intercepts = ratio1_beta0,
       ratio1_slopes = ratio1_beta1,
       #ratio2_intercepts = ratio2_beta0,
       ratio2_slopes = ratio2_beta1,
       ratio1_sigma = ratio1_sigma,
       ratio2_sigma = ratio2_sigma,
       ratio1_sigma_pc = ratio1_sigma_pc,
       ratio2_sigma_pc = ratio2_sigma_pc,
       ratio1_sigma_log = ratio1_sigma_log,
       ratio2_sigma_log = ratio2_sigma_log,
       ratio1_sigma_log_pc = ratio1_sigma_log_pc,
       ratio2_sigma_log_pc = ratio2_sigma_log_pc,
       res_check_jags1 = res_check_jags1,
       res_check_jags2 = res_check_jags2,
       n_levels = paste0(n, "-", levels_factor))
}


 ### --- 3. Calling the function --- ####
n <- c(50, 100)
levels_factor <- c(25, NA)

levels_factor <- c(5, 10, 25, NA)

levels_factor <- c(5, 10, 25, NA)
n <- c(50, 100, 500)
#levels_factor <- c(5, 10, 25, NA)
# levels_factor <- c(25, NA)

arguments <- expand.grid(n, levels_factor)
n <- arguments[,1]
levels_factor <- arguments[,2]

a <- mapply(simulations_with_slopes_iid,
       n = n,
       levels_factor = levels_factor)

a <- simulations_with_slopes_iid(50, 50)


# a[,1]
# a <- parallel::mcmapply(simulations_with_slopes_iid,
#                         n = n,
#                         levels_factor = levels_factor,
#                         mc.cores = 3)
#
# a <- parallel::mclapply(n, simulations_with_slopes_iid, level_factor = 2,
#                         mc.cores = 3)
#n <- c(50, 100)
#a <- lapply(n, simulations_with_slopes_iid, levels_factor = 10)
names(a) <- paste0("n", n)
saveRDS(a, file = "simulation4_50-500.RDS")
a <- readRDS(file = "simulation4_50-500.RDS")

levels_factor[c(10,11,12)] <- c(50, 100, 500)
colnames(a) <- paste0(n, "-", levels_factor)
a[, c("50-5")]

a[, c("50-5")]
a$n50$times
a$n50$ratio1_slopes
a

paste0()

n <- c(1000, 10000)
a <- parallel::mclapply(n, simulations_with_slopes_iid,
                        mc.cores = 2)
a <- lapply(n, simulations_with_slopes_iid)
names(a) <- paste0("n", n)

saveRDS(a, file = "simulation4_1000-10000.RDS")
b <- readRDS(file = "simulation4_1000-10000.RDS")
results <- c(a,b)

### --- 4. Extracting tables for the paper --- ####
# results <- readRDS(file = "simulation4_50-500.RDS")
results <- a
results$n50$times
results$n50$intercepts
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

result_ratio2

#Latex
library(xtable)
xtable(result_time, digits = 4)
xtable(result_ratio1, digits = 4)
xtable(result_ratio2, digits = 4)









