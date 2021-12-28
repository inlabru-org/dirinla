### In this script simulations with N=100 and d = 5, 10, 15, 20, 30           #
#  are conducted in order to check                                            #
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

### --- 2. Function for simulation --- ####
simulations_just_intercepts_nd <- function(n, d)
{
  ### --- 2. Simulation data --- ####
  set.seed(1000)

  cat("n = ", n,  ", d = ", d, " -----> Simulating data \n")
  #Covariates
  V <- as.data.frame(matrix(runif((10)*n, 0, 1), ncol=10))
  names(V) <- paste0('v', 1:(10))

  # Formula that we want to fit
  b1 <- paste0(paste0((rep("1 | ", d - 1)) , collapse = ""), 1, collapse = " ")
  formula <- eval(parse(text = paste0("y ~ ", b1)))

  names_cat <- formula_list(formula, y= NULL)

  # Parameters to fit
  if(d <4){
    print("d should be bigger than 4")
  }else if(d == 4){
    x <- c(-2.4, 1.2, -3.1, 1.3)
  }else{
    x <- c(-2.4, 1.2, -3.1, 1.3, runif(d - 4, -2, 2))
  }

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
  cat(paste0("n = ", n, ", d = ", d, " -----> Fitting using SHORT JAGS \n"))
  if(file.exists(paste0("model_jags_", n, "_", d, ".RDS"))){
    model.jags <- readRDS(paste0("model_jags_", n, "_", d, ".RDS"))
    if(d < 50)
    {
      simulation <- readRDS("simulation3_n_d.RDS")
      t_jags <- simulation[[paste0("n",n)]][[paste0("d",d)]]$times[1]
    }else{
      simulation <- readRDS("simulation3_n_d2.RDS")
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
    parameters <- c('beta0')

    cat("
    model {
    #model
    for (i in 1:N){
    y[i, ] ~ ddirch(alpha[i,])
      for(j in 1:d)
        {
          log(alpha[i,j]) <- beta0[j]
        }
    }

    #priors
    for (j in 1:d)
    {
      beta0[j]  ~ dnorm(0, 0.0001)
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

  saveRDS(file = paste0("model_jags_", n,"_", d, ".RDS"), model.jags)



  ### ----- 3.2. Fitting the model with INLA --- ####
  cat(paste0("n = ", n,  ", d = ", d, " -----> Fitting using INLA \n"))
  if(file.exists(paste0("model_inla_", n, "_", d, ".RDS"))){
    model.inla <- readRDS(paste0("model_inla_", n, "_", d, ".RDS"))
    if(d < 50)
    {
      simulation <- readRDS("simulation3_n_d.RDS")
      t_inla <- simulation[[paste0("n",n)]][[paste0("d",d)]]$times[2]
    }else{
      simulation <- readRDS("simulation3_n_d2.RDS")
      t_inla <- simulation[[paste0("n",n)]][[paste0("d",d)]]$times[2]
    }
  }else{
    t <- proc.time() # Measure the time
    model.inla <- dirinlareg( formula  = formula ,
                              y        = y,
                              data.cov = V,
                              prec     = 0.0001,
                              verbose  = FALSE,
                              sim       = 10) #As we are not interested in the linear predictor just a few

    t_inla <- proc.time()-t    # Stop the time
    t_inla <- t_inla[3]
    summary(model.inla)
  }
  saveRDS(file = paste0("model_inla_", n, "_", d,".RDS"), model.inla)

  ### ----- 3.3. Fitting the model with long jags --- ####
  cat(paste0("n = ", n, ", d = ", d,  " -----> Fitting using long JAGS \n"))
  if(file.exists(paste0("model_jags_long_", n, "_", d, ".RDS"))){
    model.jags.2 <- readRDS(paste0("model_jags_long_", n, "_", d, ".RDS"))
    if(n< 1000)
    {
      simulation <- readRDS("simulation3_n_d.RDS")
      t_jags_2 <- simulation[[paste0("n",n)]][[paste0("d",d)]]$times[3]
    }else{
      simulation <- readRDS("simulation3_n_d.RDS")
      t_jags_2 <- simulation[[paste0("n",n)]][[paste0("d",d)]]$times[3]
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
    parameters <- c('beta0')

    cat("
    model {
    #model
    for (i in 1:N){
    y[i, ] ~ ddirch(alpha[i,])
      for(j in 1:d)
        {
          log(alpha[i,j]) <- beta0[j]
        }
    }

    #priors
    for (j in 1:d)
    {
      beta0[j]  ~ dnorm(0, 0.0001)
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
  saveRDS(file = paste0("model_jags_long_", n, "_", d, ".RDS"), model.jags.2)

  t_inla
  t_jags
  t_jags_2


  ### --- 4. Comparing methodologies --- ####
  cat(paste0("n = ", n,  ", d = ", d, " -----> Comparing methodologies \n"))

  ### ----- 4.1. Computational times --- ####
  times <- c(t_jags, t_inla, t_jags_2)

  ### ----- 4.2. (E(INLA) - E(JAGS2))/SD(JAGS2) and variance ratios --- ####
  ratio1_beta0 <- ratio2_beta0 <- ratio1_mu <- ratio2_mu <-  numeric()
  for (i in 1:d)
  {
    mean_jags_2 <- mean(model.jags.2$BUGSoutput$sims.list$beta0[,i])
    sd_jags_2 <- sd(model.jags.2$BUGSoutput$sims.list$beta0[,i])
    # mean_jags_1 <- mean(model.jags$BUGSoutput$sims.list$beta0[,1])
    # sd_jags_1 <- sd(model.jags$BUGSoutput$sims.list$beta0[,1])
    mean_inla <- model.inla$summary_fixed[[i]]$mean

    ratio1_beta0 <- c(ratio1_beta0, c(mean_inla - mean_jags_2)/sd_jags_2)
    ratio2_beta0 <- c(ratio2_beta0, sd(inla.rmarginal(1000, model.inla$marginals_fixed[[i]]$intercept))^2/sd_jags_2^2)
    print(paste0(d, i))
  }



  ### ----- 4.3. Mean and sd of the posterior distributions --- ####
  ### Intercepts
  result <- numeric()
  for(i in 1:d)
  {
    result <- rbind(result,
                    t(matrix(c(model.jags$BUGSoutput$summary[paste0("beta0[", i,"]"), c("mean", "sd")],
                           model.inla$summary_fixed[[i]][,c("mean", "sd")],
                           model.jags.2$BUGSoutput$summary[paste0("beta0[", i,"]"), c("mean", "sd")]))))
    print(i)
  }
  rownames(result) <- paste0("beta0", 1:d)
  colnames(result) <- c(paste0("JAGS", c("_mean", "_sd")),
                        paste0("INLA", c("_mean", "_sd")),
                        paste0("LONG_JAGS", c("_mean", "_sd")))



  list(times = times,
       intercepts       = result,
       ratio1_intercept = ratio1_beta0,
       ratio2_intercept = ratio2_beta0)
}


### --- 3. Calling the function --- ####
result <- list()

d <- c(5, 10, 15, 20, 30)
n <- c(100)
for(i in length(n))
{
  print(paste0("n = ", n))
  result[[i]] <- parallel::mclapply(d, simulations_just_intercepts_nd,
                          mc.cores = 3, n = n[i])
  names(result[[i]]) <- paste0("d", d)
}
names(result) <- paste0("n", n)
saveRDS(result, file = "simulation3_n_d.RDS")
result <- readRDS(file = "simulation3_n_d.RDS")





### --- 4. Computing ratio1 and ratio2 for R-JAGS ####
ratios_jags <- function(n, levels_factor)
{
  print(paste0(n, "-", levels_factor))
  if(is.na(levels_factor)){
    levels_factor <- n
  }
  model.jags <- readRDS(paste0("model_jags_", n, "_", levels_factor, ".RDS"))
  model.jags.2 <- readRDS(paste0("model_jags_long_", n, "_", levels_factor, ".RDS"))

  #Beta1
  ratio1_beta1_hn_jags <-  ratio2_beta1_hn_jags <- numeric()

  for (i in 1:4)
  {
    mean_jags_2 <- mean(model.jags.2$BUGSoutput$sims.list$beta1[,i])
    sd_jags_2 <- sd(model.jags.2$BUGSoutput$sims.list$beta1[,i])
    mean_jags_1 <- mean(model.jags$BUGSoutput$sims.list$beta1[,i])
    sd_jags_1 <- sd(model.jags$BUGSoutput$sims.list$beta1[,i])


    ratio1_beta1_hn_jags <- c(ratio1_beta1_hn_jags, (mean_jags_1 - mean_jags_2)/sd_jags_2)
    ratio2_beta1_hn_jags <- c(ratio2_beta1_hn_jags, (sd_jags_1^2)/(sd_jags_2^2))

  }

  #Sigmas
  mean_jags_2_sigma <- c("sigma1", "sigma2") %>% lapply(., function(x) mean(model.jags.2$BUGSoutput$sims.list[[c(x)]])) %>% unlist(.)
  sd_jags_2_sigma <- c("sigma1", "sigma2") %>% lapply(., function(x) sd(model.jags.2$BUGSoutput$sims.list[[c(x)]])) %>% unlist(.)

  mean_jags_2_sigma_log <- c("sigma1", "sigma2") %>% lapply(., function(x) mean(log(model.jags.2$BUGSoutput$sims.list[[c(x)]]))) %>% unlist(.)
  sd_jags_2_sigma_log <- c("sigma1", "sigma2") %>% lapply(., function(x) sd(log(model.jags.2$BUGSoutput$sims.list[[c(x)]]))) %>% unlist(.)

  mean_jags_1_sigma <- c("sigma1", "sigma2") %>% lapply(., function(x) mean(model.jags$BUGSoutput$sims.list[[c(x)]])) %>% unlist(.)
  sd_jags_1_sigma <- c("sigma1", "sigma2") %>% lapply(., function(x) sd(model.jags$BUGSoutput$sims.list[[c(x)]])) %>% unlist(.)

  mean_jags_1_sigma_log <- c("sigma1", "sigma2") %>% lapply(., function(x) mean(log(model.jags$BUGSoutput$sims.list[[c(x)]]))) %>% unlist(.)
  sd_jags_1_sigma_log <- c("sigma1", "sigma2") %>% lapply(., function(x) sd(log(model.jags$BUGSoutput$sims.list[[c(x)]]))) %>% unlist(.)


  #Ratios sigma
  ratio1_sigma_hn_jags <- (mean_jags_1_sigma - mean_jags_2_sigma)/sd_jags_2_sigma
  ratio2_sigma_hn_jags <- (sd_jags_1_sigma^2)/(sd_jags_2_sigma^2)

  #Ratios logarithm
  ratio1_sigma_log_hn_jags <- (mean_jags_1_sigma_log - mean_jags_2_sigma_log)/sd_jags_2_sigma_log
  ratio2_sigma_log_hn_jags <- (sd_jags_1_sigma_log^2)/(sd_jags_2_sigma_log^2)



  #Returning
  list(ratio1_beta1_hn_jags = ratio1_beta1_hn_jags,
       ratio2_beta1_hn_jags = ratio2_beta1_hn_jags,
       ratio1_sigma_hn_jags = ratio1_sigma_hn_jags,
       ratio2_sigma_hn_jags = ratio2_sigma_hn_jags,
       ratio1_sigma_log_hn_jags = ratio1_sigma_log_hn_jags,
       ratio2_sigma_log_hn_jags = ratio2_sigma_log_hn_jags)
}
### --- 5. Extracting tables for the paper --- ####
result <- readRDS("simulation3_n_d.RDS")

#Computational times
result_time <- rbind(result$n100$d5$times,
                     result$n100$d10$times,
                     result$n100$d15$times,
                     result$n100$d20$times,
                     result$n100$d30$times)
colnames(result_time) <- c("R-JAGS", "dirinla", "long R-JAGS")
rownames(result_time) <- paste0( c(5, 10, 15, 20, 30))
result_time

#Time
xtable(result_time, digits = 4)

result_ratio <- data.frame(ratio1 = round(result$n100$d5$ratio1_intercept, 4),
                                  ratio2 = round(sqrt(result$n100$d5$ratio2_intercept), 4),
                                  label = "C = 5")
colnames(result_ratio) <- c("ratio1", "ratio2", "dimension")
for(i in 2:length(d))
{
    result_ratio <- rbind(result_ratio,
          data.frame(ratio1 = round(result$n100[[i]]$ratio1_intercept, 4),
                     ratio2 = round(sqrt(result$n100[[i]]$ratio2_intercept), 4),
                     dimension = paste0("C = ", d[i])))
}

result_ratio[,3] <- ordered(as.factor(result_ratio[,3]),
                            c("C = 5", "C = 10", "C = 15", "C = 20", "C = 30"))

#Geom boxplot
pdf("boxplot_ratio1.pdf", width = 7, height = 4)
ggplot(result_ratio, aes(x= dimension, y = ratio1)) +
  geom_boxplot() +
  #ylim(c(0, 0.2)) +
  xlab("")+
  ylab(expression('ratio'[1])) +
  theme_bw()
dev.off()

pdf("boxplot_ratio2.pdf", width = 7, height = 4)
ggplot(result_ratio, aes(x= dimension, y = ratio2)) +
  geom_boxplot() +
  #ylim(c(0.90, 1.1)) +
  xlab("")+
  ylab(expression('ratio'[2])) +
  theme_bw()
dev.off()
