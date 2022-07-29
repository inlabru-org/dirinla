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
library(dplyr)

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
    parameters <- c('beta0', 'alpha[1,1]', 'alpha[1,2]', 'alpha[1,3]', 'alpha[1,4]', 'mu')

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
    alpha0 <- sum(alpha[1,])
    for(i in 1:d){
          mu[i] <- alpha[1,i]/alpha0
        }


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
    model.inla <- dirinlareg( formula  = formula  ,
                              y        = y,
                              data.cov = V,
                              prec     = 0.0001,
                              verbose  = TRUE,
                              sim       = 4000)

    t_inla <- proc.time()-t    # Stop the time
    t_inla <- t_inla[3]
    summary(model.inla)

    saveRDS(file = paste0("model_inla_", n, "_", d,".RDS"), model.inla)

  }
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
    parameters <- c('beta0', 'alpha[1,1]', 'alpha[1,2]', 'alpha[1,3]', 'alpha[1,4]', 'mu')

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
    alpha0 <- sum(alpha[1,])
    for(i in 1:d){
          mu[i] <- alpha[1,i]/alpha0
        }


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
  saveRDS(file = paste0("model_jags_long_", n, "_", d, ".RDS"), model.jags.2)


  t_inla
  t_jags
  t_jags_2


  ### --- 4. Comparing methodologies --- ####
  cat(paste0("n = ", n, " -----> Comparing methodologies \n"))

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
    ratio2_beta0 <- c(ratio2_beta0, sd(inla.rmarginal(10000, model.inla$marginals_fixed[[i]]$intercept))^2/sd_jags_2^2)
  }

  for (i in 1:d)
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

  list(times = times,
       intercepts       = result,
       mus              = result_mus,
       ratio1_intercept = ratio1_beta0,
       ratio2_intercept = ratio2_beta0,
       ratio1_mu        = ratio1_mu,
       ratio2_mu        = ratio2_mu)
}


### --- 3. Calling the function --- ####
result <- list()

d <- c(5, 10, 15, 20, 30)

n <- c(100)
for(i in length(n))
{
  print(paste0("n = ", n))
  result[[i]] <- lapply(d, simulations_just_intercepts_nd,
                           n = n[i])
  names(result[[i]]) <- paste0("d", d)
}
names(result) <- paste0("n", n)
saveRDS(result, file = "simulation3_n_d.RDS")
result <- readRDS(file = "simulation3_n_d.RDS")



### --- 4. Computing ratio1 and ratio2 for R-JAGS ####
ratios_jags <- function(n = 100, d)
{
  print(paste0(n, "-", d))

  model.jags <- readRDS(paste0("model_jags_", n, "_", d, ".RDS"))
  model.jags.2 <- readRDS(paste0("model_jags_long_", n, "_", d, ".RDS"))

  #Beta1
  ratio1_beta0_jags <-  ratio2_beta0_jags <- numeric()

  for (i in 1:d)
  {
    mean_jags_2 <- mean(model.jags.2$BUGSoutput$sims.list$mu[,i])
    sd_jags_2 <- sd(model.jags.2$BUGSoutput$sims.list$mu[,i])
    mean_jags_1 <- mean(model.jags$BUGSoutput$sims.list$mu[,i])
    sd_jags_1 <- sd(model.jags$BUGSoutput$sims.list$mu[,i])

    ratio1_beta0_jags <- c(ratio1_beta0_jags, (mean_jags_1 - mean_jags_2)/sd_jags_2)
    ratio2_beta0_jags <- c(ratio2_beta0_jags, (sd_jags_1^2)/(sd_jags_2^2))
  }

  #Returning
  list(ratio1_beta0_jags = ratio1_beta0_jags,
       ratio2_beta0_jags = ratio2_beta0_jags)
}

ratios_jags(n = 100, d = 5)


res_ratios <- d %>% lapply(., ratios_jags, n = 100)
names(res_ratios) <- paste0("d", d)
res_ratios

saveRDS(res_ratios, file = "simulation3_ratios_jags.RDS")


### --- 5. Extracting tables for the paper --- ####
result <- readRDS("simulation3_n_d.RDS")
res_ratio <- readRDS("simulation3_ratios_jags.RDS")
d <- c(5, 10, 15, 20, 30)
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


### Ratios dirinla
result_ratio <- data.frame(ratio1 = round(result$n100$d5$ratio1_mu, 4),
                                  ratio2 = round(sqrt(result$n100$d5$ratio2_mu), 4),
                                  label = "C = 5",
                                  clas  = "dirinla")
colnames(result_ratio) <- c("ratio1", "ratio2", "dimension", "clas")
for(i in 2:length(d))
{
    result_ratio <- rbind(result_ratio,
          data.frame(ratio1 = round(result$n100[[i]]$ratio1_mu, 4),
                     ratio2 = round(sqrt(result$n100[[i]]$ratio2_mu), 4),
                     dimension = paste0("C = ", d[i]),
                     clas      = "dirinla"))
}



### Ratios R-JAGS
for(i in 1:length(d))
{
  result_ratio <- rbind(result_ratio,
                        data.frame(ratio1 = round(res_ratio[[i]]$ratio1_beta0_jags, 4),
                                   ratio2 = round(sqrt(res_ratio[[i]]$ratio2_beta0_jags), 4),
                                   dimension = paste0("C = ", d[i]),
                                   clas      = "R-JAGS"))
}

# result_ratio[,3] <- ordered(as.factor(result_ratio[,3]),
#                             c("C = 5", "C = 10", "C = 15", "C = 20", "C = 30"))

result_ratio[,3] <- ordered(as.factor(result_ratio[,3]),
                            c("C = 5", "C = 10", "C = 15", "C = 20", "C = 30"))

result_ratio$clas <- as.factor(result_ratio$clas)

# result_ratio_total1 <- result_ratio %>% dplyr::select(ratio1, dimension, clas)
# result_ratio_total2 <- result_ratio %>% dplyr::select(ratio2, dimension, clas)
# result_ratio_total1 <- cbind(result_ratio_total1, ratio = "ratio1")


#Geom boxplot
a <- ggplot(result_ratio, aes(x= dimension, y = ratio1, fill = clas)) +
  geom_boxplot() +
  scale_fill_manual(values = c("dirinla" = "white",
                               "R-JAGS" = "gray")) +
  #ylim(c(0, 0.2)) +
  xlab("")+
  ylab(expression('ratio'[1])) +
  theme_bw() +
  theme(legend.title=element_blank(),
        legend.position = "top")

b <- ggplot(result_ratio, aes(x= dimension, y = ratio2, fill = clas)) +
  geom_boxplot() +
  scale_fill_manual(values = c("dirinla" = "white",
                               "R-JAGS" = "gray")) +
  #ylim(c(0.90, 1.1)) +
  xlab("")+
  ylab(expression('ratio'[2])) +
  theme_bw()


library(ggpubr)
pdf("boxplot_ratios.pdf", width = 8, height = 5)
ggpubr::ggarrange(a, b, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()
### Boxplot

### --- 6. Plotting results --- ####


### ----- 6.1. Times --- ####
result_time2 <- cbind(result_time, C = as.numeric(rownames(result_time))) %>% as.data.frame(.)
result_time2 %>%
  tidyr::pivot_longer(1:3, names_to = "method") -> result_time2
result_time2$method <- ordered(result_time2$method, levels = c("R-JAGS", "dirinla", "long R-JAGS"))
#result_time2$N <- factor(result_time2$N)
result_time

times2 <- result_time2 %>%
  ggplot(data = .) +
  geom_point(aes(x = C, y = value, col = method, shape = method), size = 3) +
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

pdf(paste0("simulation3_times", ".pdf"), width = 6, height = 4)
times2
dev.off()





