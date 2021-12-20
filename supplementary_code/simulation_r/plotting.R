


### --- 1. Loading packages --- ####
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


#### In this script, we make beautiful plots for the results obtained in simulation 4

setwd("supplementary_code/simulation_r")
file_list <- list.files(pattern = ".RDS")


### Reading computational times and accuracy of all the models
res <- readRDS("simulation4_50-500.RDS")
res1 <- readRDS("simulation4_1000.RDS")
res <- cbind(res, res1)

### --- 1. Loading the data simulated --- ###
### Parameters fitted
tau0 <- 1
x <- c(-1.5, 2,
       1, -3)
levels_factor <- c(2, 5, 10)

#random effect
prec_w <- c(4, 9)
(sd_w <- 1/sqrt(prec_w))





### --- Function to plot all the simulations done in simulation 4 --- ####
plotting_all <- function(n, levels_factor)
{
  if(is.na(levels_factor)){
    levels_factor <- n
  }
  cat(n, "-", levels_factor, "\n")
  ### --- 2. Reading models --- ####
  model.jags <- readRDS(paste0("model_jags_", n, "_", levels_factor, ".RDS"))
  model.jags.2 <- readRDS(paste0("model_jags_long_", n, "_", levels_factor, ".RDS"))
  model.inla <- readRDS(paste0("model_inla_", n, "_", levels_factor, ".RDS"))
  model.inla.2 <- readRDS(paste0("model_inla_2_", n, "_", levels_factor, ".RDS"))

  ### --- 3. Plotting slopes --- ####
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
                         labels = c("R-JAGS", "dirinla pc", "long R-JAGS", "dirinla hn"))

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
      scale_linetype_manual(labels=c("R-JAGS", "dirinla pc", "long R-JAGS", "dirinla hn"),
                            values=c("dotted", "twodash",  "solid", "dashed"))

    #values=c("dotted", "twodash",  "solid", "longdash"))

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

  pdf(paste0("examples_simulation4_slopes_", n, "_", levels_factor, ".pdf"), width = 15, height = 4)
  gridExtra::grid.arrange(p2[[1]], p2[[2]], p2[[3]], p2[[4]],ncol = 4)
  dev.off()

  ### --- 4. Hyperparemeters --- ####
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


  #jags1
  #logSigma1 and logSigma2
  lapply(c("sigma1", "sigma2"), function(x){
    model.jags$BUGSoutput$sims.matrix[,c(x)] %>%
      #.^(-1/2) %>%
      log(.) %>%
      density(., adjust = 2) %>%
      .[1:2] %>%
      do.call(data.frame, .) -> dens
    colnames(dens) <- c("x", "y")
    dens
  }) -> dens_log
  names(dens_log) <- c("log(sigma1)", "log(sigma2)")



  lapply(c("sigma1", "sigma2"), function(x){
    model.jags$BUGSoutput$sims.matrix[,c(x)] %>%
      #.^(-1/2) %>%
      log(.)%>%
      density(., adjust = 2) %>%
      .[1:2] %>%
      inla.tmarginal(fun = function(x)exp(x), .) %>%
      do.call(data.frame, .) -> dens
    colnames(dens) <- c("x", "y")
    dens
  }) -> dens
  names(dens) <- c("sigma1", "sigma2")


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

  dens <- c(dens, dens_log)
  dens2 <- c(dens2, dens_log2)
  inla_sigma <- c(inla_sigma, inla_sigma_log)
  inla_sigma_2 <- c(inla_sigma_2, inla_sigma_log_2)


  sd_name <- list(expression(paste("p(", sigma[1], "|", "y)")),
                  expression(paste("p(", sigma[2], "|", "y)")),
                  expression(paste("p(", "log(", sigma[1],")", "|", "y)")),
                  expression(paste("p(", "log(", sigma[2],")", "|", "y)")))
  sd_name2 <- list(expression(sigma[1]),
                   expression(sigma[2]),
                   expression(paste("log(", sigma[1], ")")),
                   expression(paste("log(", sigma[2], ")")))


  p3 <- list()
  for(i in 1:length(dens))
  {
    #Data combining jags (1) and inla (2)
    dens_sigma <- rbind(cbind(dens[[i]], group = 1),
                        cbind(as.data.frame(inla_sigma[[i]]), group = 2),
                        cbind(dens2[[i]], group = 3),
                        cbind(as.data.frame(inla_sigma_2[[i]]), group = 4))

    #Adding priors on standar deviation
    hn_inla <- function(sigma)
    {
      tau0 = tau0
      log_dens = log(2) - 0.5 * log(2 * pi) + 0.5 * log(tau0)
      log_dens = log_dens - 0.5 * tau0 * sigma^2
      return(exp(log_dens))
    }



    y_pc_log_sigma <- log(1/sqrt(inla.pc.rprec(100000000, 10, 0.01))) %>% density(.)
    y_pc_log_sigma <- data.frame(x = y_pc_log_sigma$x, y = y_pc_log_sigma$y)
    y_pc_sigma <- inla.tmarginal(function(x)exp(x), y_pc_log_sigma)



    priors_df <- data.frame(x = c(y_pc_sigma$x, y_pc_sigma$x),
                            y = c(y_pc_sigma$y, hn_inla(y_pc_sigma$x)),
                            group = as.factor(rep(c(5,6),
                                                  c(length(y_pc_sigma$x), length(y_pc_sigma$x)))))

    if(i <=2)
    {
      dens_sigma <- rbind(dens_sigma, priors_df)

      ### Assigning names
      dens_sigma$group <- factor(dens_sigma$group,
                                 labels = c("R-JAGS", "dirinla pc", "long R-JAGS", "dirinla hn", "pc-prior", "hn-prior"))
    }else{ #Not the optimum way
      priors_df_log <- lapply(levels(priors_df$group), function(b){
        priors_sim <- priors_df %>% dplyr::filter(group == b) %>% dplyr::select(x,y) %>% as.matrix(.) %>%
          inla.rmarginal(1000000, .) %>% log(.)
        priors_sim <- priors_sim[which(!is.na(priors_sim))] %>% density(., adjust = 2)
        priors_sim <- data.frame(x = priors_sim$x, y = priors_sim$y, group = as.numeric(b)) %>% as.matrix(.)
        priors_sim
      }) %>% do.call(rbind, .)

      dens_sigma <- rbind(dens_sigma, priors_df_log)
      dens_sigma$group <- factor(dens_sigma$group,
                                 labels = c("R-JAGS", "dirinla pc", "long R-JAGS", "dirinla hn", "pc-prior", "hn-prior"))

    }



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
    #xmax <- 2


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
        values= c("darkgreen", "red4", "blue4", "orange2", "deeppink", "magenta2" )) +
      scale_linetype_manual(labels=c("R-JAGS", "dirinla pc", "long R-JAGS", "dirinla hn", "pc-prior", "hn-prior"),
                            values=c("dotted", "twodash",  "solid", "dashed", "longdash", "12345678"))

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


  pdf(paste0("examples_simulation4_sigma_", n, "_", levels_factor, ".pdf"), width = 15, height = 4)
  gridExtra::grid.arrange(p3[[1]], p3[[2]], p3[[3]], p3[[4]],ncol = 4)
  dev.off()

}



### --- 2. Applying function to plot --- ####
n <- c(50, 100, 500, 1000)
levels_factor <- c(2, 5, 10, 25, NA)

arguments <- expand.grid(n, levels_factor)
n <- arguments[,1]
levels_factor <- arguments[,2]

a <- mapply(plotting_all,
            n = n,
            levels_factor = levels_factor)





###################################################################

### --- 3. Tables for the paper --- ####
extract_tables <- function(levels_factor)
{
  n <- c(50, 100, 500, 1000)


  pos <- paste0(n, "-", levels_factor)


  #Computational times

  res[,pos[1]]$times
  res[,pos[2]]$times


  1:4 %>% sapply(., function(x){res[,pos[x]]$times}) %>% t(.) -> result_time
  colnames(result_time) <- c("R-JAGS", "dirinla pc", "long R-JAGS", "dirinla hn")
  rownames(result_time) <- paste0( c(50, 100, 500, 1000))
  result_time



  ### --- Ratio 1 --- ####
  ### Betas for pc and hn
  1:4 %>% sapply(., function(x){round(res[,pos[x]]$ratio1_beta1_pc,4)}) %>% t(.) -> ratio1_beta1_pc
  colnames(ratio1_beta1_pc) <- paste0("beta1-pc", 1:4)
  rownames(ratio1_beta1_pc) <- paste0( c(50, 100, 500, 1000))

  1:4 %>% sapply(., function(x){round(res[,pos[x]]$ratio1_beta1_hn,4)}) %>% t(.) -> ratio1_beta1_hn
  colnames(ratio1_beta1_hn) <- paste0("beta1-hn", 1:4)
  rownames(ratio1_beta1_hn) <- paste0( c(50, 100, 500, 1000))

  ### Sigma and log_sigma for pc and hn
  1:4 %>% sapply(., function(x){round(res[,pos[x]]$ratio1_sigma_pc, 4)}) %>% t(.) -> ratio1_sigma_pc
  colnames(ratio1_sigma_pc) <- paste0("sigma-pc", 1:2)
  rownames(ratio1_sigma_pc) <- paste0( c(50, 100, 500, 1000))

  1:4 %>% sapply(., function(x){round(res[,pos[x]]$ratio1_sigma_hn, 4)}) %>% t(.) -> ratio1_sigma_hn
  colnames(ratio1_sigma_hn) <- paste0("sigma-hn", 1:2)
  rownames(ratio1_sigma_hn) <- paste0( c(50, 100, 500, 1000))

  ### ratio1 pc
  result_ratio1_pc <- cbind(ratio1_beta1_pc,
                            ratio1_sigma_pc)

  result_ratio1_pc

  ### ratio1 hn
  result_ratio1_hn <- cbind(ratio1_beta1_hn,
                            ratio1_sigma_hn)


  ### --- Ratio 2 --- ####
  ### Betas for pc and hn
  1:4 %>% sapply(., function(x){round(res[,pos[x]]$ratio2_beta1_pc,4)}) %>% t(.) -> ratio2_beta1_pc
  colnames(ratio2_beta1_pc) <- paste0("beta1-pc", 1:4)
  rownames(ratio2_beta1_pc) <- paste0( c(50, 100, 500, 1000))

  1:4 %>% sapply(., function(x){round(res[,pos[x]]$ratio2_beta1_hn,4)}) %>% t(.) -> ratio2_beta1_hn
  colnames(ratio2_beta1_hn) <- paste0("beta1-hn", 1:4)
  rownames(ratio2_beta1_hn) <- paste0( c(50, 100, 500, 1000))

  ### Sigma and log_sigma for pc and hn
  1:4 %>% sapply(., function(x){round(res[,pos[x]]$ratio2_sigma_pc, 4)}) %>% t(.) -> ratio2_sigma_pc
  colnames(ratio2_sigma_pc) <- paste0("sigma-pc", 1:2)
  rownames(ratio2_sigma_pc) <- paste0( c(50, 100, 500, 1000))

  1:4 %>% sapply(., function(x){round(res[,pos[x]]$ratio2_sigma_hn, 4)}) %>% t(.) -> ratio2_sigma_hn
  colnames(ratio2_sigma_hn) <- paste0("sigma-hn", 1:2)
  rownames(ratio2_sigma_hn) <- paste0( c(50, 100, 500, 1000))

  ### ratio2 pc
  result_ratio2_pc <- cbind(ratio2_beta1_pc,
                            ratio2_sigma_pc)

  ### ratio2 hn
  result_ratio2_hn <- cbind(ratio2_beta1_hn,
                            ratio2_sigma_hn)



  #Latex
  library(xtable)
  list(times = xtable(result_time, digits = 4),
  ratio1_pc = xtable(result_ratio1_pc, digits = 4),
  ratio1_hn = xtable(result_ratio1_hn, digits = 4),
  ratio2_pc = xtable(result_ratio2_pc, digits = 4),
  ratio2_hn = xtable(result_ratio2_hn, digits = 4))

}
extract_tables(2)


