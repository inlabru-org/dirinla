% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.dirinlaregmodel.R
\name{predict.dirinlaregmodel}
\alias{predict.dirinlaregmodel}
\title{Finding the mode of the full posterior distribution}
\usage{
\method{predict}{dirinlaregmodel}(object, data.pred.cov, ...)
}
\arguments{
\item{object}{dirinlaregmodel object.}

\item{data.pred.cov}{Data.frame with the covariate values for the variables to predict.}

\item{...}{Other arguments.}
}
\value{
model dirinlaregmodel object
}
\description{
`predict.dirinlaregmodel` computes the posterior predictive distribution for some given values of the covariates
}
\examples{
if (dirinla_safe_inla() &&
    requireNamespace("DirichletReg", quietly = TRUE)) {
### In this example, we show how to fit a model using the dirinla package ###
### --- 1. Loading the libraries --- ####
library(INLA)
library(DirichletReg)


### --- 2. Simulating from a Dirichlet likelihood --- ####
set.seed(1000)
N <- 50 #number of data
V <- as.data.frame(matrix(runif((4) * N, 0, 1), ncol = 4)) #Covariates
names(V) <- paste0('v', 1:4)

formula <- y ~ 1 + v1 | 1 + v2 | 1 + v3 | 1 + v4
(names_cat <- formula_list(formula))

x <- c(-1.5, 1, -3, 1.5,
       2, -3 , -1, 5)

mus <- exp(x) / sum(exp(x))
C <- length(names_cat)
data_stack_construct <-
  data_stack_dirich(y = as.vector(rep(NA, N * C)),
                    covariates = names_cat,
                    data       = V,
                    d          = C,
                    n          = N)

A_construct <- data_stack_construct
A_construct[1:8, ]

eta <- A_construct \%*\% x
alpha <- exp(eta)
alpha <- matrix(alpha,
                ncol  = C,
                byrow = TRUE)
y_o <- rdirichlet(N, alpha)
colnames(y_o) <- paste0("y", 1:C)
head(y_o)


### --- 3. Fitting the model --- ####
y <- y_o
model.inla <- dirinlareg(
  formula  = y ~ 1 + v1 | 1 + v2 | 1 + v3 | 1 + v4,
  y        = y,
  data.cov = V,
  prec     = 0.0001,
  verbose  = FALSE)


summary(model.inla)
### --- 4. Predicting for v1 = 0.25, v2 = 0.5, v3 = 0.5, v4 = 0.1 --- ####
model.prediction <- predict(model.inla,
                data.pred.cov= data.frame(v1 = 0.25,
                                       v2 = 0.5,
                                       v3 = 0.5,
                                       v4 = 0.1))
model.prediction$summary_predictive_means
}
}
\author{
Joaqu??n Mart??nez-Minaya <\email{jomarminaya@gmail.com}>
}
