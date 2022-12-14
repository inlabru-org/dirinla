% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/H_matrix_eta_diag.R
\name{H_matrix_eta_diag}
\alias{H_matrix_eta_diag}
\title{Computing additional diagonal part for the real Hessian H = H0 + diag}
\usage{
H_matrix_eta_diag(eta, d, y)
}
\arguments{
\item{eta}{eta vector to compute the expected Hessian.}

\item{d}{Dimension}

\item{y}{Data corresponding to the i-individual}
}
\value{
Elements of the diagonal  such as H = H0 + diag
}
\description{
`H_matrix_eta_diag` computes the expected Hessian in eta of -loglikelihood
}
\author{
Joaquín Martínez-Minaya <\email{jomarminaya@gmail.com}>
}
