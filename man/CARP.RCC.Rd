% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CARP_RCC.R
\name{CARP.RCC}
\alias{CARP.RCC}
\title{Robust Convex Clustering using CARP}
\usage{
CARP.RCC(
  X,
  phi,
  method,
  lam.begin,
  lam.step,
  rho,
  tau,
  cl_true,
  randmode,
  max.log = 100
)
}
\arguments{
\item{X}{The n-by-p data matrix whose rows are being clustered.}

\item{phi}{The parameter in the Gaussian kernel weights.}

\item{method}{The method to be used. 
Choices are \code{cvx_uni}, \code{cvx_GK}, \code{Rcvx_uni} and \code{Rcvx_GK}.}

\item{lam.begin}{The starting value of \code{lambda}.}

\item{lam.step}{The increment step of \code{lambda}. Each iteration will mutiply \code{lambda} by this step value.}

\item{rho}{Augmented Lagrangian penalty parameter.}

\item{tau}{The robustification parameter in huber loss.}

\item{cl_true}{The true clustering results. Used for rand index calculation.}

\item{randmode}{The rand index mode. See \code{adjustedRand} in package \code{clues} for details.}

\item{max.log}{The maximal number of iterations. 
The algorithm also stops when the present iteration gives out the result where all data points are classified in the same cluster.}
}
\value{
\code{method} The method used.

\code{rand} The best rand index obtained.

\code{lam} The best \code{lambda} value, which reaches the best rand index.

\code{path} A matrix, of which each row represents the clustering result for each iteration.

\code{cl_est} The result of cluster estimation which produces the best rand index.
}
\description{
\code{CARP_RCC} performs convex clustering with the CARP algorithm. 
The starting value and increment step of \code{lambda} can be set. 
\itemize{
\item{n is the number of data observations}
\item{p is the number of features}
\item{nK is the number non-zero weights.}
}
}
