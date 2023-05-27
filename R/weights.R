#' Compute weights
#'
#' \code{robust_weights} computes the robust weights given a data matrix \code{X}, 
#' a scale parameter \code{zeta} and a parameter that controls the weights \code{delta}.
#' Namely, the lth weight \code{w[l]} is given by 
#' \deqn{
#' w[l] = exp(-\zeta{\sum_{j\in D_{1}}(X_{i'j}-X_{ij})^2+\sum_{j\in D_{2}}{\delta}^2})
#' , where the lth pair of nodes is (\code{i},\code{i'})
#' and \code{D1={j:|X_{ij}-X_{i'j}|<delta}}, \code{D2={j:|X_{ij}-X_{i'j}|>delta}}.
#' }
#' @param X The data matrix to be clustered. The rows are observations, and the columns 
#' are features.
#' @param delta The nonnegative parameter that controls the scale of robust weights 
#' when there is outlier(s) in the data.
#' @param zeta The nonnegative parameter that controls the scale of robust weights.
#' @author Chenyu Liu, Qiang Sun, Kean Ming Tan
#' @useDynLib Rcvxclustr
#' @import gdata
#' @export
#' @return A vector \cite{wt} of weights for robust convex clustering.
robust_weights <- function(X, delta, zeta){
    sol <- robustweights(X=X,delta=delta,zeta=zeta)
    weights <- lowerTriangle(sol)
    return(weights/max(weights))
}

#' \code{GKernel_weights} computes the Gaussian kernel weights, 
#' using the distance matrix as the input. 
#' 
#' @param phi The parameter in the Gaussian kernel
#' @param distance The distance matrix, with dimension n*n, 
#' where n is the row count of the data matrix. The distance matrix can be generated 
#' by function \code{distance_matrix} from the data matrix. 
#' @export
GKernel_weights <- function(phi,distance){
  wt <- exp(-phi*distance^2)
  w<-lowerTriangle(wt)/max(lowerTriangle(wt))
  return(w)
}

#' \code{uni_weights} computes the uniform weights. 
#' 
#' @param n The row count, or amount of data points, of the data matrix. 
#' @param p The column count, or amount of features, of the data matrix. 
#' @export
uni_weights <- function(n,p){
  wt <- rep(1,n*(n-1)/2)
  return (wt)
}