#' Robust Convex Clustering
#'
#' \code{robustcvxclust} performs robust convex clustering using ADMM. This is an R wrapper function around C code.
#' Dimensions of various arguments are as follows:
#' \itemize{
#' \item{n is the number of data observations}
#' \item{p is the number of features}
#' \item{nK is the number non-zero weights.}
#' }
#'
#' @param X The n-by-p data matrix whose rows are being clustered.
#' @param W The n-by-p matrix of Lagrange multipliers.
#' @param V The centroid difference matrix.
#' @param Y The nK-by-p matrix of Lagrange multipliers.
#' @param Z The n-by-p matrix of Lagrange multipliers.
#' @param max_iter The maximum number of iterations. The default value is 1e5.
#' @param rho Augmented Lagrangian penalty parameter.
#' @param tau The robustification parameter in huber loss.
#' @param lambda The regularization parameter controlling the amount of shrinkage.
#' @param wt A vector of nK positive weights.
#' @param tol_abs The convergence tolerance. The default value is 1e-05.
#'
#' @return \code{U} The centroid matrix.
#' @return \code{W} The centroid matrix.
#' @return \code{V} The centroid difference matrix.
#' @return \code{Y} The Lagrange multiplier matrix.
#' @return \code{Z} The Lagrange multiplier matrix.
#' @return \code{iter} The number of iterations taken.
#' @return \code{tol} The residual tolerances.
#' @import MASS
#' @export
#' @author  Chenyu Liu, Qiang Sun, Kean Ming Tan
#' @useDynLib Rcvxclustr

robustcvxclust <- function(X,W=NULL,V=NULL,Y=NULL,Z=NULL,max_iter=1e5,
                           rho,tau,lambda,wt,tol_abs=1e-05){
    n <- as.integer(nrow(X))
    p <- as.integer(ncol(X))
    E <- create_E_matrix(n)$E
    nK <- nrow(E)
    if(is.null(V)==TRUE){
        V <- matrix(1,nrow=nK,ncol=p)
        #V <- matrix(rnorm(nK*p),nrow=nK,ncol=p)
    }
    if(is.null(W)==TRUE){
        W <- matrix(1,nrow=n,ncol=p)
        #W <- matrix(rnorm(n*p),nrow=n,ncol=p)
    }
    if(is.null(Y)==TRUE){
        Y <- matrix(0,nrow=nK,ncol=p)
        #Y <- matrix(rnorm(nK*p),nrow=nK,ncol=p)
    }
    if(is.null(Z)==TRUE){
        Z <- matrix(0,nrow=n,ncol=p)
        #Z <- matrix(rnorm(n*p),nrow=n,ncol=p)
    }
    sol = robust_convex_cluster(X=X,W=W,V=V,Y=Y,Z=Z,E=E,max_iter=max_iter,
                                tol_abs=tol_abs,lambda=lambda,rho=rho,tau=tau,wt=wt)
    return(list(U=sol$U,W=sol$W,V=sol$V,Y=sol$Y,Z=sol$Z,iter=sol$iter,tol=sol$tol))
}


