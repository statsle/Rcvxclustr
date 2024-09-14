#' Robust Convex Clustering
#'
#' \code{ADMM_RCC} performs convex clustering using the ADMM algorithm over an increasing sequence of \code{lambda} with warm start. 
#' The starting value and increment step of \code{lambda} can be set. 
#' \itemize{
#' \item{n is the number of data observations}
#' \item{p is the number of features}
#' \item{nK is the number non-zero weights.}
#' }
#'
#' @param X The n-by-p data matrix whose rows are being clustered.
#' @param phi The parameter in the Gaussian kernel weights. 
#' @param method The method to be used. 
#' Choices are \code{cvx_uni}, \code{cvx_GK}, \code{Rcvx_uni} and \code{Rcvx_GK}. 
#' @param lam.begin The starting value of \code{lambda}. 
#' @param lam.step The increment step of \code{lambda}. Each iteration will mutiply \code{lambda} by this step value. 
#' @param rho Augmented Lagrangian penalty parameter.
#' @param tau The robustification parameter in huber loss.
#' @param cl_true The true clustering results. Used for rand index calculation. 
#' @param randmode The rand index mode. See \code{adjustedRand} in package \code{clues} for details. 
#' @param max.log The maximal number of iterations for \code{lambda}. 
#' @param max.iter_cvx The maximal number of iterations for the cvxclustr method for each \code{lambda}. 
#' @param max.iter_Rcvx The maximal number of iterations for the Rcvxclustr method for each \code{lambda}. 
#' The algorithm also stops when the present iteration gives out the result where all data points are classified in the same cluster. 
#'
#' @return \code{method} The method used. 
#' @return \code{lam} The best \code{lambda} value, which reaches the best rand index. 
#' @return \code{rand} The rand index obtained using the best \code{lambda}.
#' @return \code{U_est} The estimated centroid matrix obtained using the best \code{lambda}.
#' @return \code{path} A matrix, of which each row represents the clustering result for each iteration. 
#' @return \code{cl_est} The result of cluster estimation which produces the best rand index.  
#' @import cvxclustr
#' @import clues
#' @export
#' @useDynLib Rcvxclustr

ADMM_RCC <- function(X,phi,method,lam.begin,lam.step,rho,tau,cl_true,randmode,max.log,max.iter_cvx,max.iter_Rcvx){
  n <- dim(X)[1]
  p <- dim(X)[2]
  lams <- vector(length=max.log)
  rands <- vector(length=max.log)
  cl.matrix <- matrix(0,ncol=n,nrow=length(lams))
  U <- vector("list", length=max.log)
  d <- distance_matrix(X)
  wt.GK <- GKernel_weights(phi = phi,distance=d)
  wt.uni <- uni_weights(n,p)
  if (method=='cvx_uni' | method=='cvx_GK'){
    temp <- create_clustering_problem(p,n,method = "admm")
    ix<-temp$ix;M1<-temp$M1;M2<-temp$M2;s1<-temp$s1;s2<-temp$s2
    k <- length(wt.uni)
    Lambda<-matrix(0,p,k)
    if (method=='cvx_uni'){
      present.weight <- wt.uni
    }
    else if (method=='cvx_GK'){
      present.weight <- wt.GK
    }
    Eric.init <- cvxclust_admm(X=t(X),Lambda=Lambda,ix,M1,M2,s1,s2,w=present.weight,
                               gamma=lam.begin,nu=rho,accelerate = FALSE,max_iter=max.iter_cvx)
    B <- create_adjacency(Eric.init$V,present.weight,n,method="admm")
    present.cl <- find_clusters(B)$cluster
    cl.matrix[1,] <- present.cl
    rands[1] <- adjustedRand(cl_true,present.cl,randmode)
    ind <- 2
    present.lam <- lam.begin
    U[[1]] <- t(Eric.init$U)
    lams[1] <- present.lam
    while(1){
      present.lam <- present.lam * lam.step
      lams[ind] <- present.lam
      Eric.res <- cvxclust_admm(X=t(X),Lambda=Eric.init$Lambda,ix,M1,M2,s1,s2,w=present.weight,
                                gamma=present.lam,nu=rho,accelerate = FALSE,max_iter=max.iter_cvx)
      U[[ind]] <- t(Eric.res$U)
      B <- create_adjacency(Eric.res$V,present.weight,n,method="admm")
      present.cl <- find_clusters(B)$cluster
      cl.matrix[ind,] <- present.cl
      rands[ind] <- adjustedRand(cl_true,present.cl,randmode)
      Eric.init <- Eric.res
      if (length(unique(present.cl))==1 | ind == max.log){
        break
      }
      ind <- ind + 1
    }
    cls_num <- apply(cl.matrix,1,function(x) length(unique(x)))
  }
  else if (method=='Rcvx_uni' | method=='Rcvx_GK'){
    present.lam <- lam.begin
    if (method=='Rcvx_uni'){
      present.weight <- wt.uni
    }
    else if (method=='Rcvx_GK'){
      present.weight <- wt.GK
    }
    res.init <- robustcvxclust(X,lambda=present.lam,tau=tau,wt=present.weight,rho=rho,max_iter=max.iter_Rcvx)
    A <- create_adjacency(t(res.init$V),present.weight,n,method="admm")
    present.cl <- find_clusters(A)$cluster
    cl.matrix[1,] <- present.cl
    rands[1] <- adjustedRand(cl_true,present.cl,randmode)
    U[[1]] <- res.init$U
    lams[1] <- present.lam
    ind <- 2
    while(1){
      present.lam <- present.lam * lam.step
      lams[ind] <- present.lam
      res <- robustcvxclust(X,V=res.init$V,Y=res.init$Y,W=res.init$W,Z=res.init$Z, 
                            lambda=present.lam,tau=tau,wt=present.weight,rho=rho,max_iter=max.iter_Rcvx) 
      U[[ind]] <- res$U
      A <- create_adjacency(t(res$V),present.weight,n,method="admm")
      present.cl <- find_clusters(A)$cluster
      cl.matrix[ind,] <- present.cl
      rands[ind] <- adjustedRand(cl_true,present.cl,randmode)
      res.init <- res
      if (length(unique(present.cl))==1 | ind == max.log){
        break
      }
      ind <- ind+1
    }
    cls_num <- apply(cl.matrix,1,function(x) length(unique(x)))
  }
  # Final presentation
  fake_num <- ifelse(cls_num < 2, n, cls_num)
  true_ind <- which(fake_num==min(fake_num))
  best_rand <- max(rands)
  best_ind <- which(rands==best_rand)
  best_lam <- lams[best_ind]
  max_best_ind <- which(lams==max(best_lam))
  max_best_lam <- max(best_lam)
  max_best_U <- U[[max_best_ind]]
  return(list(method=method,lam=max_best_lam,rand=best_rand,U_est=max_best_U,
              path=cl.matrix,cl_est=cl.matrix[max_best_ind,]))
}


