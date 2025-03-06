## Rcvxclustr

**Robust Convex Clustering Algorithm Implemented in R**

***

### Description

This R package `Rcvxclustr` implements the robust convex clustering problem proposed by [Sun et al. (2025)](https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-19/issue-1/Resistant-convex-clustering--How-does-the-fusion-penalty-enhance/10.1214/25-EJS2359.full). Classical approaches towards convex clustering solves a convex optimization problem with the cost function being a squared loss plus a fusion penalty that encourages the estimated centroids for observations in the same cluster to be identical. Those approaches are not robust to arbitrary outliers, and when data are contaminated, they fail to identify the correct cluster relationships. This proposed robust convex clustering algorithm, applying Huber loss and a modified weight function, performs well in cases with outliers. It does not break down until more than half of the observations are arbitrary outliers. 

In a convex clustering problem, the input is a data matrix *X* of dimension *n*&times;*p*, with *n* is the number of samples and *p* the number of features. Convex clustering algorithm estimates a centroid matrix *U* of the same size as *X*, with *X<sub>i</sub>* and *X<sub>j</sub>* being in the same cluster if and only if *U<sub>i</sub>=U<sub>j</sub>*. 

This proposed algorithm is a modified version of that proposed by [Chi and Lange (2015)](https://arxiv.org/abs/1304.0499), implemented in the R package `cvxclustr`. The two methods are compared in [our paper](https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-19/issue-1/Resistant-convex-clustering--How-does-the-fusion-penalty-enhance/10.1214/25-EJS2359.full) and can be reproduced. 

### Installation

Install `Rcvxclustr` from GitHub: 
```r
install.packages("devtools")
library(devtools)
devtools::install_github("statsle/Rcvxclustr")
library(Rcvxclustr)
```

Helper functions can be accessed by typing `help(function name)` in R command. 

### Functions

Two main functions are implemented, and other functions in the package are dependencies. 

- `robust_weights`: This function implements the new weight function applied in fusion penalty, proposed in [our paper]([https://arxiv.org/abs/1906.09581v2](https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-19/issue-1/Resistant-convex-clustering--How-does-the-fusion-penalty-enhance/10.1214/25-EJS2359.full)). It returns the lower triangle part (in the form of vector) of the *n &times;n* weight matrix for data *X* with *n* rows (number of samples). 

- `robustcvxclust`: This function implements the alternating direction method of multipliers algorithm for solving the convex clustering objective function. It returns the matrix of centroid differences, i.e. matrix *V*, with *U<sub>i</sub>-U<sub>i'</sub>=V<sub>ii'</sub>* for all *i<i'*. 

- `create_adjacency_matrix`: This function creates the adjacency matrix *A* using matrix *V* and the weight vector. 

- `find_clusters_from_adjacency`: This function generates clustering results from the adjacency matrix *A*. 

### A simple example

We first library all packages we need: 
```r
library(clustRviz)
library(MASS)
library(Matrix)
library(igraph)
library(gdata)
library(Rcpp)
library(clues)
library(Rcvxclustr)
library(ggplot2)
library(RSNNS)
```

In the example, there are *n=50* observations that belong to two distinct non-overlapping clusters. The data matrix *X* is generated according to the model *X<sub>i</sub>=U<sub>1</sub>+&epsilon;<sub>i</sub>* when *i* belongs to the first cluster, and *X<sub>i</sub>=U<sub>2</sub>+&epsilon;<sub>i</sub>* when *i* belongs to the second cluster. Here *U<sub>1</sub>* and *U<sub>2</sub>* subject to different *20*-dimensional multivariate Gaussian distribution with identical covariance matrix and different means. We also some outliers to the data, making some dimensions in some samples deviate a lot from its neighbours. The data is generated as follows: 

```r
set.seed(1234)
N <- 50/2
p <- 20
mu1 <- mvrnorm(mu=rep(0,p),Sigma=diag(1,p))
X1 <- mvrnorm(N,mu=mu1,Sigma=diag(1,p))
mu2 <- mvrnorm(mu=rep(5,p),Sigma=diag(1,p))
X2 <- mvrnorm(N,mu=mu2,Sigma=diag(1,p))
X <- rbind(X1, X2)
n = dim(X)[1]
outliers <- sample(c(runif(n=5,min=20,max=50),runif(n=5,min=-50,max=-20)))
X[sample(1:(n*p),10)] <- outliers
```

Then create the weight vector: 
```r
wt.vec <- uni_weights(n,p)
```

Solve the convex clustering objective function: 
```r
H <- robustcvxclust(X,rho=1,tau=3,lambda=0.3,wt=wt.vec)
```

The returned `H` contains the final outcome of matrix *U*, *W*, *Y*, *V*, *Z*, the iteration time, and the final tolerance level *&epsilon;*. 

We use the output *V* (first dimension equals to the length of weight vector, and second dimension equals to *p*) to create the adjacency matrix: 
```r
A <- create_adjacency_matrix(t(H$V),wt.vec,n)
```

Finally we can obtain the clustering results by
```r
cl <- find_clusters_from_adjacency(A)$cluster
```

Using function `adjustedRand` in package `clues`, we can evaluate the performance of clustering:
```r
cl_true <- c(rep(1,N),rep(2,N))
adjustedRand(cl_true,cl)
```

In this case, the algorithm gets all the results correctly, that is, `cl` equals `cl_true`, and the `adjustedRand` all equal to 1. 

```r
> cl
 [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
[26] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
> adjustedRand(cl_true,cl)
   Rand      HA      MA      FM Jaccard 
      1       1       1       1       1 
```

### The ADMM algorithm with parameter tuning for $`\lambda`$

We use the ADMM algorithm to solve the optimization problem over an increasing sequence of $`\lambda`$ values until all data points are clustered into one single cluster. To enhance computational efficiency, we employ a warm-start strategy, initializing the ADMM algorithm for a new $`\lambda`$ with the solution obtained from the preceding $`\lambda`$. For each value of $`\lambda`$, we perform multiple iterations of the ADMM algorithm. 
We report the best solution obtained along the path of increasing $`\lambda`$’s, measured by Hubert and Arabie’s adjusted Rand index (HA Rand index). 
We use the function `ADMM_RCC` in the package to conduct the numerical experiments in our paper. 

<!---
Adopting the CARP algorithm of [Waylandt et al (2020)](https://www.tandfonline.com/doi/full/10.1080/10618600.2019.1629943), we run the clustering method while incrementing the lambda parameter. Each lambda is only used for one iteration. This algorithm can speed up the process greatly. We use the function `CARP_RCC` in the package to conduct the numerical experiments in our paper. 
-->

### Synthetic data

In our numerical studies, we show the robustness of our proposed method using the synthetic data generated as follows: 

```r
data.gen.mixed <- function(seed,N,p,out_form='entry',out_entry_prop,
                           noise_type='normal',outlier_type='uniform',df=NULL){
  set.seed(seed)
  
  # normal noise with uniform outliers
  if (noise_type=='normal' & outlier_type=='uniform'){
    mu1 <- rnorm(p,0,1)
    X1 <- mvrnorm(N,mu=mu1,Sigma=diag(1,p))
    mu2 <- c(rnorm(p/2,3,1),rnorm(p/2,-3,1))
    X2 <- mvrnorm(N,mu=mu2,Sigma=diag(1,p))
    X <- rbind(X1, X2)
    X0<-X #clean data
    n = dim(X)[1]
    if (out_form=='entry'){
      out_num <- as.integer(n*p*out_entry_prop)
      outliers <- runif(n=out_num,min=10,max=20)
      if (out_num > 0){
        X[sample(1:(n*p),out_num)] <- outliers
      }
    }
    else if (out_form=='row'){
      out_row <- as.integer(n*out_entry_prop)
      out_rows <- sample(1:n,out_row)
      out_num <- as.integer(p*0.2)
      for (row in out_rows){
        out_cols<-sample(1:p,out_num)
        outliers <- runif(n=out_num,min=10,max=20)
        X[row,out_cols]<-outliers
      }
    } 
    else {
      stop("The contamination form needs to be either 'entry' or 'row'.")
    }
  }
  # t noise with uniform outliers
  else if (noise_type=='t' & outlier_type=='uniform'){
    mu1 <- rnorm(p,0,1)
    X1 <- matrix(rep(mu1,each=N),nrow=N,ncol=p)+matrix(rt(n=N*p,df=df,ncp=0),nrow=N,ncol=p)
    mu2 <- c(rnorm(p/2,3,1),rnorm(p/2,-3,1))
    X2 <- matrix(rep(mu2,each=N),nrow=N,ncol=p)+matrix(rt(n=N*p,df=df,ncp=0),nrow=N,ncol=p)
    X <- rbind(X1, X2)
    X0<-X #clean data
    n = dim(X)[1]
    if (out_form=='entry'){
      out_num <- as.integer(n*p*out_entry_prop)
      outliers <- runif(n=out_num,min=10,max=20)
      if (out_num > 0){
        X[sample(1:(n*p),out_num)] <- outliers
      }
    }
    else if (out_form=='row'){
      out_row <- as.integer(n*out_entry_prop)
      out_rows <- sample(1:n,out_row)
      out_num <- as.integer(p*0.2)
      for (row in out_rows){
        out_cols<-sample(1:p,out_num)
        outliers <- runif(n=out_num,min=10,max=20)
        X[row,out_cols]<-outliers
      }
    } 
    else {
      stop("The contamination form needs to be either 'entry' or 'row'.")
    }
  }
  # normal noise with t outliers
  else if (noise_type=='normal' & outlier_type=='t'){
    mu1 <- rnorm(p,0,1)
    X1 <- mvrnorm(N,mu=mu1,Sigma=diag(1,p))
    mu2 <- c(rnorm(p/2,3,1),rnorm(p/2,-3,1))
    X2 <- mvrnorm(N,mu=mu2,Sigma=diag(1,p))
    X <- rbind(X1, X2)
    X0<-X #clean data
    n = dim(X)[1]
    if (out_form=='entry'){
      out_num <- as.integer(n*p*out_entry_prop)
      outliers <- rt(n=out_num,df=df,ncp=0)
      if (out_num > 0){
        X[sample(1:(n*p),out_num)] <- outliers
      }
    }
    else if (out_form=='row'){
      out_row <- as.integer(n*out_entry_prop)
      out_rows <- sample(1:n,out_row)
      out_num <- as.integer(p*0.2)
      for (row in out_rows){
        out_cols<-sample(1:p,out_num)
        outliers <- rt(n=out_num,df=df,ncp=0)
        X[row,out_cols]<-outliers
      }
    } 
    else {
      stop("The contamination form needs to be either 'entry' or 'row'.")
    }
  } 
  else{
    stop("The noise type or outlier type is not recognized.")
  }
  
  cl_true <- c(rep(1,N),rep(2,N))
  return (list(X0=X0,X=X,cl_true=cl_true))
}
```

To generate synthetic data with Gaussian noise and uniform outliers with 2\% entry-wise contamination: 
```r
gen <- data.gen.mixed(seed=2024,N=20,p=20,out_form='entry',out_entry_prop=0.02,
                      noise_type='normal',outlier_type='uniform')
X <- gen$X
cl_true <- gen$cl_true
```

To generate synthetic data with t noise (with df=5) and uniform outliers with 10\% entry-wise contamination: 
```r
gen <- data.gen.mixed(seed=2024,N=20,p=20,out_form='entry',out_entry_prop=0.1,
                      noise_type='t',outlier_type='uniform',df=5)
X <- gen$X
cl_true <- gen$cl_true
```

To generate synthetic data with Gaussian noise and t outliers (with df=1) with 50\% row-wise contamination: 
```r
gen <- data.gen.mixed(seed=2024,N=20,p=20,out_form='row',out_entry_prop=0.5,
                      noise_type='normal',outlier_type='t',df=1)
X <- gen$X
cl_true <- gen$cl_true
```


### Results

Using the cvxclustr method with uniform weights, we obtain:

```r
gen <- data.gen.mixed(seed=2024,N=10,p=20,out_form='entry',out_entry_prop=0.02,
                      noise_type='normal',outlier_type='uniform')
X <- gen$X
cl_true <- gen$cl_true
result <- ADMM_RCC(X=X,phi=.1,method='cvx_uni',lam.begin=0.01,lam.step=1.05,
                   rho=1,tau=NULL,cl_true=cl_true,randmode='HA',
                   max.log=200,max.iter_cvx=100,max.iter_Rcvx=NULL)
result$cl_est
[1]  1  1  1  2  3  1  1  1  1  1  4  5  6  7  8  9 10 11 12 13
```

Using the proposed Rcvxclustr method with uniform weights and $`\tau=0.1`$, we obtain: 

```r
gen <- data.gen.mixed(seed=2024,N=10,p=20,out_form='entry',out_entry_prop=0.02,
                      noise_type='normal',outlier_type='uniform')
X <- gen$X
cl_true <- gen$cl_true
result <- ADMM_RCC(X=X,phi=.1,method='Rcvx_uni',lam.begin=0.01,lam.step=1.05,
                   rho=1,tau=0.1,cl_true=cl_true,randmode='HA',
                   max.log=200,max.iter_cvx=NULL,max.iter_Rcvx=100)
result$cl_est
[1] 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2
```

The path graph can be obtained from `result$path`, and some pictures are shown in the paper. 

<!---
When trying to add student t distribution noise, the data generation function ```data.gen.mixed``` could be replaced with the below:

```r
data.gen.mixed <- function(seed,N,p,out_entry_prop,df,out_form='arbitrary',noise_type='uniform'){
  set.seed(seed)
  mu1 <- rnorm(p,0,1)
  X1 <- mvrnorm(N,mu=mu1,Sigma=diag(1,p))
  mu2 <- c(rnorm(p/2,3,1),rnorm(p/2,-3,1))
  X2 <- mvrnorm(N,mu=mu2,Sigma=diag(1,p))
  X <- rbind(X1, X2)
  X0<-X
  n = dim(X)[1]
  if (out_form == 'arbitrary'){
    out_num <- as.integer(n*p*out_entry_prop)
    outliers <- rt(n=out_num,df=df,ncp=0)
    if (out_num > 0){
      X[sample(1:(n*p),out_num)] <- outliers
    }
  }
  else{
    out_row <- as.integer(n*out_entry_prop)
    out_rows <- sample(1:n,out_row)
    out_sample <- c()
    for (row in out_rows){
      out_num <- as.integer(p*0.2)
      out_cols<-sample(1:p,out_num)
      out_sample<-c(out_sample,out_cols+(row-1)*p)
      outliers <- rt(n=length(out_sample),df=df,ncp=0)
      X[out_sample]<-outliers
    }
  }
  cl_true <- c(rep(1,N),rep(2,N))
  return (list(X0=X0,X=X,cl_true=cl_true))
}
```
-->











