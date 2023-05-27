## Rcvxclustr

**Robust Convex Clustering Algorithm Implemented in R**

***

### Description

This package implements the robust convex clustering problem proposed by [Liu, Sun and Tan (2019)](https://arxiv.org/abs/1906.09581v2). Classical approaches towards convex clustering solves a convex optimization problem with the cost function being a squared loss plus a fusion penalty that encourages the estimated centroids for observations in the sames cluster to be identical. Those approaches are not robust to arbitrary outliers, and when data are contaminated, they fail to identify the correct cluster relationships. This proposed robust convex clustering algorithm, applying Huber loss and a modified weight function, performs well in cases with outliers. It does not break down until more than half of the observations are arbitrary outliers. 

In a convex clustering problem, the input is a data matrix *X* of dimention *n*&times;*p*, with *n* is the number of samples and *p* the number of features. Convex clustering algorithm estimates a centroid matrix *U* of the same size as *X*, with *X<sub>i</sub>* and *X<sub>j</sub>* being in the same cluster if and only if *U<sub>i</sub>=U<sub>j</sub>*. 

This proposed algorithm is a modified version of that proposed by [Chi and Lange (2015)](https://arxiv.org/abs/1304.0499), implemented in package `cvxclustr`. The two packages are compared in [this paper](https://arxiv.org/abs/1906.09581v2) and can be reproduced. 

### Installation

Install `Rcvxclustr` from GitHub: 
```r
install.packages("devtools")
library(devtools)
devtools::install_github("JhZhang-1999/Rcvxclustr")
library(Rcvxclustr)
```

Helper functions can be accessed by typing `help(function name)` in R command. 

### Functions

Two main functions are implemented, and other functions in the package are dependencies. 

- `robust_weights`: This function implements the new weight function applied in fusion penalty, proposed in [the paper](https://arxiv.org/abs/1906.09581v2). It returns the lower triangle part (in the form of vector) of the *n &times;n* weight matrix for data *X* with *n* rows (number of samples). 

- `robustcvxclust`: This function implements the alternating direction method of multipliers algorithm for solving the convex clustering objective function. It returns the matrix of centriod differences, i.e. matrix *V*, with *U<sub>i</sub>-U<sub>i'</sub>=V<sub>ii'</sub>* for all *i<i'*. 

- `create_adjacency_matrix`: This function creates the adjacency matrix *A* using matrix *V* and the weight vector. 

- `find_clusters_from_adjacency`: This function generates clustering results from the adjacency matrix *A*. 

### A Simple Example

We first library all packages we need: 
```r
library(MASS)
library(Matrix)
library(igraph)
library(gdata)
library(Rcpp)
library(clues)
library(Rcvxclustr)
```

In the simulation, there are *n=20* observations that belong to two distinct non-overlapping clusters. The data matrix *X* is generated according to the model *X<sub>i</sub>=U<sub>1</sub>+&epsilon;<sub>i</sub>* when *i* belongs to the first cluster, and *X<sub>i</sub>=U<sub>2</sub>+&epsilon;<sub>i</sub>* when *i* belongs to the second cluster. Here *U<sub>1</sub>* and *U<sub>2</sub>* subject to different *20*-dimensional multivariate Gaussian distribution with identical covariance matrix and different means. We also add outliers to the data, making some dimensions in some samples deviate a lot from its neighbors. The data is generated as follows: 

```r
set.seed(1234)
N <- 25
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
wt.vec <- uni_weights(N,p)
```

Solve the convex clustering objective function: 
```r
H <- robustcvxclust(X,rho=1,tau=3,lambda=0.7,wt=wt.vec)
```

The returned `H` contains the final outcome of matrix *U*, *W*, *Y*, *V*, *Z*, the iteration time, and the final tolerance level *&epsilon;*. 

We use the output *V* (first dimension equals to the length of weight vector, and second dimension equals to *p*) to create the adjacency matrix: 
```r
A <- create_adjacency_matrix(t(H$V),wt,n)
```

Finally we can obtain the clustering results by
```r
cl <- find_clusters_from_adjacency(A)$cluster
```

Using function `adjustedRand` in package `clues`, we can evaluate the performance of clustering:
```r
cl_true <- c(rep(1,25),rep(2,25))
adjustedRand(cl_true,cl)
```

In this case, the algorithm gets all the results correctly, that is, `cl` equals `cl_true`, and the `adjustedRand` all equal to 1. 

```r
> cl
 [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2
[28] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
> adjustedRand(cl_true,cl)
   Rand      HA      MA      FM Jaccard 
      1       1       1       1       1 
```

### The CARP Algorithm

Adopting the CARP alrogithm of [Waylandt et al (2020)](https://www.tandfonline.com/doi/full/10.1080/10618600.2019.1629943), we run the clustering method while incrementing the lambda paramter. Each lambda is only used for one iteration. This algorithm can speed up the process greatly. We use the function `CARP_RCC` in the package to conduct the numerical experiments in our paper. 

### Artificial Demonstration Data

In the numerical section, we show the robustness of our proposed method using an artificial datag generated as follows: 

```r
data.gen.mixed <- function(seed,N,p,out_entry_prop,out_form='arbitrary'){
  set.seed(seed)
  mu1 <- rnorm(p,0,1)
  X1 <- mvrnorm(N,mu=mu1,Sigma=diag(1,p))
  mu2 <- c(rnorm(p/2,3,1),rnorm(p/2,-3,1))
  X2 <- mvrnorm(N,mu=mu2,Sigma=diag(1,p))
  X <- rbind(X1, X2)
  n = dim(X)[1]
  if (out_form == 'arbitrary'){
    out_num <- as.integer(n*p*out_entry_prop)
    outliers <- runif(n=out_num,min=10,max=20)
    if (out_num > 0){
      X[sample(1:(n*p),out_num)] <- outliers
    }
  }
  cl_true <- c(rep(1,N),rep(2,N))
  return (list(X=X,cl_true=cl_true))
}
seed <- 1262
gen <- data.gen.mixed(seed,N=12,p=20,out_entry_prop = 0.01)
X <- gen$X
cl_true <- gen$cl_true
```

Using the Eric method, we obtain:

```r
result <- CARP_RCC(X,delta=15,zeta=0.1,phi=.1,method='cvx_uni',
                        lam.begin =0.01,lam.step=1.05,rho=1,tau=3,cl_true=cl_true,randmode='HA',max.log=200)
result$cl_est
[1] 1 1 1 2 1 1 1 1 1 1 3 1 4 4 4 4 4 4 5 4 4 4 4 4
```

Using the proposed method, we obtain: 

```r
result <- CARP_RCC(X,delta=15,zeta=0.1,phi=.1,method='Rcvx_uni',
                        lam.begin =0.01,lam.step=1.05,rho=1,tau=3,cl_true=cl_true,randmode='HA',max.log=200)
result$cl_est
[1] 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2
```

The path graph can be obtained from the `cl_matrix`, and the pictures are shown in the paper. 











