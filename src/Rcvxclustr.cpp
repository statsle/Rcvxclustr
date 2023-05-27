#include <iostream>
#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;
using namespace arma;

//[[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
NumericVector robustweights(NumericMatrix X, double delta=15.0, double zeta=0.1) {
    int n = X.nrow(), p=X.ncol();
    mat X1(X.begin(),n,p,false);
    
    vec d(p,fill::zeros);
    mat w(n,n,fill::zeros);
    
    for (int i=0;i<n-1;i++){
        for (int j=i+1;j<n;j++){
            //double num1=0.;
            double num2=0.;
            for (int k=0;k<p;k++){
                d(k) = fabs(X1(i,k)-X1(j,k));
                num2=num2+pow(d(k),2.0);
            }
            num2 = sqrt(num2);
            if (num2 == 0){
            	w(i,j)=0.;
            	w(j,i)=0.;
			}
            else if(num2 < delta){
            	w(i,j)=exp(-(zeta)*(0.5*num2));
            	w(j,i)=exp(-(zeta)*(0.5*num2));
			} 
			else {
				w(i,j)=exp(-(zeta)*(delta)*0.5);
            	w(j,i)=exp(-(zeta)*(delta)*0.5);
			}            
        }
    }
    return wrap(w);
}


vec prox_L2(vec v, double sigma){
     int n= v.size();
     int i;
     vec prox(n,fill::zeros);
     double sum=0.0;
     for (i=0; i<n; i++){
         sum = sum + pow(v(i),2.0);
     }
     double l2norm=sqrt(sum);
     if (sum==0.0){
         for (i=0; i<n; i++)
             prox(i)=v(i);
     }
     else{
         for(i=0; i<n; i++)
             prox(i)=fmax(0.0, 1.0-(sigma/l2norm))*v(i);
     }
     return prox;
}
int Sign(double a){
     if(a>0){return 1;}
     else if (a<0){return -1;}
     return 0;
 }
double soft_thresholding(double a, double b){
    if (b<0){
        return 0.0;
    }else{
        return Sign(a)*fmax(0.0,fabs(a)-b);
     }
}

mat update_U(mat V, mat Y, mat Z, mat W, mat E){
    int n = W.n_rows;
    mat I(n,n,fill::eye);
    mat Et = trans(E);
    mat sumEtEI = Et*E+I;
    mat invEtEI = inv(sumEtEI);
    // get the final update for U
    mat U = invEtEI*(Et*(V+Y)+W+Z);
    // update U
    return U;
}
// update W
mat update_W(mat X, mat U, mat Z, double tau, double rho){
    int n = X.n_rows, p=X.n_cols;
    int i,j;
    mat W(n,p,fill::zeros);
    double w1, w2;
    for (i=0; i<n; i++){
        for (j=0; j<p; j++){
            if (fabs(rho*(X(i,j)-U(i,j)+Z(i,j))/(1+rho))<=tau){//fix bug
                w1=(X(i,j)+rho*(U(i,j)-Z(i,j)))/(1+rho);
                W(i,j)=w1;
             }
            else{
                w2=X(i,j)-soft_thresholding(X(i,j)-U(i,j)+Z(i,j),tau/rho);
                W(i,j)=w2;
            }
        }
    }
    return W;
}
// update V
mat update_V(mat U, mat Y, mat E, vec wt, double lambda, double rho){
    int p=U.n_cols;
    int nK = E.n_rows;
    int j, k;
    mat EU = E*U;
    mat V(nK,p,fill::zeros);
    vec a(p,fill::zeros);
    vec aa(p,fill::zeros);
    for (k=0;k<nK;k++){
        for (j=0;j<p;j++){
            a(j)= EU(k,j)-Y(k,j);
        }
        aa=prox_L2(a,wt(k)*lambda/rho);
        for (j=0;j<p;j++){
            V(k,j)=aa(j);
        }
    }
    return V;
}

//update Z
mat update_Z(mat Z_old,mat U, mat W, double rho){
    int n = U.n_rows, p=U.n_cols;
    int i,j;
    double Z_temp;
    mat Z(n,p,fill::zeros);
    for (i=0; i<n; i++){
        for (j=0;j<p;j++){
            Z_temp = Z_old(i,j)-rho*(U(i,j)-W(i,j));
            Z(i,j) = Z_temp;
        }
    }
    return Z;
}
// update Y
mat update_Y(mat Y_old, mat U, mat V, mat E,double rho){
    int p=U.n_cols, nK = E.n_rows;
    int i, j;
    double y;
    mat EU = E*U;
    mat Y(nK,p,fill::zeros);
    for (i=0; i<nK; i++){
        for (j=0; j<p; j++){
            y=Y_old(i,j);
            y=y-rho*(EU(i,j)-V(i,j));
            Y(i,j)=y;
        }
    }
    return Y;
}
double tolerance(mat W, mat W_old){
    int n = W.n_rows, p=W.n_cols;
    int i,j;
    double temp;
    double sum=0.;
    for (i=0;i<n;i++){
        for (j=0;j<p;j++){
            temp = fabs(W(i,j)-W_old(i,j));
            sum = sum + pow(temp,2.0);
        }
    }
    return sqrt(sum);
 }

//[[Rcpp::export]]
List robust_convex_cluster(NumericMatrix X, NumericMatrix W, NumericMatrix V,NumericMatrix Y,
                           NumericMatrix Z, NumericMatrix E, int max_iter, double tol_abs,
                           double lambda, double rho, double tau, NumericVector wt){
    int n = X.nrow(), p=X.ncol();
    int it, iter;
    double tol;
    mat X1 = as<arma::mat>(X);
    mat W1 = as<arma::mat>(W);
    mat V1 = as<arma::mat>(V);
    mat Y1 = as<arma::mat>(Y);
    mat Z1 = as<arma::mat>(Z);
    mat E1 = as<arma::mat>(E);
    vec wt1 = as<arma::vec>(wt);
    mat U1(n,p,fill::zeros);
    
    for (it=0; it<max_iter; it++){
        mat W_old=W1;
        mat Y_old=Y1;
        mat Z_old=Z1;
        U1=update_U(V1,Y1,Z1,W1,E1);
        W1=update_W(X1,U1,Z1,tau,rho);
        V1=update_V(U1,Y1,E1,wt1,lambda,rho);
        Y1=update_Y(Y_old,U1,V1,E1,rho);
        Z1=update_Z(Z_old,U1,W1,rho);
        tol = tolerance(W1,W_old);
        if (tol<tol_abs){
            break;
        }
    }
    if (it > max_iter){
    iter = max_iter;
    }else{
    iter = it;
    }
List a= List::create(_["U"]=wrap(U1),
                     _["W"]=wrap(W1),
                     _["V"]=wrap(V1),
                     _["Y"]=wrap(Y1),
                     _["Z"]=wrap(Z1),
                     _["iter"]=it,
                     _["tol"]=tol);
return a;
}
