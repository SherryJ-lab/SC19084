#' @title A function to estimate coefficient vector beta  
#' @description min_beta  (2 Lam)^(-1) (z - beta)^2 + p_lambda(|beta|) for L1 penalty
#' @param XXmat X is columnwise rescaled to make each column have L2 norm root-n and XXmat=n^(-1)t(X)X
#' @param cvec n^(-1)t(X)y
#' @param inival initial value of beta
#' @param lambda ridge parameter
#' @param varset variance set
#' @param maxiter maximum iteration time
#' @param tol limitation to end the iteration
#' @return etimator of beta
#' @import Rcpp microbenchmark knitr rmarkdown boot bootstrap DAAG Ball GeneralizedHyperbolic parallel
#' @example 
#' \dontrun{
#' n<-10;p<-15
#' b0<-matrix(0,p,1)
#' b0[1:7]<-c(1,-0.5,0.7,-1.2,-0.9,0.3,0.55)
#' epsi <-0.01
#' cormat=matrix(1,p,p)
#' corval = 0.5
#' for (i in 1:(p-1)){
#'   for (j in (i+1):p){
#'     cormat[i,j]=corval^abs(i-j)
#'     cormat[j,i]=cormat[i,j]
#'  }
#' }
#' sigmat=cormat^0.5
#' X<-matrix(rnorm(n*p),n,p)%*%sigmat
#' y<-X%*%b0+0.4*matrix(rnorm(n),n,1)
#' Xsca <- sqrt(colSums(X^2))/sqrt(n)
#' X <- X/(matrix(1,n,1)%*%Xsca)
#' XXmat <- n^(-1)*t(X)%*%X
#' cvec <- n^(-1)*t(X)%*%y
#' lammax = max(abs(cvec))
#' lammin = epsi*lammax
#' lamvec = lammax - (lammax - lammin)/(50 - 1)*c(0:(50 - 1))
#' inival<-rep(0,p)
#' lambda<-lamvec[50]
#' varset<-vector(mode="numeric")
#' ol1(XXmat, cvec, inival, lambda, varset, 50, 1e-4)
#' }
#' @export
ol1<-function(XXmat, cvec, inival, lambda, varset, maxiter, tol){
  p<-length(cvec)
  
  #first round of iteration
  beta<-inival
  iter<-1
  update<-1
  ind<-c(1:p)
  varset<-union(which(inival!=0),varset)
  while((iter <= maxiter) & (update > tol)){
    iter <- iter + 1
    betaold <- beta
    for(k in 1:length(ind)){
      setr <- ind
      I <- setr[k]
      setr<- setr[-k]
      
      #solve min_beta  (2 Lam)^(-1) (z - beta)^2 + p_lambda(|beta|) for scalar beta
      z1 <- XXmat[I, I]
      Lam <- 1/z1
      if(is.null(setr)){
        z<-cvec[I]/z1
      }else{
        z <- (cvec[I] - XXmat[I, setr]%*%beta[setr])/z1
      }
      beta[I]<-sign(z)*max(0, abs(z) - Lam*lambda)
    }
    ind<-which(beta!=0)
    update <-sum(abs(beta - betaold))
    setr <- setdiff(1:p, ind)
    if(length(setr)==0){
      resc<-0
    }else{
      resc<-abs(cvec[setr] - XXmat[setr, ind]%*%beta[ind])
    }
    indm<-which(resc>lambda)
    ind <- union(c(setr[indm],ind), varset)
  }
  return(beta)
}