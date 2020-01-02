#' @title uhard
#' @description min_beta  (2 Lam)^(-1) (z - beta)^2 + lambda0 |beta| + p_lambda(|beta|) for hard-thresholding penalty
#' @param z wavelet transform coefficient
#' @param Lam control coefficieent
#' @param Lamda0 for initial value of z
#' @param Lamda threshold value
#' @return beta for hard-thresholding penalty
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
#' z1 = XXmat[1, 1]
#' Lam<-1/z1
#' cvec <- n^(-1)*t(X)%*%y
#' z<-cvec[1]/z1
#' lammax = max(abs(cvec))
#' lammin = epsi*lammax
#' lamvec = lammax - (lammax - lammin)/(50 - 1)*c(0:(50 - 1))
#' lambda0<-lamvec[49]
#' lambda<-lamvec[50]
#' uhard(z, Lam, lambda0, lambda)
#' }
#' @export
uhard<-function (z, Lam, lambda0, lambda){
  z <- sign(z)*max(0, abs(z) - Lam*lambda0)
  if(Lam == 1){
    beta <- z*(abs(z) > lambda)
  }else{
    z0 <- sign(z)*max(0, abs(z) - Lam*lambda)/(1 - Lam)
  if(abs(z) <= lambda){
    beta <- z0
  }else if((abs(z) <= Lam*lambda) && ((z - z0)^2/2 + Lam*((lamda^2 - (max(0, lamda - abs(z0)))^2)/2) <= Lam*((lamda^2 - (max(0, lamda - abs(z)))^2)/2))){
    beta = z0
  }else{
    beta = z
  }
  }
  return(beta)
}
