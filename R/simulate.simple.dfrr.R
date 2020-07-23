#' Simulating a Simple \code{dfrr} Model
#'
#' Simulation from a simple dfrr model:
#'  \deqn{Y_{i}(t)=I(\beta_0(t)+\beta_1(t)*x_{i}+\varepsilon_{i}(t)+\epsilon_{i}(t)\times\sigma^2>0),}
#'  where \eqn{I(.)} is the indicator function, and \eqn{\epsilon_{i}(t)} is iid standard normal for each \eqn{i} and \eqn{t}.
#'  For demonstration purpose only.
#'
#'@param beta0,beta1 (optional) functional intercept and regression coefficients
#'@param X an (optional) vector consists of scalar covariate
#'@param time an (optional) vector of time point for which, each sample curve is observed at.
#'@param sigma2 variance of the measurement error in the \code{dfrr} model.
#' @examples
#' N<-50;M<-24
#' X<-rnorm(N,mean=0)
#' time<-seq(0,1,length.out=M)
#' Y<-simulate.simple.dfrr(beta0=function(t){cos(pi*t+pi)},
#'                         beta1=function(t){2*t},
#'                         X=X,time=time)
#'
#' @export simulate.simple.dfrr
simulate.simple.dfrr <-
function(beta0=function(t){cos(pi*t+pi)},
            beta1=function(t){2*t},X=rnorm(50),time=seq(0,1,length.out=24),sigma2=0.2){
  N<-length(X)
  c_f<-function(t,s){
    sig2<-1;rho<-0.5
    d<-abs(outer(t,s,"-"))
    tmp2<-sig2*(1+sqrt(3)*d/rho)*exp(-sqrt(3)*d/rho)
  }
  sigma<-c_f(time,time)

  mu<-beta0(time)
  beta<-beta1(time)
  Z<-MASS::mvrnorm(N,mu,sigma)+X%*%t(beta)+MASS::mvrnorm(N,mu*0,sigma2*diag(nrow=length(mu)))
  Y<-matrix(Z>0,nrow=N)
  Y*1
}
