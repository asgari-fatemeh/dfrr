#' Get the basis functions from a dfrr-object
#'
#'Returns the basis functions employed in fitting a dfrr-object.
#'
#'@inheritParams fitted.dfrr
#'
#'@examples
#' set.seed(2000)
#' \donttest{N<-50;M<-24}
#' \dontshow{N<-30;M<-12}
#' X<-rnorm(N,mean=0)
#' time<-seq(0,1,length.out=M)
#' Y<-simulate_simple_dfrr(beta0=function(t){cos(pi*t+pi)},
#'                         beta1=function(t){2*t},
#'                         X=X,time=time)
#' \donttest{dfrr_fit<-dfrr(Y~X,yind=time)}
#' \dontshow{dfrr_fit<-dfrr(Y~X,yind=time,T_E=1)}
#' coefs<-coef(dfrr_fit,return.fourier.coefs=TRUE)
#'
#' basis<-basis(dfrr_fit)
#' evaluated_coefs<-coefs%*%t(fda::eval.basis(time,basis))
#'
#' #Plotting the regression coefficients
#' par(mfrow=c(1,2))
#' plot(time,evaluated_coefs[1,],'l',main="Intercept")
#' plot(time,evaluated_coefs[2,],'l',main="X")
#'
#' @export
basis <-
function(object){
  dfrr_fit<-object
  dfrr_fit$basis
}
