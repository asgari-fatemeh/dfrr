#'QQ-plot for dfrr residuals
#'
#'The output gives the qq-plot of estimated measurment error.
#'
#'@param residuals.dfrr a \code{residuals.dfrr}-object.
#'@param ... graphical parameters passed to \code{car::\link[car]{qqPlot}}
#'@inheritParams fitted.dfrr
#'
#'@examples
#'\donttest{
#' \donttest{N<-50;M<-24}
#' \dontshow{N<-30;M<-12}
#' X<-rnorm(N,mean=0)
#' time<-seq(0,1,length.out=M)
#' Y<-simulate_simple_dfrr(beta0=function(t){cos(pi*t+pi)},
#'                         beta1=function(t){2*t},
#'                         X=X,time=time)
#' \donttest{dfrr_fit<-dfrr(Y~X,yind=time)}
#' \dontshow{dfrr_fit<-dfrr(Y~X,yind=time,T_E=1)}
#' resid<-residuals(dfrr_fit)
#' \donttest{plot(resid)}
#' #qq(dfrr_fit)
#'}
#' @method plot residuals.dfrr
#'
#' @export
#'
plot.residuals.dfrr<-function(residuals.dfrr,...){

  dfrr_fit<-attr(residuals.dfrr,"dfrr_fit")
  standardized<-attr(residuals.dfrr,"standardized")
  print(standardized)
  if(!standardized)
    residuals.dfrr<-residuals.dfrr(dfrr_fit,standardized=TRUE)

  if(is.data.frame(residuals.dfrr))
  {
    resids<-residuals.dfrr$residual
    names(resids)<-paste0(residuals.dfrr$.obs,",",residuals.dfrr$.index)
    resids<-resids[!is.na(resids)]
  }
  else{
    resids<- c(residuals.dfrr)
    idss<-rep(dfrr_fit$ids_rows,ncol(residuals.dfrr))
    timess<-rep(dfrr_fit$yind,each=nrow(residuals.dfrr))
    names(resids)<-paste0(idss,",",timess)
    resids<-resids[!is.na(resids)]
  }
print(4)
  if(requireNamespace("car",quietly = TRUE))
    car::qqPlot(resids,...)
  else
    warning("Package 'car' is not installed.")
}

#'qq-plot Generic function
#'
#'This is a generic function for qq() method.
#'
#'@param x an object
#'@param ... extra parameters passed to S3 methods
#'
#' @export
 qq<-function(x,...){UseMethod("qq")}

#' @rdname plot.residuals.dfrr
#' @method qq dfrr
#' @export
qq.dfrr<-function(dfrr_fit,...){
  resids<-residuals.dfrr(dfrr_fit)
  plot.residuals.dfrr(resids,...)

}








