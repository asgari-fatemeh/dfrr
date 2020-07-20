#'QQ-plot for dfrr residuals
#'
#'The output gives the qq-plot of estimated measurment error.
#'
#'@param residuals.dfrr a \code{residuals.dfrr}-object.
#'@param ... graphical parameters passed to \code{car::\link[car]{qqPlot}}
#'@inheritParams fitted.dfrr
#'
#'@examples
#' N<-50;M<-24
#' X<-rnorm(N,mean=0)
#' time<-seq(0,1,length.out=M)
#' Y<-simulate.simple.dfrr(beta0=function(t){cos(pi*t+pi)},
#'                         beta1=function(t){2*t},
#'                         X=X,time=time)
#' dfrr_fit<-dfrr(Y~X,yind=time)
#' resid<-residuals(dfrr_fit)
#' plot(resid)
#' #qq(dfrr_fit)
#'
#' @method plot residuals.dfrr
#'
#' @export
#'
plot.residuals.dfrr<-function(residuals.dfrr,...){
  dfrr_fit<-attr(residuals.dfrr,"dfrr_fit")
  if(attr(residuals.dfrr,"standardized"))
    residuals.dfrr<-residuals(dfrr_fit,standardized=TRUE)

  if( is.data.frame(residuals.dfrr))
  {
    resids<-residuals.dfrr$residuals
    names(resids)<-paste0(residuals.dfrr$.obs,",",residuals.dfrr$.index)
    resids<-resids[!is.na(resids)]
  }
  else{
    resids<- c(residuals.dfrr)
    idss<-rep(dfrr_fit$ids_rows,ncol(residuals.dfrr))
    timess<-rep(dfrr_fit$yind,each=nrow(residuals.dfrr))
    names(resids)<-paste0(idss,",",timess)
    resids<-resid[!is.na(resids)]
  }
  car::qqPlot(resids,...)
}

#' @export
qq<-function(x,...){UseMethod("qq")}

#' @rdname plot.residuals.dfrr
#' @export
qq.dfrr<-function(dfrr_fit,...){
  resids<-residuals(dfrr_fit)
  plot.residuals.dfrr(resids,...)

}







