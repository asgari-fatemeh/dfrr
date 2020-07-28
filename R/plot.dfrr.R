#' Plot a dfrr fit
#'
#' Plot the regression coefficients, principal components, kernel function and residuals of a \code{dfrr}-object.
#'
#'@details
#' The contour plot of the kernel function is produced if the package \code{ggplot2} is installed.
#' Plotting the 3d surface  of the kernel function is also depends on the package \code{plotly}.
#'To produce the qq-plot, the package \code{car} must be installed.
#'
#'@inheritParams plot.fitted.dfrr
#'@param ... graphical parameters passed to \code{\link{plot.coef.dfrr}}
#'@param plot.kernel a boolean indicating whether plots the kernel function or not.
#'\code{ggplot2}-package and \code{plotly}-package is required to plot contour and 3d surface of kernel function.
#'@examples
#'set.seed(2000)
#' \donttest{N<-50;M<-24}
#' \dontshow{N<-30;M<-12}
#' X<-rnorm(N,mean=0)
#' time<-seq(0,1,length.out=M)
#' Y<-simulate_simple_dfrr(beta0=function(t){cos(pi*t+pi)},
#'                         beta1=function(t){2*t},
#'                         X=X,time=time)
#' \donttest{dfrr_fit<-dfrr(Y~X,yind=time)}
#' \dontshow{dfrr_fit<-dfrr(Y~X,yind=time,T_E=1)}
#' \donttest{plot(dfrr_fit)}
#' \dontshow{plot(dfrr_fit,plot.kernel=FALSE)}
#'
#'@method plot dfrr
#'@importFrom graphics plot.new
#'@export

plot.dfrr <-
function(x,plot.kernel=TRUE,...){
  dfrr_fit<-x
 #Plotting regression coefficients
coefs<-coef.dfrr(dfrr_fit)
p<-nrow(coefs)
mfrow<-c(1,1)
if(p==1)
  mfrow<-c(1,1)
else if(p==2)
  mfrow<-c(1,2)
else if(p %in% c(3,4,5,6,7,8))
  mfrow<-c(2,2)
else
  mfrow<-c(3,3)

ppg<-mfrow[1]*mfrow[2]
nographs<-ceiling(p/ppg)
par(mfrow=mfrow)
for (i in 1:nographs) {
  ind1<-(i-1)*ppg+1
  ind2<-min(p,ppg*i)
  plot.coef.dfrr(coefs,select=ind1:ind2,ask.hit.return=FALSE,...)
  if(ind2<(ppg*i))
    for(ii in (ind2+1):(ppg*i))
      plot.new()
  invisible(readline(prompt="Hit <Returen> to see next plot:"))

}

 #Plotting principal components
pcs<-fpca(dfrr_fit)
p<-ncol(pcs$vectors)
percents<-pcs$values/sum(pcs$values)
cper<-cumsum(percents)
p<-which(cper>=0.99)[1]
mfrow<-c(1,1)
if(p==1)
  mfrow<-c(1,1)
else if(p==2)
  mfrow<-c(1,2)
else if(p %in% c(3,4,5,6,7,8))
  mfrow<-c(2,2)
else
  mfrow<-c(3,3)

ppg<-mfrow[1]*mfrow[2]
nographs<-ceiling(p/ppg)
par(mfrow=mfrow)
for (i in 1:nographs) {
  ind1<-(i-1)*ppg+1
  ind2<-min(p,ppg*i)
  plot.fpca.dfrr(pcs,select=ind1:ind2,ask.hit.return=FALSE)
  if(ind2<(ppg*i))
    for(ii in (ind2+1):(ppg*i))
      plot.new()
  invisible(readline(prompt="Hit <Returen> to see next plot:"))

}

dots<-list(...)
if(plot.kernel){
  #plotting contour of kernel function
  plot.fpca.dfrr(pcs,plot.contour=TRUE,plot.eigen.functions=FALSE)
  invisible(readline(prompt="Hit <Returen> to see next plot:"))


  #Plotting 3d surface of kernel function
  plot.fpca.dfrr(pcs,plot.3dsurface=TRUE,plot.eigen.functions=FALSE)
  invisible(readline(prompt="Hit <Returen> to see next plot:"))
}


  #Plotting residual functions
resids<-residuals.dfrr(dfrr_fit)
par(mfrow=c(1,1))
plot.residuals.dfrr(resids)

}
