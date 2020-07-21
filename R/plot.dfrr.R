#' Plot a dfrr fit
#'
#' Plot the regression coefficients, principal components, kernel function and residuals of a \code{dfrr}-object.
#'
#'@details
#' The contour plot of the kernel function is produced if the package \code{\link[ggplot2]{ggplot2}} is installed.
#' Plotting the 3d surface  of the kernel function is also depends on the package \code{\link[plotly]{plotly}}.
#'To produce the qq-plot, the package \code{\link[car]{car}} must be installed.
#'
#'@inheritParams fpca
#'
#'@examples
#'set.seed(2000)
#' N<-50;M<-24
#' X<-rnorm(N,mean=0)
#' time<-seq(0,1,length.out=M)
#' Y<-simulate.simple.dfrr(beta0=function(t){cos(pi*t+pi)},
#'                         beta1=function(t){2*t},
#'                         X=X,time=time)
#' dfrr_fit<-dfrr(Y~X,yind=time)
#' plot(dfrr_fit)
#'
#'@method plot dfrr
#'
#'@export

plot.dfrr <-
function(dfrr_fit){
 #Plotting regression coefficients
coefs<-coef(dfrr_fit)
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
  plot(coefs,select=ind1:ind2,ask.hit.return=FALSE)
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
  plot(pcs,select=ind1:ind2,ask.hit.return=FALSE)
  invisible(readline(prompt="Hit <Returen> to see next plot:"))

}
  #plotting contour of kernel function
plot(pcs,plot.contour=TRUE,plot.eigen.functions=FALSE)
invisible(readline(prompt="Hit <Returen> to see next plot:"))

  #Plotting 3d surface of kernel function
plot(plot.3dsurface=TRUE,plot.eigen.functions=FALSE)
invisible(readline(prompt="Hit <Returen> to see next plot:"))

  #Plotting residual functions
resids<-residuals(dfrr_fit)
plot(resids)

}
