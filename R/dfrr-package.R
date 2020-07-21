#' @details
#' Implementing Function-on-Scalar Regression model in which the response function is dichotomized and observed sparsely. This package provides smooth estimations of functional regression coefficients and principal components for the dfrr model.
#' The main function in the dfrr-package is \link{dfrr}().
#' @keywords Dichotomized Functional, Dichotomized Functional, FDA, Functional Regression
#' @aliases dfrr-package
#' @references Fatemeh Asgari, Alamatsaz Mohammad Hossein, Hayati Saeed (2021).
#' Dichotomized Functional Response Regression Model.
#' <http://arxive.org/adress_to_paper>
#' @examples
#' set.seed(2000)
#' N<-50;M<-24
#' X<-rnorm(N,mean=0)
#' Y<-simulate.simple.dfrr(beta0=function(t){cos(pi*t+pi)},
#'                         beta1=function(t){2*t},
#'                         X=X,M=M)
#'
#' dfrr_fit<-dfrr(Y~X,yind=time)
#'
#' coefs<-coef(dfrr_fit)
#'   plot(coefs)
#'
#' fitteds<-fitted(dfrr_fit)
#'   plot(fitteds)
#'
#' resids<-residuals(dfrr_fit)
#'
#' eig<-fpca(dfrr_fit)
#' plot(eig,plot.contour=TRUE,plot.3d.surface = TRUE)
#'
#' newdata<-data.frame(X=c(1,0))
#'   preds<-predict(dfrr_fit,newdata=newdata)
#'   plot(preds,conf.level = 0.9)
#'
#' newdata<-data.frame(X=c(1,0))
#'   preds<-predict(dfrr_fit,newdata=newdata)
#'   plot(preds,conf.level = 0.9)
#'
#' newdata<-data.frame(X=c(1,0))
#' newydata<-data.frame(.obs=rep(1,5),.index=c(0.0,0.1,0.2,0.3,0.7),.value=c(1,1,1,0,0))
#' preds<-predict(madras_dfrr,newdata=newdata,newydata = newydata)
#' plot(preds,conf.level = 0.9)
#'
"_PACKAGE"

