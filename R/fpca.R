#'Functional principal component analysis of a dfrr fit
#'
#'\code{fpca()}  returns estimations of the smooth principal components/eigen-functions
#' and the corresponding eigen-values of the residual function in the \code{dfrr} model.
#' The result is a named list containing  the vector of eigen-values and the matrix of Fourier coefficients. See Details.
#'
#'@details  Fourier coefficients which are reported are
#' based on the a set of basis which can be determined by \code{\link{basis}(dfrr_fit)}.
#' Thus the evaluation of pricipal component/eigen-function on the set of time points specified by vector \code{time},
#' equals to \code{fpca(dfrr_fit)\%*\%t(\link[fda]{eval.basis}(time,\link{basis}(dfrr_fit)))}.
#'
#' Consider that the unstandardized estimations are not identifiable. So, it is recommended to
#' extract and report the standardized estimations.
#'
#'
#'
#'@return
#' \code{fpca(dfrr_fit)} returns a list containtng the following components:
#' \item{values}{a vector containing the eigen-values of the standaridized/unstandardized covariance operator of
#' the residual function term in \code{dfrr} model,
#' sorted in decreasing order.}
#' \item{vectors}{a matrix whose columns contain the Fourier coefficients of the
#'  principal components/eigen-functions of the standaridized/unstandardized covariance operator of
#' the residual function term in \code{dfrr} model,
#' sorted based on the corresponding eigen-values.}
#'
#'
#'@inheritParams fitted.dfrr
#'@param standardized,unstandardized a \code{boolean} indicating whether stanadrdized/unstandardized pricipal components/eigen-functions are reported.
#' Only standardized pricipal components/eigen-functions are identifiable, thus the arugment is defaults to \code{standardized=TRUE}.
#'
#'
#'@seealso \code{\link{plot.fpca.dfrr}}
#'
#'@examples
#' set.seed(2000)
#' N<-50;M<-24
#' X<-rnorm(N,mean=0)
#' time<-seq(0,1,length.out=M)
#' Y<-simulate.simple.dfrr(beta0=function(t){cos(pi*t+pi)},
#'                         beta1=function(t){2*t},
#'                         X=X,time=time)
#' dfrr_fit<-dfrr(Y~X,yind=time)
#' fpcs<-fpca(dfrr_fit)
#' plot(fpcs)
#'
#'@export
fpca <-
function(dfrr_fit,standardized=NULL,unstandardized=!standardized){
  standardized<-paired.args.check(standardized,
                                  ifelse(missing(unstandardized),NA,unstandardized),
                                  "Please specify 'standardized' or 'unstandardizedd' coefficients must be reported",
                                  TRUE)
  if(standardized)
    res<-list(values=dfrr_fit$nus_std,vectors=t(dfrr_fit$Theta_std))
  else
    res<-list(values=dfrr_fit$nus,vectors=t(dfrr_fit$Theta))

  class(res)<-"fpca.dfrr"
  attr(res,"standardized")<-standardized
  attr(res,"dfrr_fit")<-dfrr_fit

  res
}
