#'Get estimated coefficients from a dfrr fit
#'
#'Returns estimations of the smooth functional regression coefficients \eqn{\beta(t)}.
#' The result is a matrix of either Fourier coefficients or evaluations. See Details.
#'
#'@details This function will return either the Fourier coefficients or the evaluation of
#' estimated coefficients. Fourier coefficients which are reported are
#' based on the a set of basis which can be determined by \code{\link{basis}(dfrr_fit)}.
#' Thus the evaluation of regression coefficients on the set of time points specified by vector \code{time},
#' equals to \code{fitted(dfrr_fit)\%*\%t(\link[fda]{eval.basis}(time,\link{basis}(dfrr_fit)))}.
#'
#' Consider that the unstandardized estimations are not identifiable. So, it is recommended to
#' extract and report the standardized estimations.
#'
#'@inheritParams summary.dfrr
#'@param standardized,unstandardized a \code{boolean} indicating whether stanadrdized/unstandardized regression coefficients are reported.
#' Only standardized regression coefficients are identifiable, thus the arugment is defaults to \code{standardized=TRUE}.
#'@param return.fourier.coefs,return.evaluations a \code{boolean} indicating whether the Fourier coefficients of regression coefficients are returned
#'              (\code{return.fourier.coefs=TRUE}), or evaluations of the regression coefficients (\code{return.evaluations=TRUE}).
#'              Defaults to \code{return.fourier.coefs=TRUE}.
#'@param time_to_evaluate a numeric vector indicating the set of time points for evaluating the functional regression  coefficients, for the case of \code{return.evaluations=TRUE}.

#'
#'@seealso \code{\link{plot.coef.dfrr}}
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
#' coefs<-coef(dfrr_fit)
#' plot(coefs)
#'@method coef dfrr
#'@export




coef.dfrr <-
function(object,standardized=NULL,unstandardized=!standardized,
                    return.fourier.coefs=NULL,
                    return.evaluations=!return.fourier.coefs,
                    time_to_evaluate=NULL,...){
  dfrr_fit<-object

  return.principal.components<-FALSE

  standardized<-paired.args.check(standardized,
                              ifelse(missing(unstandardized),NA,unstandardized),
        "Please specify 'standardized' or 'unstandardizedd' coefficients must be reported",
        TRUE)
  # return.principal.components<-paired.args.check(return.principal.components,
  #                                 ifelse(missing(return.regression.coefficients),NA,return.regression.coefficients),
  #                                 "Please specify only on of the 'return.regression.coefficients' or 'return.principal.components'",
  #                                 FALSE)
  return.fourier.coefs<-paired.args.check(return.fourier.coefs,
                                  ifelse(missing(return.evaluations),NA,return.evaluations),
                                  "Please specify only on of the 'return.fourier.coefs' or 'return.evaluations'",
                                  TRUE)


  if(return.fourier.coefs)
    if(return.principal.components){
      if(standardized)
        coefs<-dfrr_fit$Theta_std
      else
        coefs<-dfrr_fit$Theta
    }else{
      if(standardized)
        coefs<-dfrr_fit$B_std
      else
        coefs<-dfrr_fit$B
    }

  if(is.null(time_to_evaluate))
    time_to_evaluate<-seq(dfrr_fit$range[1],dfrr_fit$range[2],length.out=100)

  E<-t(fda::eval.basis(time_to_evaluate,dfrr_fit$basis))
  if(!return.fourier.coefs)
  if(return.principal.components){
    if(standardized)
      coefs<-dfrr_fit$Theta_std%*%E
    else
      coefs<-dfrr_fit$Theta%*%E
  }else{
    if(standardized)
      coefs<-dfrr_fit$B_std%*%E
    else
      coefs<-dfrr_fit$B%*%E
  }

  if(!return.principal.components)
    rownames(coefs)<-dfrr_fit$varnames

  class(coefs)<-c("coef.dfrr",class(coefs))
  attr(coefs,"dfrr_fit")<-dfrr_fit
  attr(coefs,"standardized")<-standardized
  attr(coefs,"pc.coefs")<-return.principal.components

  coefs
}
