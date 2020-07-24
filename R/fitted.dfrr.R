#' Obtain fitted curves for a dfrr model
#'
#' Fitted curves refer to the estimations of latent functional response curves.
#' The results can be either the Fourier coefficients or evaluation of the
#' fitted functions. See Details.
#'
#' @details  This function will return either the Fourier coefficients or the evaluation of
#' fitted curves to the binary sequences. Fourier coefficients which are reported are
#' based on the a set of basis which can be determined by \code{\link{basis}(dfrr_fit)}.
#' Thus the evaluation of fitted latent curves on the set of time points specified by vector \code{time},
#' equals to \code{fitted(dfrr_fit)\%*\%t(\link[fda]{eval.basis}(time,\link{basis}(dfrr_fit)))}.
#'
#' Consider that the unstandardized estimations are not identifiable. So, it is recommended to
#' extract and report the standardized estimations.
#'
#' @param object a fitted \code{dfrr}-object obtained from invoking the function \code{\link{dfrr}}.
#' @param return.fourier.coefs,return.evaluations a \code{boolean} indicating whether the Fourier coefficients of the fitted curves are returned
#'              (\code{return.fourier.coefs=TRUE}), or evaluations of the fitted curves (\code{return.evaluations=TRUE}).
#'              Defaults to \code{return.fourier.coefs=TRUE}.
#' @param time_to_evaluate a numeric vector indicating the set of time points for evaluating the fitted latent functions, for the case of \code{return.evaluations=TRUE}.
#' @param standardized,unstandardized a \code{boolean} indicating whether stanadrdized/unstandardized fitted latent curves is reported.
#' Only standardized fitted curves are identifiable, thus the arugment is defaults to \code{standardized=TRUE}.
#' @param ... dot argument, just for consistency with the generic function

#'
#' @seealso \code{\link{plot.fitted.dfrr}}
#'
#' @examples
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
#' fitteds<-fitted(dfrr_fit)
#' plot(fitteds)
#'
#' @export
fitted.dfrr <-
function(object,return.fourier.coefs=NULL,return.evaluations=!return.fourier.coefs,
         time_to_evaluate=NULL,standardized=NULL,unstandardized=!standardized,...){
  dfrr_fit<-object

  standardized<-paired.args.check(standardized,
                                  ifelse(missing(unstandardized),NA,unstandardized),
                                  "Please specify 'standardized' or 'unstandardizedd' coefficients must be reported",
                                  TRUE)

  return.fourier.coefs<-paired.args.check(return.fourier.coefs,
                                          ifelse(missing(return.evaluations),NA,return.evaluations),
                                          "Please specify only one of the 'return.fourier.coefs' or 'return.evaluations'",
                                          TRUE)

  if(standardized)
    fitted<-dfrr_fit$fitted_coefs_std
  else
    fitted<-dfrr_fit$fitted_coefs

  if(return.fourier.coefs){
    class(fitted)<-"fitted.dfrr"
    attr(fitted,"dfrr_fit")<-dfrr_fit
    attr(fitted,"standardized")<-standardized
    return(fitted)
  }


  if(is.null(time_to_evaluate))
    time_to_evaluate<-seq(dfrr_fit$range[1],dfrr_fit$range[2],length.out=100)

  E<-t(fda::eval.basis(time_to_evaluate,dfrr_fit$basis))
  fitted<-fitted%*%E
  if(!is.null(dfrr_fit$ids))
    rownames(fitted)<-dfrr_fit$ids

  fitted

  class(fitted)<-"fitted.dfrr"

  attr(fitted,"dfrr_fit")<-dfrr_fit
  attr(fitted,"standardized")<-standardized

  fitted
}
