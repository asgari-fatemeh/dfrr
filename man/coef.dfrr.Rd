% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coef.dfrr.R
\name{coef.dfrr}
\alias{coef.dfrr}
\title{Get estimated coefficients from a dfrr fit}
\usage{
\method{coef}{dfrr}(
  object,
  standardized = NULL,
  unstandardized = !standardized,
  return.fourier.coefs = NULL,
  return.evaluations = !return.fourier.coefs,
  time_to_evaluate = NULL,
  ...
)
}
\arguments{
\item{object}{a \code{dfrr}-object}

\item{standardized, unstandardized}{a \code{boolean} indicating whether stanadrdized/unstandardized regression coefficients are reported.
Only standardized regression coefficients are identifiable, thus the arugment is defaults to \code{standardized=TRUE}.}

\item{return.fourier.coefs, return.evaluations}{a \code{boolean} indicating whether the Fourier coefficients of regression coefficients are returned
(\code{return.fourier.coefs=TRUE}), or evaluations of the regression coefficients (\code{return.evaluations=TRUE}).
Defaults to \code{return.fourier.coefs=TRUE}.}

\item{time_to_evaluate}{a numeric vector indicating the set of time points for evaluating the functional regression  coefficients, for the case of \code{return.evaluations=TRUE}.}

\item{...}{dot argument, just for consistency with the generic function}
}
\description{
Returns estimations of the smooth functional regression coefficients \eqn{\beta(t)}.
The result is a matrix of either Fourier coefficients or evaluations. See Details.
}
\details{
This function will return either the Fourier coefficients or the evaluation of
estimated coefficients. Fourier coefficients which are reported are
based on the a set of basis which can be determined by \code{\link{basis}(dfrr_fit)}.
Thus the evaluation of regression coefficients on the set of time points specified by vector \code{time},
equals to \code{fitted(dfrr_fit)\%*\%t(\link[fda]{eval.basis}(time,\link{basis}(dfrr_fit)))}.

Consider that the unstandardized estimations are not identifiable. So, it is recommended to
extract and report the standardized estimations.
}
\examples{
set.seed(2000)
\donttest{N<-50;M<-24}
\dontshow{N<-30;M<-12}
X<-rnorm(N,mean=0)
time<-seq(0,1,length.out=M)
Y<-simulate_simple_dfrr(beta0=function(t){cos(pi*t+pi)},
                        beta1=function(t){2*t},
                        X=X,time=time)
\donttest{dfrr_fit<-dfrr(Y~X,yind=time)}
\dontshow{dfrr_fit<-dfrr(Y~X,yind=time,T_E=1)}
coefs<-coef(dfrr_fit)
plot(coefs)
}
\seealso{
\code{\link{plot.coef.dfrr}}
}
