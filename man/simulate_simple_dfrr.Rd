% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_simple_dfrr.R
\name{simulate_simple_dfrr}
\alias{simulate_simple_dfrr}
\title{Simulating a Simple \code{dfrr} Model}
\usage{
simulate_simple_dfrr(
  beta0 = function(t) {     cos(pi * t + pi) },
  beta1 = function(t) {     2 * t },
  X = rnorm(50),
  time = seq(0, 1, length.out = 24),
  sigma2 = 0.2
)
}
\arguments{
\item{beta0, beta1}{(optional) functional intercept and slope parameters}

\item{X}{an (optional) vector consists of scalar covariate}

\item{time}{an (optional) vector of time points for which, each sample curve is observed at.}

\item{sigma2}{variance of the measurement error in the \code{dfrr} model}
}
\description{
Simulation from a simple dfrr model:
 \deqn{Y_{i}(t)=I(\beta_0(t)+\beta_1(t)*x_{i}+\varepsilon_{i}(t)+\epsilon_{i}(t)\times\sigma^2>0),}
 where \eqn{I(.)} is the indicator function, \eqn{\varepsilon_{i}} is a Gaussian random function, and \eqn{\epsilon_{i}(t)} are iid standard normal for each \eqn{i} and \eqn{t} independent of \eqn{\varepsilon_{i}}.
 For demonstration purpose only.
}
\examples{
N<-50;M<-24
X<-rnorm(N,mean=0)
time<-seq(0,1,length.out=M)
Y<-simulate_simple_dfrr(beta0=function(t){cos(pi*t+pi)},
                        beta1=function(t){2*t},
                        X=X,time=time)
}
