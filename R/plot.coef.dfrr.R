#'Plot dfrr coefficients
#'
#'Plot a  \code{coef.dfrr} object. The output is the plot of regression coefficients.
#'
#'@param x a \code{coef.dfrr}-object.
#'@param select a vector of length one or more of indices of regression
#' coefficients to plot.
#'@param ... graphical parameters passed to \code{plot}.
#'@param ask.hit.return a boolean indicating whether to wait for interaction of the user between any two plots.
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
#'
#'@method plot coef.dfrr
#'
#'@export


plot.coef.dfrr<-function(x,select=NULL,ask.hit.return=TRUE,...){
  coefs<-x

  attr(coefs,"dfrr_fit")->dfrr_fit
  attr(coefs,"standardized")->standardized
  attr(coefs,"pc.coefs")->return.principal.components

  time<-seq(dfrr_fit$range[1],dfrr_fit$range[2],length.out = 100)
  basis<-basis(dfrr_fit)
  E<-t(fda::eval.basis(time,basis))

  if(return.principal.components){
    if(standardized)
      yval<-dfrr_fit$Theta_std%*%E
    else
      yval<-dfrr_fit$Theta%*%E
  }else{
    if(standardized)
      yval<-dfrr_fit$B_std%*%E
    else
      yval<-dfrr_fit$B%*%E
  }

  if(is.null(select))
    select<-1:nrow(yval)

  if(standardized)
    nus<-dfrr_fit$nus_std/sum(dfrr_fit$nus_std)
  else
    nus<-dfrr_fit$nus/sum(dfrr_fit$nus)

  plotnames<-if(return.principal.components)
                paste0("PC ",select)
              else
                paste0("Reg.Coef. ",dfrr_fit$varnames[select])

  for(i in 1:length(select)){
    if(return.principal.components){
      variance_explained<-round(nus[select[i]]*100,digits = 2)
      if(variance_explained<=0)
        return()
    }
    if(ask.hit.return)
      invisible(readline(prompt="Hit <Returen> to see next plot:"))

    lbl<-plotnames[select[i]]
    if(standardized)
      lbl<-paste0("Standardized ",lbl)

    if(return.principal.components){
      variance_explained<-round(nus[select[i]]*100,digits = 2)
      if(variance_explained<=0)
        return()

      lbl<-paste0(lbl," (",round(nus[select[i]]*100,digits = 2),"%)")
    }

    value<-yval[select[i],]
    plot(time,value,'l',main=lbl,...)
  }


}
