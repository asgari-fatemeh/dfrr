#'Plot dfrr fitted latent functions
#'
#'Plot a  \code{fitted.dfrr} object.
#'
#'The output is the plot of latent curves over the observed binary sequence.
#'The binary sequence is illustrated with circles and filled circles for the values
#'of zero and one, respectively.
#'
#'
#'@param fitted.dfrr the output of the function \link{fitted.dfrr}
#'@param id a vector of length one or more containing subject ids to plot. Must be matched with
#' \code{rownames(<response>)} or the \code{.obs} column of \code{ydata}. Defaults to
#'  all  subject ids.
#'@param main a vector of length one or \code{length(id)} containing the title of
#' plots.
#'@param col,lwd,lty,... graphical parameters passed to \code{\link{plot}}
#'@param cex.circle,col.circle size and color of circles and filled circles.
#'@param ylim a vector of length two indicating the range of y-axis of the plot.
#'
#'@importFrom graphics lines par points
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
#' fitteds<-fitted(dfrr_fit)
#' plot(fitteds)
#'
#'@method plot fitted.dfrr
#'@export
plot.fitted.dfrr <-
function(fitted.dfrr,id=NULL,main=NULL,
                           col='blue',lwd=2,lty="solid",cex.circle=1,col.circle='black',ylim=NULL,...)
{
  if(!is.null(main))
    if(is.na(main[1]))
      main<-""

  dfrr_fit<-attr(fitted.dfrr,"dfrr_fit")
  standardized<-attr(fitted.dfrr,"standardized")


  if(is.null(id)){
    id<-dfrr_fit$ids
  }else{
    id<-intersect(id,dfrr_fit$ids)
  }

  if(length(id)==0)
    stop("There is no sample with the specified id(s)")

  if(is.null(main))
    main<-paste0("Id: ",id)

  basis<-basis(dfrr_fit)
  time<-seq(dfrr_fit$range[1],dfrr_fit$range[2],length.out=100)
  E2<-t(fda::eval.basis(time,basis))

  if(length(main)!=length(id))
    main<-rep(main[1],length(id))

  if(standardized){
    coefs<-dfrr_fit$fitted_coefs_std
  }else{
    coefs<-dfrr_fit$fitted_coefs
  }


  for(i in 1:length(id)){
    ind<-which(dfrr_fit$ids==id[i])

    fitted.dfrr_value<-t(E2)%*%t(t(coefs[ind,]))

    lbl<-paste0("Id: ",id[i])
    mx<-max(abs(min(fitted.dfrr_value)),abs(max(fitted.dfrr_value)))*1.05
    ylim2<-c(-mx,mx)
    if(!is.null(ylim))
      ylim2<-ylim

    if(is.null(main[i]))
      plot(time,fitted.dfrr_value,'l',main=lbl,ylim=ylim2,col=col,lwd=lwd,lty=lty,...)
    else if(main[i]=="")
      plot(time,fitted.dfrr_value,'l',ylim=ylim2,col=col,lwd=lwd,lty=lty,...)
    else
      plot(time,fitted.dfrr_value,'l',main=main[i],ylim=ylim2,col=col,lwd=lwd,lty=lty,...)

    points(dfrr_fit$data$time[[ind]],dfrr_fit$data$Y[[ind]]*0,pch=1,cex=cex.circle)
    tme1<-dfrr_fit$data$time[[ind]][dfrr_fit$data$Y[[ind]]==1]
    points(tme1,tme1*0,pch=16,cex=cex.circle,col=col.circle)

    if(i<length(id))
        invisible(readline(prompt="Hit <Returen> to see next plot:"))
  }

}
