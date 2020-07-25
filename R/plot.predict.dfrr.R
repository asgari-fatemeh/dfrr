#'Plot  dfrr predictions
#'
#'Plot a  \code{predict.dfrr} object.
#'
#'The output is the plot of predictions of latent functions given the new covariates.
#' For the case in which \code{newydata} is also given, the predictions are plotted
#'  over the observed binary sequence.
#'The binary sequence is illustrated with circles and filled circles for the values
#'of zero and one, respectively.
#'
#'@param x a \code{predict.dfrr}-object
#'@param id a vector of length one or more containing subject ids to plot. Must be matched with
#' \code{rownames(newdata)}. Defaults to
#'  all subject ids.
#'@param main a vector of length one or \code{length(id)} containing the title of plots.
#'@param col,lwd,lty,... graphical parameters passed to \code{\link{plot}}
#'@param cex.circle,col.circle size and color of circles and filled circles.
#'@param ylim a vector of length two indicating the range of y-axis of the plot.
#'
#'@importFrom graphics lines par points
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
#'
#' newdata<-data.frame(X=c(1,0))
#'   preds<-predict(dfrr_fit,newdata=newdata)
#'   plot(preds)
#'
#' newdata<-data.frame(X=c(1,0))
#' newydata<-data.frame(.obs=rep(1,5),.index=c(0.0,0.1,0.2,0.3,0.7),.value=c(1,1,1,0,0))
#' preds<-predict(dfrr_fit,newdata=newdata,newydata = newydata)
#' plot(preds)
#'
#'@method plot predict.dfrr
#'
#'@export
#'
#' @references Choi, H., & Reimherr, M.  A geometric approach to confidence regions and bands for functional parameters .
#'  \emph{Journal of the Royal Statistical Society, Series B Statistical methodology} 2018; 80:239-260.
#'
plot.predict.dfrr <-
  function(x,id=NULL,
           main=id,col='blue',lwd=2,lty="solid",cex.circle=1,col.circle='black',ylim=NULL,...){
    predict.dfrr<-x





    if(!is.null(main))
      if(is.na(main))
        main<-""

    dfrr_fit<-attr(predict.dfrr,"dfrr_fit")

    coefs<-dfrr_fit$pred_data$coefs
    ids<-dfrr_fit$pred_data$ids
    ydata<-dfrr_fit$pred_data$ydata
    zzt<-dfrr_fit$pred_data$zzt
    standardized<-dfrr_fit$pred_data$standardized

    time<-seq(dfrr_fit$range[1],dfrr_fit$range[2],length.out = 100)
    basis<-basis(dfrr_fit)
    E2<-t(fda::eval.basis(time,basis))

    if(standardized)
      cv<-dfrr_fit$sigma_2*diag(nrow = length(time))
    else
      cv<-dfrr_fit$sigma_2*diag(diag(t(E2)%*%dfrr_fit$sigma_theta%*%E2))


    if(is.null(id))
      id<-ids
    else
      id<-intersect(id,ids)

    if(length(id)==0)
      stop("There is no sample with the specified id(s)")


    if(length(main)!=length(id))
      main<-rep(main[1],length(id))
    
      lbl<-paste0("Prediction (Id: ",id[i],")")

    for(i in 1:length(id)){
      ind<-which(ids==id[i])

      pred<-t(E2)%*%t(t(coefs[i,]))


      mx<-max(abs(min(pred)),abs(max(pred)))

      mx<-mx*1.05



      ylim2<-c(-mx,mx)

      if(!is.null(ylim))
        ylim2<-ylim

      if(is.null(main[i]))
        plot(time,pred,'l',main=lbl,ylim=ylim2,col=col,lwd=lwd,lty=lty,...)
      else if(main[i]=="")
        plot(time,pred,'l',ylim=ylim2,col=col,lwd=lwd,lty=lty,...)
      else
        plot(time,pred,'l',main=main[i],ylim=ylim2,col=col,lwd=lwd,lty=lty,...)


      if(!is.null(ydata$M))
        if(ydata$M[ind]>0){
          points(ydata$time[[ind]],ydata$Y[[ind]]*0,pch=1,cex=cex.circle)
          tme1<-ydata$time[[ind]][ydata$Y[[ind]]==1]
          points(tme1,tme1*0,pch=16,cex=cex.circle,col=col.circle)
        }


      if(i<length(id))
        invisible(readline(prompt="Hit <Returen> to see next plot:"))
    }
  }
