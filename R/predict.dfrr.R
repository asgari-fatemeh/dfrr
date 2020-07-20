#'Prediction for dichotomized function-on-scalar regression
#'
#'Takes a \code{dfrr}-object created by \code{\link{dfrr}()} and returns predictions
#' given a new set of values for a model covariates and an optional \code{ydata}-like
#' \code{data.frame} of observations for the dichotomized response.
#'
#' @details This function will return either the Fourier coefficients or the evaluation of
#' predictions. Fourier coefficients which are reported are
#' based on the a set of basis which can be determined by \code{\link{basis}(dfrr_fit)}.
#' Thus the evaluation of predictions on the set of time points specified by vector \code{time},
#' equals to \code{fitted(dfrr_fit,return.fourier.coefs=T)\%*\%t(\link[fda]{eval.basis}(time,\link{basis}(dfrr_fit)))}.
#'
#'
#'@inheritParams fitted.dfrr
#'@param newdata a \code{data.frame} containing the values of all of the
#' model covariates at which the latent functional response is going  to be
#'  predicted.
#'@param newydata (optional) a \code{ydata}-like \code{data.frame} containing
#' the values of dichotomized response sparsly observed in the domain of function.
#'@param return.fourier.coefs,return.evaluations a \code{boolean} indicating whether the Fourier coefficients of predictions are returned
#'              (\code{return.fourier.coefs=TRUE}), or evaluations of the predictions (\code{return.evaluations=TRUE}).
#'              Defaults to \code{return.evaluations=TRUE}.
#'@param standardized,unstandardized a \code{boolean} indicating whether stanadrdized/unstandardized predictions are reported.
#' Defaults to \code{standardized=TRUE}.
#' @param time_to_evaluate a numeric vector indicating the set of time points for evaluating the predictions, for the case of \code{return.evaluations=TRUE}.
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
#' preds<-predict(dfrr_fit)
#' plot(preds)
#'
#'@seealso \code{\link{plot.predict.dfrr}}
#'
#'@export
predict.dfrr <-
function(dfrr_fit,newdata,newydata=NULL,standardized=NULL,unstandardized=!standardized,
         return.fourier.coefs=NULL,
         return.evaluations=!return.fourier.coefs,
                       time_to_evaluate=NULL){
  standardized<-paired.args.check(standardized,
                                  ifelse(missing(unstandardized),NA,unstandardized),
                                  "Please specify 'standardized' or 'unstandardizedd' coefficients must be reported",
                                  TRUE)
  return.fourier.coefs<-paired.args.check(return.fourier.coefs,
                                          ifelse(missing(return.evaluations),NA,return.evaluations),
                                          "Please specify only on of the 'return.fourier.coefs' or 'return.evaluations'",
                                          FALSE)
  if(!is.null(newydata))
    if(!all(c(".obs",".index",".value") %in% colnames(newydata)))
      stop("newydata is not of the expected structure. See the help for more details")

  if(nrow(newdata)==0)
    stop("newdata mus be a nonempty data.frame or matrix")

  newdata_<-newdata
  ncols<-ncol(newdata)
  na_inds<-sapply(1:nrow(newdata),function(i1){any(is.na(newdata[i1,]))})
  newdata<-newdata[which(na_inds==FALSE),]

  if(is.null(dim(newdata))){
    if(ncols==1)
    {
      newdata<-newdata_
      newdata[,1]<-newdata[which(na_inds==FALSE),1]
    }
  }

  if(nrow(newdata)==0)
    stop("newdata mus be a nonempty data.frame or matrix")

  if(is.null(rownames(newdata)))
    ids<-1:nrow(newdata)
  else
    ids<-rownames(newdata)

  formula2<-attr(dfrr_fit,"formula")
  xData<-model.matrix(formula2,data=newdata)

  N<-nrow(newdata)
  J<-dfrr_fit$basis$nbasis

  X<-lapply(1:N, function(i){kronecker(t(xData[i,]),diag(nrow = J))})

  basis<-basis.dfrr(dfrr_fit)

  if(standrdized)
    b<-t(t(c(t(dfrr_fit$B_std))))
  else
    b<-t(t(c(t(dfrr_fit$B))))


  Coefs<-matrix(0,N,J)
  zzt<-list()

  for(i in 1:N)
    Coefs[i,]<-X[[i]]%*%b


  if(standrdized)
    zzt[[i]]<-dfrr_fit$sigma_theta_std
  else
    zzt[[i]]<-dfrr_fit$sigma_theta



  if(!is.null(ids))
    rownames(Coefs)<-ids

  if(is.null(rownames(newdata))){
    ids<-1:N
  }else{
    ids<-rownames(newdata)
  }

  if(!is.null(newydata))
    if(length(interaction(ids,unique(newydata$.obs)))==0)
      stop("newydata .obs column does not match with the rownames of newdata")

  Ys<-list()
  times<-list()
  Ms<-c()

  if(!is.null(newydata))
    for(i in 1:N){
      Ms[i]<-0
      ind<-which(newydata$.obs==ids[i])
      if(length(ind)==0)
        next

      ys<-newydata$.value[ind]
      time<-newydata$.index[ind]

      ind<-!is.na(ys) & !is.na(time)

      time<-time[ind]
      ys<-ys[ind]

      if(length(ys)==0)
        next

      M<-length(time)
      T_G<-500


      Ys[[i]]<-ys
      times[[i]]<-time
      Ms[i]<-M

      if(standrdized)
        sigma0<-dfrr_fit$sigma_theta_std
      else
        sigma0<-dfrr_fit$sigma_theta

      Ei<-t(fda::eval.basis(time,basis))

      if(M==1)
        kttt<-t(Ei)%*%sigma0%*%Ei
      else
        kttt<-diag(diag(t(Ei)%*%sigma0%*%Ei))

      if(standrdized)
        kttt<-diag(nrow=M[i])

      cv<-dfrr_fit$sigma_2*kttt
      b0<-b

      vnu0<-t(Ei)%*%X[[i]]%*%b0
      sigmai<-sigma0-sigma0%*%Ei%*%
        solve(t(Ei)%*%sigma0%*%Ei+cv)%*%
        t(Ei)%*%sigma0

        mu1<-X[[i]]%*%b0-sigma0%*%Ei%*%
          solve(t(Ei)%*%sigma0%*%Ei+cv)%*%
          t(Ei)%*%X[[i]]%*%b0

        mu2<-sigma0%*%Ei%*%
          solve(t(Ei)%*%sigma0%*%Ei+cv)

        sigma<-t(Ei)%*%sigma0%*%Ei+cv

        lb<-rep(-Inf,M)
        ub<-rep(Inf,M)
        lb[ys==1]<-0
        ub[ys==0]<-0

      if(M==1)
        zprimes<-tmvtnorm::rtmvnorm(T_G,mean=c(vnu0),sigma=sigma,algorithm = "rejection",
                          lower=lb,upper=ub)
      else
        zprimes<-tmvtnorm::rtmvnorm(T_G,mean=c(vnu0),sigma=sigma,algorithm = "gibbs",
                          lower=lb,upper=ub,burn.in.samples=100,thin=10,
                          start.value=NULL)

      if(M==1)
        zprimes<-matrix(c(zprimes),ncol=M)

      z<-t(sapply(1:nrow(zprimes), function(i1){
          mui<-c(mu1+mu2%*%t(t(zprimes[i1,])))
          MASS::mvrnorm(1,mui,sigmai)
      }))

      Coefs[i,]<-colMeans(z)
      zzt[[i]]<-t(z)%*%z/T_G
    }


  dfrr_fit$pred_data<-list(coefs=Coefs,zzt=zzt,ids=ids,standrdized=standrdized,
                           ydata=list(Y=Ys,time=times,M=Ms))

  if(return.fourier.coefs){
    class(Coefs)<-"predict.dfrr"
    attr(Coefs,"dfrr_fit")<-dfrr_fit
    return(Coefs)
  }



  if(is.null(time_to_evaluate))
    time_to_evaluate<-seq(dfrr_fit$range[1],dfrr_fit$range[2],length.out=100)
  E<-t(fda::eval.basis(time_to_evaluate,dfrr_fit$basis))

  preds<-Coefs%*%E
  if(!is.null(ids))
    rownames(preds)<-ids

    class(preds)<-"dfrr_pred"
    attr(preds,"dfrr_fit")<-dfrr_fit

  return(preds)
}
