#' Dichotomized Functional Response Regression
#'
#' Implementing Function-on-Scalar Regression model, in which the response function
#'  is dichotomized and observed sparsely.
#'
#' The output is a \code{dfrr}-object, which then can be injected into other methods/functions
#'  to postprocess the fitted model, including:
#'  \code{\link{coefs.dfrr}},\code{\link{fitted.dfrr}}, \code{\link{residuals.dfrr}}, \code{\link{predict.dfrr}},
#'   \code{\link{eigen.dfrr}}, \code{\link{summary.dfrr}}, \code{\link{qq.dfrr}}, \code{\link{model.matrix.dfrr}},
#'   \code{\link{plot.coefs.dfrr}}, \code{\link{plot.fitted.dfrr}}, \code{\link{plot.residuals.dfrr}},
#'   \code{\link{plot.predict.dfrr}}, \code{\link{plot.eigen.dfrr}}, \code{\link{plot.residuals.dfrr}}
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
#' summary(dfrr_fit)
#'
#' ##### Fitting dfrr model to the Madras Longitudinal Schizophrenia data
#' data(madras)
#' ydata<-data.frame(.obs=madras$id,.index=madras$month,.value=madras$y)
#' ids<-unique(madras$id)
#' q<-4
#' N<-length(ids)
#' xData<-data.frame(Age=rep(NA,N),Gender=rep(NA,N))
#' for(i in 1:N){
#'   dt<-madras[madras$id==ids[i],]
#'   xData[i,]<-c(dt$age[1],dt$gender[1])
#' }
#' rownames(xData)<-ids
#'
#' madras_dfrr<-dfrr(Y~Age+Gender+Age*Gender, data=xData, ydata=ydata, J=11,T_E=5)
#' coefs<-coef(madras_dfrr)
#' plot(coefs)
#'
#' fpcs<-fpca(madras_dfrr)
#' plot(fpcs)
#'
#' @param formula an object of class "\code{\link[stats]{formula}}" (or one that can be coerced to that class with \code{\link[stats]{as.formula}}:
#'  a symbolic description of the model to be fitted.
#' @param yind a vector with length equal to the number of columns of the matrix of functional
#'  responses giving the vector of evaluation points \eqn{(t_1,...,t_{G})}.
#'  If not supplied, \code{yind} is set to \code{1:ncol(<response>)}.
#' @param data an (optional) \code{data.frame} containing the covariate data.
#' the variable terms will be searched from the columns of \code{data},
#' covariates also can be read from the workspace if it is not available in \code{data}.
#' @param ydata an (optional) \code{data.frame} consists of three columns \code{.obs}, \code{.index} and \code{.value},
#'  supplying the functional responses that are not observed on a regular grid.
#'  ydata must be provided if the sampling design is irregular.
#' @param method detrmines the estimation method of functional parameters. Defaults to "\code{REML}" estimation.
#' @param rangeval an (optional) vector of length two, indicating the lower and upper limit of the domain of
#' latent functional response. If not specified, it will set by minimum and maximum of \code{yind} or \code{.index} column of \code{ydata}.
#' @param ... other arguments that can be passed to the inner function \code{AMCEM}.
#' @param basis an (optional) object of class \code{'basisfd'}. Defaults to cubic bspline basis.
#' @export
dfrr <-
function(formula, yind=NULL, data = NULL, ydata = NULL,
               method = c("REML","ML"),rangeval=NULL,basis=NULL,...){

    formula<-as.formula(formula)


  method<-match.arg(method)
  times_to_evaluate<-NULL
  is.regular<-is.null(ydata)

  if(is.regular){
    if(!is.null(ydata))
    {
      if(!all(c(".obs",".index",".value") %in% colnames(ydata)))
        stop("Argument 'ydata' is not of the expected structure. See the help for more details")
      if(!all(unique(ydata$.index) %in% yind))
        stop("'yind' does not contain unique values of '.index' column of ydata")
    }

  }else{
    if(is.null(ydata))
      stop("User must provide either ydata or yind")
    if(!all(c(".obs",".index",".value") %in% colnames(ydata)))
      stop("Argument 'ydata' is not of the expected structure. See the help for more details")
  }

  if(!is.null(data))
    if(is.null(rownames(data)))
      rownames(data)<-1:N

  data_ids<-rownames(data)

  if(!is.null(ydata)){
    if(!is.null(data))
    {
      int<-intersect(unique(ydata$.obs) ,data_ids)
      q<-length(attr(terms(formula),"term.labels"))
      if(length(int)==0){
        stop("Rownames of argument 'data' does not match with the unique values of column '.obs' in 'ydata'")
      }else if(length(int)<=q){
        stop("Rownames of argument 'data' does not match with the unique values of column '.obs' in 'ydata'.\\r\\n Not enough sample to fit the model.")
      }else if(length(int)!=nrow(data)){
        data<-data[int,]
        N<-nrow(data)
        excluded_ids<-setdiff(unique(ydata$.obs),int)
        excluded_ids<-paste0(excluded_ids,", ",collapse = "")
        excluded_ids<-substr(excluded_ids,1,nchar(excluded_ids)-2)
        warning(paste0("Rownames of argument 'data' does not completely match with the unique values of column '.obs' in 'ydata'.\\r\\n",
                       "The following ", nrow(data)-length(int)," cases excluded from ydata. \\r\\n .obs=",excluded_ids))
      }

    }
  }

  #Analyze formula
  # Specify response variable
  # specify coariates' varialbes
  # check for variables' environments

  #Identifying response variable
  res<-paste0(formula)
  resp_varname<-""
  if(length(res)==3)
    resp_varname<-res[2]

  if(resp_varname=="" & is.null(ydata))
    stop("response variable must be specified")

  check_wsp<-FALSE
  y_matrix<-NULL
  if(is.null(ydata)){
    check_wsp<-TRUE
    y_matrix<-eval({as.symbol(resp_varname)},envir = parent.frame())
  }else{
    # if(!is.na(resp_varname))
    #   if(resp_varname!="")
    #     invisible(cat(paste0("Argument 'ydata' is provided. The variable name '",resp_varname,"' from the left-hand side of the formula is is read from '.value' column of 'ydata'\\r\\n")))
  }

  tof<-terms(formula)
  var_names<-rownames(attr(tof,"factors"))
  var_names<-setdiff(var_names,resp_varname)

  if(!is.null(ydata))
    N<-length(unique(ydata$.obs))
  else
    N<-nrow(y_matrix)

  if(is.null(data))
    data<-data.frame(.int=rep(1,N))


  if(!is.null(data))
  {
    if(nrow(data)!=N)
      stop(paste0("The number of rows in data frame is different from number of samples in the response (N=",N,")"))

  }

  if(!is.null(var_names)){
    wsp_var_names<-setdiff(var_names,colnames(data))

    if(!is.null(wsp_var_names) & length(wsp_var_names)>0){
      for(i1 in 1:length(wsp_var_names)){
        tmp_<-eval({as.symbol(wsp_var_names[i1])},envir = parent.frame())
        if(length(tmp_)!=N)
          stop(paste0("Variable '",wsp_var_names[i1],"' has a length different from number of samples in the response (N=",N,")"))
        data[,wsp_var_names[i1]]<-tmp_
      }

      if(!is.null(ydata)){
        ids<-intersect(1:N,unique(ydata$.obs))
        if(length(ids)!=N)
          stop(paste0("The ids in '.obs' column in ydata must be the sequence 1,...,N; N=",N, " is thenmber of samples"))
      }

    }
  }

  # if(!is.null(ydata)){
  #   ids<-unique(ydata$.obs)
  #   data<-data[ids,]
  # }else
  if(is.null(ydata)){
    if(is.null(rownames(y_matrix))){
      ids<-1:N
      rownames(y_matrix)<-1:nrow(y_matrix)
    }else{
      ids<-rownames(y_matrix)
      u<-intersect(rownames(y_matrix),data_ids)
      if(length(u)!=length(rownames(y_matrix)))
        stop("Rownames of response matrix does not meatch with rownames of 'data' matrix")

      #data<-data[rownames(y_matrix),]
    }
  }


  formula2<-formula
  if(length(res)==3)
  {
    formula_string<-paste0("~",res[3])
    formula2<-formula(formula_string)
  }

  xData<-model.matrix(formula2,data=data)
  if(ncol(xData)==0)
    xData<-NULL


  time<-list()
  Y<-list()
  ids_nna<-c()
  i_nna<-c()

  kk<-0
  if(!is.null(ydata)){

    for(i in 1:N){
      if(!is.null(xData) & all(!is.na(xData[i,]))){
        ys<-ydata[ydata$.obs==ids[i],]
        ys<-ys[!is.na(ys$.value) & !is.na(ys$.index),]
        if(length(ys)==0)
          next
        kk<-kk+1
        time[[kk]]<-ys$.index
        Y[[kk]]<-ys$.value
        ids_nna[kk]<-ids[i]
        i_nna[kk]<-i
      }
    }

  }else{

    for(i in 1:N){
      if(!is.null(xData) & all(!is.na(xData[i,]))){
        ys<-y_matrix[i,]
        ty<-yind[!is.na(ys)]
        if(length(ty)==0)
          next
        kk<-kk+1
        time[[kk]]<-ty
        Y[[kk]]<-ys[!is.na(ys)]
        ids_nna[kk]<-ids[i]
        i_nna[kk]<-i
      }
    }

  }


  if(kk==0)
    stop("There is no sample with complete with data")

  if(length(i_nna)!=N)
    if(!is.null(xData))
      xData<-xData[i_nna,]

  if(!is.null(xData))
    if(qr(xData)$rank<ncol(xData))
      stop("xData does not have full column rank")

  timeT<-sort(unique(unlist(time)),decreasing = FALSE)
  if(is.null(rangeval)){
    rangeval<-c(min(timeT),max(timeT))
  }

  if(length(rangeval)!=2)
    stop("rangeval must be a vector of length 2")

  if(rangeval[1]>min(timeT) | rangeval[2]<max(timeT))
    stop("rangeval does not include some of time points")

  if(is.null(times_to_evaluate))
    times_to_evaluate<-seq(rangeval[1],rangeval[2],length.out = 100)

  #Pass Y, time, xData, ids_nna to AMCEM3


  # if(!is.null(ydata)){
  #   timeT<-sort(unique(unlist(time)),decreasing = FALSE)
  #   if(all(sapply(1:N, function(i){length(intersect(tim[[i]],timeT))==length(timeT)})))
  #     design<-"R"
  # }

  res<-AMCEM(Y,time,xData,ids = ids_nna,method = method,
             rangeval=rangeval,times_to_evaluate=times_to_evaluate,basis=basis,...)

  if(!is.null(ydata))
    res$ydata<-ydata
  else{
    res$y_matrix<-y_matrix
    res$yind<-yind
    res$ids_rows<-ids
  }
  class(res)<-"dfrr"
  attr(res,"formula")<-formula2

  res

}
