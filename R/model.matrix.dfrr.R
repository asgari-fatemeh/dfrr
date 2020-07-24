#' Obtain model matrix for a dfrr fit
#'
#' Obtain model matrix for a dfrr fit
#'
#'
#' @inheritParams summary.dfrr
#'
#' @method model.matrix dfrr
#'
#' @export

model.matrix.dfrr<-function(object,...){
  dfrr_fit<-object
    dfrr_fit$modelMatrix
}
