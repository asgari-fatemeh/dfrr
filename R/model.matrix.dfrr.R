#' Obtain model matrix for a dfrr fit
#'
#' Obtain model matrix for a dfrr fit
#'
#' @inheritParams coef.dfrr
#'
#' @method model.matrix dfrr
#' @export

model.matrix.dfrr<-function(dfrr_fit){
    dfrr_fit$modelMatrix
}
