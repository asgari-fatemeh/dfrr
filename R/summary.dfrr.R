#'Summary for a dfrr fit
#'
#'Summarise a fitted \code{dfrr}-object. Not implemented.
#'
#'@param object a \code{dfrr}-object
#'@param ... dot argument, just for consistency with the generic function
#'
#'@method summary dfrr
#'
#'@export
#'
summary.dfrr <-
function(object,..){
  dfrr_fit<-object
  warning("Not yet implemented...\r\n Use the functions coef(), fpca(), fitted(), and residuals() to see the outputs.")
}
