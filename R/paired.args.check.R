paired.args.check <-
function(a,b,conflict.msg,default){
  if(is.null(a))
  {
    if(is.na(b))
      a<-default
    else
      a<-!as.logical(b)
  }else{
    if(!is.na(b))
      if(as.logical(b)==a)
        stop(paste0("(Conflict in arguments)\\r\\n ",conflict.msg))
  }
  a
}
