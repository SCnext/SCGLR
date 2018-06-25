#@title Function that computes M^(-1/2)
#@export
#@param dt a data.frame containing  covariates
metric <- function(dt){
  dt <- as.data.frame(dt)
  mdt <- stats::model.frame(dt,data=dt)
  nm <- names(mdt)
  out <- lapply(1:ncol(mdt),function(k){
    x <- model.matrix(as.formula(paste0("~",nm[k])),data=mdt)[,-1,drop=FALSE]
    x <- scale(x[,,drop=FALSE],TRUE,FALSE)
    x <- crossprod(x)/nrow(x)
  })
  out <- bdiag(out)
  out <- svd(out)
  if(length(out$d)==1){
    out <- out$u%*%as.matrix(1/sqrt(out$d),1,1)%*%t(out$v)
  }else{
    out <- out$u%*%diag(1/sqrt(out$d))%*%t(out$v)
  }
  return(out)
}
