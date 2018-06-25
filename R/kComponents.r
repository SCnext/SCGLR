# @title Estimation algorithm for K components
# @description calculates the K components by iteratively calling function oneComponent
# @export
# @param X matrix (n*p) containing the standardized covariates
# @param Y matrix (n*q) containing dependent variables
# @param AX matrix of additional covariates used in the generalized regression but not entering the linear
# combinations giving components
# @param K integer specifying the number of components
# @param family a vector of the same length as the number of responses containing characters
# identifying the distribution families of the dependent variables.
# "bernoulli", "binomial", "poisson" or "gaussian" are allowed.
# @param size matrix of size statistical units * number of binomial responses, giving the number of trials
# for binomial dependent variables.
# @param offset used for the poisson dependent variables.
# A vector or a matrix of size: number of observations * number of Poisson dependent variables is expected
# @param crit a list of maxit and tol, default is 50 and 10e-6. If responses are bernoulli variables only, tol should generally be increased
# @param method optimization algorithm used for loadings and components. Object of class "method.SCGLR"
# built by \code{\link{methodEigen}} or \code{\link{methodIng}}
# @return a list
# @return \item{u}{matrix of size (number of regressors * number of components), contains the component-loadings,
# i.e. the coefficients of the regressors in the linear combination giving each component}
# @return \item{comp}{matrix of size (number of statistical units * number of components) having the components as column vectors}
# @return \item{compr}{matrix of size (number of statistical units * number of components) having the standardized components as column vectors}
# @return \item{ds}{the final value of the regularization degree}
kComponents <- function(X,Y,AX,K,family,size=NULL,offset=NULL,crit,method)
{
  n<-dim(X)[1]
  #p<-dim(X)[2]
  q<-dim(Y)[2]
  if (is.null(q)) q<-1
  # on stocke la matrice X et on la centre-reduit
  #Xorig<-X
  F <- NULL
  out <- oneComponent(X,Y,AX,F=F,family=family,size=size,
                      offset=offset,crit=crit,method=method)
  if(is.logical(out))
    custom_stop("convergence_failed")

  u <- out$u
  f <- X%*%u
  F <- f
  Fr <- f/sqrt(c(crossprod(f,f))/n)

  #Composantes suivantes si K>1
  if(K>1) {
    #keep_tol <- crit$tol
    for (k in 2:K) {
      #crit$tol <- keep_tol
      FAX <- cbind(F,AX)
      out <- oneComponent(X,Y,FAX,F,family=family,size=size,
                          offset=offset,crit=crit,method=method)
      if(is.logical(out))
        custom_stop("convergence_failed")

      f <- X%*%out$u
      Fr <- cbind(Fr,f/sqrt(c(crossprod(f,f))/n))
      u <- cbind(u,out$u)
      F <- cbind(F,f)
    }
  }
  u <- as.matrix(u)
  colnames(u) <- paste("u",1:ncol(u),sep="")
  colnames(F) <- paste("sc",1:ncol(F),sep="")
  colnames(Fr) <- paste("sc",1:ncol(Fr),sep="")
  gamma <- out$gamma
  if(is.null(AX)){
    rownames(gamma) <- c("(intercept)",colnames(Fr))
  }else{
    gamma <- gamma[c(1:K,(nrow(gamma)),(K+1):(nrow(gamma)-1)),,drop=FALSE]
    rownames(gamma) <- c("(intercept)",colnames(Fr),colnames(AX))
  }
  return(list(u=u,F=F,Fr=Fr,gamma=gamma))
}


