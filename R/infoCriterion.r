#' @title Function that calculates cross-validation selection criteria
#' @export
#' @keywords internal
#' @importFrom stats dbinom dpois dnorm
#' @description  Function that calculates cross-validation selection criteria
#' @param ynew data matrix corresponding to the observations used as test sample.
#' @param pred predicted value of the linear predictor obtained from Xnew and the estimated parameters.
#' @param family a vector of the same length as the number of responses containing characters
#' identifying the distribution families of the dependent variables.
#' "bernoulli", "binomial", "poisson" or "gaussian" are allowed.
#' @param type information criterion used. Likelihood, aic, bic, aicc or
#'  Mean Square Prediction Error (mspe) are defined. Area Under ROC Curve (auc) also defined for Bernoulli cases only.
#' @param size describes the number of trials for the binomial dependent variables.
#' A (number of statistical units * number of binomial dependent variables) matrix is expected.
#' @param npar number of parameters used for penalisation.
#' @return a matrix containing the criterion value for each dependent variable (row)
#' and each number of components (column).
infoCriterion <- function(ynew,pred,family,type,size=NULL,npar=0) {
  if(type=="auc")
    stop("auc not allowed as loss function!")
  checkLossFunction(type)

  nobs <- nrow(ynew)
  ny <- ncol(ynew)

  pen <- switch(type,likelihood=0,aic=npar*2,bic=npar*log(nobs),aicc=2*npar(npar+1)/(nobs-npar-1),mspe=0)

  res <- rep(0,ny)
 # success <- ynew
  if(sum("bernoulli" %in% family)>0) {
    tmpy <- ynew[,which(family %in% "bernoulli"),drop=FALSE]
    tmpp <- pred[,which(family %in% "bernoulli"),drop=FALSE]
    if(type=="mspe") {
      ldber <- apply((tmpy-tmpp)^2/(tmpp*(1-tmpp)),2,mean)
    } else {
      ldber <- -2*apply(dbinom(tmpy,1,tmpp,log=TRUE),2,sum)
    }
    res[which(family=="bernoulli")] <- ldber
    success <- NULL
  }
  if(sum("binomial" %in% family)>0) {
    tmpy <- ynew[,which(family %in% "binomial"),drop=FALSE]*size
    tmpp <- pred[,which(family %in% "binomial"),drop=FALSE]

    if(type=="mspe"){
      ldbin <-  apply((tmpy-tmpp)^2/(tmpp*(1-tmpp)*size),2,mean)
    } else {
      ldbin <- -2*apply(dbinom(tmpy,size,tmpp,log=TRUE),2,sum)
    }
    success <- NULL
    res[which(family=="binomial")] <- ldbin
  }

  if(sum("poisson"%in%family)>0){
    tmpy <- ynew[,which(family%in%"poisson"),drop=FALSE]
    tmpp <- pred[,which(family%in%"poisson"),drop=FALSE]
    if(type=="mspe"){
      ldpois <- apply((tmpy-tmpp)^2/tmpp,2,mean)#/tmpp
    } else {
      ldpois <- -2*apply(dpois(tmpy,tmpp,log=TRUE),2,sum)
    }

    #lower <- (tmpp-alpha*sqrt(tmpp))
    #lower[lower<0]<-0
    #upper <- (tmpp+alpha*sqrt(tmpp))
    #success <- (tmpy>=lower)&(tmpy<=upper)

    res[which(family=="poisson")] <- ldpois
  }

  if(sum("gaussian"%in%family)>0){
    tmpy <- ynew[,which(family%in%"gaussian"),drop=FALSE]
    tmpp <- pred[,which(family%in%"gaussian"),drop=FALSE]
    if(type=="mspe"){
      v <-apply(tmpy,2,stats::var)#apply(tmpp,2,var)
      ldgaus <- apply((tmpy-tmpp)^2,2,mean)/v#TODA check equation
    } else {
      sd <- matrix(sqrt(apply((tmpy-tmpp)^2,2,mean)),nobs,ny,byrow=TRUE)
      ldgaus <- -2*apply(dnorm(tmpy,tmpp,sd,log=TRUE),2,sum)## TODO a corriger
    }
    #sd <- matrix(sqrt(apply((tmpy-tmpp)^2,2,mean)),nobs,ny,byrow=TRUE)
    #success[,which(family%in%"gaussian")] <- (tmpy>=(tmpp-alpha*sqrt(sd)))&(tmpy<=(tmpp+alpha*sqrt(sd)))
    res[which(family=="gaussian")] <- ldgaus
  }

  #success <- apply(success,2,mean)
  #return(list(out=res+pen,success=success))
 return(res+pen)
}
