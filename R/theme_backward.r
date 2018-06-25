#' @title Theme Backward selection
#' @description Perform component selection by cross-validation backward approach
#' @export
#' @param formula an object of class "\code{\link{Formula}}" (or one that can be coerced
#' to that class): a symbolic description of the model to be fitted. The details
#' of model specification are given under Details.
#' @param data data frame.
#' @param H vector of R integer. Number of components to keep for each theme
#' @param kfolds number of folds - default is 10. Although nfolds can be as large as the sample size (leave-one-out CV),
#' it is not recommended for large datasets. Smallest value allowable is nfolds=2
#' Models for theme are specified symbolically. A model as the form \code{response ~ terms} where \code{response}
#' is the numeric response vector and terms is a series of R themes composed of
#' predictors. Themes are separated by  "|" (pipe) and are composed.... Y1+Y2+...
#' ~ X11+X12+...+X1_  | X21+X22+... | ...+X1_+...  | A1+A2+... See \code{\link{multivariateFormula}}.
#' @param family a vector of character of the same length as the number of dependent variables:
#' "bernoulli", "binomial", "poisson" or "gaussian" is allowed.
#' @param size describes the number of trials for the binomial dependent variables.
#' A (number of statistical units * number of binomial dependent variables) matrix is expected.
#' @param weights weights on individuals (not available for now)
#' @param offset used for the poisson dependent variables.
#' A vector or a matrix of size: number of observations * number of Poisson dependent variables is expected.
#' @param na.action a function which indicates what should happen when the data contain NAs. The default is set to \code{na.omit}.
#' @param crit a list of two elements : maxit and tol, describing respectively the maximum number of iterations and
#' the tolerance convergence criterion for the Fisher scoring algorithm. Default is set to 50 and 10e-6 respectively.
#' @param method structural relevance criterion. Object of class "method.SCGLR"
#' built by  \code{\link{methodSR}} for Structural Relevance.
#' @param type loss function to use for cross-validation.
#' Currently six options are available depending on whether the responses are of the same distribution family.
#' If the responses are all bernoulli distributed, then the prediction performance may be measured
#' through the area under the ROC curve: type = "auc"
#' In any other case one can choose among the following five options ("likelihood","aic","aicc","bic","mspe").
#' @param st logical (FALSE) theme build and fit order. TRUE means random, FALSE means sequential (T1, ..., Tr)

#' @return  a list containing the path followed along the selection process, the associated mean square predictor error and the best configuration.
#' @examples \dontrun{
#' library(SCGLR)
#'
#' # load sample data
#' data(genus)
#'
#' # get variable names from dataset
#' n <- names(genus)
#' n <-n[!n%in%c("geology","surface","lon","lat","forest","altitude")]
#' ny <- n[grep("^gen",n)]    # Y <- names that begins with "gen"
#' nx1 <- n[grep("^evi",n)]   # X <- remaining names
#' nx2 <- n[-c(grep("^evi",n),grep("^gen",n))]
#'
#'
#'form <- multivariateFormula(ny,nx1,nx2,A=c("geology"))
#'fam <- rep("poisson",length(ny))
#'testcv <- scglrThemeBackward(form,data=genus,H=c(2,2),family=fam,offset = genus$surface,kfolds=2)
#'Cross-validation pathway
#'testcv$H_path
#'Plot criterion
#'plot(testcv$cv_path)
#'Best combination
#'testcv$H_best
#' }
scglrThemeBackward <- function(formula, data, H, family, size = NULL, weights = NULL,
                  offset = NULL, na.action = na.omit, crit = list(), method = methodSR(), kfolds=10,type="mspe",st=FALSE){


  if(!inherits(formula,"MultivariateFormula"))
    formula <- multivariateFormula(formula,data=data)

  additional=formula$additional

  # check family
  Y_vars <- attr(formula,"Y_vars")
  if(!is.character(family)||any(!(family %in% c("gaussian","poisson","bernoulli","binomial"))))
    stop("Expecting character vector of gaussian, poisson, bernoulli or binomial for family")
  if(!(length(family) %in% c(1,length(Y_vars))))
    stop("Length of family must be equal to one or number of Y variables")
  if(length(family)==1)
    family <- rep(family, length(Y_vars))

  if(kfolds<2) stop("kfolds must be at least equal to two")

  X_expand <- lapply(formula$X,function(X){
    f <- as.formula(paste0("~",paste0(trim(deparse(X)),collapse=" ")))
    return(model.matrix(f,data)[,-1])
  })

  invsqrtm <- lapply(formula$X, function(X){
    nms <- all.vars(as.formula(paste0("~",paste0(trim(deparse(X)),collapse=" "))))
    metric(data[,nms])
  })
  X_expand <- lapply(1:length(H),function(l) scale(X_expand[[l]],TRUE,FALSE)%*%invsqrtm[[l]])
  if(formula$additional)
    A_expand <- model.matrix(as.formula(paste0("~",paste0(trim(deparse(formula$A)),collapse=" "))),data)[,-1,drop=F]
  else
    A_expand <- NULL


  folds <- sample(1:kfolds, nrow(data),replace = TRUE)
#Full model evaluation
  message("full model")
  thm <- lapply(1:kfolds,function(k) {
    out <- scglrTheme(
      formula=formula,
      data=data,
      H=H,
      family=family,
      size=size,
      weights=weights,
      offset=offset,
      subset=(1:nrow(data))[folds!=k],
      na.action = na.action,
      crit=crit,
      method=method,
      st=st
    )
  result <- list(
    u=lapply(out$themes,function(t) as.matrix(t$u)),
    gamma=out$gamma
  )
  return (result)
  })

  cv <- lapply(1:kfolds,function(k){
    XU_new <- as.matrix(do.call(cbind,lapply(which(H>0),function(l) as.matrix(X_expand[[l]][folds==k,,drop=F])%*%thm[[k]]$u[[l]])))
    if(!is.null(A_expand))
      A_new <- A_expand[folds==k,,drop=F]
    else
      A_new <- NULL
    X_new <- cbind(1,XU_new,A_new)
    pred <- multivariatePredictGlm(X_new,family=family,beta=thm[[k]]$gamma,offset = offset[folds==k])
    if(!(type%in%c("mspe","auc"))) npar <- ncol(X_new) else npar=0
    qual <- infoCriterion(ynew=as.matrix(data[folds==k,formula$Y_vars]),pred=pred,family=family,type=type,size=size[folds==k,,drop=FALSE],npar=npar)
  })
  cv <- mean(log(Reduce("+",cv)/kfolds))
  message("[",paste(H,collapse=","),"] = ",cv)

  #backward evaluation
  message("backward")
  H_cur <- H
  cv_path <- cv
  H_path <- list(H)
  while(sum(H_cur)>1) {
    H_new <- lapply(which(H_cur>0),function(i) {h <- H_cur; h[i]<-h[i]-1;h})
    cv_new <- lapply(H_new,function(h){
      cv <- lapply(1:kfolds, function(k){
        XU_fit <- as.matrix(do.call(cbind,lapply(which(h>0),function(l) as.matrix(X_expand[[l]][folds!=k,,drop=F])%*%thm[[k]]$u[[l]][,1:h[l],drop=FALSE])))
        colnames(XU_fit) <- paste0("c",1:ncol(XU_fit))
        if(!is.null(A_expand))
          A_fit <- A_expand[folds!=k,,drop=FALSE]
        else
          A_fit <- NULL
        X_fit <- cbind(XU_fit,A_fit)
        fit <- multivariateGlm.fit(Y=data[folds!=k,formula$Y_vars],comp=X_fit,family=family,offset=offset[folds!=k],size=size[folds!=k,,drop=FALSE])
        gamma <- sapply(fit, coef)

        XU_new <- as.matrix(do.call(cbind,lapply(which(h>0),function(l) as.matrix(X_expand[[l]][folds==k,,drop=F])%*%thm[[k]]$u[[l]][,1:h[l],drop=FALSE])))
        colnames(XU_new) <- paste0("c",1:ncol(XU_new))
        if(!is.null(A_expand))
          A_new <- A_expand[folds==k,,drop=FALSE]
        else
          A_new <- NULL
        X_new <- cbind(1,XU_new,A_new)

        pred <- multivariatePredictGlm(X_new,family=family,beta=gamma,offset = offset[folds==k])
        if(!(type%in%c("mspe","auc"))) npar <- ncol(X_new) else npar=0

        qual <- infoCriterion(ynew=as.matrix(data[folds==k,formula$Y_vars]),pred=pred,family=family,type=type,size=size[folds==k,,drop=FALSE],npar=npar)
      })
      return(mean(log(Reduce("+",cv)/kfolds)))
    })
    cv_new <- unlist(cv_new)
    H_cur <- H_new[[which.min(cv_new)]]

    cv_path <- c(cv_path,min(cv_new))
    H_path <- c(H_path,list(H_cur))
    message("[",paste(H_cur,collapse=","),"] = ",min(cv_new))
  }

  message("NULL model")
  cvNull <- lapply(1:kfolds,function(k){
        if(!is.null(A_expand))
          A_fit <- A_expand[folds!=k,,drop=FALSE]
        else
          A_fit <- NULL
        fit <- multivariateGlm.fit(Y=data[folds!=k,formula$Y_vars],comp=A_fit,family=family,offset=offset[folds!=k],size=size[folds!=k,,drop=FALSE])
        gamma <- sapply(fit, coef)
        if(!is.null(A_expand))
          A_new <- A_expand[folds==k,,drop=FALSE]
        else
          A_new <- NULL
        X_new <- cbind(rep(1,sum(folds==k)),A_new)
        pred <- multivariatePredictGlm(X_new,family=family,beta=gamma,offset = offset[folds==k])
        if(!(type%in%c("mspe","auc"))) npar <- ncol(X_new) else npar=0
        qual <- infoCriterion(ynew=as.matrix(data[folds==k,formula$Y_vars]),pred=pred,family=family,type=type,size=size[folds==k,,drop=FALSE],npar=npar)
  })
  message("[",paste(rep(0,length(H)),collapse=","),"] = ",mean(log(Reduce("+",cvNull)/kfolds)))
  cv_path <- c(cv_path,mean(log(Reduce("+",cvNull)/kfolds)))
  H_path <- c(H_path,list(rep(0,length(H))))
  H_path=do.call(rbind,H_path)
  colnames(H_path) <- paste0("thm",1:length(H))
  return(list(H_path=H_path,cv_path=cv_path,H_best=H_path[which.min(cv_path),]))
}
