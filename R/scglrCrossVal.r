if(getRversion()>="2.15.1") {
  utils::globalVariables(c("na.omit","coef"))
}
#' Function that fits and selects the number of component by cross-validation.
#' @export
#' @importFrom stats model.matrix model.extract coef cor
#' @importFrom pROC auc roc
#' @param formula an object of class "Formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data the data frame to be modeled.
#' @param family a vector of character of length q specifying the distributions of the responses. Bernoulli, binomial, poisson and gaussian are allowed.
#' @param K number of components, default is one.
#' @param folds number of folds, default is 10.
#' Although folds can be as large as the sample size (leave-one-out CV),
#' it is not recommended for large datasets.
#' folds can also be provided as a vector (same length as data) of fold identifiers.
#' @param type loss function to use for cross-validation.
#' Currently six options are available depending on whether the responses are of the same distribution family.
#' If the responses are all bernoulli distributed, then the prediction performance may be measured
#' through the area under the ROC curve: type = "auc"
#' In any other case one can choose among the following five options ("likelihood","aic","aicc","bic","mspe").
#' @param size specifies the number of trials of the binomial variables included in the model.  A (n*qb) matrix is expected
#' for qb binomial variables.
#' @param offset used for the poisson dependent variables.
#' A vector or a matrix of size: number of observations * number of Poisson dependent variables is expected.
#' @param na.action a function which indicates what should happen when the data contain NAs. The default is set to the \code{na.omit}.
#' @param crit a list of two elements : maxit and tol, describing respectively the maximum number of iterations and
#' the tolerance convergence criterion for the Fisher scoring algorithm. Default is set to 50 and 10e-6 respectively.
#' @param mc.cores deprecated
#' @param method Regularization criterion type. Object of class "method.SCGLR"
#' built by \code{\link{methodSR}} for Structural Relevance.
#' @param nfolds deprecated. Use \code{fold} parameter instead.
#' @return  a matrix containing the criterion values for each response (rows) and each number of components (columns).
#' @references Bry X., Trottier C., Verron T. and Mortier F. (2013) Supervised Component Generalized Linear Regression using a PLS-extension of the Fisher scoring algorithm. \emph{Journal of Multivariate Analysis}, 119, 47-60.
#' @examples \dontrun{
#' library(SCGLR)
#'
#' # load sample data
#' data(genus)
#'
#' # get variable names from dataset
#' n <- names(genus)
#' ny <- n[grep("^gen",n)]    # Y <- names that begins with "gen"
#' nx <- n[-grep("^gen",n)]   # X <- remaining names
#'
#' # remove "geology" and "surface" from nx
#' # as surface is offset and we want to use geology as additional covariate
#' nx <-nx[!nx%in%c("geology","surface")]
#'
#' # build multivariate formula
#' # we also add "lat*lon" as computed covariate
#' form <- multivariateFormula(ny,c(nx,"I(lat*lon)"),A=c("geology"))
#'
#' # define family
#' fam <- rep("poisson",length(ny))
#'
#' # cross validation
#' genus.cv <- scglrCrossVal(formula=form, data=genus, family=fam, K=12,
#'  offset=genus$surface)
#'
#' # find best K
#' mean.crit <- colMeans(log(genus.cv))
#'
#' #plot(mean.crit, type="l")
#' }
scglrCrossVal <-  function(formula,data,family,K=1,folds=10,type="mspe",size=NULL,offset=NULL,
                           na.action=na.omit,crit=list(), method=methodSR(), nfolds, mc.cores) {
  
  if(!missing(nfolds)) {
    .Deprecated("fold", msg="'nfolds' parameter as been renamed to 'folds'. I'll use provided value but update your code!")
    folds <- nfolds
  }
  
  if(!missing(mc.cores)) {
    .Deprecated("plan",msg="mc.cores is now deprecated. Use plan function from future package instead. Value ignored.")
  }
  
  checkLossFunction(type)

  if((type=="auc") && (prod(family=="bernoulli")==0))
    stop("auc loss function only when all bernoulli!")

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula","data","size","offset","na.action"),names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE

  form <- as.Formula(formula)
  mf$formula <- form
  if(!is.null(size))  size <- as.matrix(size)
  mf$size <- size
  if(!is.null(offset)) {
    if(is.vector(offset)) {
      offset <- matrix(offset,nrow(data), sum(family=="poisson"))
    } else {
      offset <- as.matrix(offset)
    }
  }

  mf$offset <- offset
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  crit <- do.call("critConvergence", crit)

  y <- as.matrix(model.part(form,data=mf,lhs=1))
  x <- model.part(form, data=mf, rhs = 1)

  # if additional variable
  if(sum(length(form))==3){
    AX <- model.part(form, data=mf, rhs = 2)
    namesAx <- names(AX)
    AX <- model.matrix(form,data=mf,rhs=2)[,-1]
    if(is.vector(AX)) {
      AX <- matrix(AX,ncol=1)
      colnames(AX) <- namesAx
    }
  }else{
    AX <- NULL
  }

  invsqrtm <- metric(as.data.frame(x))
  xdesign <- model.matrix(form,data=mf)[,-1]
  centerx <- apply(xdesign,2,mean)
  nms <- colnames(xdesign)
  xcr <- scale(xdesign,center=TRUE,scale=FALSE)
  xcr <- xcr%*%invsqrtm
  colnames(xcr) <- nms

  x <- model.matrix(form,data=mf)[,-1]
  ### Controls of  dimension between Y and Size, weights and offsets
  ## number of columns in Y
  ncy <- ncol(y)
  if(length(family)!=ncy){
    stop("number of dependent variables and family attributes are different!")
  }

  if("binomial"%in%family){
    if(is.null(size)){
      stop("Number of trials is unknown for binomial variables!")
    }else{
      if(ncol(size)!=sum("binomial"==family)){
        stop("Number of trials is different from number of binomial variables!")
      }else{
        y[,family%in%"binomial"] <- y[,family%in%"binomial"]/size
      }
    }
  }

  if(!is.null(model.extract(mf,"offset"))){
    if(ncol(offset)!=sum("poisson"==family)){
      stop("Number of offset and poisson variables are different!")
    }
  }


  ###compute the K components of scglr
  size <- model.extract(mf,"size")
  offset <- model.extract(mf,"offset")
  nobs <- nrow(y)
  ny <- ncol(y)
  
  # fold provided as a vector of user groups
  if(length(folds)>1) {
    if(length(folds)!=nobs)
      stop("length of folds must be the same as the number of observations!")
    folds <- as.factor(folds)
    kfolds <- length(levels(folds))
    foldid <- as.integer(folds)
  } else {
    kfolds <- folds
    foldid <- sample(rep(seq(kfolds), length = nobs))
  }
  
  if(kfolds<2) stop("kfolds must be at least equal to two")
  
  cv <- array(0,c(ny,K,kfolds))
  cvNull <- matrix(0,ny,kfolds)

  mainFolds <- function(nf) {
    estid <- which(foldid!=nf)
    valid <- which(foldid==nf)

    cv <- array(0,c(ny,K))
    cvNull <- matrix(0,ny)
    
    ### NULL model
    if(is.null(AX)){
      gamma.fit <- multivariateGlm.fit(Y=y[estid,,drop=FALSE],comp=NULL,family=family,
        offset=offset[estid,,drop=FALSE],
        size=size[estid,,drop=FALSE])
      xnew <- as.matrix(rep(1,length(valid)),ncol=1)
      beta.coefs <- matrix(sapply(gamma.fit, coef),nrow=1,ncol=ny)
    }else{
      gamma.fit <- multivariateGlm.fit(Y=y[estid,,drop=FALSE],comp=AX[estid,,drop=FALSE],family=family,
        offset=offset[estid,,drop=FALSE],
        size=size[estid,,drop=FALSE])
      xnew <- cbind(1,AX[valid,,drop=FALSE])
      beta.coefs <- sapply(gamma.fit, coef)
    }
    
    predict <- multivariatePredictGlm(Xnew=xnew,family,
      beta.coefs,offset[valid,,drop=FALSE])
    
    if(type=="auc"){
      cvNull[1:ny] <- auc(roc(y[valid,,drop=FALSE],predict, quiet=TRUE))
    } else if(type%in%c("likelihood","aic","bic","aicc","mspe")){
      cvNull[1:ny] <- infoCriterion(ynew=y[valid,,drop=FALSE],predict,family,
        type=type,size=size[valid,,drop=FALSE],npar=nrow(gamma.coefs))
    }
    
    try_result <- try({
    kComponent.fit <- kComponents(X=xcr[estid,,drop=FALSE],Y=y[estid,,drop=FALSE],AX=AX[estid,,drop=FALSE],K=K,
                                  family=family,size=size[estid,,drop=FALSE],
                                  offset=offset[estid,,drop=FALSE],crit=crit,method=method)
    },silent=TRUE)
    
    if(inherits(try_result,"try-error")) {
      # return error as a warning
      warning(attr(try_result,"condition")$message,", in fold ",nf, call. = FALSE)
      ..progressor..()
      return(list(cv=array(NA,c(ny,K)),cvNull=cvNull))
    }

    for(kk in seq(K)){

      if(is.null(AX)){
        gamma.fit <- suppressWarnings(multivariateGlm.fit(Y=y[estid,,drop=FALSE],comp=kComponent.fit$F[,1:kk,drop=FALSE],
                                         family=family,offset=offset[estid,,drop=FALSE],size=size[estid,,drop=FALSE]))
      }else{
        gamma.fit <- suppressWarnings(multivariateGlm.fit(Y=y[estid,,drop=FALSE],comp=cbind(kComponent.fit$F[,1:kk,drop=FALSE],AX[estid,,drop=FALSE]),
                                         family=family,offset=offset[estid,,drop=FALSE],size=size[estid,,drop=FALSE]))
      }
      gamma.coefs <- sapply(gamma.fit, coef)

      beta0 <- gamma.coefs[1,] - t(centerx)%*%invsqrtm%*%kComponent.fit$u[,1:kk,drop=FALSE]%*%gamma.coefs[2:(kk+1),]
      beta <- invsqrtm%*%kComponent.fit$u[,1:kk,drop=FALSE]%*%gamma.coefs[2:(kk+1),]
      beta.coefs <- rbind(beta0,beta,as.matrix(gamma.coefs[-c(1:(kk+1)),,drop=F]))
      beta <- as.data.frame(beta.coefs)
      rownames(beta) <- c("Intercept",colnames(xcr),colnames(AX))
      colnames(beta.coefs) <- colnames(y)


      if(is.null(AX)){
        predict <- multivariatePredictGlm(Xnew=as.matrix(cbind(rep(1,length(valid)),x[valid,,drop=FALSE])),family,
                                          beta.coefs[,,drop=FALSE],offset[valid,,drop=FALSE])
      }else{
        predict <- multivariatePredictGlm(Xnew=as.matrix(cbind(rep(1,length(valid)),x[valid,,drop=FALSE],AX[valid,,drop=FALSE])),family,
                                          beta.coefs[,,drop=FALSE],offset[valid,,drop=FALSE])
      }
      
      if(type=="auc"){
        cv[1:ny,kk] <- auc(roc(y[valid,,drop=FALSE],predict, quiet=TRUE))
      } else if(type%in%c("likelihood","aic","bic","aicc","mspe")){
        cv[1:ny,kk]<- infoCriterion(ynew=y[valid,,drop=FALSE],predict,family,
                                    type=type,size=size[valid,,drop=FALSE],npar=nrow(gamma.coefs))
      }

    }
    ..progressor..()
    return(list(cv=cv,cvNull=cvNull))
  }

  ..progressor.. <- getProgressor(kfolds)
  result <- getParallel("lapply", seq(kfolds), mainFolds)
  
  cv <- sapply(result,function(x) x$cv,simplify=FALSE)
  cv <- simplify2array(cv)
  cv <- apply(cv,c(1,2),mean, na.rm=TRUE)

  cvNull <- sapply(result,function(x) x$cvNull,simplify=FALSE)
  cvNull <- simplify2array(cvNull)
  if(is.vector(cvNull)) {
    cvNull <- mean(cvNull, na.rm=TRUE)
  } else {
    cvNull <- apply(cvNull,1,mean, na.rm=TRUE)
  }
  
  cv <- cbind(cvNull,cv)
  colnames(cv) <- c("null model",paste("nc",1:K,sep=""))
  rownames(cv) <- colnames(y)
  
  #class(cv) <- "SCGLRCV"
  return(cv)
}

summary.SCGLRCV <- function(object, ...) {}
print.SCGLRCV <- function(x, ...) {invisible(x)}
plot.SCGLRCV <- function(x, ...) {}
