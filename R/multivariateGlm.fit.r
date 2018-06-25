#' @title Multivariate generalized linear regression
#' @description \code{multivariateGlm} is used to fit multivariate generalized linear models
#' specified by a symbolic formula together with the distributions of the responses.
#' This function performs a simple GLM fit for each dependent variable with the associated distribution.
#' @export
#' @importFrom stats glm
#' @param Y matrix of dependent variables.
#' @param comp matrix of covariates.
#' @param family a vector of character giving the family distribution of each response.
#' @param size a matrix giving the number of trials for each Binomial dependent variable
#' ncol(size) must be equal to the number of Binomial variables.
#' @param offset used for the poisson dependent variables.
#' A vector or a matrix of size: number of observations * number of Poisson dependent variables is expected.
#' @return the list, each item of which is the glm object associated with each response.
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
#' # remove "geology" and "surface" from nx as surface
#' # is offset and we want to use geology as additional covariate
#' nx <-nx[!nx%in%c("geology","surface")]
#'
#' # build multivariate formula
#' # we also add "lat*lon" as computed covariate
#' form <- multivariateFormula(ny,c(nx,"I(lat*lon)"),c("geology"))
#'
#' # split genus dataset
#' sub <- sample(1:nrow(genus),100,replace=FALSE)
#' sub_fit <- (1:nrow(genus))[-sub]
#'
#' # define family
#' fam <- rep("poisson",length(ny))
#'
#' # fit the model
#' genus.scglr <- scglr(formula=form, data=genus, family=fam, K=4,
#'  offset=genus$surface, subset=sub_fit)
#'
#' # xnew, the design matrix associated to sub-sample used for prediction
#' # note rhs parameter is introduced to take into account that the
#' # covariate part of the formula is composed of two differents sets
#' xnew <- model.matrix(form, data=genus[sub,], rhs=1:2)[,-1]
#'
#' # prediction based on the scglr approch
#' pred.scglr <- multivariatePredictGlm(xnew,family=fam,
#'  beta=genus.scglr$beta, offset=genus$surface[sub])
#' cor.scglr <-diag(cor(pred.scglr,genus[sub,ny]))
#' plot(cor.scglr, col="red",ylim=c(-1,1))
#'
#' # prediction based on classical poisson glm
#' X <- model.matrix(form, data=genus)[,-1]
#' Y <- genus[,ny]
#' genus.glm <- multivariateGlm.fit(Y[sub_fit,,drop=FALSE],X[sub_fit,,drop=FALSE],
#'              family=fam, offset=matrix(genus$surface[sub_fit],
#'              length(sub_fit),length(ny)),size=NULL)
#' coefs <- sapply(genus.glm,coef)
#'
#' # rhs parameter is introduced to take into account that the
#' # covariate part of the formula is composed of two differents sets
#' pred.glm <- multivariatePredictGlm(xnew,family=fam,beta=coefs,
#'  offset=genus$surface[sub])
#' cor.glm <- diag(cor(pred.glm,genus[sub,ny]))
#'
#' points(cor.glm, col="blue")
#' }
multivariateGlm.fit <- function(Y,comp,family,offset,size){
  q <- ncol(Y)
  out <- list()
  if(!is.null(comp)){
    if (is.matrix(comp)|| is.data.frame(comp))
    {
      if(is.null(colnames(comp))){
      colnames(comp) <- paste("c",1:ncol(comp),sep="")
      }
      form <- as.formula(paste("obs~",paste(colnames(comp),collapse="+")))
      ##out <- matrix(0,ncol(comp)+1,ncol(Y))
    }else
    {
      form <- as.formula("obs~comp")
      ##out <- matrix(0,2,ncol(Y))
    }
  }else{
    form <- as.formula("obs~1")
    comp <- rep(1,nrow(Y))
  }

  kpois <- kbern <-  1
  for(i in seq(q)){
    if(family[i]=="poisson"){
      if(is.null(offset)){
        out[[i]] <- glm(form,data=data.frame(obs=Y[,i],comp),family="poisson")
      }else{
        if(is.vector(offset))
          offset <- matrix(offset,length(offset),ncol(Y))
        out[[i]] <- glm(form,data=data.frame(obs=Y[,i],comp),family=family[i],
                        offset=log(offset[,kpois]))
        kpois <- kpois+1
      }
    }
    if(family[i]=="bernoulli"){
      out[[i]] <- glm(form,data=data.frame(obs=Y[,i],comp),family="binomial")
    }
    if(family[i]=="binomial"){
      out[[i]] <- glm(form,data=data.frame(obs=Y[,i],comp),family="binomial",
                      weights=size[,kbern])
      kbern <- kbern+1
    }
    if(family[i]=="gaussian"){
      out[[i]] <- glm(form,data=data.frame(obs=Y[,i],comp),family="gaussian")
    }
  }

  return(out)
}

