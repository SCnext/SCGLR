#######################################################################
#   Functions to compute k components via the PING algorithm in the   #
#        special case of multivariate grouped count data              #
#######################################################################





# Auxiliary function (computation of 1 component)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Y:       matrix of random responses
# X:       matrix of the normalized explanatory variables
# AX:      matrix of the additional covariates
# u:       vector of initial loadings
# fmat:    matrix of the previously calculated components
# random:  vector giving the group of each unit (factor)
# loffset: matrix of size (nrow(Y), ncol(Y)) giving the log of the offset
#          relative to each unit
# method:  structural relevance criterion
oneCompRand <- function(Y, family, size=NULL,
                        X, AX=NULL, u=NULL, fMat=NULL, random, loffset=NULL,
                        init.sigma = rep(1, ncol(Y)), resACP,
                        method=methodSR("vpi", l=4, s=1/2,
                                        maxiter=1000, epsilon=10^-6, bailout=1000)){
  
  # control
  muinf <- 1e-5
  if (is.null(q)) q<-1
  
  # Useful dimensions
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]
  
  # Design random effects (grouped data)
  designXi <- model.matrix(~factor(random)-1)#,data=random)
  R <- table(random)
  N <- length(R)
  
  # initialisation working variables Z, weights W, and variance components sigma
  mu0 <- apply(Y, 2, mean)
  mu0 <- matrix(mu0, n, q, byrow = TRUE)
  Z <- Y
  
  if("bernoulli"%in%family)
  {
    tmu0 <- mu0[,family=="bernoulli"]
    Z[,family=="bernoulli"] <-  log(tmu0/(1-tmu0))+(Y[,family=="bernoulli"]-tmu0)/(tmu0*(1-tmu0))
  }
  if("binomial"%in%family)
  {
    tmu0 <- mu0[,family=="binomial"]
    Z[,family=="binomial"] <-   log(tmu0/(size-tmu0))+(Y[,family=="binomial"]-tmu0)/((tmu0*(size-tmu0))/size)
  }
  if("poisson"%in%family)
  {
    tmu0 <- mu0[,family=="poisson"]
    if(is.null(loffset)){
      Z[,family=="poisson"]<- log(tmu0)+(Y[,family=="poisson"]-tmu0)/tmu0
    }else{
      Z[,family=="poisson"]<- log(tmu0)-loffset+(Y[,family=="poisson"]-tmu0)/tmu0
    }
  }
  if("gaussian"%in%family){
    Z[,family=="gaussian"] <- Y[,family=="gaussian"]
  }
  
  # Z <- log(mu0/(50-mu0)) + (Y - mu0)/((mu0*(50-mu0))/50)
  # Z <- log(mu0) - loffset + (Y - mu0)/mu0
  
  W <- matrix(1,sum(R),q)#/sum(R)
  sigma <- init.sigma#rep(1,q)#rgamma(q,1,1)
  eta.old <- 1e10
  
  
  
  
  
  
  
  # INITIAL LOOP
  # ------------------------------------------------------------
  
  # 1. PING
  #########
  # resACP argument PING
  u.new = mixed.ping(Z=Z, X=X, AX=AX, W=W, F=fMat, u=u, method=method, resACP=resACP)
  if(p>n){
    # utilPCA = eigen(crossprod(x = X)/n,symmetric = TRUE)
    # pourcent.variance = utilPCA$values/sum(utilPCA$values)
    # nbcomp = max( which(c(pourcent.variance > (1/p))==TRUE ))
    # Xpc = X%*%utilPCA$vectors[,1:nbcomp]
    # f = Xpc %*% u.new
    utilPCA = resACP$utilPCA
    pourcent.variance = resACP$pourcent.variance
    nbcomp = resACP$nbcomp
    Xpc = resACP$Xpc
    f = Xpc %*% u.new
  }else{
    f = X %*% u.new
  }
  T2 = cbind(1,AX,fMat,f)
  
  # 2. Henderson's systems
  ########################
  m1 = crossprod(T2, T2)
  m2 = crossprod(designXi, T2)
  m3 = t(m2)#crossprod(T2,designXi)
  m4 = crossprod(designXi, designXi) + 1/sigma[1] * diag(N)
  b1 = crossprod(T2, Z)
  b2 = crossprod(designXi, Z)
  A <- rbind(cbind(m1, m3), cbind(m2,m4))
  b = rbind(b1, b2)
  coeffs.courants <-  solve(A,b)
  
  # 3. Variance components
  ########################
  for(k in 1:q){
    sigma[k] = (crossprod(coeffs.courants[(ncol(T2)+1):(ncol(T2)+N),k]))/
      (N - 1/sigma[k] * sum(diag(solve(crossprod(designXi, W[,k]*designXi) +
                                         1/as.vector(sigma[k]) * diag(N)))))
  }
  
  # 4. Update
  ###########
  eta=cbind(T2,designXi)%*%coeffs.courants#[1:(ncol(T)+2)]
  
  # etainf <- log(muinf)
  # indinf <- 1 * (eta < etainf)
  # eta <- eta * (1 - indinf) + etainf * indinf
  # mu = exp(eta+loffset)
  # Z = eta + (Y - mu)/mu
  # W <- mu
  if("poisson"%in%family)
  {
    etainf <- log(muinf)
    indinf<-1*(eta[,family=="poisson"]<etainf)
    eta[,family=="poisson"]<-eta[,family=="poisson"]*(1-indinf)+etainf*indinf
    if(is.null(loffset))
    {
      mu <- exp(eta[,family=="poisson"])
    }else{
      mu <- exp(eta[,family=="poisson"]+loffset)
    }
    Z[,family=="poisson"]<- eta[,family=="poisson"]+(Y[,family=="poisson"]-mu)/mu
    W[,family=="poisson"] <- mu
  }
  if("bernoulli"%in%family)
  {
    etainf <- log(muinf/(1-muinf))
    indinf<-1*(eta[,family=="bernoulli"]<etainf)
    eta[,family=="bernoulli"]<-eta[,family=="bernoulli"]*(1-indinf)+etainf*indinf
    indsup<-1*(eta[,family=="bernoulli"]>-etainf)
    eta[,family=="bernoulli"]<-eta[,family=="bernoulli"]*(1-indsup)-etainf*indsup
    mu <- exp(eta[,family=="bernoulli"])/(1+exp(eta[,family=="bernoulli"]))
    Z[,family=="bernoulli"] <-  eta[,family=="bernoulli"] + (Y[,family=="bernoulli"]-mu)/(mu*(1-mu))
    W[,family=="bernoulli"] <- mu*(1-mu)
  }
  if("binomial"%in%family)
  {
    etainf <- log(muinf/(1-muinf))
    indinf<-1*(eta[,family=="binomial"]<etainf)
    eta[,family=="binomial"]<-eta[,family=="binomial"]*(1-indinf)+etainf*indinf
    indsup<-1*(eta[,family=="binomial"]>-etainf)
    eta[,family=="binomial"]<-eta[,family=="binomial"]*(1-indsup)-etainf*indsup
    mu <- size*exp(eta[,family=="binomial"])/(1+exp(eta[,family=="binomial"]))
    Z[,family=="binomial"] <-  eta[,family=="binomial"] + (Y[,family=="binomial"]-mu)/((mu*(size-mu))/size)
    W[,family=="binomial"] <- (mu*(size-mu))/size
  }
  if("gaussian"%in%family){
    Z[,family=="gaussian"]<-Y[,family=="gaussian"]
  }
  
  
  
  
  u.old  <-  u.new
  sigma.old = sigma
  coeffs.courants.old = coeffs.courants
  
  deltaEta = 10
  deltaU <- 10
  deltaSigma <-10
  # ------------------------------------------------------------
  
  
  
  # WHILE LOOP
  # ------------------------------------------------------------
  compt <- 0
  while( ((deltaU > 10^(-6)) || (deltaSigma > 10^(-6)) || (deltaEta> 10^(-6)))
         && (compt<100)
  ){
    
    # 1. PING
    #########
    u.new = mixed.ping(Z=Z, X=X, AX=AX, W=W, F=fMat, u=u.old, method=method, resACP=resACP)
    if(p>n){
      f = Xpc %*% u.new
    }else{
      f = X %*% u.new
    }
    T2 = cbind(1,AX,fMat,f)
    
    # 2. Henderson's systems
    ########################
    for(k in 1:q){
      m1 = crossprod(T2, W[,k]*T2)
      m2 = crossprod(designXi, W[,k]*T2)
      m3 = t(m2)#crossprod(T2, W[,k]*designXi)
      m4 = crossprod(designXi, W[,k]*designXi) + 1/sigma[k] * diag(N)
      b1 = crossprod(T2, W[,k]*Z[,k])
      b2 = crossprod(designXi, W[,k]*Z[,k])
      A <- rbind(cbind(m1, m3), cbind(m2,m4))
      b = rbind(b1, b2)
      coeffs.courants[,k] = solve(A,b)
      
      # 3. Variance components
      ########################
      sigma[k] = (crossprod(coeffs.courants[(ncol(T2)+1):(ncol(T2)+N),k]))/
        (N - 1/sigma[k] * sum(diag(solve(crossprod(designXi, W[,k]*designXi) +
                                           1/as.vector(sigma[k]) * diag(N)))))
    }
    
    # 4. Update
    ###########
    eta=cbind(T2,designXi)%*%coeffs.courants
    
    # etainf <- log(muinf)
    # indinf <- 1 * (eta < etainf)
    # eta <- eta * (1 - indinf) + etainf * indinf
    # mu = exp(eta+loffset)
    # Z = eta + (Y - mu)/mu
    # W <- mu
    if("poisson"%in%family)
    {
      etainf <- log(muinf)
      indinf<-1*(eta[,family=="poisson"]<etainf)
      eta[,family=="poisson"]<-eta[,family=="poisson"]*(1-indinf)+etainf*indinf
      if(is.null(loffset))
      {
        mu <- exp(eta[,family=="poisson"])
      }else{
        mu <- exp(eta[,family=="poisson"]+loffset)
      }
      Z[,family=="poisson"]<- eta[,family=="poisson"]+(Y[,family=="poisson"]-mu)/mu
      W[,family=="poisson"] <- mu
    }
    if("bernoulli"%in%family)
    {
      etainf <- log(muinf/(1-muinf))
      indinf<-1*(eta[,family=="bernoulli"]<etainf)
      eta[,family=="bernoulli"]<-eta[,family=="bernoulli"]*(1-indinf)+etainf*indinf
      indsup<-1*(eta[,family=="bernoulli"]>-etainf)
      eta[,family=="bernoulli"]<-eta[,family=="bernoulli"]*(1-indsup)-etainf*indsup
      mu <- exp(eta[,family=="bernoulli"])/(1+exp(eta[,family=="bernoulli"]))
      Z[,family=="bernoulli"] <-  eta[,family=="bernoulli"] + (Y[,family=="bernoulli"]-mu)/(mu*(1-mu))
      W[,family=="bernoulli"] <- mu*(1-mu)
    }
    if("binomial"%in%family)
    {
      etainf <- log(muinf/(1-muinf))
      indinf<-1*(eta[,family=="binomial"]<etainf)
      eta[,family=="binomial"]<-eta[,family=="binomial"]*(1-indinf)+etainf*indinf
      indsup<-1*(eta[,family=="binomial"]>-etainf)
      eta[,family=="binomial"]<-eta[,family=="binomial"]*(1-indsup)-etainf*indsup
      mu <- size*exp(eta[,family=="binomial"])/(1+exp(eta[,family=="binomial"]))
      Z[,family=="binomial"] <-  eta[,family=="binomial"] + (Y[,family=="binomial"]-mu)/((mu*(size-mu))/size)
      W[,family=="binomial"] <- (mu*(size-mu))/size
    }
    if("gaussian"%in%family){
      Z[,family=="gaussian"]<-Y[,family=="gaussian"]
    }
    
    
    
    
    deltaEta = sum(colSums((eta.old - eta)^2))
    eta.old <- eta
    deltaU <- sum((u.new-u.old)^2)
    u.old = u.new
    deltaSigma <-sum((sigma-sigma.old)^2)
    sigma.old = sigma
    compt <- compt+1
  }
  # ------------------------------------------------------------
  
  return(list(Z=Z,W=W,u=u.new,coefs=coeffs.courants,sigma=sigma,f=f))
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





# General function for k >= 1 component(s)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Y:       matrix of random responses
# X:       matrix of the normalized explanatory variables
# AX:      matrix of the additional covariates
# random:  vector giving the group of each unit (factor)
# loffset: matrix of size (nrow(Y), ncol(Y)) giving the log of the offset
#          relative to each unit
#k:        Number ok components to compute
# method:  structural relevance criterion
# return:  - u:        data frame of the loading vectors
#          - comp:     data frame of the computed components
#          - beta:     estimations of fixed effects
#          - sigma:    variance components
#          - lin.pred: linear predictors
#          - blup:     predictions for random effects
#          - other things useful for plots!




#' @title Function that fits the mixed-SCGLR model
#' @description Calculates the components to predict all the response variables.
#' @export kCompRand
#' @importFrom plsdepot plsreg2
#' @importFrom stats model.matrix cor var
#' @param Y the matrix of random responses
#' @param family a vector of character of the same length as the number of response variables: "bernoulli", "binomial", "poisson" or "gaussian" is allowed.
#' @param size describes the number of trials for the binomial dependent variables: a (number of observations * number of binomial response variables) matrix is expected.
#' @param X the matrix of the standardised explanatory variables
#' @param AX the matrix of the additional explanatory variables
#' @param random the vector giving the group of each unit (factor)
#' @param loffset a matrix of size (number of observations * number of Poisson response variables) giving the log of the offset associated with each observation
#' @param k number of components, default is one
#' @param init.sigma a vector giving the initial values of the variance components, default is rep(1, ncol(Y))
#' @param init.comp a character describing how the components (loadings-vectors) are inisialised in the PING algorithm: "pca" or "pls" is allowed
#' @param method Regularization criterion type: object of class "method.SCGLR"
#' built by function \code{\link{methodSR}}.
#' @return an object of the SCGLR class.
#' @examples \dontrun{
#' library(SCGLR)
#' # load sample data
#' data(dataGen)
#' k.opt=4
#' s.opt=0.1
#' l.opt=10
#' withRandom.opt=kCompRand(Y=dataGen$Y, family=rep("poisson", ncol(dataGen$Y)),
#'                         X=dataGen$X, AX=dataGen$AX,
#'                         random=dataGen$random, loffset=log(dataGen$offset), k=k.opt,
#'                         init.sigma = rep(1, ncol(dataGen$Y)), init.comp = "pca",
#'                         method=methodSR("vpi", l=l.opt, s=s.opt,
#'                                         maxiter=1000, epsilon=10^-6, bailout=1000))
#' plot(withRandom.opt, pred=TRUE, plane=c(1,2), title="Component plane (1,2)",
#'      threshold=0.7, covariates.alpha=0.4, predictors.labels.size=6)
#'}
kCompRand <- function(Y, family, size=NULL,
                      X, AX=NULL, random, loffset=NULL, k,
                      init.sigma = rep(1, ncol(Y)), init.comp = "pca",
                      method=methodSR("vpi", l=4, s=1/2,
                                      maxiter=1000, epsilon=10^-6, bailout=1000)){
  
  # Control
  if(is.null(colnames(Y))){
    colnames(Y) <- paste("y",1:ncol(Y),sep="")
  }
  if(is.null(colnames(X))){
    colnames(X) <- paste("x",1:ncol(X),sep="")
  }
  if(is.null(AX)){
    ncolAX <- 0
    colnamesAX <- NULL
  }else{
    ncolAX <- ncol(AX)
    if(is.null(colnames(AX))){
      colnamesAX <- paste("ax",1:ncol(AX),sep="")
    }else{
      colnamesAX <- colnames(AX)
    }
  }
  
  
  
  # First component
  no <- nrow(X)
  po <- ncol(X)
  # Init resACP
  resACP = NULL
  if(po>no){
    utilPCA = eigen(crossprod(x = X)/no,symmetric = TRUE)
    pourcent.variance = utilPCA$values/sum(utilPCA$values)
    nbcomp = max( which(c(pourcent.variance > (1/po))==TRUE ))
    Xpc = X%*%utilPCA$vectors[,1:nbcomp]
    #Stock
    resACP = list(utilPCA=utilPCA, pourcent.variance=pourcent.variance, nbcomp=nbcomp, Xpc=Xpc)
    if(init.comp=="pca"){u <- eigen(crossprod(x = Xpc)/no,symmetric = TRUE)$vector[,1]}
    if(init.comp=="pls"){u <- as.vector(plsdepot::plsreg2(Xpc, Y, comps = 2, crosval = FALSE)$raw.wgs[,1])}
  }else{
  if(init.comp=="pca"){u <- eigen(crossprod(x = X)/no,symmetric = TRUE)$vector[,1]}
  if(init.comp=="pls"){u <- as.vector(plsdepot::plsreg2(X, Y, comps = 2, crosval = FALSE)$raw.wgs[,1])}
  }
  # resACP in oneCompRand
  tmp <- oneCompRand(Y=Y, family=family, size=size,
                     X=X,AX=AX,u=u,random=random,loffset=loffset,
                     init.sigma = init.sigma, resACP=resACP, method=method)
  uMat <- matrix(tmp$u,ncol=1)
  fMat <- matrix(tmp$f,ncol=1)
  
  
  # Higher rank component
  if(k>1){
    for(kk in seq(k-1)){
      if(po>no){
        tmpX <- Xpc-fMat%*%(solve(crossprod(fMat),crossprod(fMat,Xpc)))
        if(init.comp=="pca"){
          u <- eigen(crossprod(x = tmpX)/no,symmetric = TRUE)$vector[,1]
          out_h <- mixed.hFunct(Z=tmp$Z,X=X,AX=AX,W=tmp$W,u=tmp$u,method=method,resACP=resACP)
        }
        if(init.comp=="pls"){
          u.initial <- as.vector(plsdepot::plsreg2(tmpX, Y, comps = 2, crosval = FALSE)$raw.wgs[,1])
          out_h <- mixed.hFunct(Z=tmp$Z,X=X,AX=AX,W=tmp$W,u=u.initial,method=method, resACP=resACP)
        }
      }else{
        tmpX <- X-fMat%*%(solve(crossprod(fMat),crossprod(fMat,X)))
      if(init.comp=="pca"){
        u <- eigen(crossprod(x = tmpX)/no,symmetric = TRUE)$vector[,1]
        out_h <- mixed.hFunct(Z=tmp$Z,X=X,AX=AX,W=tmp$W,u=tmp$u,method=method, resACP=resACP)
      }
      if(init.comp=="pls"){
        u.initial <- as.vector(plsdepot::plsreg2(tmpX, Y, comps = 2, crosval = FALSE)$raw.wgs[,1])
        out_h <- mixed.hFunct(Z=tmp$Z,X=X,AX=AX,W=tmp$W,u=u.initial,method=method, resACP=resACP)
      }
      }
      
      if(po>no){
        C <- (crossprod(Xpc,fMat))
        proj_C_ortho <- out_h$gradh - C%*%solve(crossprod(C),crossprod(C,out_h$gradh))
        u <- c(proj_C_ortho / sqrt(sum(proj_C_ortho^2)))
      }else{
        C <- (crossprod(X,fMat))#/nrow(X))
        proj_C_ortho <- out_h$gradh - C%*%solve(crossprod(C),crossprod(C,out_h$gradh))
        u <- c(proj_C_ortho / sqrt(sum(proj_C_ortho^2)))
      }
      # resACP in oneCompRand
      tmp <- oneCompRand(Y=Y,family=family,size=size,
                         X=X,AX=AX,u=u,fMat=fMat,random=random,loffset=loffset,
                         init.sigma=init.sigma, method=method, resACP=resACP)
      fMat <- cbind(fMat,tmp$f)
      uMat <- cbind(uMat,tmp$u)
    }
  }
  
  colnames(fMat) <- paste("scr",1:ncol(fMat),sep="")
  colnames(uMat) <- paste("u",1:ncol(uMat),sep="")
  coefs <- tmp$coefs
  rownames(coefs)[c(1,(ncolAX+(2:(k+1))))] <- c("intercept",colnames(fMat))
  sigma <- tmp$sigma
  names(sigma) <- colnames(Y)
  
  
  # BLUE, BLUP, and LINEAR PREDICTOR
  designXi <- as.matrix(model.matrix(~ factor(random)-1))
  invsqrtm <- diag(apply(X,2,var))*(no-1)/no
  centerx <- apply(X,2,mean)
  if(po>no){
    beta0 <- coefs[1,] - t(centerx)%*%invsqrtm%*%utilPCA$vectors[,1:nbcomp]%*%uMat%*%coefs[(ncolAX+(2:(k+1))),]
    beta <- invsqrtm%*%utilPCA$vectors[,1:nbcomp]%*%uMat%*%coefs[(ncolAX+(2:(k+1))),]
  }else{
    beta0 <- coefs[1,] - t(centerx)%*%invsqrtm%*%uMat%*%coefs[(ncolAX+(2:(k+1))),]
    beta <- invsqrtm%*%uMat%*%coefs[(ncolAX+(2:(k+1))),]
  }
  beta.coefs <- rbind(beta0,beta)
  #browser()
  if(!is.null(AX)){
    beta <- as.data.frame(rbind(beta.coefs, as.matrix(coefs[2:(ncolAX+1),,drop=F])))
  }else{
    beta <- as.data.frame(beta.coefs)
  }
  rownames(beta) <- c("Intercept",colnames(X),colnamesAX)
  xi <- coefs[(ncolAX+(k+2)):nrow(coefs),]
  if(!is.null(AX)){
    lin.pred <- cbind(1,as.matrix(X),as.matrix(AX))%*%as.matrix(beta)
    lin.pred.cond <- cbind(1,as.matrix(X),as.matrix(AX),as.matrix(designXi))%*%as.matrix(rbind(beta,xi))
  }else{
    lin.pred <- cbind(1,as.matrix(X))%*%as.matrix(beta)
    lin.pred.cond <- cbind(1,as.matrix(X),as.matrix(designXi))%*%as.matrix(rbind(beta,xi))
  }
  inertia <- cor(as.matrix(X), fMat)
  inertia <- inertia^2
  inertia <- colMeans(inertia)
  names(inertia) <- colnames(fMat)
  
  
  # OUTPUTS
  out <- list(
    call=NULL,
    u=as.data.frame(uMat),
    comp=as.data.frame(fMat),
    compr=as.data.frame(fMat),
    gamma=coefs,#[1:(ncolAX+k+1),],
    beta=beta,
    sigma=sigma,
    lin.pred=lin.pred,
    lin.pred.cond=lin.pred.cond,
    blup = xi,
    xFactors=NULL,
    xNumeric=X,
    inertia=inertia,
    invsqrtm=invsqrtm,
    centerx=centerx,
    Z=tmp$Z,
    W=tmp$W
  )
  class(out) <- "SCGLR"
  return(out)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


