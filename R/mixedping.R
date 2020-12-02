###################################################################
#            The PROJECTED ITERATED NORMED GRADIENT               #
#         F. mortier, C. Trottier, G. Cornu and X. Bry            #
###################################################################



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Z:      matrix of the working variables
# X:      matrix of the normalized explanatory variables
# AX:     matrix of the additional covariates
# W:      matrix of weights
# u:      vector of initial loadings
# method: structural relevance criterion
# return: the FINAL updated loading vectors
mixed.ping <- function(Z,X,AX,W,F,u,method,resACP) {
  u_cur <- u
  h_cur <- mixed.hFunct(Z=Z,X=X,AX=AX,W=W,u=u_cur,method=method,resACP=resACP)$h
  u_new <- mixed.update_u(Z=Z,X=X,AX=AX,W=W,F=F,u=u_cur,method=method,resACP=resACP)
  h_new <- mixed.hFunct(Z=Z,X=X,AX=AX,W=W,u=u_new,method=method,resACP=resACP)$h
  ing_iter <- 0
  while((abs((h_new-h_cur)/h_cur)>method$epsilon)&(ing_iter<=method$maxiter)){
    u_cur <- u_new
    h_cur <- h_new
    u_new <- mixed.update_u(Z=Z,X=X,AX=AX,W=W,F=F,u=u_cur,method=method,resACP=resACP)
    h_new <- mixed.hFunct(Z=Z,X=X,AX=AX,W=W,u=u_new,method=method,resACP=resACP)$h
    ing_iter <- ing_iter+1
  }
  return(u_new)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Z:      matrix of the working variables
# X:      matrix of the normalized explanatory variables
# AX:     matrix of the additional covariates
# W:      matrix of weights
# F:      matrix of the previously calculated components
# u:      vector of initial loadings
# method: structural relevance criterion
# return: the CURRENT updated loading vectors
mixed.update_u <- function(Z,X,AX,W,F,u,method,resACP){
  out_h <- mixed.hFunct(Z=Z,X=X,AX=AX,W=W,u=u,method=method,resACP=resACP)
  if(is.null(F)) {
    m <- out_h$gradh/sqrt(sum(out_h$gradh^2))
  } else {
    no <- nrow(X)
    po <- ncol(X)
    if(po>no){
      utilPCA = resACP$utilPCA#eigen(crossprod(x = X)/no,symmetric = TRUE)
      pourcent.variance = resACP$pourcent.variance#utilPCA$values/sum(utilPCA$values)
      nbcomp = resACP$nbcomp#max( which(c(pourcent.variance > (1/po))==TRUE ))
      Xpc = resACP$Xpc#X%*%utilPCA$vectors[,1:nbcomp]
      
      C <- (crossprod(Xpc,F))#/nrow(X))
      proj_C_ortho <- out_h$gradh - C%*%solve(crossprod(C),crossprod(C,out_h$gradh))
      m <- c(proj_C_ortho / sqrt(sum(proj_C_ortho^2)))
    }else{
      C <- (crossprod(X,F))#/nrow(X))
      proj_C_ortho <- out_h$gradh - C%*%solve(crossprod(C),crossprod(C,out_h$gradh))
      m <- c(proj_C_ortho / sqrt(sum(proj_C_ortho^2)))
    }
  }
  h_m <- mixed.hFunct(Z=Z,X=X,AX=AX,W=W,u=m,method=method,resACP=resACP)$h
  h_u <- out_h$h
  k <- 1
  while((h_m<h_u)&(k<method$bailout)){
    m <- c(u)+m
    m <- m/sqrt(sum(m^2))
    h_m <- mixed.hFunct(Z=Z,X=X,AX=AX,W=W,u=m,method=method,resACP=resACP)$h
    k <- k+1
  }
  u_new <- m
  return(u_new)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Z:      matrix of the working variables
# X:      matrix of the normalized explanatory variables
# AX:     matrix of the additional covariates
# W:      matrix of weights
# u:      vector of initial loadings
# method: structural relevance criterion
# return: a list containing
#             - the value of h = (1-s)*log(psi) + s*log(phi)
#               with psi = GOF = Goodness-Of-Fit and phi = SR = Structural Relevance
#             - the values of phi, grad.phi, psi, grad.psi, and grad.h)
mixed.hFunct<- function(Z,X,AX,W,u,method,resACP){
  no <- nrow(X)
  po <- ncol(X)
  if(po>no){
    # utilPCA = eigen(crossprod(x = X)/no,symmetric = TRUE)
    # pourcent.variance = utilPCA$values/sum(utilPCA$values)
    # nbcomp = max( which(c(pourcent.variance > (1/po))==TRUE ))
    # Xpc = X%*%utilPCA$vectors[,1:nbcomp]
    utilPCA = resACP$utilPCA
    pourcent.variance = resACP$pourcent.variance
    nbcomp = resACP$nbcomp
    Xpc = resACP$Xpc
    
    f <- c(Xpc%*%u)
  }else{
    f <- c(X%*%u)
  }
  psi <- 0
  gradpsi <- rep(0,length(u))
  if(!is.null(AX)){
    for(k in 1:ncol(Z)){
      AXtWkAX <- crossprod(AX,W[,k]*AX)
      projWkfAX <- c(AX%*%solve(AXtWkAX,crossprod(AX,W[,k]*f)))
      projWkforthoAX <- f - projWkfAX
      Zk <- mixed.wtScale(Z[,k],W[,k])#Z[,k] - sum(W[,k]*Z[,k])
      WZk <- W[,k]*Zk
      projWkZAX <- AX%*%solve(AXtWkAX,crossprod(AX,WZk))
      # calculate  psi
      scalsqpfZ <- sum(c(projWkforthoAX)*WZk)^2
      scalsqpfpf <- sum(c(projWkforthoAX)^2*W[,k])
      term1psi <- sum(scalsqpfZ/(scalsqpfpf))
      term2psi <- sum(WZk*projWkZAX)
      psi <- psi+term1psi+term2psi
      #calculate grad.psi
      PiorthoPrimeWkZ <- WZk -  W[,k]*AX%*%solve(AXtWkAX,crossprod(AX,WZk))
      if(po>no){
        XprimeprojorthoWZ <- crossprod(Xpc,PiorthoPrimeWkZ)
        term1 <- c(XprimeprojorthoWZ%*%crossprod(XprimeprojorthoWZ,u))/(scalsqpfpf)
        WprojWkOrthof <- W[,k]*projWkforthoAX
        term2 <-  scalsqpfZ*c(crossprod(Xpc,WprojWkOrthof))/(scalsqpfpf^2)
        gradpsi <- gradpsi +(term1-term2)
      }else{
        XprimeprojorthoWZ <- crossprod(X,PiorthoPrimeWkZ)
        term1 <- c(XprimeprojorthoWZ%*%crossprod(XprimeprojorthoWZ,u))/(scalsqpfpf)
        WprojWkOrthof <- W[,k]*projWkforthoAX
        term2 <-  scalsqpfZ*c(crossprod(X,WprojWkOrthof))/(scalsqpfpf^2)
        gradpsi <- gradpsi +(term1-term2)
      }
    }
    gradpsi <- 2*gradpsi
  }else{
    for(k in 1:ncol(Z)){
      Zk <- mixed.wtScale(Z[,k],W[,k])#Z[,k] - sum(W[,k]*Z[,k])
      WZk <- W[,k]*Zk
      scalsqpfZ <- sum(c(f)*WZk)^2
      scalsqpfpf <- sum(c(f)^2*W[,k])
      # calculate  psi
      psi <- psi+sum(scalsqpfZ/(scalsqpfpf))
      #calculate grad.psi
      if(po>no){
        XprimeWZ <- crossprod(Xpc,WZk) #X'W_k Z_k
        term1 <- c(XprimeWZ%*%crossprod(XprimeWZ,u))/(scalsqpfpf)
        term2 <-  scalsqpfZ*c(crossprod(Xpc,W[,k]*f))/(scalsqpfpf^2)
        gradpsi <- gradpsi +(term1-term2)
      }else{
        XprimeWZ <- crossprod(X,WZk) #X'W_k Z_k
        term1 <- c(XprimeWZ%*%crossprod(XprimeWZ,u))/(scalsqpfpf)
        term2 <-  scalsqpfZ*c(crossprod(X,W[,k]*f))/(scalsqpfpf^2)
        gradpsi <- gradpsi +(term1-term2)
      }
    }
    gradpsi <- 2*gradpsi
  }
  n <- nrow(X)
  # calculate phi when method="cv" (Component Variance)
  if(method$phi=="cv") {
    phi <- c(crossprod(f))/n
    # calculate grad.phi
    if(po>no){
      gradphi <- c(2*crossprod(Xpc,f/n))
    }else{
      gradphi <- c(2*crossprod(X,f/n))
    }
  } else {
    ### calculate phi when method="vpi" (Variable Powered Inertia) with l>=1
    scalsqfX <- colSums(f*X/n)
    phi <- (sum((scalsqfX^2)^method$l))^(1/method$l)
    # calculate grad.phi
    if(po>no){
      XpcWX <- crossprod(Xpc, X)
      gradphi <- 2*phi^(1-method$l)*rowSums(XpcWX%*%diag( (scalsqfX^2)^method$l * scalsqfX^(-1) ))
    }else{
      XtWX <- crossprod(X)/n
      gradphi <- 2*phi^(1-method$l)*rowSums(XtWX%*%diag( (scalsqfX^2)^method$l * scalsqfX^(-1) ))
      #gradphi <- 2*phi^(1-method$l)*rowSums(XtWX%*%diag(scalsqfX)^(2*method$l-1))
    }
  }
  # s in [0..1]
  # h = log(psi)+method$s*log(phi)
  # gradh=gradpsi/psi+method$s*gradphi/phi
  h = (1-method$s)*log(psi)+method$s*log(phi)
  gradh=(1-method$s)*gradpsi/psi+method$s*gradphi/phi
  return(list(h=h, gradh=gradh,psi=psi,gradpsi=gradpsi,phi=phi,gradphi=gradphi))
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # method for STRUCTURAL RELEVANCE criterion
# # phi: character string describing structural relevance
# #        used in the regularization process
# #        Allowed values are "vpi" for Variable Powered Inertia
# #        and "cv" for Component Variance. Default to "vpi".
# # l:   numeric argument (>1) tuning the importance of variable bundle locality.
# # s:   numeric argument (in [0,1]) tuning the strength of structural relevance
# #        with respect to goodness of fit.
# # maxiter:  integer for maximum number of iterations of \code{SR} function
# # epsilon:  positive convergence threshold
# # bailout:  integer argument
# methodSR <- function(phi="vpi",l=1,s=1/2,maxiter=1000,epsilon=1e-6,bailout=10){
#   # check arguments
#   if(!(phi %in% c("vpi","cv")))
#     stop("phi should be \"vpi\" or \"cv\"")
#   if(!is.numeric(l) || l<1)
#     stop("l must be greater than 1")
#   if(!is.numeric(s) || s<0 || s>1)
#     stop("s must be between 0 and 1")
#   if(!is.numeric(maxiter) || maxiter<1)
#     stop("maxiter must be an integer greater than 1")
#   if(!is.numeric(epsilon) || epsilon<=0)
#     stop("epsilon must be a positive numeric")
#   if(!is.numeric(bailout) || bailout<1)
#     stop("bailout must be an integer greater than 1")
#
#   structure(list(
#     method="sr",
#     phi=phi,
#     l=l,
#     s=s,
#     maxiter=maxiter,
#     epsilon=epsilon,
#     bailout=bailout#1000
#   ),
#   class="method.SCGLR",
#   description="Method iterative normed gradient (ING) for Structural Relevance"
#   )
# }
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# auxiliary functions

############################
# weighted scale function
# x:      vector to scale
# w:      weight
# return: scaled vector
mixed.wtScale <-function(x,w){
  xc=x-sum(w*x)
  v=sum(xc^2*w)
  xcr=xc/sqrt(v)
  return(xcr)
}
############################


############################
# x:      vector to center
# w:      weight
# return: centered vector
mixed.wtCenter=function(x,w){
  xc=x-sum(w*x)
  return(xc)
}
############################


############################
mixed.checkLossFunction <- function(type){
  if(!type %in% c("auc","likelihood","aic","aicc","bic","mspe"))
    stop("Unknown loss function!")
}

# returns string w/o leading or trailing whitespace
mixed.trim <- function(x){
  gsub("^\\s+|\\s+$", "", x)}
############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

