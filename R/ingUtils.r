# @title Projected Iterated Normed Gradient
# @author F. mortier, C. Trottier, G. Cornu and X. Bry
# @param Z matrix of the working variables
# @param X matrix of the normalized covariates
# @param AX matrix of the suplementary covariates
# @param W matrix of weights
# @param u vector of loadings
# @param method structural relevance criterion. Object of class "method.SCGLR"  built by  \code{\link{methodSR}} for Structural Relevance.
# @return u_new the updated loading vectors

ping <- function(Z,X,AX,W,F,u,method) {
  u_cur <- u
  h_cur <- hFunct(Z=Z,X=X,AX=AX,W=W,u=u_cur,method=method)$h
  u_new <- update_u(Z=Z,X=X,AX=AX,W=W,F=F,u=u_cur,method=method)
  h_new <- hFunct(Z=Z,X=X,AX=AX,W=W,u=u_new,method=method)$h
  ing_iter <- 0
  while((abs((h_new-h_cur)/h_cur)>method$epsilon)&(ing_iter<=method$maxiter)){
    u_cur <- u_new
    h_cur <- h_new
    u_new <- update_u(Z=Z,X=X,AX=AX,W=W,F=F,u=u_cur,method=method)
    h_new <- hFunct(Z=Z,X=X,AX=AX,W=W,u=u_new,method=method)$h
    ing_iter <- ing_iter+1
  }
  return(u_new)
}


# ping <- function(Z,X,AX,W,F,u,method) {
#   u_cur <- u
#   h_cur <- hFunct(Z=Z,X=X,AX=AX,W=W,u=u_cur,method=method)$h
#
#   u_tmp <- update_u(Z=Z,X=X,AX=AX,W=W,F=F,u=u_cur,method=method)
#   h_tmp <- hFunct(Z=Z,X=X,AX=AX,W=W,u=u_tmp,method=method)$h
#
#   ing_iter <- 0
#
#   browser()
#   while((abs((h_tmp-h_cur)/h_cur)>method$epsilon)&(ing_iter<=method$maxiter)){
#     u_cur <- u_tmp
#     h_cur <- h_tmp
#     u_tmp <- update_u(Z=Z,X=X,AX=AX,W=W,F=F,u=u_cur,method=method)
#     h_tmp <- hFunct(Z=Z,X=X,AX=AX,W=W,u=u_new,method=method)$h
#     ing_iter <- ing_iter+1
#     print("j'y rentre")
#   }
#   print(c(h_cur,h_new))
#   return(u_new)
# }


update_u <- function(Z,X,AX,W,F,u,method){
  out_h <- hFunct(Z=Z,X=X,AX=AX,W=W,u=u,method=method)
  if(is.null(F)) {
    m <- out_h$gradh/sqrt(sum(out_h$gradh^2))
  } else {
    C <- (crossprod(X,F))#/nrow(X))
    proj_C_ortho <- out_h$gradh - C%*%solve(crossprod(C),crossprod(C,out_h$gradh))
    m <- c(proj_C_ortho / sqrt(sum(proj_C_ortho^2)))
  }
  h_m <- hFunct(Z=Z,X=X,AX=AX,W=W,u=m,method=method)$h
  h_u <- out_h$h
  k <- 1
  while((h_m<h_u)&(k<method$bailout)){
    m <- c(u)+m
    m <- m/sqrt(sum(m^2))
    h_m <- hFunct(Z=Z,X=X,AX=AX,W=W,u=m,method=method)$h
    k <- k+1
  }
  if(k>method$bailout) print("ARGH")
  u_new <- m
  return(u_new)
}

# updateGamma <- function(Z,X,AX,W,gamma,method){
#   g <- hFunct(Z=Z,X=X,AX=AX,W=W,u=gamma,method=method)
#   newgamma <- g$gradh/sqrt(sum(g$gradh^2))
#   tmp <- hFunct(Z=Z,X=X,AX=AX,W=W,u=newgamma,method=method)
#   k <- 1
#   while((tmp$h<g$h)&(k<method$bailout)){
#     newgamma <- gamma+newgamma
#     newgamma <- newgamma/sqrt(sum(newgamma^2))
#     tmp <- hFunct(Z=Z,X=X,AX=AX,W=W,u=newgamma,method=method)
#
#     k <- k+1
#   }
#   return(newgamma)
# }


hFunct<- function(Z,X,AX,W,u,method)
{
  f <- c(X%*%u)
  psi <- 0
  gradpsi <- rep(0,length(u))
  if(!is.null(AX)){
    for(k in 1:ncol(Z)){
      AXtWkAX <- crossprod(AX,W[,k]*AX)
      projWkfAX <- c(AX%*%solve(AXtWkAX,crossprod(AX,W[,k]*f)))
      projWkforthoAX <- f - projWkfAX
      Zk <- wtScale(Z[,k],W[,k])#Z[,k] - sum(W[,k]*Z[,k])
      WZk <- W[,k]*Zk
      projWkZAX <- AX%*%solve(AXtWkAX,crossprod(AX,WZk))
      #calcul de psi
      scalsqpfZ <- sum(c(projWkforthoAX)*WZk)^2
      scalsqpfpf <- sum(c(projWkforthoAX)^2*W[,k])
      term1psi <- sum(scalsqpfZ/(scalsqpfpf))
      term2psi <- sum(WZk*projWkZAX)
      psi <- psi+term1psi+term2psi
      #calcul de grad de psi
      PiorthoPrimeWkZ <- WZk -  W[,k]*AX%*%solve(AXtWkAX,crossprod(AX,WZk))
      XprimeprojorthoWZ <- crossprod(X,PiorthoPrimeWkZ)
      term1 <- c(XprimeprojorthoWZ%*%crossprod(XprimeprojorthoWZ,u))/(scalsqpfpf)

      WprojWkOrthof <- W[,k]*projWkforthoAX
      term2 <-  scalsqpfZ*c(crossprod(X,WprojWkOrthof))/(scalsqpfpf^2)
      gradpsi <- gradpsi +(term1-term2)
    }
    gradpsi <- 2*gradpsi
  }else{
    for(k in 1:ncol(Z)){
      Zk <- wtScale(Z[,k],W[,k])#Z[,k] - sum(W[,k]*Z[,k])
      WZk <- W[,k]*Zk
      scalsqpfZ <- sum(c(f)*WZk)^2
      scalsqpfpf <- sum(c(f)^2*W[,k])
      #calcul de psi
      psi <- psi+sum(scalsqpfZ/(scalsqpfpf))
      #calcul de grad de psi
      XprimeWZ <- crossprod(X,WZk) #X'W_k Z_k
      term1 <- c(XprimeWZ%*%crossprod(XprimeWZ,u))/(scalsqpfpf)
      term2 <-  scalsqpfZ*c(crossprod(X,W[,k]*f))/(scalsqpfpf^2)
      gradpsi <- gradpsi +(term1-term2)
    }
    gradpsi <- 2*gradpsi
  }
  n <- nrow(X)
  # calcul phi Component Variance: cv
  if(method$phi=="cv") {
    phi <- c(crossprod(f))/n
    # calcul grad phi
    gradphi <- c(2*crossprod(X,f/n))
  } else {
    ### autre calcul de phi avec l>=1 : vpi: Variable Powered Inertia
    scalsqfX <- colSums(f*X/n)
    XtWX <- crossprod(X)/n
    phi <- (sum((scalsqfX^2)^method$l))^(1/method$l)
    # calcul de grad phi
    gradphi <- 2*phi^(1-method$l)*rowSums(XtWX%*%diag(scalsqfX)^(2*method$l-1))
  }
  # calcul de h (s in R+)
  #h = log(psi)+method$s*log(phi)
  #gradh=gradpsi/psi+method$s*gradphi/phi
  # calcul de h (s in [0..1])
  h = (1-method$s)*log(psi)+method$s*log(phi)
  gradh=(1-method$s)*gradpsi/psi+method$s*gradphi/phi
  return(list(h=h, gradh=gradh,psi=psi,gradpsi=gradpsi,phi=phi,gradphi=gradphi))
}

