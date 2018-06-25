# @title Function that calculates the current scglr-component
# @export
# @param X a matrix (n*p) containing covariates standardized
# @param Y a matrix (n*q) containing dependent variables
# @param AX a matrix containing additional covariates
# @param family a vector of charater of length q specifying the distributions of the responses.
# Bernoulli, binomial, Poisson and gaussian are allowed.
# @param size specifies the numbers of trials of the binomial responses.  A (n*qb) matrix is expected
#  for qb binomial variables.
# @param offset used for the poisson dependent variables.
# A vector or a matrix of size: number of observations * number of Poisson dependent variables is expected
# when calculating the component
# @param crit a list of two elements : maxit and tol, describing respectively the maximum number of iterations and
# the tolerance convergence criterion for the Fisher scoring algorithm. Default is set to 50 and 10e-6 respectively.
# @param method optimization algorithm used for loadings and components. Object of class "method.SCGLR"
# built by \code{\link{methodEigen}} or \code{\link{methodIng}}
# @return the unit vector of loadings associated with the current component,
# i.e. the coefficients of the regressors in the linear combination giving each component


# NB AX = cbind(F,AX)
# TODO renommer AX en FAX
oneComponent <- function(X, Y, AX, F,family, size = NULL, offset = NULL, crit, method) {
  ## cst control and iteration
  muinf <- 1e-05
  maxit <- crit$maxit
  tol <- crit$tol
  tol1 <- 1
  tol2 <- 1e10
  iter <- 0

  ## dimension
  n <- dim(X)[1]
  #p <- dim(X)[2]
  q <- dim(Y)[2]
  if (is.null(q))
    q <- 1

  # initialisation de Z et de W, de alpha, de u et de a
  if (!is.null(offset))
    loffset <- log(offset)

  ### Initialization, Z working variables
  mu0 <- apply(Y, 2, mean)
  mu0 <- matrix(mu0, n, q, byrow = TRUE)
  Z <- Y

  if ("bernoulli" %in% family) {
    tmu0 <- mu0[, family == "bernoulli"]
    Z[, family == "bernoulli"] <- log(tmu0/(1 - tmu0)) + (Y[, family == "bernoulli"] - tmu0)/(tmu0 *
                                                                                              (1 - tmu0))
  }
  if ("binomial" %in% family) {
    tmu0 <- mu0[, family == "binomial"]
    Z[, family == "binomial"] <- log(tmu0/(1 - tmu0)) + (Y[, family == "binomial"] - tmu0)/(tmu0 *
                                                                                            (1 - tmu0))
  }
  if ("poisson" %in% family) {
    tmu0 <- mu0[, family == "poisson"]
    if (is.null(offset)) {
      Z[, family == "poisson"] <- log(tmu0) + (Y[, family == "poisson"] - tmu0)/tmu0
    } else {
      Z[, family == "poisson"] <- log(tmu0) - loffset + (Y[, family == "poisson"] - tmu0)/tmu0
    }
  }
  if ("gaussian" %in% family) {
    Z[, family == "gaussian"] <- Y[, family == "gaussian"]
  }

  #prohection ortho de X sur F Control sur singularite du systeme

  if(is.null(F)){
    u <- eigen(crossprod(x = X)/n,symmetric = TRUE)$vector[,1]
  }else{
    tmpX <- X-F%*%(solve(crossprod(F),crossprod(F,X)))
    u <- eigen(crossprod(x = tmpX)/n,symmetric = TRUE)$vector[,1]
  }


  # calcul du premier f
  f <- X %*% u

  # initialisation des eta
  if (is.null(AX)) {
    reg <- cbind(1, f)
    sol <- solve(crossprod(reg), crossprod(reg, Z))
    eta <- apply(sol, 2, function(x) x[1] + X %*% (u * x[2]))
  } else {
    r <- dim(AX)[2]
    reg <- cbind(1, AX, f)
    sol <- solve(crossprod(reg), crossprod(reg, Z))
    eta <- apply(sol, 2, function(x) x[1] + AX %*% x[2:(r + 1)] + X %*% (u * x[r + 2]))
  }

  # Update initialization of Z and initialization of W Z <- eta
  W <- matrix(1, n, q)

  if ("bernoulli" %in% family) {
    etainf <- log(muinf/(1 - muinf))
    indinf <- 1 * (eta[, family == "bernoulli"] < etainf)
    eta[, family == "bernoulli"] <- eta[, family == "bernoulli"] * (1 - indinf) + etainf * indinf
    indsup <- 1 * (eta[, family == "bernoulli"] > -etainf)
    eta[, family == "bernoulli"] <- eta[, family == "bernoulli"] * (1 - indsup) - etainf * indsup
    mu <- exp(eta[, family == "bernoulli"])/(1 + exp(eta[, family == "bernoulli"]))
    Z[, family == "bernoulli"] <- eta[, family == "bernoulli"] + (Y[, family == "bernoulli"] - mu)/(mu *
                                                                                                    (1 - mu))
    W[, family == "bernoulli"] <- mu * (1 - mu)
  }
  if ("binomial" %in% family) {
    etainf <- log(muinf/(1 - muinf))
    indinf <- 1 * (eta[, family == "binomial"] < etainf)
    eta[, family == "binomial"] <- eta[, family == "binomial"] * (1 - indinf) + etainf * indinf
    indsup <- 1 * (eta[, family == "binomial"] > -etainf)
    eta[, family == "binomial"] <- eta[, family == "binomial"] * (1 - indsup) - etainf * indsup
    mu <- exp(eta[, family == "binomial"])/(1 + exp(eta[, family == "binomial"]))
    Z[, family == "binomial"] <- eta[, family == "binomial"] + (Y[, family == "binomial"] - mu)/(mu *
                                                                                                 (1 - mu))
    W[, family == "binomial"] <- mu * (1 - mu) * size
  }
  if ("poisson" %in% family) {
    etainf <- log(muinf)
    indinf <- 1 * (eta[, family == "poisson"] < etainf)
    eta[, family == "poisson"] <- eta[, family == "poisson"] * (1 - indinf) + etainf * indinf
    if (is.null(offset)) {
      mu <- exp(eta[, family == "poisson"])
    } else {
      mu <- exp(eta[, family == "poisson"] + loffset)
    }
    Z[, family == "poisson"] <- eta[, family == "poisson"] + (Y[, family == "poisson"] - mu)/mu
    W[, family == "poisson"] <- mu
  }
  if ("gaussian" %in% family) {
    Z[, family == "gaussian"] <- Y[, family == "gaussian"]
  }
  W <- apply(W, 2, function(x) x/sum(x))

  if(!is.null(F)){
    out_h <- hFunct(Z=Z,X=X,AX=AX,W=W,u=u,method=method)
    C <- (crossprod(X,F))#/nrow(X))
    proj_C_ortho <- out_h$gradh - C%*%solve(crossprod(C),crossprod(C,out_h$gradh))
    u <- c(proj_C_ortho / sqrt(sum(proj_C_ortho^2)))
  }


### Loop
  sol_old <- sol

  while (((tol1 > tol) || (tol2>tol)) && (iter < maxit)) {

    unew <- ping(Z = Z, X = X, AX = AX, W = W, F = F, u = u, method = method)
    fnew <- X %*% unew
     if (is.null(AX)) {
      reg <- cbind(1, fnew)
      for (j in seq(q)) {
        sol[,j] <- solve(crossprod(reg, W[, j] * reg), crossprod(reg, W[, j] * Z[, j]))
        eta[, j] <- sol[1,j] + fnew %*% sol[2,j]
      }
    } else {
      r <- dim(AX)[2]
      reg <- cbind(rep(1, n), AX, fnew)
      for (j in seq(q)) {
        sol[,j] <- solve(crossprod(reg, W[, j] * reg), crossprod(reg, W[, j] * Z[, j]))
        eta[, j] <- sol[1,j] + AX %*% sol[2:(r + 1),j] + fnew %*% sol[r + 2,j]
      }
    }

   # Update of Z and W Z <- eta
    if ("bernoulli" %in% family) {
      etainf <- log(muinf/(1 - muinf))
      indinf <- 1 * (eta[, family == "bernoulli"] < etainf)
      eta[, family == "bernoulli"] <- eta[, family == "bernoulli"] * (1 - indinf) + etainf * indinf
      indsup <- 1 * (eta[, family == "bernoulli"] > -etainf)
      eta[, family == "bernoulli"] <- eta[, family == "bernoulli"] * (1 - indsup) - etainf * indsup
      mu <- exp(eta[, family == "bernoulli"])/(1 + exp(eta[, family == "bernoulli"]))
      Z[, family == "bernoulli"] <- eta[, family == "bernoulli"] + (Y[, family == "bernoulli"] -
                                                                    mu)/(mu * (1 - mu))
      W[, family == "bernoulli"] <- mu * (1 - mu)
    }
    if ("binomial" %in% family) {
      etainf <- log(muinf/(1 - muinf))
      indinf <- 1 * (eta[, family == "binomial"] < etainf)
      eta[, family == "binomial"] <- eta[, family == "binomial"] * (1 - indinf) + etainf * indinf
      indsup <- 1 * (eta[, family == "binomial"] > -etainf)
      eta[, family == "binomial"] <- eta[, family == "binomial"] * (1 - indsup) - etainf * indsup
      mu <- exp(eta[, family == "binomial"])/(1 + exp(eta[, family == "binomial"]))
      Z[, family == "binomial"] <- eta[, family == "binomial"] + (Y[, family == "binomial"] - mu)/(mu *
                                                                                                   (1 - mu))
      W[, family == "binomial"] <- mu * (1 - mu) * size
    }
    if ("poisson" %in% family) {
      etainf <- log(muinf)
      indinf <- 1 * (eta[, family == "poisson"] < etainf)
      eta[, family == "poisson"] <- eta[, family == "poisson"] * (1 - indinf) + etainf * indinf
      if (is.null(offset)) {
        mu <- exp(eta[, family == "poisson"])
      } else {
        mu <- exp(eta[, family == "poisson"] + loffset)
      }
      Z[, family == "poisson"] <- eta[, family == "poisson"] + (Y[, family == "poisson"] - mu)/mu
      W[, family == "poisson"] <- mu
    }
    if ("gaussian" %in% family) {
      Z[, family == "gaussian"] <- Y[, family == "gaussian"]
    }

    W <- apply(W, 2, function(x) x/sum(x))
    f <- c(X %*% u)
    fnew <- c(fnew)
    # NB les 1/n et 1/srqt(n) se simplifient
    f <- f/sqrt(c(crossprod(f))) # divise par 1/sqrt(n)
    fnew <- fnew/sqrt(c(crossprod(fnew))) # divise par 1/sqrt(n)
    tol1 <- 1 - crossprod(f, fnew)^2 # divise par 1/n
    u <- unew
    tol2 <- mean((sol_old-sol)^2)
    sol_old<-sol
    iter <- iter + 1
  }

  if (iter == maxit) {
    # TODO verbose option
    if (FALSE)
      message("  Warning !!! max number of iterations in oneComponent tol=", tol1)
    return(FALSE)
  }
  if (tol1 <= tol) {
    # TODO verbose option
    if (FALSE)
      message("  Convergence in ", iter, " and tol=", tol)
  }
  return(list(u = u, gamma = sol))
}

