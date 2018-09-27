#' @title Function that fits the theme model
#' @description Calculates the components to predict all the dependent variables.
#' @export
#' @param formula an object of class "\code{\link[=multivariateFormula]{MultivariateFormula}}" (or one that can be coerced
#' to that class): a symbolic description of the model to be fitted. The details
#' of model specification are given under Details.
#' @param data data frame.
#' @param H vector of R integer. Number of components to keep for each theme
#' @param family a vector of character of the same length as the number of dependent variables:
#' "bernoulli", "binomial", "poisson" or "gaussian" is allowed.
#' @param size describes the number of trials for the binomial dependent variables.
#' A (number of statistical units * number of binomial dependent variables) matrix is expected.
#' @param weights weights on individuals (not available for now)
#' @param offset used for the poisson dependent variables.
#' A vector or a matrix of size: number of observations * number of Poisson dependent variables is expected.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param na.action a function which indicates what should happen when the data contain NAs. The default is set to \code{na.omit}.
#' @param crit a list of two elements : maxit and tol, describing respectively the maximum number of iterations and
#' the tolerance convergence criterion for the Fisher scoring algorithm. Default is set to 50 and 10e-6 respectively.
#' @param method structural relevance criterion. Object of class "method.SCGLR"
#' built by  \code{\link{methodSR}} for Structural Relevance.
#' @param st logical (FALSE) theme build and fit order. TRUE means random, FALSE means sequential (T1, ..., Tr)
#' @details
#' Models for theme are specified symbolically. A model as the form \code{response ~ terms} where \code{response}
#' is the numeric response vector and terms is a series of R themes composed of
#' predictors. Themes are separated by  "|" (pipe) and are composed. ... Y1+Y2+...
#' ~ X11+X12+...+X1_  | X21+X22+... | ...+X1_+...  | A1+A2+... See \code{\link{multivariateFormula}}.
#' @return a list of SCGLRTHM class. Each element is a SCGLR object
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
#' form <- multivariateFormula(ny,nx1,nx2,A=c("geology"))
#' fam <- rep("poisson",length(ny))
#' testthm <-scglrTheme(form,data=genus,H=c(2,2),family=fam,offset = genus$surface)
#' plot(testthm)
#' }

scglrTheme <- function(formula, data, H, family, size = NULL, weights = NULL,
                  offset = NULL, subset = NULL, na.action = na.omit, crit = list(), method = methodSR(),st=FALSE){
  # TODO data vs subset, na.action, drop.levels, ...

  if(!inherits(formula,"MultivariateFormula"))
    formula <- multivariateFormula(formula,data=data)

  additional=formula$additional

  # check data
  if(!inherits(data, "data.frame"))
    stop("data must be compatible with a data.frame!")
  data_size <- nrow(data)

  # check crit
  crit <- do.call("critConvergence", crit)

  # extract left-hand side (Y)
  # Y is a formula of the form ~...
  if(length(formula)[1] != 1)
    stop("Left hand side part of formula (Y) must have ONE part!")
  theme_Y <- stats::terms(formula, lhs=1, rhs=0)
  Y_vars <- all.vars(theme_Y)

  # check family
  if(!is.character(family)||any(!(family%in%c("gaussian","poisson","bernoulli","binomial"))))
    stop("Expecting character vector of gaussian, poisson, bernoulli or binomial for family")
  if(!(length(family)%in%c(1,length(Y_vars))))
    stop("Length of family must be equal to one or number of Y variables")
  if(length(family)==1)
    family <- rep(family,length(Y_vars))

  # check and preprocess size parameter
  binomial_count <- sum(family=="binomial")
  if(!is.null(size)&&(binomial_count>0)) {
    if(is.vector(size)) {
      if(sum(family=="binomial")>1)
        message("Assuming that each binomial variable has same size!")
      size <- matrix(size, data_size, binomial_count)
    } else {
      size <- as.matrix(size)
    }
  }

  # check and preprocess offset parameter
  poisson_count <- sum(family=="poisson")
  if(!is.null(offset)&&(poisson_count>0)) {
    if(is.vector(offset)) {
      offset <- matrix(offset, data_size, poisson_count)
    } else {
      offset <- as.matrix(offset)
    }
  }

  # TODO check and preprocess weights parameter
  if(!is.null(weights))
    weights <- as.matrix(weights)

  # check and process themes #######################################################################

  # check part counts
  if(length(formula)[2] < 1+additional)
    if(additional) {
      stop("Right hand side part of formula with additional variables must have at least TWO parts!")
    } else {
      stop("Right hand side part of formula must have at least ONE part!")
    }

  # theme count
  theme_R <- length(formula)[2] - additional

  # check H (number of components to keep per theme)
  H <- as.integer(H)
  if(length(H) != theme_R)
    stop(sprintf("H must be a vector of R integers! (R=number of themes=%i for this call)", theme_R))
  if(any(H<0))
    stop("H must be positive integers")
  if(sum(H)==0)
    stop("H must contain at least one non null value")

  # extract right-hand side (X)
  # X is a list of R formulae of the form ~...
  # As a convenience we also keep value of current 'r' and 'Hr' as attributes
  theme_X <- lapply(1:theme_R, function(r) structure(stats::terms(formula, lhs=0, rhs=r), oldr=r, label=sprintf("T%s",r)))

  # name themes
  theme_labels <- sprintf("T%s", 1:theme_R)

  # removes themes when no components are requested (Hr=0)
  theme_X[which(H==0)] <- NULL
  theme_labels <- theme_labels[which(H>0)]
  H <- H[which(H>0)]
  theme_R = length(H)
  for(r in 1:theme_R) {
    attr(theme_X[[r]],"r") <- r
    attr(theme_X[[r]],"Hr") <- H[r]
  }

  # extract additional variables (A)
  # A is a formula of the form ~...
  if(additional) {
    theme_A <- stats::terms(formula, lhs=0, rhs=length(formula)[[2]])
  } else {
    theme_A <- NULL
  }

  # extract vars
  data_vars <- names(data)
  X_vars <- unique(unlist(lapply(theme_X, all.vars)))
  A_vars <- all.vars(theme_A)
  YXA_vars <-unique(c(Y_vars, X_vars, A_vars))

  # check if all variables can be found within data
  missing_vars <- YXA_vars[!YXA_vars %in% data_vars]
  if(length(missing_vars))
    stop("Some variable(s) where not found in data! '", paste(missing_vars, collapse="', '"),"'")

  # check that Y and X+A variable do not overlap
  error_vars <- Y_vars[Y_vars %in% c(X_vars, A_vars)]
  if(length(error_vars))
    stop("LHS and RHS variables must be different! '", paste(error_vars, collapse="', '"),"'")

  # check that A variables do not overlap with X
  error_vars <- A_vars[A_vars %in% X_vars]
  if(length(error_vars))
    stop("Additional variables must be different of theme variables! '", paste(error_vars, collapse="', '"),"'")

  # build data

  data <- data[,YXA_vars]

  # subset it if needed
  if(!is.null(subset)) {
    data <- data[subset,]
    offset <- offset[subset,]
    size <- size[subset,]
    weights <- weights[subset,]
    data_size <- nrow(data)
  }

  # build T (list of R Tr)
  theme_T <- lapply(theme_X, function(Xr) {
    r <- attr(Xr, "r")
    Hr <- attr(Xr, "Hr")
    oldr <- attr(Xr, "oldr")

    # build model frame
    Tr <- stats::model.frame(Xr, data)

    # get beta vars

    beta_vars <- sprintf("%s_%s", theme_labels[r], colnames(model.matrix(Xr, data[1,]))[-1])
    # check Hr value
    if(Hr > ncol(Tr))
      stop(sprintf("Number of required components (%i) for theme must be at most number of predictors (%i)!", Hr, ncol(Tr)),
           "\n   in theme '", Xr, "'",
           "\n   that expand to '", paste(names(Tr), collapse=" + "), "'")

        # extract Hr principal components
    if(has.categorical(Tr)) {
      Tr <- ade4::dudi.mix(Tr, scannf=FALSE, nf=Hr)$l1
    } else {
      Tr <- ade4::dudi.pca(Tr, scannf=FALSE, nf=Hr)$l1
    }

    # rename
    names(Tr) <- sprintf("%s_sc%i", theme_labels[r], 1:ncol(Tr))

    return(structure(Tr, r=r, beta_vars=beta_vars, oldr=oldr))
  })

  # get all beta vars
  TA_beta_vars <- unlist(lapply(theme_T, function(t) attr(t, "beta_vars")))
  if(additional) {
    A_beta_vars <- sprintf("A_%s", colnames(model.matrix(theme_A, data[1,]))[-1])
    TA_beta_vars <- c(TA_beta_vars, A_beta_vars)
  } else {
    A_beta_vars <- NULL
  }

  # build work dataset from data and Tr
  work_data <- do.call("cbind", list(data, theme_T))

  # init stop values
  tol1 <- rep(Inf, theme_R)
  tol2 <- Inf
  gamma_old <- matrix(Inf, 1+sum(H)+length(A_beta_vars), length(Y_vars))

  # init result
  out <-list(themes=list())

  # main optimization loop ####
  iter <- 0
  repeat {
    iter <- iter+1

    index <- NULL

    # main loop over themes ####

    if(st) {
      rand_themeX <- sample(theme_X)
    } else {
      rand_themeX <- theme_X
    }

    for(Xr in rand_themeX) {
      r <- attr(Xr, "r")
      Hr <- attr(Xr, "Hr")

      # get additional variables and T (except Tr) component variables
      TA_vars <- c(do.call("c", lapply(theme_T[-r], names)), all.vars(theme_A))
      if(length(TA_vars)==0)
        TA_vars <- NULL
      #TA <- formula(paste("~", paste(TA_vars, collapse="+")))

      # build formula
      scglr_formula <- multivariateFormula(Y_vars, attr(Xr, "term.labels"),A=TA_vars)
      #scglr_formula <- build_YXA(theme_Y, Xr, TA)

      tryCatch({
          # call scglr

          result <- scglr(scglr_formula, data=work_data, K=Hr, family=family,
                        offset=offset, size=size, weights=weights,
                        na.action = na.action, crit=crit, method=method)
        },
        convergence_failed = function(e) stop(
          "Convergence failed in theme ", Xr, " during ",iter," iteration",
          "\n  Try to increase 'l' value in methodSR (current=", method$l, ")",
          "\n  original message was: ", e,
          call.=FALSE)
      )
      # rename component names

      colnames(result$comp) <- sprintf("%s_%s", theme_labels[r], colnames(result$comp))
      colnames(result$compr) <- sprintf("%s_%s", theme_labels[r], colnames(result$compr))
      result$label <- theme_labels[r]

      # update result
      out$themes[[r]] <- result

      # update stop condition over components
      if(Hr>1){
        f_old <- as.matrix(work_data[, names(theme_T[[r]])])
        f_old <- apply(f_old,2,function(x)x/(sqrt(c(crossprod(x)/data_size))))
        f_new <- as.matrix(result$compr[, 1:Hr])
       } else {
        f_old <- c(work_data[, names(theme_T[[r]])])
        f_old <- f_old/(sqrt(c(crossprod(f_old)/data_size)))
        f_new <- c(result$compr[, 1:Hr])
       }
      tol1[r] <- sum(1 - diag(crossprod(f_old, f_new)/data_size)^2)
      # update Tr with Hr first components of scglr result
      work_data[, names(theme_T[[r]])] <- result$comp[, 1:Hr]#normalisation
     # index <- c(index, rep(r,Hr))

    } # end of main loop over themes

    # stop condition over gamma
    tol2 <- sum((gamma_old-result$gamma)^2)

    gamma_old <- result$gamma

    # check convergence
    if(all(c(tol1, tol2)<crit$tol) || (iter>crit$maxit))
      break;
  } # end of main optimization loop


  #  gamma <- result$gamma
  #
  # # reorder rows
  # if(length(H)>1) {
  #   if(additional) {
  #     gamma <- gamma[c(1,(H[length(H)]+2):(sum(H)+1), 2:(H[length(H)]+1), (sum(H)+2):nrow(gamma)),]
  #   } else {
  #     gamma <- gamma[c(1,(H[length(H)]+2):(sum(H)+1), 2:(H[length(H)]+1)),]
  #   }
  # }
  # rename rows
  # sc <- str_detect(rownames(gamma),"^sc")

  #refit model on components non centre
  comp <- do.call("cbind",lapply(out$themes, function(x) x$comp))
  x <- data.frame(comp,data[,A_vars])
  colnames(x) <- c(colnames(comp),A_vars)
  fit <- multivariateGlm.fit(Y=data[,Y_vars],comp=x,family=family, size = size, offset=offset)
  gamma <- sapply(fit,function(x) coef(x))

  #rownames(gamma)[sc] <- sprintf("%s_sc%s", theme_labels[r], 1:Hr)
  # if(additional) {
  #  rownames(gamma)[(sum(H)+2):nrow(gamma)] <- sprintf("A_%s", rownames(gamma)[(sum(H)+2):nrow(gamma)])
  # }
  # rownames(test) <- rownames(gamma)
  # gamma <- test


  index <- rep(1:theme_R,times=H)
  gammacr <- by(gamma[2:(sum(H)+1),,drop=FALSE],index,list)
  beta0 <- gamma[1,] - colSums(matrix(unlist(lapply(1:length(H),function(k) t(unlist(out$themes[[k]]$centerx))%*%as.matrix(out$themes[[k]]$invsqrtm)%*%as.matrix(out$themes[[k]]$u)%*%as.matrix(gammacr[[k]]))),ncol=length(Y_vars),byrow=T))
  beta <- do.call(rbind,lapply(1:theme_R,function(k) out$themes[[k]]$invsqrtm%*%as.matrix(out$themes[[k]]$u)%*%as.matrix(gammacr[[k]])))


  beta <- rbind(as.vector(beta0),as.matrix(beta))
  if(additional){
    beta <- rbind(beta,as.matrix(gamma[-c(1:(sum(H)+1)),,drop=FALSE]))
  }
  rownames(beta) <- c("(intercept)", TA_beta_vars)

  # compute global linear predictors
  comp <- cbind(1,as.matrix(do.call("cbind",lapply(out$themes, function(x) x$comp))))
  if(additional) {
    comp <- cbind(comp, model.matrix(theme_A, data)[,-1])
  }
  lin.pred.global <- multivariatePredictGlm(comp, family=family, beta=gamma, offset=offset)

  # update individual theme gamma and beta
  for(t in seq(1, length(out))) {
    par <- out$themes[[t]]$gamma[1:(H[t]+1),]
    out$themes[[t]]$lin.pred <- multivariatePredictGlm(cbind(1,as.matrix(out$themes[[t]]$compr)), family=family, beta=par, offset=offset)
    out$themes[[t]]$gamma <- gamma
    out$themes[[t]]$beta <- beta
  }
  # add them to result
  out$gamma <- gamma
  out$beta <- beta
  out$lin.pred.global <- lin.pred.global
  class(out) <- "SCGLRTHM"
  return(out)
}
