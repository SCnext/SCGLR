
#' @title Formula construction
#' @description Helper function for building multivariate scglr formula.
#'
#' NOTE: Interactions involving factors are not allowed for now.
#' For interactions between two quantitative variables, use \code{I(x*y)} as usual.
#' @export
#' @param Y a formula or a vector of character containing the names of the dependent variables.
#' @param X a vector of character containing the names of the covariates (X) involved in the components or a list of it.
#' @param ... additional groups of covariates (theme)
#' @param A a vector of character containing the names of the additional covariates.
#' @param additional logical (if A is not provided, should we consider last X to be additional covariates)
#' @param data a data frame against which formula's variable will be checked
#' @return an object of class \code{MultivariateFormula, Formula, formula} with additional attributes: Y, X, A, X_vars, Y_vars,A_vars,XA_vars, YXA_vars, additional
#' @details
#' If Y is given as a formula, groups of covariates must be separated by \code{|} (pipes). To declare last
#' group as additional covariates, one can use \code{||} (double pipes) as last group separator or set
#' \code{additional} parameter as \code{TRUE}.
#' @examples \dontrun{
#' # build multivariate formula
#' ny <- c("y1","y2")
#' nx1 <- c("x11","x12")
#' nx2 <- c("x21","x22")
#' nadd <- c("add1","add2")
#' form <- multivariateFormula(ny,nx1,nx2,nadd,additional=T)
#' form2 <- multivariateFormula(ny,list(nx1,nx2,nadd),additional=T)
#' form3 <- multivariateFormula(ny,list(nx1,nx2),A=nadd)
#' form4 <- multivariateFormula(y1+y2~x11+x12|x21+x22||add1+add2)
#' # Print formulas
#' form
#' form2
#' form3
#' }
multivariateFormula <- function(Y, X=NULL, ..., A=NULL, additional=NULL, data=NULL) {
   if(inherits(Y,"formula")) {
    if(any(!is.null(X),!is.null(A))) {
      X <- NULL
      A <- NULL
      warning("As Y is a formula, I'm ignoring X and A")
    }

    # try to detect sugar operator || for additionals
    sugar_found <- FALSE
    walker <- function(expr) {
      if(is.name(expr)) {
        if(identical(expr,quote(`||`))) {
          if(sugar_found)
            stop("Only one || is allowed to mark additional covariate group!")
          sugar_found <<- TRUE
          expr <- quote(`|`)
        }
      } else if(is.call(expr)) {
        for (cc in seq_along(expr)) {
          if (is.name(expr[[cc]]) && expr[[cc]] == "")
            next
          expr[[cc]] <- walker(expr[[cc]])
        }
      }
      expr
    }
    Y <- walker(Y)

    if(sugar_found) {
      if(is.logical(additional)&&!additional)
        warning("|| found so ignoring additional parameter")
      additional <- TRUE
    }

    if(is.null(additional))
      additional=FALSE

    # split formula and check parts
    formula <- Formula(Y) # to handle multiple |
    l <- length(formula)
    if(l[1] != 1)
      stop("Left hand side part of formula (Y) must have ONE part!")

    # check part counts
    if(l[2] < 1+additional)
      if(additional) {
        stop("Right hand side part of formula with additional variables must have at least TWO parts!")
      } else {
        stop("Right hand side part of formula must have at least ONE part!")
      }

    Y <- stats::terms(formula,lhs=1,rhs=0)[[2]]
    X <- lapply(1:l[2], function(i) stats::terms(formula,lhs=0,rhs=i)[[2]])

  } else {

    if(!is.vector(Y)||!is.character(Y))
      stop("Y must be provided as vectors of response names")
    Y <- as.Formula(paste0("~",paste(Y,collapse="+")))[[2]]

    # parts are given as vector of names
    if(!is.list(X))
      X <- list(X)
    X <- c(X, list(...))

    # it must be a list of character vectors
    if(any(sapply(X,function(x) !is.vector(x)||!is.character(x))))
      stop("X and ... must be provided as vectors of covariate names")

    X <- lapply(X, function(x) as.Formula(paste0("~",paste(x,collapse="+")))[[2]])

    if(!is.null(A)) {
      A <- as.Formula(paste0("~",paste(A,collapse="+")))[[2]]
      additional <- TRUE
    }
  }

  # handle additional covariates removing them from X
  if(is.null(A)&&is.logical(additional)&&additional) {
    A <- unlist(X[[length(X)]])
    X <- X[-length(X)]
  }

  # give name to covariate groups
  names(X) <- paste0("T",1:length(X))

  # extract var names
  Y_vars <- all.vars(Y)
  X_vars <- unique(unlist(lapply(X, all.vars)))
  A_vars <- all.vars(A)
  XA_vars <- unique(c(X_vars,A_vars))
  YXA_vars <- unique(c(Y_vars,X_vars,A_vars))

  ## check consistency with data if provided
  # check if all variables can be found within data
  if(!is.null(data)) {
    data_vars <- names(data)
    missing_vars <- setdiff(YXA_vars,data_vars)
    if(length(missing_vars))
      stop("Some variable(s) where not found in data! '", paste(missing_vars, collapse="', '"),"'")
  }

  # check that Y and X+A variable do not overlap
  error_vars <- intersect(Y_vars,XA_vars)
  if(length(error_vars))
    stop("LHS and RHS variables must be different! '", paste(error_vars, collapse="', '"),"'")

  # check that Xs and A vars do not overlap with each other
  error_vars <- intersect(X_vars,A_vars)
  if(length(error_vars))
    stop("X and A variables must be different! '", paste(error_vars, collapse="', '"),"'")

  # formula builder from parts
  # first build formula . ~ 1 | 2 ...... | n
  # then replace numbered placeholders with corresponding parts and front dot with response expr
  cov <- c(X,A)
  walker2 <- function(expr) {
    if(is.atomic(expr)) {
      expr <- cov[[expr]]
    } else if(is.call(expr)) {
      for (cc in seq_along(expr)) {
        if (is.name(expr[[cc]]) && expr[[cc]] == "")
          next
        expr[[cc]] <- walker2(expr[[cc]])
      }
    }
    expr
  }
  formula <- as.formula(paste0(".~",paste0(1:length(cov),collapse="|")))
  formula <- walker2(formula)
  formula[[2]] <- Y
  formula <- Formula(formula)

  # document formula with collected metadata
  structure(
    formula,
    class=c("MultivariateFormula","Formula","formula"),
    Y=Y,
    X=X,
    A=A,
    additional=!is.null(A),
    Y_vars=Y_vars,
    X_vars=X_vars,
    A_vars=A_vars,
    XA_vars=XA_vars,
    YXA_vars=YXA_vars
  )
}

#' @title $ operator for multivariate formula
#' @export
#' @keywords internal
#' @param f formula
#' @param a attribute
#' @description
#' S3 helper function to retrieve attributes as if it was named values
#' @examples
#' frm <- multivariateFormula(y~x+z)
#' print(frm$Y)
'$.MultivariateFormula' <- function(f,a) {
  attr(f,a,exact = TRUE)
}

#' @title print a multivariate formula
#' @export
#' @keywords internal
#' @description
#' S3 helper function to print a multivariate formula
#' NB use $ semantic to retrieve metadata
#' @param x a formula
#' @param ... unused
print.MultivariateFormula <- function(x, ...) {

  deparse <- function(x) {
    paste0(base::deparse(x,60),collapse="\n")
  }
  
  cat("Multivariate formula \n   ",deparse(x),"\n")

  # response part
  cat("  Response: \n     Y   = ",deparse(x$Y),"\n")

  # covariates
  cat("  Covariates:\n")
  for(i in seq_along(x$X)) {
    cat("    ",paste0("T",i)," = ",deparse(x$X[[i]]),"\n")
  }

  # additional
  if(!is.null(x$A))
    cat("    ","A","  = ",deparse(x$A),"\n")

  cat("\n")
  invisible(x)
}

# print(multivariateFormula("y","1"))
# print(multivariateFormula(y1+y2~x1+x2|x3+x4|x5+x6*x7||a1+a2))
# print(multivariateFormula(c("y1","y2"),list(c("x1","x2"),c("x3","x4")),c("x5","x6*x7"),additional = TRUE))
# print(multivariateFormula(y1+y2~x1+x2|x3*x4|I(lon*lat)))
# print(multivariateFormula(Sepal.Length+Species~Petal.Length||Species,data=iris))

