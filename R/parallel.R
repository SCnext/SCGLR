##%######################################################%##
#                                                          #
####             PARALLEL HELPER FUNCTIONS              ####
#                                                          #
##%######################################################%##

quietRequire <- function(package) {
  suppressPackageStartupMessages(
    suppressMessages(suppressWarnings(
      require(package, character.only = TRUE)
    )))
}

# check if parallelization is available
hasParallel <- function() {
  res <- isNamespaceLoaded("future") &&
    !inherits(future::plan(), "uniprocess")
  if(res && !quietRequire("future.apply")) {
    warning("Good start but future.apply package is also required to allow parallelization!\n  install.packages(\"future.apply\")")
    res <- FALSE
  }
  res
}

# return a progressor depending of whether progressr package is loaded or not
getProgressor <- function(...) {
  if(isNamespaceLoaded("progressr")) {
    progressr::progressor(..., envir=parent.frame(1))
  } else {
    function(...) {}
  }
}

# get parallel or not version of apply like function
# if only fun is provided it will return corresponding function
# if no parallel version is found then fallback to base one
# if additionnal parameters are given then the call will be also performed
# eg:  res <- getParallel("lapply", 1:10, function(x) {x+1})
getParallel <- function(fun, ...) {
  #fun <- as.character(substitute(fun))
  dots <- list(...)
  
  # resolve function
  base_fun <- get(fun, envir = asNamespace("base"))
  if(hasParallel()) {
    fun <- get0(
      paste0("future_", fun), 
      envir = asNamespace("future.apply")
    )
    if(is.null(fun)) {
      warning("Parallel version of ", fun, " was not found! Fallback to base")
      fun <- base_fun
    }
    dots <- c(dots, future.seed=TRUE)
  } else { 
    fun <- base_fun
  }
  
  # apply it if requested otherwise return it
  if(...length()) {
    do.call(fun, dots)
  } else {
    invisible(fun)
  }
}
